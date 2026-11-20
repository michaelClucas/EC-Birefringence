"""Reproduce entry point (from the repo root, pinned, 1 thread):
    .venv/Scripts/python tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --log <file>

Steps:
  1. likelihood validation (task 010): published beta significances and SPT-3G results.
  2. claim numbers, regeneration (tasks 014, 017): calls the work/ scripts (src/claims_numbers.py).
  3. candidate fits (task 014): EC-03', KK-02, generic ALP vs null on isotropic beta and SPT A_CB. Output csv goes to
     artifacts/reproduce/, not to work/. The fit numbers quoted in CLAIMS.md are merged into the claim numbers.
  4. claim numbers vs CLAIMS.md: compared with src/claims_expected.csv, whose printed text is also checked to still
     occur in CLAIMS.md.
  4b. draft-only numbers (task 033): src/draft_expected.csv, numbers due in the RESULT draft but not in CLAIMS.md
     (029 Pauli bound, Cosmoglobe penalties). A wrong or missing value fails; text not (yet) in
     work/20261120-018-result-draft/RESULT_DRAFT.md is reported as "not in draft" and does not fail.
  5. fit summary table (model, beta_max, penalty vs free ALP, Delta chi2 vs null) for ACT DR6, ACT+WMAP+Planck (EK22)
     and ACT+Cosmoglobe WMAP+Planck.
  5b. only with --checks (task 027, skeptic 020 MINOR 12): the derivation/symbolic check scripts behind the
     no-number verdicts (tasks 003, 004, 009, 012, 016, 019, 022, and 025 if present), each in its own pinned
     process; fails if any exits nonzero or prints FAIL. Output files they write into work/ are restored.
  6. path scan: the log (--log) and the text outputs in artifacts/reproduce/ are scanned for absolute local paths;
     the repo root becomes <repo>, the home folder <home>, anything else absolute <abs-path>. stdout/stderr are
     scrubbed the same way as they are written, so a shell redirect is also path-free.
Options: --skip-validation (steps 2-6 only), --fast-scan (EC-03' scan on the 17-mass sub-grid instead of all 45),
         --log FILE (tee the output to FILE, relative to the current directory), --checks (step 5b).
Exit code is nonzero if any step fails or any live claim number does not match.
"""
import warnings

# candl has an invalid escape sequence; the SyntaxWarning text contains the absolute site-packages path
warnings.filterwarnings("ignore", category=SyntaxWarning)

import csv  # noqa: E402
import importlib.util  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

from . import claims_numbers as CN  # noqa: E402
from .data import ROOT  # noqa: E402

VALIDATE = ROOT / "work" / "20261120-010-data-likelihood" / "validate.py"
FIT = ROOT / "work" / "20261120-014-data-reproduce-fit" / "fit_candidates.py"
FIT_DATASETS = {"act_dr6": "ACT DR6", "act_wmap_planck": "ACT+WMAP+Planck",
                "act_cosmoglobe": "ACT+Cosmoglobe WMAP+Planck"}  # Cosmoglobe added by task 033 (skeptic 020 S4)


# ------------------------------------------------------------------------------------------ path scrubbing

def _path_regex(p):
    """Regex for path p with either separator, any case, also in /c/... form."""
    p = Path(p).resolve()
    parts = [re.escape(x) for x in p.parts[1:]]
    drive = p.drive[:1]
    head = rf"(?:{drive}:|/{drive})" if drive else ""
    return re.compile(head + r"[\\/]+" + r"[\\/]+".join(parts), re.IGNORECASE)


_SCRUB = [(_path_regex(ROOT), "<repo>"), (_path_regex(Path.home()), "<home>"),
          (re.compile(r"(?<![\w<])[A-Za-z]:[\\/][^\s'\"<>|,;)]*"), "<abs-path>"),
          (re.compile(r"(?<![\w<.])/(?:Users|home|[a-z]/Users)/[^\s'\"<>|,;)]*"), "<abs-path>")]


def scrub(text):
    n = 0
    for rx, rep in _SCRUB:
        text, k = rx.subn(rep, text)
        n += k
    return text, n


class _Scrubbed:
    """Text stream wrapper: scrubs absolute paths on write, optionally tees to a log file."""

    def __init__(self, stream, log=None):
        self._s, self._log, self.hits = stream, log, 0

    def write(self, text):
        text, n = scrub(text)
        self.hits += n
        if self._log is not None:
            self._log.write(text)
        return self._s.write(text)

    def flush(self):
        self._s.flush()
        if self._log is not None:
            self._log.flush()

    def __getattr__(self, name):
        return getattr(self._s, name)


def scan_files(paths):
    """Rewrite each text file with absolute paths replaced; returns (files scanned, replacements)."""
    nfiles = nrep = 0
    for p in paths:
        if not p.is_file():
            continue
        text = p.read_text(encoding="utf-8", errors="replace")
        new, n = scrub(text)
        nfiles += 1
        if n:
            p.write_text(new, encoding="utf-8")
            nrep += n
    return nfiles, nrep


# ------------------------------------------------------------------------------------------ steps

def _load(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def regenerate(full_scan):
    vals, times = CN.regenerate(full_scan=full_scan)
    for k, t in times.items():
        print(f"  {k}: {t:.1f} s")
    npts, diff = vals["_EC-03'.scan_points"], vals["_EC-03'.scan_vs_committed_maxreldiff"]
    print(f"  [INFO] EC-03' scan: {npts} masses, max rel. diff of beta/(g phi_i) vs committed "
          f"work/20261120-008 scan.csv = {diff:.1e}")
    return vals


def fits(vals):
    """Run 014's fits with output redirected to artifacts/reproduce; returns (failed, fit rows, claim values)."""
    fit = _load(FIT)
    CN.OUTDIR.mkdir(parents=True, exist_ok=True)
    fit.HERE = CN.OUTDIR
    fit.DATASETS = list(FIT_DATASETS)  # 014 fits ACT DR6 and EK22 only; the Cosmoglobe rows feed the draft numbers
    failed = bool(fit.main(vals))
    with open(CN.OUTDIR / "fit_results.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    caps = fit.ultralight_caps(fit.cosmo_inputs()["omch2"])
    return failed, rows, CN.fit_values(rows, caps)


def claim_table(vals):
    flags = {k: bool(v) for k, v in vals.items() if k.startswith("_") and str(v) in ("True", "False")}
    for k, v in flags.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k[1:]}")
    bad = CN.print_table(CN.compare(vals))
    print("### draft-only numbers vs the RESULT draft (src/draft_expected.csv; task 033)")
    bad_draft = CN.print_draft_table(CN.compare_draft(vals))
    return bool(bad or bad_draft or [k for k, v in flags.items() if not v])


def fit_summary(rows):
    for ds, label in FIT_DATASETS.items():
        sub = [r for r in rows if r["dataset"] == ds]
        print(f"{label}:")
        print(f"  {'model':28s} | {'beta_max (deg)':>14s} | {'penalty vs ALP':>14s} | {'dchi2 vs null':>13s}")
        for r in sub:
            free = r["model"].startswith("ALP")
            bm = "free" if free else f"{float(r['beta_max_deg']):.3f}"
            print(f"  {r['model']:28s} | {bm:>14s} | {float(r['penalty_vs_ALP']):14.2f} | "
                  f"{float(r['dchi2_vs_null']):13.2f}")


# symbolic/derivation checks behind the no-number verdicts (--checks). Task 025 is picked up if it has a run_checks.py.
CHECKS = ["20261120-003-theory-ec-photon/run_checks.py", "20261120-004-kk-reduction/run_checks.py",
          "20261120-009-kk-bounds/run_checks.py", "20261120-012-theory-shift-uholst/run_checks.py",
          "20261120-016-kk05-edm-cs/run_checks.py", "20261120-019-adjudicate-kk05-edm/adjudicate_edm.py",
          "20261120-022-theory-ec06/run_checks.py"]
OPTIONAL_CHECKS = ["20261120-025-*/run_checks.py"]
# outputs the check scripts write into their own work folder (run.log, csv, png); restored afterwards so --checks
# leaves work/ unchanged. Scripts and READMEs are never touched.
_OUTPUT_SUFFIXES = (".log", ".csv", ".png", ".json", ".txt")
_FAIL_RX = re.compile(r"^\s*(\[FAIL\]|FAIL\b)|\b[1-9]\d* FAIL|SOME FAIL", re.M)


def _snapshot(folder):
    return {p: p.read_bytes() for p in folder.rglob("*") if p.is_file() and p.suffix in _OUTPUT_SUFFIXES
            and "__pycache__" not in p.parts}


def _restore(folder, snap):
    n = 0
    for p in folder.rglob("*"):
        if not p.is_file() or p.suffix not in _OUTPUT_SUFFIXES or "__pycache__" in p.parts:
            continue
        if p not in snap:
            p.unlink()
            n += 1
        elif p.read_bytes() != snap[p]:
            p.write_bytes(snap[p])
            n += 1
    return n


def run_checks():
    """Run each check script in its own pinned process (tools/pin.py, 1 thread) from the repo root. Full output goes
    to artifacts/reproduce/checks_<task>.log. Fails if a script exits nonzero or prints a FAIL line/count."""
    work = ROOT / "work"
    scripts = [work / s for s in CHECKS] + sorted(p for g in OPTIONAL_CHECKS for p in work.glob(g))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", PYTHONIOENCODING="utf-8")
    py = [sys.executable, "-W", "ignore::SyntaxWarning"]
    CN.OUTDIR.mkdir(parents=True, exist_ok=True)
    bad, t_all = [], time.time()
    for s in scripts:
        rel = s.relative_to(ROOT).as_posix()
        if not s.is_file():
            print(f"  [FAIL] {rel}: missing")
            bad.append(rel)
            continue
        snap = _snapshot(s.parent)
        t0 = time.time()
        r = subprocess.run(py + [str(ROOT / "tools" / "pin.py")] + py + [str(s)], cwd=ROOT, env=env,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        dt = time.time() - t0
        restored = _restore(s.parent, snap)
        text, _ = scrub(r.stdout + r.stderr)
        (CN.OUTDIR / f"checks_{s.parent.name.split('-')[1]}.log").write_text(text, encoding="utf-8")
        npass = len(re.findall(r"^\s*(\[PASS\]|PASS\b)", text, re.M))
        fails = [ln.strip() for ln in text.splitlines() if _FAIL_RX.search(ln)]
        ok = r.returncode == 0 and not fails
        print(f"  [{'PASS' if ok else 'FAIL'}] {rel}: exit {r.returncode}, {npass} PASS lines, {len(fails)} FAIL lines, "
              f"{dt:.1f} s; {restored} output file(s) in its folder restored")
        for ln in fails[:10]:
            print(f"      {ln}")
        if not ok:
            bad.append(rel)
    print(f"  checks: {len(scripts) - len(bad)}/{len(scripts)} scripts pass ({time.time() - t_all:.0f} s)")
    return bad


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    args = sys.argv[1:]
    log_path = Path(args[args.index("--log") + 1]) if "--log" in args else None
    log = open(log_path, "w", encoding="utf-8") if log_path else None
    sys.stdout = out = _Scrubbed(sys.stdout, log)
    sys.stderr = err = _Scrubbed(sys.stderr, log)
    failed = []
    t0 = time.time()
    try:
        if "--skip-validation" not in args:
            print("### likelihood validation (20261120-010)")
            if _load(VALIDATE).main():
                failed.append("likelihood validation")
        print("### claim numbers: regeneration (20261120-014, -017)")
        vals = regenerate(full_scan="--fast-scan" not in args)
        print("### candidate fits (20261120-014; csv to artifacts/reproduce/)")
        fit_failed, rows, fit_vals = fits(vals)
        if fit_failed:
            failed.append("candidate fits")
        vals.update(fit_vals)
        print("### claim numbers vs CLAIMS.md (20261120-014, -017)")
        if claim_table(vals):
            failed.append("claim numbers")
        print("### fit summary (penalty = chi2_min(model) - chi2_min(free ALP); dchi2 = chi2_null - chi2_min)")
        fit_summary(rows)
        if "--checks" in args:
            print("### derivation checks behind the no-number verdicts (--checks; logs in artifacts/reproduce/)")
            if run_checks():
                failed.append("derivation checks")
    except Exception:
        # task 033: an exception used to print "reproduce: OK" and its traceback bypassed the path scrubbing
        import traceback
        traceback.print_exc()  # sys.stderr is the scrubbed wrapper here
        failed.append("exception")
    finally:
        print("### path scan")
        sys.stdout.flush()
        files = sorted(p for p in CN.OUTDIR.glob("*") if p.suffix in (".csv", ".txt", ".log", ".json"))
        nf, nr = scan_files(files)
        print(f"  scrubbed on write: {out.hits + err.hits} absolute path(s); artifacts/reproduce: {nf} text files, "
              f"{nr} replaced")
        print(f"### reproduce: {'OK' if not failed else f'FAILED: {failed}'} ({time.time() - t0:.0f} s)")
        sys.stdout.flush()
        if log is not None:
            log.close()
            sys.stdout, sys.stderr = out._s, err._s
            _, n = scan_files([log_path])
            print(f"  log rescan: {n} absolute path(s) replaced")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
