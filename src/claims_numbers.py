"""Regenerate the numbers in CLAIMS.md ("predicted beta" and "needed parameter vs bound" columns) by calling the
work/ scripts that produced them. Every measured input is read from data/ by those scripts.

Expected values live in src/claims_expected.csv (claim, key, the text as printed in CLAIMS.md, numeric value,
comparison rule). `compare()` checks (1) that the printed text still occurs in that claim's CLAIMS.md row, so
drift of CLAIMS.md is flagged, and (2) that the regenerated number rounds to the printed one.

Sources (scripts are imported, not copied; their stdout is captured and parsed where they only print):
  EC-01         work/20261120-007 common.Cosmo (data/-driven cosmology): |xi S| = 2 beta / int dt (LSS -> today)
  EC-02         work/20261120-007 ec02_nogo.main()
  EC-03, EC-03' work/20261120-008 variants.qcd_variants(); beta_scan.main() + variants.main() on a sub-grid
  EC-04         work/20261120-007 ec04_estimates.main()
  EC-05         work/20261120-007 ec05_loop.main()
  KK-02/03/04   work/20261120-009 kk_bounds.main() (results kept in memory, written to artifacts/ not to work/)
Files the scripts write (plots, csv) are redirected to artifacts/reproduce/ by patching the module's HERE.
"""
import contextlib
import csv
import importlib
import io
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

from .data import ARTIFACTS, ROOT

WORK = ROOT / "work"
OUTDIR = ARTIFACTS / "reproduce"
EXPECTED = Path(__file__).resolve().parent / "claims_expected.csv"
T007 = "20261120-007-theory-nogo-minimal"
T008 = "20261120-008-cosmo-torsion-axion"
T009 = "20261120-009-kk-bounds"


# ------------------------------------------------------------------------------------------ module loading

def load(folder, name):
    """Import work/<folder>/<name>.py. Work folders share module names (common.py), so modules from other work
    folders are dropped from sys.modules first; modules already imported keep their own references."""
    d = (WORK / folder).resolve()
    for k, m in list(sys.modules.items()):
        f = getattr(m, "__file__", None)
        if f:
            p = Path(f).resolve().parent
            if p.parent == WORK.resolve() and p != d:
                del sys.modules[k]
    old = list(sys.path)
    sys.path.insert(0, str(d))
    try:
        return importlib.import_module(name)
    finally:
        sys.path[:] = old


def run_main(mod, *args):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ret = mod.main(*args)
    return buf.getvalue(), ret


def num(pattern, text, group=1):
    m = re.search(pattern, text)
    if not m:
        raise ValueError(f"pattern not found in script output: {pattern}")
    return float(m.group(group))


# ------------------------------------------------------------------------------------------ per task

def ec01_ec02_ec04_ec05():
    out = {}
    common = load(T007, "common")
    c = common.load()
    cos = common.Cosmo(c)
    beta = common.beta_rad(c)
    out["EC-01.xiS_needed_eV"] = 2 * beta / cos.t_int(lambda a: 1.0)
    out["EC-01.S0_lab_eV"] = c["S0_lab_max"]

    txt, ok = run_main(load(T007, "ec02_nogo"))
    row = re.search(r"nu-A massless, xi95\s*\|[^|]*\|[^|]*\|\s*([\d.eE+-]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.eE+-]+)", txt)
    out["EC-02.beta_max_deg"] = float(row.group(1))
    out["EC-02.orders_short"] = float(row.group(2))
    out["EC-02.xi_needed"] = float(row.group(3))
    out["EC-02.xi_max_helioscope"] = num(r"helioscope bound via dual axion: xi < ([\d.eE+-]+)", txt)
    out["_EC-02.checks_ok"] = bool(ok)

    txt, ok = run_main(load(T007, "ec04_estimates"))
    out["EC-04.xi_max_CAST"] = num(r"xi < ([\d.eE+-]+)", txt)
    out["EC-04.xi_min_lab"] = num(r"xi > ~([\d.eE+-]+)", txt)
    out["EC-04.w_tuning"] = num(r"w = A/m\^2 >~ ([\d.eE+-]+)", txt)
    out["_EC-04.checks_ok"] = bool(ok)

    txt, ok = run_main(load(T007, "ec05_loop"))
    gen = re.search(r"nu-C generous N_eff\s+JK \(max\)\s*: .*?max beta = ([\d.eE+-]+) deg, shortfall ([\d.]+) orders", txt)
    out["EC-05.beta_max_deg"] = float(gen.group(1))
    out["EC-05.orders_short"] = float(gen.group(2))
    out["_EC-05.checks_ok"] = bool(ok)
    return out


# masses of beta_scan.M kept for EC-03': the plateau maximum of beta under the 5% cap lies at 5.6e-33 eV; for
# m > 1e-30 eV the cap is 1% (Rogers) and beta_max is lower, for m < 1e-34 eV F << 1. `full=True` runs all 45.
SCAN_RANGE = (1e-34, 1.01e-30)


def load_scan_csv(path):
    with open(path, newline="") as f:
        return {float(r["m"]): float(r["beta_maxabs"]) for r in csv.DictReader(f)}


def ec03(full=False):
    out = {}
    OUTDIR.mkdir(parents=True, exist_ok=True)
    variants = load(T008, "variants")
    with contextlib.redirect_stdout(io.StringIO()):
        rows, b_univ = variants.qcd_variants()
    masses = [r["m_eV"] for r in rows]
    out["EC-03.m_QCD_min_eV"] = min(masses)
    out["EC-03.m_QCD_max_eV"] = max(masses)
    out["EC-03.beta_osc_deg"] = max(r["beta_loc_amp_deg"] for r in rows)
    out["EC-03.orders_short"] = min(r["orders_short"] for r in rows)
    out["_EC-03.f_a_independent"] = max(r["beta_loc_amp_deg"] for r in rows) / min(r["beta_loc_amp_deg"] for r in rows) - 1 < 1e-9

    beta_scan = load(T008, "beta_scan")
    M = np.asarray(beta_scan.M)
    if not full:
        M = M[(M >= SCAN_RANGE[0]) & (M <= SCAN_RANGE[1])]
    beta_scan.M = M
    beta_scan.HERE = str(OUTDIR)
    variants.HERE = str(OUTDIR)
    t0 = time.time()
    run_main(beta_scan)
    txt, _ = run_main(variants)
    out["_EC-03'.scan_seconds"] = time.time() - t0
    out["_EC-03'.scan_points"] = len(M)
    committed = load_scan_csv(WORK / T008 / "scan.csv")
    regen = load_scan_csv(OUTDIR / "scan.csv")
    common_m = sorted(set(committed) & set(regen))
    out["_EC-03'.scan_vs_committed_maxreldiff"] = max(
        abs(regen[m] - committed[m]) / max(abs(committed[m]), 1e-30) for m in common_m)
    out["EC-03'.g_T_GeVinv"] = num(r"g_T = ([\d.eE+-]+) GeV\^-1, phi_i <= pi f/N_q", txt)
    m_best = num(r"phi_i <= pi f/N_q: max beta = [\d.]+ deg at m = ([\d.eE+-]+) eV", txt)
    # the printed beta_max is rounded to 3 decimals; recompute it unrounded from the same arrays variants.main uses
    sc = variants.load_scan()
    m, F, cap, phiE, _, _ = variants.required_g(sc)
    ok = m <= 1.5e-27
    g_T = variants.ALPHA / (2 * np.pi * (variants.F_CF / variants.NQ)) * variants.EN
    bT = 0.5 * g_T * F * np.minimum(np.pi * variants.F_CF / variants.NQ, phiE)
    j = int(np.argmax(np.where(ok, bT, -1)))
    assert abs(m[j] / m_best - 1) < 1e-2
    bmax = bT[j] / variants.DEG
    out["EC-03'.beta_max_deg"] = bmax
    out["EC-03'.m_best_eV"] = m[j]
    out["EC-03'.cap_at_best"] = cap[j]
    out["EC-03'.cap_needed_pct"] = 100 * cap[j] * (variants.BETA_DEG / bmax) ** 2
    # inputs for the fits (task 014): beta scales as sqrt(cap) while phi_E < pi f/N_q
    out["_EC-03'.phiE_GeV"] = phiE[j] / variants.GEV
    out["_EC-03'.phi_excursion_cap_GeV"] = np.pi * variants.F_CF / variants.NQ / variants.GEV
    out["_EC-03'.F"] = F[j]
    out["_EC-03'.g_T_eVinv"] = g_T
    return out


def kk():
    out = {}
    OUTDIR.mkdir(parents=True, exist_ok=True)
    kb = load(T009, "kk_bounds")
    kb.HERE = str(OUTDIR)
    kb.OUT.clear()
    txt, _ = run_main(kb)
    r = {k: v for k, v, _ in kb.OUT}
    out["_KK.checks_ok"] = "FAIL" not in txt
    out["KK-02.beta_max_k1_deg"] = r["KK02_k1_beta_max_deg"]
    out["KK-02.invL_min_GeV"] = r["KK02_k1_invL_min_GeV"]
    out["KK-02.invL_max_GeV"] = r["KK02_invL_max_GeV"]
    xs = []
    for (k, v, note) in kb.OUT:
        if k.startswith("KK02_xMl_window"):
            xs += [v, float(re.search(r"to ([\d.]+)", note).group(1))]
    out["KK-02.Ml_min"] = min(xs)
    out["KK-02.Ml_max"] = max(xs)
    out["KK-03.g_min_GeVinv"] = r["KK03_g_min_GeVinv"]
    out["KK-03.g_max_GeVinv"] = r["KK03_g_max_GeVinv"]
    out["KK-03.natural_excluded"] = float(r["KK03_natural_dNnu_BBN"] > 1e3 * r["dNnu_max_BBN"])
    out["_KK-03.natural_dNnu"] = r["KK03_natural_dNnu_BBN"]
    out["KK-04.gamma_minus_1"] = r["radion_gamma_minus_1"]
    out["KK-04.m_ISL_eV"] = r["radion_m_min_ISL_eV"]
    out["_KK-04.mISL_over_Hstar"] = r["radion_m_min_ISL_eV"] / r["H_star_eV"]
    return out


def regenerate(full_scan=True):
    vals, times = {}, {}
    for label, fn in (("007 (EC-01/02/04/05)", ec01_ec02_ec04_ec05), ("008 (EC-03, EC-03')", lambda: ec03(full_scan)),
                      ("009 (KK-02/03/04)", kk)):
        t0 = time.time()
        vals.update(fn())
        times[label] = time.time() - t0
    return vals, times


# ------------------------------------------------------------------------------------------ comparison

def claims_rows():
    rows = {}
    for line in (ROOT / "CLAIMS.md").read_text(encoding="utf-8").splitlines():
        if line.startswith("| ") and not line.startswith("| id") and not line.startswith("|--"):
            cid = line.split("|")[1].strip()
            rows[cid] = line
    return rows


def expected():
    with open(EXPECTED, newline="", encoding="utf-8") as f:
        return [r for r in csv.DictReader(f) if not r["claim"].startswith("#")]


def _match(v, e, rule, digits):
    """True if regenerated v agrees with printed e under the stated rounding."""
    if rule == "sig":  # same value at `digits` significant figures
        r = lambda x: float(f"{x:.{digits - 1}e}")
        return r(v) == r(e)
    if rule == "dec":  # same value at `digits` decimals
        return round(v, digits) == round(e, digits)
    if rule == "oom":  # same power of ten after rounding log10 ("~10^120")
        return round(math.log10(abs(v))) == round(math.log10(abs(e)))
    if rule == "floor":  # printed as ">= e"
        return math.floor(v) == e
    if rule == "exact":
        return v == e
    raise ValueError(rule)


def compare(vals):
    rows = claims_rows()
    table = []
    for x in expected():
        printed = x["claims_text"]
        in_claims = printed in rows.get(x["claim"], "")
        v = vals.get(f"{x['claim']}.{x['key']}")
        e = float(x["expected"])
        if v is None:
            status = "NOT REGENERATED"
        else:
            status = "match" if _match(v, e, x["rule"], int(x["digits"] or 0)) else "MISMATCH"
        if not in_claims:
            status += " (text not in CLAIMS.md)"
        superseded = x.get("status", "live").startswith("superseded")
        if superseded:  # still regenerated and reported, but not a live claim number: never a failure
            status = f"{status} [superseded]"
        table.append(dict(claim=x["claim"], quantity=x["quantity"], value=v, printed=printed, rule=f"{x['rule']}{x['digits']}",
                          status=status, superseded=superseded))
    return table


def fmt(v):
    if v is None:
        return "-"
    if v == 0 or 1e-3 <= abs(v) < 1e4:
        return f"{v:.4g}"
    return f"{v:.3e}"


def print_table(table):
    w = max(len(t["quantity"]) for t in table)
    print(f"{'claim':7s} | {'quantity':{w}s} | {'regenerated':>11s} | {'in CLAIMS.md':24s} | rule  | match")
    for t in table:
        print(f"{t['claim']:7s} | {t['quantity']:{w}s} | {fmt(t['value']):>11s} | {t['printed']:24s} | {t['rule']:5s} | {t['status']}")
    live = [t for t in table if not t["superseded"]]
    bad = [t for t in live if t["status"] != "match"]
    sup = [t for t in table if t["superseded"]]
    print(f"claim numbers (live): {len(live) - len(bad)} match, {len(bad)} not matching; "
          f"superseded (EC-03 QCD-mass reading, skeptic 011): {sum(t['status'].startswith('match') for t in sup)}/{len(sup)} "
          "still reproduced by 008's scripts")
    return bad
