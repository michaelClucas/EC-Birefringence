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
  EC-03 (photon decoupling), KK-05   work/20261120-012 check_formfactor, bounds_uholst (top-level scripts: imported
                with stdout captured, values read from module attributes or parsed from the printed [INFO] lines)
  fit numbers (EC-03' 2.7% beta_max, penalties)   work/20261120-014 fit_candidates.main(), merged by src.reproduce
                via fit_values()
  EC-06, EC-06a-d   work/20261120-022 ec06_numbers (needed xi, GW bound, f(T) numbers, SNY 0+ form), -023 ec06_common
                (xi eps for p tracking H) and ec06_aniso (A_CB), -024 kaf_lowz (low-z xi eps limit, run from the repo
                root), -026 skep_beta (skeptic corrections: MZ GW side, visibility-averaged
                reionization share), -025 ec06c_numbers (EC-06c photon-mass window, EFT cutoff), sny_background
                (EC-06d xi and H(z*) on SNY cases III-V), check_sny_limit (EC-06d GR-regime xi). Table EC06 below,
                one line per number (tasks 027, 030).
  iter-8 skeptic numbers (task 033)   work/20261120-028 window_cutoff (EC-06c Lambda_eta, white-dwarf lower edge),
                025 sny_background (EC-06d H/H_LCDM at z = 10), 025 check_sny_limit (xi for 10% deviations),
                032 kk04_shortfall (KK-04 shortfall vs Cassini), 032 hz_bounds (EC-06d H(z*) ratios with radiation and
                shortfall vs the cited H(z*) bounds; replaces 025's no-radiation ratios); in EC06.
  iter-9/10 numbers (task 038)   work/20261120-031 mz_numbers (EC-06a Lambda_SC 3.6 keV, 9.4 orders below 10 TeV,
                xi >= 4.7e16 from k*), -034 skeptic_mz (EC-06a no-GW tracker <= 9 MeV, f(Lambda)FFtilde xi ratio 21),
                -035 ec06d_endpoint (EC-06d ground/L2 ratio), 032 kk04_shortfall (KK-04 brane / 2 sigma 3.6-4.1); in EC06.
  draft-only numbers (task 033)   src/draft_expected.csv: numbers quoted in the RESULT draft
                (work/20261120-018-result-draft/RESULT_DRAFT.md) but not in CLAIMS.md. Table DRAFT_NUMBERS below
                (029 pauli_bound) plus the Cosmoglobe fit penalties from fit_values(). compare_draft() checks each
                value and whether the printed text occurs in the draft ("not in draft" is reported, never a failure).
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
T012 = "20261120-012-theory-shift-uholst"
T016 = "20261120-016-kk05-edm-cs"
T022 = "20261120-022-theory-ec06"
T023 = "20261120-023-cosmo-ec06"
T024 = "20261120-024-lit-ec06-novelty"
T025 = "20261120-025-theory-ec06c-costs"
T026 = "20261120-026-skeptic-ec06"
T028 = "20261120-028-skeptic-ec06c-costs"
T029 = "20261120-029-final-skeptic-rev2"
T031 = "20261120-031-theory-mz-verdict"
T032 = "20261120-032-lit-bounds-fixes"
T034 = "20261120-034-skeptic-mz"
T035 = "20261120-035-final-skeptic-rev3"
DRAFT_EXPECTED = Path(__file__).resolve().parent / "draft_expected.csv"
DRAFT = WORK / "20261120-018-result-draft" / "RESULT_DRAFT.md"

# KK-05 numbers: one line each (task 016 may revise them; change the source here and the row in claims_expected.csv).
# key -> (work folder, module, getter on the imported module)
KK05 = {
    "c_min_MICROSCOPE": (T012, "bounds_uholst", lambda m: m.c_ep[0.39]),
    "c_min_clocks": (T012, "bounds_uholst", lambda m: m.c_clock_min),
    "c_max_EDM": (T016, "edm_bound", lambda m: m.c_edm),  # 012 value superseded (EDM/magnetic pieces swapped)
}


# ------------------------------------------------------------------------------------------ module loading

def _drop_other_work_modules(d):
    for k, m in list(sys.modules.items()):
        f = getattr(m, "__file__", None)
        if f:
            p = Path(f).resolve().parent
            if p.parent == WORK.resolve() and p != d:
                del sys.modules[k]


def load(folder, name):
    """Import work/<folder>/<name>.py. Work folders share module names (common.py), so modules from other work
    folders are dropped from sys.modules first; modules already imported keep their own references."""
    d = (WORK / folder).resolve()
    _drop_other_work_modules(d)
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


def exec_script(folder, name, here=None, cwd=None):
    """Run the top-level script work/<folder>/<name>.py in-process and return (namespace, stdout, exit code).
    `here`: directory the script sees as its own (__file__), so files it writes next to itself land there instead
    of in work/. `cwd`: working directory during the run (024's kaf_lowz.py reads data/ relative to the repo root)."""
    d = (WORK / folder).resolve()
    _drop_other_work_modules(d)
    src = d / f"{name}.py"
    g = {"__name__": f"_reproduce_{name}", "__file__": str(Path(here or d) / f"{name}.py")}
    old_path, old_argv = list(sys.path), list(sys.argv)
    sys.path.insert(0, str(d))
    sys.argv = [src.name]
    buf, code = io.StringIO(), 0
    try:
        with contextlib.redirect_stdout(buf), contextlib.chdir(cwd or Path.cwd()):
            try:
                exec(compile(src.read_text(encoding="utf-8"), str(src), "exec"), g)
            except SystemExit as e:
                code = e.code or 0
    finally:
        sys.path[:], sys.argv = old_path, old_argv
    return g, buf.getvalue(), code


def load_script(folder, name):
    """Import a work/ script that runs at import time; returns (module, captured stdout)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mod = load(folder, name)
    return mod, buf.getvalue()


def ec03_photon():
    """012 F5: g_eff/g_anom for a derivative coupling. Range = min/max of the values 012 prints; coefficient
    (m_phi/m_e)^2/coeff from its ratio() at a test q."""
    out = {}
    ff, txt = load_script(T012, "check_formfactor")
    r = [float(x) for x in re.findall(r"g_eff/g_anomaly at [^:\n]*: ([\d.eE+-]+)", txt)]
    if not r:
        raise ValueError("012 check_formfactor printed no g_eff/g_anomaly values")
    out["EC-03.g_ratio_min"] = min(r)
    out["EC-03.g_ratio_max"] = max(r)
    q = 1e-30
    out["EC-03.g_ratio_coeff"] = (q / ff.M_E) ** 2 / ff.ratio(q)
    out["_EC-03.012_formfactor_checks_ok"] = "[FAIL]" not in txt
    return out


def kk05():
    out, mods = {}, {}
    for key, (folder, name, get) in KK05.items():
        if (folder, name) not in mods:
            mod, txt = load_script(folder, name)
            mods[(folder, name)] = mod
            out[f"_KK-05.{name}_checks_ok"] = "[FAIL]" not in txt
        out[f"KK-05.{key}"] = float(get(mods[(folder, name)]))
    return out


# ------------------------------------------------------------------------------------------ EC-06 (tasks 022-024)
# Sources: each runs once per reproduce; getters below read from its namespace. Files the scripts write
# (ec06_backgrounds.csv, acb_estimate.csv) go to artifacts/reproduce/.

def _ec06_022():
    g, txt, _ = exec_script(T022, "ec06_numbers", here=OUTDIR)
    g["_checks_ok"] = "FAIL:" not in txt and re.search(r"== \d+ PASS, 0 FAIL ==", txt) is not None
    return g


def _ec06_023():
    c = load(T023, "ec06_common")
    return {"xieps_n0": c.BETA_ACT * c.DEG / c.beta_unit(c.ZSTAR, 0.0), "_checks_ok": True}


def _ec06_023_aniso():
    g, txt, _ = exec_script(T023, "ec06_aniso", here=OUTDIR)
    g["_checks_ok"] = all(g["checks"]) and "FAIL:" not in txt
    return g


def _ec06_024():
    g, txt, code = exec_script(T024, "kaf_lowz", cwd=ROOT)
    g["_txt"] = txt
    g["_checks_ok"] = code == 0 and txt.rstrip().endswith("PASS")
    return g


def _ec06_026():
    g, txt, code = exec_script(T026, "skep_beta", cwd=ROOT)  # read-only; loads data/ relative to the repo root
    g["_checks_ok"] = code == 0
    return g


def _ec06_025(name):
    """025 scripts share one common.py PASS/FAIL counter, so judge each by its own FAIL: lines and exit code.
    here=OUTDIR: sny_background.csv is written to artifacts/reproduce/, not work/."""
    def f():
        g, txt, code = exec_script(T025, name, here=OUTDIR, cwd=ROOT)
        g["_checks_ok"] = code == 0 and "FAIL:" not in txt and "PASS:" in txt
        return g
    return f


def _ec06_028():
    """028 window_cutoff.py: reads its own inputs_028.csv next to the script and writes no files, so it runs with
    its own folder as __file__ (no `here` redirect)."""
    g, txt, code = exec_script(T028, "window_cutoff", cwd=ROOT)
    g["_checks_ok"] = code == 0 and re.search(r"== \d+ PASS, 0 FAIL ==", txt) is not None
    return g


def _pauli_029():
    """029 pauli_bound.py: top-level asserts; read-only."""
    g, txt, code = exec_script(T029, "pauli_bound", cwd=ROOT)
    g["_checks_ok"] = code == 0 and txt.count("PASS:") == 2
    return g


def _kk04_032():
    """032 kk04_shortfall.py: |gamma - 1| of the radion (bulk matter) vs the Cassini 1 sigma error
    (data/kk_bounds_inputs.csv); read-only. 032's summary line has no closing '=='."""
    g, txt, code = exec_script(T032, "kk04_shortfall", cwd=ROOT)
    g["_checks_ok"] = code == 0 and "FAIL:" not in txt and re.search(r"== \d+ PASS, 0 FAIL", txt) is not None
    return g


def _hz_032():
    """032 hz_bounds.py: SNY H(z)/H_LCDM with radiation included, from 025's committed sny_background.csv (read
    only), and the shortfall in orders of delta H/H against the cited 95% bounds (data/expansion_history_bounds.csv)."""
    g, txt, code = exec_script(T032, "hz_bounds", cwd=ROOT)
    g["_checks_ok"] = code == 0 and "FAIL:" not in txt and "PASS:" in txt
    return g


def _mz_031():
    """031 mz_numbers.py: EC-06a numbers (tracker, GW and k* bounds, Lambda_SC); read-only, data/ via 031 common."""
    g, txt, code = exec_script(T031, "mz_numbers", cwd=ROOT)
    g["_checks_ok"] = code == 0 and "FAIL:" not in txt and re.search(r"== \d+ PASS, 0 FAIL ==", txt) is not None
    return g


def _mz_034():
    """034 skeptic_mz.py: independent EC-06a checks (own loader); read-only."""
    g, txt, code = exec_script(T034, "skeptic_mz", cwd=ROOT)
    g["_checks_ok"] = code == 0 and "FAIL:" not in txt and re.search(r"== \d+ PASS, 0 FAIL ==", txt) is not None
    return g


def _endpoint_035():
    """035 ec06d_endpoint.py: EC-06d observer-endpoint dependence; reads local_inputs.csv next to itself; read-only.
    It prints PASS/FAIL lines but always exits 0, so both PASS lines are required."""
    g, txt, code = exec_script(T035, "ec06d_endpoint", cwd=ROOT)
    g["_checks_ok"] = code == 0 and "FAIL" not in txt and txt.count("PASS:") == 2
    return g


def _kk04_alt(g):
    """032 kk04_shortfall: orders short for brane matter (1 sigma, 2 sigma) and bulk matter at 2 sigma."""
    lg = g["math"].log10
    bulk, brane = g["pred"]["bulk matter (omega = 0)"], g["pred"]["brane matter (alpha0^2 = 1/6)"]
    return [lg(bulk / g["b2"]), lg(brane / g["b1"]), lg(brane / g["b2"])]


def _hz_star(g):
    """032 rows at z*: (case, z, H/H_LCDM with radiation, bound, orders in H, orders in rho); 3 cases x 2 bounds."""
    rows = [o for o in g["out"] if abs(o[1] - g["zstar"]) < 1]
    if len({o[0] for o in rows}) != 3:
        raise ValueError(f"032 hz_bounds: {len(rows)} rows at z*, need all three SNY cases")
    return rows


def _sny_at(g, z):
    """H_SNY/H_LCDM at redshift z for every SNY case (025 sny_background.py rows; all three must reach z)."""
    rows = [r for r in g["res"] if abs(r["z"] - z) < 1e-9]
    if len(rows) != 3:
        raise ValueError(f"025 sny_background: {len(rows)} of 3 cases reached z = {z}")
    return [r["H_over_LCDM"] for r in rows]


def _sny_star(g):
    """SNY cases III-V rows at z* from 025 sny_background.py (all three must reach z*)."""
    rows = g["reached_star"]
    if len(rows) != 3:
        raise ValueError(f"025 sny_background: {len(rows)} of 3 cases reached z*")
    return rows


def _sny_coeff(g):
    """Coefficient c in Phi = c d ln(R + 6 mu/b)/dt from 022's sympy solve of SNY Eq. 46."""
    import sympy as sp
    R, t, mu, b = g["R"], g["t"], g["mu"], g["b"]
    return float(sp.simplify(g["Phi"] * (R + 6 * mu / b) / sp.diff(R, t)))


def _sny_range(g, i):
    """Phi/H today for SNY cases III-V: transcribed in 022's background table (SNY 0805.3834 Table 1), not computed."""
    row = next(r for r in g["rows"] if r[0].startswith("PGT 0+"))
    return float(re.search(r"([\d.]+) to ([\d.]+)", row[3]).group(i))


def _kaf_xieps(g, rid):
    return float(re.search(rf"^{re.escape(rid)}\s.*xi\*eps<\s*([\d.eE+-]+)", g["_txt"], re.M).group(1))


EC06_SOURCES = {"022": _ec06_022, "023": _ec06_023, "023-aniso": _ec06_023_aniso, "024": _ec06_024, "026": _ec06_026,
                "025-num": _ec06_025("ec06c_numbers"), "025-sny": _ec06_025("sny_background"),
                "025-lim": _ec06_025("check_sny_limit"), "028-win": _ec06_028, "029-pauli": _pauli_029,
                "032-kk04": _kk04_032, "032-hz": _hz_032, "031-mz": _mz_031, "034-mz": _mz_034,
                "035-end": _endpoint_035}

# One line per CLAIMS.md number: "claim.key" -> (source, getter on that source's namespace). Task 025 may revise
# EC-06c/EC-06d: change the getter here and the row in claims_expected.csv.
EC06 = {
    "EC-06.xieps_needed": ("023", lambda g: g["xieps_n0"]),
    "EC-06.xieps_max_lowz": ("024", lambda g: _kaf_xieps(g, "L-AFG10-UV")),
    "EC-06.ACB_max": ("023-aniso", lambda g: max(g["res"]["A"].max(), g["res"]["B"].max())),
    # skeptic 026: GW170817 applied on the c_T > 1 side (+7e-16), eps_S = 2 gamma/3 (022 used the -3e-15 side)
    "EC-06a.epsS_max_GW": ("026", lambda g: 2 * g["g_hi"] / 3),
    "EC-06a.xi_min_GW": ("026", lambda g: g["xieps"] / (2 * g["g_hi"] / 3)),
    "EC-06c.eps_T": ("022", lambda g: g["eps_T"](0.0)),
    "EC-06c.post_reion_min_pct": ("026", lambda g: 100 * g["wavg"] / g["math"].log(1 + g["zs"])),  # visibility-avg
    "EC-06c.post_reion_max_pct": ("022", lambda g: 100 * g["fr_re"]),  # at z_re
    "EC-06c.xi_needed": ("022", lambda g: g["xi_T"]),
    "EC-06c.xi_min_k1": ("022", lambda g: g["xi_k"]),
    "EC-06d.Phi_coeff": ("022", _sny_coeff),
    "EC-06d.PhiH_min": ("022", lambda g: _sny_range(g, 1)),
    "EC-06d.PhiH_max": ("022", lambda g: _sny_range(g, 2)),
    # task 025 (EC-06c repair (a) Stueckelberg photon; EC-06d at z*), added by task 030
    "EC-06c.mgamma_min_eV": ("025-num", lambda g: g["m_low"]),  # max(conversion A < 0.1, CFJ stability)
    "EC-06c.mgamma_max_eV": ("025-num", lambda g: g["m_up_proca"]),  # data/ec06c_costs_inputs.csv (hep-ph/0306245)
    "EC-06c.cutoff_min_GeV": ("025-num", lambda g: (g["m_low"] * g["MP"]**2 / g["xi"])**(1 / 3) / 1e9),
    "EC-06c.cutoff_max_GeV": ("025-num", lambda g: (g["m_up_proca"] * g["MP"]**2 / g["xi"])**(1 / 3) / 1e9),
    "EC-06d.xi_SNY_zstar": ("025-sny", lambda g: sum(r["xi_needed"] for r in _sny_star(g)) / 3),  # mean of III-V
    # task 032 (iter 9) re-sources the z* ratios with radiation in H_LCDM (025's 33-66 omitted it)
    "EC-06d.Hratio_min_zstar": ("032-hz", lambda g: min(o[2] for o in _hz_star(g))),
    "EC-06d.Hratio_max_zstar": ("032-hz", lambda g: max(o[2] for o in _hz_star(g))),
    "EC-06d.orders_short_Hzstar_min": ("032-hz", lambda g: min(o[4] for o in _hz_star(g))),  # weakest: dNeff, case V
    "EC-06d.orders_short_Hzstar_max": ("032-hz", lambda g: max(o[4] for o in _hz_star(g))),  # Omega_e, case III
    "EC-06d.xi_min_GRlike": ("025-lim", lambda g: g["xi_per"](1e-2)),  # |delta H^2/H^2 (z*)| <= 1e-2, mu = 1
    # skeptic 028 (iter 8), added by task 033. Lambda_eta = (64 pi alpha M_P^2 m^2/(3 xi^2))^(1/4) at alpha = 1, on
    # 025's window [m_low, m_up_proca] and at the PDG photon-mass limit (data/ec06c_costs_inputs.csv)
    "EC-06c.Lambda_eta_min_keV": ("028-win", lambda g: g["lam_eta"](g["_src"]("025-num")["m_low"], 1.0) / 1e3),
    "EC-06c.Lambda_eta_max_keV": ("028-win", lambda g: g["lam_eta"](g["m_proca"], 1.0) / 1e3),
    "EC-06c.Lambda_eta_PDG_MeV": ("028-win", lambda g: g["lam_eta"](g["m_pdg"], 1.0) / 1e6),
    # typical 0.6 M_sun white dwarf, A < 0.1 (025's criterion); Sirius B gives 1e-23 eV, A < 1 gives 2.5e-25 eV
    "EC-06c.mgamma_min_WD_eV": ("028-win", lambda g: g["res"]["WD 0.6 Msun"][0.1]),
    "EC-06d.Hratio_max_z10": ("025-sny", lambda g: max(_sny_at(g, 10.0))),  # case III
    "EC-06d.xi_min_10pct": ("025-lim", lambda g: g["xi_per"](1e-1)),  # |delta H^2/H^2 (z*)| <= 0.1, mu = 1
    # skeptic 029 item 3, computed by task 032: orders between the bulk-matter radion |gamma - 1| = 1/2 and Cassini
    "KK-04.orders_short_Cassini": ("032-kk04", lambda g: g["math"].log10(g["pred"]["bulk matter (omega = 0)"] / g["b1"])),
    # task 038 (iter 11): numbers 035 listed as not regenerated
    "KK-04.orders_short_alt_min": ("032-kk04", lambda g: min(_kk04_alt(g))),  # brane matter, 2 sigma (3.63)
    "KK-04.orders_short_alt_max": ("032-kk04", lambda g: max(_kk04_alt(g))),  # brane matter, 1 sigma (4.09)
    # EC-06a (031): theta-exchange strong-coupling scale at xi(GW speed), PDG photon mass, 1 - c_s = Omega_m today
    "EC-06a.LambdaSC_GW_keV": ("031-mz", lambda g: g["lam_sc"](g["xi_gw"], g["ecc"]["m_gamma_PDG"], g["one_m_cs"]) / 1e3),
    "EC-06a.orders_below_Esafe": ("031-mz", lambda g: g["math"].log10(
        g["E_safe"] / g["lam_sc"](g["xi_gw"], g["ecc"]["m_gamma_PDG"], g["one_m_cs"]))),
    "EC-06a.xi_min_kstar": ("031-mz", lambda g: g["xi_k"]),  # k* = H/eps_S > k_GW(24 Hz) on the tracker
    # skeptic 034: plain tracker at the N_eff limit (no GW bound), PDG mass; f(Lambda)FFtilde needs xi x J/J_alt
    "EC-06a.LambdaSC_noGW_MeV": ("034-mz", lambda g: g["lam_sc"](g["xi_noGW"], g["mP"], g["omc_mean"]) / 1e6),
    "EC-06a.xi_ratio_vac": ("034-mz", lambda g: g["J"] / g["Jalt"]),
    # skeptic 035: beta/xi ratio, ground (air at 5.2 km) vs L2 observer endpoint, heavy-0+ branch (mu = 1, delta = 1%)
    "EC-06d.endpoint_ratio_ground_L2": ("035-end", lambda g: g["r"]),
}

# Numbers due in the RESULT draft but not in CLAIMS.md: "claim.key" -> (source, getter). Expected values in
# src/draft_expected.csv. The Cosmoglobe penalties come from fit_values() (src.reproduce merges them).
DRAFT_NUMBERS = {
    # 029 Pauli bound: any fermionic axial source in minimal EC, g = 100 helicity states carrying all of rho_LCDM
    "EC-02.orders_short_pauli": ("029-pauli", lambda g: g["math"].log10(g["beta_obs"] / g["beta_max_deg"](100, 1.0))),
    "EC-02.c_vector_needed": ("029-pauli", lambda g: g["c_need"]),  # vector-current (baryon) source, xi = 1
}
# EC-03'.m_best_eV (draft only, task 038) comes from ec03(): grid mass of beta_max under the 5% cap.


def _run_table(table, ns, out):
    """Evaluate {key: (source, getter)}; sources run once and are cached in ns. Each namespace gets `_src(name)` so a
    getter can read another source (028's Lambda_eta uses 025's lower mass)."""
    def src_ns(src):
        if src not in ns:
            ns[src] = EC06_SOURCES[src]()
            ns[src]["_src"] = src_ns
            out[f"_EC-06.{src}_checks_ok"] = bool(ns[src]["_checks_ok"])
        return ns[src]
    for key, (src, get) in table.items():
        out[key] = float(get(src_ns(src)))
    return out


def ec06():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ns, out = {}, {}
    _run_table(EC06, ns, out)
    _run_table(DRAFT_NUMBERS, ns, out)
    return out


def fit_values(rows, caps):
    """Claim numbers from the 014 fits. rows: fit_results.csv rows; caps: fit_candidates.ultralight_caps()."""
    def row(ds, pred):
        hit = [r for r in rows if r["dataset"] == ds and pred(r["model"])]
        assert len(hit) == 1, (ds, [r["model"] for r in hit])
        return hit[0]
    ec27 = lambda m: m.startswith("EC-03'") and "P18+BOSS" in m
    return {
        "EC-03'.beta_max_cap27_deg": float(row("act_wmap_planck", ec27)["beta_max_deg"]),
        "EC-03'.penalty_cap27_AWP": float(row("act_wmap_planck", ec27)["penalty_vs_ALP"]),
        "KK-02.penalty_k1_AWP": float(row("act_wmap_planck", lambda m: m == "KK-02 k=1")["penalty_vs_ALP"]),
        "EC-03'.cap_P18BOSS_pct": 100 * caps["ULF-R23-P18BOSS-m1e-32"],
        "EC-03'.cap_P18_pct": 100 * caps["ULF-R23-P18-m1e-32"],
        # draft only (task 033; skeptic 020 S4, 029 item 13): penalties with the ACT + Cosmoglobe combination
        "EC-03'.penalty_cap27_CG": float(row("act_cosmoglobe", ec27)["penalty_vs_ALP"]),
        "KK-02.penalty_k1_CG": float(row("act_cosmoglobe", lambda m: m == "KK-02 k=1")["penalty_vs_ALP"]),
    }


def regenerate(full_scan=True):
    vals, times = {}, {}
    for label, fn in (("007 (EC-01/02/04/05)", ec01_ec02_ec04_ec05), ("008 (EC-03, EC-03')", lambda: ec03(full_scan)),
                      ("009 (KK-02/03/04)", kk), ("012 (EC-03 photon, KK-05)", lambda: {**ec03_photon(), **kk05()}),
                      ("022-035 (EC-06, KK-04, skeptics)", ec06)):
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


def expected(path=EXPECTED):
    with open(path, newline="", encoding="utf-8") as f:
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
    if rule == "exp":  # same decimal exponent in scientific notation ("4e-70" printed as "10^-70")
        return math.floor(math.log10(abs(v))) == math.floor(math.log10(abs(e)))
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


def compare_draft(vals):
    """Draft-only numbers (src/draft_expected.csv). A value mismatch or a missing value is a failure; text absent from
    the draft is reported as "not in draft" and is not."""
    text = DRAFT.read_text(encoding="utf-8") if DRAFT.is_file() else ""
    table = []
    for x in expected(DRAFT_EXPECTED):
        v = vals.get(f"{x['claim']}.{x['key']}")
        if v is None:
            status = "NOT REGENERATED"
        else:
            status = "match" if _match(v, float(x["expected"]), x["rule"], int(x["digits"] or 0)) else "MISMATCH"
        in_draft = x["draft_text"] in text
        table.append(dict(claim=x["claim"], quantity=x["quantity"], value=v, printed=x["draft_text"],
                          rule=f"{x['rule']}{x['digits']}", status=status if in_draft else f"{status} (not in draft)",
                          ok=status == "match", in_draft=in_draft))
    return table


def print_draft_table(table):
    w = max(len(t["quantity"]) for t in table)
    print(f"{'claim':7s} | {'quantity':{w}s} | {'regenerated':>11s} | {'text':24s} | rule  | match")
    for t in table:
        print(f"{t['claim']:7s} | {t['quantity']:{w}s} | {fmt(t['value']):>11s} | {t['printed']:24s} | {t['rule']:5s} | "
              f"{t['status']}")
    bad = [t for t in table if not t["ok"]]
    print(f"draft numbers: {len(table) - len(bad)}/{len(table)} regenerated and matching; "
          f"{sum(not t['in_draft'] for t in table)} not (yet) in the draft (non-fatal)")
    return bad
