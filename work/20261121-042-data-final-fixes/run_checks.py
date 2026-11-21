"""Task 20261121-042: numbers for skeptic 040's data items. Read-only; run from the repo root, pinned, 1 thread:
    .venv/Scripts/python tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261121-042-data-final-fixes/run_checks.py

C1  EC-04 upper coupling edge from Chandra NGC 1275 (data/bounds.csv A-Chandra, m < 1e-12 eV, which covers the
    m ~ H0 mode) instead of CAST: g_phigamma = xi/sqrt(A), A = 3/(2 kappa) (007 ec04_estimates), i.e. g = 0.82 xi/M_P.
C2  KK-05 charge-to-mass: m >= e M_P/sqrt(2) vs m_e, in orders (020 S2 arithmetic, independent recomputation).
C3  020 S3 hardcodes xi_needed(EC-02) = 2.9e28; compare with 007's regenerated value; Lambda scales as xi^-1/2.
C4  Strong-coupling comparisons against the highest observed SM diphoton mass (ATLAS 2.36 TeV,
    data/collider_photon_bounds.csv) instead of the 10 TeV convention: EC-06a Lambda_SC (031).
C5  2608.06480 Gaussian penalties computed by hand vs the src.reproduce fit rows (artifacts/reproduce/fit_results.csv,
    only if present, i.e. after a reproduce run; otherwise INFO).
"""
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src import claims_numbers as CN  # noqa: E402
from src import likelihood as L  # noqa: E402

NP = NF = 0


def check(name, ok, msg=""):
    global NP, NF
    NP += bool(ok)
    NF += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {msg}")


def rows(name):
    with open(ROOT / "data" / name, newline="", encoding="utf-8-sig") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


cos = {k: float(r["value"]) for k, r in rows("cosmo_inputs.csv").items()}
bnd = rows("bounds.csv")
col = rows("collider_photon_bounds.csv")
MP, alpha, me = cos["MP_reduced_eV"], cos["alpha_em"], cos["m_e_eV"]
gval = lambda rid: float(bnd[rid]["value_eV"].replace("<", "").strip())  # eV^-1

# C1 -------------------------------------------------------------------------------------------------------------
sqrtA = math.sqrt(1.5) * MP  # A = 3/(2 kappa) = 1.5 M_P^2
print(f"C1 g_phigamma = xi/sqrt(A) = {MP / sqrtA:.4f} xi/M_P")
g_ch, g_cast = gval("A-Chandra"), gval("A-CAST")
xi_ch, xi_cast = g_ch * sqrtA, g_cast * sqrtA
print(f"C1 Chandra g < {g_ch * 1e9:.1e} GeV^-1 ({bnd['A-Chandra']['experiment']}) -> xi < {xi_ch:.3e}")
print(f"C1 CAST    g < {g_cast * 1e9:.1e} GeV^-1 -> xi < {xi_cast:.3e}")
v007 = CN.ec01_ec02_ec04_ec05()
# 007 prints xi to 2 significant figures (1.7e+08), so agreement is to that rounding
check("C1a CAST edge equals 007's regenerated EC-04 xi_max_CAST (2 s.f.)", f"{xi_cast:.1e}" == f"{v007['EC-04.xi_max_CAST']:.1e}",
      f"({v007['EC-04.xi_max_CAST']:.3e})")
check("C1b Chandra edge = CAST edge x g_Chandra/g_CAST", abs(xi_ch / (xi_cast * g_ch / g_cast) - 1) < 1e-12)
check("C1c EC-04 window stays open: lab lower end < Chandra upper end", v007["EC-04.xi_min_lab"] < xi_ch,
      f"({v007['EC-04.xi_min_lab']:.0e} < xi < {xi_ch:.1e})")
print(f"C1 RESULT: EC-04 10^{math.log10(v007['EC-04.xi_min_lab']):.0f} <~ xi < {xi_ch:.2e} (Chandra; CAST gave "
      f"{xi_cast:.1e}, {math.log10(xi_cast / xi_ch):.2f} orders weaker)")

# C2 -------------------------------------------------------------------------------------------------------------
e = math.sqrt(4 * math.pi * alpha)
m_min = e * MP / math.sqrt(2)
orders = math.log10(m_min / me)
print(f"C2 KK-05: m >= e M_P/sqrt2 = {m_min:.3e} eV = {m_min / me:.3e} m_e ({orders:.2f} orders)")
check("C2 KK-05 mass gap ~21 orders", round(orders) == 21)

# C3 -------------------------------------------------------------------------------------------------------------
xi020, xi007 = 2.9e28, v007["EC-02.xi_needed"]
shift = 0.5 * math.log10(xi007 / xi020)
print(f"C3 020 hardcoded xi = {xi020:.2e}; 007 regenerated xi_needed = {xi007:.3e}; Lambda shifts by {shift:+.3f} dex")
check("C3 020's hardcoded xi agrees with 007 to 5% (Lambda to 2.5%)", abs(xi007 / xi020 - 1) < 0.05)

# C4 -------------------------------------------------------------------------------------------------------------
E_atlas = float(col["ATLAS21-mgg-max"]["value"]) * 1e12  # TeV -> eV
g31 = CN._mz_031()
lam = g31["lam_sc"](g31["xi_gw"], g31["ecc"]["m_gamma_PDG"], g31["one_m_cs"])
o10, oat = math.log10(g31["E_safe"] / lam), math.log10(E_atlas / lam)
print(f"C4 EC-06a Lambda_SC = {lam / 1e3:.2f} keV: {o10:.2f} orders below 10 TeV (convention), {oat:.2f} below the "
      f"ATLAS {E_atlas / 1e12:.2f} TeV diphoton maximum; BaBar sqrt(s) = {col['BABAR17-sqrts-4S']['value']} GeV")
check("C4 031 namespace ran clean", bool(g31["_checks_ok"]))
check("C4 the 10 TeV -> 2.36 TeV change moves shortfalls by < 1 order", 0 < o10 - oat < 1, f"({o10 - oat:.2f})")

# C5 -------------------------------------------------------------------------------------------------------------
fr = ROOT / "artifacts" / "reproduce" / "fit_results.csv"
if fr.is_file():
    with open(fr, newline="", encoding="utf-8") as f:
        fit = list(csv.DictReader(f))
    for ds in ("act_planck_esk26", "act_planck_esk26_dust"):
        sub = {r["model"]: r for r in fit if r["dataset"] == ds}
        if not sub:
            print(f"[INFO] C5 {ds}: no rows in fit_results.csv (reproduce not yet run with task 042's src)")
            continue
        lk = L.IsotropicBeta(ds)
        for m, r in sub.items():
            B = float(r["beta_max_deg"])
            pen = ((lk.mu - B) / lk.sigma) ** 2 if B < lk.mu else 0.0
            d0 = (lk.mu / lk.sigma) ** 2 - pen
            ok = abs(pen - float(r["penalty_vs_ALP"])) < 1e-3 and abs(d0 - float(r["dchi2_vs_null"])) < 1e-3
            check(f"C5 {ds} {m}: hand penalty {pen:.3f}, dchi2 vs null {d0:.3f}", ok)
else:
    print("[INFO] C5 skipped: artifacts/reproduce/fit_results.csv absent (run src.reproduce first)")
try:
    L.joint("act_planck_esk26", "act_wmap_planck")
    check("C5 joint() refuses ESK26 x DPK25 (shared ACT DR6 + Planck PR4)", False)
except Exception as ex:  # noqa: BLE001
    check("C5 joint() refuses ESK26 x DPK25 (shared ACT DR6 + Planck PR4)", True, f"({type(ex).__name__})")

print(f"== {NP} PASS, {NF} FAIL ==")
sys.exit(1 if NF else 0)
