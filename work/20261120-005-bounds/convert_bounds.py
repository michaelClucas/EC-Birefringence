"""Convert and derive bounds for task 20261120-005-bounds.

Reads data/bounds.csv, data/bounds_inputs.csv and data/birefringence_measurements.csv.
No measured number is hardcoded here. Writes derived_bounds.csv next to this script
and prints a summary. Single-threaded, runs in well under a second.

Run from the repo root:
    .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-005-bounds/convert_bounds.py
"""
import csv
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "derived_bounds.csv"

NUM = r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?"


def load(name, key="id"):
    with open(ROOT / "data" / name, encoding="utf-8") as f:
        return {r[key]: r for r in csv.DictReader(f)}


B = load("bounds.csv")
I = load("bounds_inputs.csv")
M = load("birefringence_measurements.csv")
inp = {k: float(v["value"]) for k, v in I.items()}


def nums(s):
    return [float(x) for x in re.findall(NUM, s, flags=re.I)]


# ---------------------------------------------------------------- unit check
# Recompute value_eV from value_published for rows with a plain unit factor
# and compare with what was typed into the CSV.
factor = {"GeV": inp["GeV"], "GeV^-1": 1 / inp["GeV"], "eV": 1.0, "-": 1.0,
          "m^-1": inp["hbarc"]}
checked, mismatches = 0, []
for k, r in B.items():
    u = r["units_published"]
    if u not in factor:
        continue
    pub, conv = nums(r["value_published"]), nums(r["value_eV"])
    if not pub or len(pub) != len(conv) or "(" in r["value_published"]:
        continue
    for a, b in zip(pub, conv):
        if not math.isclose(a * factor[u], b, rel_tol=0.02):
            mismatches.append((k, a * factor[u], b))
    checked += 1
print(f"unit check: {checked} rows recomputed, mismatches: {mismatches}")

# ------------------------------------------------ fermion b -> axial torsion S
# Minimal coupling (KRT arXiv:0712.4393 Eq. (2)-(3), xi_4^(4) = 3/4): b_mu = -(3/4) S_mu
# with KRT's A_mu = (1/6) eps_{mu a b c} T^{abc}, same contraction as our S^mu.
MIN = 3.0 / 4.0


def two_sigma(row_id):
    c, s = nums(B[row_id]["value_eV"])[:2]
    return abs(c) + 2 * s


derived = []


def add(id_, quantity, value, note):
    derived.append({"id": id_, "quantity": quantity, "value_eV": f"{value:.3g}", "note": note})


b_perp = {
    "electron": (max(two_sigma("F-e-bX"), two_sigma("F-e-bY")), "|c|+2sigma of Heckel 2008 b~_X,Y"),
    "neutron": (nums(B["F-n-bperp-Allm"]["value_eV"])[0], "Allmendinger 2014 |b~_perp| (68% CL)"),
    "proton": (nums(B["F-p-bperp"]["value_eV"])[0], "Brown 2010 via Data Tables D10"),
}
b_T = {
    "electron": (2 * two_sigma("F-e-bT"), "2 x (|c|+2sigma) of Heckel 2008 half-combination, torsion-only"),
    "neutron": (nums(B["F-n-bT"]["value_eV"])[0], "Data Tables S2/D15 order of magnitude"),
    "proton": (nums(B["F-p-bT"]["value_eV"])[0], "Data Tables S2/D10 order of magnitude"),
}

# Sun velocity relative to the CMB, projected on the equatorial (X-Y) plane of the
# Sun-centred frame. A pure CMB-frame S_0 appears in the Sun frame as
# S_J ~ -gamma v_J S_0, so |S_perp| = v_perp |S_0| to first order in v.
l, b = math.radians(inp["l_sun_cmb"]), math.radians(inp["b_sun_cmb"])
ra_p, dec_p, l_ncp = (math.radians(inp[k]) for k in ("gal_ngp_ra", "gal_ngp_dec", "gal_l_ncp"))
sin_dec = math.sin(dec_p) * math.sin(b) + math.cos(dec_p) * math.cos(b) * math.cos(l_ncp - l)
dec = math.asin(sin_dec)
v = inp["v_sun_cmb"] / inp["c_kms"]
v_perp = v * math.cos(dec)
add("v-perp", "Sun-CMB velocity component in equatorial plane (units of c)", v_perp,
    f"v = {v:.4e}, dipole Dec = {math.degrees(dec):.2f} deg")

for sp in ("electron", "neutron", "proton"):
    bp, note = b_perp[sp]
    bt, note_t = b_T[sp]
    S_perp = bp / MIN
    add(f"S-perp-{sp}", f"|S_X|,|S_Y| bound from {sp} (Sun frame), minimal coupling", S_perp, note)
    add(f"S0-boost-{sp}", f"|S_0| (CMB frame) from {sp} b_perp via Sun-CMB boost", S_perp / v_perp,
        "assumes S_mu purely timelike in CMB frame (FRW-compatible) and constant over the solar system")
    add(f"S0-direct-{sp}", f"|S_0| from {sp} b_T directly", bt / MIN, note_t)

# ------------------------------------------- beta -> constant CFJ k_AF and p_0
# Lookback time to recombination, flat LCDM, radiation neglected; integrate to z = inf.
H0_eV = inp["H0"] / (inp["Mpc"] / 1e3) * inp["hbar"]  # km/s/Mpc -> s^-1 -> eV
Om = inp["Omega_m"]
OL = 1 - Om
# closed form: t(inf) = 2/(3 H0 sqrt(OL)) asinh(sqrt(OL/Om))
t_L = 2 / (3 * H0_eV * math.sqrt(OL)) * math.asinh(math.sqrt(OL / Om))
t_L_Gyr = t_L * inp["hbar"] / inp["julian_year"] / 1e9
add("t-L", "lookback time to z=inf (eV^-1)", t_L, f"= {t_L_Gyr:.2f} Gyr; H0 = {H0_eV:.3e} eV")

beta = M["DPK25-beta-baseline"]
beta_rad = math.radians(float(beta["value"]))
sig_rad = math.radians(float(beta["err_plus"]))
# KM 0905.0031 Eq. (158), d=3, isotropic: delta psi = Y_00 k_(V)00 t_L = -(k_AF)_T t_L
kAF = beta_rad / t_L
add("kAF-ACT", "|(k_AF)_T| equivalent to ACT DR6 beta (constant CFJ, KM Eq. 158)", kAF,
    f"+- {sig_rad / t_L:.2e} eV; beta = {beta['value']} deg")
add("kV00-ACT", "|k^(3)_(V)00| = sqrt(4pi)|k_AF| for ACT beta", math.sqrt(4 * math.pi) * kAF, "")
add("p0-ACT", "|p_0| in L = -(1/2) p_mu A_nu Ftilde^{mu nu} for ACT beta (p = -2 k_AF)", 2 * kAF,
    "beta = (1/2) Int p_0 dt; matches KM: beta = |k_AF| t_L")

# Data Tables pairing beta(MK) = 0.35 deg <-> k_(V)00 = 7.32e-45 GeV (refs [142],[143])
kV_MK = nums(B["P-kV00-MK"]["value_eV"])[0]
kV_MK_KM = math.sqrt(4 * math.pi) * math.radians(0.35) / t_L
add("kV00-MK-check", "k_(V)00 for beta=0.35 deg from KM Eq. 158 (compare Data Tables row P-kV00-MK)",
    kV_MK_KM, f"Data Tables value {kV_MK:.3g} eV; ratio {kV_MK_KM / kV_MK:.2f} (unresolved)")

# ---------------------------------------------------------------- axion g
for k in ("A-CAST", "A-HB-Ayala", "A-HB-Dolan", "A-SN1987A", "A-Chandra", "A-MWD"):
    add(f"g-{k}", f"g_agg bound ({k}) in eV^-1", nums(B[k]["value_eV"])[0], B[k]["experiment"])

# ------------------------------------------------------------- KK-specific
yr_eV = inp["hbar"] / inp["julian_year"]
ad, ad_s = nums(B["K-alphadot-clock"]["value_published"])[:2]
add("alphadot-eV", "alpha_dot/alpha today in eV", ad * yr_eV, f"+- {ad_s * yr_eV:.2e} eV")
add("alphadot-H0", "alpha_dot/alpha today in units of H0", ad * yr_eV / H0_eV,
    f"+- {ad_s * yr_eV / H0_eV:.2e}")
lam = nums(B["K-ISL-EotWash"]["value_published"])[0] * 1e-6
add("m-yukawa", "mass of gravitational-strength scalar must exceed (eV)", inp["hbarc"] / lam, "hbar c/lambda")
invR = nums(B["K-mUED-LHC"]["value_published"])[0] * 1e12
add("R-mUED", "R upper bound if SM in bulk (eV^-1)", 1 / invR, f"= {inp['hbarc'] / invR:.2e} m")

# ------------------------------------------------- BBN trace torsion (Kranas)
lo, hi = nums(B["C-BBN-Kranas"]["value_published"])[:2]
T = inp["T_BBN"]
H_bbn = math.sqrt(math.pi ** 2 * inp["g_star_BBN"] / 90) * T ** 2 / inp["MP_reduced"]
add("H-BBN", "H at T = 1 MeV (eV)", H_bbn, "H = sqrt(pi^2 g*/90) T^2 / M_P")
add("T0-BBN", "|T_0| (trace torsion) bound at T ~ 1 MeV (eV)", 6 * max(abs(lo), abs(hi)) * H_bbn,
    "T_mu = -6 phi u_mu; trace part only, sign depends on their signature")

with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["id", "quantity", "value_eV", "note"])
    w.writeheader()
    w.writerows(derived)

for r in derived:
    print(f"{r['id']:<22} {r['value_eV']:>10}  {r['quantity']}  [{r['note']}]")
