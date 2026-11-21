"""Skeptic 039, item A2-A4: 036's endpoint argument, measured ratio, and HB size shortfall.

A2  Measured ratio beta_ACT/beta_EK22 recomputed from data/.
A3  036 puts the Minami-Komatsu source endpoint at the diffuse-ISM mean density. But the calibrating photons are
    emitted inside dust grains, and the 0+ field at an emitting dipole includes the grain's own coherent field.
    For a uniform sphere (radius a, density rho_g) the Yukawa-smoothed density at its centre, for lam >> a, is
      rho_bar_grain = m^2 * 3M/(8 pi a) = rho_g a^2 /(2 lam^2)   (exact interior solution; checked numerically below).
    Compare with rho_*: below lam_cross = a sqrt(rho_g/(2 rho_*)) the MK endpoint is the grain, so
    beta_MK ~ (rho_grain - rho_*) > 0 and the predicted beta_ACT/beta_MK is POSITIVE (036 claims negative).
A4  Size at the HB photon bound, |beta| = sqrt(2 alpha) g rho_bar /(4 M_P m^2) (036 Eq. B10b, re-derived in README),
    for the ACT observer piece versus lam (ground observer with rock half-space, which 036 omitted), and for
    the grain-endpoint MK piece. Where does the size kill stop, and does it overlap the endpoint kill?
"""
import csv
import math
import os

import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
npass = nfail = 0


def check(name, ok):
    global npass, nfail
    npass += bool(ok)
    nfail += (not ok)
    print(("PASS" if ok else "FAIL") + ": " + name)


def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


D = lambda n: load(os.path.join(ROOT, "data", n))
cos, pin, meas, bnd = D("cosmo_inputs.csv"), D("pgt0plus_inputs.csv"), D("birefringence_measurements.csv"), D("bounds.csv")
loc = load(os.path.join(HERE, "inputs_039.csv"))
v = lambda d, k: float(d[k]["value"])

# ---- A2 measured ratio
b_act, s_act = v(meas, "DPK25-beta-baseline"), float(meas["DPK25-beta-baseline"]["err_plus"])
b_ek, s_ekp, s_ekm = v(meas, "EK22-beta"), float(meas["EK22-beta"]["err_plus"]), float(meas["EK22-beta"]["err_minus"])
b_mk20, s_mk20 = v(meas, "MK20-beta"), float(meas["MK20-beta"]["err_plus"])
r = b_act / b_ek
sr = r * math.hypot(s_act / b_act, 0.5 * (s_ekp + s_ekm) / b_ek)
print(f"A2 beta_ACT/beta_EK22 = {r:.3f} +- {sr:.3f} (linear propagation); ACT/MK20 = {b_act / b_mk20:.2f}")
# Monte Carlo of the ratio (Gaussian, independent) to check the linear error and P(ratio < 0)
rng = np.random.default_rng(39)
ra = rng.normal(b_act, s_act, 400000) / rng.normal(b_ek, 0.5 * (s_ekp + s_ekm), 400000)
print(f"   MC: median {np.median(ra):.2f}, 16-84% [{np.percentile(ra, 16):.2f}, {np.percentile(ra, 84):.2f}], "
      f"P(ratio < 0) = {np.mean(ra < 0):.4f}")
check("A2 measured ratio 0.63 +- 0.28 reproduced from data/", abs(r - 0.63) < 0.01 and abs(sr - 0.28) < 0.01)

# ---- units
hbarc = v(cos, "hbar_c_eV_m")
H0 = v(cos, "H0") / (v(pin, "Mpc_m") / 1e3) * v(cos, "hbar_eV_s")
MP, kg, mp = v(pin, "M_P_red_eV"), v(pin, "eV_per_kg"), v(pin, "m_p_eV")
KGM3 = kg * hbarc**3
rho_star = v(cos, "Omega_m") * 3 * H0**2 * MP**2 * (1 + v(cos, "z_star"))**3
rho_ism = v(pin, "n_ism_diffuse") * 1e6 * hbarc**3 * mp
rho_g = v(loc, "rho_grain") * KGM3
print(f"   rho_* = {rho_star / KGM3:.2e} kg/m^3, rho_ISM = {rho_ism / KGM3:.1e} kg/m^3")

# ---- A3 grain endpoint
# numerical check of the interior formula: rho_bar(0) = m^2 int_0^a rho_g e^{-m r} r dr  -> rho_g a^2/(2 lam^2) for lam >> a
for a_m in (1e-7,):
    for lam in (1e-5, 1e-3):
        num = (1 / lam**2) * quad(lambda rr: math.exp(-rr / lam) * rr, 0, a_m)[0]
        check(f"A3a grain-centre smoothed density at lam = {lam:g} m equals a^2/(2 lam^2) (ratio {num / (a_m**2 / (2 * lam**2)):.4f})",
              abs(num / (a_m**2 / (2 * lam**2)) - 1) < 2e-2)


def rho_grain(lam, a):
    """centre of a uniform sphere, valid for all lam (closed form of m^2 int_0^a e^{-mr} r dr) in eV^4."""
    x = a / lam
    f = x**2 / 2 - x**3 / 3 + x**4 / 8 if x < 1e-3 else 1 - (1 + x) * math.exp(-x)  # series avoids cancellation
    return rho_g * f


for key in ("a_grain_small", "a_grain_large"):
    a = v(loc, key)
    lc = a * math.sqrt(rho_g / (2 * rho_star))
    print(f"A3 a = {a * 1e6:.2f} um: grain endpoint beats rho_* for lam < {lc / 1e3:.1f} km "
          f"(rho_grain/rho_* at 100 um: {rho_grain(1e-4, a) / rho_star:.1e})")
lc_small = v(loc, "a_grain_small") * math.sqrt(rho_g / (2 * rho_star))
lc_large = v(loc, "a_grain_large") * math.sqrt(rho_g / (2 * rho_star))
check("A3b for lam below ~2-5 km the MK source endpoint is the grain, not the ISM mean: predicted ratio sign is + there",
      rho_grain(1e3, v(loc, "a_grain_small")) > rho_star)

# ---- A4 size at the HB bound
g_hb = float(bnd["A-HB-Ayala"]["value_eV"].replace("<", "").strip())
alpha = 1 / 3
rho_film, t_film = v(pin, "rho_film") * KGM3, v(pin, "t_film_assumed")
rho_air, rho_rock = v(pin, "rho_air_5200m") * KGM3, v(loc, "rho_rock") * KGM3
beta_need = math.radians(b_act)


def rho_ground(lam):
    """lower estimate at an analyzer on the ground: film for lam < 1 m, then max(air/2, rock half-space/2 e^{-h/lam})."""
    if lam < 1.0:
        return rho_film * min(1.0, t_film / (2 * lam))
    h = 10.0  # receiver ~10 m above rock (choice)
    rock = 0.5 * rho_rock * math.exp(-h / lam) if lam < 1e6 else 0.5 * rho_rock  # half-space, lam << R_E
    return max(0.5 * rho_air, rock)


def beta_size(rho_bar, lam):
    m = hbarc / lam
    return math.sqrt(2 * alpha) * g_hb * rho_bar / (4 * MP * m**2)


b100 = beta_size(rho_ground(1e-4), 1e-4)
print(f"A4 ACT observer piece at 100 um: {b100:.2e} rad, short by {math.log10(beta_need / b100):.1f} orders (036: 7e-27, 23.7)")
check("A4a reproduces 036's 7e-27 rad / 23.7 orders", abs(math.log10(b100 / 7e-27)) < 0.05)
bmk100 = beta_size(rho_grain(1e-4, v(loc, "a_grain_large")), 1e-4)
print(f"   MK piece at 100 um with the grain endpoint (a = 0.25 um): {bmk100:.1e} rad, short by "
      f"{math.log10(beta_need / bmk100):.1f} orders (036's ISM-mean 'universal piece': 41.4 orders)")
check("A4b the MK-measured piece is set by the grain, ~15 orders larger than 036's universal piece",
      math.log10(beta_need / bmk100) < 30)
print("   lam [m] | rho_ground [kg/m3] | beta_ACT size [rad] | short [orders]")
lam_reach = None
for lam in np.logspace(-4, 6.5, 22):
    b = beta_size(rho_ground(lam), lam)
    if lam_reach is None and b > beta_need:
        lam_reach = lam
    print(f"   {lam:9.2e} | {rho_ground(lam) / KGM3:9.2e} | {b:9.2e} | {math.log10(beta_need / b):6.1f}")
print(f"A4 size kill (HB g, ground observer, half-space) holds up to lam ~ {lam_reach:.1e} m")
check("A4c size kill and the (corrected) endpoint-sign kill overlap: size holds to >> 5 km",
      lam_reach is not None and lam_reach > 10 * lc_large)
print(f"== {npass} PASS, {nfail} FAIL ==")
