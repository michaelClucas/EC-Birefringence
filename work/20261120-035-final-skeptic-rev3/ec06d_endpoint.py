"""Skeptic 035: observer-endpoint dependence of the EC-06d (SNY 0+) GR-compatible branch.

For closed p, CFJ gives beta = (3 xi/4) [ln Y(obs) - ln Y(src)], Y = R + Yc, Yc = 6 mu/b (025 derivation Sec. 4,
check_sny_limit.py Q1-Q2). 025 used the cosmological observer (R0 = mean rho today). In the heavy-mode
(GR-compatible) regime R tracks the LOCAL density (R = 8 pi G rho for dust), so the observer value is set by the
environment over the 0+ Compton length. We take Yc = R(z*)/s*, s* = 4 delta/(9 mu) (025/028), mu = 1, delta = 1%.
Densities are in units of rho_crit; only ratios enter.
"""
import csv, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.data import measurements  # noqa: E402


def rows(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def find(name, key):
    for f in (ROOT / "data").glob("*.csv"):
        for r in rows(f):
            if r.get("id") == key:
                return float(r["value"]), f.name
    raise KeyError(key)


Om, _ = find("", "Omega_m")
zs, _ = find("", "z_star")
rc_h2_GeVcm3, _ = find("", "rho_crit_h2")
h = 0.674  # only converts kg/m^3 to rho_crit; enters inside a log (0.3 per factor 2)
loc = {r["id"]: r for r in rows(Path(__file__).with_name("local_inputs.csv"))}
rho_crit_kgm3 = 1.87834e-26 * h**2  # PDG value in data/ note (= rho_crit_h2 converted)
m_p = 1.6726e-27  # kg
dens = {  # in units of rho_crit
    "cosmological observer (mean rho_m today)": Om,
    "local ISM (1 p/cm^3)": float(loc["n_ism"]["value"]) * 1e6 * m_p / rho_crit_kgm3,
    "L2 / solar wind (Planck, if Compton length >> spacecraft)": float(loc["n_sw_1AU"]["value"]) * 1e6 * m_p / rho_crit_kgm3,
    "ground, air at 5.2 km (ACT, Compton length >> instrument)": float(loc["rho_air_5200m"]["value"]) / rho_crit_kgm3,
    "inside a silicon detector (Compton length << mm)": float(loc["rho_si"]["value"]) / rho_crit_kgm3,
}
beta = measurements()["DPK25-beta-baseline"]["value"] * math.pi / 180
mu, delta = 1.0, 0.01
s_star = 4 * delta / (9 * mu)
R_star = Om * (1 + zs) ** 3
Yc = R_star / s_star
Y_star = Yc + R_star
print(f"s* = {s_star:.3e}, R(z*)/rho_crit = {R_star:.3e}, Yc/rho_crit = {Yc:.3e}")
out = []
for k, rho in dens.items():
    bx = 0.75 * math.log((Yc + rho) / Y_star)  # beta/xi in rad
    xi = beta / abs(bx)
    out.append((k, bx, xi))
    print(f"{k:60s} beta/xi = {bx:+.3e} rad  |xi| needed = {xi:.2e}  sign {'+' if bx > 0 else '-'}")
r = out[3][1] / out[2][1]
print(f"ratio ground/L2 = {r:+.2e}")
ok = abs(out[0][2] - 3 * mu * beta / delta) / (3 * mu * beta / delta) < 0.05
print(("PASS" if ok else "FAIL") + ": cosmological-observer xi reproduces 025's xi = 3 mu beta/delta")
print(("PASS" if abs(r) > 100 else "FAIL") + ": ground vs space endpoints differ by > 100x (environment dependence)")
