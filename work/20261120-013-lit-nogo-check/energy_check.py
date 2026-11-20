"""Independent check of the EC-03' energy requirement against data/ultralight_fraction_bounds.csv.

Task 20261120-013. Spectator field, quadratic potential, flat LCDM background,
field frozen at phi_i until recombination. For each mass, find phi_i giving the
ACT beta with the torsion-fixed coupling g_T = 8 alpha/(pi f_CF), then report
rho_phi today / rho_DM. Independent of task 008 code (written from scratch).

Run: .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-013-lit-nogo-check/energy_check.py
"""
import csv
import os

import numpy as np
from scipy.integrate import solve_ivp

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def load(name):
    with open(os.path.join(ROOT, "data", name), newline="", encoding="utf-8-sig") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


CI = load("cosmo_torsion_axion_inputs.csv")
BI = load("bounds_inputs.csv")
ME = load("birefringence_measurements.csv")
v = lambda d, k: float(d[k]["value"])

h = v(CI, "H0") / 100
H0 = v(CI, "H0") / (v(BI, "Mpc") / 1e3) * v(BI, "hbar")     # eV
MP = v(BI, "MP_reduced")                                     # eV
Oc = v(CI, "omch2") / h**2
Om = (v(CI, "omch2") + v(CI, "ombh2")) / h**2
TK = v(CI, "T_CMB") * v(CI, "k_B")                           # eV
rhoc = 3 * H0**2 * MP**2
Og = (np.pi**2 / 15) * TK**4 / rhoc
Or = Og * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * v(CI, "N_eff"))
OL = 1 - Om - Or
zs = v(CI, "z_star")
alpha = v(CI, "alpha_em")
gT = 8 * alpha / (np.pi * v(CI, "f_CF"))                     # eV^-1
beta = v(ME, "DPK25-beta-baseline") * np.pi / 180


def E(N):
    a = np.exp(N)
    return np.sqrt(Or / a**4 + Om / a**3 + OL)


def dlnE(N):
    a = np.exp(N)
    return (-4 * Or / a**4 - 3 * Om / a**3) / (2 * E(N) ** 2)


def need(m):
    mu = m / H0
    rhs = lambda N, y: [y[1], -(3 + dlnE(N)) * y[1] - (mu / E(N)) ** 2 * y[0]]
    s = solve_ivp(rhs, [-np.log(1 + zs), 0.0], [1.0, 0.0], rtol=1e-10, atol=1e-13)
    p0, dp0 = s.y[0, -1], s.y[1, -1]            # phi/phi_i and dphi/dN today
    dphi = abs(p0 - 1.0)
    phii = 2 * beta / (gT * dphi)
    rho0 = 0.5 * (H0 * dp0 * phii) ** 2 + 0.5 * (m * p0 * phii) ** 2
    return phii, rho0 / (Oc * rhoc)


print(f"g_T = {gT*1e9:.3e} GeV^-1, beta = {np.degrees(beta):.3f} deg, H0 = {H0:.4e} eV")
print("m [eV]     phi_i/M_P   rho_phi,0/rho_DM  (needed for beta = ACT central)")
rows = []
for m in [2e-33, 3e-33, 4e-33, 5e-33, 6e-33, 7e-33, 8e-33, 1e-32, 2e-32]:
    phii, fr = need(m)
    rows.append((m, fr))
    print(f"{m:.1e}   {phii/MP:8.3f}   {fr:10.3f}")

# PASS checks
ok = []
fmin = min(r[1] for r in rows)
mmin = [r[0] for r in rows if r[1] == fmin][0]
ok.append(("minimum required fraction within 0.09-0.14 (task 008: ~12% at ~6e-33 eV)", 0.09 <= fmin <= 0.14))
ok.append(("minimum lies in 4e-33..8e-33 eV", 4e-33 <= mmin <= 8e-33))
UB = load("ultralight_fraction_bounds.csv")
r23 = float(UB["ULF-R23-P18-m1e-32"]["value"].strip("< ")) / 0.12
f32 = [r[1] for r in rows if r[0] == 1e-32][0]
ok.append((f"at 1e-32 eV required {f32:.3f} > Planck18 cap {r23:.3f} (excluded there)", f32 > r23))
for name, flag in ok:
    print("PASS" if flag else "FAIL", name)
