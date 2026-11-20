"""Skeptic 011: independent recomputation of the EC-02 shortfall (task 007 numbers).

beta = (xi kappa / 4) * int n5 dt, LSS -> today, n5 = net helicity density x speed.
Massless-neutrino envelope (speed = 1), which is an upper bound on |n5|.
Own LCDM integration; inputs from data/cosmo_inputs.csv only. Does not import task-007 code.
"""
import csv
import os

import numpy as np
from scipy.integrate import quad

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
with open(os.path.join(ROOT, "data", "cosmo_inputs.csv"), newline="", encoding="utf-8-sig") as fh:
    D = {r["id"]: r for r in csv.DictReader(fh)}
with open(os.path.join(ROOT, "data", "birefringence_measurements.csv"), newline="", encoding="utf-8-sig") as fh:
    M = {r[next(iter(r))]: r for r in csv.DictReader(fh)}


def v(k):
    return float(D[k]["value"])


HBAR = v("hbar_eV_s")
H0 = v("H0") / v("Mpc_km") * HBAR               # eV
MP = v("MP_reduced_eV")
KAPPA = 1 / MP**2
T0 = v("T_cmb_K") * v("k_B_eV_K")
h = v("H0") / 100
OM = v("Omega_m")
OG = (np.pi**2 / 15) * T0**4 / (3 * H0**2 * MP**2)
OR = OG * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * v("N_eff_SM"))
OL = 1 - OM - OR
ZS = v("z_star")
BETA = float(M["DPK25-beta-baseline"]["value"]) * np.pi / 180
XI95 = v("xi_nu") + 2 * float(D["xi_nu"]["err"])
Tnu0 = (4 / 11) ** (1 / 3) * T0


def H(a):
    return H0 * np.sqrt(OR / a**4 + OM / a**3 + OL)


def dn(xi):
    """net number density (particles - antiparticles) of one helicity-state species / T^3."""
    return (np.pi**2 * xi + xi**3) / (6 * np.pi**2)


# int_{LSS}^{today} a^-3 dt  (n5 = n5_0 a^-3 while relativistic-counting)
I3 = quad(lambda a: a**-3 / (a * H(a)), 1 / (1 + ZS), 1, epsrel=1e-10, limit=400)[0]

out = []
for lab, nflav, xi in (("xi_nu<=%.3f, 3 flavours" % XI95, 3, XI95),
                       ("xi_mu,tau<=0.57 (N_eff only), 2 flavours", 2, 0.57)):
    n5_0 = nflav * dn(xi) * Tnu0**3
    beta = KAPPA / 4 * n5_0 * I3
    out.append((lab, n5_0 * (1 + ZS) ** 3, beta, np.log10(BETA / beta), BETA / beta))
    print(f"{lab}: n5(LSS) = {n5_0*(1+ZS)**3:.3e} eV^3, beta(xi=1) = {beta*180/np.pi:.3e} deg, "
          f"shortfall 10^{np.log10(BETA/beta):.2f}, xi needed = {BETA/beta:.2e}, "
          f"Lambda_eff = M_P sqrt(2/xi) = {MP*np.sqrt(2*beta/BETA)/1e9:.2e} GeV")

# Reig dual-axion completion: g_ag g_ae = xi kappa m_e  -> xi_max
xi_max = v("helioscope_gag_gae_GeV") * 1e-9 / (KAPPA * v("m_e_eV"))
print(f"Reig/helioscope: xi < {xi_max:.2e}")


def report(name, ok, det=""):
    print(("PASS" if ok else "FAIL") + ": " + name + (" -- " + det if det else ""))


report("EC-02 007 value 7.5e-30 deg reproduced within 10%", abs(out[0][2] * 180 / np.pi / 7.5e-30 - 1) < 0.1,
       f"{out[0][2]*180/np.pi:.3e}")
report("EC-02 007 'xi needed 2.9e28' within 10%", abs(out[0][4] / 2.9e28 - 1) < 0.1, f"{out[0][4]:.3e}")
report("generous case 6.3e-29 deg within 10%", abs(out[1][2] * 180 / np.pi / 6.3e-29 - 1) < 0.1,
       f"{out[1][2]*180/np.pi:.3e}")
report("Reig xi bound 9.4e17 within 5%", abs(xi_max / 9.4e17 - 1) < 0.05, f"{xi_max:.3e}")
