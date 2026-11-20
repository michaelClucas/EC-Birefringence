"""Independent check of g_eff/g_anom = (m_phi/m_e)^2/96 for a pure derivative coupling (d phi/f) J5.

Route (different from 012's closed form vs double Feynman integral):
  pseudoscalar triangle P(r) = 2 int_simplex dx dy 1/(1 - x y r - i0), r = q^2/m^2.
  Its absorptive part follows from delta(1 - x y r) on the simplex: Im P(s) = (2 pi/s) ln((1+v)/(1-v)), v = sqrt(1-4/s).
  P(r) is then rebuilt from an unsubtracted dispersion integral P(r) = (1/pi) int_4^inf Im P(s)/(s - r) ds.
  B1 = 1 - P (the "1" is the anomaly; P(0) = 1 is the heavy-fermion low-energy theorem, Quevillon-Smith 1903.12559:
  photon couplings come entirely from non-anomalous loops, so a heavy vector-like fermion decouples).
  dB1/dr at r = 0 from the dispersion integral: -(1/pi) int Im P(s)/s^2 ds -> expect -1/12.
  Normalising to g_anom (all SM fermions, sum N_c Q^2 = 8, universal axial coupling) gives 1/(12*8) = 1/96.
Also: sensitivity of 96 to light quarks (perturbative current masses vs hadronic, where the u,d,s pieces go
through pi0/eta mixing and scale as (m_phi/m_pi)^2).
"""
import numpy as np
from scipy.integrate import quad

from common15 import check, info, val, NPASS


def imP(s):
    v = np.sqrt(1 - 4 / s)
    return 2 * np.pi / s * np.log((1 + v) / (1 - v))


def P_disp(r, power=1):
    # substitution s = 4/(1-u^2)... use t = 4/s in (0,1] for a finite range
    f = lambda t: imP(4 / t) / (4 / t - r) ** power * 4 / t ** 2
    return quad(f, 0, 1, epsabs=1e-14, epsrel=1e-12, limit=400)[0] / np.pi


def P_closed(r):
    return 4 / r * np.arcsin(np.sqrt(r) / 2) ** 2


check("D1 dispersion P(0) = 1 (heavy-fermion pseudoscalar loop cancels the anomaly)", abs(P_disp(0.0) - 1) < 1e-9,
      f"{P_disp(0.0):.12f}")
ok = all(abs(P_disp(r) / P_closed(r) - 1) < 1e-8 for r in (0.01, 0.5, 2.0, 3.6))
check("D2 dispersion P(r) = (4/r) arcsin^2(sqrt(r)/2) for r = 0.01..3.6", ok)
slope = P_disp(0.0, power=2)      # dP/dr at 0 = (1/pi) int ImP/s^2
check("D3 dB1/dr|_0 = -1/12 from the dispersion integral", abs(slope - 1 / 12) < 1e-9, f"dP/dr = {slope:.12f}")
curv = P_disp(0.0, power=3)       # (1/2) d2P/dr2 at 0 = (1/pi) int ImP/s^3
check("D4 next term: B1 = -r/12 - r^2/90", abs(curv - 1 / 90) < 1e-9, f"{curv:.12f}")
C = 1 / (slope / 8)
check("D5 coefficient C in g_eff/g_anom = (m_phi/m_e)^2/C with the electron only", abs(C - 96) < 1e-6, f"C = {C:.6f}")

m_e = val("cosmo_torsion_axion_inputs.csv", "m_e")
for m in (6e-33, 1e-28):
    info(f"g_eff/g_anom at m_phi = {m:.0e} eV: {(m / m_e) ** 2 / 96:.2e}")

# Light quarks: quark masses are not in data/, so no number here. Treating u,d perturbatively with current masses
# m_q would give C = 96/(1 + (4/3)(m_e/m_u)^2 + (1/3)(m_e/m_d)^2), a few-% shift; see README.
info("hadronic (physical) treatment: u,d pieces enter via pi0/eta mixing, ~(m_phi/m_pi)^2, so C = 96 to < 1e-3")
print(f"formfactor_dispersion: {NPASS[0]} PASS, {NPASS[1]} FAIL")
