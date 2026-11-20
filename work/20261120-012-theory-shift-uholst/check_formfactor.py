"""Question A(2): photon coupling of a pseudoscalar with only a derivative coupling (d phi/f) c_f J5_f.

Effective coupling (1708.00443 Eq. for C_gamma gamma^eff; Quevillon-Smith 1903.12559):
  g_eff = (alpha/(pi f)) sum_f N_c Q_f^2 c_f B1(tau_f) + (explicit phi F Ftilde term, absent here),
  B1(tau) = 1 - tau f(tau)^2, tau = 4 m_f^2/q^2, f(tau) = arcsin(1/sqrt(tau)) for tau >= 1.
The "1" is the anomaly; tau f^2 is the pseudoscalar (2 i m psibar gamma5 psi) triangle.
F1  tau f(tau)^2 = 2 int_0^1 dx int_0^{1-x} dy m^2/(m^2 - x y q^2)  (independent Feynman-parameter form)
F2  series: B1 = -q^2/(12 m^2) - q^4/(90 m^4) + ...  (sympy)
F3  heavy-fermion limit: B1 -> 0, light-fermion limit (q >> m): B1 -> 1 (anomaly restored)
F4  sum_f N_c Q_f^2 = 8 for unit charges on all SM fermions (E of 1502.03694 / 008)
F5  suppression g_eff/g_anomaly for cosmological q (m_phi or H at recombination); medium estimate
On-shell, q^2 -> m_phi^2 (box phi = m_phi^2 phi for a free field), so the relevant suppression is (m_phi/m_e)^2/12.
"""
import numpy as np
import sympy as sp
from scipy.integrate import dblquad

from common import ALPHA, M_E, H0, OM, MP, cin, check, info

# F1
def tauf2(r):
    """tau f(tau)^2 for r = q^2/m^2 < 4."""
    tau = 4 / r
    return tau * np.arcsin(1 / np.sqrt(tau)) ** 2


ok = True
for r in (0.01, 0.5, 2.0, 3.6):
    num = 2 * dblquad(lambda y, x: 1 / (1 - x * y * r), 0, 1, 0, lambda x: 1 - x, epsabs=1e-13, epsrel=1e-12)[0]
    ok &= abs(num / tauf2(r) - 1) < 1e-9
check("F1 pseudoscalar triangle: closed form = Feynman-parameter integral (q^2/m^2 = 0.01..3.6)", ok)

# F2
z = sp.symbols("z", positive=True)   # z = q^2/m^2
B1 = 1 - (4 / z) * sp.asin(sp.sqrt(z) / 2) ** 2
ser = sp.series(B1, z, 0, 3).removeO()
check("F2 B1 = -z/12 - z^2/90 + O(z^3), z = q^2/m^2", sp.simplify(ser + z / 12 + z ** 2 / 90) == 0, str(ser))

# F3: large-q limit (tau -> 0, analytic continuation f = -(1/4)[ln((1+s)/(1-s)) - i pi]^2, s = sqrt(1-tau))
def B1_large(tau):
    s = np.sqrt(1 - tau)
    fz = -0.25 * (np.log((1 + s) / (1 - s)) - 1j * np.pi) ** 2
    return 1 - tau * fz


check("F3 heavy fermion: B1 -> 0; light fermion (q >> m): |B1 - 1| -> 0",
      abs(1 - tauf2(1e-6)) < 1e-6 and abs(B1_large(1e-8) - 1) < 1e-5,
      f"B1(q/m=1e-3) = {1 - tauf2(1e-6):.2e}, B1(m/q=5e-5) = {B1_large(1e-8):.6f}")

# F4
QSUM = 3 * (1 * 1) + 3 * 3 * ((2 / 3) ** 2 + (1 / 3) ** 2)
check("F4 sum_f N_c Q_f^2 = 8 over three generations", abs(QSUM - 8) < 1e-12)

# F5: numbers.  Leptons dominate (lightest charged fermions); electron term ~ (q/m_e)^2/12
masses = {"e": (M_E, 1, 1)}   # only the electron is needed for the leading term; heavier ones are smaller


def ratio(q):
    return sum(Nc * Q ** 2 * (q / m) ** 2 / 12 for m, Nc, Q in masses.values()) / QSUM


z_star = cin("z_star")
H_star = H0 * np.sqrt(OM * (1 + z_star) ** 3 * (1 + (1 + z_star) / 3400))   # radiation share approx (order of magnitude)
for lab, q in (("m_phi = 6e-33 eV (008 plateau)", 6e-33), ("m_phi = 1e-28 eV (washout edge)", 1e-28),
               ("q ~ H(z*)", H_star)):
    info(f"g_eff/g_anomaly at {lab}: {ratio(q):.2e}")
g_T = 8 * ALPHA / (np.pi * np.sqrt(8 / 3) * MP) * 1e9    # GeV^-1, EC-03' torsion-scale anomaly value
info(f"EC-03' anomaly value g_T = {g_T:.3e} GeV^-1 (008: 4.67e-21); derivative coupling gives g_eff <= {g_T * ratio(1e-28):.1e} GeV^-1")
check("F5 g_eff/g_anomaly < 1e-60 for every m_phi in the beta window (<= 1e-27 eV)", ratio(1e-27) < 1e-60,
      f"{ratio(1e-27):.1e} at 1e-27 eV")

# Medium (plasma) estimate: dimensionless density n_e/m_e^3 at recombination (crude upper scale of any
# thermal correction to the pseudoscalar triangle; not a computed amplitude)
from common import binp
n_b0 = cin("ombh2") * cin("rho_crit_h2") * 1e9 / cin("m_H")        # cm^-3
n_e = (1 - cin("Y_P")) * n_b0 * (1 + z_star) ** 3                  # fully ionised upper value
hbarc_cm = binp("hbarc") * 100
info(f"n_e(z*) <= {n_e:.0f} cm^-3, n_e/m_e^3 = {n_e * hbarc_cm ** 3 / M_E ** 3:.1e} (crude upper scale of medium terms)")
check("F6 medium-term scale n_e/m_e^3 at z* below 1e-25", n_e * hbarc_cm ** 3 / M_E ** 3 < 1e-25)
