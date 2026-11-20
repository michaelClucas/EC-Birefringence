"""Checks of the coupling/mass algebra and reproduction of published numbers (PASS/FAIL)."""
import numpy as np
import sympy as sp

from common import inp, MP, GEV, H_ev, ALPHA
from variants import soz_window


def report(name, ok, detail):
    print(f"{'PASS' if ok else 'FAIL'}: {name} -- {detail}")


# 1. 1502.03694 Eq. (24): f = kappa^-1 sqrt(8/3), kappa^-1 = reduced Planck mass
f = np.sqrt(8 / 3) * MP
report("f_CF = sqrt(8/3) M_P matches data value (~4e18 GeV)", abs(f / inp("f_CF") - 1) < 1e-3,
       f"{f/GEV:.4e} GeV")

# 2. E/N for unit axial charges on all SM fermions = 2 Qbar^2/N_q = 8/3 (sympy, exact)
Qu, Qd, Qe = sp.Rational(2, 3), sp.Rational(-1, 3), -1
Nc, ng = 3, 3
E_ = ng * (Nc * Qu**2 + Nc * Qd**2 + Qe**2)          # sum_f N_c Q_f^2 (Dirac fermions)
N_ = ng * (sp.Rational(1, 2) + sp.Rational(1, 2))   # sum over quark flavours of T(R)=1/2
Qbar2 = E_                                          # 1502.03694: Qbar^2 = sum_f Q_f^2 (colour counted)
report("E/N = 8/3 for universal unit charges", sp.simplify(E_ / N_ - sp.Rational(8, 3)) == 0,
       f"E = {E_}, N = {N_}")
report("2 Qbar^2 / N_q = E/N (1502.03694 normalisation f_a = f/N_q)",
       sp.simplify(2 * Qbar2 / 6 - E_ / N_) == 0, "N_q = 6")
report("E/N stored in data file", abs(inp("EoverN_universal") - 8 / 3) < 1e-6, "")

# 3. local-amplitude beta for a QCD-mass axion is independent of f_a (sympy)
fa, a_, c, rho, K = sp.symbols("f_a alpha c rho K", positive=True)
g = a_ * c / (2 * sp.pi * fa)
m = K / fa
beta = g * sp.sqrt(2 * rho) / (2 * m)
report("beta_amp = g sqrt(2 rho)/(2m) independent of f_a when m f_a = const",
       sp.diff(sp.simplify(beta), fa) == 0, f"beta = {sp.simplify(beta)}")

# 4. QCD mass vs Lattanzi-Mercuri quoted value (LO formula vs NLO 5.70 ueV)
m_lm = inp("maxfa_QCD") / inp("f_LM")
report("m(f_LM) from 1511.02867 Eq. (22) vs 0911.2698 Eq. (5) quoted 3.04e-12 eV",
       abs(m_lm / inp("m_LM_quoted") - 1) < 0.1, f"{m_lm:.3e} eV")
m_cf = inp("maxfa_QCD") / inp("f_CF")
report("m(f_CF) is O(1e-12 eV) as quoted in 1502.03694 Eq. (27)",
       0.3 < m_cf / inp("m_CF_quoted") < 3, f"{m_cf:.3e} eV")

# 5. reproduce 2512.21888 Eqs. (4), (6), (12) at m = 1e-24 eV, kappa = f_phi = 1
phi0, g_req, g_wash, x = soz_window()
report("2512.21888 Eq. (4): phi_0 = 2.5e12 GeV", abs(phi0 / GEV / 2.5e12 - 1) < 0.03, f"{phi0/GEV:.3e} GeV")
report("2512.21888 Eq. (6): g_min = 2.16e-15 GeV^-1", abs(g_req * GEV / 2.16e-15 - 1) < 0.03,
       f"{g_req*GEV:.3e}")
report("2512.21888 Eq. (12): washout g_max = 9.6e-16 GeV^-1", abs(g_wash * GEV / 9.6e-16 - 1) < 0.03,
       f"{g_wash*GEV:.3e} (Omega_c h^2 = {inp('omch2')} vs their 0.11933)")

# 6. heavy-fermion decoupling of a derivative coupling: B1(tau) = 1 - tau arcsin^2(1/sqrt tau)
#    -> -1/(3 tau) = -m_a^2/(12 m_f^2) for tau = 4 m_f^2/m_a^2 -> infinity (1708.00443 loop function)
t = sp.symbols("t", positive=True)       # t = 1/tau
B1 = 1 - sp.asin(sp.sqrt(t)) ** 2 / t
ser = sp.series(B1, t, 0, 2).removeO()
report("B1 -> -1/(3 tau) for heavy fermions", sp.simplify(ser + t / 3) == 0, f"series {ser}")
supp = (H_ev(1 / (1 + inp("z_star"))) / inp("m_e")) ** 2 / 12
report("derivative-coupling photon term suppressed by (H_LSS/m_e)^2/12 < 1e-60", supp < 1e-60,
       f"{supp:.1e}")
