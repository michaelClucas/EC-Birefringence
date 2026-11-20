"""Skeptic 028, item 1 and item 4 (Lorentz-mode exchange).

E1  Compare 025's zeta solution with Golovnev-Koivisto 1808.05565 Eq. (6) (conformal H = a'/a):
      Delta zeta = -3 (psi' + H phi - (H' - H^2)/H psi)
E2  Rebuild the f(T) scalar boost sector (tetrad = boost(eta) x diag(a(1+Phi), a(1-Psi)...)) independently of
    025's action assembly and extract the eta-eta part of the second-order Lagrangian e f(T).
    This is the 'stiffness' of the Lorentz mode eta that the CFJ term pushes on.
E3  Add the CFJ term -1/2 xi T_mu K^mu and get the photon source in the eta equation.
Uses only the generic torsion helpers of task 025's common.py (Weitzenboeck torsion series, torsion scalar),
which were validated in tasks 022/025; the action, the eta-eta extraction and the comparisons are new here.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "20261120-025-theory-ec06c-costs"))
from common import EPS, series_inv, torsion_scalar, torsion_series, trunc  # noqa: E402

npass = nfail = 0


def check(name, ok, msg=""):
    global npass, nfail
    npass += bool(ok)
    nfail += (not ok)
    print(f"{'PASS' if ok else 'FAIL'}: {name} {msg}".rstrip())


tau, x, y, z = sp.symbols("tau x y z", real=True)
X = [tau, x, y, z]
a = sp.Function("a", positive=True)(tau)
ap, app = sp.diff(a, tau), sp.diff(a, tau, 2)
Hc = ap / a
Phi, Psi, eta = (sp.Function(n)(tau, x) for n in ("Phi", "Psi", "eta"))
zeta = sp.Function("zeta")(tau, x)

# E1: 025 form vs GK Eq. 6
f025 = -3 * sp.diff(Psi, tau) + 3 * app / ap * Psi - 3 * Hc * (Phi + 2 * Psi)
Hp = sp.diff(Hc, tau)
fGK = -3 * (sp.diff(Psi, tau) + Hc * Phi - (Hp - Hc**2) / Hc * Psi)
check("E1 025's zeta solution equals Golovnev-Koivisto 1808.05565 Eq. (6)", sp.simplify(f025 - fGK) == 0)

# E2: boosted tetrad, second order
eta_m = sp.diag(-1, 1, 1, 1)
base = sp.diag(a * (1 + EPS * Phi), a * (1 - EPS * Psi), a * (1 - EPS * Psi), a * (1 - EPS * Psi))
Lb = sp.eye(4)
Lb[0, 0] = Lb[1, 1] = 1 + EPS**2 * eta**2 / 2
Lb[0, 1] = Lb[1, 0] = EPS * eta
e = (Lb * base).applyfunc(lambda v: trunc(v))
g = (e.T * eta_m * e).applyfunc(lambda v: trunc(v))
gi = series_inv(g)
Tsc, v = torsion_scalar(torsion_series(e, X), g, gi)
T1 = sp.expand(Tsc).coeff(EPS, 1)
T2 = sp.expand(Tsc).coeff(EPS, 2)
T1_eta = sp.simplify(T1.subs({Phi: 0, Psi: 0}).doit())
print("   delta_1 T_sc (Phi = Psi = 0) =", T1_eta)
check("E2a first-order torsion scalar contains the boost: delta_1 T_sc = -4 a'/a^3 eta_x (up to sign)",
      sp.simplify(T1_eta**2 - 16 * ap**2 / a**6 * sp.diff(eta, x)**2) == 0)
T2_eta = sp.simplify(T2.subs({Phi: 0, Psi: 0}).doit())
print("   delta_2 T_sc (Phi = Psi = 0) =", T2_eta)
# e f(T) to second order: a^4 [f_T T2 + f_TT T1^2/2] + e1 f_T T1 + e2 f ; e has no eta dependence (boost det = 1)
F1, F2 = sp.Function("F1")(tau), sp.Function("F2")(tau)
L_ee = sp.expand(a**4 * (F1 * T2_eta + F2 * T1_eta**2 / 2))
# TEGR part a^4 F1 T2_eta must be a total derivative when F1 is constant
check("E2b delta_2 T_sc has no eta-eta term, so TEGR has no eta stiffness (frame is pure gauge in TEGR)", T2_eta == 0)
# eta-eta part with time-dependent F1: integrate by parts -> only F1' = F2 Tbar' pieces + F2 T1^2/2
Tbar = 6 * ap**2 / a**4
EL_ee = sp.euler_equations(L_ee, [eta], [tau, x])[0].lhs
EL_ee = sp.simplify(EL_ee.subs(sp.Derivative(F1, tau), F2 * sp.diff(Tbar, tau)))
print("   eta-eta part of the eta Euler-Lagrange expression:", sp.factor(EL_ee))
c_xx = sp.simplify(EL_ee.coeff(sp.diff(eta, x, 2)))
c_tt = sp.simplify(EL_ee.coeff(sp.diff(eta, tau, 2)))
print("   coefficient of eta_xx:", sp.factor(c_xx), "; of eta_tautau:", c_tt)
# physical units, a = 1 today: c_xx -> 16 f_TT H^2 (expected from the F2 T1^2/2 term)
check("E2c eta-eta Lagrangian is 8 f_TT (a'/a)^2 eta_x^2 in e f(T): coefficient of eta_xx = -16 f_TT a'^2/a^2, "
      "no eta_tautau (eta is a constraint field)",
      sp.simplify(c_xx + 16 * F2 * ap**2 / a**2) == 0 and c_tt == 0)

# E3: CFJ source. L_CFJ = -(xi/2) e T_mu K^mu with K^mu = A_nu Ft^{mu nu}, here as densities eK^mu
xi = sp.Symbol("xi")
Kt, Kx = sp.Function("Kt")(tau, x), sp.Function("Kx")(tau, x)
v1 = [sp.expand(c).coeff(EPS, 1) for c in v]
Lcfj = -xi / 2 * (v1[0] * Kt + v1[1] * Kx)
EL_cfj = sp.simplify(sp.euler_equations(Lcfj, [eta], [tau, x])[0].lhs)
print("   CFJ source in the eta equation:", EL_cfj)
check("E3 CFJ source = -(xi/2)(d_x K^tau + d_tau K^x)",
      sp.simplify(EL_cfj + xi / 2 * (sp.diff(Kt, x) + sp.diff(Kx, tau))) == 0)
print(f"== {npass} PASS, {nfail} FAIL ==")
