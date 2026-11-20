"""Question A(1): does a torsion pseudoscalar coupled as (d phi/f) c_q J5_q get a QCD potential?

S1  theta-bar bookkeeping under a field-dependent chiral rotation (sympy), four cases:
    (a) derivative coupling only (KR dual axion, BI x Nieh-Yan, PGT 0-)
    (b) "reading A": anomaly shift of theta_GG kept, quark-mass phases dropped
    (c) derivative coupling + Green-Schwarz term (explicit lambda phi G Gtilde)
    (d) KSVZ control: phi in the heavy-quark Yukawa phase
S2  KR dualization with and without a GS-modified Bianchi identity (sympy, algebraic).
S3  LO two-flavour chiral Lagrangian vacuum energy V(phi) for cases (a)-(c) (numerical
    minimisation over the meson phase), curvature compared with the textbook axion formula.
Toy parameters only (K, m_u, m_d are generic numbers); no measured input is needed.
"""
import numpy as np
import sympy as sp
from scipy.optimize import minimize_scalar

from common import check, info

# ---------------- S1: theta-bar bookkeeping ----------------
phi, f, th0, lam = sp.symbols("phi f theta0 lambda", real=True)
nq = 3
c = sp.symbols(f"c1:{nq + 1}", real=True)
rho0 = sp.symbols(f"rho1:{nq + 1}", real=True)
alpha = [ci * phi / f for ci in c]          # q -> exp(i alpha_q gamma5) q removes (d phi/f) c_q J5_q
# conventions: mass phase rho_q -> rho_q + 2 alpha_q ; Fujikawa Jacobian theta -> theta - 2 alpha_q (fundamental)
theta_a = th0 - 2 * sum(alpha)
rho_a = [r + 2 * a for r, a in zip(rho0, alpha)]
tb_a = sp.simplify(theta_a + sum(rho_a))
check("S1a derivative coupling: theta-bar = theta + arg det M is phi-independent",
      sp.diff(tb_a, phi) == 0, f"theta-bar = {tb_a}")
tb_b = sp.simplify(theta_a + sum(rho0))
check("S1b reading A (anomaly kept, mass phases dropped): theta-bar shifts by -2 sum(c) phi/f",
      sp.simplify(sp.diff(tb_b, phi) + 2 * sum(c) / f) == 0)
tb_c = sp.simplify(theta_a + lam * phi / f + sum(rho_a))
check("S1c derivative + GS term lambda phi GGtilde/f: theta-bar = theta-bar0 + lambda phi/f (shift broken)",
      sp.simplify(sp.diff(tb_c, phi) - lam / f) == 0)
rhoQ = sp.Symbol("rhoQ", real=True)
tb_d = th0 + sum(rho0) + rhoQ + phi / f    # KSVZ: Yukawa y Phi Qbar_L Q_R, Phi ~ exp(i phi/f)
check("S1d KSVZ control (phi in a mass phase, non-derivative): theta-bar depends on phi",
      sp.diff(tb_d, phi) != 0)

# ---------------- S2: KR dualization with a GS-modified Bianchi identity ----------------
# h^mu = (1/6) eps^{mu nu rho sigma} H_{nu rho sigma}; dH <-> d_mu h^mu.  First-order action
#   L = (1/2) s h.h + a h.J + phi (d.h - cGS X),   X = F Ftilde (+ G Gtilde) density.
# phi d.h = -h.d(phi) + total derivative.  s = +1 gives a healthy -(1/2)(d phi)^2 in (-,+,+,+).
eta = sp.diag(-1, 1, 1, 1)
h = sp.Matrix(sp.symbols("h0:4", real=True))
p = sp.Matrix(sp.symbols("p0:4", real=True))      # d_mu phi
J = sp.Matrix(sp.symbols("J0:4", real=True))
a, cGS, X = sp.symbols("a c_GS X", real=True)
dot = lambda u, v: (u.T * eta * v)[0]
L1 = sp.Rational(1, 2) * dot(h, h) + a * dot(h, J) - dot(h, p) - cGS * phi * X
sol = sp.solve([sp.diff(L1, hi) for hi in h], list(h), dict=True)[0]
Leff = sp.expand(L1.subs(sol))
target = sp.expand(-sp.Rational(1, 2) * dot(p, p) + a * dot(p, J) - sp.Rational(1, 2) * a ** 2 * dot(J, J)
                   - cGS * phi * X)
check("S2 KR dual: L = -(1/2)(d phi)^2 + a d phi.J5 - (a^2/2) J5^2 - c_GS phi X",
      sp.simplify(Leff - target) == 0)
check("S2b without GS (c_GS = 0) phi enters only through d phi (exact shift symmetry)",
      sp.diff(Leff.subs(cGS, 0), phi) == 0)
check("S2c with GS the action changes by c_GS * const * int X under phi -> phi + const",
      sp.simplify(sp.diff(Leff, phi) + cGS * X) == 0)

# ---------------- S3: LO chiral Lagrangian ----------------
K, mu_, md_ = 1.0, 0.47, 1.0            # generic toy values


def vac(theta, ru, rd):
    """min over a1 of -K[m_u cos(a1+ru) + m_d cos(theta-a1+rd)] (det U fixed by theta)."""
    fun = lambda a1: -K * (mu_ * np.cos(a1 + ru) + md_ * np.cos(theta - a1 + rd))
    grid = np.linspace(-np.pi, np.pi, 721)
    s = grid[np.argmin(fun(grid))]
    r = minimize_scalar(fun, bounds=(s - 0.02, s + 0.02), method="bounded", options={"xatol": 1e-13})
    return r.fun


cu, cd, fval, lamv = 0.7, 1.3, 1.0, 0.9
phis = np.linspace(-0.4, 0.4, 41)


def V(case, ph):
    al_u, al_d = cu * ph / fval, cd * ph / fval
    th = -2 * (al_u + al_d)
    if case == "a":
        return vac(th, 2 * al_u, 2 * al_d)
    if case == "b":
        return vac(th, 0.0, 0.0)
    return vac(th + lamv * ph / fval, 2 * al_u, 2 * al_d)


Va = np.array([V("a", x) for x in phis])
check("S3a chiral Lagrangian: V(phi) flat for pure derivative coupling", np.ptp(Va) < 1e-12,
      f"max-min = {np.ptp(Va):.1e}")
closed = lambda tb: -K * np.sqrt(mu_ ** 2 + md_ ** 2 + 2 * mu_ * md_ * np.cos(tb))
Vb = np.array([V("b", x) for x in phis])
check("S3b reading A: V = -K sqrt(mu^2+md^2+2 mu md cos theta-bar(phi))",
      np.max(np.abs(Vb - closed(-2 * (cu + cd) * phis / fval))) < 1e-10)
dph = 1e-3
for case, slope in (("b", -2 * (cu + cd) / fval), ("c", lamv / fval)):
    curv = (V(case, dph) - 2 * V(case, 0.0) + V(case, -dph)) / dph ** 2
    pred = K * mu_ * md_ / (mu_ + md_) * slope ** 2
    check(f"S3{case} curvature = K m_u m_d/(m_u+m_d) (d theta-bar/d phi)^2 (axion mass formula)",
          abs(curv / pred - 1) < 1e-4, f"{curv:.6f} vs {pred:.6f}")
info("Only an explicit, non-derivative phi G Gtilde (GS term, or phi in a quark-mass phase) gives a QCD mass.")
