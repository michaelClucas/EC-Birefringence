"""036 step B: SNY 0+ PGT reduced to a scalar-tensor theory; 0+ mass, matter coupling, photon coupling.

SNY (0805.3834) Lagrangian: L = -(a0/2) R + (b/24) R^2 + (a1/8)(torsion^2 combination), R the Riemann-Cartan
scalar, with only the 0+ mode (pure-trace torsion) excited. From check_rc_trace.py (A3-A5):
    R = R_LC + 2 nabla.T - (2/3) T^2,   (a1/8)(combination) = -(a1/3) T^2.
Auxiliary field: (b/24) R^2 = (b/12) chi R - (b/24) chi^2 (chi = R on shell). Omega = -a0/2 + b chi/12.
Units kappa = 1/M_P^2 = 1 where stated (then a0 = -1); a0, a1 carry dimension mass^2, b is dimensionless.

B1  torsion equation: T = -3 dOmega/(2 Omega + a1) = -(3/2) d ln Y, Y = chi + 6 mu/b, mu = a1 - a0 (= SNY Eq. 20).
B2  on-shell torsion Lagrangian = (b/8)(d chi)^2/Y  (Jordan frame: (1/2)[F R_LC - K (d chi)^2] - V,
    F = 2 Omega, K = -(b/4)/Y, V = (b/24) chi^2).
B3  Einstein-frame kinetic coefficient Z = K/F + (3/2)(F'/F)^2 = a1 b/(4 F^2 Y) > 0 iff a1 b > 0 (SNY no-ghost).
B4  mass from Einstein frame (a0 = -1): m^2 = U''(0)/Z(0) = 2 mu/(a1 b).
B5  mass from the linearized Jordan-frame trace + chi equations (general a0): box chi - m^2 chi = (2mu/(a1 b)) T_m,
    m^2 = -2 a0 mu/(a1 b); static heavy limit chi = -T_m/(-a0) = kappa rho (GR value R = 8 pi G rho).
B6  SNY FRW equations (025 sny_background.py docstring, their Eqs. 43-47) linearized at H = 0:
    R'' = -m^2 R - (2 mu/(a1 b)) T_m with the same m^2; SNY Eq. 53 period 2 pi sqrt(-a1 b/(2 a0 mu)) = 2 pi/m.
B7  f(R) limit a1 -> infinity: m^2 -> -2 a0/b = Starobinsky scalaron mass for (M_P^2/2)(R + R^2/(6 m^2)).
B8  Brans-Dicke form: omega = K F/F'^2 = -3/(2 mu) (a0 = -1); Yukawa strength alpha = 1/(3 + 2 omega) = mu/(3 a1)
    >= 1/3; Einstein-frame check alpha = 2 Q^2 with Q = F'/(2 F sqrt Z); PPN gamma - 1 = -2mu/(4mu - 3),
    |gamma - 1| > 1/2 for every a1 > 0.
B9  photon coupling: CFJ p = xi T = -(3 xi/2) d ln Y <-> -(g/4) phi F F~ with p = -g d phi (CONVENTIONS), so
    g phi = (3 xi/2) ln Y; canonical phi = sqrt(Z) delta chi -> g = xi sqrt(3/(2 mu a1)) (units M_P = 1).
B10 heavy limit: beta = -(3 xi/4) ln(Y_obs/Y_src) -> -(3 xi/4) (chi_obs - chi_src)/Y_c, Y_c = 3 a1 m^2, and in terms
    of (g, alpha, m): |beta| = sqrt(2 alpha) g |Delta rho| /(4 M_P m^2).
"""
import sympy as sp

from common import check, summary

a0, a1, b, chi, xi, m2s = sp.symbols("a0 a1 b chi xi m2", real=True)
mu = a1 - a0
Om = -a0 / 2 + b * chi / 12
Y = chi + 6 * mu / b

# B1: L_T = -2 dOmega.T - ((2 Omega + a1)/3) T^2 (after integrating 2 Omega nabla.T by parts); vary T
dchi, T = sp.symbols("dchi T")
dOm = sp.diff(Om, chi) * dchi
LT = -2 * dOm * T - (2 * Om + a1) / 3 * T**2
Tsol = sp.solve(sp.diff(LT, T), T)[0]
check("B1 torsion eq gives T = -(3/2) d chi / Y (= -(3/2) d ln Y), Y = chi + 6(a1 - a0)/b",
      sp.simplify(Tsol + sp.Rational(3, 2) * dchi / Y) == 0)
# B2
LTon = sp.simplify(LT.subs(T, Tsol))
check("B2 on-shell torsion Lagrangian = (b/8) (d chi)^2 / Y", sp.simplify(LTon - b / 8 * dchi**2 / Y) == 0)
F = 2 * Om
Kc = -(b / 4) / Y
check("B2b Jordan form -(1/2) K (dchi)^2 reproduces it", sp.simplify(-Kc / 2 * dchi**2 - LTon) == 0)
# B3
Z = sp.simplify(Kc / F + sp.Rational(3, 2) * (sp.diff(F, chi) / F)**2)
check("B3 Z = a1 b/(4 F^2 Y)", sp.simplify(Z - a1 * b / (4 * F**2 * Y)) == 0)
# B4 Einstein frame mass, a0 = -1
V = b / 24 * chi**2
U = V / F**2
m2E = sp.simplify((sp.diff(U, chi, 2) / Z).subs(chi, 0).subs(a0, -1))
check("B4 Einstein-frame m^2 = 2 mu/(a1 b) (a0 = -1)", sp.simplify(m2E - 2 * (a1 + 1) / (a1 * b)) == 0, f"({m2E})")
# B5 linearized Jordan equations, flat background, scalar sector: unknowns chi, Rlc ; box -> s symbol
s, Tm = sp.symbols("s T_m")  # s stands for box acting on chi (linear, plane-wave)
# chi eq: Omega' R_LC + K box chi - V' = 0 ; trace eq: -Omega R_LC + 3 box Omega = T_m/2
K0 = Kc.subs(chi, 0)
Rlc = sp.Symbol("Rlc")
e1 = sp.diff(Om, chi) * Rlc + K0 * s * chi - sp.diff(V, chi)
e2 = -Om.subs(chi, 0) * Rlc + 3 * sp.diff(Om, chi) * s * chi - Tm / 2
sol = sp.solve([e1, e2], [Rlc, Tm], dict=True)[0]
trace_eq = sp.expand(sp.simplify(sol[Tm]))  # T_m expressed through (s - m^2) chi
coef = sp.simplify(trace_eq.coeff(chi))
# coef = c (s - m2)
c_s = sp.simplify(sp.diff(coef, s))
m2J = sp.simplify(-(coef - c_s * s) / c_s)
check("B5a Jordan trace route m^2 = -2 a0 mu/(a1 b)", sp.simplify(m2J + 2 * a0 * mu / (a1 * b)) == 0, f"({sp.factor(m2J)})")
check("B5b source: box chi - m^2 chi = (2 mu/(a1 b)) T_m", sp.simplify(1 / c_s - 2 * mu / (a1 * b)) == 0)
check("B5c static heavy limit chi = -T_m/(-a0) (GR: R = kappa rho)",
      sp.simplify((-(2 * mu / (a1 * b)) * Tm / m2J) - (-Tm / (-a0))) == 0)
check("B4=B5 at a0 = -1", sp.simplify(m2J.subs(a0, -1) - m2E) == 0)
# B6 SNY FRW linearization (a0 general in the Phi equation)
Rf, Ph, Tr = sp.symbols("R Phi Tr")
Rdot = -sp.Rational(2, 3) * (Rf + 6 * mu / b) * Ph
Phdot_lin = -a0 * Rf / (2 * a1) + Tr / (2 * a1)  # H = 0, Phi^2 dropped
Rddot = sp.expand(sp.diff(Rdot, Ph).subs(Rf, 0) * Phdot_lin)  # linear: R' = -(4mu/b) Phi at R ~ 0
check("B6a SNY FRW: R'' = -m^2 R - (2mu/(a1 b)) T_m", sp.simplify(Rddot - (-m2J * Rf - 2 * mu / (a1 * b) * Tr)) == 0)
per = 2 * sp.pi * sp.sqrt(-a1 * b / (2 * a0 * mu))
check("B6b SNY Eq. 53 period = 2 pi/m", sp.simplify(per**2 - 4 * sp.pi**2 / m2J) == 0)
# B7
check("B7 f(R) limit a1 -> oo: m^2 -> -2 a0/b (Starobinsky with b/24 = M_P^2/(12 m^2), M_P^2 = -a0)",
      sp.simplify(sp.limit(m2J, a1, sp.oo) - (-2 * a0 / b)) == 0)
# B8 BD, a0 = -1
Fp = sp.diff(F, chi)
omega = sp.simplify((Kc * F / Fp**2).subs(chi, 0).subs(a0, -1))
check("B8a omega_BD = -3/(2 mu)", sp.simplify(omega + sp.Rational(3, 2) / (a1 + 1)) == 0, f"({omega})")
alpha = sp.simplify(1 / (3 + 2 * omega))
check("B8b alpha = mu/(3 a1)", sp.simplify(alpha - (a1 + 1) / (3 * a1)) == 0)
Q = (Fp / (2 * F * sp.sqrt(Z))).subs(chi, 0).subs(a0, -1)
check("B8c Einstein-frame alpha = 2 Q^2", sp.simplify(2 * Q**2 - alpha) == 0)
gm1 = sp.simplify(-1 / (2 + omega))
check("B8d gamma - 1 = -2mu/(4mu - 3)", sp.simplify(gm1 + 2 * (a1 + 1) / (4 * (a1 + 1) - 3)) == 0)
ap = sp.Symbol("ap", positive=True)
check("B8e |gamma - 1| - 1/2 > 0 for all a1 > 0, -> 0 as a1 -> oo",
      sp.simplify(-gm1.subs(a1, ap) - sp.Rational(1, 2)) == sp.simplify(3 / (2 * (4 * ap + 1)))
      and sp.limit(-gm1.subs(a1, ap), ap, sp.oo) == sp.Rational(1, 2))
check("B8f alpha >= 1/3 for a1 > 0", sp.simplify(alpha.subs(a1, ap) - sp.Rational(1, 3) - 1 / (3 * ap)) == 0)
# B9 photon coupling, a0 = -1
Z0 = Z.subs(chi, 0).subs(a0, -1)
Y0 = Y.subs(chi, 0).subs(a0, -1)
g = sp.Rational(3, 2) * xi / (Y0 * sp.sqrt(Z0))
check("B9 g = xi sqrt(3/(2 mu a1))", sp.simplify(g**2 - xi**2 * 3 / (2 * (a1 + 1) * a1)) == 0)
# B10
Yc = Y0
check("B10a Y_c = 3 a1 m^2 (a0 = -1)", sp.simplify(Yc - 3 * a1 * m2E) == 0)
drho = sp.Symbol("drho")
beta = sp.Rational(3, 4) * xi * drho / Yc  # magnitude, chi = kappa rho = rho (M_P = 1)
gsym = sp.Symbol("g")
beta_g = beta.subs(xi, sp.solve(sp.Eq(g, gsym), xi)[0])
check("B10b |beta| = sqrt(2 alpha) g Delta rho/(4 m^2) (M_P = 1)",
      sp.simplify(beta_g**2 - (2 * alpha * gsym**2 * drho**2 / (16 * m2E**2))) == 0)
summary()
