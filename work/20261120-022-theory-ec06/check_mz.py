"""Magueijo-Zlosnik parity-violating FRW torsion (1908.05184) in our conventions (task 022, part a).

MZ use forms: T^0 = 0, T^i = -T(t) e^0 e^i + P(t) eps_ijk e^j e^k, e^0 = dt, e^i = a dx^i, and
g = adot/a + T (their Eq. 19), so the metric Hubble rate is H = g - T.
M1  frame components -> our trace and axial vectors: T_t = -3T, S^0hat = -2P, S_t = +2P (cosmic time, lower).
M2  gamma -> infinity parity-odd branch (Eqs. 32, 34, 35 => 37): T = (g^2 - P^2)/(2g), H = (g^2 + P^2)/(2g).
M3  finite gamma (Eqs. 26, 28, 29), P != 0: T/g = gamma x (1 - x^2 - 2x/gamma)/(1 - x^2 + 2 gamma x), x = P/g;
    -> M2 as gamma -> infinity; T ~ gamma P for x << 1/gamma.
M4  tracking (Eqs. 57-59, 64-66, R_Lambda = 0): x^2 = R_c/(1+R_c); Eq. 66 gives R_c ~ gamma^2/9 for small gamma
    and R_c -> 1/2 for gamma -> infinity (MZ Fig. 1).
M5  effective Friedmann ratio H^2/(kappa rho/3) at Lambda = 0: (1 - T/g)^2/(1 - x^2) -> 1 - gamma^2/9 + O(gamma^3) (a deficit, like Delta N_eff < 0).
"""
import sympy as sp

from common import R4, check, summary

T_, P_, g_, rho, Lam, Lamdot, gam, x, k = sp.symbols("T P g rho Lambda Lambdadot gamma x kappa", real=True)

# M1: frame torsion tensor T^A_{BC} (orthonormal), eta = diag(-1,1,1,1), eps_{0123} = +1
Tf = [[[sp.Integer(0)] * 4 for _ in R4] for _ in R4]
for i in (1, 2, 3):
    Tf[i][0][i] = -T_
    Tf[i][i][0] = T_
    for j in (1, 2, 3):
        for kk in (1, 2, 3):
            Tf[i][j][kk] = 2 * P_ * sp.LeviCivita(i, j, kk)  # P eps_ijk e^j^e^k = (1/2)(2P eps_ijk) e^j^e^k
eta = sp.diag(-1, 1, 1, 1)
trace = [sum(Tf[l][m][l] for l in R4) for m in R4]
Tl = [[[eta[l, l] * Tf[l][m][n] for n in R4] for m in R4] for l in R4]
S_up = [sp.Rational(1, 6) * sum(-sp.LeviCivita(mu, n, r, s) * Tl[n][r][s] for n in R4 for r in R4 for s in R4)
        for mu in R4]
check("M1a MZ trace: T_t = -3T, T_i = 0", sp.simplify(trace[0] + 3 * T_) == 0 and all(t == 0 for t in trace[1:]))
check("M1b MZ axial: S^0hat = -2P, S^i = 0 (so S_t = +2P)", sp.simplify(S_up[0] + 2 * P_) == 0
      and all(s == 0 for s in S_up[1:]))

# M2: gamma -> infinity branch
sol37 = sp.solve(sp.Eq(1, 3 * Lamdot / Lam**2 * g_), Lamdot)[0]  # Eq. 37
T34 = (Lamdot / (2 * Lam)) * (1 + k * rho / Lam)  # Eq. 34
krho = 3 * (g_**2 - P_**2) - Lam  # Eq. 32
Tinf = sp.simplify(T34.subs(Lamdot, sol37).subs(k * rho, krho).subs(rho, krho / k))
check("M2a gamma->inf: T = (g^2 - P^2)/(2g)", sp.simplify(Tinf - (g_**2 - P_**2) / (2 * g_)) == 0, f"(T = {Tinf})")
check("M2b gamma->inf: H = g - T = (g^2 + P^2)/(2g)", sp.simplify(g_ - Tinf - (g_**2 + P_**2) / (2 * g_)) == 0)

# M3: finite gamma
E29 = sp.Eq(P_, 3 * Lamdot / Lam**2 * (g_ * P_ + (Lam + k * rho) / (6 * gam)))
ld = sp.solve(E29, Lamdot)[0]
T28 = Lamdot / (2 * Lam**2) * (Lam + k * rho - 6 * g_ * P_ / gam)
Tfin = sp.simplify(T28.subs(Lamdot, ld).subs(rho, (3 * (g_**2 - P_**2) - Lam) / k))
formula = g_ * gam * x * (1 - x**2 - 2 * x / gam) / (1 - x**2 + 2 * gam * x)
check("M3a finite gamma: T/g = gamma x (1-x^2-2x/gamma)/(1-x^2+2 gamma x)",
      sp.simplify(Tfin.subs(P_, x * g_) - formula) == 0)
check("M3b finite-gamma T -> gamma->inf result", sp.simplify(sp.limit(formula, gam, sp.oo) - g_ * (1 - x**2) / 2) == 0)
check("M3c small x: T = gamma P + O(x^2)", sp.simplify(sp.series(formula / g_, x, 0, 2).removeO() - gam * x) == 0)

# M4: Eq. 66 small gamma
Rc = sp.symbols("R_c", positive=True)
E66 = 2 * gam * Rc**2 + gam * Rc + 2 * sp.sqrt(Rc + 1) * Rc**sp.Rational(3, 2) + 3 * sp.sqrt(Rc + 1) * sp.sqrt(Rc) - gam
u = sp.symbols("u", positive=True)  # u = sqrt(R_c)
E66u = sp.simplify(E66.subs(Rc, u**2))
# leading order: 3u - gamma = 0 -> u = gamma/3
uu = sp.symbols("uu")
lead = sp.series(E66u.subs(u, gam * uu), gam, 0, 2).removeO()
usol = sp.solve(sp.simplify(lead / gam), uu)
check("M4a Eq. 66: R_c = gamma^2/9 at leading order", sp.Rational(1, 3) in usol, f"(sqrt(R_c)/gamma = {usol})")
big = sp.limit(E66 / gam, gam, sp.oo)
check("M4b Eq. 66: R_c -> 1/2 as gamma -> infinity (2R_c^2 + R_c - 1 = 0)",
      sp.Rational(1, 2) in sp.solve(big, Rc))

# M5: Friedmann ratio at Lambda = 0 along tracking, small gamma
xs = gam / 3
ratio = (1 - formula / g_) ** 2 / (1 - x**2)
ser = sp.series(ratio.subs(x, xs), gam, 0, 3).removeO()
check("M5 H^2/(kappa rho/3) = 1 - gamma^2/9 + O(gamma^3) on the small-gamma tracker",
      sp.simplify(ser - (1 - gam**2 / 9)) == 0, f"(series {sp.simplify(ser)})")
summary()
