"""Teleparallel FRW torsion (task 022, part b).

T1  flat FRW, diagonal Cartesian tetrad (omega = 0): T_tau = 3a'/a, S = 0, T closed.
T2  TEGR frame dependence: a local boost of the same tetrad (same metric) changes T_mu by a
    non-closed vector, so p = xi T_mu is not determined by the metric in TEGR.
T3  closed FRW (k = +1), tetrads from left-/right-invariant forms on S^3 (the two 'axial' branches):
    metric is a^2(-dtau^2 + unit S^3), T_tau = 3a'/a, S^0hat = a S^tau = +-2 sqrt(k)/a, closed on FRW.
T5  scalar perturbations, diagonal tetrad e^0 = a(1+Phi)dtau, e^i = a(1-Psi)dx^i (linear order):
    T_tau = d_tau(3 ln a - 3 Psi), T_i = d_i(Phi - 2 Psi), so dT = d_i d_tau(Phi + Psi) != 0 unless
    Phi + Psi is static (true in matter domination, broken in the radiation and Lambda eras: ISW).
T4  f(T) antisymmetric field equations S_[mu nu]^tau dT/dtau = 0 hold for T1 and T3 ('good' tetrads);
    torsion scalar values.
"""
import random

import sympy as sp

from common import R4, axial_vec, check, is_closed, summary, trace_vec, weitzenbock_torsion

tau, x, y, z = sp.symbols("tau x y z", real=True)
a = sp.Function("a", positive=True)(tau)
ap = sp.diff(a, tau)

# ---------------------------------------------------------------- T1 flat diagonal tetrad
X = [tau, x, y, z]
e = sp.diag(a, a, a, a)
g = sp.diag(-a**2, a**2, a**2, a**2)
T = weitzenbock_torsion(e, X)
Tv = trace_vec(T)
S = axial_vec(T, g)
check("T1a flat FRW trace T_tau = 3a'/a", sp.simplify(Tv[0] - 3 * ap / a) == 0 and all(v == 0 for v in Tv[1:]))
check("T1b flat FRW axial S = 0", all(s == 0 for s in S))
check("T1c T_mu closed on FRW (T_mu = 3 d ln a)", is_closed(Tv, X))

# ---------------------------------------------------------------- T2 local Lorentz boost (TEGR gauge)
eta = sp.Function("eta")(tau, x, y)
L = sp.eye(4)
L[0, 0] = L[1, 1] = sp.cosh(eta)
L[0, 1] = L[1, 0] = sp.sinh(eta)
e2 = L * e
g2 = sp.simplify(e2.T * sp.diag(-1, 1, 1, 1) * e2)
check("T2a boosted tetrad gives the same metric", sp.simplify(g2 - g) == sp.zeros(4))
T2 = weitzenbock_torsion(e2, X)
Tv2 = trace_vec(T2)
dT = [sp.simplify(Tv2[m] - Tv[m]) for m in R4]
print("   Delta T_mu under boost eta(tau,x,y) along x:", dT)
check("T2b TEGR: trace torsion changes under a local boost (frame dependent)", any(d != 0 for d in dT))
check("T2c the boosted p = xi T_mu is not closed (generic eta)", not is_closed(Tv2, X))
S2 = axial_vec(T2, g)
print("   axial S^mu after boost:", [sp.simplify(s) for s in S2])
check("T2d a local boost also generates axial torsion", any(sp.simplify(s) != 0 for s in S2))

# ---------------------------------------------------------------- T3 closed FRW, S^3 invariant forms
psi, th, ph = sp.symbols("psi theta phi", real=True)
Y = [tau, psi, th, ph]


def forms(left=True):
    u, w = (psi, ph) if left else (ph, psi)
    s1 = {th: sp.sin(u), w: -sp.cos(u) * sp.sin(th)}
    s2 = {th: sp.cos(u), w: sp.sin(u) * sp.sin(th)}
    s3 = {u: 1, w: sp.cos(th)}
    return [s1, s2, s3]


def tetrad(left):
    E = sp.zeros(4, 4)
    E[0, 0] = a
    for i, s in enumerate(forms(left)):
        for c, v in s.items():
            E[i + 1, Y.index(c)] = a * v / 2
    return E


g3 = None
res = {}
for left in (True, False):
    E = tetrad(left)
    gm = sp.simplify(E.T * sp.diag(-1, 1, 1, 1) * E)
    if g3 is None:
        g3 = gm
    else:
        check("T3a left- and right-invariant tetrads give the same metric", sp.simplify(gm - g3) == sp.zeros(4))
    Tt = weitzenbock_torsion(E, Y)
    tv = trace_vec(Tt)
    sv = [sp.simplify(v.subs(sp.Abs(sp.sin(th)), sp.sin(th))) for v in axial_vec(Tt, gm)]  # 0 < theta < pi
    res[left] = (Tt, tv, sv, gm, E)
    name = "left" if left else "right"
    check(f"T3b k=+1 ({name}) trace T_tau = 3a'/a, spatial T_i = 0",
          sp.simplify(tv[0] - 3 * ap / a) == 0 and all(sp.simplify(v) == 0 for v in tv[1:]))
    S0hat = sp.simplify(a * sv[0])
    print(f"   k=+1 {name}: S^0hat = a S^tau = {S0hat}; spatial S^i = {[sp.simplify(s) for s in sv[1:]]}")
    check(f"T3c k=+1 ({name}) |S^0hat| = 2/a (= 2 sqrt(k)/a, comoving radius 1), S^i = 0",
          sp.simplify(S0hat**2 - 4 / a**2) == 0 and all(sp.simplify(s) == 0 for s in sv[1:]))
    check(f"T3d k=+1 ({name}) p = xi S_mu closed on FRW",
          is_closed([sp.simplify(sum(gm[m, n] * sv[n] for n in R4)) for m in R4], Y))
check("T3e the two branches have opposite axial torsion (sign = branch choice, not fixed by the theory)",
      sp.simplify(res[True][2][0] + res[False][2][0]) == 0)

# spatial metric is the unit round S^3: Ricci scalar of h_ij = g_ij / a^2 equals 6 (k = +1)
h = sp.simplify(g3[1:, 1:] / a**2)
Z = [psi, th, ph]
hi = h.inv()
Chr = [[[sp.simplify(sum(hi[l, r] * (sp.diff(h[r, m], Z[n]) + sp.diff(h[r, n], Z[m]) - sp.diff(h[m, n], Z[r]))
                         for r in range(3)) / 2) for n in range(3)] for m in range(3)] for l in range(3)]
Ric = sp.zeros(3)
for m in range(3):
    for n in range(3):
        Ric[m, n] = sum(sp.diff(Chr[l][m][n], Z[l]) - sp.diff(Chr[l][m][l], Z[n])
                        + sum(Chr[l][l][r] * Chr[r][m][n] - Chr[l][n][r] * Chr[r][m][l] for r in range(3))
                        for l in range(3))
Rs = sp.simplify(sum(hi[m, n] * Ric[m, n] for m in range(3) for n in range(3)))
check("T3f spatial metric is the unit S^3 (3-Ricci scalar = 6k, k = 1)", sp.simplify(Rs - 6) == 0, f"(R3 = {Rs})")


# ---------------------------------------------------------------- T4 superpotential, f(T) antisymmetric eqs
def superpot(Tt, gm):
    """Aldrovandi-Pereira superpotential with T^r_{mn} of CONVENTIONS.md:
    K^{mn}_r = -1/2 (T^{mn}_r - T^{nm}_r - T_r^{mn}),  T^{mn}_r := g^{nb} T^m_{br},
    S_r^{mn} = 1/2 (K^{mn}_r + delta^m_r T^{an}_a - delta^n_r T^{am}_a),  T^{an}_a = g^{nb} T_b (trace),
    torsion scalar = S_r^{mn} T^r_{mn}.  Overall sign conventions do not affect the checks below."""
    gi = gm.inv()
    tv = trace_vec(Tt)
    tvu = [sum(gi[n, b] * tv[b] for b in R4) for n in R4]
    A = [[[sum(gi[n, b] * Tt[m][b][r] for b in R4) for r in R4] for n in R4] for m in R4]  # A[m][n][r] = T^{mn}_r
    B = [[[sum(gm[r, al] * gi[m, be] * gi[n, ga] * Tt[al][be][ga] for al in R4 for be in R4 for ga in R4)
           for n in R4] for m in R4] for r in R4]  # B[r][m][n] = T_r^{mn}
    Sp = [[[None] * 4 for _ in R4] for _ in R4]
    for r in R4:
        for m in R4:
            for n in R4:
                K = -sp.Rational(1, 2) * (A[m][n][r] - A[n][m][r] - B[r][m][n])
                Sp[r][m][n] = sp.Rational(1, 2) * (K + int(r == m) * tvu[n] - int(r == n) * tvu[m])
    Tsc = sp.simplify(sum(Sp[r][m][n] * Tt[r][m][n] for r in R4 for m in R4 for n in R4))
    return Sp, Tsc


for label, (Tt, gm) in {"flat": (T, g), "k=+1 left": (res[True][0], res[True][3]),
                        "k=+1 right": (res[False][0], res[False][3])}.items():
    Sp, Tsc = superpot(Tt, gm)
    print(f"   {label}: torsion scalar = {sp.simplify(Tsc)}")
    # antisymmetric f(T) equation: f_TT dT_scalar/dtau * Y_[mu nu] = 0 with Y_{mu nu} = g_{nu a} S_mu^{a tau}
    Ymat = sp.Matrix(4, 4, lambda m, n: sp.simplify(sum(gm[n, k] * Sp[m][k][0] for k in R4)))
    check(f"T4 {label}: antisymmetric part of S_mu^(nu tau) vanishes (good tetrad for any f(T))",
          sp.simplify(Ymat - Ymat.T) == sp.zeros(4))
# teeth: the boosted flat tetrad of T2 (same metric) is NOT a good tetrad for f(T) with f_TT != 0
Sp, Tsc = superpot(T2, g)
dTs = [sp.diff(Tsc, v) for v in X]
Yb = sp.Matrix(4, 4, lambda m, n: sp.simplify(sum(g[n, k] * Sp[m][k][r] * dTs[r] for k in R4 for r in R4)))
print(f"   boosted flat tetrad: torsion scalar = {sp.simplify(Tsc)}")
check("T4 boosted tetrad violates the f(T) antisymmetric equation (check has teeth)",
      sp.simplify(Yb - Yb.T) != sp.zeros(4))
# T5 linear scalar perturbations
epsl = sp.Symbol("epsilon")
Phi = sp.Function("Phi")(tau, x, y, z)
Psi = sp.Function("Psi")(tau, x, y, z)
ep = sp.diag(a * (1 + epsl * Phi), a * (1 - epsl * Psi), a * (1 - epsl * Psi), a * (1 - epsl * Psi))
Tp = weitzenbock_torsion(ep, X)
tvp = [sp.series(v, epsl, 0, 2).removeO() for v in trace_vec(Tp)]
exp_ = [3 * ap / a - 3 * epsl * sp.diff(Psi, tau)] + [epsl * sp.diff(Phi - 2 * Psi, X[i]) for i in (1, 2, 3)]
check("T5a perturbed trace: T_tau = d(3 ln a - 3Psi), T_i = d_i(Phi - 2Psi) (linear)",
      all(sp.simplify(tvp[m] - exp_[m]) == 0 for m in R4))
curl = sp.simplify(sp.diff(tvp[1], tau) - sp.diff(tvp[0], x))
check("T5b (dT)_{tau x} = d_x d_tau (Phi + Psi): closed only if Phi + Psi is static",
      sp.simplify(curl - epsl * sp.diff(Phi + Psi, tau, x)) == 0, f"({curl})")
summary()
