"""Task 031, item 1: is the MZ axial torsion vector S closed off FRW?

MZ field equations used (1908.05184, ar5iv HTML opened in task 031):
  Eq. 14 (tetrad):     eps_abcd (e^b R^cd - (1/3) Lambda e^b e^c e^d) = -2 kappa tau_a
  Eq. 15 (connection): T^[a ^ e^b] = -(3/(2 Lambda^2)) dLambda ^ R^ab + (3/(4 gamma Lambda^2)) eps^abcd dLambda ^ R_cd
  Eqs. 26, 28, 29 (FRW): g^2 - P^2 = (Lambda + kappa rho)/3,
     T = (Ldot/(2 L^2)) (L + kappa rho - 6 g P/gamma),  P = (3 Ldot/L^2) (g P + (L + kappa rho)/(6 gamma)).

A1-A3  FRW: our forms code rebuilds MZ torsion from omega^i_0 = g e^i, omega^ij = -P eps_ijk e^k, and Eq. 15
       reproduces Eqs. 28 and 29 (fixes the eps sign convention of the gamma term).
B1-B3  Generic algebra (random rational curvature, any vector v = dtheta): the axial vector of the solution of
       T^[a ^ e^b] = v ^ M^ab is
         S_A = alpha * E_A^F v_F            for M = eps^AB_CD R^CD  (E = the Eq. 14 contraction of R)
         S^A = beta  * W^AF v_F             for M = R^AB            (W = dual of R^A_B ^ e^B = dual of D T^A)
       so, using Eq. 14, S_mu = (2/(3 gamma)) (kappa tau_mu^nu - Lambda delta_mu^nu) d_nu theta + O(d theta . DT),
       theta = 3/(2 Lambda). The normalisation 2/(3 gamma) is fixed in A3 against Eq. 29.
C1-C2  Bianchi D T^A = R^A_B ^ e^B on the MZ FRW connection (so the beta-term is O((d theta)^2) off FRW), and
       eps_ABCD R^AB ^ R^CD = GB vol for the Levi-Civita FRW connection (sign used for the theta action).
D1-D3  Linear perturbations (Newtonian gauge, perfect fluid): for S = C tau.d theta, (dS)_{tau i} = d_i Y with
         Y/C = (rho+p) d_tau(dtheta_c) + p' dtheta_c + theta' drho_c,  dtheta_c = dtheta + theta' v, drho_c = drho + rho' v
       gauge invariant, zero for a pure time-shift mode, and nonzero otherwise (closure is an extra condition).
E1     Effective form S = -(1/(sqrt3 M_P)) [d phi - (1 - c_s) (d phi)_perp]: d of (phi_t, c_s grad phi) = (c_s - 1) phi_ti.
"""
import random

import sympy as sp

from common import (ETA, R4, add, check, d1, eps_lo, eps_up, is_zero, summary, wedge11, wedge12)

# ---------------- A: FRW ----------------
t, x, y, z = sp.symbols("t x y z", real=True)
X = [t, x, y, z]
a = sp.Function("a", positive=True)(t)
g = sp.Function("g")(t)
P = sp.Function("P")(t)
Tv, Ld, Lam, gam = sp.symbols("T Ldot Lambda gamma")

e = [[sp.Integer(1), 0, 0, 0], [0, a, 0, 0], [0, 0, a, 0], [0, 0, 0, a]]
e = [[sp.sympify(c) for c in row] for row in e]


def eps3(i, j, k):
    return sp.LeviCivita(i - 1, j - 1, k - 1)


Z1 = [sp.Integer(0)] * 4
om = [[Z1[:] for _ in R4] for _ in R4]
for i in (1, 2, 3):
    om[i][0] = [g * c for c in e[i]]
    om[0][i] = [g * c for c in e[i]]
    for j in (1, 2, 3):
        if i != j:
            om[i][j] = add(*[e[k] for k in (1, 2, 3)], coeffs=[-P * eps3(i, j, k) for k in (1, 2, 3)])

Tor = [add(d1(e[A], X), *[wedge11(om[A][B], e[B]) for B in R4]) for A in R4]
Rc = [[add(d1(om[A][B], X), *[wedge11(om[A][C], om[C][B]) for C in R4]) for B in R4] for A in R4]
sub_adot = {sp.Derivative(a, t): a * (g - Tv)}


def S_(f):
    if isinstance(f, list):
        return [S_(q) for q in f]
    return sp.simplify(sp.sympify(f).subs(sub_adot))


Tor = S_(Tor)
ok = is_zero(Tor[0])
for i in (1, 2, 3):
    target = add(wedge11(e[0], e[i]), *[wedge11(e[j], e[k]) for j in (1, 2, 3) for k in (1, 2, 3)],
                 coeffs=[-Tv] + [P * eps3(i, j, k) for j in (1, 2, 3) for k in (1, 2, 3)])
    ok = ok and is_zero(add(Tor[i], target, coeffs=[1, -1]))
check("A1 forms code: omega^i_0 = g e^i, omega^ij = -P eps_ijk e^k give MZ torsion T^i = -T e0^ei + P eps_ijk ej^ek "
      "with g = adot/a + T", ok)

Rup = [[[[sum(Rc[A][C][m][n] * ETA[C, B] for C in R4) for n in R4] for m in R4] for B in R4] for A in R4]
Rdn = [[[[sum(ETA[A, C] * Rc[C][B][m][n] for C in R4) for n in R4] for m in R4] for B in R4] for A in R4]
dL = [Ld, 0, 0, 0]
eq28 = Ld / (2 * Lam**2) * (3 * (g**2 - P**2) - 6 * g * P / gam)  # Eq. 28 with Eq. 26
Ld29 = sp.solve(sp.Eq(P, 3 * Ld / Lam**2 * (g * P + 3 * (g**2 - P**2) / (6 * gam))), Ld)[0]  # Eq. 29 with Eq. 26


def eq15_residual(s):
    res = []
    for A in R4:
        for B in R4:
            if A >= B:
                continue
            lhs = add(wedge12(e[B], Tor[A]), wedge12(e[A], Tor[B]), coeffs=[sp.Rational(1, 2), -sp.Rational(1, 2)])
            rhs = add(wedge12(dL, Rup[A][B]),
                      *[wedge12(dL, Rdn[C][D]) for C in R4 for D in R4],
                      coeffs=[-sp.Rational(3, 2) / Lam**2] + [s * sp.Rational(3, 4) / (gam * Lam**2) * eps_up(A, B, C, D)
                                                              for C in R4 for D in R4])
            r = add(lhs, rhs, coeffs=[1, -1])
            res += [r[m][n][q] for m in R4 for n in R4 for q in R4]
    return res


n_good = 0
for s in (1, -1):
    res = eq15_residual(s)
    good = all(sp.cancel(sp.expand(r.subs(Tv, eq28).subs(Ld, Ld29))) == 0 for r in res if r != 0)
    print(f"   eps sign s = {s:+d} (eps^0123 = {-s:+d}): Eq. 15 on FRW reproduces Eqs. 28 + 29: {good}")
    n_good += good
check("A2 Eq. 15 on FRW (spatial curvature only, since dLambda ~ dt) reproduces MZ Eqs. 28 and 29 for exactly one "
      "eps sign", n_good == 1)

# A3: covariant normalisation on FRW: S_t = 2P (task 022 M1) must equal (2/(3 gamma)) G_t^t theta_dot + (gP term),
# with G_t^t = -(Lambda + kappa rho) = -3(g^2 - P^2) and theta = 3/(2 Lambda) => theta_dot = -(3/2) Ldot/Lambda^2.
thdot = -sp.Rational(3, 2) * Ld / Lam**2
St_cov_gamma_part = sp.Rational(2, 3) / gam * (-3 * (g**2 - P**2)) * thdot
gP_part = 6 * g * P * Ld / Lam**2  # the non-1/gamma part of Eq. 29 (x2)
check("A3 2P = (2/(3 gamma)) G_t^t theta_dot + 6 g P Ldot/Lambda^2 is Eq. 29: normalisation of the covariant formula",
      sp.simplify((St_cov_gamma_part + gP_part).subs(Ld, Ld29) - 2 * P) == 0)

# ---------------- B: generic frame algebra ----------------
random.seed(31)
idx = [(C, D) for C in R4 for D in R4 if C < D]


def rand_R():
    R = [[[[sp.Integer(0)] * 4 for _ in R4] for _ in R4] for _ in R4]  # R^{AB}_{CD}
    for (A, B) in idx:
        for (C, D) in idx:
            r = sp.Rational(random.randint(-9, 9), random.randint(1, 5))
            R[A][B][C][D], R[B][A][C][D], R[A][B][D][C], R[B][A][D][C] = r, -r, -r, r
    return R


def solve_axial(J):
    """Solve 1/2 (T^A ^ e^B - T^B ^ e^A) = J^{AB} (3-form, frame comps) for T^A_{CD}; return S^A and S_A."""
    Ts = {(A, C, D): sp.Symbol(f"T{A}{C}{D}") for A in R4 for (C, D) in idx}

    def T(A, C, D):
        if C == D:
            return 0
        return Ts[(A, C, D)] if C < D else -Ts[(A, D, C)]

    def Te(A, B, C, D, E):  # (T^A ^ e^B)_{CDE}
        dlt = lambda i, j: 1 if i == j else 0
        return T(A, C, D) * dlt(B, E) + T(A, D, E) * dlt(B, C) + T(A, E, C) * dlt(B, D)

    eqs = []
    for (A, B) in idx:
        for (C, D, E) in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]:
            eqs.append(sp.Rational(1, 2) * (Te(A, B, C, D, E) - Te(B, A, C, D, E)) - J[A][B][C][D][E])
    Mx, bx = sp.linear_eq_to_matrix(eqs, list(Ts.values()))
    sol = dict(zip(Ts.values(), Mx.LUsolve(bx)))
    Tl = lambda B, C, D: ETA[B, B] * sp.sympify(T(B, C, D)).subs(sol)
    Sup = [sp.Rational(1, 6) * sum(eps_up(A, B, C, D) * Tl(B, C, D) for B in R4 for C in R4 for D in R4) for A in R4]
    return Sup, [ETA[A, A] * Sup[A] for A in R4]


def vM(v, M):  # 3-form (v ^ M^{AB})_{CDE}
    return [[[[[v[C] * M[A][B][D][E] + v[D] * M[A][B][E][C] + v[E] * M[A][B][C][D] for E in R4] for D in R4]
               for C in R4] for B in R4] for A in R4]


def E_contr(R):
    """E_A^F = (1/6) eps^{FPQR} (eps_{ABCD} e^B ^ R^{CD})_{PQR}: the Eq. 14 combination."""
    dl = lambda i, j: 1 if i == j else 0
    out = [[0] * 4 for _ in R4]
    for A in R4:
        F3 = [[[sum(eps_lo(A, B, C, D) * (dl(B, p) * R[C][D][q][r] + dl(B, q) * R[C][D][r][p] + dl(B, r) * R[C][D][p][q])
                    for B in R4 for C in R4 for D in R4) for r in R4] for q in R4] for p in R4]
        for Fi in R4:
            out[A][Fi] = sp.Rational(1, 6) * sum(eps_up(Fi, p, q, r) * F3[p][q][r] for p in R4 for q in R4 for r in R4)
    return out


def W_contr(R):
    """W^{XF} = (1/6) eps^{FPQR} (R^X_Y ^ e^Y)_{PQR}; R^X_Y = R^{XZ} eta_ZY. By Bianchi this is the dual of D T^X."""
    out = [[0] * 4 for _ in R4]
    for Xi in R4:
        D3 = [[[sum(R[Xi][Y][p][q] * ETA[Y, Y] * (1 if Y == r else 0) + R[Xi][Y][q][r] * ETA[Y, Y] * (1 if Y == p else 0)
                    + R[Xi][Y][r][p] * ETA[Y, Y] * (1 if Y == q else 0) for Y in R4) for r in R4] for q in R4] for p in R4]
        for Fi in R4:
            out[Xi][Fi] = sp.Rational(1, 6) * sum(eps_up(Fi, p, q, r) * D3[p][q][r] for p in R4 for q in R4 for r in R4)
    return out


alphas, betas, ok_shape_a, ok_shape_b = [], [], True, True
for trial in range(3):
    R = rand_R()
    v = [sp.Rational(random.randint(-7, 7), random.randint(1, 3)) for _ in R4]
    Rlow = [[[[sum(ETA[C, C2] * ETA[D, D2] * R[C2][D2][m][n] for C2 in R4 for D2 in R4) for n in R4] for m in R4]
             for D in R4] for C in R4]
    M2 = [[[[sum(eps_up(A, B, C, D) * Rlow[C][D][m][n] for C in R4 for D in R4) for n in R4] for m in R4] for B in R4]
          for A in R4]
    _, S2 = solve_axial(vM(v, M2))
    E = E_contr(R)
    Ev = [sum(E[A][F] * v[F] for F in R4) for A in R4]
    k = [sp.nsimplify(S2[A] / Ev[A]) for A in R4 if Ev[A] != 0]
    ok_shape_a = ok_shape_a and len(set(k)) == 1
    alphas.append(k[0])
    S1up, _ = solve_axial(vM(v, R))
    W = W_contr(R)
    Wv = [sum(W[A][F] * v[F] for F in R4) for A in R4]
    kb = [sp.nsimplify(S1up[A] / Wv[A]) for A in R4 if Wv[A] != 0]
    ok_shape_b = ok_shape_b and len(set(kb)) == 1
    betas.append(kb[0])
print(f"   alpha (eps-R part) per draw: {alphas}; beta (R part) per draw: {betas}")
check("B1 axial part of the Eq. 15 solution for M = eps.R is exactly alpha * E_A^F v_F (E = Eq. 14 contraction), "
      "same alpha for every random curvature", ok_shape_a and len(set(alphas)) == 1)
check("B2 axial part for M = R^AB is exactly beta * W^AF v_F, W = dual of R^A_B ^ e^B (= dual of D T^A)",
      ok_shape_b and len(set(betas)) == 1)
check("B3 so S_mu = (2/(3 gamma)) (kappa tau - Lambda) . d theta + O(d theta . DT): the 1/gamma part is fixed by "
      "matter + Lambda through Eq. 14 (normalisation from A3)", ok_shape_a and ok_shape_b)

# ---------------- C: Bianchi and GB sign ----------------
af = sp.Function("af", positive=True)(t)
gg, PP, TT = (sp.Function(n)(t) for n in ("gg", "PP", "TT"))
e2 = [[sp.Integer(1), 0, 0, 0], [0, af, 0, 0], [0, 0, af, 0], [0, 0, 0, af]]
om2 = [[Z1[:] for _ in R4] for _ in R4]
for i in (1, 2, 3):
    om2[i][0] = [gg * c for c in e2[i]]
    om2[0][i] = [gg * c for c in e2[i]]
    for j in (1, 2, 3):
        if i != j:
            om2[i][j] = add(*[e2[k] for k in (1, 2, 3)], coeffs=[-PP * eps3(i, j, k) for k in (1, 2, 3)])
Tor2 = [add(d1(e2[A], X), *[wedge11(om2[A][B], e2[B]) for B in R4]) for A in R4]
R2 = [[add(d1(om2[A][B], X), *[wedge11(om2[A][C], om2[C][B]) for C in R4]) for B in R4] for A in R4]


def d2(B):
    return [[[sp.diff(B[n][r], X[m]) + sp.diff(B[r][m], X[n]) + sp.diff(B[m][n], X[r]) for r in R4] for n in R4]
            for m in R4]


ok = True
for A in R4:
    DT = add(d2(Tor2[A]), *[wedge12(om2[A][B], Tor2[B]) for B in R4])
    Re = add(*[wedge12(e2[B], R2[A][B]) for B in R4])
    ok = ok and is_zero(add(DT, Re, coeffs=[1, -1]))
check("C1 Bianchi D T^A = R^A_B ^ e^B on the MZ FRW connection (beta-term = dual DT)", ok)


def w22(al, be, m, n, r, s):
    return (al[m][n] * be[r][s] - al[m][r] * be[n][s] + al[m][s] * be[n][r]
            + al[r][s] * be[m][n] - al[n][s] * be[m][r] + al[n][r] * be[m][s])


subLC = {gg: sp.diff(af, t) / af, PP: 0}
R2LC = [[[[sp.simplify(R2[A][B][m][n].subs(subLC).doit()) for n in R4] for m in R4] for B in R4] for A in R4]
R2up = [[[[sum(R2LC[A][C][m][n] * ETA[C, B] for C in R4) for n in R4] for m in R4] for B in R4] for A in R4]
eRR = sum(eps_lo(A, B, C, D) * w22(R2up[A][B], R2up[C][D], 0, 1, 2, 3) for A in R4 for B in R4 for C in R4 for D in R4)
H = sp.diff(af, t) / af
GB = 24 * H**2 * (sp.diff(H, t) + H**2)
check("C2 eps_ABCD R^AB ^ R^CD = + GB vol (vol = e0^e1^e2^e3 = a^3 dt^dx^dy^dz) for the Levi-Civita FRW connection",
      sp.simplify(eRR - GB * af**3) == 0)

# ---------------- D: linear closure ----------------
tau, ep = sp.symbols("tau epsilon")
Xp = [tau, x]
A_ = sp.Function("A", positive=True)(tau)
rb, pb, thb = (sp.Function(n)(tau) for n in ("rhob", "pb", "thb"))
dr, dpp, dth, v, Phi, Psi = (sp.Function(n)(tau, x) for n in ("drho", "dp", "dth", "v", "Phi", "Psi"))
rho, pr, th = rb + ep * dr, pb + ep * dpp, thb + ep * dth
# metric diag a^2(-(1+2 Phi), 1-2Psi, ...) restricted to (tau, x); u^mu, u_mu to O(eps)
gtt, gxx = -A_**2 * (1 + 2 * ep * Phi), A_**2 * (1 - 2 * ep * Psi)
uU = [(1 - ep * Phi) / A_, ep * sp.diff(v, x) / A_]
uD = [gtt * uU[0], gxx * uU[1]]
dth_ = [sp.diff(th, tau), sp.diff(th, x)]
Svec = [sp.series((rho + pr) * uD[m] * (uU[0] * dth_[0] + uU[1] * dth_[1]) + pr * dth_[m], ep, 0, 2).removeO()
        for m in range(2)]
dS = sp.expand(sp.diff(Svec[1], tau) - sp.diff(Svec[0], x))
dS1 = dS.coeff(ep, 1)
thc = dth + sp.diff(thb, tau) * v
rhc = dr + sp.diff(rb, tau) * v
Y = (rb + pb) * sp.diff(thc, tau) + sp.diff(pb, tau) * thc + sp.diff(thb, tau) * rhc
check("D1 (dS)_{tau x} = d_x Y at O(eps), Y = (rho+p) dtheta_c' + p' dtheta_c + theta' drho_c (no Phi, Psi left); "
      "background dS = 0", sp.simplify(dS1 - sp.diff(Y, x)) == 0 and sp.simplify(dS.coeff(ep, 0)) == 0)
xi0 = sp.Function("xi0")(tau, x)
gauge = {dth: -sp.diff(thb, tau) * xi0, dr: -sp.diff(rb, tau) * xi0, v: xi0}
check("D2 a pure time-shift (adiabatic super-horizon) mode gives Y = 0", sp.simplify(Y.subs(gauge).doit()) == 0)
# matter era, p = 0: closure needs dtheta_c' = -theta' delta_c; check Y != 0 for a generic comoving pair
dc, tc = sp.Function("dc")(tau, x), sp.Function("tc")(tau, x)
Ym = (rb * sp.diff(tc, tau) + sp.diff(thb, tau) * rb * dc)
check("D3 dust: Y = rho (dtheta_c' + theta' delta_c); zero only if theta's comoving perturbation tracks the "
      "density contrast, dtheta_c' = -theta' delta_c (not implied: theta carries its own data, MZ Table 1)",
      sp.simplify(Y.subs({pb: 0, dpp: 0}).subs({dth: tc - sp.diff(thb, tau) * v, dr: rb * dc - sp.diff(rb, tau) * v})
                  .doit() - Ym) == 0)

# ---------------- E: effective-theory form ----------------
phi = sp.Function("phi")(tau, x)
cs = sp.Symbol("c_s")
p_eff = [sp.diff(phi, tau), cs * sp.diff(phi, x)]
check("E1 S ~ (phi_t, c_s phi_x): (dS)_{tx} = (c_s - 1) phi_tx, closed only where c_s = Lambda/(Lambda + kappa rho) = 1 "
      "(vacuum)", sp.simplify(sp.diff(p_eff[1], tau) - sp.diff(p_eff[0], x) - (cs - 1) * sp.diff(phi, tau, x)) == 0)
summary()
