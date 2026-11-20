"""EC-04: the PGT 0- mode as a pseudoscalar (derivation Sec. C). Conventions of CONVENTIONS.md.

1. Pure axial torsion on flat space, T^l_{mn} = eta^{lr} eps_{rmnk} S^k with S^k(x) = s^k + M^k_a x^a:
   the pseudoscalar curvature X = eps^{mnrs} R_{mnrs}(Gamma) equals c * d_k S^k with no S^2 term.
   (So a PGT "E^2" term, Chen et al. 0908.3323 Eq. (7), is (d.S)^2 up to normalisation for axial torsion.)
2. Auxiliary-field reduction of L = (A/2) S.S + (w/2)(d.S)^2 + c_A S.J:
   introduce Phi: (w/2)X^2 = Phi X - Phi^2/(2w); eliminating S gives
   L = -(1/2A)(dPhi - c_A J)^2 - Phi^2/(2w): one scalar, healthy iff A > 0, non-tachyonic iff w > 0, m^2 = A/w.
3. w = 0 reproduces minimal EC: S = -(kappa/2) J for A = 3/(2 kappa), c_A = 3/4 (003: R ⊃ (3/2) S.S).
4. d_m (A_n Ftilde^{mn}) = (1/2) F_{mn} Ftilde^{mn}: the CFJ term with S = d theta is (xi/4) theta F Ftilde.
5. The 0- mode seen in momentum space: the S propagator (A eta + w k k)^-1 has a single pole at k^2 = -A/w.
"""
import itertools
import random
import sys

import sympy as sp

from common import report

R4 = range(4)
ETA = sp.diag(-1, 1, 1, 1)


def eps_l(a, b, c, d):  # eps_{abcd}, flat, [0123] = +1
    return sp.LeviCivita(a, b, c, d)


def eps_u(a, b, c, d):  # eps^{abcd} = -[abcd]
    return -sp.LeviCivita(a, b, c, d)


def main():
    ok = True
    x = sp.symbols("x0:4")
    s = sp.symbols("s0:4")
    M = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"M{i}{j}"))
    S = [s[k] + sum(M[k, a] * x[a] for a in R4) for k in R4]
    # T_{lmn} (all lower) = eps_{lmnk} S^k ; T^l_{mn} = eta^{ll} T_{lmn}
    Tl = [[[sum(eps_l(l, m, n, k) * S[k] for k in R4) for n in R4] for m in R4] for l in R4]
    # axial vector back: S^m = (1/6) eps^{mnrs} T_{nrs}
    Sback = [sp.Rational(1, 6) * sum(eps_u(m, n, r, q) * Tl[n][r][q] for n, r, q in itertools.product(R4, R4, R4)) for m in R4]
    ok &= report("axial torsion ansatz returns S^m = (1/6) eps^{mnrs} T_nrs", all(sp.expand(Sback[m] - S[m]) == 0 for m in R4))
    K = [[[sp.Rational(1, 2) * (Tl[l][m][n] + Tl[n][l][m] - Tl[m][n][l]) for n in R4] for m in R4] for l in R4]
    G = [[[ETA[l, l] * K[l][m][n] for n in R4] for m in R4] for l in R4]  # Gamma^l_{mn}, flat metric

    def Riem(l, sg, m, n):  # R^l_{sigma m n}
        v = sp.diff(G[l][n][sg], x[m]) - sp.diff(G[l][m][sg], x[n])
        v += sum(G[l][m][r] * G[r][n][sg] - G[l][n][r] * G[r][m][sg] for r in R4)
        return v

    X = 0
    for l, sg, m, n in itertools.product(R4, R4, R4, R4):
        e = eps_u(l, sg, m, n)
        if e != 0:
            X += e * ETA[l, l] * Riem(l, sg, m, n)
    X = sp.expand(X.subs({xi: 0 for xi in x}))
    divS = sum(M[k, k] for k in R4)
    c = sp.simplify(X / divS) if sp.expand(X - (X.coeff(M[0, 0])) * divS) == 0 else None
    print(f"X = eps^{{mnrs}} R_mnrs = {X}")
    ok &= report("X = c d.S with no S^2 term", c is not None and all(sp.diff(X, si) == 0 for si in s))
    print(f"   c = {c}")

    # 2. auxiliary-field reduction (component symbols, signature -+++)
    A, w, cA = sp.symbols("A w c_A", nonzero=True)
    Sv = sp.symbols("S0:4")
    dP = sp.symbols("dPhi0:4")
    J = sp.symbols("J0:4")
    Phi = sp.Symbol("Phi")
    dot = lambda u, v: sum(ETA[i, i] * u[i] * v[i] for i in R4)
    Laux = A / 2 * dot(Sv, Sv) - dot(dP, Sv) - Phi**2 / (2 * w) + cA * dot(J, Sv)  # Phi d.S -> -dPhi.S by parts
    sol = sp.solve([sp.diff(Laux, si) for si in Sv], Sv, dict=True)[0]
    Lred = sp.expand(Laux.subs(sol))
    target = sp.expand(-dot([dP[i] - cA * J[i] for i in R4], [dP[i] - cA * J[i] for i in R4]) / (2 * A) - Phi**2 / (2 * w))
    ok &= report("eliminating S: L = -(1/2A)(dPhi - c_A J)^2 - Phi^2/(2w)", sp.simplify(Lred - target) == 0)
    ok &= report("S = (dPhi - c_A J)/A on shell", all(sp.simplify(sol[Sv[i]] - (dP[i] - cA * J[i]) / A) == 0 for i in R4))
    Xs = sp.Symbol("X")
    Pst = sp.solve(sp.diff(Xs * Phi - Phi**2 / (2 * w), Phi), Phi)[0]
    ok &= report("(w/2) X^2 = max_Phi [Phi X - Phi^2/(2w)] (Phi = w X)", sp.simplify((Xs * Phi - Phi**2 / (2 * w)).subs(Phi, Pst) - w * Xs**2 / 2) == 0)
    print("   canonical phi = Phi/sqrt(A): kinetic -(1/2)(d phi)^2 needs A > 0; mass^2 = A/w needs w > 0")

    # 3. minimal-EC limit
    kap = sp.Symbol("kappa", positive=True)
    LEC = sp.Rational(3, 4) / kap * dot(Sv, Sv) + sp.Rational(3, 4) * dot(J, Sv)
    solEC = sp.solve([sp.diff(LEC, si) for si in Sv], Sv, dict=True)[0]
    ok &= report("w = 0, A = 3/(2 kappa), c_A = 3/4 gives S = -(kappa/2) J (003)",
                 all(sp.simplify(solEC[Sv[i]] + kap / 2 * J[i]) == 0 for i in R4))

    # 4. d_m K^m = (1/2) F Ftilde with K^m = A_n Ftilde^{mn}
    Af = [sp.Function(f"A{i}")(*x) for i in R4]
    Al = Af  # treat as A_mu (lower); flat
    F = [[sp.diff(Al[n], x[m]) - sp.diff(Al[m], x[n]) for n in R4] for m in R4]
    Fu = [[ETA[m, m] * ETA[n, n] * F[m][n] for n in R4] for m in R4]
    Ft = [[sp.Rational(1, 2) * sum(eps_u(m, n, r, q) * F[r][q] for r in R4 for q in R4) for n in R4] for m in R4]
    Kv = [sum(Al[n] * Ft[m][n] for n in R4) for m in R4]
    lhs = sp.expand(sum(sp.diff(Kv[m], x[m]) for m in R4))
    rhs = sp.expand(sp.Rational(1, 2) * sum(F[m][n] * Ft[m][n] for m in R4 for n in R4))
    ok &= report("d_m (A_n Ftilde^{mn}) = (1/2) F_mn Ftilde^mn", sp.simplify(lhs - rhs) == 0)
    print("   => -(1/2) xi (d_m theta) K^m = (xi/4) theta F Ftilde + total derivative: g_theta = -xi (L = -(g/4) theta F Ftilde)")

    # 5. propagator pole
    k = sp.symbols("k0:4")
    Ap, wp = sp.symbols("A_p w_p", positive=True)
    kl = [ETA[i, i] * k[i] for i in R4]  # k_mu
    Mmat = sp.Matrix(4, 4, lambda i, j: Ap * ETA[i, j] + wp * kl[i] * kl[j])
    det = sp.factor(Mmat.det())
    k2 = sum(ETA[i, i] * k[i] ** 2 for i in R4)
    ok &= report("det(A eta + w k k) = -A^3 (A + w k^2): single pole k^2 = -A/w", sp.simplify(det + Ap**3 * (Ap + wp * k2)) == 0)
    return ok


if __name__ == "__main__":
    random.seed(1)
    sys.exit(0 if main() else 1)
