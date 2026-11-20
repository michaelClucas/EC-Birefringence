"""Candidate (iii): Nieh-Yan density in our conventions.

 1. Identity check (first-order jets at a point, random metric with g(0)=eta and random
    x-dependent torsion): d_m([m n r s] T_{n r s}) = a [m n r s] T^l_{mn} T_{l r s} + b [m n r s] R_{m n r s}
    with constant a, b, where R_{m n r s} = g_{m l} R^l_{n r s}(Gamma) and
    R^r_{s m n} = d_m Gam^r_{n s} - d_n Gam^r_{m s} + Gam^r_{m l} Gam^l_{n s} - Gam^r_{n l} Gam^l_{m s}.
    The fitted (a, b) must be the same for independent random configurations.
 2. [m n r s] T_{n r s} = -6 sqrt(-g) S^m (random metric and torsion), so
    int theta * d_m([mnrs] T_nrs) d^4x = +6 int sqrt(-g) S^m d_m theta d^4x (after integration by parts).
"""
import itertools
import sys

import sympy as sp

from common import (ETA, R4, axial_vector, contortion_correct, lower1, perm, random_metric_point,
                    random_torsion, report, rng, zeros3)

X = sp.symbols("x0:4")
PERMS = [(p, perm(*p)) for p in itertools.permutations(range(4))]


def trunc(e, deg=1):
    e = sp.expand(e)
    if e == 0:
        return e
    P = sp.Poly(e, *X)
    return sum((c * sp.prod([X[i]**k for i, k in enumerate(mon)])
                for mon, c in P.terms() if sum(mon) <= deg), sp.Integer(0))


def at0(e):
    return e.subs({x: 0 for x in X})


def config(seed):
    r = rng(seed)
    # metric g = eta + h1 + h2, h1 linear, h2 quadratic, vanishing at x = 0
    g = sp.Matrix(ETA)
    for a in R4:
        for b in range(a, 4):
            v = sum(sp.Rational(r.randint(-5, 5), 10) * x for x in X)
            v += sum(sp.Rational(r.randint(-5, 5), 10) * X[i] * X[j] for i in R4 for j in range(i, 4))
            g[a, b] += v
            if a != b:
                g[b, a] += v
    T = zeros3()
    for l in R4:
        for m in R4:
            for n in range(m + 1, 4):
                v = sp.Rational(r.randint(-5, 5), 3) + sum(sp.Rational(r.randint(-5, 5), 4) * x for x in X)
                T[l][m][n], T[l][n][m] = v, -v
    return g, T


def nieh_yan_terms(g, T):
    h1 = (g - ETA).applyfunc(lambda e: trunc(e, 1))
    gi = (ETA - ETA * h1 * ETA).applyfunc(lambda e: trunc(e, 1))  # first-order inverse at g(0)=eta
    dg = [[[sp.diff(g[a, b], X[k]) for b in R4] for a in R4] for k in R4]
    Chr = [[[trunc(sp.Rational(1, 2) * sum(gi[l, rr] * (dg[m][rr][n] + dg[n][rr][m] - dg[rr][m][n])
                                            for rr in R4)) for n in R4] for m in R4] for l in R4]
    Tl = [[[trunc(sum(g[l, rr] * T[rr][m][n] for rr in R4)) for n in R4] for m in R4] for l in R4]
    Kl = contortion_correct(Tl)
    K = [[[trunc(sum(gi[l, rr] * Kl[rr][m][n] for rr in R4)) for n in R4] for m in R4] for l in R4]
    G = [[[Chr[l][m][n] + K[l][m][n] for n in R4] for m in R4] for l in R4]
    G0 = [[[at0(G[l][m][n]) for n in R4] for m in R4] for l in R4]
    dG = [[[[at0(sp.diff(G[l][m][n], X[k])) for n in R4] for m in R4] for l in R4] for k in R4]

    def Riem(r_, s, m, n):
        return (dG[m][r_][n][s] - dG[n][r_][m][s]
                + sum(G0[r_][m][l] * G0[l][n][s] - G0[r_][n][l] * G0[l][m][s] for l in R4))

    lhs = at0(sum(sgn * sp.diff(Tl[p[1]][p[2]][p[3]], X[p[0]]) for p, sgn in PERMS))
    T0 = [[[at0(T[l][m][n]) for n in R4] for m in R4] for l in R4]
    Tl0 = [[[at0(Tl[l][m][n]) for n in R4] for m in R4] for l in R4]
    A = sum(sgn * sum(T0[l][p[0]][p[1]] * Tl0[l][p[2]][p[3]] for l in R4) for p, sgn in PERMS)
    B = sum(sgn * sum(ETA[p[0], l] * Riem(l, p[1], p[2], p[3]) for l in R4) for p, sgn in PERMS)
    return sp.nsimplify(lhs), sp.nsimplify(A), sp.nsimplify(B)


def main():
    ok = True
    a, b = sp.symbols("a b")
    rows = [nieh_yan_terms(*config(s)) for s in (5, 6, 7)]
    sol = sp.solve([L - a * A_ - b * B_ for L, A_, B_ in rows[:2]], [a, b], dict=True)
    ok &= report("Nieh-Yan fit solvable", len(sol) == 1)
    L3, A3, B3 = rows[2]
    ok &= report("same (a, b) holds for an independent third configuration",
                 sp.simplify(L3 - sol[0][a] * A3 - sol[0][b] * B3) == 0)
    print(f"INFO: d_m([mnrs] T_nrs) = {sol[0][a]} [mnrs] T^l_mn T_lrs + {sol[0][b]} [mnrs] R_mnrs")
    r = rng(9)
    g, _ = random_metric_point(r)
    T = random_torsion(r)
    Tl = lower1(g, T)
    S = axial_vector(g, T)
    sq = sp.sqrt(-g.det())
    ok &= report("[mnrs] T_nrs = -6 sqrt(-g) S^m",
                 all(sp.simplify(sum(sgn * Tl[p[1]][p[2]][p[3]] for p, sgn in PERMS if p[0] == m)
                                 + 6 * sq * S[m]) == 0 for m in R4))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
