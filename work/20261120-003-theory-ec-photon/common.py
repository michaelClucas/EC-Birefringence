"""Shared helpers for the symbolic checks (conventions of CONVENTIONS.md).

Index conventions used everywhere:
  Gam[l][m][n] = Gamma^l_{m n}, derivative index m (nabla_m V^l = d_m V^l + Gamma^l_{m n} V^n)
  T[l][m][n]   = T^l_{m n} = Gamma^l_{m n} - Gamma^l_{n m}
  eps_lower(g)[a][b][c][d] = epsilon_{abcd} = +sqrt(-g) [abcd],  [0123] = +1
"""
import itertools
import random

import sympy as sp

N = 4
R4 = range(N)
ETA = sp.diag(-1, 1, 1, 1)


def perm(*idx):
    return sp.LeviCivita(*idx)


def rr(rng, lo=-9, hi=9):
    return sp.Rational(rng.randint(lo, hi), rng.randint(1, 7))


def zeros3():
    return [[[sp.Integer(0)] * N for _ in R4] for _ in R4]


def random_metric_point(rng):
    """Random Lorentzian metric value g and first derivatives dg[k][a][b] at a point."""
    while True:
        h = sp.Matrix(4, 4, lambda i, j: 0)
        for i in R4:
            for j in range(i, N):
                v = rr(rng) / 10
                h[i, j] = h[j, i] = v
        g = ETA + h
        if g.det() < 0:
            break
    dg = [[[None] * N for _ in R4] for _ in R4]
    for k in R4:
        for a in R4:
            for b in range(a, N):
                v = rr(rng)
                dg[k][a][b] = dg[k][b][a] = v
    return g, dg


def random_torsion(rng):
    T = zeros3()
    for l in R4:
        for m in R4:
            for n in range(m + 1, N):
                v = rr(rng)
                T[l][m][n] = v
                T[l][n][m] = -v
    return T


def christoffel(g, dg):
    gi = g.inv()
    C = zeros3()
    for l, m, n in itertools.product(R4, R4, R4):
        C[l][m][n] = sp.Rational(1, 2) * sum(
            gi[l, r] * (dg[m][r][n] + dg[n][r][m] - dg[r][m][n]) for r in R4)
    return C


def lower1(g, X):
    """X^l_{mn} -> X_{lmn}."""
    return [[[sum(g[l, r] * X[r][m][n] for r in R4) for n in R4] for m in R4] for l in R4]


def raise1(g, X):
    gi = g.inv()
    return [[[sum(gi[l, r] * X[r][m][n] for r in R4) for n in R4] for m in R4] for l in R4]


def contortion_correct(Tl):
    """K_{lmn} = 1/2 (T_{lmn} + T_{nlm} - T_{mnl}) (derived in derivation_ec_photon.md, Sec. 1)."""
    return [[[sp.Rational(1, 2) * (Tl[l][m][n] + Tl[n][l][m] - Tl[m][n][l])
              for n in R4] for m in R4] for l in R4]


def contortion_conventions_md(Tl):
    """Formula as printed in CONVENTIONS.md iter 1: 1/2 (T_{lmn} - T_{mln} - T_{nlm})."""
    return [[[sp.Rational(1, 2) * (Tl[l][m][n] - Tl[m][l][n] - Tl[n][l][m])
              for n in R4] for m in R4] for l in R4]


def eps_lower(g):
    s = sp.sqrt(-g.det())
    return lambda a, b, c, d: s * perm(a, b, c, d)


def eps_upper(g):
    s = sp.sqrt(-g.det())
    return lambda a, b, c, d: -perm(a, b, c, d) / s


def axial_vector(g, T):
    """S^m = (1/6) eps^{m n r s} T_{n r s}."""
    Tl = lower1(g, T)
    eu = eps_upper(g)
    return [sp.Rational(1, 6) * sum(eu(m, n, r, s) * Tl[n][r][s]
                                    for n, r, s in itertools.product(R4, R4, R4)) for m in R4]


def trace_vector(T):
    """T_m = T^l_{m l}."""
    return [sum(T[l][m][l] for l in R4) for m in R4]


def torsion_from_parts(g, Tvec, Svec, q=None):
    """T^l_{mn} = 1/3 (delta^l_n T_m - delta^l_m T_n) + g^{lr} eps_{r m n k} S^k + q."""
    gi = g.inv()
    el = eps_lower(g)
    T = zeros3()
    for l, m, n in itertools.product(R4, R4, R4):
        v = sp.Rational(1, 3) * ((1 if l == n else 0) * Tvec[m] - (1 if l == m else 0) * Tvec[n])
        v += sum(gi[l, r] * el(r, m, n, k) * Svec[k] for r in R4 for k in R4)
        if q is not None:
            v += q[l][m][n]
        T[l][m][n] = sp.expand(v)
    return T


def report(name, ok):
    print(("PASS" if ok else "FAIL") + ": " + name)
    return bool(ok)


def rng(seed):
    return random.Random(seed)
