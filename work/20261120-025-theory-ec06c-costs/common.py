"""Shared helpers for task 20261120-025 (copied from task 022, plus perturbative tools at the end).

Gam[l][m][n] = Gamma^l_{mn}, derivative index m first; T^l_{mn} = Gam^l_{mn} - Gam^l_{nm};
trace T_m = T^l_{ml}; axial S^m = (1/6) eps^{mnrs} T_{nrs}, eps_{0123} = +sqrt(-g) so eps^{0123} = -1/sqrt(-g).
"""
import csv
import os

import sympy as sp

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
R4 = range(4)
_count = {"PASS": 0, "FAIL": 0}


def check(name, ok, msg=""):
    tag = "PASS" if ok else "FAIL"
    _count[tag] += 1
    print(f"{tag}: {name} {msg}".rstrip())


def summary():
    print(f"== {_count['PASS']} PASS, {_count['FAIL']} FAIL ==")


def load(name):
    with open(os.path.join(ROOT, "data", name), newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def weitzenbock_torsion(e, X):
    """e[A, mu] = e^A_mu (frame index row). Weitzenboeck Gamma^l_{mn} = e_A^l d_m e^A_n (omega = 0)."""
    einv = e.inv()  # einv[mu, A] = e_A^mu
    G = [[[sum(einv[l, A] * sp.diff(e[A, n], X[m]) for A in R4) for n in R4] for m in R4] for l in R4]
    return [[[sp.simplify(G[l][m][n] - G[l][n][m]) for n in R4] for m in R4] for l in R4]


def trace_vec(T):
    return [sp.simplify(sum(T[l][m][l] for l in R4)) for m in R4]


def axial_vec(T, g):
    """S^mu = (1/6) eps^{mu nu rho sigma} T_{nu rho sigma}, eps^{0123} = -1/sqrt(-g)."""
    sq = sp.sqrt(-g.det())
    Tl = [[[sum(g[l, k] * T[k][m][n] for k in R4) for n in R4] for m in R4] for l in R4]
    out = []
    for mu in R4:
        s = 0
        for n in R4:
            for r in R4:
                for q in R4:
                    lc = sp.LeviCivita(mu, n, r, q)
                    if lc:
                        s += -lc / sq * Tl[n][r][q]
        out.append(sp.simplify(s / 6))
    return out


def is_closed(p, X):
    return all(sp.simplify(sp.diff(p[m], X[n]) - sp.diff(p[n], X[m])) == 0 for m in R4 for n in R4)


# ---------------------------------------------------------------- perturbative tools (task 025)
EPS = sp.Symbol("epsilon")


def trunc(expr, order=2):
    """Keep powers of EPS up to `order`."""
    ex = sp.expand(expr)
    return sum(ex.coeff(EPS, k) * EPS**k for k in range(order + 1))


def series_inv(e, order=2):
    """Inverse of a matrix polynomial in EPS, truncated at EPS^order."""
    e0 = e.applyfunc(lambda v: sp.expand(v).coeff(EPS, 0))
    i0 = e0.inv()
    D = (i0 * (e - e0)).applyfunc(lambda v: trunc(v, order))
    out = sp.eye(4)
    term = sp.eye(4)
    for _ in range(order):
        term = (-term * D).applyfunc(lambda v: trunc(v, order))
        out += term
    return (out * i0).applyfunc(lambda v: trunc(v, order))


def torsion_series(e, X, order=2):
    """Weitzenboeck torsion T^l_{mn} of tetrad e[A, mu], truncated at EPS^order."""
    einv = series_inv(e, order)
    G = [[[trunc(sum(einv[l, A] * sp.diff(e[A, n], X[m]) for A in R4), order) for n in R4] for m in R4]
         for l in R4]
    return [[[G[l][m][n] - G[l][n][m] for n in R4] for m in R4] for l in R4]


def torsion_scalar(T, g, gi, order=2):
    """T_sc = 1/4 T^r_{mn} T_r^{mn} + 1/2 T^r_{mn} T^{nm}_r - v_m v^m with v_m = T^l_{ml}.
    Equals +6H^2 on flat FRW with the diagonal tetrad. Returns (T_sc, v)."""
    v = [sum(T[l][m][l] for l in R4) for m in R4]
    s = 0
    for r in R4:
        for m in R4:
            for n in R4:
                if T[r][m][n] == 0:
                    continue
                low = sum(g[r, q] * gi[m, b] * gi[n, c] * T[q][b][c] for q in R4 for b in R4 for c in R4)
                mix = sum(gi[m, b] * T[n][b][r] for b in R4)  # T^{n m}_r := g^{mb} T^n_{br}
                s += sp.Rational(1, 4) * T[r][m][n] * low + sp.Rational(1, 2) * T[r][m][n] * mix
    s -= sum(v[m] * gi[m, n] * v[n] for m in R4 for n in R4)
    return trunc(s, order), [trunc(c, order) for c in v]
