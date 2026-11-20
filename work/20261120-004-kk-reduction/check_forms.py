"""Check 4: the only epsilon-free (parity-odd) 5-form built from e, T, R is exact.

Identity in any dimension (Bianchi DT^A = R^A_B ^ e^B):
    d(T_A ^ T^A) = 2 R_AB ^ e^B ^ T^A .
In 5D, e_A^R^AB^T_B is the unique Lorentz-scalar 5-form built from (e, T, R) without eps, so a
"5D Holst / Nieh-Yan" term is a total derivative. Verified with an explicit exterior algebra on
random polynomial vielbein e^A_M(x) and an independent (non-Levi-Civita) connection w^AB_M(x),
x in R^5, evaluated at random points.
"""
import itertools
import numpy as np
import sympy as sp

rng = np.random.default_rng(11)
N = 5
X = sp.symbols("x0:5")
ETA = np.diag([-1, 1, 1, 1, 1])


def wedge(f, g):
    out = {}
    for I, a in f.items():
        for J, b in g.items():
            if set(I) & set(J):
                continue
            idx = I + J
            perm = sorted(range(len(idx)), key=lambda k: idx[k])
            # sign of sorting permutation
            p = list(perm); s = 1
            for i in range(len(p)):
                while p[i] != i:
                    j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
            key = tuple(sorted(idx))
            out[key] = out.get(key, 0) + s * a * b
    return out


def add(*fs, coef=None):
    out = {}
    coef = coef or [1] * len(fs)
    for c, f in zip(coef, fs):
        for k, v in f.items():
            out[k] = out.get(k, 0) + c * v
    return out


def d(f):
    out = {}
    for I, a in f.items():
        for i in range(N):
            if i in I:
                continue
            key = tuple(sorted((i,) + I))
            s = (-1) ** sum(1 for j in I if j < i)
            out[key] = out.get(key, 0) + s * sp.diff(a, X[i])
    return out


def rpoly():
    c = rng.integers(-3, 4, size=1 + N)
    return int(c[0]) + sum(int(ci) * x for ci, x in zip(c[1:], X))


def run():
    e = [{(M,): (1 if A == M else 0) + rpoly() for M in range(N)} for A in range(N)]
    w = [[None] * N for _ in range(N)]          # w^A_B with w_AB antisymmetric
    for A in range(N):
        w[A][A] = {}
        for B in range(A + 1, N):
            comp = {(M,): rpoly() for M in range(N)}
            # w^A_B = eta^AA w_AB ; w_BA = -w_AB
            w[A][B] = {k: ETA[A, A] * v for k, v in comp.items()}
            w[B][A] = {k: -ETA[B, B] * v for k, v in comp.items()}
    T = [add(d(e[A]), *[wedge(w[A][B], e[B]) for B in range(N)]) for A in range(N)]
    R = [[add(d(w[A][B]), *[wedge(w[A][C], w[C][B]) for C in range(N)]) for B in range(N)] for A in range(N)]
    TT = add(*[wedge({k: ETA[A, A] * v for k, v in T[A].items()}, T[A]) for A in range(N)])
    lhs = d(TT)
    rhs = add(*[wedge(wedge({k: ETA[A, A] * v for k, v in R[A][B].items()}, e[B]), T[A])
                for A in range(N) for B in range(N)], coef=[2] * (N * N))
    diff = add(lhs, rhs, coef=[1, -1])
    worst = 0
    for _ in range(2):
        pt = {x: sp.Rational(int(rng.integers(-5, 6)), 3) for x in X}
        for v in diff.values():
            worst = max(worst, abs(sp.sympify(v).xreplace(pt)))
    scale = max(abs(sp.sympify(v).xreplace({x: sp.Integer(1) for x in X})) for v in lhs.values()) if lhs else 0
    return [("5D: d(T_A^T^A) = 2 R_AB^e^B^T^A, so e^R^T (unique eps-free 5-form) is exact",
             worst == 0 and scale != 0, f"max |difference| at random rational points = {worst}, |lhs| sample {scale}")]


if __name__ == "__main__":
    for name, ok, det in run():
        print(("PASS" if ok else "FAIL"), name, "|", det)
