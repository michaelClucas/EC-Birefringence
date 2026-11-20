"""Skeptic 011: independent check of the P' = diag(1,-1,-1,-1,-1) argument (task 009, KK-01).

1. P' is in the identity component of SO(4,1): it is exp of a boost-free generator, det = +1,
   preserves time orientation, and has a spinor lift S with S^-1 Gamma^M S = P'^M_N Gamma^N.
2. Pull back the KK metric g_{mu5} = sqrt2 phi^2 A_mu, g_55 = phi^2 under (x^i, y) -> (-x^i, -y)
   with generic functions; read off phi', A'; show phi'^p E'.B' = -phi^p E.B at the image point
   (odd) and phi'^p (E'^2 - B'^2) = + (even). A bulk-U(1) Wilson line theta = a_5 is odd, so
   theta f f~ is even (allowed): the check distinguishes the two.
3. Loophole probe: one-loop KK tower of a 5D Dirac fermion with mass M and a Wilson-line phase w.
   Schematic graviphoton FF~ coefficient ~ sum_n k^2 arg (graviphoton charge k = n + w)(M + i (n + w)/R). Zero for w = 0 (P'-even
   vacuum), nonzero for w != 0, pi (P'-odd background) -- consistent with the theorem.
"""
import itertools

import numpy as np
import sympy as sp
from scipy.linalg import expm


def report(name, ok, det=""):
    print(("PASS" if ok else "FAIL") + ": " + name + (" -- " + det if det else ""))


# ---- 1 ------------------------------------------------------------------------------
def J(i, j, n=5):
    m = np.zeros((n, n))
    m[i, j], m[j, i] = -1.0, 1.0          # rotation generator in the (i, j) spatial plane
    return m


P = expm(np.pi * (J(1, 2) + J(3, 4)))
Pexp = np.diag([1, -1, -1, -1, -1.0])
report("exp(pi (J12 + J3y)) = diag(1,-1,-1,-1,-1)", np.allclose(P, Pexp, atol=1e-12))
eta = np.diag([-1, 1, 1, 1, 1.0])
report("P' in SO(4,1), det +1, time orientation preserved",
       np.allclose(P.T @ eta @ P, eta) and abs(np.linalg.det(P) - 1) < 1e-12 and P[0, 0] > 0)

# 5D Clifford algebra (signature -++++): Gamma^0 = i gamma^0_E ... build from Pauli matrices
s0 = np.eye(2); sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1])
E = [np.kron(sx, s) for s in (sx, sy, sz)] + [np.kron(sy, s0), np.kron(sz, s0)]   # Euclidean 5D
G = [1j * E[0]] + E[1:]                                   # Gamma^0 squares to -1
okc = all(np.allclose(G[a] @ G[b] + G[b] @ G[a], 2 * eta[a, b] * np.eye(4)) for a in range(5) for b in range(5))
S = expm(np.pi / 2 * (G[1] @ G[2] + G[3] @ G[4]))        # exp((theta/2) G_a G_b), theta = pi per plane
Si = np.linalg.inv(S)
oks = all(np.allclose(Si @ G[m] @ S, sum(P[m, n] * G[n] for n in range(5))) for m in range(5))
report("Clifford algebra built", okc)
report("spinor lift S exists with S^-1 Gamma^M S = P'^M_N Gamma^N (5D fermions incl. masses covariant)",
       oks)
report("S^2 = +1 on spinors (two rotations by pi)", np.allclose(S @ S, np.eye(4)))

# ---- 2 ------------------------------------------------------------------------------
t, x1, x2, x3, y = X = sp.symbols("t x1 x2 x3 y", real=True)
phi = sp.Function("phi")(t, x1, x2, x3)
A = [sp.Function(f"A{i}")(t, x1, x2, x3) for i in range(4)]
g4 = sp.diag(-1, 1, 1, 1)
G5 = sp.zeros(5, 5)
for m in range(4):
    for n in range(4):
        G5[m, n] = g4[m, n] + 2 * phi**2 * A[m] * A[n]
    G5[m, 4] = G5[4, m] = sp.sqrt(2) * phi**2 * A[m]
G5[4, 4] = phi**2
# diffeo x' = P x (P is its own inverse): g'_{MN}(x') = P^A_M P^B_N g_AB(P x')
Pm = sp.diag(1, -1, -1, -1, -1)
sub = {t: t, x1: -x1, x2: -x2, x3: -x3}
Gp = (Pm.T * G5 * Pm).subs(sub, simultaneous=True)
report("g'_55(x) = phi(Px)^2", sp.simplify(Gp[4, 4] - phi.subs(sub, simultaneous=True)**2) == 0)
phip = phi.subs(sub, simultaneous=True)   # phi > 0
Ap = [sp.simplify(Gp[m, 4] / (sp.sqrt(2) * Gp[4, 4])) for m in range(4)]
report("A'_0(x) = -A_0(Px), A'_i(x) = +A_i(Px)",
       sp.simplify(Ap[0] + A[0].subs(sub, simultaneous=True)) == 0 and
       all(sp.simplify(Ap[i] - A[i].subs(sub, simultaneous=True)) == 0 for i in (1, 2, 3)))

xs = [t, x1, x2, x3]


def F(Avec):
    return sp.Matrix(4, 4, lambda m, n: sp.diff(Avec[n], xs[m]) - sp.diff(Avec[m], xs[n]))


def invariants(Avec):
    Fm = F(Avec)
    Ev = [Fm[0, i] for i in (1, 2, 3)]
    Bv = [Fm[2, 3], Fm[3, 1], Fm[1, 2]]
    EB = sum(e * b for e, b in zip(Ev, Bv))
    E2B2 = sum(e**2 for e in Ev) - sum(b**2 for b in Bv)
    return EB, E2B2


EB, E2B2 = invariants(A)
EBp, E2B2p = invariants(Ap)
p = sp.Symbol("p")
odd = sp.simplify(phip**p * EBp + (phi**p * EB).subs(sub, simultaneous=True))
even = sp.simplify(phip**p * E2B2p - (phi**p * E2B2).subs(sub, simultaneous=True))
report("phi^p F_KK F~_KK is P'-odd (forbidden in a P'-invariant vacuum)", odd == 0)
report("phi^p F_KK^2 is P'-even", even == 0)
# negative control: the same transformation on an ordinary bulk-U(1) zero mode a_mu (no g_{mu5} sign)
a = [sp.Function(f"a{i}")(t, x1, x2, x3) for i in range(4)]
a5 = sp.Function("a5")(t, x1, x2, x3)
ap = [a[0].subs(sub, simultaneous=True)] + [-a[i].subs(sub, simultaneous=True) for i in (1, 2, 3)]
a5p = -a5.subs(sub, simultaneous=True)
eba, _ = invariants(a)
ebap, _ = invariants(ap)
report("bulk U(1): theta f f~ (theta = a_5) is P'-even (KK-02 allowed, check discriminates)",
       sp.simplify(a5p * ebap - (a5 * eba).subs(sub, simultaneous=True)) == 0)

# ---- 3 ------------------------------------------------------------------------------
def tower(w, MR=0.7, N=4000):
    n = np.arange(-N, N + 1)
    k = n + w
    return np.sum(k**2 * np.arctan2(k, MR) * np.exp(-(k / 400.0) ** 2))   # smooth regulator, even


c0 = tower(0.0)
c1 = tower(0.13)
report("KK fermion tower: graviphoton FF~ coefficient = 0 with no Wilson line", abs(c0) < 1e-6 * abs(c1),
       f"w=0: {c0:.2e}, w=0.13: {c1:.2e}")
report("nonzero with a P'-odd Wilson line (w = 0.13)", abs(c1) > 1e-4)
report("w = 1/2 (theta = pi, P'-even point) gives zero again", abs(tower(0.5)) < 1e-6 * abs(c1),
       f"{tower(0.5):.2e}")
