"""Checks P1-P8: the proper-rotation parity P' and the 2608.09299 term.

P' = diag(1,-1,-1,-1,-1) on (x0,x1,x2,x3,y). It is a rotation by pi in the
(x1,x2) and (x3,y) planes, so it lies in the identity component of SO(4,1) and
is a symmetry of any diffeomorphism-invariant 5D action on M4 x S1 (y -> -y maps
the circle to itself). Fields are random exact polynomials; densities are
compared at a point p and at its image P'p.
"""
import itertools
import random

import numpy as np
import sympy as sp
from scipy.linalg import expm

from common import check, info

random.seed(9)
X = sp.symbols("x0 x1 x2 x3", real=True)
Y = sp.Symbol("y", real=True)
C5 = list(X) + [Y]
SQ2 = sp.sqrt(2)


def rpoly(deg=2):
    mons = [sp.Integer(1)] + list(X) + [a * b for a, b in itertools.combinations_with_replacement(X, 2)]
    if deg >= 3:
        mons += [X[0] * X[1] * X[2], X[3] ** 3]
    return sum(sp.Rational(random.randint(-9, 9), random.randint(1, 5)) * m for m in mons)


def mirror(expr):
    """expr(P x): x_i -> -x_i (i = 1,2,3); fields are y-independent."""
    return expr.subs({X[1]: -X[1], X[2]: -X[2], X[3]: -X[3]}, simultaneous=True)


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


PERM5 = [(p, perm_sign(p)) for p in itertools.permutations(range(5))]
PERM4 = [(p, perm_sign(p)) for p in itertools.permutations(range(4))]


def _jet(fields, pt):
    """Values and first derivatives (5D, d_y = 0) of a list of fields at pt, exact."""
    sub = dict(zip(X, pt))
    val = [sp.nsimplify(f.subs(sub)) for f in fields]
    der = [[sp.nsimplify(sp.diff(f, C).subs(sub)) for f in fields] for C in C5]
    return val, der


def cs5(n, pt):
    """eps~^{MNPQR} n_M d_N n_P d_Q n_R at pt (symbol eps~^{0123y} = +1)."""
    v, d = _jet(n, pt)
    return sp.expand(sum(s * v[p[0]] * d[p[1]][p[2]] * d[p[3]][p[4]] for p, s in PERM5))


def ffd(A, pt):
    """eps~^{mnrs} d_m A_n d_r A_s = (1/4) eps~ F F at pt."""
    _, d = _jet(A, pt)
    return sp.expand(sum(s * d[p[0]][p[1]] * d[p[2]][p[3]] for p, s in PERM4))


def ff(A, pt):
    """F_mn F^mn with flat eta = diag(-1,1,1,1) at pt."""
    eta = [-1, 1, 1, 1]
    _, d = _jet(A, pt)
    F = [[d[m][n] - d[n][m] for n in range(4)] for m in range(4)]
    return sp.expand(sum(eta[m] * eta[n] * F[m][n] ** 2 for m in range(4) for n in range(4)))


def at(expr, pt):
    return sp.nsimplify(expr.subs(dict(zip(X, pt))))


def main():
    # P1: P' in identity component, det +1
    J = np.zeros((5, 5))
    J[1, 2], J[2, 1] = -1, 1
    J[3, 4], J[4, 3] = -1, 1
    R = expm(np.pi * J)
    check("P1 P' = exp(pi(J12 + J3y)) = diag(1,-1,-1,-1,-1), det +1",
          np.allclose(R, np.diag([1, -1, -1, -1, -1]), atol=1e-12) and abs(np.linalg.det(R) - 1) < 1e-12)

    phi = 3 + rpoly() / 10
    A = [rpoly(3) for _ in range(4)]
    pts = [(sp.Rational(1, 3), sp.Rational(-2, 7), sp.Rational(5, 11), sp.Rational(3, 4)),
           (sp.Rational(-4, 5), sp.Rational(1, 6), sp.Rational(2, 3), sp.Rational(-1, 9)),
           (sp.Rational(7, 8), sp.Rational(3, 5), sp.Rational(-5, 6), sp.Rational(2, 7))]
    pt = pts[0]
    mpt = (pt[0], -pt[1], -pt[2], -pt[3])

    # transformed KK fields: g'_{MN}(x) = L^A_M L^B_N g_AB(P'x), L = P'
    phiT = mirror(phi)
    AT = [-mirror(A[0])] + [mirror(A[i]) for i in (1, 2, 3)]
    n = [SQ2 * phi * A[m] for m in range(4)] + [phi]

    # P2: 2608 Eq. (41) EM term: n ^ dn ^ dn = 2 phi^3 (dx^5 ^ dA ^ dA); d phi terms cancel
    sgn = set()
    for q in pts:
        S1, S2, ph = cs5(n, q), ffd(A, q), at(phi, q)
        if sp.simplify(S1 - 2 * ph**3 * S2) == 0:
            sgn.add("+")
        elif sp.simplify(S1 + 2 * ph**3 * S2) == 0:
            sgn.add("-")
        else:
            sgn.add("x")
    check("P2 eps n dn dn = 2 phi^3 eps dA dA at 3 random points (2608 Eq. 41 EM term)",
          len(sgn) == 1 and "x" not in sgn, "sign " + "".join(sgn))

    # P3: teleparallel A = a part, e^5 ^ T^a ^ T_a: Weitzenboeck torsion of e^a = (e^a_mu(x), 0)
    eta4 = [-1, 1, 1, 1]
    tet = [[(1 if a == m else 0) + rpoly() / 20 for m in range(4)] for a in range(4)]
    sub = dict(zip(X, pt))
    Ta = [[[sp.nsimplify((sp.diff(tet[a][nn], C5[mm]) - sp.diff(tet[a][mm], C5[nn])).subs(sub))
            if (mm < 4 and nn < 4) else 0 for nn in range(5)] for mm in range(5)] for a in range(4)]

    def nTT(nvec):
        v = [sp.nsimplify(c.subs(sub)) for c in nvec]
        return sp.expand(sum(s * v[p[0]] * sum(eta4[a] * Ta[a][p[1]][p[2]] * Ta[a][p[3]][p[4]] for a in range(4))
                             for p, s in PERM5))
    n2 = [SQ2 * phi * (A[m] + rpoly()) for m in range(4)] + [phi]
    v1, v2 = nTT(n), nTT(n2)
    check("P3 e^5 ^ T^a ^ T_a independent of A_mu (the Nieh-Yan part carries no photon)",
          sp.simplify(v1 - v2) == 0 and v1 != 0)

    # P4/P5: radion x F_KK Ftilde_KK odd under P'; Maxwell phi^3 F^2 even
    D_m = at(phi, mpt) ** 3 * ffd(A, mpt)
    D_T = at(phiT, pt) ** 3 * ffd(AT, pt)
    K_m = at(phi, mpt) ** 3 * ff(A, mpt)
    K_T = at(phiT, pt) ** 3 * ff(AT, pt)
    check("P4 phi^3 F_KK Ftilde_KK odd under P' (forbidden in a P'-invariant 5D theory, all orders)",
          sp.simplify(D_T + D_m) == 0 and D_m != 0)
    check("P5 phi^3 F_KK^2 even under P'", sp.simplify(K_T - K_m) == 0)

    # P6: KK-02 Wilson line: theta = Acal_5 (odd under P'), bulk photon a_mu (polar vector)
    th = rpoly()
    a = [rpoly(3) for _ in range(4)]
    thT = -mirror(th)
    aT = [mirror(a[0])] + [-mirror(a[i]) for i in (1, 2, 3)]
    check("P6 theta f ftilde (KK-02) even under P' (allowed)",
          sp.simplify(at(thT, pt) * ffd(aT, pt) - at(th, mpt) * ffd(a, mpt)) == 0)

    # P7: the 2608 functional with n fixed along +d_y flips sign
    nT = [SQ2 * phiT * AT[m] for m in range(4)] + [phiT]
    check("P7 n^dn^dn with n oriented along +d_y is odd under P' (needs a P'-odd background)",
          sp.simplify(cs5(nT, pt) + cs5(n, mpt)) == 0 and cs5(n, mpt) != 0)

    # P8: covariant pullback P'^* n = -n[phi',A']
    Lam = [1, -1, -1, -1, -1]
    npull = [Lam[M] * mirror(n[M]) for M in range(5)]
    check("P8 pullback P'^* n = -n[phi',A']: the sign of u relative to d_y is the P'-odd datum",
          all(sp.expand(npull[M] + nT[M]) == 0 for M in range(5)))
    info("P'-odd backgrounds that allow a P'-even scalar x F_KK Ftilde_KK: an oriented unit vector u "
         "(2608 Stueckelberg vacuum with det(d psi) > 0), FRW axial torsion S_0, a fermion chiral asymmetry.")


if __name__ == "__main__":
    main()
