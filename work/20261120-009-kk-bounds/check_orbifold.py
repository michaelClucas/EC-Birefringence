"""Checks O1-O5: Chern-Simons densities on S1/Z2 (y -> -y), goal 4.

A term is kept on the orbifold only if its density L(x, y) (coefficient of
d^4x dy in a fixed orientation) is even in y; an odd density integrates to zero
over the circle and needs an explicit eps(y) coefficient. Fields with parity +
are cos-modes plus a y-independent zero mode, parity - are sin-modes.
"""
import itertools
import random

import sympy as sp

from common import check

random.seed(4)
X = sp.symbols("x0 x1 x2 x3", real=True)
Y = sp.Symbol("y", real=True)
C5 = list(X) + [Y]


def rp():
    mons = [sp.Integer(1)] + list(X) + [a * b for a, b in itertools.combinations(X, 2)]
    return sum(sp.Rational(random.randint(-7, 7), random.randint(1, 4)) * m for m in mons)


def field(par):
    if par > 0:
        return rp() + sp.cos(Y) * rp() + sp.cos(2 * Y) * rp()
    return sp.sin(Y) * rp() + sp.sin(2 * Y) * rp()


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


PERM5 = [(p, perm_sign(p)) for p in itertools.permutations(range(5))]


def cs_density(A1, A2, A3, pt):
    """eps~^{MNPQR} A1_M d_N A2_P d_Q A3_R at pt = (x0..x3, y)."""
    sub = dict(zip(C5, pt))
    v = [sp.nsimplify(c.subs(sub)) for c in A1]
    d2 = [[sp.nsimplify(sp.diff(c, C).subs(sub)) for c in A2] for C in C5]
    d3 = [[sp.nsimplify(sp.diff(c, C).subs(sub)) for c in A3] for C in C5]
    return sp.nsimplify(sum(s * v[p[0]] * d2[p[1]][p[2]] * d3[p[3]][p[4]] for p, s in PERM5))


def parity_of(A1, A2, A3):
    x = [sp.Rational(1, 3), sp.Rational(-1, 2), sp.Rational(2, 5), sp.Rational(3, 7)]
    yv = sp.pi / 5
    lp = cs_density(A1, A2, A3, x + [yv])
    lm = cs_density(A1, A2, A3, x + [-yv])
    if sp.simplify(lp - lm) == 0 and lp != 0:
        return +1
    if sp.simplify(lp + lm) == 0 and lp != 0:
        return -1
    return 0


def vec(mu_par, five_par):
    return [field(mu_par) for _ in range(4)] + [field(five_par)]


def main():
    A_odd = vec(-1, +1)   # A_mu odd, A_5 even (A_5 zero mode survives)
    A_even = vec(+1, -1)  # A_mu even (photon zero mode survives), A_5 odd
    check("O1 A^F^F with A_mu odd, A_5 even: density even (allowed; Wilson line kept, bulk photon zero mode lost)",
          parity_of(A_odd, A_odd, A_odd) == +1)
    check("O2 A^F^F with A_mu even, A_5 odd: density odd (needs eps(y); Wilson line projected out)",
          parity_of(A_even, A_even, A_even) == -1)
    B = vec(-1, +1)
    Yc = vec(+1, -1)
    check("O3 mixed A'^F_Y^F_Y (A'_mu odd, A'_5 even; A_Y mu even): density even (allowed)",
          parity_of(B, Yc, Yc) == +1)
    # zero-mode reduction of O3: only A'_5 = theta(x), A_Y mu = a_mu(x)
    th = rp()
    a = [rp() for _ in range(4)]
    Bz = [sp.Integer(0)] * 4 + [th]
    Yz = a + [sp.Integer(0)]
    x = [sp.Rational(1, 3), sp.Rational(-1, 2), sp.Rational(2, 5), sp.Rational(3, 7)]
    val = cs_density(Bz, Yz, Yz, x + [0])
    d = [[sp.diff(a[n], X[m]).subs(dict(zip(X, x))) for n in range(4)] for m in range(4)]
    ffd = sum(perm_sign(p) * d[p[0]][p[1]] * d[p[2]][p[3]] for p in itertools.permutations(range(4)))
    check("O4 O3 zero modes give theta' f_Y ftilde_Y (nonzero), density = +-theta eps da da",
          sp.simplify(val - th.subs(dict(zip(X, x))) * ffd) == 0 or sp.simplify(val + th.subs(dict(zip(X, x))) * ffd) == 0)
    # brane-localized theta(0) F Ftilde is not gauge invariant: Lambda odd => d_y Lambda(0) != 0
    lam = sp.sin(Y) * rp()
    check("O5 brane term A'_5(y=0) F Ftilde not gauge invariant (delta A'_5(0) = d_y Lambda(0) != 0 for odd Lambda)",
          sp.simplify(sp.diff(lam, Y).subs(Y, 0)) != 0 and lam.subs(Y, 0) == 0)


if __name__ == "__main__":
    main()
