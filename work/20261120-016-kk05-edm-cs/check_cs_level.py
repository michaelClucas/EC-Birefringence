"""C1-C6: quantisation of the 5D U(1) Chern-Simons level (KK-02) and beta_max.

S_CS = kappa_CS int A^F^F with kappa_CS = k/(24 pi^2), charges normalised so that int F/2pi = c_1 in Z.
exp(i S) must not depend on the 6D extension: S_6 = kappa_CS int_6 F^F^F must lie in 2 pi Z on closed spin 6-manifolds.
"""
import sympy as sp

from common import ALPHA, BETA_ACT, DEG, check, info

h, n, k, pi = sp.symbols("h n k pi", positive=True)
c1 = sp.Symbol("c1")

# C1: 6D extension of the CS term in units of 2 pi
F = 2 * sp.pi * c1
S6 = sp.Rational(1, 24) / sp.pi ** 2 * F ** 3
check("C1 (1/(24 pi^2)) int F^3 = 2 pi * (1/6) int c_1^3", sp.simplify(S6 - 2 * sp.pi * c1 ** 3 / 6) == 0)


def top(expr):
    """coefficient of h^3 (int_{CP^3} h^3 = 1)."""
    return sp.expand(sp.series(expr, h, 0, 4).removeO()).coeff(h, 3)


# C2: CP^3 is spin (c_1(T) = 4h even). Dirac index twisted by O(n): A-hat route vs Hirzebruch-Riemann-Roch route
p1 = 4 * h ** 2                                       # p_1(CP^3) = (1 + h^2)^4 -> 4 h^2
idx_ahat = top(sp.exp(n * h) * (1 - p1 / 24))
todd = (h / (1 - sp.exp(-h))) ** 4
idx_hrr = top(sp.exp((n - 2) * h) * todd)             # spinor bundle = K^{1/2} x Lambda^{0,*}, K^{1/2} = O(-2)
check("C2 CP^3: index(D x O(n)) = (n^3 - n)/6 from A-hat, equals chi(O(n-2)) = binom(n+1,3) from HRR",
      sp.simplify(idx_ahat - (n ** 3 - n) / 6) == 0 and sp.simplify(idx_hrr - idx_ahat) == 0
      and sp.simplify(sp.expand_func(sp.binomial(n + 1, 3)) - idx_ahat) == 0)

# C3: without the gravitational term the minimal level on spin manifolds is 6; with it, 1
ints = [((n_ ** 3 - n_) % 6 == 0) for n_ in range(-20, 21)]
check("C3 (k/6) int c1^3 with c1^3 = 1 on (CP^3, O(1)) needs k in 6Z; (k/6)c1^3 - (k/24)c1 p1 = k index is integral for k in Z",
      not sp.Rational(1, 6).is_integer and all(ints))
# product of spheres: int (a x + b y + c z)^3 = 6abc, p_1 = 0: no obstruction there (shows CP^3 is needed)
a, b, c = sp.symbols("a b c", integer=True)
x, y, z = sp.symbols("x y z")
cube = sp.Poly(sp.expand((a * x + b * y + c * z) ** 3), x, y, z).coeff_monomial(x * y * z)
check("C3b on S^2 x S^2 x S^2 int c1^3 = 6abc (always 6Z), so the k in 6Z condition comes from manifolds like CP^3",
      sp.simplify(cube - 6 * a * b * c) == 0)

# C4: the gravitational partner in 5D: 2 pi (k/24) c1 p1 with p1 = -(1/8 pi^2) tr R^R
trRR = sp.Symbol("trRR")
grav = 2 * sp.pi * k / 24 * (F / (2 * sp.pi)) * (-trRR / (8 * sp.pi ** 2))
info(f"C4 6D gravitational partner = {sp.simplify(grav)}  ->  5D term -(k/(192 pi^2)) A^tr(R^R) (sign convention dependent)")
check("C4 6D partner = -(k/(192 pi^2)) F^tr(R^R), i.e. 5D term -(k/(192 pi^2)) A^tr(R^R)",
      sp.simplify(grav + k * F * trRR / (192 * sp.pi ** 2)) == 0)

# C5: 4D theta-term from the Wilson line: (k/8 pi^2) theta f^f; shift theta -> theta + 2 pi on spin 4-manifolds
q = sp.Symbol("q2", integer=True)                      # int c1^2, even on spin 4-manifolds
shift = (k / (8 * sp.pi ** 2)) * 2 * sp.pi * (4 * sp.pi ** 2) * (2 * q)
check("C5 theta -> theta + 2 pi shifts S by 2 pi k q (int c1^2 = 2q on spin 4-manifolds): period 2 pi for |k| = 1",
      sp.simplify(shift - 2 * sp.pi * k * q) == 0)

# C6: beta_max = k alpha / 2 (Delta theta <= pi), from data alpha
bmax1 = ALPHA / 2 / DEG
bmax6 = 6 * bmax1
need1 = BETA_ACT / (ALPHA / 2) * 180
need6 = BETA_ACT / (6 * ALPHA / 2) * 180
info(f"C6 beta_max: |k| = 1 -> {bmax1:.3f} deg (needs Delta theta = {need1 / 180:.2f} pi); "
     f"k in 6Z -> {bmax6:.3f} deg (k = 6 needs Delta theta = {need6 / 180:.3f} pi)")
check("C6 minimal |k| = 1 (spin + gravitational CS term): beta_max = 0.209 deg", abs(bmax1 - 0.209) < 1e-3, f"{bmax1:.4f}")
