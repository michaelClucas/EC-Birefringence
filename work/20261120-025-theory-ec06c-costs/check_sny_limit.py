"""EC-06d (SNY 0+ PGT) closure and the GR-like regime (task 025).

Q1  SNY Eq. 20 (covariant, valid off FRW when P = Q = 0): grad R = -(2/3) Y T, Y = R + 6mu/b
    => T_mu = -(3/2) d_mu ln|Y|, closed for any Y(x): the CFJ term is gauge invariant (an ALP-like coupling
    -(3 xi/4) ln|Y| F F~ up to a boundary term, with g_phi-gamma <-> -xi via CONVENTIONS p = -g d phi).
Q2  GR-like regime (heavy 0+ mode, R << Yc = 6mu/b, R = rho for dust with 8 pi G = 1, a0 = -1):
    eps = Phi/H = (9/2) R/(Yc + R),  beta/xi = (3/4) ln(Y*/Y0) ~ eps*/6.
Q3  Friedmann constraint (SNY Eq. 47) in that regime: H^2 = (rho/3)(1 + delta), delta = c1 * R b/mu + ...;
    so eps* and the fractional H(z*) shift are tied: xi ~ beta / (eps*/6) = O(1e-2/delta) for a
    GR-compatible recombination (delta <~ 1e-2 gives xi >~ O(1)).
"""
import sympy as sp

from common import check, is_closed, summary

t, x, y, z = sp.symbols("t x y z", real=True)
X = [t, x, y, z]
Yf = sp.Function("Y")(*X)
T = [-sp.Rational(3, 2) * sp.diff(sp.log(Yf), v) for v in X]
check("Q1 T_mu = -(3/2) d ln Y is closed off FRW (any Y(x))", is_closed(T, X))
# and it solves SNY Eq. 20 with R = Y - 6mu/b
mu, b = sp.symbols("mu b", positive=True)
R = Yf - 6 * mu / b
check("Q1b it solves grad R + (2/3) Y T = 0 identically",
      all(sp.simplify(sp.diff(R, v) + sp.Rational(2, 3) * Yf * T[i]) == 0 for i, v in enumerate(X)))

# Q2: dust GR-like history R = rho = r0 a^-3, a(t), H = a'/a
a = sp.Function("a", positive=True)(t)
r0, Yc = sp.symbols("r0 Y_c", positive=True)
Rt = r0 / a**3
Phi = -sp.Rational(3, 2) * sp.diff(Rt, t) / (Yc + Rt)
eps = sp.simplify(Phi / (sp.diff(a, t) / a))
check("Q2a eps = Phi/H = (9/2) R/(Yc + R)", sp.simplify(eps - sp.Rational(9, 2) * Rt / (Yc + Rt)) == 0)
s, s0 = sp.symbols("s s0", positive=True)  # s = R/Yc at the source, s0 today
bx = sp.Rational(3, 4) * sp.log((1 + s) / (1 + s0))
ser = sp.series(bx.subs(s0, 0), s, 0, 2).removeO()
check("Q2b beta/xi = (3/4) ln(Y*/Y0) -> (3/4) s* = eps*/6 for s* << 1, s0 -> 0",
      sp.simplify(ser - sp.Rational(9, 2) * s / 6) == 0)

# Q3: Friedmann constraint (b/18) Y (3H - Phi)^2 - (b/24) R^2 - 3 a1 H^2 = rho, a1 = mu + a0 = mu - 1
H2, rho = sp.symbols("H2 rho", positive=True)
a1 = mu - 1
Ycv = 6 * mu / b
Rv = rho  # GR leading order
epsv = sp.Rational(9, 2) * Rv / (Ycv + Rv)
con = b / 18 * (Ycv + Rv) * (3 - epsv)**2 * H2 - b / 24 * Rv**2 - 3 * a1 * H2 - rho
H2sol = sp.solve(con, H2)[0]
lam = sp.Symbol("lambda", positive=True)
delta = sp.series(sp.simplify(H2sol.subs(rho, lam * rho) / (lam * rho / 3)) - 1, lam, 0, 2).removeO().subs(lam, 1)
delta = sp.simplify(delta)
print("   delta = H^2/(rho/3) - 1 to first order:", sp.factor(delta))
c1 = sp.simplify(delta / (rho * b / mu))
print("   delta / (rho b / mu) =", c1)
check("Q3a GR limit: H^2 -> rho/3 as rho b/mu -> 0", sp.limit(H2sol / (rho / 3), rho, 0) == 1)
check("Q3b first-order shift delta is linear in s = rho b/(6 mu) with an O(1) coefficient",
      sp.simplify(sp.diff(c1, rho)) == 0 and c1 != 0)
# beta/xi = (3/4) s*, delta = c1 * 6 s*  =>  xi = beta / ((3/4) delta /(6 c1))
c1n = complex(c1.subs(mu, 1)).real
xi_per = lambda d, beta=3.75e-3: beta / (0.75 * d / (6 * abs(c1n)))
print(f"   with mu = 1: |delta| = {6 * abs(c1n):.2f} s*; xi needed for |delta(z*)| = 1e-2: {xi_per(1e-2):.2f}")
check("Q3c GR-compatible recombination (|delta| <= 1e-2) needs xi >= 0.1", xi_per(1e-2) >= 0.1)
summary()
