"""Symbolic checks (sympy) of the beta bookkeeping for EC-06 (CONVENTIONS.md: beta = +1/2 int p_tau dtau).

S1 p_tau = a p_0hat for the FRW tetrad e^0hat = a dtau (lower indices); raising flips the sign (eta = diag(-1,1,1,1)).
S2 p_0hat = eps H (H = a'/a^2 in conformal time): 1/2 int p_tau dtau = (eps/2) ln(a_0/a_*) = (eps/2) ln(1+z*).
S3 p_0hat = eps H (1+z)^n: 1/2 int p_0hat dt = eps ((1+z*)^n - 1)/(2n), and -> (eps/2) ln(1+z*) as n -> 0.
S4 constant p_0hat: 1/2 int p_0hat dt = (p/2) t_L (lookback time), matching KM 2009 Eq. 158 with p = -2 k_AF (|beta| = |k_AF| t_L).
S5 a CFJ vector p_mu = (p_tau(tau), 0, 0, 0) is closed (dp = 0), so -1/2 p A F~ is gauge invariant up to a total derivative.
"""
import sympy as sp

checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}: {name}")


tau, eps, n = sp.symbols("tau epsilon n", real=True)
a = sp.Function("a", positive=True)(tau)
# S1: tetrad e^A_mu = a delta^A_mu; p_mu = e^A_mu p_A
p0hat = sp.Symbol("p0hat")
e = sp.diag(a, a, a, a)
p_hat = sp.Matrix([p0hat, 0, 0, 0])
p_coord = e.T * p_hat
check("S1a p_tau = a p_0hat (lower index, tetrad e^0hat = a dtau)", sp.simplify(p_coord[0] - a * p0hat) == 0)
eta = sp.diag(-1, 1, 1, 1)
check("S1b p^0hat = -p_0hat (signature -+++)", (eta * p_hat)[0] == -p0hat)

# S2: H = (da/dtau)/a^2, p_tau = a eps H = eps a'/a
H = sp.diff(a, tau) / a**2
integrand = sp.Rational(1, 2) * a * eps * H
x = sp.Symbol("x", positive=True)  # x = a
# change variable tau -> a: dtau = da / a'
integrand_a = sp.simplify(integrand / sp.diff(a, tau)).subs(a, x)
astar, a0 = sp.symbols("a_star a_0", positive=True)
I2 = sp.integrate(integrand_a, (x, astar, a0))
zs = sp.Symbol("z_star", positive=True)
check("S2 1/2 int p_tau dtau = (eps/2) ln(1+z*) for p_0hat = eps H",
      sp.simplify(I2.subs({a0: 1, astar: 1 / (1 + zs)}) - eps / 2 * sp.log(1 + zs)) == 0)

# S3: dt = -dz/((1+z)H) -> 1/2 int eps H (1+z)^n dt = (eps/2) int_0^z* (1+z)^(n-1) dz
z = sp.Symbol("z", positive=True)
I3 = sp.integrate(sp.Rational(1, 2) * eps * (1 + z) ** (n - 1), (z, 0, zs), conds="none")
target = eps * ((1 + zs) ** n - 1) / (2 * n)
check("S3a n != 0: (eps/2n)((1+z*)^n - 1)", sp.simplify(sp.piecewise_fold(I3).rewrite(sp.Piecewise).args[0][0] - target) == 0
      if isinstance(I3, sp.Piecewise) else sp.simplify(I3 - target) == 0)
check("S3b n -> 0 limit gives (eps/2) ln(1+z*)", sp.simplify(sp.limit(target, n, 0) - eps / 2 * sp.log(1 + zs)) == 0)

# S4: constant p: 1/2 int p dt = p t_L/2; KM: |delta psi| = |k_AF| t_L with |p| = 2|k_AF|
p, tL, k = sp.symbols("p t_L k_AF", positive=True)
check("S4 |beta| = |k_AF| t_L for p = 2 k_AF", sp.simplify((p * tL / 2).subs(p, 2 * k) - k * tL) == 0)

# S5: closedness of p_mu = (f(tau), 0, 0, 0)
X = sp.symbols("tau x y w")
f = sp.Function("f")(X[0])
pv = [f, 0, 0, 0]
dp = [[sp.diff(pv[j], X[i]) - sp.diff(pv[i], X[j]) for j in range(4)] for i in range(4)]
check("S5 dp = 0 for a homogeneous timelike CFJ vector", all(sp.simplify(v) == 0 for row in dp for v in row))

print(f"\n{sum(checks)} PASS, {len(checks) - sum(checks)} FAIL")
