"""Gauge invariance of the CFJ term L = -1/2 p_mu A_nu Ft^{mu nu} for a non-constant p (task 022, Q2).

G1  delta_chi L = d_nu(...) + (1/4) chi (dp)_{nu mu} Ft^{mu nu} with (dp)_{nu mu} = d_nu p_mu - d_mu p_nu:
    gauge invariant (up to a boundary term) iff dp = 0.
G2  Euler-Lagrange expressions E^nu of -F^2/4 + L_CFJ satisfy identically
    d_nu E^nu = -(1/4)(dp)_{nu mu} Ft^{mu nu}; with a conserved current J the equations E^nu = -J^nu then force
    (dp).Ft = 0 on every solution: the system is inconsistent unless dp = 0.
G3  FRW background p = p_tau(tau) dtau is closed; a perturbed p is closed iff it is a gradient
    (delta p_tau = d_tau s, delta p_i = d_i s); a transverse vector piece or an independent delta p_tau breaks it.
G4  Stueckelberg completion A -> A + d sigma: the CFJ term produces sigma (dp).Ft / 4, a mixing of the
    longitudinal mode with the photon proportional to dp (no kinetic term for sigma without a photon mass).
Flat metric suffices (the CFJ term is metric independent: p_mu A_nu Ft^{mu nu} sqrt(-g) has no metric).
"""
import sympy as sp

from common import R4, check, summary

X = sp.symbols("t x y z", real=True)
A = [sp.Function(f"A{m}")(*X) for m in R4]
p = [sp.Function(f"p{m}")(*X) for m in R4]
chi = sp.Function("chi")(*X)
eta = sp.diag(-1, 1, 1, 1)
EPS_UP = lambda a, b, c, d: -sp.LeviCivita(a, b, c, d)  # eps^{0123} = -1 (flat, CONVENTIONS.md)


def Fl(Av):
    return [[sp.diff(Av[n], X[m]) - sp.diff(Av[m], X[n]) for n in R4] for m in R4]


def Ftu(Av):
    F = Fl(Av)
    return [[sp.Rational(1, 2) * sum(EPS_UP(m, n, r, s) * F[r][s] for r in R4 for s in R4) for n in R4] for m in R4]


def Lcfj(Av, pv):
    Ft = Ftu(Av)
    return -sp.Rational(1, 2) * sum(pv[m] * Av[n] * Ft[m][n] for m in R4 for n in R4)


# G1
Ft = Ftu(A)
dL = sp.expand(Lcfj([A[m] + sp.diff(chi, X[m]) for m in R4], p) - Lcfj(A, p))
bound = sum(sp.diff(-sp.Rational(1, 2) * sum(p[m] * chi * Ft[m][n] for m in R4), X[n]) for n in R4)
dp = [[sp.diff(p[m], X[n]) - sp.diff(p[n], X[m]) for m in R4] for n in R4]  # dp[n][m] = d_n p_m - d_m p_n
bulk = sp.Rational(1, 4) * chi * sum(dp[n][m] * Ft[m][n] for m in R4 for n in R4)
check("G1 delta_chi L_CFJ = total derivative + (1/4) chi (dp)_{nu mu} Ft^{mu nu}", sp.simplify(dL - bound - bulk) == 0)

# G2 Euler-Lagrange of L = -1/4 F^2 + L_CFJ
F = Fl(A)
Fu = [[sum(eta[m, a] * eta[n, b] * F[a][b] for a in R4 for b in R4) for n in R4] for m in R4]
L = -sp.Rational(1, 4) * sum(F[m][n] * Fu[m][n] for m in R4 for n in R4) + Lcfj(A, p)
E = []
for nu in R4:
    dLdA = sp.diff(L, A[nu])
    dLddA = sum(sp.diff(sp.diff(L, sp.Derivative(A[nu], X[m])), X[m]) for m in R4)
    E.append(sp.expand(dLdA - dLddA))
divE = sp.expand(sum(sp.diff(E[n], X[n]) for n in R4))
target = sp.expand(-sp.Rational(1, 2) * sum(dp[n][m] * Ft[m][n] for m in R4 for n in R4))
ratio = sp.simplify(divE / target) if target != 0 else None
print(f"   d_nu E^nu / [-(1/2)(dp).Ft] = {ratio}")
check("G2 divergence of the field equations is proportional to (dp).Ft (consistency needs dp = 0)",
      ratio is not None and ratio.is_number and ratio != 0)
pconst = {p[m]: sp.Symbol(f"c{m}") for m in R4}
check("G2b for constant p the divergence vanishes identically",
      sp.simplify(divE.subs(pconst).doit()) == 0)

# G3 FRW background and perturbations
tau = X[0]
pb = [sp.Function("pt")(tau), 0, 0, 0]
closed = lambda q: all(sp.simplify(sp.diff(q[m], X[n]) - sp.diff(q[n], X[m])) == 0 for m in R4 for n in R4)
check("G3a FRW background p = p_tau(tau) dtau is closed", closed(pb))
s = sp.Function("s")(*X)
grad = [pb[0] + sp.diff(s, tau)] + [sp.diff(s, X[i]) for i in (1, 2, 3)]
check("G3b gradient perturbation (scalar potential s) keeps dp = 0", closed(grad))
q = sp.Function("q")(*X)
indep = [pb[0] + q] + [sp.diff(s, X[i]) for i in (1, 2, 3)]
check("G3c independent scalar perturbation delta p_tau = q != d_tau s breaks dp = 0", not closed(indep))
v = sp.Function("v")(tau, X[3])
curl = [pb[0], v, 0, 0]  # transverse vector mode v_x(tau, z), k along z
check("G3d transverse vector perturbation breaks dp = 0", not closed(curl))

# G4 Stueckelberg
sig = sp.Function("sigma")(*X)
dLs = sp.expand(Lcfj([A[m] + sp.diff(sig, X[m]) for m in R4], p) - Lcfj(A, p))
bounds = sum(sp.diff(-sp.Rational(1, 2) * sum(p[m] * sig * Ft[m][n] for m in R4), X[n]) for n in R4)
mix = sp.Rational(1, 4) * sig * sum(dp[n][m] * Ft[m][n] for m in R4 for n in R4)
check("G4 Stueckelberg: CFJ(A + d sigma) - CFJ(A) = boundary + sigma (dp).Ft / 4", sp.simplify(dLs - bounds - mix) == 0)
summary()
