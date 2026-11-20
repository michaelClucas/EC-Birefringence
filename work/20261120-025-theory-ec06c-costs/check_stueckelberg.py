"""Stueckelberg / Proca repair of the CFJ term with a non-closed p (task 025, repair a). Flat space,
eps^{0123} = -1 (CONVENTIONS.md), F~^{mn} = 1/2 eps^{mnrs} F_rs.

S1  L_CFJ(A - d sigma) = L_CFJ(A) - 1/4 sigma (dp)_{nm} F~^{mn} + total derivative, for arbitrary p(x).
    So the gauge-invariant completion couples the Stueckelberg scalar to (dp).F~; with canonical
    pi = m sigma the coupling is -(1/4m) pi (dp).F~, which diverges as m -> 0.
S2  Proca + CFJ with constant timelike p = (p0, 0, 0, 0): transverse modes obey
    (w^2 - k^2 - m^2)^2 = p0^2 k^2, so w^2 >= m^2 - p0^2/4: a photon mass m > |p0|/2 removes the
    tachyonic band k < |p0| of massless timelike CFJ (Carroll-Field-Jackiw 1990; Adam-Klinkhamer hep-ph/0101087).
S3  Mixing: for (dp)_{tau x} = D constant and a photon moving along z, the pi-A_y off-diagonal element is
    D k/(2m) (up to sign); the two modes are degenerate (w^2 = k^2 + m^2), so the converted amplitude
    grows as |D| L/(4m) over a coherent path L (checked numerically on the 2-mode system).
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from common import check, summary

t, x, y, z = sp.symbols("t x y z", real=True)
X = [t, x, y, z]
eta = sp.diag(-1, 1, 1, 1)
A = [sp.Function(f"A{i}")(*X) for i in range(4)]  # lower components A_mu
p = [sp.Function(f"p{i}")(*X) for i in range(4)]  # lower components p_mu
sig = sp.Function("sigma")(*X)


def eps_up(*idx):
    return -sp.LeviCivita(*idx)  # eps^{0123} = -1


def Fdual_up(Al):
    F = [[sp.diff(Al[n], X[m]) - sp.diff(Al[m], X[n]) for n in range(4)] for m in range(4)]
    return [[sum(eps_up(m, n, r, s) * F[r][s] for r in range(4) for s in range(4)) / 2 for n in range(4)]
            for m in range(4)]


def L_cfj(Al):
    Fd = Fdual_up(Al)
    return -sp.Rational(1, 2) * sum(p[m] * Al[n] * Fd[m][n] for m in range(4) for n in range(4))


As = [A[m] - sp.diff(sig, X[m]) for m in range(4)]
Fd = Fdual_up(A)
dp = [[sp.diff(p[m], X[n]) - sp.diff(p[n], X[m]) for m in range(4)] for n in range(4)]  # dp[n][m] = (dp)_{nm}
extra = -sp.Rational(1, 4) * sig * sum(dp[n][m] * Fd[m][n] for m in range(4) for n in range(4))
diff = sp.expand(L_cfj(As) - L_cfj(A) - extra)
EL = sp.euler_equations(diff, A + [sig], X)
check("S1 L(A - d sigma) - L(A) + 1/4 sigma (dp).F~ is a total derivative (all Euler-Lagrange terms vanish)",
      all(sp.expand(e.lhs) == 0 for e in EL))

# ---------------------------------------------------------------- S2 dispersion (Proca + constant timelike CFJ)
w, k, m, p0 = sp.symbols("omega k m p0", positive=True)
eps_ = [sp.Symbol(f"e{i}") for i in range(4)]
# plane wave A_mu = e_mu exp(i(k z - w t)); field equation from L = -1/4 F^2 - 1/2 m^2 A^2 - 1/2 p.A F~
Ar = [sp.Function(f"B{i}")(*X) for i in range(4)]
Fr = [[sp.diff(Ar[n], X[mm]) - sp.diff(Ar[mm], X[n]) for n in range(4)] for mm in range(4)]
F2 = sum(eta[a_, a_] * eta[b, b] * Fr[a_][b]**2 for a_ in range(4) for b in range(4))
A2 = sum(eta[a_, a_] * Ar[a_]**2 for a_ in range(4))
Fdr = Fdual_up(Ar)
Lp = -F2 / 4 - m**2 * A2 / 2 - sp.Rational(1, 2) * p0 * sum(Ar[n] * Fdr[0][n] for n in range(4))  # p_mu = (p0,0,0,0)
eqs = [e.lhs for e in sp.euler_equations(Lp, Ar, X)]
phase = sp.exp(sp.I * (k * z - w * t))
sub = {Ar[i]: eps_[i] * phase for i in range(4)}
M = sp.Matrix(4, 4, lambda r, c: sp.simplify(sp.expand(sp.simplify(eqs[r].subs(sub).doit() / phase)).coeff(eps_[c])))
det = sp.factor(M.det())
print("   det M =", det)
target = (w**2 - k**2 - m**2)**2 - p0**2 * k**2
check("S2a det contains the transverse factor (w^2 - k^2 - m^2)^2 - p0^2 k^2", sp.simplify(sp.rem(sp.expand(det), sp.expand(target), w)) == 0)
w2min = sp.minimum(k**2 + m**2 - p0 * k, k, sp.Interval(0, sp.oo))
check("S2b min_k w^2 = m^2 - p0^2/4 (stable iff m > |p0|/2; m = 0 gives a tachyonic band 0 < k < p0)",
      sp.simplify(w2min - (m**2 - p0**2 / 4)) == 0, f"({w2min})")

# ---------------------------------------------------------------- S3 mixing element and conversion
D = sp.Symbol("D", real=True)
Ay = sp.Function("Ay")(t, z)
pi_ = sp.Function("pi")(t, z)
Al = [0, 0, Ay, 0]
Fd3 = Fdual_up(Al)
# (dp)_{tau x} = D = -(dp)_{x tau}; mixing term -(1/4m) pi (dp)_{nm} F~^{mn}
mix = -1 / (4 * m) * pi_ * (D * Fd3[1][0] - D * Fd3[0][1])
mix = sp.simplify(mix)
print("   mixing Lagrangian:", mix)
coef = sp.simplify(mix / (pi_ * sp.diff(Ay, z)))
check("S3a mixing = c * pi * d_z A_y with |c| = |D|/(2m) (off-diagonal element D k/(2m) for a mode e^{ikz})",
      sp.simplify(sp.Abs(coef) - sp.Abs(D) / (2 * m)) == 0, f"(c = {coef})")
# numeric: degenerate 2-mode system in the slowly varying envelope approx: da/dt = -i mu b, db/dt = -i mu a,
# mu = (D k/2m)/(2w) ~ D/(4m) for k ~ w >> m.  Check |b(L)| = sin(mu L) ~ D L/(4m).
Dn, mn, kn = 1e-3, 1.0, 50.0
wn = np.sqrt(kn**2 + mn**2)
off = Dn * kn / (2 * mn)


def rhs(tt, s):
    # full second-order system: (d_t^2 + w^2) A = -off * pi, (d_t^2 + w^2) pi = -off * A
    Aa, dA, P, dP = s
    return [dA, -wn**2 * Aa - off * P, dP, -wn**2 * P - off * Aa]


Ln = 400.0
sol = solve_ivp(rhs, (0, Ln), [1.0, 0.0, 0.0, 0.0], rtol=1e-10, atol=1e-12, dense_output=True)
ts = np.linspace(Ln - 2 * np.pi / wn, Ln, 400)
amp = np.max(np.abs(sol.sol(ts)[2]))
pred = abs(np.sin(off / (2 * wn) * Ln))
print(f"   numeric converted amplitude {amp:.5f}, analytic sin(D k L/(4 m w)) = {pred:.5f}, D L/(4m) = {Dn * Ln / (4 * mn):.5f}")
check("S3b degenerate mixing: converted amplitude = sin(D k L/(4 m w)) ~ |D| L/(4m) (numeric, 1%)",
      abs(amp - pred) < 0.01 * pred)
summary()
