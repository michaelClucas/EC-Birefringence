"""Skeptic 026: Weitzenboeck torsion of flat FRW tetrads, written independently of 022.

Conventions (CONVENTIONS.md): Gamma^l_{mu nu} = e_A^l d_mu e^A_nu (derivative index first), so with zero spin
connection T^A_{mu nu} = d_mu e^A_nu - d_nu e^A_mu. Work in the orthonormal frame: T^A_{BC} = e_B^mu e_C^nu T^A_{mu nu},
trace T_C = T^A_{CA}, axial S^D = (1/6) eps^{DABC} T_{ABC} with frame eps^{0123} = -1 (eps_{0123} = +1), then
T_mu = e^C_mu T_C and S^mu = e_D^mu S^D. Coordinates (tau, x, y, z). Perturbations to first order in eps.
"""
import itertools

import sympy as sp
from sympy.combinatorics import Permutation

tau, x, y, z, eps = sp.symbols("tau x y z epsilon")
X = [tau, x, y, z]
a = sp.Function("a")(tau)
etaM = sp.diag(-1, 1, 1, 1)


def trunc(q):
    q = sp.expand(q)
    return q.coeff(eps, 0) + eps * q.coeff(eps, 1)


def torsion(e):
    """e[A, mu] linear in eps. Returns metric, T_mu (coordinate, lower), S^mu (coordinate, upper)."""
    i0 = e.subs(eps, 0).inv()
    einv = (i0 - eps * i0 * sp.diff(e, eps) * i0).applyfunc(sp.expand)      # einv[mu, A] to O(eps)
    TA = [[[sp.diff(e[A, n], X[m]) - sp.diff(e[A, m], X[n]) for n in range(4)] for m in range(4)] for A in range(4)]
    TF = [[[trunc(sum(einv[m, B] * einv[n, C] * TA[A][m][n] for m in range(4) for n in range(4)))
            for C in range(4)] for B in range(4)] for A in range(4)]
    Tlow = [[[sum(etaM[D, A] * TF[A][B][C] for A in range(4)) for C in range(4)] for B in range(4)] for D in range(4)]
    Tc = [sum(TF[A][C][A] for A in range(4)) for C in range(4)]
    Sf = []
    for D in range(4):
        s = 0
        for p in itertools.permutations(range(4)):
            if p[0] == D:
                s += -Permutation(list(p)).signature() * Tlow[p[1]][p[2]][p[3]]
        Sf.append(s / 6)
    Tmu = [sp.simplify(trunc(sum(e[C, m] * Tc[C] for C in range(4)))) for m in range(4)]
    Smu = [sp.simplify(trunc(sum(einv[m, D] * Sf[D] for D in range(4)))) for m in range(4)]
    g = (e.T * etaM * e).applyfunc(trunc)
    return g, Tmu, Smu


def curl(V):
    return {(m, n): sp.simplify(sp.diff(V[n], X[m]) - sp.diff(V[m], X[n])) for m in range(4) for n in range(m + 1, 4)}


res = []


def chk(name, ok):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + ": " + name)


def first(v):
    return [sp.simplify(sp.diff(q, eps)) for q in v]


# 1. diagonal flat FRW tetrad
e0 = sp.diag(a, a, a, a)
g0, T0, S0 = torsion(e0)
print("diag tetrad: T_mu =", T0, " S^mu =", S0)
chk("T_tau = 3a'/a, T_i = 0 (T_t = T_tau/a = 3H, eps_T = 3)",
    sp.simplify(T0[0] - 3 * sp.diff(a, tau) / a) == 0 and all(c == 0 for c in T0[1:]))
chk("S = 0 on the diagonal tetrad", all(c == 0 for c in S0))

# 2. local boost in the tau-x plane, rapidity eps*eta(tau, x, y)
eta = sp.Function("eta")(tau, x, y)
L = sp.eye(4); L[0, 1] = eps * eta; L[1, 0] = eps * eta
g1, T1, S1 = torsion(L * e0)
chk("metric unchanged to O(eps)", (g1 - g0).applyfunc(sp.simplify) == sp.zeros(4))
dT, dS = first(T1), first(S1)
print("boost: delta T_mu =", dT, "; delta S^mu =", dS)
c = curl(dT)
print("d(delta T) nonzero components:", {k: v for k, v in c.items() if v != 0})
chk("delta T_mu = (-eta_x, -eta_tau, 0, 0)",
    sp.simplify(dT[0] + sp.diff(eta, x)) == 0 and sp.simplify(dT[1] + sp.diff(eta, tau)) == 0 and dT[2] == 0 and dT[3] == 0)
chk("delta T not closed for generic eta", any(v != 0 for v in c.values()))
chk("delta S^z = 2 eta_y/(3a^2), other components 0 (022's T2)",
    sp.simplify(dS[3] - 2 * sp.diff(eta, y) / (3 * a**2)) == 0 and all(dS[i] == 0 for i in range(3)))

# 2b. time-only boost: rapidity h(tau): changes the spatial trace only
ht = sp.Function("h")(tau)
L2 = sp.eye(4); L2[0, 1] = eps * ht; L2[1, 0] = eps * ht
_, T2, S2 = torsion(L2 * e0)
d2 = first(T2)
print("time-only boost: delta T_mu =", d2, "; delta S =", first(S2))
chk("time-only boost: T_tau unchanged, T_x = -h': line-of-sight beta shift (xi/2) n^x [h] is a pure dipole",
    d2[0] == 0 and sp.simplify(d2[1] + sp.diff(ht, tau)) == 0)

# 3. scalar perturbations, diagonal tetrad (Newtonian gauge, ds^2 = a^2[-(1+2Phi)dtau^2 + (1-2Psi)dx^2])
Phi = sp.Function("Phi")(tau, x, y, z); Psi = sp.Function("Psi")(tau, x, y, z)
ep = sp.diag(a * (1 + eps * Phi), a * (1 - eps * Psi), a * (1 - eps * Psi), a * (1 - eps * Psi))
_, Tp, Sp_ = torsion(ep)
Tp1 = first(Tp)
print("perturbed diagonal tetrad: delta T_mu =", Tp1, "; delta S =", first(Sp_))
chk("delta T_tau = -3 Psi', delta T_i = d_i(Phi - 2 Psi)",
    sp.simplify(Tp1[0] + 3 * sp.diff(Psi, tau)) == 0 and sp.simplify(Tp1[1] - sp.diff(Phi - 2 * Psi, x)) == 0)
chk("(dT)_{tau x} = d_x d_tau (Phi + Psi) (022's T5)", sp.simplify(curl(Tp1)[(0, 1)] - sp.diff(Phi + Psi, x, tau)) == 0)

# 4. add the extra Lorentz scalar zeta of f(T) perturbation theory: boost with rapidity eps * d_i zeta
zeta = sp.Function("zeta")(tau, x, y, z)
Lz = sp.eye(4)
for i in (1, 2, 3):
    Lz[0, i] = eps * sp.diff(zeta, X[i]); Lz[i, 0] = eps * sp.diff(zeta, X[i])
_, Tz, Sz = torsion(Lz * ep)
Tz1 = first(Tz)
extra = sp.simplify(curl(Tz1)[(0, 1)] - sp.diff(Phi + Psi, x, tau))
print("with zeta: delta T_mu =", Tz1, "; delta S =", first(Sz))
print("extra (dT)_{tau x} from zeta:", extra)
chk("zeta adds a non-zero piece to (dT)_{tau i}: 022's formula is incomplete unless zeta = 0", extra != 0)
print(f"{sum(res)} PASS, {len(res) - sum(res)} FAIL")
