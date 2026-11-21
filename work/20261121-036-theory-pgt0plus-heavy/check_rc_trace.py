"""036 step A: Riemann-Cartan scalar curvature with pure-trace torsion, in CONVENTIONS.md conventions.

Gamma^l_{mn} = Christoffel + K^l_{mn}, derivative index m first; T^l_{mn} = Gamma^l_{mn} - Gamma^l_{nm};
K_{lmn} = (1/2)(T_{lmn} + T_{nlm} - T_{mnl}); R^r_{smn} = d_m G^r_{ns} - d_n G^r_{ms} + G^r_{ml} G^l_{ns} - G^r_{nl} G^l_{ms};
Ric_{sn} = R^m_{smn}; R = g^{sn} Ric_{sn}. Pure trace: T^l_{mn} = (1/3)(delta^l_n T_m - delta^l_m T_n).

A1  FRW Levi-Civita check: R = 6(a''/a + a'^2/a^2) (cosmic time, -+++).
A2  T^l_{ml} = T_m for the pure-trace ansatz, and the connection reproduces it.
A3  Flat space, generic T_m(x): R = 2 d_m T^m - (2/3) T_m T^m.
A4  FRW with trace torsion T_t = Phi(t): R = R_LC + 2 nabla_LC.T - (2/3) T^2 (covariance spot check).
A5  SNY torsion-square combination (their index order T_{ij}^k, antisymmetric in the first two:
    SNY T_{nsm} = our T_{mns}) on pure trace: T.T + 2 T_{nsm}T^{msn} - 4 T_m T^m = -(8/3) T^2,
    so (a1/8)(...) = -(a1/3) T^2.
"""
import sympy as sp

from common import check, summary

t, x, y, z = sp.symbols("t x y z", real=True)
X = [t, x, y, z]
R4 = range(4)


def christoffel(g):
    gi = g.inv()
    return [[[sp.simplify(sum(gi[l, r] * (sp.diff(g[r, m], X[n]) + sp.diff(g[r, n], X[m]) - sp.diff(g[m, n], X[r]))
                              for r in R4) / 2) for n in R4] for m in R4] for l in R4]


def trace_torsion(Tv, g):
    """T^l_{mn} for pure trace with covector Tv."""
    d = sp.eye(4)
    return [[[sp.Rational(1, 3) * (d[l, n] * Tv[m] - d[l, m] * Tv[n]) for n in R4] for m in R4] for l in R4]


def contortion(Tup, g):
    gi = g.inv()
    Tlo = [[[sum(g[l, r] * Tup[r][m][n] for r in R4) for n in R4] for m in R4] for l in R4]
    Klo = [[[(Tlo[l][m][n] + Tlo[n][l][m] - Tlo[m][n][l]) / 2 for n in R4] for m in R4] for l in R4]
    return [[[sum(gi[l, r] * Klo[r][m][n] for r in R4) for n in R4] for m in R4] for l in R4]


def scalar_curv(G, g):
    gi = g.inv()
    Rs = 0
    for s in R4:
        for n in R4:
            if gi[s, n] == 0:
                continue
            ric = 0
            for m in R4:
                ric += (sp.diff(G[m][n][s], X[m]) - sp.diff(G[m][m][s], X[n])
                        + sum(G[m][m][l] * G[l][n][s] - G[m][n][l] * G[l][m][s] for l in R4))
            Rs += gi[s, n] * ric
    return sp.simplify(Rs)


# A1
a = sp.Function("a", positive=True)(t)
gF = sp.diag(-1, a**2, a**2, a**2)
GF = christoffel(gF)
RF = scalar_curv(GF, gF)
check("A1 FRW Levi-Civita R = 6(a''/a + (a'/a)^2)",
      sp.simplify(RF - 6 * (sp.diff(a, t, 2) / a + sp.diff(a, t)**2 / a**2)) == 0)

# A2, A3 flat
eta = sp.diag(-1, 1, 1, 1)
Tv = [sp.Function(f"T{i}")(*X) for i in R4]
Tup = trace_torsion(Tv, eta)
check("A2a trace of ansatz T^l_{ml} = T_m", all(sp.simplify(sum(Tup[l][m][l] for l in R4) - Tv[m]) == 0 for m in R4))
K = contortion(Tup, eta)
G = K
check("A2b connection torsion G^l_{mn} - G^l_{nm} = T^l_{mn}",
      all(sp.simplify(G[l][m][n] - G[l][n][m] - Tup[l][m][n]) == 0 for l in R4 for m in R4 for n in R4))
Kl = [[[sum(eta[l, r] * K[r][m][n] for r in R4) for n in R4] for m in R4] for l in R4]
check("A2c metric compatibility K_{lmn} = -K_{nml}",
      all(sp.simplify(Kl[l][m][n] + Kl[n][m][l]) == 0 for l in R4 for m in R4 for n in R4))
Rflat = scalar_curv(G, eta)
divT = sum(eta.inv()[m, m] * sp.diff(Tv[m], X[m]) for m in R4)
T2 = sum(eta.inv()[m, m] * Tv[m]**2 for m in R4)
check("A3 flat: R = 2 d.T - (2/3) T^2", sp.simplify(Rflat - (2 * divT - sp.Rational(2, 3) * T2)) == 0)

# A4 FRW + trace torsion along t
Phi = sp.Function("Phi")(t)
TvF = [Phi, 0, 0, 0]
TupF = trace_torsion(TvF, gF)
KF = contortion(TupF, gF)
GtF = [[[GF[l][m][n] + KF[l][m][n] for n in R4] for m in R4] for l in R4]
RtF = scalar_curv(GtF, gF)
sqg = a**3
divLC = sp.diff(sqg * (-Phi), t) / sqg  # nabla_mu T^mu, T^t = -Phi
check("A4 FRW: R = R_LC + 2 nabla.T - (2/3) T^2 (T^2 = -Phi^2)",
      sp.simplify(RtF - (RF + 2 * divLC - sp.Rational(2, 3) * (-Phi**2))) == 0)

# A5 SNY combination
Tlo = [[[sum(eta[l, r] * Tup[r][m][n] for r in R4) for n in R4] for m in R4] for l in R4]
ei = eta.inv()


def up3(A, i, j, k):
    return ei[i, i] * ei[j, j] * ei[k, k] * A[i][j][k]


S = lambda n, s, m: Tlo[m][n][s]  # SNY T_{nsm} = our T_{mns}
t1 = sum(S(n, s, m) * ei[n, n] * ei[s, s] * ei[m, m] * S(n, s, m) for n in R4 for s in R4 for m in R4)
t2 = sum(S(n, s, m) * ei[n, n] * ei[s, s] * ei[m, m] * S(m, s, n) for n in R4 for s in R4 for m in R4)
comb = sp.simplify(t1 + 2 * t2 - 4 * T2)
check("A5 SNY combination on pure trace = -(8/3) T^2", sp.simplify(comb + sp.Rational(8, 3) * T2) == 0, f"({comb})")
summary()
