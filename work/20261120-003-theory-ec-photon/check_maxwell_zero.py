"""Task 2(a,b): minimally coupled Maxwell does not see torsion; 'torsionful minimal coupling' breaks U(1).

(a) F = dA: d_m A_n - d_n A_m contains no connection. With the Levi-Civita connection the
    covariant curl equals dA (Christoffel symmetric). With the torsionful Gamma:
    nabla_m A_n - nabla_n A_m = F_mn - T^l_mn A_l.
(b) Under A -> A + d chi, F' = F - T A changes by -T^l_mn d_l chi: not gauge invariant.
(c) FRW: sqrt(-g) g^{ma} g^{nb} = eta^{ma} eta^{nb} for g = a^2 eta (Maxwell is conformally
    invariant in 4D), and sqrt(-g) eps^{mnrs} = -[mnrs] is metric independent, so the comoving
    Maxwell (and CFJ) equations are the flat ones. Torsion never enters.
(d) With FRW torsion (axial S^0, trace h), the gauge-breaking Lagrangian -1/4 F'F' equals
    -1/4 FF - 1/2 p_m A_n Ftilde^{mn} with p_m = 2 S_m, plus derivative-free A^2 terms
    a parity-even gauge-breaking term (1/3) T_m A_n F^{mn}, and (h, S)-dependent
    photon 'mass' terms. Checked in flat space (a = 1).
"""
import itertools
import sys

import sympy as sp

from common import ETA, R4, axial_vector, perm, report, trace_vector, zeros3

X = sp.symbols("x0:4")


def main():
    ok = True
    A = [sp.Function(f"A{i}")(*X) for i in R4]
    chi = sp.Function("chi")(*X)
    Gam = [[[sp.Function(f"G{l}{m}{n}")(*X) for n in R4] for m in R4] for l in R4]
    T = [[[Gam[l][m][n] - Gam[l][n][m] for n in R4] for m in R4] for l in R4]

    def covcurl(Av, G):
        # nabla_m A_n = d_m A_n - Gamma^l_{mn} A_l
        return [[sp.diff(Av[n], X[m]) - sum(G[l][m][n] * Av[l] for l in R4)
                 - (sp.diff(Av[m], X[n]) - sum(G[l][n][m] * Av[l] for l in R4))
                 for n in R4] for m in R4]

    F = [[sp.diff(A[n], X[m]) - sp.diff(A[m], X[n]) for n in R4] for m in R4]
    Fp = covcurl(A, Gam)
    ok &= report("(a) nabla_m A_n - nabla_n A_m = F_mn - T^l_mn A_l (generic Gamma)",
                 all(sp.expand(Fp[m][n] - (F[m][n] - sum(T[l][m][n] * A[l] for l in R4))) == 0
                     for m in R4 for n in R4))
    Gsym = [[[sp.Function(f"C{l}{min(m, n)}{max(m, n)}")(*X) for n in R4] for m in R4] for l in R4]
    Fs = covcurl(A, Gsym)
    ok &= report("(a) symmetric (Levi-Civita) connection drops out: curl = dA",
                 all(sp.expand(Fs[m][n] - F[m][n]) == 0 for m in R4 for n in R4))
    Ag = [A[i] + sp.diff(chi, X[i]) for i in R4]
    Fpg = covcurl(Ag, Gam)
    dF = [[sp.expand(Fpg[m][n] - Fp[m][n] + sum(T[l][m][n] * sp.diff(chi, X[l]) for l in R4))
           for n in R4] for m in R4]
    ok &= report("(b) gauge variation of F' is -T^l_mn d_l chi (nonzero for T != 0)",
                 all(v == 0 for row in dF for v in row))
    Fg = [[sp.diff(Ag[n], X[m]) - sp.diff(Ag[m], X[n]) for n in R4] for m in R4]
    ok &= report("(b) F = dA is gauge invariant",
                 all(sp.expand(Fg[m][n] - F[m][n]) == 0 for m in R4 for n in R4))

    # (c) conformal invariance in FRW
    a = sp.Function("a", positive=True)(X[0])
    g = a**2 * ETA
    gi = g.inv()
    sq = sp.sqrt(-g.det())
    ok &= report("(c) sqrt(-g) g^ma g^nb = eta^ma eta^nb for g = a^2 eta",
                 all(sp.simplify(sq * gi[m, al] * gi[n, be] - ETA[m, al] * ETA[n, be]) == 0
                     for m, n, al, be in itertools.product(R4, R4, R4, R4)))
    ok &= report("(c) sqrt(-g) eps^{0123} = -1 independent of a",
                 sp.simplify(sq * (-perm(0, 1, 2, 3) / sq) + 1) == 0)

    # (d) FRW torsion in flat space (a = 1): S^0 = s, trace T^i_{0j} = h delta_ij
    s, h = sp.symbols("s h")
    Tf = zeros3()
    for i in (1, 2, 3):
        Tf[i][0][i] = h
        Tf[i][i][0] = -h
    for i, j, k in itertools.product((1, 2, 3), repeat=3):
        # T_{ijk} = eps_{ijk0} S^0 = -[0ijk] s ; index i raised with eta_ii = +1
        Tf[i][j][k] = -perm(0, i, j, k) * s
    Sv = [sp.simplify(c) for c in axial_vector(ETA, Tf)]
    ok &= report("(d) axial vector of the test torsion is S^m = (s,0,0,0)", Sv == [s, 0, 0, 0])
    print("INFO: trace vector T_m =", trace_vector(Tf))
    Fpf = [[F[m][n] - sum(Tf[l][m][n] * A[l] for l in R4) for n in R4] for m in R4]

    def sq_(Fm):
        return sum(ETA[m, m] * ETA[n, n] * Fm[m][n] ** 2 for m in R4 for n in R4)

    Ldiff = sp.expand(-sp.Rational(1, 4) * sq_(Fpf) + sp.Rational(1, 4) * sq_(F))
    Slow = [ETA[m, m] * Sv[m] for m in R4]
    p = [2 * Slow[m] for m in R4]
    # Ftilde^{mn} = 1/2 eps^{mnrs} F_rs, eps^{mnrs} = -[mnrs] (flat)
    Ft = [[sp.Rational(1, 2) * sum(-perm(m, n, r, ss) * F[r][ss] for r in R4 for ss in R4)
           for n in R4] for m in R4]
    LCFJ = sp.expand(-sp.Rational(1, 2) * sum(p[m] * A[n] * Ft[m][n] for m in R4 for n in R4))
    Tv = trace_vector(Tf)
    Fup = [[ETA[m, m] * ETA[n, n] * F[m][n] for n in R4] for m in R4]
    Ltr = sp.expand(sp.Rational(1, 3) * sum(Tv[m] * A[n] * Fup[m][n] for m in R4 for n in R4))
    rest = sp.expand(Ldiff - LCFJ - Ltr)
    has_deriv = rest.has(sp.Derivative)
    ok &= report("(d) -1/4 F'F' + 1/4 FF = CFJ(p_m = 2 S_m) + (1/3) T_m A_n F^mn "
                 "+ derivative-free A^2 terms", not has_deriv)
    print("INFO: derivative-free remainder (photon 'mass' terms):", sp.factor(rest))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
