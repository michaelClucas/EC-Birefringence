"""Task 2(c): Einstein-Cartan + Dirac at tree level.

Dirac conventions (Srednicki, mostly plus): {gamma^a, gamma^b} = -2 eta^{ab},
L_D = (i/2)(psibar gamma^m nabla_m psi - (nabla_m psibar) gamma^m psi) - m psibar psi,
psibar = psi^dagger beta, beta = gamma^0, gamma_5 = i gamma^0 gamma^1 gamma^2 gamma^3,
J_5^m = psibar gamma^m gamma_5 psi, J_5^0 = psi^dag gamma_5 psi (= n_R - n_L).
nabla_m psi = d_m psi + Omega_m psi, Omega_m = c omega_{m a b} gamma^a gamma^b, with c fixed by
[gamma^a, Omega_m] = omega_m^a_b gamma^b (so psibar gamma^a psi is a Lorentz vector).
The torsion part of omega is omega_m^a_b = K^a_{m b} (flat tetrad at the point).

Checks:
 1. c solved, and beta Omega^dag beta^-1 = -Omega (so nabla psibar = d psibar - psibar Omega).
 2. (i/2){gamma^m, Omega^K_m} = c_A S_m gamma^m gamma_5 for random full torsion: only the axial
    vector couples (trace and tensor parts drop out); c_A printed.
 3. Torsion-squared part of R(Gamma) for constant torsion in flat space:
    R_Q = b_T T_m T^m + b_S S_m S^m + (pure tensor part), no cross terms; b_T, b_S printed.
 4. Eliminating algebraic S from (1/2kappa) R + c_A S.J5: S_m = -(kappa c_A / b_S) J5_m and
    L_eff = -(kappa c_A^2 / (2 b_S)) J5.J5. No A_m appears: tree-level EC with Dirac fermions
    has no photon-torsion coupling.
"""
import itertools
import sys

import sympy as sp

from common import (ETA, R4, axial_vector, contortion_correct, lower1, perm, random_torsion,
                    report, rng, torsion_from_parts, trace_vector, zeros3)

I2 = sp.eye(2)
Z2 = sp.zeros(2)
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]


def blk(a, b, c, d):
    return sp.Matrix(sp.BlockMatrix([[a, b], [c, d]]))


GAM = [blk(Z2, I2, I2, Z2)] + [blk(Z2, s, -s, Z2) for s in sig]
G5 = sp.I * GAM[0] * GAM[1] * GAM[2] * GAM[3]
BETA = GAM[0]


def ricci_scalar_const(T):
    """R for constant torsion, flat metric: Gamma = K, no derivatives."""
    Tl = lower1(ETA, T)
    Kl = contortion_correct(Tl)
    G = [[[ETA[l, l] * Kl[l][m][n] for n in R4] for m in R4] for l in R4]
    # R^r_{s m n} = d_m G^r_{n s} - d_n G^r_{m s} + G^r_{m l} G^l_{n s} - G^r_{n l} G^l_{m s}
    R = 0
    for s_, n in itertools.product(R4, R4):
        if ETA[s_, n] == 0:
            continue
        Ric = sum(G[r][r][l] * G[l][n][s_] - G[r][n][l] * G[l][r][s_] for r in R4 for l in R4)
        R += ETA[s_, n] * Ric
    return sp.expand(R)


def main():
    ok = True
    ok &= report("Clifford {g^a,g^b} = -2 eta^ab",
                 all(GAM[a] * GAM[b] + GAM[b] * GAM[a] == -2 * ETA[a, b] * sp.eye(4)
                     for a in R4 for b in R4))
    print("INFO: gamma_5 =", list(G5.diagonal()), "(diag; upper = left-handed)")
    r = rng(11)
    # 1. solve for c with random antisymmetric omega_{ab}
    c = sp.Symbol("c")
    om = sp.zeros(4)
    for a in R4:
        for b in range(a + 1, 4):
            v = sp.Rational(r.randint(-5, 5), r.randint(1, 4))
            om[a, b], om[b, a] = v, -v
    Om = c * sum((om[a, b] * GAM[a] * GAM[b] for a in R4 for b in R4), sp.zeros(4))
    eqs = []
    for a in R4:
        lhs = GAM[a] * Om - Om * GAM[a]
        rhs = sum((ETA[a, a] * om[a, b] * GAM[b] for b in R4), sp.zeros(4))  # omega^a_b
        eqs += list(lhs - rhs)
    sol = sp.solve([e for e in eqs if e != 0], c, dict=True)
    ok &= report("unique c for spinor spin connection", len(sol) == 1)
    cval = sol[0][c]
    print("INFO: Omega_m = c omega_{m ab} gamma^a gamma^b with c =", cval)
    Omv = Om.subs(c, cval)
    ok &= report("beta Omega^dag beta^-1 = -Omega", sp.simplify(BETA * Omv.H * BETA.inv() + Omv) == sp.zeros(4))

    # 2. axial coupling coefficient c_A with random full torsion
    T = random_torsion(r)
    Kl = contortion_correct(lower1(ETA, T))
    M = sp.zeros(4)
    for m in R4:
        Omm = cval * sum((Kl[a][m][b] * GAM[a] * GAM[b] for a in R4 for b in R4), sp.zeros(4))
        gm = GAM[m]  # flat tetrad: gamma^m = gamma^a delta
        M += sp.I / 2 * (gm * Omm + Omm * gm)
    S = axial_vector(ETA, T)
    Sl = [ETA[m, m] * S[m] for m in R4]
    cA = sp.Symbol("c_A")
    target = cA * sum((Sl[m] * GAM[m] * G5 for m in R4), sp.zeros(4))
    sA = sp.solve(list(M - target), cA, dict=True)
    ok &= report("torsion enters the Dirac action only via c_A S_m J5^m (trace/tensor drop out)",
                 len(sA) == 1)
    cAval = sA[0][cA]
    print("INFO: c_A =", cAval, " (L_D contains c_A S_m psibar gamma^m gamma_5 psi)")

    # 3. torsion-squared part of R
    t0, t1, s0, s1 = sp.symbols("t0 t1 s0 s1")
    Tt = torsion_from_parts(ETA, [t0, t1, 0, 0], [0, 0, 0, 0])
    Ts = torsion_from_parts(ETA, [0, 0, 0, 0], [s0, s1, 0, 0])
    ok &= report("decomposition: trace part has T_m = (t0,t1,0,0), S = 0",
                 [sp.simplify(v) for v in trace_vector(Tt)] == [t0, t1, 0, 0]
                 and all(sp.simplify(v) == 0 for v in axial_vector(ETA, Tt)))
    ok &= report("decomposition: axial part T_lmn = eps_lmnk S^k has S = (s0,s1,0,0), T_m = 0",
                 [sp.simplify(v) for v in axial_vector(ETA, Ts)] == [s0, s1, 0, 0]
                 and all(sp.simplify(v) == 0 for v in trace_vector(Ts)))
    Rt = ricci_scalar_const(Tt)
    Rs = ricci_scalar_const(Ts)
    TT = -t0**2 + t1**2
    SS = -s0**2 + s1**2
    bT = sp.simplify(Rt / TT)
    bS = sp.simplify(Rs / SS)
    ok &= report("R_Q(trace) = b_T T.T, R_Q(axial) = b_S S.S with constant b",
                 bT.is_number and bS.is_number)
    print(f"INFO: b_T = {bT}, b_S = {bS}")
    Tsum = [[[Tt[l][m][n] + Ts[l][m][n] for n in R4] for m in R4] for l in R4]
    ok &= report("no trace-axial cross term in R_Q", sp.expand(ricci_scalar_const(Tsum) - Rt - Rs) == 0)
    # tensor part: q = T_random - trace part - axial part; check no cross terms with q
    Tr = random_torsion(r)
    tv = trace_vector(Tr)
    sv = axial_vector(ETA, Tr)
    base = torsion_from_parts(ETA, tv, sv)
    q = [[[Tr[l][m][n] - base[l][m][n] for n in R4] for m in R4] for l in R4]
    ok &= report("q has zero trace and zero axial part",
                 all(sp.simplify(v) == 0 for v in trace_vector(q) + axial_vector(ETA, q)))
    Tts = torsion_from_parts(ETA, [t0, t1, 0, 0], [s0, s1, 0, 0])
    Tall = [[[Tts[l][m][n] + q[l][m][n] for n in R4] for m in R4] for l in R4]
    cross = sp.expand(ricci_scalar_const(Tall) - ricci_scalar_const(Tts) - ricci_scalar_const(q))
    ok &= report("no cross terms between tensor part q and (trace, axial)", cross == 0)

    # 4. eliminate algebraic S
    kap = sp.Symbol("kappa", positive=True)
    Sx = sp.symbols("S0:4")
    J = sp.symbols("J0:4")
    L = bS / (2 * kap) * sum(ETA[m, m] * Sx[m]**2 for m in R4) + cAval * sum(Sx[m] * ETA[m, m] * J[m] for m in R4)
    solS = sp.solve([sp.diff(L, v) for v in Sx], Sx, dict=True)[0]
    Leff = sp.simplify(L.subs(solS))
    JJ = sum(ETA[m, m] * J[m]**2 for m in R4)
    print("INFO: S^m =", [solS[v] for v in Sx], " (J are J5^m upper components)")
    print("INFO: L_eff =", sp.factor(Leff), " = ", sp.simplify(Leff / JJ), "* J5.J5")
    ok &= report("L_eff = -(kappa c_A^2/(2 b_S)) J5.J5",
                 sp.simplify(Leff + kap * cAval**2 / (2 * bS) * JJ) == 0)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
