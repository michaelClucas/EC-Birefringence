"""Checks D1-D9: fermion dipole couplings of the KK-05 candidate mechanism.

5D Einstein-Cartan, orthonormal frame, eta = diag(-1,1,1,1,1), index 4 = 5^ (fifth direction).
KK vielbein e^4 = dy + k A (radion stabilised, phi = 1); W = de^4 = (1/2) W_ab e^a e^b, W = k F.
Action (kappa_5 = 1):  S = int [ (1/2) R(w) + c *(u ^ e^A ^ e^B ^ R_AB(w)) ],  u = e^4 (oriented: P'-odd).
w = w0 (Levi-Civita) + K (contortion 1-form K_AB = K[A,B,C] e^C).
K-dependent part (after the u-Holst D0 K term is integrated by parts, D0 e = 0):
  L_K = (1/2) Q(K) + c odd(K) + mix(K; X),   X = d(c u) = c W + c' dchi ^ u,
  Q   = scalar from K^K, odd = *(e^4 ^ e^A ^ e^B ^ (K^K)_AB), mix = *(X ^ e^A ^ e^B ^ K_AB).
Fermion: L = i Psibar Gamma^C D_C Psi - M Psibar Psi, {Gamma^A, Gamma^B} = -2 eta^AB, Psibar = Psi^dag gamma^0,
D_C = d_C + (1/2) w_ABC Sigma^AB, hermitian part kept. Zero mode: Psi(x, y) = psi(x)/sqrt(2 pi R).

This is an independent implementation (exact sympy linear algebra for K, explicit gamma matrices for the
fermion); it reproduces the 009/012 photon result (Z = 1 + 12c^2, Y = 16c^3) as a normalisation check and
then adds the dipole couplings, which 012 only estimated.
"""
import itertools

import numpy as np
import sympy as sp

from common import check, info

D = 5
ETA = [-1, 1, 1, 1, 1]


def perm_sign(p):
    p, s = list(p), 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


EPS5 = {p: perm_sign(p) for p in itertools.permutations(range(5))}
EPS4 = {p: perm_sign(p) for p in itertools.permutations(range(4))}

# ---------------------------------------------------------------- contortion sector (sympy, exact)
IDX = [(A, B, C) for A in range(D) for B in range(A + 1, D) for C in range(D)]
NB = len(IDX)
KS = sp.symbols(f"k0:{NB}")


def Kten(kv):
    K = [[[0] * D for _ in range(D)] for _ in range(D)]
    for (A, B, C), v in zip(IDX, kv):
        K[A][B][C] = v
        K[B][A][C] = -v
    return K


def riemann_KK(K):
    """R^A_{B C D} from w = K with constant components: K^A_{E C} K^E_{B D} - (C <-> D)."""
    Ku = [[[ETA[A] * K[A][B][C] for C in range(D)] for B in range(D)] for A in range(D)]
    R = {}
    for A in range(D):
        for B in range(D):
            for C in range(D):
                for Dd in range(D):
                    R[A, B, C, Dd] = sum(Ku[A][E][C] * Ku[E][B][Dd] - Ku[A][E][Dd] * Ku[E][B][C] for E in range(D))
    return R


def L_K(kv, c, X):
    K = Kten(kv)
    R = riemann_KK(K)
    Q = sum(ETA[B] * R[A, B, A, B] for A in range(D) for B in range(D))
    odd = sp.Rational(1, 2) * sum(s * ETA[p[0]] * R[p[0], p[1], p[2], p[3]] for p, s in EPS4.items())
    mix = sp.Rational(1, 2) * sum(s * X[p[0]][p[1]] * K[p[2]][p[3]][p[4]] for p, s in EPS5.items())
    return sp.Rational(1, 2) * Q + c * odd + mix


def solve_K(c, X):
    L = sp.expand(L_K(KS, c, X))
    M = sp.hessian(L, KS)
    v = sp.Matrix([sp.diff(L, k).subs({kk: 0 for kk in KS}) for k in KS])
    sol = M.LUsolve(-v)
    Leff = sp.Rational(1, 2) * (v.T * sol)[0]
    return [sp.nsimplify(s) for s in sol], sp.simplify(Leff), M


def two_form(pairs):
    X = [[0] * D for _ in range(D)]
    for (a, b), val in pairs.items():
        X[a][b] = val
        X[b][a] = -val
    return X


def alt3(K, A, B, C):
    return sum(s * K[[A, B, C][p[0]]][[A, B, C][p[1]]][[A, B, C][p[2]]]
               for p, s in {q: perm_sign(q) for q in itertools.permutations(range(3))}.items()) / 6


# ---------------------------------------------------------------- spin connection of the KK vielbein
def lc_connection(W):
    """w0_{ABC} for e^4 = dy + kA with constant W: from de^A = -w^A_B ^ e^B."""
    w = np.zeros((D, D, D))
    for b in range(4):
        for cc in range(4):
            w[4, b, cc] = 0.5 * W[b][cc]
            w[b, 4, cc] = -0.5 * W[b][cc]
    for a in range(4):
        for b in range(4):
            w[a, b, 4] = -0.5 * W[a][b]
    return w


def de_minus_wedge(w, W):
    """Residual of de^A + w^A_B ^ e^B = 0 as 2-form components [A, C, D]."""
    res = np.zeros((D, D, D))
    for C in range(D):
        for Dd in range(D):
            res[4, C, Dd] += W[C][Dd] if C < 4 and Dd < 4 else 0.0
    for A in range(D):
        for B in range(D):
            wu = ETA[A] * w[A, B]           # w^A_B components (1-form)
            for C in range(D):
                for Dd in range(D):
                    res[A, C, Dd] += wu[C] * (1.0 if B == Dd else 0.0) - wu[Dd] * (1.0 if B == C else 0.0)
    return res


def ricci_scalar_lc(w, W):
    """R0 = eta^{BD} R^A_{BAD}, R = dw + w^w, dw^A_B = w^A_{B4} W (constant frame components)."""
    wu = np.einsum("a,abc->abc", np.array(ETA, float), w)
    R = np.zeros((D, D, D, D))
    Wf = np.zeros((D, D))
    Wf[:4, :4] = np.array(W, float)[:4, :4]
    for A in range(D):
        for B in range(D):
            R[A, B] += wu[A, B, 4] * Wf
            R[A, B] += np.einsum("ec,ed->cd", wu[A, :, :], wu[:, B, :]) - np.einsum("ed,ec->cd", wu[A, :, :], wu[:, B, :])
    return sum(ETA[B] * R[A, B, A, B] for A in range(D) for B in range(D))


# ---------------------------------------------------------------- gamma matrices (Dirac rep, Srednicki signs)
s0 = np.eye(2)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
PAULI = [sx, sy, sz]
Z2 = np.zeros((2, 2))
g0 = np.block([[s0, Z2], [Z2, -s0]]).astype(complex)
gi = [np.block([[Z2, s], [-s, Z2]]) for s in PAULI]
g5 = 1j * g0 @ gi[0] @ gi[1] @ gi[2]
GU = [g0] + gi + [1j * g5]                         # Gamma^A, A = 0..4
GL = [ETA[A] * GU[A] for A in range(D)]           # Gamma_A


def sigma_spin(sgn):
    return {(A, B): sgn * 0.25 * (GU[A] @ GU[B] - GU[B] @ GU[A]) for A in range(D) for B in range(D)}


def herm(X):
    return 0.5 * (X + g0 @ X.conj().T @ g0)


def spin_term(w, SIG):
    """Hermitian part of i Gamma^C (1/2) w_ABC Sigma^AB."""
    M = np.zeros((4, 4), complex)
    for C in range(D):
        Om = sum(0.5 * w[A, B, C] * SIG[A, B] for A in range(D) for B in range(D))
        M += 1j * GU[C] @ Om
    return herm(M)


def rest_frame(H):
    blk = (g0 @ H)[:2, :2]
    return np.array([0.5 * np.trace(s @ blk).real for s in PAULI]), blk


def main():
    rng = np.random.default_rng(16)

    # D0: spinor generator sign so that Gamma is covariantly constant
    wr = rng.normal(size=(D, D))
    wr = wr - wr.T
    ok_sgn = None
    for sgn in (1, -1):
        SIG = sigma_spin(sgn)
        Om = sum(0.5 * wr[A, B] * SIG[A, B] for A in range(D) for B in range(D))
        good = all(np.allclose(Om @ GU[A] - GU[A] @ Om, -sum(ETA[A] * wr[A, B] * GU[B] for B in range(D))) for A in range(D))
        if good:
            ok_sgn = sgn
    check("D0 Clifford {G,G} = -2 eta and a spinor generator sign making Gamma covariantly constant exist",
          ok_sgn is not None and all(np.allclose(GU[A] @ GU[B] + GU[B] @ GU[A], -2 * ETA[A] * (A == B) * np.eye(4))
                                     for A in range(D) for B in range(D)), f"sign {ok_sgn}")
    SIG = sigma_spin(ok_sgn)

    # D1: Levi-Civita connection of the KK vielbein, totally antisymmetric part = W/6
    Wn = rng.normal(size=(4, 4))
    Wn = Wn - Wn.T
    Wfull = np.zeros((D, D))
    Wfull[:4, :4] = Wn
    w0 = lc_connection(Wfull)
    res = de_minus_wedge(w0, Wfull)
    a0 = np.array([[alt3(w0, a, b, 4) for b in range(4)] for a in range(4)])
    check("D1 KK Levi-Civita connection solves de = -w^e; w0_[ab5] = W_ab/6 (task 004 value)",
          np.allclose(res, 0) and np.allclose(a0, Wn / 6), f"max resid {abs(res).max():.1e}")

    # D2: Maxwell normalisation from (1/2) R0
    R0 = ricci_scalar_lc(w0, Wfull)
    WW = sum(ETA[a] * ETA[b] * Wn[a, b] ** 2 for a in range(4) for b in range(4))
    check("D2 (1/2) R0 = -(1/8) W.W (Maxwell term of the KK photon, healthy sign)", abs(0.5 * R0 + WW / 8) < 1e-12,
          f"{0.5 * R0:.6f} vs {-WW / 8:.6f}")

    # D3/D4: exact contortion elimination for rational c, generic rational W
    Wq = two_form({(0, 1): sp.Rational(3, 7), (0, 2): sp.Rational(-2, 5), (0, 3): sp.Rational(1, 3),
                   (1, 2): sp.Rational(5, 11), (1, 3): sp.Rational(-1, 2), (2, 3): sp.Rational(2, 9)})
    WWq = sum(ETA[a] * ETA[b] * Wq[a][b] ** 2 for a in range(4) for b in range(4))
    WWt = sp.Rational(1, 4) * sum(s * Wq[p[0]][p[1]] * Wq[p[2]][p[3]] for p, s in EPS4.items())
    Wtil = [[sp.Rational(1, 2) * sum(EPS4.get((a, b, cc, d), 0) * ETA[cc] * ETA[d] * Wq[cc][d]
                                     for cc in range(4) for d in range(4)) for b in range(4)] for a in range(4)]
    ok3, ok4, ok4b = True, True, True
    Kdat = {}
    for c in (sp.Rational(1, 3), sp.Integer(2), sp.Integer(7)):
        X = [[c * Wq[a][b] for b in range(D)] for a in range(D)]
        sol, Leff, M = solve_K(c, X)
        ok3 &= sp.simplify(Leff - (-sp.Rational(3, 2) * c ** 2 * WWq - 4 * c ** 3 * WWt)) == 0
        K = Kten(sol)
        Kab4 = [[sp.nsimplify(alt3(K, a, b, 4)) for b in range(4)] for a in range(4)]
        fits = [all(sp.simplify(Kab4[a][b] - (-sp.Rational(4, 3) * c ** 2 * Wq[a][b] + sg * c * Wtil[a][b])) == 0
                    for a in range(4) for b in range(4)) for sg in (1, -1)]
        ok4 &= any(fits)
        # no other totally antisymmetric component sourced by W
        others = [alt3(K, A, B, C) for A, B, C in itertools.combinations(range(4), 3)]
        ok4b &= all(sp.simplify(o) == 0 for o in others)
        Kdat[c] = (fits, M)
        info(f"c={c}: L_eff matches; K_[ab5] = -(4/3)c^2 W + ({'+' if fits[0] else '-'})c W~")
    check("D3 independent exact elimination: L_eff = -(3c^2/2) W.W - 4c^3 *(u^W^W)  => Z = 1+12c^2, Y = 16c^3 (009/012)", ok3)
    check("D4 K_[ab5] = -(4/3)c^2 W_ab +- c W~_ab exactly (012 U3); K_[abc] (4D indices) = 0", ok4 and ok4b)

    # D5: fermion dipole couplings in the rest frame (numerical, float), as functions of c
    # numeric K response: use the exact formula verified in D4 (sign of the W~ term from D4)
    sg = 1 if Kdat[sp.Integer(2)][0][0] else -1

    def total_w(c, W):
        Wt = np.zeros((4, 4))
        for a in range(4):
            for b in range(4):
                Wt[a, b] = 0.5 * sum(EPS4.get((a, b, cc, d), 0) * ETA[cc] * ETA[d] * W[cc, d] for cc in range(4) for d in range(4))
        Wf = np.zeros((D, D))
        Wf[:4, :4] = W
        w = lc_connection(Wf)
        Ka = -4.0 / 3.0 * c ** 2 * W + sg * c * Wt          # totally antisymmetric K_[ab4]
        Kt = np.zeros((D, D, D))
        for a in range(4):
            for b in range(4):
                for (i, j, k) in itertools.permutations((0, 1, 2)):
                    idx = [a, b, 4]
                    Kt[idx[i], idx[j], idx[k]] += perm_sign((i, j, k)) * Ka[a, b] / 2.0   # sum over perms of alt part
        # Kt is totally antisymmetric with Kt[a,b,4] = Ka[a,b] (each ordered pair counted once per orientation)
        return w + Kt

    def dipoles(c):
        E1 = np.zeros((4, 4))
        E1[0, 1], E1[1, 0] = -1.0, 1.0          # W_01 = -E_x (E_x = 1)
        B1 = np.zeros((4, 4))
        B1[2, 3], B1[3, 2] = 1.0, -1.0          # W_23 = B_x = 1
        vE, _ = rest_frame(spin_term(total_w(c, E1), SIG))
        vB, _ = rest_frame(spin_term(total_w(c, B1), SIG))
        return vE, vB

    # sanity: the constructed Kt has the intended alt part
    Wtest = Wn
    wt = total_w(0.7, Wtest) - lc_connection(Wfull)
    alt_ok = np.allclose([[alt3(wt, a, b, 4) for b in range(4)] for a in range(4)],
                         -4 / 3 * 0.49 * Wtest + sg * 0.7 * np.array([[0.5 * sum(EPS4.get((a, b, cc, d), 0) * ETA[cc] * ETA[d] * Wtest[cc, d]
                                                                                      for cc in range(4) for d in range(4)) for b in range(4)] for a in range(4)]))
    rows = []
    for c in (0.0, 0.1, 0.3, 1.0, 10.0, 100.0, 1000.0):
        vE, vB = dipoles(c)
        rows.append((c, vE, vB))
        info(f"c={c:g}: rest-frame L = phi^dag (D_E sigma.E + D_B sigma.B) phi per unit W: "
             f"D_E={vE[0]:+.6g} (sigma_x), D_B={vB[0]:+.6g}; off-axis {max(abs(vE[1:]).max(), abs(vB[1:]).max()):.1e}")
    dE = np.array([abs(r[1][0]) for r in rows])
    dB = np.array([abs(r[2][0]) for r in rows])
    cs = np.array([r[0] for r in rows])
    check("D5 c = 0: pure KK Pauli term is an EDM of size W/4 (D_E = 1/4, D_B = 0) for a real-mass zero mode",
          abs(dE[0] - 0.25) < 1e-12 and dB[0] < 1e-12 and alt_ok)
    check("D6 D_E(c) = (3/2)|1/6 - (4/3)c^2| (P'-even part -> EDM), D_B(c) = (3/2)c (P'-odd part -> MDM)",
          np.allclose(dE, 1.5 * abs(1 / 6 - 4 / 3 * cs ** 2), rtol=1e-10) and np.allclose(dB, 1.5 * cs, rtol=1e-10, atol=1e-12),
          "012 assigned the c-odd piece to the EDM; it is the magnetic one")

    # D7: hermitian projection removes the non-totally-antisymmetric parts; chirality
    wr3 = rng.normal(size=(D, D, D))
    wr3 = wr3 - wr3.transpose(1, 0, 2)
    walt = np.zeros_like(wr3)
    for p in itertools.permutations(range(3)):
        walt += perm_sign(p) * wr3.transpose(p)
    walt /= 6
    check("D7a only w_[ABC] survives the hermitian projection (vector parts drop)",
          np.allclose(spin_term(wr3, SIG), spin_term(walt, SIG)))
    PL = 0.5 * (np.eye(4) - g5)
    H = spin_term(total_w(30.0, Wn), SIG)
    check("D7b a single-chirality zero mode has no Pauli coupling: P_L^dag gamma0 H_dipole P_L = 0",
          np.allclose(PL.conj().T @ g0 @ H @ PL, 0))

    # D8: rolling chi, source X = c' chi-dot e^0 ^ e^4 (W = 0): background contortion and fermion term
    BASIS16 = {"1": np.eye(4), "g5": g5}
    for m in range(4):
        BASIS16[f"g{m}"] = GU[m]
        BASIS16[f"g{m}g5"] = GU[m] @ g5
    for m, n in itertools.combinations(range(4), 2):
        BASIS16[f"s{m}{n}"] = GU[m] @ GU[n]
    Xq = two_form({(0, 4): sp.Integer(1)})
    ok8 = True
    for c in (sp.Integer(0), sp.Integer(1), sp.Integer(30)):
        sol, Leff, _ = solve_K(c, Xq)
        K = np.array(Kten([float(s) for s in sol]), float)
        H = spin_term(K, SIG)
        coefs = {n: np.trace(np.linalg.inv(Bm) @ g0 @ H).real / 4 + 1j * np.trace(np.linalg.inv(Bm) @ g0 @ H).imag / 4
                 for n, Bm in BASIS16.items()}
        nz = {n: v for n, v in coefs.items() if abs(v) > 1e-12}
        K123 = alt3(K, 1, 2, 3)
        info(f"c={c}: chi-dot source -> K_[123] = {K123:+.6g} (units c' chi-dot); L_eff = {Leff}; "
             f"Psibar(...)Psi = Psi^dag(...)Psi with nonzero {nz}")
        ok8 &= (set(nz) == {"g5"} and abs(abs(K123) - 1) < 1e-12 and abs(abs(list(nz.values())[0]) - 1.5) < 1e-12
                and sp.simplify(Leff - 3) == 0)
    check("D8 rolling chi: K_[123] = -+c' chi-dot (c-independent), L_eff = 3 c'^2 chi-dot^2 (= 012's -3c'^2 (dchi)^2), "
          "fermion term 1.5 c' chi-dot Psi^dag gamma5 Psi = axial b_0 = 1.5 c' chi-dot (psibar gamma^0 gamma5 psi); no dipole", ok8)
    return rows


if __name__ == "__main__":
    main()
