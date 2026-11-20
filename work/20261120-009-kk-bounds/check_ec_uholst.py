"""Checks E1-E7: Einstein-Cartan analogue of the 2608.09299 term.

Pointwise algebra in an orthonormal 5D frame (eta = diag(-1,1,1,1,1), index 4 = 5^).
Contortion 1-form K_AB = K[A,B,C] e^C, K[A,B,C] = -K[B,A,C] (50 components).
Torsion T^A = K^A_B ^ e^B (Levi-Civita torsion is zero). R(w+K) = R_LC + D_LC K + K^K.
u is a unit 1-form (u = e^4 here, the zero-mode value of the 2608 Stueckelberg vector),
W = du is a background 2-form (for u = e^5 of the KK vielbein, W_ab = e^{beta phi} k F_ab).

Terms (kappa = 1):
  L_EC   = (1/2) Q(K)                                   [K-part of R(w), plus a total derivative]
  mixing = a * *( du ^ e^A ^ e^B ^ K_AB )               [from u ^ e^A^e^B ^ D_LC K_AB, after parts]
  odd    = b * *( u ^ T^A ^ T_A )                        [= b * *(u ^ e^A^e^B^(K^K)_AB), checked]
u-Holst  c * *(u ^ e^A ^ e^B ^ R_AB(w)) has a = b = c (the R_LC part vanishes by Bianchi).
u ^ Nieh-Yan = u ^ d(e^A ^ T_A) has a = -c, b = 0.
2608 transplanted, c * u ^ T^A ^ T_A with EC torsion, has a = 0, b = c.
Integrating out K gives L_eff(W) = A W.W + B_star *(u ^ W ^ W).
"""
import itertools
import math

import numpy as np

from common import check, info

D = 5
ETA = np.diag([-1.0, 1, 1, 1, 1])
rng = np.random.default_rng(20261120)


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


_PERMS = {}


def alt(t):
    """Antisymmetrize all indices (average over permutations)."""
    r = t.ndim
    if r <= 1:
        return t
    if r not in _PERMS:
        _PERMS[r] = [(p, perm_sign(p)) for p in itertools.permutations(range(r))]
    out = np.zeros_like(t)
    for p, s in _PERMS[r]:
        out += s * np.transpose(t, p)
    return out / math.factorial(r)


def wedge(a, b):
    p, q = a.ndim, b.ndim
    return math.comb(p + q, p) * alt(np.multiply.outer(a, b))


def top(w):
    return w[0, 1, 2, 3, 4]


E1F = np.eye(D)  # e^A as 1-forms: E1F[A] has components delta^A_C


def kform(K):
    """K^A_B 1-forms: array [A,B,C] with first index raised."""
    return np.einsum("aa,abc->abc", ETA, K)


def torsion(K):
    Ku = kform(K)
    return np.array([sum(wedge(Ku[A, B], E1F[B]) for B in range(D)) for A in range(D)])


def kk2(K):
    """(K^K)_AB 2-forms, lowered A: K_AC ^ K^C_B."""
    Ku = kform(K)
    return np.array([[sum(wedge(K[A, C], Ku[C, B]) for C in range(D)) for B in range(D)] for A in range(D)])


def Qscalar(K):
    """K-quadratic part of the scalar curvature: eta^{BD} (K^K)^A_{B,AD}."""
    KK = np.einsum("aa,abcd->abcd", ETA, kk2(K))  # raise first index
    return sum(ETA[B, B] * KK[A, B, A, B] for A in range(D) for B in range(D))


BASIS = []
for A in range(D):
    for B in range(A + 1, D):
        for C in range(D):
            k = np.zeros((D, D, D))
            k[A, B, C], k[B, A, C] = 1.0, -1.0
            BASIS.append(k)
NB = len(BASIS)


def ee():
    return np.array([[wedge(E1F[A], E1F[B]) for B in range(D)] for A in range(D)])


EE = ee()
U = E1F[4]


def L_parts(K, W):
    """(Q, mixing, odd) scalars for contortion K and background W = du."""
    T = torsion(K)
    Tl = np.einsum("ab,bcd->acd", ETA, T)
    TT = sum(wedge(T[A], Tl[A]) for A in range(D))
    odd = top(wedge(U, TT))
    eeK = sum(wedge(EE[A, B], K[A, B]) for A in range(D) for B in range(D))
    mix = top(wedge(W, eeK))
    return Qscalar(K), mix, odd


def levi(n):
    e = np.zeros((n,) * n)
    for p in itertools.permutations(range(n)):
        e[p] = perm_sign(p)
    return e


EPS5 = levi(5)
EPS4 = levi(4)
BST = np.array(BASIS)  # (50, 5, 5, 5)


def Q_bilinear(K1, K2):
    """Symmetric bilinear form with Q(K) = Q_bilinear(K, K); K stacks (n,5,5,5)."""
    K2u = np.einsum("cc,ncbe->ncbe", ETA, K2)   # K^C_{B E}
    # (K_AC ^ K^C_B)_{DE} = K_ACD K^C_BE - K_ACE K^C_BD ; Q = eta^AA eta^BD (.)_{A D}
    t1 = np.einsum("aa,bb,maca,ncbb->mn", ETA, ETA, K1, K2u)
    t2 = np.einsum("aa,bb,macb,ncba->mn", ETA, ETA, K1, K2u)
    M = t1 - t2
    return 0.5 * (M + M.T) if K1 is K2 else M


def odd_bilinear(K1, K2):
    """u ^ T^A ^ T_A with u = e^4: (T^A ^ T_A)_{0123} = 1/4 eps~^{abcd} T^A_ab T_A cd."""
    def T(K):  # T^A_{CD} = K^A_{B C} delta^B_D - K^A_{B D} delta^B_C = K^A_{DC} - K^A_{CD}
        Ku = np.einsum("aa,nabc->nabc", ETA, K)
        return np.transpose(Ku, (0, 1, 3, 2)) - Ku
    T1, T2 = T(K1), T(K2)
    T2l = np.einsum("aa,nacd->nacd", ETA, T2)
    M = 0.25 * np.einsum("abcd,mxab,nxcd->mn", EPS4, T1[:, :, :4, :4], T2l[:, :, :4, :4])
    return 0.5 * (M + M.T)


def mix_linear(Wstack):
    """*(W ^ e^A^e^B^K_AB) = (1/12) eps~^{abcde} W_ab (6 alt K)_cde, for each W and basis K."""
    Y = np.array([6 * alt(k) for k in BASIS])
    return np.einsum("abcde,wab,ncde->wn", EPS5, Wstack, Y) / 12.0


def main():
    # E1: identity T^A ^ T_A = e^A ^ e^B ^ (K^K)_AB for random K
    K = sum(rng.normal() * b for b in BASIS)
    T = torsion(K)
    Tl = np.einsum("ab,bcd->acd", ETA, T)
    TT = sum(wedge(T[A], Tl[A]) for A in range(D))
    eeKK = sum(wedge(EE[A, B], kk2(K)[A, B]) for A in range(D) for B in range(D))
    check("E1 T^A^T_A = e^A^e^B^(K^K)_AB (so u^Nieh-Yan has no K^2 piece)", np.allclose(TT, eeKK, atol=1e-12))

    # E2: Q for totally antisymmetric K equals -K_ABC K^ABC (matches task 004)
    Ka = alt(K)
    Kup = np.einsum("aa,bb,cc,abc->abc", ETA, ETA, ETA, Ka)
    check("E2 Q(K_[ABC]) = -K_ABC K^ABC (same normalization as task 004)",
          abs(Qscalar(Ka) + np.sum(Ka * Kup)) < 1e-10)

    # E3: algebraic Bianchi: e^A^e^B^R_AB = 0 for a torsion-free Riemann tensor (2608 P1 vanishes)
    Rm = rng.normal(size=(D, D, D, D))
    Rm = Rm - Rm.transpose(1, 0, 2, 3)
    Rm = Rm - Rm.transpose(0, 1, 3, 2)
    Rm = Rm + Rm.transpose(2, 3, 0, 1)
    Rm = Rm - alt(Rm) * 1.0  # remove totally antisymmetric part -> first Bianchi holds
    Rm = 0.5 * (Rm + Rm.transpose(2, 3, 0, 1))
    eeR = sum(wedge(EE[A, B], Rm[A, B]) for A in range(D) for B in range(D))
    check("E3 e^A^e^B^R_AB = 0 for Levi-Civita curvature (no EC term without torsion)",
          np.allclose(eeR, 0, atol=1e-12))

    # matrices (einsum forms), cross-checked against the generic wedge code
    MQ = 2 * Q_bilinear(BST, BST)       # Q(K) = 1/2 k.MQ.k
    MO = 2 * odd_bilinear(BST, BST)     # odd(K) = 1/2 k.MO.k
    pairs = [(a, b) for a in range(4) for b in range(a + 1, 4)]
    Wb = []
    for a_, b_ in pairs:
        w = np.zeros((D, D))
        w[a_, b_], w[b_, a_] = 1.0, -1.0
        Wb.append(w)
    Wb = np.array(Wb)
    V = mix_linear(Wb)                   # mix(K) = V.k
    kv = rng.normal(size=NB)
    Kr = np.einsum("n,nabc->abc", kv, BST)
    Wr = np.einsum("w,wab->ab", rng.normal(size=6), Wb)
    q_g, m_g, o_g = L_parts(Kr, Wr)
    wv = np.linalg.lstsq(Wb.reshape(6, -1).T, Wr.ravel(), rcond=None)[0]
    ok = (abs(0.5 * kv @ MQ @ kv - q_g) < 1e-10 and abs(0.5 * kv @ MO @ kv - o_g) < 1e-10
          and abs(wv @ V @ kv - m_g) < 1e-10)
    check("E4a einsum forms agree with generic wedge algebra at a random point", ok)
    check("E4 Q non-degenerate on 50 contortion components", np.linalg.matrix_rank(MQ, tol=1e-9) == NB)

    WW = np.einsum("iab,jcd,ac,bd->ij", Wb[:, :4, :4], Wb[:, :4, :4], ETA[:4, :4], ETA[:4, :4])
    WtW = 0.25 * np.einsum("abcd,iab,jcd->ij", EPS4, Wb[:, :4, :4], Wb[:, :4, :4])
    WtW = 0.5 * (WtW + WtW.T)

    def leff(a, b):
        """G with L_eff(W) = w.G.w; L_K = 1/2 Q + b odd + a mix (kappa = 1)."""
        Mt = 0.5 * MQ + b * MO
        lin = a * V
        G = -0.5 * lin @ np.linalg.solve(Mt, lin.T)
        return 0.5 * (G + G.T)

    def fit(G):
        X = np.stack([WW.ravel(), WtW.ravel()], axis=1)
        coef, res, *_ = np.linalg.lstsq(X, G.ravel(), rcond=None)
        resid = np.linalg.norm(X @ coef - G.ravel()) / max(np.linalg.norm(G), 1e-300)
        return coef, resid

    # E5: u ^ Nieh-Yan (a = -c, b = 0): parity-even only
    (A_ny, B_ny), r_ny = fit(leff(-0.3, 0.0))
    check("E5 EC u^Nieh-Yan: L_eff = A W.W only (no W Wtilde)", abs(B_ny) < 1e-12 and r_ny < 1e-10,
          f"A={A_ny:.4g}, B*={B_ny:.2e}")
    # E6: 2608 transplanted (a = 0): nothing
    G0 = leff(0.0, 0.3)
    check("E6 EC u^T^A^T_A alone (2608 transplanted): no photon term", np.allclose(G0, 0, atol=1e-14))

    # E7: u-Holst (a = b = c): parity-odd W Wtilde appears at O(c^3)
    rows = []
    for c in [1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0]:
        (Ac, Bc), r = fit(leff(c, c))
        Z = 1 - 8 * Ac
        codd = -4 * Bc / Z
        rows.append((c, Ac, Bc, r, codd))
        info(f"u-Holst c={c:g}: A={Ac:.6g}  B*={Bc:.6g}  fit resid={r:.1e}  c_odd={codd:.6g}")
    small = [r for r in rows if r[0] <= 1e-2]
    p_B = np.polyfit(np.log([r[0] for r in small]), np.log([abs(r[2]) for r in small]), 1)[0]
    p_A = np.polyfit(np.log([r[0] for r in small]), np.log([abs(r[1]) for r in small]), 1)[0]
    check("E7 u-Holst: fits close (W.W, W Wtilde only) and B* ~ c^3, A ~ c^2 at small c",
          all(r[3] < 1e-8 for r in rows) and abs(p_B - 3) < 0.02 and abs(p_A - 2) < 0.02,
          f"slopes A {p_A:.3f}, B {p_B:.3f}")
    check("E8 u-Holst exact: A = -3c^2/2, B* = -4c^3, c_odd = 16c^3/(1+12c^2) for all scanned c",
          all(abs(r[1] + 1.5 * r[0] ** 2) < 1e-9 * max(1, r[0] ** 2) and abs(r[2] + 4 * r[0] ** 3) < 1e-9 * max(1, r[0] ** 3)
              and abs(r[4] - 16 * r[0] ** 3 / (1 + 12 * r[0] ** 2)) < 1e-9 * max(1, r[0]) for r in rows))
    cmax = max(abs(r[4]) for r in rows)
    info(f"max |c_odd| over scanned c: {cmax:.4g}")
    return rows


if __name__ == "__main__":
    main()
