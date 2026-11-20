"""Question A, variant (a): dynamical Barbero-Immirzi field in the Holst term with minimally
coupled fermions.  Does the scalar keep a shift symmetry?

Pointwise algebra in a 4D orthonormal frame, eta = diag(-1,1,1,1), eps^{0123} = -1.
Contortion K[a,b,c] = -K[b,a,c] (24 components), R^{ab}(w+K) = R_LC + D_LC K + K^K.
  L = (1/2kappa)[ Q(K) + gb * H(K) ] - (1/2kappa) eps^{abcd} v_c K_{abd} + lam eps^{abcd} K_{abc} J_d
with gb = 1/gamma (Barbero-Immirzi field), v = d(gb),
  Q = (K^K)^{ab}_{ab} (scalar curvature), H = (1/2) eps^{abcd} (K^K)_{abcd} (Holst density),
the linear term from gb * (1/2) eps (D K) after integration by parts, and the minimal Dirac coupling
to the totally antisymmetric contortion (lam fixed by the repo result -(3 kappa/16) J5^2 at gb = 0).
Eliminating K gives L_eff = P v.v + X v.J + Y J.J.
H1  lam^2 fixed; Y(gb=0) = -3/16 (kappa = 1)
H2  Y(gb) = -(3/16) gamma^2/(1+gamma^2)   [Perez-Rovelli gr-qc/0505081 Eq. (15)]
H3  only the three invariants appear (fit residual)
H4  P(gb) proportional to 1/(1+gb^2)  [Taveras-Yunes 0807.2652 structure]; print P, X
H5  in the canonical field Y depends on phi: shift symmetry broken by a parity-even J5^2 term
H6  the phi-dependence of Y and P is even in gb, so no gb-odd (theta-bar) coupling arises at tree level
"""
import itertools

import numpy as np

from common import check, info

D = 4
ETA = np.diag([-1.0, 1, 1, 1])
rng = np.random.default_rng(12)


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


EPS_LO = np.zeros((4,) * 4)
for p in itertools.permutations(range(4)):
    EPS_LO[p] = perm_sign(p)          # eps_{0123} = +1
EPS_UP = -EPS_LO                       # eps^{0123} = -1 (orthonormal frame)

BASIS = []
for A in range(D):
    for B in range(A + 1, D):
        for C in range(D):
            k = np.zeros((D, D, D))
            k[A, B, C], k[B, A, C] = 1.0, -1.0
            BASIS.append(k)
BST = np.array(BASIS)
NB = len(BST)


def KK(K1, K2):
    """(K1^K2)_{ab,cd} = K1_{ae c} eta^{ee} K2_{e b d} - (c<->d), all indices down."""
    t = np.einsum("aec,ee,ebd->abcd", K1, ETA, K2)
    return t - t.transpose(0, 1, 3, 2)


def Qf(K1, K2):
    return np.einsum("aa,bb,abab->", ETA, ETA, KK(K1, K2))


def Hf(K1, K2):
    return 0.5 * np.einsum("abcd,abcd->", EPS_UP, KK(K1, K2))


MQ = np.array([[Qf(a, b) + Qf(b, a) for b in BST] for a in BST])   # Q(k) = 1/2 k.MQ.k
MH = np.array([[Hf(a, b) + Hf(b, a) for b in BST] for a in BST])
LV = np.array([[-np.einsum("abd,abd->", EPS_UP[:, :, c, :], k) for k in BST] for c in range(D)])  # v_c coeff
LJ = np.array([[np.einsum("abc,abc->", EPS_UP[:, :, :, d], k) for k in BST] for d in range(D)])   # J_d coeff


def leff_matrix(gb, lam):
    """8x8 G with L_eff = s.G.s, s = (v_0..3, J_0..3); kappa = 1."""
    M = 0.5 * (MQ + gb * MH)                 # (1/2kappa) * (1/2) k.(MQ+gb MH).k = 1/2 k.M.k
    lin = np.vstack([0.5 * LV, lam * LJ])    # b = lin^T s
    G = -0.5 * lin @ np.linalg.solve(M, lin.T)
    return 0.5 * (G + G.T)


def invariants(G):
    """Fit G blocks to P eta, X/2 eta, Y eta; return P, X, Y, residual."""
    Gvv, Gvj, Gjj = G[:4, :4], G[:4, 4:], G[4:, 4:]
    P = np.trace(ETA @ Gvv) / 4
    Xh = np.trace(ETA @ Gvj) / 4
    Y = np.trace(ETA @ Gjj) / 4
    res = (np.linalg.norm(Gvv - P * ETA) + np.linalg.norm(Gvj - Xh * ETA) + np.linalg.norm(Gjj - Y * ETA))
    return P, 2 * Xh, Y, res / max(np.linalg.norm(G), 1e-300)


check("H0 EC quadratic form non-degenerate on 24 contortion components", np.linalg.matrix_rank(MQ) == NB)
P0, X0, Y0, r0 = invariants(leff_matrix(0.0, 1.0))
lam = np.sqrt((3 / 16) / -Y0) if Y0 < 0 else np.nan
info(f"gb = 0, lam = 1: Y = {Y0:.6g}  ->  lam = {lam:.6g} reproduces -(3/16) J5^2 (normalization of eps K J sums over all a,b,c)")
check("H1 minimal EC: Y(0) < 0, so lam is real (normalization fixed by the repo result)", np.isfinite(lam))

rows = []
for gb in [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, -1.3]:
    P, X, Y, r = invariants(leff_matrix(gb, lam))
    rows.append((gb, P, X, Y, r))
    info(f"gb={gb:+.2f}: P={P:.8f}  X={X:.8f}  Y={Y:.8f}  resid={r:.1e}")
check("H2 Y = -(3/16) gamma^2/(1+gamma^2) = -(3/16)/(1+gb^2) (Perez-Rovelli Eq. 15)",
      all(abs(Y + 3 / 16 / (1 + gb ** 2)) < 1e-12 for gb, P, X, Y, r in rows))
check("H3 only v.v, v.J, J.J invariants (no other structure)", all(r < 1e-12 for *_, r in rows))
Pfit = [P * (1 + gb ** 2) for gb, P, X, Y, r in rows]
check("H4 P(gb) = P0/(1+gb^2) (Taveras-Yunes structure)", np.ptp(Pfit) < 1e-12, f"P0 = {Pfit[0]:.6f}")
Xfit = [X * (1 + gb ** 2) for gb, P, X, Y, r in rows]
info(f"X(gb)(1+gb^2) values: {[round(x, 10) for x in Xfit]}")
check("H4b X(gb) = X0/(1+gb^2): cross term X(gb) d(gb).J = d G(gb).J, a total-derivative-type coupling",
      np.ptp(Xfit) < 1e-12, f"X0 = {Xfit[0]:.6f}")

# Canonical field: P0/(1+gb^2) (d gb)^2 = -(1/2)(d phi)^2  ->  gb = sinh(phi/F), F = sqrt(-2 P0)/kappa^{1/2}
P0v = Pfit[0]
F = np.sqrt(-2 * P0v)
info(f"P0 = {P0v:.6f} -> healthy (P0 < 0): {P0v < 0}; gb = sinh(phi/F), F = {F:.6f}/sqrt(kappa) = {F:.4f} M_P")
phis = np.linspace(-2, 2, 9) * F
Ycan = [-(3 / 16) / np.cosh(x / F) ** 2 for x in phis]
check("H5 canonical field: Y(phi) = -(3 kappa/16)/cosh^2(phi/F) is phi-dependent (shift symmetry broken)",
      P0v < 0 and np.ptp(Ycan) > 0.1)
P_, X_, Y_, _ = invariants(leff_matrix(0.7, lam))
check("H4c complete square: L_eff = -(3/4kappa)(v + J/2)^2/(1+gb^2)",
      abs(X_ - P_) < 1e-12 and abs(Y_ - P_ / 4) < 1e-12)
check("H6 Y and P even in gb: no gb-odd non-derivative coupling at tree level (no theta-bar shift)",
      all(abs(invariants(leff_matrix(g, lam))[2] - invariants(leff_matrix(-g, lam))[2]) < 1e-13 for g in (0.3, 1.7)))
