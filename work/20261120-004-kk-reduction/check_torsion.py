"""Check 2: torsion/contortion algebra in 5D and its 4D decomposition.

(a) contortion-from-torsion formula consistent with CONVENTIONS.md index order
    (Gam^lam_{mu nu}, derivative index mu). Tests the CONVENTIONS.md formula and ours.
(b) counting: T^A_{BC} in 5D = 50 = 5 (trace) + 10 (3-form) + 35 (rest); 4D zero-mode split.
(c) R(Gam_LC + K) = R_LC + nabla_mu V^mu + Q(K) at random points (generic metric, generic K
    and dK). F lives inside Gam_LC, so this shows K couples to F only through a total derivative.
(d) Q restricted to totally antisymmetric K equals -K_{ABC}K^{ABC}; Q is non-degenerate on the
    50-dim space of metric-compatible K, so vacuum EC (no spin source) forces K = 0.
(e) Dirac fermion: {G^C, G^{AB}} = 2 G^{CAB} (only K_[ABC] couples), and in 5D G^{ABC} is
    proportional to eps^{ABCDE} G_{DE} (the source is a 2-form bilinear).
(f) 4D limit: eps_{abcd} eps^{abce} = -6 delta_d^e, which gives the 4D four-fermion coefficient 3k/16.
(g) Levi-Civita spin connection of the KK vielbein: its totally antisymmetric part (the only part
    a Dirac fermion sees) is (k/6) F_ab in the [a b 5hat] slot: a Pauli term, linear in F.
"""
import itertools
import numpy as np

from geom import christoffel, ricci_scalar

rng = np.random.default_rng(7)
ETA5 = np.diag([-1.0, 1, 1, 1, 1])
D = 5


def perm_sign(p):
    p = list(p); s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s


def levi(n):
    e = np.zeros((n,) * n)
    for p in itertools.permutations(range(n)):
        e[p] = perm_sign(p)
    return e


def rand_T(n, g):
    """T_{lam mu nu} (all lower) antisymmetric in mu nu."""
    t = rng.normal(size=(n, n, n))
    return t - t.transpose(0, 2, 1)


def K_ours(Tl):
    # K_{lmn} = 1/2 (T_{lmn} + T_{nlm} - T_{mnl})
    return 0.5 * (Tl + Tl.transpose(1, 2, 0) - Tl.transpose(2, 0, 1))  # numpy: x.transpose(1,2,0)[l,m,n] = x[n,l,m]


def K_conv(Tl):
    # CONVENTIONS.md: K_{lmn} = 1/2 (T_{lmn} - T_{mln} - T_{nlm})
    return 0.5 * (Tl - Tl.transpose(1, 0, 2) - Tl.transpose(1, 2, 0))


def check_K_formula():
    g = ETA5 + 0.1 * np.diag(rng.normal(size=5))
    Tl = rand_T(D, g)
    out = []
    for name, f in (("ours  K=1/2(T_lmn+T_nlm-T_mnl)", K_ours), ("CONVENTIONS K=1/2(T_lmn-T_mln-T_nlm)", K_conv)):
        Kl = f(Tl)
        compat = np.max(np.abs(Kl + Kl.transpose(2, 1, 0)))       # K_{lmn} = -K_{nml}
        tors = np.max(np.abs(Kl - Kl.transpose(0, 2, 1) - Tl))      # K_lmn - K_lnm = T_lmn
        out.append((name, compat, tors))
    return out


def check_counting():
    # basis of T^A_{BC} antisym in BC
    basis = []
    for a in range(D):
        for b in range(D):
            for c in range(b + 1, D):
                t = np.zeros((D, D, D)); t[a, b, c] = 1; t[a, c, b] = -1
                basis.append(t)
    n = len(basis)
    eta = ETA5
    def trace(t): return np.einsum("aba->b", t)
    def axial(t):
        tl = np.einsum("ad,dbc->abc", eta, t)
        return (tl + tl.transpose(1, 2, 0) + tl.transpose(2, 0, 1)).reshape(-1)
    Mtr = np.array([trace(t) for t in basis]).T
    Max = np.array([axial(t) for t in basis]).T
    r_tr = np.linalg.matrix_rank(Mtr)
    r_ax = np.linalg.matrix_rank(Max)
    r_both = np.linalg.matrix_rank(np.vstack([Mtr, Max]))
    # 4D split of zero-mode components (tangent indices a=0..3, 5hat=4)
    blocks = {"T^a_bc": 4 * 6, "T^a_b5": 4 * 4, "T^5_ab": 6, "T^5_a5": 4}
    return n, r_tr, r_ax, r_both, blocks


def quad_Q(gi, K):
    """Q = g^{sn}(K^m_{ml}K^l_{ns} - K^m_{nl}K^l_{ms}); K[r,m,n] = K^r_{mn}."""
    return (np.einsum("sn,mml,lns->", gi, K, K) - np.einsum("sn,mnl,lms->", gi, K, K))


def check_RC_identity(n=5, trials=5):
    worst = 0.0
    for _ in range(trials):
        h = rng.normal(size=(n, n)) * 0.1
        g = np.diag([-1.0] + [1.0] * (n - 1)) + h + h.T
        dg = rng.normal(size=(n, n, n)) * 0.3; dg = dg + dg.transpose(0, 2, 1)
        ddg = rng.normal(size=(n, n, n, n)) * 0.3
        ddg = ddg + ddg.transpose(1, 0, 2, 3); ddg = ddg + ddg.transpose(0, 1, 3, 2)
        gi, G, dG = christoffel(g, dg, ddg)
        K = rng.normal(size=(n, n, n))           # generic (1,2) tensor
        dK = rng.normal(size=(n, n, n, n))       # dK[a] = d_a K
        R_full = ricci_scalar(gi, G + K, dG + dK)
        R_lc = ricci_scalar(gi, G, dG)
        dgi = -np.einsum("im,amn,nj->aij", gi, dg, gi)
        # V^m = g^{sn} K^m_{ns} - g^{ms} K^n_{ns}
        dV = (np.einsum("asn,mns->am", dgi, K) + np.einsum("sn,amns->am", gi, dK)
              - np.einsum("ams,nns->am", dgi, K) - np.einsum("ms,anns->am", gi, dK))
        V = np.einsum("sn,mns->m", gi, K) - np.einsum("ms,nns->m", gi, K)
        divV = np.trace(dV) + np.einsum("mml,l->", G, V)
        res = R_full - R_lc - divV - quad_Q(gi, K)
        worst = max(worst, abs(res) / (1 + abs(R_full)))
    return worst


def check_Q_antisym_and_rank():
    g = ETA5
    a = rng.normal(size=(D, D, D))
    Kl = sum(perm_sign(p) * a.transpose(p) for p in itertools.permutations(range(3)))
    K = np.einsum("rs,smn->rmn", g, Kl)
    q = quad_Q(g, K)
    q_expect = -np.einsum("abc,ad,be,cf,def->", Kl, g, g, g, Kl)
    # Q as a quadratic form on metric-compatible K parametrised by T (50 dims)
    basis = []
    for l in range(D):
        for m in range(D):
            for nn in range(m + 1, D):
                t = np.zeros((D, D, D)); t[l, m, nn] = 1; t[l, nn, m] = -1
                basis.append(np.einsum("rs,smn->rmn", g, K_ours(t)))
    N = len(basis)
    M = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            M[i, j] = 0.5 * (quad_Q(g, basis[i] + basis[j]) - quad_Q(g, basis[i]) - quad_Q(g, basis[j]))
    rank5 = np.linalg.matrix_rank(M)
    return abs(q - q_expect), rank5, N


def gammas5():
    I2 = np.eye(2); sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1])
    Z = np.zeros((2, 2))
    g0 = np.block([[I2, Z], [Z, -I2]]).astype(complex)                 # Dirac rep, mostly-minus
    gi = [np.block([[Z, s], [-s, Z]]).astype(complex) for s in (sx, sy, sz)]
    g5 = 1j * g0 @ gi[0] @ gi[1] @ gi[2]
    G = [1j * g0] + [1j * x for x in gi] + [g5]                          # mostly-plus 5D
    return G


def check_gamma():
    G = gammas5()
    eta = ETA5
    ok_cliff = all(np.allclose(G[a] @ G[b] + G[b] @ G[a], 2 * eta[a, b] * np.eye(4))
                   for a in range(5) for b in range(5))
    import math
    def antisym(idx):
        out = np.zeros((4, 4), complex)
        for p in itertools.permutations(range(len(idx))):
            m = np.eye(4, dtype=complex)
            for k in p:
                m = m @ G[idx[k]]
            out += perm_sign(p) * m
        return out / math.factorial(len(idx))
    err_ac = 0.0
    for c in range(5):
        for a in range(5):
            for b in range(5):
                lhs = G[c] @ antisym((a, b)) + antisym((a, b)) @ G[c]
                err_ac = max(err_ac, np.max(np.abs(lhs - 2 * antisym((c, a, b)))))
    # duality G^{ABC} = c eps^{ABCDE} G_{DE}/2 with eps^{01234} = -1 (eps_{01234}=+1, det eta=-1)
    eps = levi(5) * -1.0
    Glow = [eta[a, a] * G[a] for a in range(5)]
    ratios = []
    for (a, b, c) in itertools.combinations(range(5), 3):
        rhs = np.zeros((4, 4), complex)
        for d in range(5):
            for e in range(5):
                if eps[a, b, c, d, e] != 0:
                    rhs += 0.5 * eps[a, b, c, d, e] * (Glow[d] @ Glow[e] - Glow[e] @ Glow[d]) / 2
        lhs = antisym((a, b, c))
        k = np.trace(np.linalg.pinv(rhs.reshape(-1, 1)) @ lhs.reshape(-1, 1))
        ratios.append((k, np.max(np.abs(lhs - k * rhs))))
    ks = np.array([r[0] for r in ratios])
    dual_ok = np.allclose(ks, ks[0]) and max(r[1] for r in ratios) < 1e-12
    return ok_cliff, err_ac, dual_ok, ks[0]


def check_eps4():
    eps = levi(4)  # eps_{0123}=+1 (flat); eps^{0123} = -1
    eps_up = -eps
    return np.allclose(np.einsum("abcd,abce->de", eps, eps_up), -6 * np.eye(4))


def check_pauli():
    """omega_{[CAB]} for e^a = dx^a, e^5 = dy + k A, flat g, phi = 0, at a point with
    A = c (random), dA = s + f/2 (s symmetric random, f antisymmetric)."""
    k = 0.7
    c = rng.normal(size=4)
    s = rng.normal(size=(4, 4)); s = s + s.T
    f = rng.normal(size=(4, 4)); f = f - f.T
    dA = s + 0.5 * f                       # dA[m, n] = d_m A_n
    # vielbein E[A, M]; derivative dE[M, A, N]
    E = np.eye(5); E[4, :4] = k * c
    dE = np.zeros((5, 5, 5)); dE[:4, 4, :4] = k * dA
    g = E.T @ ETA5 @ E
    dg = np.einsum("mAn,AB,Bp->mnp", dE, ETA5, E) + np.einsum("An,AB,mBp->mnp", E, ETA5, dE)
    ddg = np.zeros((5, 5, 5, 5))
    gi, Gam, _ = christoffel(g, dg, ddg)
    Einv = np.linalg.inv(E)                # Einv[N, A] = e_A^N
    El = np.einsum("AB,BN->AN", ETA5, E)   # e_{A N}
    dEl = np.einsum("AB,mBN->mAN", ETA5, dE)
    # omega_{M A B} = - e_B^N (d_M e_{A N} - Gam^P_{M N} e_{A P})
    om = -np.einsum("NB,mAN->mAB", Einv, dEl) + np.einsum("NB,pmN,Ap->mAB", Einv, Gam, El)
    omT = np.einsum("mC,mAB->CAB", Einv, om)    # tangent first index
    tot = sum(perm_sign(p) * omT.transpose(p) for p in itertools.permutations(range(3))) / 6
    P = np.zeros((5, 5, 5)); P[:4, :4, 4] = (k / 6) * f
    pred = sum(perm_sign(p) * P.transpose(p) for p in itertools.permutations(range(3))) / 2
    # f with tangent indices: e^a = dx^a so f_ab = f_mn here
    r = np.max(np.abs(tot - pred)); r_minus = np.max(np.abs(tot + pred))
    only_one5 = np.max(np.abs(tot[:4, :4, :4]))
    return min(r, r_minus), ("+" if r < r_minus else "-"), only_one5, np.max(np.abs(om + om.transpose(0, 2, 1)))


def run():
    res = []
    kf = check_K_formula()
    ours, conv = kf
    res.append(("contortion formula K_lmn = 1/2(T_lmn + T_nlm - T_mnl) is metric-compatible and reproduces T",
                ours[1] < 1e-12 and ours[2] < 1e-12, f"compat {ours[1]:.1e}, torsion {ours[2]:.1e}"))
    res.append(("INFO CONVENTIONS.md formula 1/2(T_lmn - T_mln - T_nlm) with Gam^lam_{mu nu}, mu=derivative index",
                True, f"metric-compat violation {conv[1]:.2f}, torsion mismatch {conv[2]:.2f} (nonzero => formula is for the other index order; reported in README)"))
    n, r_tr, r_ax, r_both, blocks = check_counting()
    res.append(("torsion counting 50 = 5 trace + 10 three-form + 35 tensor; 4D split 24+16+6+4",
                n == 50 and r_tr == 5 and r_ax == 10 and r_both == 15 and sum(blocks.values()) == 50,
                f"dim {n}, trace rank {r_tr}, 3-form rank {r_ax}, joint {r_both}, remainder {n - r_both}, blocks {blocks}"))
    w5 = check_RC_identity(5); w4 = check_RC_identity(4)
    res.append(("R(LC+K) = R_LC + div V + Q(K) (generic metric, K, dK; D=5 and D=4)",
                w5 < 1e-11 and w4 < 1e-11, f"max rel residual {max(w5, w4):.1e}"))
    dq, rank5, N = check_Q_antisym_and_rank()
    res.append(("Q(K) = -K_ABC K^ABC for totally antisymmetric K; Q non-degenerate on 50 compatible K",
                dq < 1e-10 and rank5 == N == 50, f"|dQ| {dq:.1e}, rank {rank5}/{N}"))
    ok_cliff, err_ac, dual_ok, kdual = check_gamma()
    res.append(("5D Clifford algebra; {G^C,G^AB} = 2 G^CAB; G^ABC = c eps^ABCDE G_DE/2 with one constant c",
                ok_cliff and err_ac < 1e-12 and dual_ok, f"anticomm err {err_ac:.1e}, c = {complex(kdual):.3f}"))
    res.append(("eps_abcd eps^abce = -6 delta (4D Lorentzian), gives J_abc J^abc = -6 A^2",
                check_eps4(), ""))
    r, sgn, only5, anti = check_pauli()
    res.append(("KK Levi-Civita spin connection: omega_[CAB] = (sign) (k/6) F_ab in [a b 5hat], zero for [abc]",
                r < 1e-12 and only5 < 1e-12 and anti < 1e-12, f"sign {sgn}, residual {r:.1e}, [abc] part {only5:.1e}"))
    return res


if __name__ == "__main__":
    for name, ok, det in run():
        print(("PASS" if ok else "FAIL"), name, "|", det)
