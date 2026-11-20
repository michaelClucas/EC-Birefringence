"""Check 3: reduction of the 5D Chern-Simons term and related 4D photon algebra.

(a) 5D CS density  L5 = epsT^{MNPQR} Acal_M d_N Acal_P d_Q Acal_R  (epsT = permutation symbol,
    epsT^{01234}=+1, coordinate 4 = y), i.e. Acal^dAcal^dAcal = L5 d^4x^dy.
    Zero-mode ansatz Acal = a_mu dx^mu + theta (dy + c A_mu dx^mu), all fields y-independent
    (A = KK photon, c = its KK normalisation; a = 4D part of the bulk vector).
    Claim (mod total derivatives):
      L4 = epsT^{mnrs} [ 3 theta d_m a_n d_r a_s + 3 c theta^2 d_m a_n d_r A_s + c^2 theta^3 d_m A_n d_r A_s ]
    Verified by equality of Euler-Lagrange expressions of L5 - L4 for all 9 fields
    (total derivatives have identically zero EL expressions).
(b) f^f = -(1/2) sqrt(-g) F_{mn} Ft^{mn} d^4x with eps_{0123} = +sqrt(-g), Ft^{mn} = eps^{mnrs}F_rs/2.
(c) A purely mixed coupling -(g/2) b G_mn Ft^mn between the photon F and a second massless
    vector G gives no polarization rotation for a photon that starts with no G admixture.
"""
import itertools
import numpy as np
import sympy as sp
from scipy.linalg import expm
from sympy.calculus.euler import euler_equations


def perm_sign(p):
    p = list(p); s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s


def check_cs_reduction():
    X = sp.symbols("x0:4")
    c = sp.Symbol("c")
    a = [sp.Function(f"a{i}")(*X) for i in range(4)]
    A = [sp.Function(f"A{i}")(*X) for i in range(4)]
    th = sp.Function("theta")(*X)
    cal = [a[m] + c * th * A[m] for m in range(4)] + [th]
    d = lambda f, N: sp.diff(f, X[N]) if N < 4 else 0
    L5 = 0
    for p in itertools.permutations(range(5)):
        M, N, P, Q, R = p
        if N == 4 or Q == 4:
            continue  # d_y = 0
        L5 += perm_sign(p) * cal[M] * d(cal[P], N) * d(cal[R], Q)
    L4 = 0
    for p in itertools.permutations(range(4)):
        m, n, r, s = p
        sg = perm_sign(p)
        L4 += sg * (3 * th * d(a[n], m) * d(a[s], r)
                    + 3 * c * th**2 * d(a[n], m) * d(A[s], r)
                    + c**2 * th**3 * d(A[n], m) * d(A[s], r))
    diff = sp.expand(L5 - L4)
    eqs = euler_equations(diff, a + A + [th], X)
    bad = [e for e in eqs if sp.expand(e.lhs) != 0]
    # control: dropping the factor 3 on theta f^f must fail
    L4w = L4 - sum(perm_sign(p) * th * d(a[p[1]], p[0]) * d(a[p[3]], p[2]) for p in itertools.permutations(range(4)))
    eqw = euler_equations(sp.expand(L5 - L4w), [th], X)
    control_fails = any(sp.expand(e.lhs) != 0 for e in eqw)
    return len(bad) == 0, control_fails


def check_ffdual():
    rng = np.random.default_rng(3)
    F = rng.normal(size=(4, 4)); F = F - F.T
    epsT = np.zeros((4,) * 4)
    for p in itertools.permutations(range(4)):
        epsT[p] = perm_sign(p)
    eps_up = -epsT                                   # flat: eps^{0123} = -1
    Ft = 0.5 * np.einsum("mnrs,rs->mn", eps_up, F)
    FFt = np.einsum("mn,mn->", Ft, F)              # Ft^{mn} F_{mn}
    wedge = 0.25 * np.einsum("mnrs,mn,rs->", epsT, F, F)   # f^f coefficient of dx0123
    return np.isclose(wedge, -0.5 * FFt)


def check_mixed_no_rotation():
    # helicity +/- evolution of (A, G) amplitudes: exp(+- i Phi M), M = [[0,1],[1,0]]
    Phi = 0.37
    M = np.array([[0, 1.0], [1.0, 0]])
    psi0 = np.array([1.0, 0.0])
    ap = (expm(1j * Phi * M) @ psi0)[0]
    am = (expm(-1j * Phi * M) @ psi0)[0]
    rot = 0.5 * (np.angle(ap) - np.angle(am))
    # diagonal control: M = diag(1, 0) gives rotation Phi
    Md = np.diag([1.0, 0.0])
    rd = 0.5 * (np.angle((expm(1j * Phi * Md) @ psi0)[0]) - np.angle((expm(-1j * Phi * Md) @ psi0)[0]))
    return abs(rot) < 1e-14 and np.isclose(abs(ap), np.cos(Phi)) and np.isclose(rd, Phi), rot, abs(ap)


def run():
    ok_cs, ctrl = check_cs_reduction()
    ok_mix, rot, amp = check_mixed_no_rotation()
    return [
        ("5D CS A^F^F on S1 -> 3 theta f^f + 3c theta^2 f^F_KK + c^2 theta^3 F_KK^F_KK (EL-equivalent)",
         ok_cs and ctrl, f"EL match {ok_cs}; control with coefficient 2 instead of 3 fails: {ctrl}"),
        ("f^f = -(1/2) sqrt(-g) F Ft d^4x  (eps_0123 = +sqrt(-g))", bool(check_ffdual()), ""),
        ("mixed coupling b G Ft: no rotation for an initial pure photon (amplitude cos Phi only)",
         bool(ok_mix), f"rotation {rot:.1e}, photon amplitude {amp:.4f} = cos(0.37)"),
    ]


if __name__ == "__main__":
    for name, ok, det in run():
        print(("PASS" if ok else "FAIL"), name, "|", det)
