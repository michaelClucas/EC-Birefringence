"""Skeptic, independent check of the CFJ dispersion relation and the sign of beta.

Written without importing any code from work/20261120-003 or -004.

Route:
 1. Build L = -1/4 F_mn F^mn - (g/4) phi(t) F_mn Ft^mn in flat space, signature (-+++),
    Ft^mn = 1/2 eps^mnrs F_rs with eps^{0123} = -1 (our eps_{0123} = +sqrt(-g)).
    Euler-Lagrange equations via sympy.
 2. Plane wave A_j = a_j(t) exp(i k.x) along an ARBITRARY direction k (not the z axis):
    extract the 3x3 matrix ODE a'' = -M a, check the helicity eigenvalues
    omega^2 = k^2 -/+ p_t k (p_t = -g phidot) for helicity +/- about k.
 3. Integrate the ODE numerically (no helicity algebra), read the polarization at the
    observer in the HEALPix/COSMO basis (e_theta, e_phi) of the line of sight n = -k,
    with Q = |a_th|^2 - |a_ph|^2, U = 2 Re(a_th conj(a_ph)) (COSMO: U = -U_IAU).
    beta_cosmo = change of psi = (1/2) atan2(U, Q). 2202.13919: Q +- iU = P e^{+-2i beta},
    and 2509.13654 Eq. 4 rotates (E + iB) by e^{2i beta}, i.e. beta = delta psi_COSMO.
 4. Check separately that psi_COSMO increasing is clockwise as seen by the observer
    (up = north = -e_theta, observer's right = n x up).
 5. Control: flipping the sign of eps (Komatsu 2202.13919: eps^{0123} = +1/sqrt(-g))
    must flip the sign of beta; this is how the check could fail.
"""
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from sympy.calculus.euler import euler_equations

t, x, y, z = X = sp.symbols("t x y z", real=True)
ETA = np.diag([-1, 1, 1, 1])
results = []


def rep(name, ok, info=""):
    results.append(ok)
    print(("PASS" if ok else "FAIL") + ": " + name + (("  | " + info) if info else ""))


def eom_matrix(eps_sign):
    A = [sp.Function(f"A{i}")(*X) for i in range(4)]
    g, w = sp.symbols("g phidot", real=True)
    phi = w * t
    F = [[sp.diff(A[n], X[m]) - sp.diff(A[m], X[n]) for n in range(4)] for m in range(4)]
    Fup = [[ETA[m, m] * ETA[n, n] * F[m][n] for n in range(4)] for m in range(4)]
    epsU = lambda *i: eps_sign * (-1) * sp.LeviCivita(*i)  # eps_sign=+1: ours, eps^{0123} = -1
    Ft = [[sp.Rational(1, 2) * sum(epsU(m, n, r, s) * F[r][s] for r in range(4) for s in range(4))
           for n in range(4)] for m in range(4)]
    L = (-sp.Rational(1, 4) * sum(F[m][n] * Fup[m][n] for m in range(4) for n in range(4))
         - g / 4 * phi * sum(F[m][n] * Ft[m][n] for m in range(4) for n in range(4)))
    eqs = euler_equations(L, A, X)
    kx, ky, kz = sp.symbols("k_x k_y k_z", real=True)
    a = [sp.Function(f"a{i}")(t) for i in range(1, 4)]
    ph = sp.exp(sp.I * (kx * x + ky * y + kz * z))
    sub = {A[0]: 0, A[1]: a[0] * ph, A[2]: a[1] * ph, A[3]: a[2] * ph}
    rows = []
    for j in (1, 2, 3):
        e = sp.expand(sp.simplify(eqs[j].lhs.subs(sub).doit() / ph))
        # e = c2 * a_j'' + sum_l c_l a_l  (+ terms)
        c2 = e.coeff(sp.diff(a[j - 1], t, 2))
        row = [sp.simplify(-e.coeff(a[l]) / c2) for l in range(3)]  # a_j'' = -sum row_l a_l  ->  M = -row
        rest = sp.simplify(e - c2 * sp.diff(a[j - 1], t, 2) - sum(e.coeff(a[l]) * a[l] for l in range(3)))
        assert rest == 0, rest
        rows.append([-r for r in row])
    return sp.Matrix(rows), (g, w, kx, ky, kz)


def main():
    M, (g, w, kx, ky, kz) = eom_matrix(+1)
    k = sp.sqrt(kx**2 + ky**2 + kz**2)
    p_t = -g * w  # CFJ vector p_mu = -g d_mu phi (claim in 003)
    # helicity vectors about k: e_pm = (u1 +- i u2)/sqrt2 with (u1, u2, khat) right-handed
    rng = np.random.default_rng(7)
    kv = rng.normal(size=3)
    kv /= np.linalg.norm(kv)
    kv *= 25.0  # k >> |p| so the helicity splitting is linear in p
    subsk = {kx: kv[0], ky: kv[1], kz: kv[2]}
    khat = kv / np.linalg.norm(kv)
    u1 = np.cross(khat, [0.3, -0.5, 0.8]); u1 /= np.linalg.norm(u1)
    u2 = np.cross(khat, u1)
    Mn = sp.lambdify((g, w), M.subs(subsk), "numpy")
    gv, wv = 0.8, 0.0625  # p_t = -0.05
    Mv = np.array(Mn(gv, wv), dtype=complex)
    kn = np.linalg.norm(kv)
    for sgn, lab in ((1, "+"), (-1, "-")):
        e = (u1 + sgn * 1j * u2) / np.sqrt(2)
        lam = (e.conj() @ Mv @ e)
        resid = np.linalg.norm(Mv @ e - lam * e)
        pred = kn**2 - sgn * (-gv * wv) * kn
        rep(f"helicity {lab} about k is an eigenmode with omega^2 = k^2 {'-' if sgn > 0 else '+'} p_t k",
            resid < 1e-12 and abs(lam - pred) < 1e-12, f"omega^2 = {lam.real:.6f}, predicted {pred:.6f}")
    # symbolic: transverse eigenvalues
    rep("longitudinal mode (A_0 = 0 gauge) is pure gauge: M khat = 0",
        np.linalg.norm(Mv @ khat) < 1e-10)

    # numerical evolution, positive-frequency wave along +k, initial linear polarization u1
    T_end = 10.0

    def run(Mmat, pol0):
        om = kn
        y0 = np.concatenate([pol0, -1j * om * pol0]).astype(complex)

        def rhs(_, Y):
            aa, vv = Y[:3], Y[3:]
            return np.concatenate([vv, -Mmat @ aa])
        s = solve_ivp(rhs, (0, T_end), y0, rtol=1e-11, atol=1e-13, method="DOP853")
        return s.y[:3, -1]

    # line of sight n = -khat; HEALPix basis at n
    n = -khat
    th, ph_ = np.arccos(n[2]), np.arctan2(n[1], n[0])
    e_th = np.array([np.cos(th) * np.cos(ph_), np.cos(th) * np.sin(ph_), -np.sin(th)])
    e_ph = np.array([-np.sin(ph_), np.cos(ph_), 0.0])

    def psi_cosmo(avec):
        # strip global phase: polarization of a linearly polarized complex vector
        at, ap = e_th @ avec, e_ph @ avec
        Q = abs(at)**2 - abs(ap)**2
        U = 2 * np.real(at * np.conj(ap))
        return 0.5 * np.arctan2(U, Q)

    beta_num = psi_cosmo(run(Mv, u1)) - psi_cosmo(u1.astype(complex))
    dphi = wv * T_end
    beta_claim = -(gv / 2) * dphi  # 003 claim: beta = -(g/2)(phi_0 - phi_LSS)
    rep("beta (COSMO psi change, numerical) = -(g/2) Delta phi  [003 claim]",
        abs(beta_num - beta_claim) < 2e-3 * abs(beta_claim), f"beta_num = {beta_num:.6f}, claim = {beta_claim:.6f} rad")
    rep("equivalently beta = +1/2 int p_t dt with p = -g dphi",
        abs(beta_num - 0.5 * (-gv * wv) * T_end) < 2e-3 * abs(beta_claim))

    # psi_COSMO increasing = clockwise as seen by the observer
    up = -e_th  # north
    right = np.cross(n, up)  # observer faces n
    # angle measured from up toward right is clockwise on a clock face
    ang = lambda v: np.arctan2(v @ right, v @ up)
    psi = 0.2
    v0 = np.cos(0) * e_th + np.sin(0) * e_ph
    v1 = np.cos(psi) * e_th + np.sin(psi) * e_ph
    d = np.angle(np.exp(1j * (ang(v1) - ang(v0))))
    rep("psi_COSMO increasing is clockwise as seen by the observer (so beta > 0 clockwise)", d > 0,
        f"clock-face angle change {d:.3f} for dpsi = {psi}")
    rep("east (e_phi) is on the observer's left with north up", np.dot(e_ph, right) < -0.99)

    # control: Komatsu eps sign flips beta
    Mk, _ = eom_matrix(-1)
    Mkv = np.array(sp.lambdify((g, w), Mk.subs(subsk), "numpy")(gv, wv), dtype=complex)
    beta_k = psi_cosmo(run(Mkv, u1)) - psi_cosmo(u1.astype(complex))
    rep("control: with eps^{0123} = +1/sqrt(-g) (2202.13919) beta = +(g/2) Delta phi, as Komatsu states",
        abs(beta_k - (gv / 2) * dphi) < 2e-3 * abs(beta_claim), f"beta = {beta_k:.6f}")
    # control: a wrong sky mapping would fail (beta in IAU = -beta_COSMO)
    rep("control: the check distinguishes signs (|beta| >> tolerance)", abs(beta_claim) > 1e-2)
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
