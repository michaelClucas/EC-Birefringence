"""Task 4: CFJ / phi F Ftilde equivalence, sign, dispersion relation, rotation, gauge variation.

Flat comoving coordinates (tau, x, y, z) right-handed; by check_maxwell_zero (c) the FRW
equations are these with p_m the conformal-coordinate covariant components.
 1. -(g/4) phi F Ftilde and -(1/2) p_m A_n Ftilde^{mn} give identical Euler-Lagrange equations
    iff p_m = -g d_m phi (CONVENTIONS.md states p_m = +g d_m phi: sign error, see README).
 2. EOM for constant p: d_m F^{mn} = p_m Ftilde^{mn}; for p = (p0,0,0,0), A_0 = 0, div A = 0:
    A'' - lap A = p0 curl A.
 3. Helicity modes along +z: omega_pm^2 = k^2 -+ p0 k, (e_x +- i e_y).
 4. Linear polarization: e_+ e^{i th} + e_- e^{-i th} = 2(cos th e_x - sin th e_y), th = p0 dtau/2
    -> polarization angle about +k (right-handed) changes by chi = -(1/2) int p_tau dtau.
 5. Independent numerical check: evolve the EL equations (derived in step 2, not the helicity
    algebra) for a linearly polarized packet and measure chi.
 6. Gauge variation of -(1/2) xi S_m A_n Ftilde^{mn}: equals a total derivative plus
    (1/2) xi chi d_n S_m Ftilde^{mn}, which vanishes iff dS = 0 (true for FRW S_m = (S_tau(tau),0,0,0)).
"""
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from sympy.calculus.euler import euler_equations

from common import ETA, R4, perm, report

X = sp.symbols("x0:4")


def fields():
    return [sp.Function(f"A{i}")(*X) for i in R4]


def Fdown(A):
    return [[sp.diff(A[n], X[m]) - sp.diff(A[m], X[n]) for n in R4] for m in R4]


def Ftilde_up(F):
    return [[sp.Rational(1, 2) * sum(-perm(m, n, r, s) * F[r][s] for r in R4 for s in R4)
             for n in R4] for m in R4]


def maxwell(F):
    return -sp.Rational(1, 4) * sum(ETA[m, m] * ETA[n, n] * F[m][n]**2 for m in R4 for n in R4)


def FFt(F):
    Ft = Ftilde_up(F)
    return sum(F[m][n] * Ft[m][n] for m in R4 for n in R4)


def cfj(A, p):
    Ft = Ftilde_up(Fdown(A))
    return -sp.Rational(1, 2) * sum(p[m] * A[n] * Ft[m][n] for m in R4 for n in R4)


def el(L, A):
    return [sp.expand(e.lhs) for e in euler_equations(L, A, X)]


def main():
    ok = True
    A = fields()
    F = Fdown(A)
    gg = sp.Symbol("g")
    phi = sp.Function("phi")(*X)
    E1 = el(maxwell(F) - gg / 4 * phi * FFt(F), A)
    dphi = [sp.diff(phi, x) for x in X]
    E2 = el(maxwell(F) + cfj(A, [-gg * d for d in dphi]), A)
    E3 = el(maxwell(F) + cfj(A, [gg * d for d in dphi]), A)
    ok &= report("1. -(g/4) phi FFt == CFJ with p_m = -g d_m phi (same EL equations)",
                 all(sp.expand(a - b) == 0 for a, b in zip(E1, E2)))
    ok &= report("1. ... and NOT with p_m = +g d_m phi (CONVENTIONS.md sign)",
                 not all(sp.expand(a - b) == 0 for a, b in zip(E1, E3)))
    # FFt = 2 d_m (A_n Ftilde^{mn}) identity
    Ft = Ftilde_up(F)
    ok &= report("1. F_mn Ftilde^mn = 2 d_m(A_n Ftilde^mn)",
                 sp.expand(FFt(F) - 2 * sum(sp.diff(A[n] * Ft[m][n], X[m]) for m in R4 for n in R4)) == 0)

    # 2. EOM with constant p
    p = sp.symbols("p0:4")
    E = el(maxwell(F) + cfj(A, list(p)), A)
    Fup = [[ETA[m, m] * ETA[n, n] * F[m][n] for n in R4] for m in R4]
    eom = [sum(sp.diff(Fup[m][n], X[m]) for m in R4) - sum(p[m] * Ft[m][n] for m in R4) for n in R4]
    # euler_equations returns dL/dA - d(dL/ddA) = 0 ; compare up to overall factor
    lam = [sp.simplify(E[n] / eom[n]) for n in R4]
    ok &= report("2. EL equations are d_m F^mn = p_m Ftilde^mn (same overall factor for all n)",
                 all(l_ == lam[0] for l_ in lam) and lam[0].is_number)

    # restrict to A = (0, Ax(t,z), Ay(t,z), 0), p = (p0,0,0,0)
    t, z, p0 = X[0], X[3], sp.Symbol("p_0")
    Ax, Ay = sp.Function("Ax")(t, z), sp.Function("Ay")(t, z)
    sub = {A[0]: 0, A[3]: 0, A[1]: Ax, A[2]: Ay}
    eqx = sp.expand(eom[1].subs({p[0]: p0, p[1]: 0, p[2]: 0, p[3]: 0}).subs(sub).doit())
    eqy = sp.expand(eom[2].subs({p[0]: p0, p[1]: 0, p[2]: 0, p[3]: 0}).subs(sub).doit())
    curl_x, curl_y = -sp.diff(Ay, z), sp.diff(Ax, z)
    wave = lambda f: -sp.diff(f, t, 2) + sp.diff(f, z, 2)
    ok &= report("2. 1D reduction: A'' - A_zz = p0 curl A",
                 sp.expand(eqx - (wave(Ax) + p0 * curl_x)) == 0 and sp.expand(eqy - (wave(Ay) + p0 * curl_y)) == 0)

    # 3. dispersion for e_pm = (1, +-i) e^{i(kz - w t)}
    k, w = sp.symbols("k omega", positive=True)
    res = {}
    for sgn in (1, -1):
        ph = sp.exp(sp.I * (k * z - w * t))
        s = {Ax: ph, Ay: sgn * sp.I * ph}
        ex = sp.simplify(eqx.subs(s).doit() / ph)
        ey = sp.simplify(eqy.subs(s).doit() / ph)
        # ex = 0 is the dispersion relation; ey must be i*sgn*ex (same mode, consistent)
        res[sgn] = (sp.expand(ex), sp.simplify(ey - sgn * sp.I * ex))
    print("INFO: dispersion e_+:", res[1][0], "= 0 ; e_-:", res[-1][0], "= 0")
    ok &= report("3. omega_+^2 = k^2 - p0 k, omega_-^2 = k^2 + p0 k",
                 sp.expand(res[1][0] - (w**2 - k**2 + p0 * k)) == 0
                 and sp.expand(res[-1][0] - (w**2 - k**2 - p0 * k)) == 0
                 and res[1][1] == 0 and res[-1][1] == 0)
    print("INFO: omega_pm ~ k -+ p0/2 for k >> |p0|")

    # 4. superposition algebra
    th = sp.Symbol("theta", real=True)
    v = sp.Matrix([1, sp.I]) * sp.exp(sp.I * th) + sp.Matrix([1, -sp.I]) * sp.exp(-sp.I * th)
    ok &= report("4. e_+ e^{i th} + e_- e^{-i th} = 2(cos th, -sin th)",
                 sp.simplify((v - 2 * sp.Matrix([sp.cos(th), -sp.sin(th)])).applyfunc(lambda e: e.rewrite(sp.cos))) == sp.zeros(2, 1))

    # 5. numerical evolution of the EL equations (spectral in z)
    Np, Lbox = 1024, 200.0
    zz = np.linspace(0, Lbox, Np, endpoint=False)
    kk = 2 * np.pi * np.fft.fftfreq(Np, d=Lbox / Np)
    p0v, Tend, k0 = 0.05, 10.0, 2 * np.pi * 8 / 1.0 / 8  # carrier k0 = 2 pi
    env = np.exp(-((zz - 60.0) / 10.0) ** 2)
    f0 = env * np.cos(k0 * zz)
    df0 = np.real(np.fft.ifft(1j * kk * np.fft.fft(f0)))

    def dz(fv):
        return np.real(np.fft.ifft(1j * kk * np.fft.fft(fv)))

    def rhs(_, y):
        ax, ay, vx, vy = np.split(y, 4)
        # from eqx = 0: Ax_tt = Ax_zz + p0 curl_x = Ax_zz - p0 Ay_z ; Ay_tt = Ay_zz + p0 Ax_z
        axx, ayy = dz(dz(ax)), dz(dz(ay))
        return np.concatenate([vx, vy, axx - p0v * dz(ay), ayy + p0v * dz(ax)])

    y0 = np.concatenate([f0, 0 * f0, -df0, 0 * f0])  # right-moving A_x = f(z - t)
    sol = solve_ivp(rhs, (0, Tend), y0, rtol=1e-9, atol=1e-11, method="DOP853")
    ax, ay = np.split(sol.y[:, -1], 4)[:2]
    chi_num = np.arctan2(np.sum(ax * ay), np.sum(ax * ax))
    chi_pred = -p0v * Tend / 2
    print(f"INFO: numerical chi = {chi_num:.6f} rad, predicted -p0 dtau/2 = {chi_pred:.6f} rad")
    ok &= report("5. numerical rotation matches chi = -(1/2) p0 dtau (1e-3 rel)",
                 abs(chi_num - chi_pred) < 1e-3 * abs(chi_pred) + 1e-5)
    # the CFJ sign flip check for the derived EOM signs:
    ok &= report("5. RHS used in simulation equals the symbolic EOM",
                 sp.expand(eqx.subs({sp.diff(Ax, t, 2): 0}) - (sp.diff(Ax, z, 2) - p0 * sp.diff(Ay, z))) == 0
                 and sp.expand(eqy.subs({sp.diff(Ay, t, 2): 0}) - (sp.diff(Ay, z, 2) + p0 * sp.diff(Ax, z))) == 0)

    # 6. gauge variation
    xi = sp.Symbol("xi")
    chi = sp.Function("chi")(*X)
    S = [sp.Function(f"S{i}")(*X) for i in R4]
    Lc = cfj(A, [xi * s for s in S])
    Ag = [A[i] + sp.diff(chi, X[i]) for i in R4]
    dL = sp.expand(cfj(Ag, [xi * s for s in S]) - Lc)
    div = sum(sp.diff(sp.Rational(1, 2) * xi * chi * sum(S[m] * Ft[m][n] for m in R4), X[n]) for n in R4)
    rem = sp.expand(dL + div - sp.Rational(1, 2) * xi * chi *
                    sum(sp.diff(S[m], X[n]) * Ft[m][n] for m in R4 for n in R4))
    ok &= report("6. delta L_CFJ = -d_n(xi chi S_m Ft^mn / 2) + (xi/2) chi d_n S_m Ft^mn", rem == 0)
    ph_ = sp.Function("vphi")(*X)
    resid = sp.expand(sum(sp.diff(sp.diff(ph_, X[m]), X[n]) * Ft[m][n] for m in R4 for n in R4))
    ok &= report("6. remainder vanishes for S_m = d_m vphi (curl-free)", resid == 0)
    Sfrw = [sp.Function("Stau")(X[0]), 0, 0, 0]
    resid2 = sp.expand(sum(sp.diff(Sfrw[m], X[n]) * Ft[m][n] for m in R4 for n in R4))
    ok &= report("6. remainder vanishes for FRW S_m = (S_tau(tau),0,0,0)", resid2 == 0)

    # 7. map to the sky and to the paper's beta (literature/conventions_beta.md):
    # photon moves along +z toward the observer, who faces f = -z. With up u = +y the observer's
    # right is r = f x u. If r = +x, then x -> y (chi > 0) runs right -> up: counterclockwise as seen.
    f_, u_ = sp.Matrix([0, 0, -1]), sp.Matrix([0, 1, 0])
    r_ = f_.cross(u_)
    ccw = r_ == sp.Matrix([1, 0, 0])
    ok &= report("7. chi > 0 (right-handed about k) is counterclockwise on the sky as seen by the observer", ccw)
    # paper/Komatsu convention: beta > 0 clockwise as seen  ->  beta = -chi
    gsym, dphi_ = sp.symbols("g Delta_phi")
    chi_phi = -sp.Rational(1, 2) * (-gsym * dphi_)  # chi = -(1/2) int p_tau dtau, p = -g dphi
    beta = -chi_phi
    print("INFO: our conventions, L = -(g/4) phi F Ftilde: chi =", chi_phi, " beta_paper = -chi =", beta)
    ok &= report("7. beta_paper = -(g/2)(phi_obs - phi_emit) for L = -(g/4) phi F Ftilde in CONVENTIONS.md",
                 sp.simplify(beta + gsym * dphi_ / 2) == 0)
    print("INFO: CFJ form: beta_paper = +(1/2) int p_tau dtau (p_tau: conformal-coordinate covariant component)")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
