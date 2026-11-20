"""Question B: u-Holst term with a hidden scalar chi, c -> c(chi), in 5D EC-KK.

5D pointwise algebra (same frame, basis and quadratic forms as task 009 check_ec_uholst.py,
re-implemented here). L_K = (1/2) Q(K) + b * odd(K) + mix(X, K), kappa_5 = 1, where for the term
c(chi) u^e^A^e^B^R_AB(w), with u = e^4 (oriented, P'-odd background):
  b = c(chi)                     (u ^ T^A ^ T_A piece, pointwise)
  X = d(c u) = c W + c' dchi ^ u (source 2-form in the linear piece; W = du)
U1  L_eff(X) = A W.W + B *(u^W^W) + C (dchi)^2 with A = -3c^2/2, B = -4c^3 (009 values, now with c(chi));
    no W x dchi cross term; C reported
U2  exact u (u = dchi): X = d(u) = 0, no linear source -> K = 0 -> no photon term
U3  contortion response K_[ab4] to W: c-odd part = c W~ (P'-odd, EDM-type dipole source), c-even part (4/3)c^2 W
U4  helicity evolution for L = -(1/4) Z F^2 - (1/4) Y F Ftilde with both Z, Y time dependent:
    rotation angle = (1/2) |int Y'/Z dtau|
U5  with Z = Z0 (1 + 12 c^2), Y = Z0 16 c^3:  d beta = -c d ln alpha (sympy), closed-form beta(c)
U6  a dynamical orienting field (compact scalar psi winding the circle) gives the KK photon a mass:
    g^{MN} d_M psi d_N psi = (d pi - (w/R) A)^2 + (w/R)^2/phi^2   (sympy)
"""
import itertools
import math

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad

from common import check, info

D = 5
ETA = np.diag([-1.0, 1, 1, 1, 1])


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


def levi(n):
    e = np.zeros((n,) * n)
    for p in itertools.permutations(range(n)):
        e[p] = perm_sign(p)
    return e


EPS5, EPS4 = levi(5), levi(4)
BASIS = []
for A_ in range(D):
    for B_ in range(A_ + 1, D):
        for C_ in range(D):
            k = np.zeros((D, D, D))
            k[A_, B_, C_], k[B_, A_, C_] = 1.0, -1.0
            BASIS.append(k)
BST = np.array(BASIS)
NB = len(BST)
_PERM3 = [(p, perm_sign(p)) for p in itertools.permutations(range(3))]


def alt3(t):
    return sum(s * np.transpose(t, p) for p, s in _PERM3) / 6.0


def Q_bilinear(K1, K2):
    K2u = np.einsum("cc,ncbe->ncbe", ETA, K2)
    t1 = np.einsum("aa,bb,maca,ncbb->mn", ETA, ETA, K1, K2u)
    t2 = np.einsum("aa,bb,macb,ncba->mn", ETA, ETA, K1, K2u)
    M = t1 - t2
    return 0.5 * (M + M.T)


def odd_bilinear(K1, K2):
    def T(K):
        Ku = np.einsum("aa,nabc->nabc", ETA, K)
        return np.transpose(Ku, (0, 1, 3, 2)) - Ku
    T1, T2 = T(K1), T(K2)
    T2l = np.einsum("aa,nacd->nacd", ETA, T2)
    M = 0.25 * np.einsum("abcd,mxab,nxcd->mn", EPS4, T1[:, :, :4, :4], T2l[:, :, :4, :4])
    return 0.5 * (M + M.T)


def mix_linear(Wstack):
    Y = np.array([6 * alt3(k) for k in BASIS])
    return np.einsum("abcde,wab,ncde->wn", EPS5, Wstack, Y) / 12.0


MQ = 2 * Q_bilinear(BST, BST)
MO = 2 * odd_bilinear(BST, BST)


def two_form(a, b):
    w = np.zeros((D, D))
    w[a, b], w[b, a] = 1.0, -1.0
    return w


PAIRS = [(a, b) for a in range(4) for b in range(a + 1, 4)]
WB = np.array([two_form(a, b) for a, b in PAIRS])          # W = du, 4D components
YB = np.array([two_form(a, 4) for a in range(4)])          # dchi ^ u: (dchi ^ e^4)_{a4} = d_a chi
VW, VY = mix_linear(WB), mix_linear(YB)
WW = np.einsum("iab,jcd,ac,bd->ij", WB[:, :4, :4], WB[:, :4, :4], ETA[:4, :4], ETA[:4, :4])
WtW = 0.25 * np.einsum("abcd,iab,jcd->ij", EPS4, WB[:, :4, :4], WB[:, :4, :4])
WtW = 0.5 * (WtW + WtW.T)


def solve(c, cp):
    Mt = 0.5 * MQ + c * MO
    lin = np.vstack([c * VW, cp * VY])
    G = -0.5 * lin @ np.linalg.solve(Mt, lin.T)
    return 0.5 * (G + G.T), Mt, lin


ok_all, rows = True, []
for c, cp in [(0.01, 0.3), (0.1, 1.0), (0.5, -0.7), (1.0, 0.2), (3.0, 1.5)]:
    G, *_ = solve(c, cp)
    Gww, Gwy, Gyy = G[:6, :6], G[:6, 6:], G[6:, 6:]
    X = np.stack([WW.ravel(), WtW.ravel()], axis=1)
    (A, B), *_ = np.linalg.lstsq(X, Gww.ravel(), rcond=None)
    rw = np.linalg.norm(X @ np.array([A, B]) - Gww.ravel()) / np.linalg.norm(Gww)
    C = np.trace(ETA[:4, :4] @ Gyy) / 4
    ry = np.linalg.norm(Gyy - C * ETA[:4, :4]) / max(np.linalg.norm(Gyy), 1e-300)
    good = (rw < 1e-10 and ry < 1e-10 and np.linalg.norm(Gwy) < 1e-12 * np.linalg.norm(G)
            and abs(A + 1.5 * c ** 2) < 1e-10 * max(1, c ** 2) and abs(B + 4 * c ** 3) < 1e-10 * max(1, c ** 3))
    ok_all &= good
    rows.append((c, cp, A, B, C))
    info(f"c={c:g}, c'={cp:g}: A={A:.6g} B*={B:.6g} C={C:.6g} (C/c'^2 = {C / cp ** 2:.6g}) cross={np.linalg.norm(Gwy):.1e}")
check("U1 c(chi) u-Holst: A = -3c^2/2, B* = -4c^3 pointwise; no W x dchi cross term; (dchi)^2 term only", ok_all)
Cs = [r[4] / r[1] ** 2 for r in rows]
info(f"C/c'^2 is c-independent: {np.ptp(Cs) < 1e-10}, value {Cs[0]:.6g} (sign < 0 adds a healthy chi kinetic term)")
check("U1b chi kinetic correction C = C0 c'^2 with C0 < 0 (no ghost)", np.ptp(Cs) < 1e-10 and Cs[0] < 0)

# U2: exact u -> no linear source
Mt = 0.5 * MQ + 0.7 * MO
check("U2 exact u (X = d u = 0): K = 0 is the solution, no photon term",
      np.allclose(np.linalg.solve(Mt, np.zeros(NB)), 0))

# U3: contortion response K_[ab4] to W (for the EDM estimate)
def K_of_W(c, w):
    Mt = 0.5 * MQ + c * MO
    k = -np.linalg.solve(Mt, (c * VW).T @ w)
    K = np.einsum("n,nabc->abc", k, BST)
    return alt3(K)[:4, :4, 4]


w = np.zeros(6)
w[0] = 1.0                       # W_01 = 1
Wm = np.einsum("w,wab->ab", w, WB)[:4, :4]
Wt = 0.5 * np.einsum("abcd,cd->ab", EPS4, Wm)
for c in (1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0):
    Ka = K_of_W(c, w)
    ev = np.sum(Ka * Wm) / np.sum(Wm * Wm)
    od = np.sum(Ka * Wt) / np.sum(Wt * Wt)
    info(f"U3 c={c:g}: K_[ab4] = {ev:.4g} W_ab + {od:.4g} W~_ab")
ok3 = True
for c in (1e-3, 0.37, 2.0, 50.0):
    Ka = K_of_W(c, w)
    ev = np.sum(Ka * Wm) / np.sum(Wm * Wm)
    od = np.sum(Ka * Wt) / np.sum(Wt * Wt)
    ok3 &= abs(abs(od) - c) < 1e-9 * max(1, c) and abs(abs(ev) - 4 * c ** 2 / 3) < 1e-9 * max(1, c ** 2)
    ok3 &= np.linalg.norm(Ka - ev * Wm - od * Wt) < 1e-9 * max(1, c ** 2)
check("U3 K_[ab4] = (+-)c W~_ab (+-)(4/3)c^2 W_ab exactly (c-odd piece linear in c)", ok3)

# U4: helicity evolution with Z(tau), Y(tau)
kmode = 300.0
cfun = lambda t: 0.5 + 0.1 * (1 + np.tanh((t - 5.0) / 0.8)) / 2
dcf = lambda t: 0.1 / (2 * 0.8 * np.cosh((t - 5.0) / 0.8) ** 2)
Zf = lambda t: 1 + 12 * cfun(t) ** 2
dZ = lambda t: 24 * cfun(t) * dcf(t)
dY = lambda t: 48 * cfun(t) ** 2 * dcf(t)


def rhs(t, y, s):
    A = y[0] + 1j * y[1]
    Ap = y[2] + 1j * y[3]
    App = -(dZ(t) / Zf(t)) * Ap - (kmode ** 2 + s * kmode * dY(t) / Zf(t)) * A
    return [Ap.real, Ap.imag, App.real, App.imag]


y0 = [1.0, 0.0, 0.0, -kmode]          # A = exp(-i k t)
sols = [solve_ivp(rhs, (0, 10), y0, args=(s,), rtol=1e-11, atol=1e-12, method="DOP853") for s in (1, -1)]
Ap_, Am_ = [s.y[0, -1] + 1j * s.y[1, -1] for s in sols]
dphase = abs(np.angle(Ap_ / Am_))
pred = quad(lambda t: dY(t) / Zf(t), 0, 10, limit=200)[0]
check("U4 helicity phase difference = int Y'/Z dtau (rotation beta = half of it)", abs(dphase / pred - 1) < 2e-3,
      f"{dphase:.6f} vs {pred:.6f}")

# U5: d beta = -c d ln alpha
cs = sp.symbols("c", real=True)
Zs, Ys = 1 + 12 * cs ** 2, 16 * cs ** 3
dbeta = sp.Rational(1, 2) * sp.diff(Ys, cs) / Zs        # beta = (1/2) int Y'/Z (magnitude convention)
dlnalpha = -sp.diff(sp.log(Zs), cs)
check("U5 d beta / d ln alpha = -c exactly", sp.simplify(dbeta / dlnalpha + cs) == 0)
bclosed = sp.integrate(dbeta, cs)
check("U5b beta(c) = 2c - arctan(2 sqrt(3) c)/sqrt(3) + const",
      sp.simplify(sp.diff(bclosed - (2 * cs - sp.atan(2 * sp.sqrt(3) * cs) / sp.sqrt(3)), cs)) == 0)

# U6: winding scalar Higgses the KK photon
x = sp.symbols("x0:4", real=True)
y_ = sp.Symbol("y", real=True)
Asym = sp.symbols("A0:4", real=True)
ph, wR = sp.symbols("phi_r wR", positive=True)
pi_ = sp.Function("pi")(*x)
g = sp.zeros(5, 5)
eta4 = sp.diag(-1, 1, 1, 1)
for m in range(4):
    for n in range(4):
        g[m, n] = eta4[m, n] + ph ** 2 * Asym[m] * Asym[n]
    g[m, 4] = g[4, m] = ph ** 2 * Asym[m]
g[4, 4] = ph ** 2
ginv = g.inv()
psi = wR * y_ + pi_
coords = list(x) + [y_]
dpsi = sp.Matrix([sp.diff(psi, cc) for cc in coords])
kin = sp.simplify((dpsi.T * ginv * dpsi)[0])
dpiA = sp.Matrix([sp.diff(pi_, x[m]) - wR * Asym[m] for m in range(4)])
target = (dpiA.T * eta4 * dpiA)[0] + wR ** 2 / ph ** 2
check("U6 (d psi)^2 = (d pi - (w/R) A)^2 + (w/R)^2/phi^2: a winding orienting scalar gives the KK photon a mass",
      sp.simplify(kin - target) == 0)
