"""Skeptic 011: does a derivative (shift-symmetric) coupling d(phi).J5/f give a QCD mass or a
low-energy photon coupling? Independent of task 007/008 code.

Part A (QCD mass). Two-flavour U(2) chiral model, U = diag(e^{i p_u}, e^{i p_d}),
  V = -L (m_u cos p_u + m_d cos p_d) + (chi/2)(p_u + p_d - theta_G)^2
(the chi term is the GG~ / eta' anomaly piece; theta_G is the GG~ angle).
  Reading B (correct): the field-dependent axial rotation q -> exp(i c phi gamma5/2f) q that removes
     the derivative coupling moves phi into BOTH the quark mass phases (arg m_q -> arg m_q - c phi/f)
     AND theta_G (theta_G -> theta_G + 2 c phi/f). Minimise V over (p_u, p_d) at fixed phi.
  Reading A (008 'as published'): keep only the anomaly piece theta_G = 2 c phi/f, mass phases untouched.
A falsifiable check: reading A must give the textbook axion mass; reading B must give V_min(phi) flat.

Part B (photon). Derivative coupling to a fermion of mass m: amplitude = anomaly (1) minus
mass-insertion triangle 2 int dx dy /(1 - r x y), r = q^2/m^2. Computed by direct numerical
integration of the Feynman-parameter integral (not from the closed form) and compared with
B1 = 1 - tau arcsin^2(1/sqrt(tau)), tau = 4/r, and with the limit -r/12.
"""
import numpy as np
import sympy as sp
from scipy.integrate import dblquad
from scipy.optimize import minimize


def report(name, ok, det=""):
    print(("PASS" if ok else "FAIL") + ": " + name + (" -- " + det if det else ""))


# ---------------- Part A -------------------------------------------------------------
L, mu, md, chi, c, f = 1.0, 1.0, 2.0, 50.0, 1.0, 10.0      # toy units


def V(p, thG, ph_u=0.0, ph_d=0.0):
    pu, pd = p
    return -L * (mu * np.cos(pu + ph_u) + md * np.cos(pd + ph_d)) + 0.5 * chi * (pu + pd - thG) ** 2


def Vmin(phi, reading):
    if reading == "B":
        a = c * phi / f
        args = (2 * a, -a, -a)        # theta_G shift 2a (two flavours x a), mass phases -a each
    else:
        args = (2 * c * phi / f, 0.0, 0.0)
    r = minimize(V, x0=[0.0, 0.0], args=args, method="BFGS", options={"gtol": 1e-12})
    return r.fun


phis = np.linspace(-0.3, 0.3, 13) * f
VB = np.array([Vmin(p, "B") for p in phis])
VA = np.array([Vmin(p, "A") for p in phis])
report("reading B (full chiral rotation): V_min(phi) flat => m_phi = 0",
       np.ptp(VB) < 1e-9, f"ptp V_B = {np.ptp(VB):.1e}")
# numeric curvature of reading A vs textbook  m^2 f^2 = (2c)^2 L mu md/(mu+md)  (chi -> inf, 2 flavours)
hstep = 1e-3 * f
curvA = (Vmin(hstep, "A") - 2 * Vmin(0, "A") + Vmin(-hstep, "A")) / hstep**2
pred = (2 * c / f) ** 2 * L * mu * md / (mu + md) / (1 + L * mu * md / (mu + md) / chi)
report("reading A reproduces the textbook QCD-axion mass (so the test can fail)",
       abs(curvA / pred - 1) < 1e-3, f"m^2 numeric {curvA:.5e}, formula {pred:.5e}")
report("reading A is not flat", np.ptp(VA) > 1e-3, f"ptp V_A = {np.ptp(VA):.2e}")

# symbolic: the invariant combination theta_bar = theta_G - sum arg(m) is phi independent
a_ = sp.symbols("a")
thG, pu, pd = sp.symbols("theta_G p_u p_d")
theta_bar_B = (thG + 2 * a_) + ((pu - a_) + (pd - a_))      # after rotation
report("theta_bar = theta_G + arg det M is phi-independent under the rotation (sympy)",
       sp.diff(theta_bar_B, a_) == 0)

# kinetic mixing route (derivative coupling to pi0 after hadronisation): mass eigenvalue of phi
lam = sp.Symbol("lambda")
eps, mpi = sp.symbols("epsilon m_pi", positive=True)
K = sp.Matrix([[1, eps], [eps, 1]])
Mm = sp.Matrix([[mpi**2, 0], [0, 0]])
roots = sp.solve(sp.det(Mm - lam * K), lam)
report("phi-pi0 kinetic mixing alone leaves a massless eigenstate", 0 in [sp.simplify(r) for r in roots],
       f"roots {roots}")

# ---------------- Part B -------------------------------------------------------------
def tri(r):
    val, _ = dblquad(lambda y, x: 1.0 / (1.0 - r * x * y), 0, 1, 0, lambda x: 1 - x,
                     epsabs=1e-13, epsrel=1e-12)
    return 2 * val


ok = True
rows = []
for r in (1e-3, 0.1, 1.0, 3.0):
    Bnum = 1 - tri(r)
    tau = 4 / r
    Bcl = 1 - tau * np.arcsin(1 / np.sqrt(tau)) ** 2
    rows.append((r, Bnum, Bcl, -r / 12))
    ok &= abs(Bnum - Bcl) < 1e-9 * max(1, abs(Bcl)) + 1e-12
print("r = q^2/m^2, B1 (Feynman-parameter integral), B1 (closed form), -r/12")
for row in rows:
    print("  %.0e  %.6e  %.6e  %.6e" % row)
report("B1 from the Feynman-parameter integral = closed form", ok)
report("B1 -> -r/12 for r -> 0 (heavy-fermion decoupling)", abs(rows[0][1] / rows[0][3] - 1) < 1e-3)
report("massless/light fermion limit is NOT suppressed (r = 3: |B1| = O(1))", abs(rows[3][1]) > 0.3,
       f"B1(3) = {rows[3][1]:.3f}")
