"""Goal 4: anisotropic rotation from torsion perturbations, order of magnitude, vs SPT-3G A_CB < 1.2e-4 rad^2.

A homogeneous p_0hat gives no anisotropy at linear order. If the torsion tracks the LOCAL expansion, perturbations of the
expansion rate rotate different lines of sight differently. Two assumptions (the theory, task 022, decides which if any):
  (A) kinematic: p_0hat = xi eps theta/3 with theta the local expansion of the matter flow. Linear continuity:
      delta theta / 3 = -f H delta_m / 3   (f = dlnD/dlna).
  (B) potential: delta p / p = Phi (Newtonian potential), Phi = -(3/2) Omega_m H0^2 (1+z) D delta_0 / k^2.
delta beta(n) = 1/2 int dt delta p_0hat = 1/2 int dchi a delta p_0hat (comoving chi). Limber:
  C_L = int dchi W(chi, k)^2 / chi^2 P_0(k = (L+1/2)/chi),  D_L = L(L+1)C_L/2pi  (A_CB is D_L for a scale-invariant
  spectrum, src/templates.py). Linear P(k): Eisenstein-Hu 1998 no-wiggle transfer, normalised with A_s, n_s from data/.
Shape: p_0hat = xi eps H (n = 0) with xi eps at the ACT best fit. Linear theory only; non-linear scales noted in README.
"""
import math

import numpy as np
from scipy.integrate import quad

from ec06_common import PL18, ZSTAR, OMM, OMB, h, E, DEG, BETA_ACT, MEAS, beta_unit, ec6

checks = []


def check(name, ok, info=""):
    checks.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}: {name} -- {info}")


C_KMS = 299792.458
H0_MPC = 100 * h / C_KMS                       # H0/c in Mpc^-1
ombh2, omch2 = float(PL18["ombh2"]["value"]), float(PL18["omch2"]["value"])
omh2 = OMM * h**2
As = math.exp(float(PL18["logA"]["value"])) * 1e-10
ns = float(PL18["ns"]["value"])
from ec06_common import c008  # noqa: E402
T27 = float(c008.INP["T_CMB"]["value"]) / 2.7   # EH98 Theta_2.7
fb = ombh2 / omh2


def T_eh_nowiggle(k):
    """EH98 Eqs. (26)-(31); k in Mpc^-1."""
    s = 44.5 * np.log(9.83 / omh2) / np.sqrt(1 + 10 * ombh2 ** 0.75)
    aG = 1 - 0.328 * np.log(431 * omh2) * fb + 0.38 * np.log(22.3 * omh2) * fb**2
    Gam = OMM * h * (aG + (1 - aG) / (1 + (0.43 * k * s) ** 4))
    q = (k / h) * T27**2 / Gam
    L0 = np.log(2 * np.e + 1.8 * q)
    C0 = 14.2 + 731 / (1 + 62.5 * q)
    return L0 / (L0 + C0 * q**2)


def E_late(a):
    return np.sqrt(OMM / a**3 + (1 - OMM))


def D_unnorm(a):
    """Linear growth for matter + Lambda, D -> a deep in matter era: D = (5 Om/2) E(a) int_0^a da'/(a' E)^3."""
    a = np.atleast_1d(a)
    return np.array([2.5 * OMM * E_late(x) * quad(lambda y: 1 / (y * E_late(y)) ** 3, 0, x, epsrel=1e-10)[0] for x in a])


D0 = D_unnorm(1.0)[0]


def P0(k):
    """Linear matter P(k) today, Mpc^3."""
    D2R = As * (k / 0.05) ** (ns - 1)
    D2d = (4.0 / 25.0) * D2R * (k / H0_MPC) ** 4 * T_eh_nowiggle(k) ** 2 * (D0 / OMM) ** 2
    return 2 * np.pi**2 / k**3 * D2d


def sigmaR(R):
    W = lambda x: 3 * (np.sin(x) - x * np.cos(x)) / x**3
    lk = np.linspace(np.log(1e-5), np.log(1e2), 20001)
    k = np.exp(lk)
    return math.sqrt(np.trapezoid(k**3 * P0(k) / (2 * np.pi**2) * W(k * R) ** 2, lk))


s8 = sigmaR(8.0 / h)
s8_pl = ec6("sigma8")
check("A1 sigma8 from A_s, n_s with EH98 no-wiggle P(k) vs Planck 0.8111 (no massive nu, no wiggles: 5%)",
      abs(s8 / s8_pl - 1) < 0.05, f"{s8:.4f}")
fD = D_unnorm(np.array([0.5]))[0]
check("A2 growth D(a=0.5)/D(1) for Om = 0.315 between 0.6 and 0.65 (standard LCDM value ~0.61)",
      0.6 < fD / D0 < 0.65, f"{fD / D0:.4f}")

# line-of-sight tables (late-time background without radiation is enough for z < 50; radiation included in E for chi)
zz = np.concatenate([np.geomspace(1e-3, 50, 3000), np.geomspace(50.01, ZSTAR, 300)])
aa = 1 / (1 + zz)
chi = np.array([quad(lambda z: 1 / (H0_MPC * E(1 / (1 + z))), 0, z, epsrel=1e-9)[0] for z in zz])
Dz = D_unnorm(aa) / D0
lnD = np.log(Dz)
fz = np.gradient(lnD, np.log(aa))
Hz = H0_MPC * E(aa)

K = beta_unit(ZSTAR, 0.0)
xieps = BETA_ACT * DEG / K


def cl_limber(L, model):
    k = (L + 0.5) / chi
    if model == "A":
        W = 0.5 * xieps * aa * Hz * fz * Dz / 3.0
        integrand = W**2 / chi**2 * P0(k)
    else:
        W = 0.5 * xieps * aa * Hz * 1.5 * OMM * H0_MPC**2 * (1 + zz) * Dz / k**2
        integrand = W**2 / chi**2 * P0(k)
    return np.trapezoid(integrand, chi)


Ls = np.unique(np.geomspace(2, 3000, 60).astype(int))
acb_lim = float(MEAS["BAL25-ACB-SPT3G"]["value"])
res = {}
for model in ("A", "B"):
    DL = np.array([L * (L + 1) * cl_limber(L, model) / (2 * np.pi) for L in Ls])
    res[model] = DL
    i = DL.argmax()
    print(f"model {model}: max D_L = {DL[i]:.2e} rad^2 at L = {Ls[i]}; D_L(L=100) = {np.interp(100, Ls, DL):.2e}; "
          f"ratio to SPT A_CB limit {acb_lim:.1e}: {DL[i] / acb_lim:.1e}")
# analytic cross-check of model A: delta beta_rms^2 ~ (1/2 xi eps)^2 int dchi (aHf/3)^2 sigma^2(R) * 2R (coherent cells)
R = 10.0 / h
sig = sigmaR(R)
m = zz < 3
rough = (0.5 * xieps) ** 2 * np.trapezoid((aa[m] * Hz[m] * fz[m] * Dz[m] / 3) ** 2 * sig**2 * 2 * R, chi[m])
DLmaxA = res["A"].max()
var_limber = np.trapezoid(res["A"] / Ls, Ls)   # sum (2L+1) C_L / 4pi ~ int D_L dL / L
check("A3 Limber variance (model A, L = 2-3000) within 2 orders of magnitude of the coherent-cell estimate",
      1e-2 < var_limber / rough < 1e2, f"cell estimate (delta beta)^2 = {rough:.2e}, Limber variance = {var_limber:.2e}")
check("A4 both assumptions give A_CB far below the SPT-3G limit (factor > 1e3)",
      max(res["A"].max(), res["B"].max()) < 1e-3 * acb_lim, f"max {max(res['A'].max(), res['B'].max()):.2e}")
# super-horizon / Sachs-Wolfe-type estimate: delta beta ~ beta * Phi
phi_rms = 0.6 * math.sqrt(As)  # matter era Phi = (3/5) R, rms per ln k
print(f"Sachs-Wolfe-type: (beta * Phi)^2 ~ ({BETA_ACT * DEG:.2e} * {phi_rms:.1e})^2 = {(BETA_ACT * DEG * phi_rms) ** 2:.1e} rad^2")
np.savetxt(__file__.replace("ec06_aniso.py", "acb_estimate.csv"), np.column_stack([Ls, res["A"], res["B"]]),
           delimiter=",", header="L,D_L_modelA_kinematic_rad2,D_L_modelB_potential_rad2", comments="", fmt="%.6g")
print(f"\n{sum(checks)} PASS, {len(checks) - sum(checks)} FAIL")
