"""EC-06 numbers per theory and the background table for the cosmology agent (task 022).

beta = +1/2 int p_tau dtau (CONVENTIONS.md) with p_mu = xi X_mu, X = S (axial) or T (trace).
Define eps_X(z) = X_t / H, X_t the LOWER cosmic-time component (X_tau dtau = X_t dt). Then
    beta = (xi/2) int_0^{z*} eps_X(z) dz/(1+z).
For the axial vector S_t = -S^0hat (S^0hat = comoving-observer upper component used in CLAIMS EC-01).
All measured numbers are read from data/ (and the derived lab bound from work/20261120-005-bounds).
"""
import csv
import math
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq

from common import ROOT, check, load, summary

HERE = os.path.dirname(os.path.abspath(__file__))
cos = {k: float(r["value"]) for k, r in load("cosmo_inputs.csv").items() if r["value"] not in ("",)}
cos_err = {k: float(r["err"]) for k, r in load("cosmo_inputs.csv").items() if r["err"] not in ("",)}
bir = load("birefringence_measurements.csv")
e6c = load("ec06_cosmo_inputs.csv")
e6t = load("ec06_theory_inputs.csv")
bnd = load("bounds.csv")
with open(os.path.join(ROOT, "work", "20261120-005-bounds", "derived_bounds.csv"), newline="", encoding="utf-8") as f:
    derived = {r["id"]: r for r in csv.DictReader(f)}

beta_deg = float(bir["DPK25-beta-baseline"]["value"])
beta = math.radians(beta_deg)
zs = cos["z_star"]
zre = float(e6c["z_re"]["value"])
h = cos["H0"] / 100
H0 = cos["H0"] / cos["Mpc_km"] * cos["hbar_eV_s"]  # eV
Om = cos["Omega_m"]
Or = 2.469e-5 / h**2 * (cos["T_cmb_K"] / 2.7255) ** 4 * (1 + 0.2271 * cos["N_eff_SM"])
OL = 1 - Om - Or
E = lambda z: math.sqrt(Om * (1 + z) ** 3 + Or * (1 + z) ** 4 + OL)
Hs = H0 * E(zs)
L1 = math.log(1 + zs)
print(f"inputs: beta = {beta_deg} deg = {beta:.3e} rad, z* = {zs}, z_re = {zre}, H0 = {H0:.3e} eV, H* = {Hs:.3e} eV")


def beta_int(eps, z1=0.0, z2=zs):
    return 0.5 * quad(lambda z: eps(z) / (1 + z), z1, z2, limit=200)[0]


rows = []

# ---------------------------------------------------------------- symbolic: closed-p beta = endpoint difference
t = sp.symbols("t")
th = sp.Function("theta")(t)
xi_s = sp.symbols("xi")
check("N0 closed p = xi d theta: beta = (xi/2)[theta_0 - theta_*] (endpoint only, achromatic)",
      sp.simplify(sp.integrate(xi_s * sp.diff(th, t) / 2, (t, sp.Symbol("ts"), sp.Symbol("t0")))
                  - xi_s / 2 * (th.subs(t, sp.Symbol("t0")) - th.subs(t, sp.Symbol("ts")))) == 0)
# SNY Eq. 46: Rdot = -(2/3)(R + 6 mu/b) Phi  =>  Phi = -(3/2) d/dt ln(R + 6 mu/b): T_t is a gradient on FRW
R = sp.Function("R")(t)
mu, b = sp.symbols("mu b", positive=True)
Phi = sp.solve(sp.Eq(sp.diff(R, t), -sp.Rational(2, 3) * (R + 6 * mu / b) * sp.Symbol("Phi")), sp.Symbol("Phi"))[0]
check("N0b SNY 0+ mode: Phi = -(3/2) d ln(R + 6 mu/b)/dt", sp.simplify(Phi + sp.Rational(3, 2) * sp.diff(sp.log(R + 6 * mu / b), t)) == 0)

# ---------------------------------------------------------------- (b1) f(T) / TEGR flat: trace T_t = 3H
eps_T = lambda z: 3.0
bT = beta_int(eps_T)
xi_T = beta / bT
fr_re = beta_int(eps_T, 0, zre) / bT
fr_10 = beta_int(eps_T, 0, 10) / bT
print(f"(b1) teleparallel flat, CFJ on trace: beta = {bT:.3f} xi rad [(3/2)ln(1+z*) = {1.5*L1:.3f}];"
      f" xi = {xi_T:.3e}; fraction after z_re={zre}: {fr_re:.3f}, after z=10: {fr_10:.3f}")
check("N1 teleparallel trace: beta = (3 xi/2) ln(1+z*)", abs(bT - 1.5 * L1) < 1e-6)
rows.append(["teleparallel f(T) flat (TEGR: frame dependent)", "trace T_t", "3 (constant)", 3.0, 3.0,
             f"xi = {xi_T:.2e}; beta_reion/beta = {fr_re:.2f} (z_re) {fr_10:.2f} (z=10)",
             "check_teleparallel T1; f(T) good tetrad T4"])

# ---------------------------------------------------------------- (b2) teleparallel k = +1 axial branch
OK = float(e6t["Omega_K"]["value"]); OKe = float(e6t["Omega_K"]["err"])
OK_closed_max = abs(OK - 1.96 * OKe)  # largest |Omega_K| with Omega_K < 0 at 95%
eps_k = lambda z, ok=OK_closed_max: 2 * math.sqrt(ok) * (1 + z) / E(z)
bk_unit = beta_int(lambda z: eps_k(z, 1.0))  # beta / (xi sqrt|Omega_K|)
chi_star = quad(lambda z: 1 / E(z), 0, zs, limit=200)[0]
check("N2 k=+1 axial: beta = xi sqrt|Omega_K| H0 chi(z*) (conformal-distance weighted)",
      abs(bk_unit - chi_star) / chi_star < 1e-6, f"(H0 chi* = {chi_star:.3f})")
xi_k = beta / (math.sqrt(OK_closed_max) * chi_star)
fr_k = quad(lambda z: 1 / E(z), 0, zre)[0] / chi_star
fr_k10 = quad(lambda z: 1 / E(z), 0, 10)[0] / chi_star
print(f"(b2) teleparallel k=+1 axial: eps(z) = +-2 sqrt|Omega_K| (1+z)/E(z); |Omega_K| <= {OK_closed_max:.4f} (95%, closed);"
      f" eps_0 <= {eps_k(0):.3f}, eps(z*) <= {eps_k(zs):.2e}; xi >= {xi_k:.3e} (xi sqrt|Omega_K| = {beta/chi_star:.2e});"
      f" fraction after z_re: {fr_k:.3f}, after z=10: {fr_k10:.3f}")
rows.append(["teleparallel f(T) closed k=+1 (axial branch)", "axial S_t (sign = branch)",
             "+-2 sqrt|Omega_K| (1+z)/E(z)", eps_k(0), eps_k(zs),
             f"xi >= {xi_k:.2e} at |Omega_K| = {OK_closed_max:.4f}; beta_reion/beta = {fr_k:.2f} (z_re) {fr_k10:.2f} (z=10)",
             "check_teleparallel T3; Omega_K from data/ec06_theory_inputs.csv"])
rows.append(["teleparallel f(T) closed k=+1 (axial branch)", "trace T_t", "3 (constant)", 3.0, 3.0,
             "as flat case", "check_teleparallel T3b"])

# ---------------------------------------------------------------- (a) Magueijo-Zlosnik tracker
def Rc_of(gam):
    # MZ Eq. 66 in u = sqrt(R_c) (robust for tiny gamma, R_c ~ gamma^2/9)
    f = lambda u: 2 * gam * u**4 + gam * u**2 + 2 * math.sqrt(u**2 + 1) * u**3 + 3 * math.sqrt(u**2 + 1) * u - gam
    u = brentq(f, 0.0, math.sqrt(0.5), xtol=1e-300, rtol=1e-14)
    return u * u


def mz(gam):
    Rc = Rc_of(gam)
    x = math.sqrt(Rc / (1 + Rc))
    Tg = gam * x * (1 - x**2 - 2 * x / gam) / (1 - x**2 + 2 * gam * x)
    epsS = 2 * x / (1 - Tg)          # S_t / H = 2P/H, H = g - T
    epsT = -3 * Tg / (1 - Tg)        # T_t / H = -3T/H
    ratio = (1 - Tg) ** 2 / (1 - x**2)  # H^2 / (kappa rho/3)
    return Rc, x, Tg, epsS, epsT, ratio


# N_eff: tracker changes H^2 by 'ratio' in the radiation era; translate to Delta N_eff (photons + nu)
Neff_lo = cos["N_eff_obs"] - 2 * cos_err["N_eff_obs"]
f_nu = 0.2271
dmax = f_nu * (cos["N_eff_SM"] - Neff_lo) / (1 + f_nu * cos["N_eff_SM"])
gmax = brentq(lambda gm: (1 - mz(gm)[5]) - dmax, 1e-4, 50)
Rc, x, Tg, eS, eT, ratio = mz(gmax)
for gm in (0.05, 0.1, 0.3, gmax, 1.0, 1e4):
    r_ = mz(gm)
    print(f"   MZ tracker gamma = {gm:.3g}: R_c = {r_[0]:.3e}, x = P/g = {r_[1]:.3f}, T/g = {r_[2]:.3f},"
          f" eps_S = {r_[3]:.3f}, eps_T = {r_[4]:.3f}, H^2/(kappa rho/3) = {r_[5]:.3f}")
inf = mz(1e6)
check("N3a MZ gamma->inf tracker: H^2/(kappa rho/3) -> 2/3 (33% deficit, excluded by N_eff)", abs(inf[5] - 2 / 3) < 1e-3)
xi_mz = beta / (eS / 2 * L1)
print(f"(a) MZ finite-gamma tracker: N_eff >= {Neff_lo:.2f} (95%) => |1 - H^2/(kappa rho/3)| <= {dmax:.3f}"
      f" => gamma <= {gmax:.3f}, eps_S <= {eS:.3f}, eps_T = {eT:.3f};"
      f" if the tracker persists to z=0: xi_S >= {xi_mz:.2e}; beta_reion/beta = {math.log(1+zre)/L1:.2f}")
check("N3b N_eff alone would allow eps_S = O(0.1-1) on the MZ tracker (superseded by the GW bound N3d)", 1e-2 < eS < 2)
# GW170817 tensor speed in the same theory (2001.06373): tracker P^2 = gamma^2 kappa rho/27 (P/H = gamma/3, agrees
# with Rc ~ gamma^2/9 above); Eq. 72 omega^2 = [1 + (1+3w) gamma^2/9] k^2 -> |c_T - 1| = gamma^2/18 (dust);
# Lambda era Eq. 65 (rho = 0): omega^2 = [1 - 8 sqrt3 gamma eps_P + 84 eps_P^2] k^2, eps_P = P/sqrt(Lambda).
dv = max(abs(float(e6t["cT_GW170817_lo"]["value"])), abs(float(e6t["cT_GW170817_hi"]["value"])))
dv_hi = abs(float(e6t["cT_GW170817_hi"]["value"]))
check("N3c MZ tracker x = P/g -> gamma/3 at small gamma (2001.06373: P^2 = gamma^2 kappa rho/27)",
      abs(mz(1e-3)[1] / (1e-3 / 3) - 1) < 1e-2)
g_gw = math.sqrt(18 * dv)
eS_gw = mz(g_gw)[3]
epsP_max = math.sqrt(dv_hi / 42)
eS0_gw = 2 * epsP_max * math.sqrt(3 * OL)  # 2P/H0 with sqrt(Lambda) = sqrt(3 Omega_L) H0
g_paper = math.sqrt(1e-15)
print(f"(a) GW170817 |dv/v| <= {dv:.0e}: tracker gamma <= {g_gw:.1e} -> eps_S <= {eS_gw:.1e} -> xi_S >= {beta/(eS_gw/2*L1):.1e};"
      f" with the paper's own gamma^2 <~ 1e-15 (Eq. 68): eps_S ~ {2*g_paper/3:.1e}, xi >~ {beta/(g_paper/3*L1):.1e};"
      f" Lambda era: eps_P <= {epsP_max:.1e} -> eps_S(today) <= {eS0_gw:.1e}; |gamma| << 1e-21 (Eq. 70) -> xi >> 1e18")
check("N3d MZ finite gamma: GW speed forces eps_S <~ 1e-7, so xi >~ 1e4 (dead as a perturbative coupling)",
      beta / (eS_gw / 2 * L1) > 1e3)
rows.append(["Magueijo-Zlosnik EC + varying Lambda, finite-gamma tracker", "axial S_t (sign = initial P)",
             "constant on tracker (R_c from MZ Eq. 66), P ~ a^-3 in Lambda era", eS_gw, eS_gw,
             f"N_eff alone: gamma <= {gmax:.2f}, xi_S >= {xi_mz:.2e}; GW170817 (2001.06373): gamma <= {g_gw:.1e}, eps_S <= {eS_gw:.1e}, xi >= {beta/(eS_gw/2*L1):.0e}: DEAD",
             "check_mz M1-M5; 1908.05184 Eqs. 26-29, 57-66"])
rows.append(["Magueijo-Zlosnik EC + varying Lambda, finite-gamma tracker", "trace T_t",
             "constant on tracker", mz(g_gw)[4], mz(g_gw)[4], f"at GW-allowed gamma = {g_gw:.1e}", "check_mz M3"])
rows.append(["Magueijo-Zlosnik, gamma -> infinity (Euler term only)", "axial + trace", "x^2 = 1/3: eps_S = 4x/(1+x^2), eps_T = -3(1-x^2)/(1+x^2)",
             4 * math.sqrt(1 / 3) / (4 / 3), 4 * math.sqrt(1 / 3) / (4 / 3),
             "DEAD background: tracker 33% H^2 deficit; MZ call it unviable; tensor modes undetermined (2001.06373)",
             "check_mz M2; N3a"])

# ---------------------------------------------------------------- (c) PGT
rows.append(["PGT 0+ scalar torsion (Shie-Nester-Yo 0805.3834)", "trace T_t = Phi",
             "Phi = -(3/2) d ln(R+6mu/b)/dt; damped oscillation", "0.8 to 1.4 (SNY cases III-V, Phi/H at t0)", "not computed (SNY early-time solutions unreliable)",
             f"if eps ~ 1 constant: xi ~ {2*beta/L1:.1e}; sign of Phi vs our T_t not verified",
             "N0b; SNY Eqs. 43-47, Table 1"])
rows.append(["PGT 0- pseudoscalar torsion (EC-04; Chen et al. 0908.3323; BHN 1009.5112)", "axial S = d phi/sqrt(A)",
             "ALP field: thawing, m ~ H0 needed", "free", "free",
             "see CLAIMS EC-04: ~1e120 tuning of the E^2 coefficient", "work/20261120-007"])


# ---------------------------------------------------------------- prior-art forms (task 024): p = c d ln(scalar)
# Yu-Yang-Qiu 2508.10532: ln Q, Q = 6H^2 (= f(T) torsion scalar on flat FRW, check_teleparallel T4);
# Li-Xia-Li hep-ph/0611192: ln R, R = 6(Hdot + 2H^2). In matter domination d ln H^2 = d ln R = -3 d ln a,
# so p = xi T = 3 xi d ln a equals p = -xi d ln Q there. beta per unit coefficient c (p = c d ln X): (c/2) ln[X0/X*].
def Rsc(z):
    Hz2 = E(z) ** 2
    dlnE2 = (3 * Om * (1 + z) ** 3 + 4 * Or * (1 + z) ** 4) / Hz2  # -d ln E^2 / d ln a
    return Hz2 * (2 - 0.5 * dlnE2)  # (Hdot + 2H^2)/H0^2 = H^2 (2 + dlnH/dlna)
bQ = 0.5 * math.log(E(0) ** 2 / E(zs) ** 2)
bR = 0.5 * math.log(Rsc(0) / Rsc(zs))
print(f"prior-art forms, beta per unit c: ln Q: {bQ:.3f}, ln R: {bR:.3f}; trace torsion p = xi T: (3/2)ln(1+z*) = {1.5*L1:.3f}"
      f" (= -ln Q form x {-1.5*L1/bQ:.3f}); post-reion share ln Q: {math.log(E(zre)**2)/math.log(E(zs)**2):.2f}")
check("N7 flat f(T) trace CFJ equals the ln Q (Q = 6H^2 = torsion scalar) coupling to within 5% in beta (matter-era identical)",
      abs(-1.5 * L1 / bQ - 1) < 0.05)

# ---------------------------------------------------------------- bounds
S0lab = float(derived["S0-boost-neutron"]["value_eV"])
kaf = float(bnd["P-kAF-CFJ90"]["value_eV"].replace("<", "").strip())
th_rg = float(e6c["RG_UV_theta_deg"]["value"]); th_rg_e = float(e6c["RG_UV_theta_deg"]["err"])
zrg = float(e6c["RG_UV_zmin"]["value"])
print(f"bounds: lab |S_0| < {S0lab:.1e} eV vs O(H0) torsion {H0:.2e} eV: ratio {H0/S0lab:.1e} (does not bind for any eps < 1e11)")
for name, eps0, xi in (("trace flat", 3.0, xi_T), ("k=+1 axial", eps_k(0), xi_k), ("MZ axial (N_eff-only case)", eS, xi_mz)):
    p0 = xi * eps0 * H0
    print(f"   {name}: p today = xi eps H0 = {p0:.2e} eV, |k_AF| ~ p/2 = {p0/2:.1e} eV vs radio-galaxy CFJ bound {kaf:.0e} eV")
    check(f"N4 {name}: today's CFJ vector is below the radio-galaxy k_AF bound", p0 / 2 < kaf)
b2 = math.degrees(xi_T * 1.5 * math.log(1 + zrg))
print(f"   trace flat: beta(z={zrg}) = {b2:.3f} deg vs radio-galaxy UV rotation {th_rg} +- {th_rg_e} deg")
check("N5 low-z rotation (z=2) is far below the radio-galaxy UV limit", abs(b2) < th_rg_e)

# gauge-invariant higher-dimension alternative: -(1/4)(div S / M^2) F Ft, S_t = eps H
# div S = a^-3 d(a^3 S^t)/dt = -eps (Hdot + 3H^2) = -eps H^2 (3 + Hdot/H^2); matter era Hdot/H^2 = -3/2
for eps in (1.0,):
    theta_star = eps * Hs**2 * 1.5
    M = math.sqrt(theta_star / (2 * beta))
    print(f"   (div S) F Ft / M^2 with eps = {eps}: needs M = {M:.1e} eV (H* = {Hs:.1e} eV): cutoff below H*")
    check("N6 dimension-6 gauge-invariant torsion-photon operators need M ~ H* (not an EFT)", M < 100 * Hs)

# ---------------------------------------------------------------- table
out = os.path.join(HERE, "ec06_backgrounds.csv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["theory", "component", "eps_z_form", "eps_today", "eps_zstar", "needed_xi_and_notes", "source"])
    for r in rows:
        w.writerow([f"{v:.4g}" if isinstance(v, float) else v for v in r])
print(f"wrote ec06_backgrounds.csv ({len(rows)} rows); eps = X_t/H, beta = (xi/2) int eps dz/(1+z)")
summary()
