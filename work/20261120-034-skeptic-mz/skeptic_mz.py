"""Task 034 skeptic of 031 (EC-06a, Magueijo-Zlosnik). Independent checks; own loader, no import from 031.

K1  Stueckelberg vertex from scratch: delta(-1/2 p.A Ftilde) under A -> A + d sigma equals
    -1/4 sigma (dp)_{mu nu} Ftilde^{mu nu} + total derivative (explicit 4D components, generic fields).
K2  Vertex coefficient for the non-gradient piece p_i = xi (1 - c_s) d_i phi/(sqrt3 M_P): c = xi (1-c_s)/(2 sqrt3 M_P m);
    power counting gamma gamma -> phi phi by pi exchange: M ~ c^2 E^4 (dimensionless check); Lambda_SC numbers.
K3  theta normalisation, independent route: MZ FRW Eqs. 26-29 (1908.05184) at O(gamma^2) give
    H^2 = (Lambda + kappa rho)/3 + 5P^2 - 2 gamma H P, i.e. -gamma^2 H^2/9 on the tracker P = gamma H/3:
    the torsion sector carries energy O(M_P^2 P^2), so S ~ d phi/M_P to O(1) with gamma cancelling.
K4  Closure in matter without free data (separate-universe, dust, Lambda = 0): from nabla.(G.S) = (gamma/4) GB and
    GB = C^2 - (4/3) kappa^2 rho^2 (dust, Einstein eq.), s = -u.S obeys s_dot = gamma kappa rho/3 along each worldline.
    The growing mode delta ~ t^(2/3) gives delta s = 3 s_bar delta + const(x): d_i s != 0 for every choice of data.
K5  Off-tracker, numerically: integrate the FRW transport law for eps = 2P/H forward (tracker vs 031 closed form) and
    backward from today with eps(0) shifted (response of eps(z*)).
K6  k* condition vs frequency: xi_k scales with f; paper's own Eq. 70 (2001.06373): |gamma| << 1e-21.
K7  Alternative gauge-invariant repair: couple only to the exact-gradient vacuum piece S_vac = (1/gamma) d ln Lambda
    (= c_s S_t on the background). beta_alt, xi needed, vacuum ALP coupling g = xi/(sqrt3 M_P) for each gamma bound.
"""
import csv
import math
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
NP = {"PASS": 0, "FAIL": 0}


def check(name, ok):
    NP["PASS" if ok else "FAIL"] += 1
    print(f"{'PASS' if ok else 'FAIL'}: {name}")


def load(name):
    with open(os.path.join(ROOT, "data", name), newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


# ---------------- K1 ----------------
X = sp.symbols("x0:4", real=True)
sig = sp.Function("sigma")(*X)
p = [sp.Function(f"p{i}")(*X) for i in range(4)]
A = [sp.Function(f"A{i}")(*X) for i in range(4)]
eup = lambda *i: -sp.LeviCivita(*i)  # eps^{0123} = -1 (flat, CONVENTIONS.md)


def Ft(Avec):
    F = [[sp.diff(Avec[n], X[m]) - sp.diff(Avec[m], X[n]) for n in range(4)] for m in range(4)]
    return [[sp.Rational(1, 2) * sum(eup(m, n, r, s) * F[r][s] for r in range(4) for s in range(4))
             for n in range(4)] for m in range(4)]


FtA = Ft(A)
L = lambda Av: -sp.Rational(1, 2) * sum(p[m] * Av[n] * FtA[m][n] for m in range(4) for n in range(4))
dsig = [sp.diff(sig, X[n]) for n in range(4)]
dL = L(dsig)  # Ftilde unchanged under gauge shift
dp = [[sp.diff(p[n], X[m]) - sp.diff(p[m], X[n]) for n in range(4)] for m in range(4)]
target = -sp.Rational(1, 4) * sig * sum(dp[m][n] * FtA[m][n] for m in range(4) for n in range(4))
Xv = [-sp.Rational(1, 2) * sig * sum(p[m] * FtA[m][n] for m in range(4)) for n in range(4)]
td = sum(sp.diff(Xv[n], X[n]) for n in range(4))
check("K1 delta(-1/2 p.A Ftilde) = -1/4 sigma (dp).Ftilde + d_nu(-1/2 sigma p_mu Ftilde^{mu nu})",
      sp.expand(dL - target - td) == 0)

# ---------------- K2 ----------------
cos = load("cosmo_inputs.csv")
c = {k: float(v["value"]) for k, v in cos.items()}
ecc = {k: float(v["value"]) for k, v in load("ec06c_costs_inputs.csv").items() if k.startswith("m_gamma")}
th = {k: float(v["value"]) for k, v in load("ec06_theory_inputs.csv").items() if k.startswith("cT")}
mzi = {k: float(v["value"]) for k, v in load("mz_inputs.csv").items()}
beta = math.radians(float(load("birefringence_measurements.csv")["DPK25-beta-baseline"]["value"]))
MP, hbar = c["MP_reduced_eV"], c["hbar_eV_s"]
H0 = c["H0"] / c["Mpc_km"] * hbar
Om = c["Omega_m"]
h = c["H0"] / 100
Or = 2.469e-5 / h**2 * (c["T_cmb_K"] / 2.7255)**4 * (1 + 0.2271 * c["N_eff_SM"])
OL = 1 - Om - Or
zs = c["z_star"]

# vertex: p_t = 0, p_i = k d_i phi, k = xi (1-c_s)/(sqrt3 M_P); (dp)_{ti} = k d_t d_i phi; -(1/4) sigma (dp).Ft
# = -(1/2) sigma k phi_ti Ft^{ti}; sigma = pi/m  =>  |c| = k/(2m)
xi, cs, m, E = sp.symbols("xi c_s m E", positive=True)
k_ = xi * (1 - cs) / (sp.sqrt(3) * sp.Symbol("M_P", positive=True))
cvert = k_ / (2 * m)
check("K2a c = xi (1-c_s)/(2 sqrt3 M_P m) (same as 031 T5, derived from K1)",
      sp.simplify(cvert - xi * (1 - cs) / (2 * sp.sqrt(3) * sp.Symbol("M_P", positive=True) * m)) == 0)
# dimensions: [c] = -2 (M_P, m mass dim 1); vertex pi phi_ti F: 1+3+2 = 6 -> c dim -2; amp c^2 E^3 E^3/E^2 dim 0
dim_c = -2
check("K2b power counting: operator pi d d phi F has dim 6, [c] = -2, amplitude c^2 E^6/E^2 = c^2 E^4 dimensionless",
      1 + (1 + 2) + 2 == 6 and 2 * dim_c + 4 == 0)


def lam_sc(xi_, m_, omc, conv=16 * math.pi):
    cc = xi_ * omc / (2 * math.sqrt(3) * MP * m_)
    return (conv / cc**2)**0.25


def Efn(a):
    return math.sqrt(Om / a**3 + Or / a**4 + OL)


Lam = 3 * OL


def eps_tr(a):  # eps_S/gamma on the tracker (031 T2a closed form; checked numerically in K5)
    I3 = quad(lambda b: b**2 / Efn(b), 0, a, limit=200)[0]
    return 2 / 3 - 2 / 9 * Lam**2 * I3 / (a**3 * Efn(a)**3)


a_s = 1 / (1 + zs)
J = quad(lambda a: eps_tr(a) / a, a_s, 1, limit=200)[0]
xig = 2 * beta / J
g_gw = math.sqrt(th["cT_GW170817_hi"] / 0.057)
xi_gw = xig / g_gw
f_nu = 0.2271
Neff_lo = c["N_eff_obs"] - 2 * float(cos["N_eff_obs"]["err"])
g_neff = math.sqrt(9 * f_nu * (c["N_eff_SM"] - Neff_lo) / (1 + f_nu * c["N_eff_SM"]))
xi_noGW = xig / g_neff
mP = ecc["m_gamma_PDG"]
omc_mean = Om / (Om + OL)
print(f"K2 xi*gamma = {xig:.3e}; xi(GW speed) = {xi_gw:.2e}; xi(tracker, N_eff only, no GW) = {xi_noGW:.2e}")
for lab, xx in (("GW-speed tracker", xi_gw), ("tracker, no GW bound", xi_noGW)):
    for omc in (omc_mean, 1.0):
        print(f"   {lab}: 1-c_s = {omc:.3f}: Lambda_SC(PDG) = {lam_sc(xx, mP, omc):.2e} eV, "
              f"(Proca gal) {lam_sc(xx, ecc['m_gamma_Proca_gal'], omc):.2e} eV; "
              f"unitarity |a0|<=1/2 with a0=M/(16pi): {lam_sc(xx, mP, omc, 8 * math.pi):.2e} eV")
L1 = lam_sc(xi_gw, mP, omc_mean)
check("K2c Lambda_SC = 3-4 keV at xi(GW speed), PDG m, 1-c_s = Omega_m (031: 3.6 keV)", 3e3 < L1 < 4.5e3)
check("K2d no-GW-bound tracker (gamma at N_eff limit) still < 20 MeV at PDG m", lam_sc(xi_noGW, mP, omc_mean) < 2e7)
E_safe = mzi["E_collider_safe"]
for xx, lab in ((xi_gw, "GW"), (xi_noGW, "noGW")):
    mneed = (E_safe / (16 * math.pi)**0.25)**2 * xx * omc_mean / (2 * math.sqrt(3) * MP)
    print(f"   m needed for 10 TeV ({lab}): {mneed:.1e} eV = 10^{math.log10(mneed / mP):.1f} x PDG; "
          f"energy shortfall {math.log10(E_safe / lam_sc(xx, mP, omc_mean)):.1f} orders")
# c_s momentum enhancement (031 did not use it): phi momenta E/c_s; vertex has one spatial derivative on phi,
# so the amplitude gains >= 1/c_s^2 (two vertices) and Lambda_SC drops by >= sqrt(c_s).
print(f"   with the phi momentum enhancement, Lambda_SC -> Lambda_SC sqrt(c_s) at least; c_s = {OL / (Om + OL):.2f} today "
      f"(cosmic mean), -> 0 wherever kappa rho >> Lambda")

# ---------------- K3 ----------------
g_, P_, T_, Ld, Lm, kr, gm, Hh = sp.symbols("g P T Ldot Lambda krho gamma H")
Ld29 = sp.solve(sp.Eq(P_, 3 * Ld / Lm**2 * (g_ * P_ + (Lm + kr) / (6 * gm))), Ld)[0]
T28 = (Ld / (2 * Lm**2) * (Lm + kr - 6 * g_ * P_ / gm)).subs(Ld, Ld29)
# g = H + T (031 A1: adot/a = g - T); solve Eq. 26 for H to O(gamma^2) with P = gamma * pt, H fixed: compute residual
pt = sp.Symbol("pt")
eps_ = sp.Symbol("e")  # bookkeeping: gamma -> e gamma
expr = ((Hh + T28.subs(g_, Hh + sp.Symbol("Tg")))**2 - P_**2 - (Lm + kr) / 3)
# iterate T at leading order: T = gamma P - 6 g P^2/(Lambda + krho) + ...
T_lead = sp.series(T28.subs({P_: eps_ * gm * pt, gm: eps_ * gm, g_: Hh}, simultaneous=True).subs(kr, 3 * Hh**2 - Lm), eps_, 0, 3).removeO()
fried = sp.expand((Hh + T_lead)**2 - (eps_ * gm * pt)**2 - Hh**2)  # = (Lambda+krho)/3 - H^2 with krho = 3H^2 - Lambda
fried2 = sp.series(fried, eps_, 0, 3).removeO().coeff(eps_, 2)
print(f"K3 O(gamma^2) of g^2 - P^2 - H^2 at P = gamma pt: {sp.factor(fried2)}  [Friedmann: H^2 = (L+kr)/3 - this]")
deficit = sp.simplify(-fried2.subs(pt, Hh / 3))
check("K3 MZ Eqs. 26-29: H^2 - (Lambda + kappa rho)/3 = 5P^2 - 2 gamma H P = -gamma^2 H^2/9 on the tracker",
      sp.simplify(-fried2 - (5 * (gm * pt)**2 - 2 * gm**2 * Hh * pt)) == 0 and sp.simplify(deficit + gm**2 * Hh**2 / 9) == 0)
# energy O(M_P^2 P^2) <-> phi_dot ~ M_P P ~ sqrt3 M_P S: S = O(1) d phi/M_P, gamma-free (031: 1/sqrt3)

# ---------------- K4 ----------------
kap, rho, Lc = sp.symbols("kappa rho Lambda_c", positive=True)
eta = sp.diag(-1, 1, 1, 1)
u = sp.Matrix([1, 0, 0, 0])
ud = eta * u
Tmn = rho * ud * ud.T
Ttr = sum(eta[i, i] * Tmn[i, i] for i in range(4))
Rmn = kap * (Tmn - Ttr * eta / 2) + Lc * eta
Rup = eta * Rmn * eta
RR = sum(Rmn[i, j] * Rup[i, j] for i in range(4) for j in range(4))
Rs = sum(eta[i, i] * Rmn[i, i] for i in range(4))
GB_noW = sp.expand(-2 * RR + sp.Rational(2, 3) * Rs**2)  # GB = C^2 - 2 R_mn R^mn + (2/3) R^2
# flat FRW check: GB = 24 H^2 (Hdot + H^2), H^2 = (krho + L)/3, Hdot + H^2 = (L - krho/2)/3
GB_frw = 24 * (kap * rho + Lc) / 3 * (Lc - kap * rho / 2) / 3
check("K4a GB = C^2 + f(rho) with f = 24 (krho+L)/3 (L - krho/2)/3: GB is local in rho at linear order (C^2 is 2nd order)",
      sp.simplify(GB_noW - GB_frw) == 0)
t = sp.symbols("t", positive=True)
d0, t0, gmm = sp.symbols("delta0 t0 gamma", positive=True)
kr_bar = 4 / (3 * t**2)  # matter era, Lambda = 0: krho = 3H^2, H = 2/(3t)
# transport: d/dt s = -(gamma/(4 krho)) GB = gamma krho/3 (Lambda = 0, C^2 dropped)
sdot = sp.simplify(-(gmm / (4 * kap * rho)) * GB_noW.subs(Lc, 0))
check("K4b s_dot = gamma kappa rho/3 along each dust worldline (Lambda = 0)", sp.simplify(sdot - gmm * kap * rho / 3) == 0)
sbar = sp.integrate(gmm * kr_bar / 3, t)  # = -4 gamma/(9t): |s| = 2 gamma H/3 = 2P (tracker)
check("K4c background s = -(2/3) gamma H, |s| = 2P with P = gamma H/3 (tracker, matches 031 T2)",
      sp.simplify(sbar + sp.Rational(2, 3) * gmm * 2 / (3 * t)) == 0)
delta = d0 * (t / t0)**sp.Rational(2, 3)
ds = sp.integrate(gmm * kr_bar * delta / 3, t)
check("K4d growing mode: delta s = 3 s_bar delta (+ x-dependent constant): grows as t^(-1/3), so no constant of "
      "integration can cancel it; d_i s != 0 => dS != 0 for every solution with a growing density mode",
      sp.simplify(ds - 3 * sbar * delta) == 0)

# ---------------- K5 ----------------
# transport law d/dt[a^3 3H^2 2P] = (gamma/4) a^3 (GB - (8/3) Lambda^2) rewritten for eps = 2P/H (gamma = 1, H0 = 1):
# Y = 3 a^3 H^3 eps, dY/dlna = (3/4) a^3 (GB - (8/3) Lambda^2)/H
def deps(lna, y):
    a = math.exp(lna)
    Ea = Efn(a)
    dlnE = a * (-3 * Om / a**4 - 4 * Or / a**5) / (2 * Ea**2)
    GB = 24 * Ea**2 * (Ea**2 * dlnE + Ea**2)
    dY = 0.25 * a**3 * (GB - 8 / 3 * Lam**2) / Ea
    return [dY / (3 * a**3 * Ea**3) - 3 * y[0] * (1 + dlnE)]


a_i = 1e-4  # start late: numerical errors grow as (aH)^-3 (the homogeneous mode)
fw = solve_ivp(deps, [math.log(a_i), 0], [eps_tr(a_i)], rtol=1e-12, atol=1e-16, method="LSODA", dense_output=True)
e_star_fw, e0_fw = fw.sol(math.log(a_s))[0], fw.sol(0.0)[0]
print(f"K5 forward-integrated tracker eps/gamma: z* {e_star_fw:.4f}, today {e0_fw:.4f}; closed form "
      f"{eps_tr(a_s):.4f}, {eps_tr(1.0):.4f}")
check("K5a numerical integration of the transport law reproduces 031's closed-form tracker to 1e-3",
      abs(e_star_fw - eps_tr(a_s)) < 1e-3 and abs(e0_fw - eps_tr(1.0)) < 1e-3)
bk = [solve_ivp(deps, [0.0, math.log(a_s)], [e0_fw + d], rtol=1e-12, atol=1e-16, method="LSODA").y[0][-1]
      for d in (0.0, 1.0)]
resp = abs(bk[1] - bk[0])
print(f"   backward from today with eps(0) shifted by 1 (gamma units): eps(z*) shifts by {resp:.2e}")
check("K5b any shift of eps(0) moves eps(z*) by < 1e-3 of it: no branch large at z* and small today "
      "except tuned cancellation at z ~ 0", resp < 1e-3)

# ---------------- K6 ----------------
eps0 = eps_tr(1.0)
print("K6 k* = H/eps_S > k_GW today -> xi_min vs frequency:")
xik = {}
for f in (10.0, mzi["GW170817_f_low"], 100.0, 1000.0):
    r = H0 / (2 * math.pi * f * hbar)
    xik[f] = xig * eps0 / r
    print(f"   f = {f:6.0f} Hz: gamma < {r / eps0:.1e}, xi >= {xik[f]:.1e}, Lambda_SC(PDG, Omega_m) = "
          f"{lam_sc(xik[f], mP, omc_mean):.1e} eV")
check("K6a for any f in 10 Hz - 1 kHz the k* condition needs xi >= 1e16 (031: 4.7e16 at 24 Hz)", min(xik.values()) > 1e16)
xi70 = xig / 1e-21
print(f"   2001.06373 Eq. 70 (|gamma| << 1e-21): xi >> {xi70:.1e}")
check("K6b paper's own Eq. 70 bound gives xi > 1e18, stronger than 031's 4.7e16", xi70 > 1e18)

# ---------------- K7 ----------------
cs_bg = lambda a: OL / (OL + Om / a**3 + Or / a**4)
Jalt = quad(lambda a: cs_bg(a) * eps_tr(a) / a, a_s, 1, limit=200)[0]
xig_alt = 2 * beta / Jalt
print(f"K7 S_vac = (1/gamma) d ln Lambda = c_s S on the background; int c_s eps dlna = {Jalt:.3f} gamma vs {J:.3f}: "
      f"xi*gamma(alt) = {xig_alt:.2e}")
for lab, gb in (("GW speed (outside k* regime)", g_gw), ("k* at 24 Hz", 1.45e-20 / eps0 * 1.0),
                ("Eq. 70 gamma = 1e-21", 1e-21), ("N_eff only", g_neff)):
    if lab.startswith("k*"):
        gb = H0 / (2 * math.pi * mzi["GW170817_f_low"] * hbar) / eps0
    xa = xig_alt / gb
    gvac = xa / (math.sqrt(3) * MP) * 1e9
    print(f"   {lab}: gamma = {gb:.1e}, xi = {xa:.1e}, vacuum ALP coupling g = xi/(sqrt3 M_P) = {gvac:.1e} GeV^-1; "
          f"in matter g_eff = c_s g -> 0")
check("K7 the vacuum-gradient-only coupling needs xi larger by J/J_alt > 3 than 031's full-S coupling", J / Jalt > 3)
print(f"== {NP['PASS']} PASS, {NP['FAIL']} FAIL ==")
