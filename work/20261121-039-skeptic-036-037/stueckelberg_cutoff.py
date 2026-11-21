"""Skeptic 039, item B: Stueckelberg strong-coupling scale of the 037 weak corner, computed independently.

Operator (037, CONVENTIONS): L = C J_mu A_nu Ft^{mu nu}, C = xi kappa/4, J = sum_f fbar gamma gamma5 f (EC: universal).
Stueckelberg photon A -> A + d pi/m (m = photon mass). Since Ft is gauge invariant, only the explicit A moves:
   L_pi = (C/m) J_mu d_nu pi Ft^{mu nu} = -(C/m) pi (d_nu J_mu) Ft^{mu nu}  (Bianchi d_nu Ft^{mu nu} = 0),
equivalently -(1/4m) pi (dp).Ft with p = xi S = -(xi kappa/2) J (the form used in 028/031).

S1  Gauge variation of the operator (sympy): delta O = total derivative - a (d_nu J_mu) Ft^{mu nu} (037's check_gauge).
S2  Tree amplitude f fbar -> gamma pi from explicit massless helicity spinors (numpy, chiral rep, +--- internally;
    |M| is convention independent). Only J = 1 contributes; partial wave a_1 from
       |a_1|^2 = (3/2) int |M_{lam mu}|^2 dcos(theta) / (48 pi)^2   (M = 16 pi sum (2J+1) a_J d^J).
    Lambda_tree = sqrt(s) where max |a_1| = 1. Compare with 037's (m/C)^(1/3).
S3  028/031-style exchange estimate: gamma gamma -> pi pi through the <J5 J5> correlator of N light Dirac fermions
    (one loop, Pi ~ N/(12 pi^2) ln), |a_0| ~ (C/m)^2 E^6 N L/(12 pi^2 * 16 pi). Order of magnitude only.
S4  Numbers at the xi each route needs for beta = 0.215 deg (Pauli envelope recomputed here; EC-02 nu from 037's
    7.5e-30 deg reproduced in 007/037), with the weakest m_gamma (Retino 7.9e-14 eV) and the PDG value.
    Compare with energies where e+e- -> gamma + X and Compton scattering are measured standard.
S5  NDA star-emission cross-check: gamma e -> pi e in a red-giant core vs the axion-Compton rate allowed by
    g_ae < 1.6e-13 (2007.03694). Order of magnitude only.
"""
import csv
import itertools
import math
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
npass = nfail = 0


def check(name, ok):
    global npass, nfail
    npass += bool(ok)
    nfail += (not ok)
    print(("PASS" if ok else "FAIL") + ": " + name)


def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


D = lambda n: load(os.path.join(ROOT, "data", n))
v = lambda d, k: float(d[k]["value"])
cos, meas, ec6, cst, bi, sb = (D("cosmo_inputs.csv"), D("birefringence_measurements.csv"), D("ec06c_costs_inputs.csv"),
                              D("compact_star_inputs.csv"), D("bounds_inputs.csv"), D("shortfall_bound_inputs.csv"))
loc = load(os.path.join(HERE, "inputs_039.csv"))

# ---------------- S1 gauge variation
X = sp.symbols("x0:4", real=True)
eta = sp.diag(1, -1, -1, -1)
al = sp.Function("al")(*X)
Jl = [sp.Function(f"J{i}")(*X) for i in range(4)]
A = [sp.Function(f"A{i}")(*X) for i in range(4)]
Fl = [[sp.diff(A[n], X[m]) - sp.diff(A[m], X[n]) for n in range(4)] for m in range(4)]
Fu = [[sum(eta[m, a] * eta[n, b] * Fl[a][b] for a in range(4) for b in range(4)) for n in range(4)] for m in range(4)]
Ft = [[sum(sp.LeviCivita(m, n, r, s) * Fl[r][s] for r in range(4) for s in range(4)) / 2 for n in range(4)] for m in range(4)]
# with lower F the dual has upper indices (sign/normalisation irrelevant for the identity)
dO = sum(Jl[m] * sp.diff(al, X[n]) * Ft[m][n] for m in range(4) for n in range(4))
rhs_td = sum(sp.diff(al * sum(Jl[m] * Ft[m][n] for m in range(4)), X[n]) for n in range(4))
rhs_rest = -al * sum(sp.diff(Jl[m], X[n]) * Ft[m][n] for m in range(4) for n in range(4))
check("S1 delta(J.A Ft) = total derivative - alpha (d_nu J_mu) Ft^{mu nu}", sp.simplify(sp.expand(dO - rhs_td - rhs_rest)) == 0)

# ---------------- S2 tree amplitude f fbar -> gamma pi
s0, s1, s2, s3 = (np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
sig = [s0, s1, s2, s3]
Z2 = np.zeros((2, 2))
gam = [np.block([[Z2, s0], [s0, Z2]])] + [np.block([[Z2, s], [-s, Z2]]) for s in (s1, s2, s3)]
g5 = np.block([[-s0, Z2], [Z2, s0]])
etan = np.diag([1., -1, -1, -1])
check("S2a Clifford algebra {g^mu, g^nu} = 2 eta^{mu nu} and g5 anticommutes",
      all(np.allclose(gam[a] @ gam[b] + gam[b] @ gam[a], 2 * etan[a, b] * np.eye(4)) for a in range(4) for b in range(4))
      and all(np.allclose(g5 @ gm + gm @ g5, 0) for gm in gam))


def msqrt(Mh):
    w, U = np.linalg.eigh(Mh)
    return U @ np.diag(np.sqrt(np.clip(w, 0, None))) @ U.conj().T


def psig(p):  # p.sigma and p.sigmabar for upper-index p
    ps = p[0] * s0 - p[1] * s1 - p[2] * s2 - p[3] * s3
    pb = p[0] * s0 + p[1] * s1 + p[2] * s2 + p[3] * s3
    return ps, pb


def u_sp(p, xi_):
    ps, pb = psig(p)
    return np.concatenate([msqrt(ps) @ xi_, msqrt(pb) @ xi_])


def v_sp(p, et):
    ps, pb = psig(p)
    return np.concatenate([msqrt(ps) @ et, -msqrt(pb) @ et])


LC = np.zeros((4, 4, 4, 4))
for perm in itertools.permutations(range(4)):
    LC[perm] = sp.LeviCivita(*perm)
lower = lambda vec: etan @ vec


def amp(E, th, xi_, et, hel):
    """M/(C/m) for f(p1) fbar(p2) -> gamma(k, hel) pi(q): M = P_nu eps^{mu nu rho sigma} k_rho eps*_sigma j_mu."""
    p1, p2 = np.array([E, 0, 0, E]), np.array([E, 0, 0, -E])
    k = np.array([E, E * np.sin(th), 0, E * np.cos(th)])
    eps = np.array([0, np.cos(th), 1j * hel, -np.sin(th)]) / np.sqrt(2)
    u, vv = u_sp(p1, xi_), v_sp(p2, et)
    vbar = vv.conj() @ gam[0]
    j_up = np.array([vbar @ gam[mu] @ g5 @ u for mu in range(4)])
    P_l, k_l, e_l, j_l = lower(p1 + p2), lower(k), lower(eps.conj()), lower(j_up)
    return np.einsum("mnrs,m,n,r,s->", LC, j_l, P_l, k_l, e_l)


E = 1.0
spins = [np.array([1, 0], complex), np.array([0, 1], complex)]
a1max = 0.0
for (i, xi_), (j, et), hel in itertools.product(enumerate(spins), enumerate(spins), (1, -1)):
    I = quad(lambda c: abs(amp(E, math.acos(c), xi_, et, hel)) ** 2, -1, 1)[0]
    a1 = math.sqrt(1.5 * I) / (48 * math.pi)
    if I > 1e-12:
        print(f"   spins ({i},{j}) photon hel {hel:+d}: int|M|^2 dcos = {I:.4f} (E=1, C/m=1) -> |a_1| = {a1:.3e} (C/m) E^3")
    a1max = max(a1max, a1)
# |a_1| = K (C/m) E^3 with E the beam energy, sqrt(s) = 2E
K = a1max
# also confirm scaling with E^3
Isc = quad(lambda c: abs(amp(2.0, math.acos(c), spins[0], spins[0], 1)) ** 2 + abs(amp(2.0, math.acos(c), spins[1], spins[1], 1)) ** 2 +
           abs(amp(2.0, math.acos(c), spins[0], spins[1], 1)) ** 2 + abs(amp(2.0, math.acos(c), spins[1], spins[0], 1)) ** 2, -1, 1)[0]
I1 = quad(lambda c: abs(amp(1.0, math.acos(c), spins[0], spins[0], 1)) ** 2 + abs(amp(1.0, math.acos(c), spins[1], spins[1], 1)) ** 2 +
          abs(amp(1.0, math.acos(c), spins[0], spins[1], 1)) ** 2 + abs(amp(1.0, math.acos(c), spins[1], spins[0], 1)) ** 2, -1, 1)[0]
check(f"S2b |M|^2 scales as E^6 (ratio {Isc / I1:.2f} for E x2, expect 64)", abs(Isc / I1 / 64 - 1) < 1e-6)
print(f"S2 max |a_1| = {K:.4e} (C/m) E^3 = {K / 8:.4e} (C/m) s^(3/2)")


def Lam_tree(Cm, amax=1.0):
    """sqrt(s) at which |a_1| = amax (single fermion species)."""
    return 2 * (amax / (K * Cm)) ** (1 / 3)


# ---------------- S3 loop exchange estimate
def Lam_loop(Cm, N, Lg=3.0):
    return (16 * math.pi * 12 * math.pi**2 / (N * Lg)) ** (1 / 6) * Cm ** (-1 / 3)


# ---------------- S4 numbers
MP = v(cos, "MP_reduced_eV")
kap = 1 / MP**2
hb = v(cos, "hbar_eV_s")
H0 = v(cos, "H0") / v(cos, "Mpc_km") * hb
T0 = v(cos, "T_cmb_K") * v(cos, "k_B_eV_K")
rho_c = 3 * H0**2 * MP**2
Or = math.pi**2 / 15 * T0**4 * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * v(cos, "N_eff_SM")) / rho_c
Om = v(cos, "Omega_m")
OL = 1 - Om - Or
Hf = lambda a: H0 * math.sqrt(Or / a**4 + Om / a**3 + OL)
a_s = 1 / (1 + v(cos, "z_star"))
# Pauli envelope (029): whole density as a maximally chiral degenerate gas, g = 100
n5 = lambda rho: 100**0.25 * (8 * math.pi**2 * rho) ** 0.75 / (6 * math.pi**2)
I_n = quad(lambda la: n5(3 * Hf(math.exp(la)) ** 2 * MP**2) / Hf(math.exp(la)), math.log(a_s), 0, limit=400, epsrel=1e-10)[0]
b_pauli = math.degrees(kap / 4 * I_n)
beta_obs = v(meas, "DPK25-beta-baseline")
xi_P = beta_obs / b_pauli
xi_nu = beta_obs / 7.5e-30  # 007's EC-02 cited-nu value (reproduced by 037)
print(f"S4 Pauli beta(xi=1) = {b_pauli:.3e} deg (029/037: 6.91e-27) -> xi needed {xi_P:.2e}; EC-02 nu: xi needed {xi_nu:.2e}")
check("S4a Pauli envelope reproduces 029/037 within 2%", abs(b_pauli / 6.91e-27 - 1) < 0.02)

mg_R = v(cst, "m_gamma_Retino")
mg_P = v(ec6, "m_gamma_PDG")
T_BBN = v(bi, "T_BBN")
E_babar = v(loc, "E_babar_monophoton") * 1e9
E_lep = v(loc, "E_lep2") * 1e9
print("\n| route | m_gamma [eV] | xi | C/m [eV^-3] | 037 (m/C)^(1/3) | tree sqrt(s), |a_1|=1 | tree, |a_1|=1/2 | loop gg->pipi (N=4) |")
res = {}
for rname, xi in (("Pauli envelope", xi_P), ("EC-02 cited nu", xi_nu)):
    for mname, mg in (("weakest (Retino)", mg_R), ("PDG", mg_P)):
        Cm = xi * kap / 4 / mg
        L037, Lt, Lt2, Ll = Cm ** (-1 / 3), Lam_tree(Cm), Lam_tree(Cm, 0.5), Lam_loop(Cm, 4)
        res[(rname, mname)] = Lt
        print(f"| {rname} | {mg:.1e} | {xi:.2e} | {Cm:.2e} | {L037 / 1e6:.3g} MeV | {Lt / 1e6:.3g} MeV | {Lt2 / 1e6:.3g} MeV | {Ll / 1e6:.3g} MeV |")

Lw = res[("Pauli envelope", "weakest (Retino)")]
print(f"\nS4 weak corner (Pauli, weakest m_gamma): Lambda = {Lw / 1e6:.2f} MeV; below BaBar mono-photon sqrt(s) by "
      f"{math.log10(E_babar / Lw):.1f} orders, below LEP2 by {math.log10(E_lep / Lw):.1f} orders")
check("S4b weak-corner cutoff (tree, |a_1| = 1) is below 10 MeV", Lw < 1e7)
check("S4c weak-corner cutoff is > 3 orders below e+e- -> gamma + invisible energies (BaBar)", E_babar / Lw > 1e3)


# xi allowed if Lambda_tree >= E_test, and the resulting Pauli shortfall
def xi_max(mg, Etest):
    Cm = (1.0 / K) / (Etest / 2) ** 3
    return Cm * mg * 4 / kap


for nm, Et in (("1 MeV (037 S-a)", T_BBN), ("BaBar 10.58 GeV", E_babar), ("LEP2 209 GeV", E_lep)):
    xm = xi_max(mg_R, Et)
    print(f"S4 Lambda_tree >= {nm}, weakest m_gamma: xi < {xm:.2e}; Pauli short by {math.log10(xi_P / xm):.1f} orders, "
          f"EC-02 nu short by {math.log10(xi_nu / xm):.1f}")
check("S4d with Lambda >= BaBar energy the Pauli corner is > 10 orders short", math.log10(xi_P / xi_max(mg_R, E_babar)) > 10)

# ---------------- S5 NDA star emission gamma e -> pi e vs axion Compton (red giant core, omega ~ 3T, T ~ 8.6 keV)
me = v(cos, "m_e_eV")
alpha_em = 1 / 137.035999
g_ae = v(sb, "RG_gae_TRGB")
w = 3 * 8.6e3
Cm = xi_P * kap / 4 / mg_R
sig_pi = Cm**2 * w**4 / (4 * math.pi)  # NDA: M ~ (C/m) w q 2 m_e, q ~ w, dsigma/dOmega = |M|^2/(64 pi^2 m_e^2)
sig_a = alpha_em * g_ae**2 * w**2 / (3 * me**4)  # axion Compton, nonrelativistic (Raffelt's textbook form, O(1))
print(f"S5 RG core, w = {w / 1e3:.0f} keV: sigma(g e -> pi e) ~ {sig_pi:.1e} eV^-2 vs allowed axion-Compton {sig_a:.1e} eV^-2 "
      f"(ratio {sig_pi / sig_a:.0e}); xi must drop by ~{0.5 * math.log10(sig_pi / sig_a):.1f} orders")
check("S5 NDA star emission excludes the weak corner by > 10 orders in xi", 0.5 * math.log10(sig_pi / sig_a) > 10)
print(f"== {npass} PASS, {nfail} FAIL ==")
