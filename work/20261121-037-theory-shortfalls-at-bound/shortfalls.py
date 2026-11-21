"""037: maximal beta for the spin-sourced routes at the largest allowed xi (or xi*c), per completion class.

Operator (CONVENTIONS, 007): -(1/2) xi S.A Ft with S = -(kappa/2) J  ->  C J.A Ft,  C = xi kappa/4  (dimension 6).
Reig 2611.07578 normalisation: -(1/2)(1/Lambda*^2) j.A Ft  ->  1/Lambda*^2 = 2C = xi kappa/2.
beta(xi) = xi * beta(xi = 1)  (linear in C).

Routes (beta at xi = 1 recomputed here, cross-checked against 007 / 029):
  EC-02   relic-nu asymmetry, cited |xi_nu| (95%) and the generous N_eff-only case        J = universal J5
  Pauli   any fermionic axial source, whole LCDM density as a maximally chiral gas (029)   J = universal J5
  vector  vector baryon current via a non-minimal coupling, S^0 = -(c kappa/2) n_B (029)   J = J_B, coefficient xi*c

Completion classes for the non-gauge-invariant operator (check_gauge.py: delta O ~ curl J):
  D  massless photon: Kalb-Ramond / dual massless axion (Reig: required for a local Lorentz-invariant completion)
       axial: g_ag g_ae = 4 m_e C = xi kappa m_e;     vector (conserved J): contact 1/(2 f^2) J^2, f = g_ag Lambda*^2
  S  Stueckelberg (massive photon): pi couples as (C/m_gamma) pi dJ Ft, Lambda_s = (m_gamma/C)^(1/3)
       S-a  Lambda_s >= T_BBN (Reig's criterion, conservative)
       S-b  NDA: pi must decouple above T = m_mu (else Delta N_eff >= 0.39); Gamma ~ T^7/Lambda_s^6 < H
All inputs from data/.
"""
import csv
import math
import os

from scipy.integrate import quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))


def rows(name):
    with open(os.path.join(ROOT, "data", name), newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def num(s):
    return float(str(s).replace("<", "").strip())


def report(name, ok):
    print(("PASS" if ok else "FAIL") + ": " + name)
    return bool(ok)


C0 = rows("cosmo_inputs.csv")
c = {k: float(r["value"]) for k, r in C0.items()}
beta_obs = float(rows("birefringence_measurements.csv")["DPK25-beta-baseline"]["value"])
B = rows("bounds.csv")
g_CAST = num(B["A-CAST"]["value_eV"])  # eV^-1
g_Chandra = num(B["A-Chandra"]["value_eV"])
mg_PDG = float(rows("ec06c_costs_inputs.csv")["m_gamma_PDG"]["value"])
mg_gal = float(rows("ec06c_costs_inputs.csv")["m_gamma_Proca_gal"]["value"])
mg_Ret = float(rows("compact_star_inputs.csv")["m_gamma_Retino"]["value"])
T_BBN = float(rows("bounds_inputs.csv")["T_BBN"]["value"])
N = {k: float(r["value"]) for k, r in rows("shortfall_bound_inputs.csv").items()}
GeV = 1e9

# ---------------- cosmology (flat LCDM, as 007/029)
hb, MP = c["hbar_eV_s"], c["MP_reduced_eV"]
kap = 1 / MP**2
H0 = c["H0"] / c["Mpc_km"] * hb
T0 = c["T_cmb_K"] * c["k_B_eV_K"]
Tnu0 = (4 / 11) ** (1 / 3) * T0
rho_c = 3 * H0**2 * MP**2
Or = math.pi**2 / 15 * T0**4 * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * c["N_eff_SM"]) / rho_c
Om = c["Omega_m"]
OL = 1 - Om - Or
a_s = 1 / (1 + c["z_star"])
H = lambda a: H0 * math.sqrt(Or / a**4 + Om / a**3 + OL)
tint = lambda f: quad(lambda la: f(math.exp(la)) / H(math.exp(la)), math.log(a_s), 0, limit=400, epsabs=0, epsrel=1e-10)[0]
deg = math.degrees
ok = True

# ---------------- beta at xi = 1 for each route: beta = (xi kappa/4) int n dt  (007 Sec. A)
dn = lambda T, x: T**3 / (6 * math.pi**2) * (math.pi**2 * x + x**3)  # massless Weyl species, n - nbar = |n5|
xi95 = c["xi_nu"] + 2 * float(C0["xi_nu"]["err"])
dN95 = c["N_eff_obs"] + 2 * float(C0["N_eff_obs"]["err"]) - c["N_eff_SM"]
dneff = lambda x: 30 / 7 * (x / math.pi) ** 2 + 15 / 7 * (x / math.pi) ** 4
xi_gen = brentq(lambda x: dneff(x) - dN95 / 2, 0, 5)
b_nu = deg(kap / 4 * tint(lambda a: 3 * dn(Tnu0 / a, xi95)))
b_nu_gen = deg(kap / 4 * tint(lambda a: dn(Tnu0 / a, xi95) + 2 * dn(Tnu0 / a, xi_gen)))


def n5_pauli(rho, g):  # 029: massless degenerate maximally chiral gas
    return g**0.25 * (8 * math.pi**2 * rho) ** 0.75 / (6 * math.pi**2)


b_pauli = deg(kap / 4 * tint(lambda a: n5_pauli(3 * H(a) ** 2 * MP**2, 100)))
nB0 = c["Omega_b_h2"] * 3 * (100 / c["Mpc_km"] * hb) ** 2 * MP**2 / c["m_p_eV"]
b_vec = deg(kap / 4 * tint(lambda a: nB0 / a**3))  # xi = c = 1

print(f"beta(xi=1): EC-02 cited xi_nu<={xi95:.3f}: {b_nu:.2e} deg; generous (xi_mu,tau={xi_gen:.2f}): {b_nu_gen:.2e}; "
      f"Pauli: {b_pauli:.2e}; vector (xi=c=1): {b_vec:.2e}")
ok &= report("EC-02 cited case reproduces 007 (7.5e-30 deg) within 3%", abs(b_nu / 7.5e-30 - 1) < 0.03)
ok &= report("EC-02 generous case reproduces 007 (6.3e-29 deg) within 5%", abs(b_nu_gen / 6.3e-29 - 1) < 0.05)
ok &= report("Pauli envelope reproduces 029 (6.91e-27 deg) within 2%", abs(b_pauli / 6.91e-27 - 1) < 0.02)
ok &= report("vector route reproduces 029 (xi c needed 1.7e36) within 5%", abs(beta_obs / b_vec / 1.7e36 - 1) < 0.05)

# ---------------- class D: dual massless axion
me = c["m_e_eV"]
xi_D_CAST = c["helioscope_gag_gae_GeV"] / GeV / (kap * me)
xi_D_ChRG = g_Chandra * N["RG_gae_TRGB"] / (kap * me)
xi_D_CaRG = g_CAST * N["RG_gae_TRGB"] / (kap * me)
ok &= report(f"CAST g_ag*g_ae bound -> xi < {xi_D_CAST:.2e} (007: 9.4e17)", abs(xi_D_CAST / 9.4e17 - 1) < 0.02)
# Reig check targets for the mapping
LamR = lambda xinu: c["Reig_Lambda_star_TeV_per_sqrt_xi"] * 1e12 * math.sqrt(xinu)  # eV
prod = me / LamR(0.05) ** 2 * GeV
ok &= report(f"Reig Eq. 11: m_e/Lambda*^2 = {prod:.2e} GeV^-1 vs {N['Reig_gag_gae_required_xi005']:.1e} (10%)",
             abs(prod / N["Reig_gag_gae_required_xi005"] - 1) < 0.1)
fR = N["ALPSII_gag"] / GeV * LamR(0.1) ** 2 / GeV
ok &= report(f"Reig Eq. 12: f = g_ALPS Lambda*^2 = {fR:.2f} GeV vs {N['Reig_f_max_xi01']} GeV (10%)", abs(fR / N["Reig_f_max_xi01"] - 1) < 0.1)
ok &= report(f"Reig Eq. S14: 1/(2f^2) = {1 / (2 * fR**2):.1f} GeV^-2 vs {N['Reig_contact_min_xi01']} (10%)",
             abs(1 / (2 * fR**2) / N["Reig_contact_min_xi01"] - 1) < 0.1)

# vector route in class D: contact (1/2f^2) J_B^2 with f = g Lambda*^2, 1/Lambda*^2 = xi c kappa/2
# CMS LL quark contact: |coef of (qbar_L gamma q_L)^2| < 2 pi/Lambda^2; J_B = (1/3) sum_q qbar gamma q -> C_B/9 < 2 pi/Lambda^2
LamCI = N["CMS_dijet_CI_LL_destructive"] * 1e12
CBmax = 18 * math.pi / LamCI**2  # eV^-2
xic_D = {nm: 2 * g * math.sqrt(2 * CBmax) / kap for nm, g in
         [("ALPS II", N["ALPSII_gag"] / GeV), ("CAST", g_CAST), ("Chandra", g_Chandra)]}

# ---------------- class S: Stueckelberg
xi_S = lambda mg, E: 4 * mg * MP**2 / E**3
mmu = N["m_mu"] * GeV
gstar = 10.75 + 7 / 8 * 4  # muons in (SM count below the QCD crossover, pions neglected)
Lam_s_NDA = (mmu**5 * MP / math.sqrt(math.pi**2 * gstar / 90)) ** (1 / 6)
dNeff_late = 4 / 7 * (10.75 / gstar) ** (4 / 3)
ok &= report(f"a thermal pi decoupling below m_mu gives Delta N_eff >= {dNeff_late:.2f} > 95% room {dN95:.2f}", dNeff_late > dN95)
print(f"S-b NDA: pi decouples above m_mu only if Lambda_s > {Lam_s_NDA:.2e} eV (rate T^7/Lambda_s^6 vs H, g* = {gstar})")

bounds = [
    ("D: CAST g_ag*g_ae (Reig; 007)", xi_D_CAST, "2611.07578 + CAST via Reig ref [45]"),
    ("D: CAST g_ag x RG g_ae", xi_D_CaRG, "2406.16840 x 2007.03694"),
    ("D: Chandra g_ag x RG g_ae", xi_D_ChRG, "1907.05475 x 2007.03694"),
    ("S-a: Lambda_s>T_BBN, m_g Retino weak end", xi_S(mg_Ret, T_BBN), "1302.6168; T_BBN choice"),
    ("S-a: Lambda_s>T_BBN, m_g PDG", xi_S(mg_PDG, T_BBN), "PDG 2024 (Ryutov 07); T_BBN choice"),
    ("S-a: Lambda_s>T_BBN, m_g Proca galactic", xi_S(mg_gal, T_BBN), "hep-ph/0306245 via PDG"),
    ("S-b: NDA decoupling, m_g Retino weak end", xi_S(mg_Ret, Lam_s_NDA), "1302.6168; NDA"),
    ("S-b: NDA decoupling, m_g PDG", xi_S(mg_PDG, Lam_s_NDA), "PDG 2024; NDA"),
]
routes = [("EC-02 nu, cited xi_nu", b_nu), ("EC-02 nu, generous N_eff", b_nu_gen), ("Pauli, any fermion", b_pauli)]
print("\n| bound on xi (J5 operator) | xi_max | source | " + " | ".join(f"{r} beta_max deg / short" for r, _ in routes) + " |")
for nm, xm, src in bounds:
    cells = [f"{b1 * xm:.1e} / {math.log10(beta_obs / (b1 * xm)):.1f}" for _, b1 in routes]
    print(f"| {nm} | {xm:.1e} | {src} | " + " | ".join(cells) + " |")

print("\n| bound on xi*c (vector J_B operator) | (xi c)_max | beta_max deg | shortfall |")
vec = [(f"D: contact 1/(2f^2) J_B^2 < CMS dijet, g = {k}", v) for k, v in xic_D.items()]
vec += [("S-a: T_BBN, m_g Retino weak end", xi_S(mg_Ret, T_BBN)), ("S-a: T_BBN, m_g PDG", xi_S(mg_PDG, T_BBN))]
for nm, v in vec:
    print(f"| {nm} | {v:.1e} | {b_vec * v:.1e} | {math.log10(beta_obs / (b_vec * v)):.1f} |")

# FMT torsion contact alone (no completion): (3 kappa/16) c_q^2 (gamma^2+1) (sum_q qbar gamma q)^2, c_B = 3 c_q
cq_max = math.sqrt(2 * math.pi / LamCI**2 / (3 * kap / 16))
cB_max = 3 * cq_max
print(f"\nFMT contact only: c (baryon normalisation) < {cB_max:.1e} / sqrt(gamma^2+1)")
for nm, xm in [("CAST product", xi_D_CAST), ("S-a PDG", xi_S(mg_PDG, T_BBN))]:
    r_need = beta_obs / b_vec / xm  # needed r = alpha/gamma (baryon normalisation) at zeta = xi_max
    gmax = math.sqrt(cB_max / r_need)  # gamma^2 r < c_max (gamma << 1)
    print(f"  at zeta = {xm:.1e} ({nm}): needed r = {r_need:.1e}; FMT contact allows it iff Immirzi gamma < {gmax:.2f}")

# ---------------- verdict checks
largest = max(xm for _, xm, _ in bounds[:5])  # weakest examined bound (excluding Proca-galactic and NDA rows)
ok &= report(f"no J5 route reaches beta at the weakest examined xi bound ({largest:.1e})",
             all(b1 * largest < beta_obs for _, b1 in routes))
ok &= report("no vector route reaches beta in either completion class", all(b_vec * v < beta_obs for _, v in vec))
ok &= report("ordering: dual-axion bounds (all three) are stronger than the T_BBN Stueckelberg bound with the PDG mass",
             max(xi_D_CAST, xi_D_ChRG, xi_D_CaRG) < xi_S(mg_PDG, T_BBN))
print("ALL PASS" if ok else "SOME FAIL")
