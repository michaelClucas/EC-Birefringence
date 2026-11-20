"""Skeptic 029: model-independent upper bound on beta from ANY fermionic axial-charge source in minimal EC.

Draft claim tested: "spin-sourced torsion (any coupling) falls >= 28 orders short". 007 computed this only for
cited sources (relic-nu asymmetry, helicity-polarised DM, baryons). A homogeneous fermion condensate or an exotic
dark fermion with a large axial charge (e.g. Poplawski 1102.5667) was not bounded. Here: Pauli bound.

For g fermionic helicity states with total energy density rho_f, each mode contributes |v| <= 1 to |n5| and costs
E >= |p| in energy, occupation <= 1. So |n5| is maximised by a massless degenerate gas:
    n = g p_F^3/(6 pi^2), rho = g p_F^4/(8 pi^2)  =>  |n5| <= g^(1/4) (8 pi^2 rho)^(3/4) / (6 pi^2).
Minimal EC: S^0 = -(kappa/2) n5 (CONVENTIONS / 003). EC-01 coupling: beta = (xi/2) Int S^0 dt (cosmic time).
Most generous choice: rho_f = total rho(z) of LCDM (the whole universe is one maximally chiral degenerate gas),
g = 100. Inputs from data/cosmo_inputs.csv only.
"""
import csv
import math
from pathlib import Path

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
D = {r["id"]: float(r["value"]) for r in csv.DictReader(open(ROOT / "data" / "cosmo_inputs.csv", encoding="utf-8"))}
H0_kms = D["H0"]
Om = D["Omega_m"]
zs = D["z_star"]
MP = D["MP_reduced_eV"]
hbar = D["hbar_eV_s"]
Mpc_km = D["Mpc_km"]
H0 = H0_kms / Mpc_km * hbar  # km/s/Mpc -> 1/s -> eV
Tcmb = D["T_cmb_K"] * D["k_B_eV_K"]
Neff = D["N_eff_SM"]
Or = (math.pi**2 / 15) * Tcmb**4 * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * Neff) / (3 * H0**2 * MP**2)
OL = 1 - Om - Or


def H(z):
    return H0 * math.sqrt(Om * (1 + z) ** 3 + Or * (1 + z) ** 4 + OL)


def n5_max(rho, g):
    return g**0.25 * (8 * math.pi**2 * rho) ** 0.75 / (6 * math.pi**2)


def beta_max_deg(g, frac=1.0):
    kappa = 1 / MP**2
    # dt = dz / ((1+z) H)
    f = lambda z: 0.5 * kappa * n5_max(frac * 3 * H(z) ** 2 * MP**2, g) / ((1 + z) * H(z))
    I = quad(f, 0, zs, limit=400, epsabs=0, epsrel=1e-8)[0]
    return math.degrees(0.5 * I)  # xi = 1


beta_obs = next(float(r["value"]) for r in csv.DictReader(open(ROOT / "data" / "birefringence_measurements.csv",
                encoding="utf-8")) if r["id"] == "DPK25-beta-baseline")
print(f"H0 = {H0:.4e} eV, Omega_r = {Or:.3e}, z* = {zs}")
for g in (2, 100):
    for frac in (1.0, 0.05):
        b = beta_max_deg(g, frac)
        print(f"g = {g:3d}, rho_f = {frac:4.2f} rho_tot: beta_max(xi=1) = {b:.2e} deg; "
              f"xi needed = {beta_obs / b:.1e}; orders short = {math.log10(beta_obs / b):.1f}")
# S/H at z* for the generous case
z = zs
S = 0.5 * n5_max(3 * H(z) ** 2 * MP**2, 100) / MP**2
print(f"S/H at z* (g=100, rho_f = rho_tot) = {S / H(z):.2e}; today: "
      f"{0.5 * n5_max(3 * H0**2 * MP**2, 100) / MP**2 / H0:.2e}")
b = beta_max_deg(100, 1.0)
assert math.log10(beta_obs / b) > 20, "Pauli bound does not close the condensate loophole"
print("PASS: any fermionic axial source in minimal EC is > 20 orders short at xi = 1 (Pauli bound)")

# Part 2: non-minimal fermion-torsion coupling (Holst + Freidel-Minic-Takeuchi hep-th/0507253 type) lets the axial
# torsion be sourced by the VECTOR (baryon) current, S^0 = -(c kappa/2) n_B, which is homogeneous and nonzero.
# Coefficient c needed at xi = 1 with the actual baryon density (Omega_b h^2 from data/cosmo_inputs.csv).
rho_c_h2 = 3 * (100 / Mpc_km * hbar) ** 2 * MP**2  # eV^4 for h = 1
nB0 = D["Omega_b_h2"] * rho_c_h2 / D["m_p_eV"]  # eV^3 (baryons as protons; He correction < 10%)
fB = lambda z: 0.5 * nB0 * (1 + z) ** 3 / MP**2 / ((1 + z) * H(z))
IB = quad(fB, 0, zs, limit=400, epsrel=1e-8)[0]
bB = math.degrees(0.5 * IB)
c_need = beta_obs / bB
print(f"n_B today = {nB0:.2e} eV^3 = {nB0 / D['hbar_c_eV_m']**3:.2f} m^-3; beta(xi=1, c=1) = {bB:.1e} deg; "
      f"c needed = {c_need:.1e} ({math.log10(c_need):.1f} orders)")
print(f"induced B-B contact scale if the four-fermion term is c*kappa*J^2: M_P/sqrt(c) = "
      f"{MP / math.sqrt(c_need) / 1e9:.1e} GeV; if c^2*kappa*J^2: M_P/c = {MP / c_need / 1e9:.1e} GeV")
assert math.log10(c_need) > 20
print("PASS: vector-current-sourced torsion (non-minimal coupling) needs c > 1e20 at xi = 1")
