"""Order-of-magnitude size estimates for the candidate mechanisms (PRELIMINARY).

Formulas (derived and checked in check_cfj_rotation.py, check_dirac_torsion.py; our conventions):
  beta_paper = +(1/2) int p_tau dtau      (beta > 0 clockwise on the sky, arXiv:2509.13654)
  phi F Ftilde baseline, L = -(g/4) phi F Ftilde:   beta = -(g/2) [phi_0 - phi_LSS]
  candidate (i), L = -(1/2) xi S_m A_n Ftilde^{mn}:  p_tau = xi S_tau = -xi a S^0hat
        => beta = -(xi/2) int S^0hat dt     (S^0hat = a S^tau, comoving-observer component)
  minimal EC + Dirac: S^m = -(kappa/2) J5^m  => beta = +(xi kappa/4) int n5 dt, n5 = J5^0hat
Only magnitudes are used below; the overall sign is a free sign of xi or of n5.
"""
import csv
import json
import math
import os
import sys

from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    with open(os.path.join(HERE, "inputs.json")) as f:
        d = json.load(f)
    out = {k: v["value"] for k, v in d.items() if isinstance(v, dict)}
    # beta from the shared data file (row DPK25-beta-baseline, arXiv:2509.13654)
    with open(os.path.join(HERE, "..", "..", "data", "birefringence_measurements.csv"), newline="") as f:
        rows = {r["id"]: r for r in csv.DictReader(f)}
    out["beta_deg"] = float(rows["DPK25-beta-baseline"]["value"])
    return out


def main():
    c = load()
    ok = True
    beta = math.radians(c["beta_deg"])
    hbar = c["hbar_eV_s"]
    H0 = c["H0_km_s_Mpc"] / c["Mpc_km"] * hbar  # eV
    MP = c["M_P_reduced_GeV"] * 1e9
    kappa = 1 / MP**2
    T = c["k_B_eV_K"] * c["T_cmb_K"]
    rho_c = 3 * H0**2 * MP**2
    rho_g = math.pi**2 / 15 * T**4
    Om_r = rho_g * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * c["N_eff"]) / rho_c
    Om_m = c["Omega_m"]
    Om_L = 1 - Om_m - Om_r
    a_star = 1 / (1 + c["z_star"])

    def H(a):
        return H0 * math.sqrt(Om_r / a**4 + Om_m / a**3 + Om_L)

    dt = quad(lambda a: 1 / (a * H(a)), a_star, 1, limit=200)[0]  # eV^-1
    I3 = quad(lambda a: a**-3 / (a * H(a)), a_star, 1, limit=200)[0]
    gyr = dt * hbar / (3.15576e16)
    print(f"beta = {c['beta_deg']} deg = {beta:.4e} rad;  H0 = {H0:.4e} eV;  kappa = {kappa:.4e} eV^-2")
    print(f"t(LSS->today) = {dt:.4e} eV^-1 = {gyr:.3f} Gyr;  int a^-3 dt = {I3:.4e} eV^-1")
    ok &= gyr > 13.0 and gyr < 14.5
    print(("PASS" if ok else "FAIL") + ": cosmic time since LSS in [13, 14.5] Gyr")

    print("\n[baseline phi F Ftilde]  |g Delta_phi| = 2|beta| = %.3e (rad)" % (2 * beta))
    print("  sign (our conventions): beta > 0 needs g (phi_0 - phi_LSS) < 0")

    xiS_const = 2 * beta / dt
    xiS_a3 = 2 * beta / I3
    print("\n[(i) CFJ, xi S]  constant S^0hat:   |xi S^0hat| = %.2e eV  (= %.2e H0)" % (xiS_const, xiS_const / H0))
    print("                 S^0hat ~ a^-3:      |xi S^0hat(today)| = %.2e eV" % xiS_a3)
    print("  gauge-breaking 'torsionful minimal coupling' is the case xi = 2 (plus photon mass terms).")

    # minimal EC: S^0hat = -(kappa/2) n5
    Tnu = (4 / 11) ** (1 / 3) * T
    xn = c["xi_nu_max"]
    dn_flavor = Tnu**3 / (6 * math.pi**2) * (math.pi**2 * xn + xn**3)
    n5_nu = 3 * dn_flavor  # |J5^0| for relativistic Weyl neutrinos = |n_nu - n_nubar|
    cm3 = c["hbar_c_eV_cm"] ** 3  # 1 cm^-3 = cm3 eV^3
    S_nu = kappa / 2 * n5_nu
    print("\n[minimal EC source, S = -(kappa/2) J5]")
    print(f"  CnuB asymmetry (|xi_nu| = {xn}, 3 flavors, relativistic, upper bound): n5 = {n5_nu:.2e} eV^3 = {n5_nu / cm3:.1f} cm^-3")
    print(f"  |S^0hat| today = {S_nu:.2e} eV")
    xi_need_const = xiS_a3 / S_nu  # a^-3 scaling of n5 (generous: ignores non-relativistic suppression)
    print(f"  needed |xi| (n5 ~ a^-3 from LSS) = {xi_need_const:.1e}  -> log10 shortfall at xi = 1: {math.log10(xi_need_const):.1f}")
    rho_dm = c["Omega_c_h2"] * (H0 / (c["H0_km_s_Mpc"] / 100)) ** 2 * 3 * MP**2
    n_dm = rho_dm / c["m_fermion_dm_eV"]
    S_dm = kappa / 2 * n_dm
    print(f"  fermion DM (m = {c['m_fermion_dm_eV']:.0f} eV, fully chiral, generous): n5 = {n_dm / cm3:.2f} cm^-3, |S^0hat| = {S_dm:.1e} eV,"
          f" needed |xi| = {xiS_a3 / S_dm:.1e}")

    # (iv) radiative CFJ from a fermion loop: |p| = 2 C e^2 |b|, |b| = (3/4)|S|
    e2 = 4 * math.pi * c["alpha_em"]
    fac = 2 * c["C_JK"] * e2 * 0.75
    print(f"\n[(iv) fermion-loop CFJ, one charged Dirac fermion, C = {c['C_JK']:.4f}]  |p|/|S| = {fac:.2e}")
    print(f"  needed |S^0hat| (constant) = {xiS_const / fac:.1e} eV;  minimal-EC S from CnuB: {S_nu:.1e} eV")

    # (ii)/(iii) Planck-suppressed pseudoscalar couplings: derivative coupling d phi J5 / f gives a
    # photon coupling only through the fermion loop, suppressed by (omega/m_f)^2 for omega << m_f.
    me = 0.51099895e6
    print(f"\n[(ii),(iii) loop-induced phi F Ftilde from d(phi) J5_e/f]  suppression ~ (H_LSS/m_e)^2 = {(H(a_star) / me) ** 2:.1e}")
    print(f"  with g added by hand, g ~ 1/M_P:  |Delta_phi| = 2|beta| M_P = {2 * beta * MP / 1e9:.2e} GeV")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
