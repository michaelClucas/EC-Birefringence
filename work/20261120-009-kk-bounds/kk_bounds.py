"""Numbers for KK-02, KK-03 and the radion x F Ftilde class (2608.09299 and its EC analogue).

All measured inputs come from data/. Writes results.csv next to this file.
"""
import csv
import os

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from common import (ALPHA, A_STAR, BETA, BETA_DEG, BETA_ERR, DEG, E, G4SQ, GEV, H0, H_STAR, HBAR, HERE,
                    MEAS, MP, NEFF_SM, OMG, RNU1, T0, YR, binp, bound_num, check, info, kk, BOUNDS)

OUT = []


def rec(key, value, note=""):
    OUT.append((key, value, note))
    info(f"{key} = {value:.4g}  {note}")


def main():
    g4 = np.sqrt(G4SQ)
    gmax = bound_num("A-Chandra")              # eV^-1, m < 1e-12 eV
    rec("beta_rad", BETA, f"ACT DR6 {BETA_DEG} deg")
    rec("H0_eV", H0)
    rec("H_star_eV", H_STAR, "H at z_star")

    # ---------------- KK-02: Wilson-line axion, quantized CS level ----------------
    # g = k alpha/(pi f_a), beta = (g/2) Delta phi = k alpha Delta theta / (2 pi): independent of L and f_a
    for k in (1, 2, 3):
        bmax = k * ALPHA / 2
        dth = 2 * np.pi * BETA / (k * ALPHA)
        lo = 2 * np.pi * (BETA - BETA_ERR) / (k * ALPHA)
        hi = 2 * np.pi * (BETA + BETA_ERR) / (k * ALPHA)
        rec(f"KK02_k{k}_beta_max_deg", bmax / DEG, "Delta theta = pi")
        rec(f"KK02_k{k}_Dtheta_over_pi", dth / np.pi, f"1 sigma range {lo / np.pi:.2f}-{hi / np.pi:.2f}")
        fmin = k * ALPHA / (np.pi * gmax)
        invL_min = 2 * np.pi * g4 * fmin
        rec(f"KK02_k{k}_fa_min_GeV", fmin / GEV, "Chandra g < 8e-13/GeV")
        rec(f"KK02_k{k}_invL_min_GeV", invL_min / GEV, "1/L = 2 pi g4 f_a")
    invL_max = MP / np.sqrt(2 * np.pi)
    rec("KK02_invL_max_GeV", invL_max / GEV, "1/L < M_5 with M_5^3 (2 pi L) = M_P^2")
    check("K1 k = 1 reaches the ACT central value with Delta theta <= 1.03 pi; k >= 2 needs Delta theta < pi",
          2 * np.pi * BETA / ALPHA / np.pi < 1.1 and 2 * np.pi * BETA / (2 * ALPHA) < np.pi)

    # Hosotani: m = q g4 sqrt(x^2+3x+3) e^{-x/2}/(pi l), l = 2 pi L; solve x = M l for m in [H0, H_star]
    for invL in (2 * np.pi * g4 * ALPHA / (np.pi * gmax), 1e12 * GEV, 1e16 * GEV, invL_max):
        l = 2 * np.pi / invL
        xs = []
        for m in (H0, H_STAR):
            f = lambda x: np.log(g4 * np.sqrt(x**2 + 3 * x + 3) / (np.pi * l)) - x / 2 - np.log(m)
            xs.append(brentq(f, 2.0, 2000.0))
        rec(f"KK02_xMl_window_at_invL_{invL / GEV:.2e}GeV", xs[0], f"to {xs[1]:.1f} (m from H0 to H_star); M = x/(2 pi L)")

    # isotropy: max anisotropic rotation if theta fluctuates by O(pi) across the sky
    acb = float(MEAS["BAL25-ACB-SPT3G"]["value"])
    for k in (1, 3):
        rec(f"KK02_k{k}_max_dbeta2_rad2", (k * ALPHA / 2) ** 2, f"vs SPT A_CB < {acb:g} rad^2 (not binding)")

    # ---------------- stiff (massless, kicked) field from z_star to today ----------------
    I_late, _ = quad(lambda a: 1.0 / (a**4 * E(a) * H0), A_STAR, 1.0, limit=400)   # eV^-1
    # Delta X = Xdot_* a_*^3 I_late ; rho_* = Xdot_*^2 / 2
    def rho_star(dX):
        return 0.5 * (dX / (A_STAR**3 * I_late)) ** 2
    rho_g_star = (np.pi**2 / 15) * (T0 / A_STAR) ** 4
    Neff_P, sNeff_P = kk("Neff_Planck_BAO")
    Nnu, sNnu = kk("Nnu_BBN")
    dN_cmb = Neff_P + 2 * sNeff_P - NEFF_SM
    dN_bbn = Nnu + 2 * sNnu - 3.0
    rec("dNeff_max_CMB", dN_cmb, "Planck+BAO, +2 sigma minus 3.046")
    rec("dNnu_max_BBN", dN_bbn, "Fields et al. +2 sigma minus 3")
    # BBN: a_BBN from entropy conservation (theory input g*s)
    gs0 = 2 + 7 / 8 * 6 * 4 / 11
    gsB = binp("g_star_BBN")
    TB = binp("T_BBN")
    a_B = (gs0 / gsB) ** (1 / 3) * T0 / TB
    rho_1nu_B = 7 / 8 * (np.pi**2 / 15) * TB**4
    # Delta X_max from each bound
    dX_cmb = np.sqrt(dN_cmb * RNU1 * rho_g_star * 2) * A_STAR**3 * I_late
    dX_bbn = np.sqrt(dN_bbn * rho_1nu_B * 2 / (A_STAR / a_B) ** 6) * A_STAR**3 * I_late
    rec("stiff_DeltaX_max_CMB_GeV", dX_cmb / GEV)
    rec("stiff_DeltaX_max_BBN_GeV", dX_bbn / GEV, "binding")
    check("K2 stiff-field round trip: Delta X_max(BBN) reproduces dN_BBN",
          abs(rho_star(dX_bbn) * (A_STAR / a_B) ** 6 / rho_1nu_B / dN_bbn - 1) < 1e-10)

    # ---------------- KK-03: KR axion, g = sqrt6 |lambda| g4^2 M_P ----------------
    gmin = 2 * BETA / min(dX_cmb, dX_bbn)
    rec("KK03_g_min_GeVinv", gmin * GEV, "massless b kicked before BBN, beta = ACT")
    rec("KK03_g_max_GeVinv", gmax * GEV, "Chandra")
    lam_k = lambda g: g * MP / (np.sqrt(6) * G4SQ)  # lambda in units of kappa_4 = 1/M_P^2
    rec("KK03_lambda_over_kappa4_min", lam_k(gmin))
    rec("KK03_lambda_over_kappa4_max", lam_k(gmax))
    rec("KK03_Ms_max_GeV", MP / np.sqrt(lam_k(gmin)) / GEV, "if lambda = 1/M_s^2")
    rec("KK03_Ms_min_GeV", MP / np.sqrt(lam_k(gmax)) / GEV)
    g_nat = np.sqrt(6) * G4SQ / MP
    dX_nat = 2 * BETA / g_nat
    dN_nat = rho_star(dX_nat) * (A_STAR / a_B) ** 6 / rho_1nu_B
    rec("KK03_g_natural_GeVinv", g_nat * GEV, "lambda = kappa_4")
    rec("KK03_natural_dNnu_BBN", dN_nat, "excluded if >> dN_BBN")
    check("K3 KK-03 window nonempty: g_min(BBN) < g_max(Chandra); natural lambda = kappa_4 excluded",
          gmin < gmax and dN_nat > 1e3 * dN_bbn)
    # axial torsion today from the KR background, S ~ g_B bdot_0, g_B = sqrt(2/3)/M_P
    bdot0 = (2 * BETA / gmin) / I_late
    S0 = np.sqrt(2 / 3) / MP * bdot0
    v = binp("v_sun_cmb") / binp("c_kms")
    b_perp_n = bound_num("F-n-bperp-Allm")
    rec("KK03_S0_today_eV_at_gmin", S0, f"S_perp ~ v S0 = {v * S0:.2e} vs neutron |b_perp| < {b_perp_n:.1e} eV")

    # ---------------- radion x F Ftilde (2608 and EC u-Holst) ----------------
    # |Delta alpha/alpha|_{LSS->0} = 2 beta/|c_odd|
    row = BOUNDS["K-alpha-CMB"]["value_eV"].split("+-")
    a_cmb, s_cmb = float(row[0]), float(row[1])
    da_max = abs(a_cmb) + 2 * s_cmb
    rec("radion_codd_min_CMB_alpha", 2 * BETA / da_max, "c_odd needed so that |Delta alpha/alpha|(LSS) < 2 sigma")
    rec("2608_c_min_CMB_alpha", BETA / da_max, "c_odd = 2c")
    # late-time variation if the radion is a massless kicked field (still moving today)
    I_q, _ = quad(lambda a: 1.0 / (a**4 * E(a) * H0), 1 / (1 + 1.15), 1.0, limit=400)
    frac_q = I_q / I_late
    row = BOUNDS["K-alpha-QSO"]["value_eV"].split("+-")
    q_max = abs(float(row[0])) + 2 * float(row[1])
    rec("radion_frac_after_z1.15", frac_q, "share of Delta ln Z after z = 1.15 (stiff field)")
    rec("radion_codd_min_QSO", 2 * BETA * frac_q / q_max)
    row = BOUNDS["K-alphadot-clock"]["value_eV"].split("+-")
    clk = abs(float(row[0])) + 2 * float(row[1])     # eV
    rate = 2 * BETA / I_late                          # c_odd * |alpha_dot/alpha| today
    rec("radion_codd_min_clock", rate / clk)
    # Cassini / ISL: massless radion is Brans-Dicke omega = 0
    alpha0sq = 1 / (2 * 0 + 3)
    gm1 = -2 * alpha0sq / (1 + alpha0sq)
    cass, scass = kk("gamma_minus_1_Cassini")
    rec("radion_gamma_minus_1", gm1, f"Cassini {cass:g} +- {scass:g}")
    m_isl = float(BOUNDS["K-ISL-EotWash"]["value_eV"].replace("m_phi", "").replace(">", "").strip())
    rec("radion_m_min_ISL_eV", m_isl, f"= {m_isl / H_STAR:.2e} H_star")
    check("K4 massless radion excluded by Cassini by > 1e4 sigma; ISL mass >> H_star (cannot roll after LSS)",
          abs(gm1 - cass) / scass > 1e4 and m_isl / H_STAR > 1e20)
    # EC u-Holst c_odd(c) = 16 c^3/(1+12 c^2): c needed for c_odd_min(clock)
    cod = rate / clk
    c_needed = brentq(lambda c: 16 * c**3 / (1 + 12 * c**2) - cod, 1e-3, 1e12)
    rec("EC_uHolst_c_needed_clock", c_needed)

    with open(os.path.join(HERE, "results.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["key", "value", "note"])
        for r in OUT:
            w.writerow([r[0], f"{r[1]:.6g}", r[2]])


if __name__ == "__main__":
    main()
