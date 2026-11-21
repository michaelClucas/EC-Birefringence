"""EC-06a remnant: Magueijo-Zlosnik Lambda-field with a hand-added f(Lambda) F Ftilde, f = (xi/gamma) ln Lambda.

Only the exact-gradient vacuum piece S_vac = (1/gamma) d ln Lambda = c_s S couples (skeptic 034, K7).
beta = (xi gamma / 2) J_alt with J_alt = int_{a*}^1 c_s eps_tr dln a, eps_tr the tracker (031 T2a closed form),
so xi gamma is fixed by beta and the vacuum ALP coupling g = xi/(sqrt3 M_P) = (xi gamma)/(gamma sqrt3 M_P) is set by
whichever upper bound on gamma holds. Every gamma bound is an upper bound, so every g is a lower limit.
Chandra: (a) naive transfer, g < 8e-13 GeV^-1; (b) estimate with the field's own dispersion omega = c_s k
(031: L = 1/2 phidot^2 - 1/2 c_s^2 (grad phi)^2): photon-phi momentum mismatch q = omega (1/c_s - 1), so the
conversion probability is at most (2 g_eff B / (2 q))^2 instead of (g_eff B L/2)^2; the bound on g_eff weakens by
q L/2. c_s <= its cosmic mean today in any overdense region (c_s ~ Lambda/(Lambda + kappa rho)).
"""
import math
import sys

from scipy.integrate import quad

from common import (A_STAR, ACT, ALPC, DEG, E, EC06T, ESK, GEV, H0, HBAR, HBARC, MP, MPC, MZIN, OML, OMM, OMR, RNU1,
                    bound_num, check, done, info, load)


def main():
    print("EC-06a remnant: MZ Lambda-field with f(Lambda) F Ftilde (candidate mechanism, ALP class)")
    Lam = 3 * OML

    def eps_tr(a):
        I3 = quad(lambda b: b**2 / E(b), 0, a, limit=200)[0]
        return 2 / 3 - 2 / 9 * Lam**2 * I3 / (a**3 * E(a) ** 3)

    cs_bg = lambda a: OML / (OML + OMM / a**3 + OMR / a**4)
    J = quad(lambda a: eps_tr(a) / a, A_STAR, 1, limit=200)[0]
    Jalt = quad(lambda a: cs_bg(a) * eps_tr(a) / a, A_STAR, 1, limit=200)[0]
    info(f"J = {J:.3f}, J_alt = {Jalt:.3f} (034: 4.584, 0.216)")
    check("M1 tracker integrals reproduce skeptic 034 (J = 4.58, J_alt = 0.216) to 2%",
          abs(J / 4.584 - 1) < 0.02 and abs(Jalt / 0.216 - 1) < 0.02)
    late = quad(lambda a: cs_bg(a) * eps_tr(a) / a, 0.5, 1, limit=200)[0] / Jalt
    info(f"share of the rotation after z = 1: {late:.2f} (all CMB photons see the full Delta, no l dependence at l > 100)")

    eps0 = eps_tr(1.0)
    cos = load("cosmo_inputs.csv")
    Nsm = float(cos["N_eff_SM"]["value"])
    Nlo = float(cos["N_eff_obs"]["value"]) - 2 * float(cos["N_eff_obs"]["err"])
    gam = {
        "N_eff only (no GW bound)": math.sqrt(9 * RNU1 * (Nsm - Nlo) / (1 + RNU1 * Nsm)),
        "GW170817 speed (formula used outside its k* regime)":
            math.sqrt(float(EC06T["cT_GW170817_hi"]["value"]) / float(ALPC["MZ-cT-coeff"]["value"])),
        "k* >> k_LIGO at 24 Hz (2001.06373 consistency)":
            H0 / (2 * math.pi * float(MZIN["GW170817_f_low"]["value"]) * HBAR) / eps0,
        "2001.06373 Eq. 70 (|gamma| << 1e-21)": float(ALPC["MZ-gamma-Eq70"]["value"]),
    }
    gch = bound_num("A-Chandra")
    res = {}
    for lab, (b, s) in (("ACT", ACT), ("2608.06480", ESK)):
        xig = 2 * b * DEG / Jalt
        info(f"{lab} beta = {b} deg: xi*gamma = {xig:.3e}")
        for k, gm in gam.items():
            g = xig / gm / (math.sqrt(3) * MP)
            res[(lab, k)] = g
            info(f"   {k}: gamma <= {gm:.2e} -> g >= {g * GEV:.2e} GeV^-1 "
                 f"({'passes' if g <= gch else 'fails'} naive Chandra 8e-13, factor {g / gch:.1e})")
    keys = list(gam)
    check("M2 g lower limits reproduce 034: 1.2e-20, 7.4e-14, 0.24, 8.2 GeV^-1 (ACT) to 5%",
          all(abs(res[("ACT", k)] * GEV / v - 1) < 0.05 for k, v in zip(keys, (1.2e-20, 7.4e-14, 0.24, 8.2))))
    check("M3 naive Chandra: N_eff and GW-speed branches pass, both 2001.06373-consistent branches fail by > 1e11",
          res[("ACT", keys[0])] < gch and res[("ACT", keys[1])] < gch and res[("ACT", keys[2])] > 1e11 * gch)
    check("M4 2608.06480 shifts every g up by 0.277/0.215 = 1.29 (linear in beta)",
          abs(res[("2608.06480", keys[0])] / res[("ACT", keys[0])] - ESK[0] / ACT[0]) < 1e-12)

    # conversion suppression estimate
    w = float(ALPC["conv_omega_min"]["value"])
    L = float(ALPC["conv_L_min"]["value"]) * MPC / 1e3 / HBARC   # kpc -> eV^-1
    cs0 = cs_bg(1.0)
    qL = w * L * (1 / cs0 - 1)
    relax = qL / 2
    info(f"c_s(today, cosmic mean) = {cs0:.3f}; q L = omega L (1/c_s - 1) = {qL:.2e} at omega = 1 keV, L = 1 kpc; "
          f"conversion suppressed by (2/qL)^2 = {(2 / qL) ** 2:.1e}; Chandra bound on g_eff relaxes to "
          f"{gch * relax * GEV:.1e} GeV^-1")
    qL_min = 1.0 * (1 / HBARC) * (1 / cs0 - 1)   # omega = 1 eV, L = 1 m
    check("M5 even at omega = 1 eV and L = 1 m, q L >> 1 (no coherent conversion; Chandra does not bound g)",
          qL_min > 1e5, f"q L = {qL_min:.1e}")
    check("M6 with the dispersion estimate the Eq. 70 branch (g >= 8 GeV^-1) is below the relaxed Chandra bound",
          res[("ACT", keys[3])] < gch * relax)
    return res


if __name__ == "__main__":
    main()
    sys.exit(done())
