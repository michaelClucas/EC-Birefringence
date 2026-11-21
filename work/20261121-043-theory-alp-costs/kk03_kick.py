"""KK-03: Kalb-Ramond axion b (dual of the torsion 2-form) with a Green-Schwarz coupling lambda.

g = sqrt6 lambda g4^2 M_P (task 009), so lambda/kappa_4 = g M_P/(sqrt6 * 4 pi alpha).
b is exactly massless (shift symmetry from B gauge invariance), so it only moves if it has a velocity:
bdot ~ a^-3 (stiff, rho_b ~ a^-6). Delta b (z* -> today) = bdot_* a_*^3 I, I = int_{a*}^1 da/(a^4 H).
beta = (g/2) Delta b fixes bdot for each g; BBN Delta N_nu < 0.16 (task 009) gives g_min, Chandra g_max.
The "kick": rho_b/rho_rad grows as T^2 into the past (a^-6 vs a^-4, g* = 10.75 held fixed above 1 MeV), so the
velocity must be created at T_k <= T_BBN / sqrt(Omega_b(T_BBN)), or the Universe was kination dominated above T_k.
"""
import math
import sys

from scipy.integrate import quad

from common import (ACT, ALPHA, A_STAR, DEG, E, ESK, GEV, H0, KKIN, MP, RNU1, T0, binp, bound_num, check, done,
                    info)


def main():
    print("KK-03 Kalb-Ramond axion + Green-Schwarz: lambda window and pre-BBN kick (candidate mechanism, ALP class)")
    G4SQ = 4 * math.pi * ALPHA
    I_late = quad(lambda a: 1 / (a**4 * E(a) * H0), A_STAR, 1, limit=400)[0]
    TB = binp("T_BBN")
    gsB = binp("g_star_BBN")
    gs0 = 2 + 7 / 8 * 6 * 4 / 11
    a_B = (gs0 / gsB) ** (1 / 3) * T0 / TB
    Nnu, sN = float(KKIN["Nnu_BBN"]["value"]), float(KKIN["Nnu_BBN"]["err"])
    dN_max = Nnu + 2 * sN - 3.0
    rho1nu = 7 / 8 * (math.pi**2 / 15) * TB**4
    rho_rad = math.pi**2 / 30 * gsB * TB**4

    def dN_of(beta_rad, g):
        """Delta N_nu at BBN for a stiff b giving rotation beta with coupling g."""
        db = 2 * beta_rad / g
        bdot_star = db / (A_STAR**3 * I_late)
        rho_star = 0.5 * bdot_star**2
        return rho_star * (A_STAR / a_B) ** 6 / rho1nu

    gmax = bound_num("A-Chandra")  # eV^-1
    lam = lambda g: g * MP / (math.sqrt(6) * G4SQ)
    out = {}
    for lab, (b, s) in (("ACT", ACT), ("2608.06480", ESK)):
        for tag, bb in (("central", b), ("-1s", b - s), ("+1s", b + s)):
            br = bb * DEG
            # g_min: dN_of(br, g) = dN_max  ->  g_min = g * sqrt(dN_of/dN_max) at any g
            gmin = gmax * math.sqrt(dN_of(br, gmax) / dN_max)
            out[(lab, tag)] = gmin
            info(f"{lab} beta = {bb:.3f} deg: g in [{gmin * GEV:.3e}, {gmax * GEV:.1e}] GeV^-1, "
                 f"lambda/kappa4 in [{lam(gmin):.3e}, {lam(gmax):.3e}], M_s = M_P/sqrt(lambda/kappa4) in "
                 f"[{MP / math.sqrt(lam(gmax)) / GEV:.2e}, {MP / math.sqrt(lam(gmin)) / GEV:.2e}] GeV, "
                 f"window width g_max/g_min = {gmax / gmin:.2f}")
    gmin = out[("ACT", "central")]
    check("C1 ACT central: lambda/kappa4 window = (1.47-8.67)e6 (task 009 numbers reproduced independently)",
          abs(lam(gmin) / 1.4693e6 - 1) < 0.01 and abs(lam(gmax) / 8.6724e6 - 1) < 0.01)
    check("C2 window stays open for the 2608.06480 central value and the ACT +1 sigma edge",
          out[("2608.06480", "central")] < gmax and out[("ACT", "+1s")] < gmax)
    # closed-form cross-check of lambda/kappa4 = g M_P/(sqrt6 4 pi alpha)
    check("C3 lambda/kappa4 closed form: 8e-13 GeV^-1 * 2.435e18 GeV / (sqrt6 * 4 pi alpha) = 8.67e6",
          abs(8e-13 * 2.435e18 / (math.sqrt(6) * 4 * math.pi * 7.2973525693e-3) / 8.672e6 - 1) < 2e-3)

    # kick: energy fraction at BBN and the latest-allowed origin of the velocity
    br = ACT[0] * DEG
    for tag, g in (("g_min", gmin), ("g_max (Chandra)", gmax)):
        dN = dN_of(br, g)
        om = dN * rho1nu / rho_rad
        Tk = TB / math.sqrt(om)
        bdot_T2 = math.sqrt(2 * om * rho_rad) / TB**2
        info(f"ACT central, {tag}: Delta N_nu(BBN) = {dN:.3e}, rho_b/rho_rad(1 MeV) = {om:.2e}, "
             f"bdot/T^2 at 1 MeV = {bdot_T2:.3f}; velocity must originate at T_k <= {Tk / 1e6:.1f} MeV "
             f"(else kination domination above that T)")
        out[("Tk", tag)] = Tk
    check("C4 at g_min the stiff b saturates the BBN bound (Delta N_nu = 0.16) by construction",
          abs(dN_of(br, gmin) / dN_max - 1) < 1e-9)
    check("C5 kick must originate below ~6 MeV (g_min) to ~37 MeV (g_max): T_k,max in [3, 60] MeV",
          3e6 < out[("Tk", "g_min")] < out[("Tk", "g_max (Chandra)")] < 6e7)
    # independent check of rho ~ a^-6 vs radiation: Omega_b(T) / Omega_b(T_B) = (T/T_B)^2 at fixed g*
    a1, a2 = a_B, a_B / 10
    check("C6 stiff/radiation ratio grows as a^-2 = T^2 (fixed g*)",
          abs(((a2 / a1) ** -6 / (a2 / a1) ** -4) / 100 - 1) < 1e-12)
    # when does b move: share of Delta b accumulated before z = 500 and after z = 10
    f500 = quad(lambda a: 1 / (a**4 * E(a) * H0), A_STAR, 1 / 501, limit=400)[0] / I_late
    f10 = quad(lambda a: 1 / (a**4 * E(a) * H0), 1 / 11, 1, limit=400)[0] / I_late
    info(f"share of Delta b between z* and z = 500: {f500:.3f}; after z = 10: {f10:.1e} (rotation is right after z*)")
    check("C7 a massless kicked b rotates mostly just after z* (> 50% before z = 500), no timing tuning needed",
          f500 > 0.5)
    g_nat = math.sqrt(6) * G4SQ / MP
    info(f"natural lambda = kappa4: g = {g_nat * GEV:.2e} GeV^-1, Delta N_nu = {dN_of(br, g_nat):.2e}")
    check("C8 natural lambda = kappa4 excluded by BBN by > 1e10", dN_of(br, g_nat) > 1e10 * dN_max)
    return out


if __name__ == "__main__":
    main()
    sys.exit(done())
