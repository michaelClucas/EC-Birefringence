"""EC-02: minimal EC spin source + CFJ term -(1/2) xi S_mu A_nu Ftilde^{mu nu}. Derivation Sec. A.

  beta = -(xi/2) int S^0hat dt,  S^0hat = -(kappa/2) n5   =>   beta = +(xi kappa/4) int n5 dt
  n5 = J5^0hat = sum over species of (helicity) x (speed) x (number density)   (check_helicity.py)

Sources of isotropic n5 (all inputs from data/cosmo_inputs.csv):
  nu-A  relic neutrino asymmetry, all 3 states massless (upper envelope), |xi_nu| <= 95% bound
  nu-B  same, normal ordering with minimal masses (m1 = 0)
  nu-C  generous: xi_e at 95% bound, xi_mu = xi_tau saturating the Planck N_eff 95% upper bound, massless
  DM-th fully helicity-polarised thermal-relic fermion DM at the Lyman-alpha mass bound
  DM-P  fully helicity-polarised degenerate (Pauli-limited) fermion DM with the same velocity bound
  e-env every electron with the same helicity (no known process does this; envelope only)
Spin polarisation by a magnetic field gives a spatial S^i (anisotropic), reported separately.
"""
import math
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from common import Cosmo, ZETA3, beta_rad, load, report


def fd(x):
    return 1.0 / (np.exp(x) + 1.0) if x < 700 else 0.0


def dn_closed(T, xi):
    """n - nbar for one Weyl species with degeneracy xi (exact for massless FD)."""
    return T**3 / (6 * math.pi**2) * (math.pi**2 * xi + xi**3)


def n5_species(T, xi, m):
    """|J5^0| for one helicity-asymmetric species: int d^3p/(2pi)^3 [f(xi) - f(-xi)] p/E."""
    g = lambda y: y**2 * (fd(y - xi) - fd(y + xi)) * (y / math.sqrt(y**2 + (m / T) ** 2))
    return T**3 / (2 * math.pi**2) * quad(g, 0, 60, limit=200, epsrel=1e-11)[0]


def dneff(xi):
    x = xi / math.pi
    return 30 / 7 * x**2 + 15 / 7 * x**4


def main():
    c = load()
    cos = Cosmo(c)
    k = cos.kappa
    ok = True
    beta = beta_rad(c)
    hb = c["hbar_eV_s"]
    gyr = cos.t_int(lambda a: 1.0) * hb / 3.15576e16
    ok &= report(f"cosmic time LSS->today = {gyr:.3f} Gyr in [13, 14.5]", 13 < gyr < 14.5)

    # --- checks of the n5 machinery
    T, xi = 1.0, 0.3
    ok &= report("massless n5 integral = closed form (pi^2 xi + xi^3) T^3/6pi^2",
                 abs(n5_species(T, xi, 0.0) / dn_closed(T, xi) - 1) < 1e-8)
    m = 1e4
    num = n5_species(T, xi, m)
    pav = quad(lambda y: y**3 * (fd(y - xi) - fd(y + xi)), 0, 60)[0] / quad(lambda y: y**2 * (fd(y - xi) - fd(y + xi)), 0, 60)[0]
    ok &= report("NR limit n5 -> Delta n <p>/m", abs(num / (dn_closed(T, xi) * pav * T / m) - 1) < 1e-6)
    f3 = lambda a: a**-3
    I_quad = cos.t_int(f3)
    la = np.linspace(math.log(cos.a_star), 0, 200001)
    y = np.array([f3(math.exp(v)) / cos.H(math.exp(v)) for v in la])
    I_trap = float(np.sum((y[1:] + y[:-1]) / 2 * np.diff(la)))
    ok &= report("int a^-3 dt: quad vs trapezoid agree to 1e-6", abs(I_quad / I_trap - 1) < 1e-6)

    # --- inputs
    xi95 = c["xi_nu"] + 2 * c["xi_nu_err"]
    dN95 = c["N_eff_obs"] + 2 * c["N_eff_obs_err"] - c["N_eff_SM"]
    xi_gen = brentq(lambda x: dneff(x) - dN95 / 2, 0, 5)  # mu and tau share the N_eff room (xi_e term negligible)
    masses = [0.0, math.sqrt(c["dm21_sq"]), math.sqrt(c["dm31_sq"])]
    print(f"\ninputs: beta = {c['beta_deg']} deg; |xi_nu| <= {xi95:.3f} (95%, 2409.08280 T4); "
          f"Delta N_eff(95%) <= {dN95:.3f} -> generous xi_mu,tau <= {xi_gen:.2f}; masses NO-min {['%.3g' % x for x in masses]} eV")

    Tnu = lambda a: cos.Tnu0 / a
    sources = {}
    sources["nu-A massless, xi95"] = lambda a: 3 * dn_closed(Tnu(a), xi95)
    sources["nu-B NO m1=0, xi95"] = lambda a: sum(n5_species(Tnu(a), xi95, mi) for mi in masses)
    sources["nu-C generous N_eff, massless"] = lambda a: dn_closed(Tnu(a), xi95) + 2 * dn_closed(Tnu(a), xi_gen)

    # thermal-relic WDM velocity at the Lyman-alpha bound (g = 2 FD relic with the observed density)
    mX = c["m_WDM_min_eV"]
    nX0 = cos.rho_dm0 / mX
    # n = (3 zeta3 / (4 pi^2)) g T^3 for fermions, g = 2  ->  T = (n 4 pi^2 / (6 zeta3))^(1/3)
    TX0 = (nX0 * 4 * math.pi**2 / (3 * ZETA3 * 2)) ** (1 / 3)
    pmean = 7 * math.pi**4 / (180 * ZETA3)  # <p>/T for FD
    v0 = pmean * TX0 / mX  # mean speed today (NR)
    vstar = v0 / cos.a_star
    ok &= report(f"WDM mean speed at LSS {vstar:.2e} << 1 (NR treatment valid)", vstar < 1e-3)
    sources["DM-th 5.3 keV relic, all same helicity"] = lambda a: cos.rho_dm0 / a**3 * (v0 / a) / mX
    # Pauli-limited: degenerate one-helicity gas, <v> = (3/4) p_F/m <= v_max, n = p_F^3/(6 pi^2)
    mmin = ((3 / 4) * (6 * math.pi**2 * cos.rho_dm0) ** (1 / 3) / v0) ** (3 / 4)
    sources["DM-P degenerate, Pauli-limited"] = lambda a: cos.rho_dm0 / a**3 * (v0 / a) / mmin
    print(f"  WDM: T_X0 = {TX0:.3e} eV, <v>_0 = {v0:.2e} (={v0 * 2.998e5:.3f} km/s), Pauli m_min = {mmin:.3e} eV")

    me, mp_ = c["m_e_eV"], c["m_p_eV"]
    ne = lambda a: cos.rho_b0 / mp_ / a**3  # n_e <= n_b (fully ionised envelope)
    ve = lambda a: math.sqrt(8 * cos.T0 / a / (math.pi * me))  # T_e = T_gamma (overestimate after z ~ 150)
    sources["e-env all electrons one helicity"] = lambda a: ne(a) * ve(a)

    print("\n source                               | n5(LSS) eV^3 | |S^0hat|(LSS) eV | max beta(xi=1) deg | shortfall (log10) | xi needed | Lambda_eff = M_P sqrt(2/xi) GeV")
    res = {}
    for name, f in sources.items():
        I = cos.t_int(f)
        bmax = k / 4 * I  # rad at xi = 1
        xi_need = beta / bmax
        lam = cos.MP * math.sqrt(2 / xi_need) / 1e9
        res[name] = (bmax, xi_need)
        print(f" {name:37s}| {f(cos.a_star):.2e}     | {k / 2 * f(cos.a_star):.2e}         | {math.degrees(bmax):.2e}           | {math.log10(xi_need):6.1f}            | {xi_need:.1e}   | {lam:.1e}")

    # sanity: dominant contribution of the a^-3 integral comes from near LSS
    late = cos.t_int(f3, a1=0.01) / I_quad
    print(f"\n fraction of int a^-3 dt from z < 99: {late:.3f}")

    # --- magnetic-field spin polarisation of free electrons (spatial S, anisotropic)
    B0T = c["B_1Mpc_nG"] * 1e-13  # nG -> T (comoving)
    for aa, lab in [(cos.a_star, "LSS")]:
        P = c["mu_B_eV_T"] * B0T / aa**2 / (cos.T0 / aa)
        Ssp = k / 2 * ne(aa) * P
        print(f" B-field electron spin polarisation at {lab}: P = {P:.1e}, |S^i| <= {Ssp:.1e} eV (random direction per coherence patch; isotropic part 0)")

    # --- cross-check with Reig 2611.07578 Eq. (3): L = -(1/2)(1/Lambda^2) j.A Ftilde, j^0 = Delta n_nu_e (massless)
    xiR = 0.05
    IR = cos.t_int(lambda a: dn_closed(Tnu(a), xiR))
    lamR = math.sqrt(IR / (2 * beta)) / 1e12
    lamReig = c["Reig_Lambda_star_TeV_per_sqrt_xi"] * math.sqrt(xiR)
    print(f"\n Reig mapping: our integral gives Lambda* = {lamR:.1f} TeV for xi_nu_e = 0.05 and beta = {c['beta_deg']} deg; Reig Eq. 3: {lamReig:.1f} TeV")
    ok &= report("our time integral reproduces Reig's Lambda* within a factor 2 (different beta target, Yin-Komatsu fit)",
                 0.5 < lamR / lamReig < 2)

    # --- what bounds xi: KR/dual-axion completion (Reig) => g_ag g_ae = xi kappa m_e
    bound = c["helioscope_gag_gae_GeV"] * 1e-9  # eV^-1
    xi_max = bound / (k * me)
    best = max(res.values(), key=lambda t: t[0])
    print(f" helioscope bound via dual axion: xi < {xi_max:.1e}; best-case needed xi = {best[1]:.1e} -> excluded by {math.log10(best[1] / xi_max):.1f} orders")
    print(f" lab: minimal-EC |S^0hat| today from nu-A = {k / 2 * sources['nu-A massless, xi95'](1.0):.1e} eV vs lab bound {c['S0_lab_max']:.1e} eV (irrelevant)")
    ok &= report("no source reaches beta at xi = 1 (no-go holds for every listed source)",
                 all(v[1] > 1 for v in res.values()))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
