"""EC-04: order-of-magnitude cosmology of the PGT 0- mode (derivation Sec. C).

From check_pgt_0minus.py: L = (A/2) S.S + (w/2)(d.S)^2 + (3/4) S.J5 reduces to a pseudoscalar
phi = Phi/sqrt(A) with mass^2 = A/w and S_mu = d_mu phi / sqrt(A) (+ contact term -(3/4) J5/A).
Minimal-EC normalisation A = 3/(2 kappa) (PGT can shift A by its axial-torsion-squared coefficient).
  fermion coupling (3/4)(d phi.J5)/sqrt(A) = (1/2f) d phi.J5 with f = (2/3) sqrt(A)  -> g_ae = m_e/f
  photon: no tree coupling; with the CFJ term (EC-01) -(1/2) xi S.K -> g_phigamma = -xi/sqrt(A), beta = -(g/2) Delta phi
FRW: phi'' + 3 H phi' + m^2 phi = source(div J5) ~ 0;  S^0hat = -phi'/sqrt(A) (sign immaterial here).
"""
import csv
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

from common import ROOT, Cosmo, beta_rad, load, report


def cast_bound():
    with open(os.path.join(ROOT, "data", "bounds.csv"), newline="", encoding="utf-8") as f:
        r = {row["id"]: row for row in csv.DictReader(f)}["A-CAST"]
    return float(r["value_eV"].replace("<", "").strip())  # eV^-1


def evolve(cos, mu, a_start=None):
    """phi(N), phi frozen = 1 at LSS; mu = m/H0. Returns (a, phi, dphi/dt in H0 units)."""
    a0 = cos.a_star if a_start is None else a_start
    E = lambda a: cos.H(a) / cos.H0

    def rhs(N, y):
        a = math.exp(N)
        e = E(a)
        return [y[1] / e, (-3 * e * y[1] - mu**2 * y[0]) / e]

    sol = solve_ivp(rhs, [math.log(a0), 0.0], [1.0, 0.0], method="DOP853", rtol=1e-10, atol=1e-13, dense_output=True)
    return sol


def main():
    c = load()
    cos = Cosmo(c)
    ok = True
    beta = beta_rad(c)
    MP, H0, me = cos.MP, cos.H0, c["m_e_eV"]
    A = 1.5 * MP**2
    f = 2 / 3 * math.sqrt(A)
    print(f"A = 3/(2 kappa) = {A:.3e} eV^2;  f = (2/3) sqrt(A) = {f / 1e9:.2e} GeV;  g_ae = m_e/f = {me / f:.1e}")
    g1 = 1 / math.sqrt(A)
    print(f"with CFJ xi-term: |g_phigamma| = xi/sqrt(A) = {g1 * 1e9:.2e} xi GeV^-1;  needed |xi Delta theta| = 2 beta = {2 * beta:.2e}")
    print(f"needed |Delta phi| = 2 beta sqrt(A)/xi = {2 * beta * math.sqrt(A) / 1e9:.2e} GeV / xi")
    gcast = cast_bound()
    print(f"CAST g < {gcast * 1e9:.1e} GeV^-1  ->  xi < {gcast / g1:.1e}")

    # massless check: a^3 phidot = const
    sol = evolve(cos, 0.0)
    ok &= report("m = 0: phi stays frozen if started at rest", abs(sol.y[0][-1] - 1) < 1e-10)
    # kinetic solution check: start with phidot, m = 0
    E = lambda a: cos.H(a) / cos.H0
    s2 = solve_ivp(lambda N, y: [y[1] / E(math.exp(N)), -3 * y[1]], [math.log(cos.a_star), 0], [0.0, 1.0],
                   method="DOP853", rtol=1e-11, atol=1e-14)
    ok &= report("m = 0 kinetic mode: a^3 phidot conserved", abs(s2.y[1][-1] / cos.a_star**3 - 1) < 1e-6)

    print("\n thawing 0- mode (frozen until H ~ m), Delta theta/theta_i by today, needed theta_i (xi = 1), Omega_phi, w = A/m^2, |S0| today")
    lab_ok = True
    for mu in (0.3, 1.0, 2.0, 5.0):
        sol = evolve(cos, mu)
        dphi = sol.y[0][-1] - 1.0
        theta_i = 2 * beta / abs(dphi)
        Om = 0.5 * (mu * H0) ** 2 * A * theta_i**2 / cos.rho_c  # frozen potential energy at start
        S0 = abs(sol.y[1][-1]) * H0 * theta_i  # |theta-dot| today
        print(f"  m = {mu:3.1f} H0: Delta/theta_i = {dphi:+.3f}, theta_i = {theta_i:.2e}, Omega_phi ~ {Om:.1e}, w = {A / (mu * H0) ** 2:.1e}, |S0| = {S0:.1e} eV")
        lab_ok &= S0 < c["S0_lab_max"]
    ok &= report("lab bound |S_0| < 7.3e-22 eV satisfied in all thawing cases at xi = 1", lab_ok)

    # kination alternative: massless, phidot ~ a^-3 from LSS
    Ia = cos.t_int(lambda a: (cos.a_star / a) ** 3)
    thdot_star = 2 * beta / Ia
    rho_k = 0.5 * A * thdot_star**2
    rho_tot = 3 * cos.H(cos.a_star) ** 2 * MP**2
    print(f"\n massless kinetic mode: theta-dot(LSS) = {thdot_star:.1e} eV, rho/rho_tot at LSS = {rho_k / rho_tot:.1e} (xi = 1), scales as 1/xi^2")
    xi_min = 2 * beta * H0 / c["S0_lab_max"]
    print(f" lower bound on xi from lab |S_0| (thawing, S0 ~ H0 Delta theta): xi > ~{xi_min:.0e}")
    print(f" naturalness: m <~ H0 needs w = A/m^2 >~ {A / H0**2:.1e} (dimensionless coefficient of (d.S)^2)")
    # photon coupling without xi: derivative coupling through charged-fermion loop, decoupled as m_phi^2/(12 m_f^2)
    print(f" no xi-term: loop g_eff ~ (alpha/pi f) (m_phi/m_e)^2/12 = {c['alpha_em'] / math.pi / f * (H0 / me) ** 2 / 12 * 1e9:.1e} GeV^-1 for m_phi = H0")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
