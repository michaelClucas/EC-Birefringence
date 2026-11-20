"""E1-E8: electron EDM bound on the KK-05 coupling c, and the other bounds it is compared with.

From check_dipole.py (D5, D6), per unit W = k F_raw, for a 5D Dirac zero mode with real bulk mass:
  EDM coefficient D_E = (3/2)|1/6 - (4/3)c^2|,   MDM coefficient D_B = (3/2)c.
Canonical photon: F_raw = F/sqrt(Z), Z = 1 + 12c^2 (D3). So in physical units
  d(c) = k (3/2)|1/6 - (4/3)c^2| / sqrt(1 + 12c^2),   k = sqrt(2 kappa_4) = sqrt(2)/M_P,
and d/e in cm uses e = sqrt(4 pi alpha). Large c: d -> (2/sqrt(12)) c k = 0.577 c k (grows linearly).
If the electron mass has a P'-odd phase theta (mass term m' Psibar u-slash Psi, allowed only by the oriented u):
  d(c, theta) = k/sqrt(Z) |D_E(c) cos(theta) +- D_B(c) sin(theta)|.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from common import (ALPHA, BETA_ACT, H0, HBAR, HBARC, MP, OM, YR, bound_pair, check, cin, edmin, info, shu)

beta = BETA_ACT
k = np.sqrt(2.0) / MP                       # eV^-1
e = np.sqrt(4 * np.pi * ALPHA)
CM = HBARC * 100.0                          # 1 eV^-1 in cm
d_JILA = shu("eEDM_JILA")                   # e cm, 90% CL


def DE(c):
    return 1.5 * abs(1 / 6 - 4 / 3 * c ** 2)


def DB(c):
    return 1.5 * c


def d_ecm(c, theta=0.0, sgn=1):
    Z = 1 + 12 * c ** 2
    return k / np.sqrt(Z) * abs(DE(c) * np.cos(theta) + sgn * DB(c) * np.sin(theta)) / e * CM


# E1: c = 0 KK Pauli EDM, and the skeptic's (c-odd only) formula for comparison
info(f"E1 pure KK (c = 0) Pauli EDM of a real-mass vector-like zero mode: d = k/4 = {d_ecm(0.0):.2e} e cm")
d_skep = lambda c: 0.75 * c * k / np.sqrt(1 + 12 * c ** 2) / e * CM
info(f"E1b c-odd piece alone (the 015 formula (3/4)ck/sqrt(Z)): {d_skep(1e3):.2e} e cm at c = 1e3 (saturates); "
     f"this piece is the MDM (D6), not the EDM")
check("E1 015 formula reproduced for the c-odd (magnetic) piece: (3/4)k/sqrt(12)/e ~ 8e-33 e cm",
      abs(d_skep(1e6) / 8.2e-33 - 1) < 0.05, f"{d_skep(1e6):.3e}")

# E2: EDM bound on c (theta = 0, minimal coupling)
c_edm = brentq(lambda c: d_ecm(c) - d_JILA, 1.0, 1e6)
slope = d_ecm(1e4) / 1e4
info(f"E2 large-c slope d/c = {slope:.3e} e cm (expected 0.577 k/e = {2 / np.sqrt(12) * k / e * CM:.3e})")
check("E2 EDM bound (JILA 4.1e-30 e cm, theta = 0): c <= c_EDM", 100 < c_edm < 300, f"c_EDM = {c_edm:.1f}")

# E3: MICROSCOPE (012 formula, |d_e| from data), field excursion in units of M_P
d_e = shu("d_e_MICROSCOPE")
c_mic = lambda x: np.sqrt(2) * beta / (d_e * x)
x_crit = np.sqrt(2) * beta / (d_e * c_edm)
for x in (0.39, 1.0, 2.0):
    info(f"E3 MICROSCOPE: Delta chi = {x} M_P -> c >= {c_mic(x):.0f}  (EDM allows c <= {c_edm:.0f})")
check("E3 EDM + MICROSCOPE close the window for Delta chi < x_crit M_P (008 cap excursion 0.39 M_P is inside)",
      0.39 < x_crit < 1.0, f"x_crit = {x_crit:.2f}")

# E4: clocks (free field on LCDM; same construction as 012 B3)
h = cin("H0") / 100
Tg = cin("T_CMB") * cin("k_B")
rho_g = np.pi ** 2 / 15 * Tg ** 4
rho_c = cin("rho_crit_h2") * 1e9 * h ** 2 * (HBARC * 100) ** 3
OR = rho_g / rho_c * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * cin("N_eff"))
OL = 1 - OM - OR
Hof = lambda a: H0 * np.sqrt(OR / a ** 4 + OM / a ** 3 + OL)
dlnH = lambda a: -0.5 * (4 * OR / a ** 4 + 3 * OM / a ** 3) / (OR / a ** 4 + OM / a ** 3 + OL)
z_star = cin("z_star")


def run(m):
    def rhs(N, y):
        a = np.exp(N)
        return [y[1], -(3 + dlnH(a)) * y[1] - (m / Hof(a)) ** 2 * y[0]]
    sol = solve_ivp(rhs, (np.log(1e-6), 0.0), [1.0, 0.0], method="DOP853", rtol=1e-10, atol=1e-13, dense_output=True)
    chi0, chip0 = sol.y[0, -1], sol.y[1, -1]
    dchi = abs(sol.sol(-np.log1p(z_star))[0] - chi0)
    env = np.sqrt((H0 * chip0) ** 2 + (m * chi0) ** 2)
    return env / (H0 * dchi), abs(H0 * chip0) / (H0 * dchi)


mc, sc = bound_pair("K-alphadot-clock")
adot_max = abs(mc) + 2 * sc
Rs = {m: run(m) for m in (1e-33, 6e-33, 1e-31, 1e-30, 1e-29)}
c_clock = min(beta * H0 * R[0] / adot_max for R in Rs.values())
info(f"E4 clocks (|alpha-dot/alpha| <= {adot_max / HBAR * YR:.1e}/yr): c >= {c_clock:.2e} over m = 1e-33..1e-29 eV")
check("E4 clock bound exceeds the EDM bound (window closed for a field still moving today)", c_clock > c_edm,
      f"ratio {c_clock / c_edm:.0f}")

# E5: theta escape (non-minimal P'-odd electron mass)
for sgn in (1, -1):
    th = brentq(lambda t: d_ecm(c_clock, t, sgn) - d_JILA, 0.0, np.pi / 2)
    info(f"E5 sign {sgn:+d}: at c = c_clock the EDM is below JILA only for |cos theta| <~ {abs(np.cos(th)):.1e} "
         f"(theta within {abs(np.pi / 2 - th):.1e} rad of pi/2)")
d_half = d_ecm(c_clock, np.pi / 2)
check("E5 at theta = pi/2 (electron mass entirely from Psibar u-slash Psi) the EDM is the c-odd piece, below JILA",
      d_half < d_JILA, f"{d_half:.1e} e cm")

# E6: rolling chi -> static axial background b_0 = 1.5 c' chi-dot (D8); c-dot_0 = (beta/2) H0 R_vel
Rvel = max(R[1] for R in Rs.values())
cdot0 = beta / 2 * H0 * Rvel                                       # eV (large-c: d beta = 2 dc)
b0 = 1.5 * cdot0
bT, sT = bound_pair("F-e-bT")
bT_lim = 2 * (abs(bT) + 2 * sT)                                    # torsion-only b~_T = 2 x listed combination
info(f"E6 b_0(electron) = 1.5 c-dot ~ {b0:.1e} eV vs Eot-Wash b~_T limit {bT_lim:.1e} eV")
check("E6 axial background from the rolling chi is unobservable (> 15 orders below the lab limit)", b0 < 1e-15 * bT_lim,
      f"{b0 / bT_lim:.0e}")

# E7: time variation of the EDM today vs the oscillating-EDM band
rate = cdot0 / c_edm / HBAR                                        # s^-1, fractional d-dot/d at c = c_EDM
fmin = edmin("osc_EDM_HfF_fmin")
info(f"E7 (1/d) dd/dt today <= {rate * YR:.1e} per yr; HfF+ oscillating-EDM band starts at {fmin:.1e} Hz")
check("E7 the EDM drift rate is far below the lowest searched frequency: static EDM for lab purposes",
      rate < 1e-6 * 2 * np.pi * fmin, f"rate/2 pi f_min = {rate / (2 * np.pi * fmin):.0e}")

info(f"SUMMARY c_EDM = {c_edm:.0f}; MICROSCOPE c >= {c_mic(0.39):.0f} (0.39 M_P), {c_mic(1.0):.0f} (1 M_P); "
     f"clocks c >= {c_clock:.2e}; excluded by factor {c_clock / c_edm:.0f} (clocks) and {c_mic(0.39) / c_edm:.1f} (MICROSCOPE at 0.39 M_P)")
