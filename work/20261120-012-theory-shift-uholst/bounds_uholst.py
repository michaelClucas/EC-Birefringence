"""Question B: bounds on the hidden-scalar u-Holst candidate mechanism (c(chi) u^e^e^R in 5D EC-KK).

Exact relation from check_uholst_chi.py U5:  d beta = -c d ln alpha  (pointwise in c).
With c roughly constant over the roll (large-c regime, Delta c = beta/2 << c):
  |Delta ln alpha|_{LSS->0} = beta/c,   d ln alpha/d chi = beta/(c Delta chi).
B1  CMB alpha (Planck z ~ 1e3):           c >= beta / (|mu| + 2 sigma)
B2  MICROSCOPE d_e (alpha(chi) = alpha(1 + d_e chi/(sqrt2 M_P)), M_P reduced):  c >= sqrt2 beta M_P/(d_e Delta chi)
B3  clocks (alpha-dot today) and quasar (z = 1.15): field ODE on LCDM for a free field of mass m
    starting frozen; R_clock = |chi-dot_0|_env/(H0 Delta chi), R_qso = |chi(z=1.15) - chi_0|/Delta chi
B4  EDM order-of-magnitude (from U3: P'-odd dipole source K_[ab5] = c W~; Pauli normalisation of task 004):
    d ~ (3/4) c k e-units with k = sqrt(2 kappa_4); only meaningful if the electron couples to the KK photon.
"""
import numpy as np
from scipy.integrate import solve_ivp

from common import (ALPHA, BETA_ACT, H0, HBAR, HBARC, MP, OM, YR, binp, bound_pair, check, cin, info, new)

beta = BETA_ACT
info(f"beta (ACT DR6 baseline) = {beta:.4e} rad")

# B1
mu, sig = bound_pair("K-alpha-CMB")
c_cmb = beta / (abs(mu) + 2 * sig)
info(f"B1 CMB alpha: |dalpha/alpha| <= {abs(mu) + 2 * sig:.2e} (2 sigma) -> c >= {c_cmb:.3f}")

# B2
d_e = new("d_e_MICROSCOPE")
c_ep = {x: np.sqrt(2) * beta / (d_e * x) for x in (0.1, 0.39, 1.0)}
for x, v in c_ep.items():
    info(f"B2 MICROSCOPE |d_e| < {d_e:g}: Delta chi = {x} M_P -> c >= {v:.3g}")

# B3: free field on LCDM
h = cin("H0") / 100
Tg = cin("T_CMB") * cin("k_B")                       # eV
rho_g = np.pi ** 2 / 15 * Tg ** 4                     # eV^4
rho_c = cin("rho_crit_h2") * 1e9 * h ** 2 * (HBARC * 100) ** 3   # GeV/cm^3 -> eV^4
OR = rho_g / rho_c * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * cin("N_eff"))
OL = 1 - OM - OR
Hof = lambda a: H0 * np.sqrt(OR / a ** 4 + OM / a ** 3 + OL)
dlnH = lambda a: -0.5 * (4 * OR / a ** 4 + 3 * OM / a ** 3) / (OR / a ** 4 + OM / a ** 3 + OL)
z_star = cin("z_star")


def run(m):
    def rhs(N, y):
        a = np.exp(N)
        return [y[1], -(3 + dlnH(a)) * y[1] - (m / Hof(a)) ** 2 * y[0]]
    N0 = np.log(1e-6)
    sol = solve_ivp(rhs, (N0, 0.0), [1.0, 0.0], method="DOP853", rtol=1e-10, atol=1e-13, dense_output=True)
    chi = lambda z: sol.sol(-np.log1p(z))[0]
    chi0, chip0 = sol.y[0, -1], sol.y[1, -1]
    dchi = abs(chi(z_star) - chi0)
    env = np.sqrt((H0 * chip0) ** 2 + (m * chi0) ** 2)   # max |chi-dot| over the phase of today's oscillation
    vel = abs(H0 * chip0)
    return dchi, env / (H0 * dchi), vel / (H0 * dchi), abs(chi(1.15) - chi0) / dchi


mc, sc = bound_pair("K-alphadot-clock")                 # eV (per-year rate converted to eV in bounds.csv)
adot_max = abs(mc) + 2 * sc                              # eV
mq, sq = bound_pair("K-alpha-QSO")
qso_max = abs(mq) + 2 * sq
info(f"clock |alpha-dot/alpha| <= {adot_max / HBAR * YR:.2e}/yr; QSO |dalpha/alpha|(z=1.15) <= {qso_max:.2e}")
rows = []
for m in (1e-33, 6e-33, 1e-31, 1e-30, 1e-29):
    dchi, Renv, Rvel, Rq = run(m)
    c_clk_env = beta * H0 * Renv / adot_max
    c_clk_now = beta * H0 * Rvel / adot_max
    c_q = beta * Rq / qso_max
    rows.append((m, Renv, c_clk_env, c_clk_now, c_q))
    info(f"m={m:.0e} eV: Delta chi/chi_i={dchi:.3f}  R_clock(env)={Renv:.3f}  c_clock >= {c_clk_env:.2e} "
         f"(this phase {c_clk_now:.2e})  R_qso={Rq:.3f}  c_qso >= {c_q:.2e}")
c_clock_min = min(r[2] for r in rows)
check("B3 clock bound on c is m-independent within an O(1) factor over 1e-33..1e-29 eV (matter-era scaling)",
      max(r[2] for r in rows) / c_clock_min < 30, f"min {c_clock_min:.2e}, max {max(r[2] for r in rows):.2e}")

# B4: EDM order of magnitude
k = np.sqrt(2) / (MP / 1e9)                              # GeV^-1, k = sqrt(2 kappa_4)
e = np.sqrt(4 * np.pi * ALPHA)
gev_cm = HBARC * 100 / 1e9                              # 1 GeV^-1 in cm
d_per_c = 0.75 * k * gev_cm / e                          # e cm per unit c
c_edm = new("eEDM_JILA") / d_per_c
info(f"B4 EDM estimate d ~ {d_per_c:.2e} c e cm -> c <~ {c_edm:.0f} (JILA 4.1e-30 e cm); order of magnitude only")

c_need = max(c_cmb, c_ep[0.39], c_clock_min)
info(f"Required c (Delta chi = 0.39 M_P, field moving today): c >= {c_need:.2e}; EDM estimate allows c <~ {c_edm:.0f}")
check("B5 window closed if the EDM estimate holds (c_needed > c_EDM)", c_need > c_edm,
      f"{c_need:.2e} vs {c_edm:.0f}")
check("B6 MICROSCOPE alone needs c > 100 at Delta chi = 0.39 M_P", c_ep[0.39] > 100, f"{c_ep[0.39]:.0f}")
