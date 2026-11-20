"""KK-05 relations and bounds, rechecked independently of work/20261120-012.

K1  Z = 1 + 12 c^2, Y = 16 c^3 (012 U1, taken as input), alpha ~ 1/Z, beta = (1/2) int dY/Z
    => d beta / d ln alpha = -c exactly; beta(c) closed form. Also: the relation only needs dY/dZ = 2c.
K2  MICROSCOPE: alpha(chi) = alpha(1 + d_e chi/(sqrt2 M_P)); |Delta ln alpha| = beta/c over Delta chi
    => c >= sqrt2 beta M_P/(d_e Delta chi).
K3  clocks, analytic matter-era scaling: chi-dot_0 ~ H0 Delta chi => c >~ beta H0 / |alpha-dot/alpha|_max (R = 1).
K4  EDM: 012 uses d = (3/4) c k (k = sqrt(2 kappa_4)) against the UNNORMALISED KK field strength.
    The canonical photon is A_c = sqrt(Z0 Z) A with Z = 1 + 12 c^2 (Z0 = 1 for k = sqrt(2 kappa_4)), and the
    measured alpha refers to A_c. So d_phys = (3/4) c k / sqrt(1 + 12 c^2), which saturates at (3/4) k/sqrt(12).
    Same rescaling for the P'-even (4/3) c^2 W piece (an anomalous magnetic moment): mu_phys = c^2 k/sqrt(1+12c^2).
"""
import numpy as np
import sympy as sp

from common15 import NPASS, check, info, pair, table, val

c = sp.symbols("c", positive=True)
Z, Y = 1 + 12 * c ** 2, 16 * c ** 3
dbeta = sp.Rational(1, 2) * sp.diff(Y, c) / Z
dlnalpha = -sp.diff(Z, c) / Z
check("K1a d beta/d ln alpha = -c (from Z = 1+12c^2, Y = 16c^3)", sp.simplify(dbeta / dlnalpha + c) == 0)
bc = sp.integrate(dbeta, c)
target = 2 * c - sp.atan(2 * sp.sqrt(3) * c) / sp.sqrt(3)
check("K1b beta(c) = 2c - arctan(2 sqrt3 c)/sqrt3", sp.simplify(sp.diff(bc - target, c)) == 0)
check("K1c the ratio only needs dY/dZ = 2c", sp.simplify(sp.diff(Y, c) / sp.diff(Z, c) - 2 * c) == 0)

DEG = np.pi / 180
beta = float(table("birefringence_measurements.csv")["DPK25-beta-baseline"]["value"]) * DEG
MP = val("bounds_inputs.csv", "MP_reduced")          # eV
d_e = val("shift_uholst_inputs.csv", "d_e_MICROSCOPE")
edm = val("shift_uholst_inputs.csv", "eEDM_JILA")    # e cm
alpha = val("cosmo_torsion_axion_inputs.csv", "alpha_em")
hbarc = val("bounds_inputs.csv", "hbarc")            # eV m
H0 = val("cosmo_torsion_axion_inputs.csv", "H0") * 1e3 / val("bounds_inputs.csv", "Mpc") * val("bounds_inputs.csv", "hbar")  # eV

# K2 (beta/c is dimensionless; M_P cancels in chi/M_P units)
for x in (1.0, 0.39, 0.1):
    info(f"K2 MICROSCOPE Delta chi = {x} M_P: c >= {np.sqrt(2) * beta / (d_e * x):.4g}")
c_mic39 = np.sqrt(2) * beta / (d_e * 0.39)
check("K2 reproduces 012: c >= 340 at 0.39 M_P", abs(c_mic39 / 340 - 1) < 0.01, f"{c_mic39:.1f}")

# K3
mu_c, s_c = pair("K-alphadot-clock")
adot = abs(mu_c) + 2 * s_c
c_clk = beta * H0 / adot
info(f"K3 H0 = {H0:.4e} eV; |alpha-dot/alpha| <= {adot:.3e} eV; c_clock(R=1) >= {c_clk:.3g}")
check("K3 analytic clock bound within a factor 2 of 012's ODE value 5.1e4", 0.5 < c_clk / 5.13e4 < 2, f"{c_clk:.3g}")

# K4
k = np.sqrt(2) / (MP / 1e9)                   # GeV^-1
e = np.sqrt(4 * np.pi * alpha)
gev_cm = hbarc * 100 / 1e9
d_unnorm = lambda cc: 0.75 * cc * k * gev_cm / e
d_phys = lambda cc: d_unnorm(cc) / np.sqrt(1 + 12 * cc ** 2)
info(f"K4 012 form d = {d_unnorm(1):.3e} c e cm -> c <~ {edm / d_unnorm(1):.0f}")
info(f"K4 normalised d_phys(c->inf) = {0.75 * k * gev_cm / e / np.sqrt(12):.3e} e cm; JILA bound {edm:.1e} e cm")
for cc in (1, 140, 340, 5e4):
    info(f"   c = {cc:g}: d_phys = {d_phys(cc):.2e} e cm  (ratio to bound {d_phys(cc) / edm:.1e})")
check("K4 with canonical normalisation the EDM estimate gives NO upper bound on c (d_phys < bound for all c)",
      0.75 * k * gev_cm / e / np.sqrt(12) < edm)
mu_phys = lambda cc: cc ** 2 * k / np.sqrt(1 + 12 * cc ** 2)     # GeV^-1 (times e-normalisation O(1))
info(f"K4b P'-even (MDM-type) piece: mu_phys ~ c k/sqrt12 = {mu_phys(5e4):.2e} GeV^-1 at c = 5e4 "
     f"(Delta a_e ~ 2 m_e mu/e ~ {2 * val("cosmo_torsion_axion_inputs.csv", "m_e") * 1e-9 * mu_phys(5e4) / e:.1e}); grows with c")
print(f"kk05_check: {NPASS[0]} PASS, {NPASS[1]} FAIL")
