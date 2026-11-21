"""036 step C: endpoints, fifth-force and photon-coupling bounds for the heavy SNY 0+ branch (EC-06d).

Uses check_scalar_mode.py results (units M_P = 1 restored explicitly here):
  m^2 = 2 mu M_P^2/(a1 b), Yukawa strength alpha = mu/(3 a1) >= 1/3, gamma - 1 = -2 mu/(4 mu - 3) (|.| >= 1/2),
  Y_c = 3 (a1/M_P^2) m^2, heavy limit chi = kappa * (Yukawa-smoothed rho over the Compton length lam = 1/m),
  beta = -(3 xi/4) ln(Y_obs/Y_src) ~ -(3 xi/4) kappa (rho_obs - rho_src)/Y_c,
  |beta| = sqrt(2 alpha) g |rho_obs - rho_src| /(4 M_P m^2), g = photon coupling of the canonical 0+ field.
Calibration endpoints: ACT DR6 fixes the instrument angle with an optics-model prior, so its beta contains the
observer endpoint; Planck/WMAP Minami-Komatsu (EK22) fix it with Galactic-foreground EB, so their beta is
rho_dust - rho_* and the observer endpoint cancels.
"""
import math

from common import check, load, summary, val

cos = load("cosmo_inputs.csv")
pi_ = load("pgt0plus_inputs.csv")
meas = load("birefringence_measurements.csv")
bnd = load("bounds.csv")
kkb = load("kk_bounds_inputs.csv")

hbarc = val(cos, "hbar_c_eV_m")
hbar = val(cos, "hbar_eV_s")
H0 = val(cos, "H0") / (val(pi_, "Mpc_m") / 1e3) * hbar  # eV
Om, zs = val(cos, "Omega_m"), val(cos, "z_star")
MP = val(pi_, "M_P_red_eV")
kg = val(pi_, "eV_per_kg")
mp = val(pi_, "m_p_eV")
KGM3 = kg * hbarc**3  # 1 kg/m^3 in eV^4
rho_crit = 3 * H0**2 * MP**2
rho_star = Om * rho_crit * (1 + zs)**3
rho_mean = Om * rho_crit
cm3 = 1e6 * hbarc**3  # 1 cm^-3 in eV^3
rho_dust = val(pi_, "n_ism_diffuse") * cm3 * mp
rho_sw = val(pi_, "n_sw_1AU") * cm3 * mp
rho_air = val(pi_, "rho_air_5200m") * KGM3
rho_film, t_film = val(pi_, "rho_film") * KGM3, val(pi_, "t_film_assumed")
Msun, ME, RE, AU, dL2, Msc = (val(pi_, k) * f for k, f in
                              [("M_sun", kg), ("M_earth", kg), ("R_earth", 1), ("AU", 1), ("d_L2", 1), ("M_planck_sc", kg)])
beta_act = math.radians(val(meas, "DPK25-beta-baseline"))
beta_mk = math.radians(val(meas, "EK22-beta"))
sig_act = math.radians(float(meas["DPK25-beta-baseline"]["err_plus"]))
sig_mk = math.radians(float(meas["EK22-beta"]["err_plus"]))
g_hb = float(bnd["A-HB-Ayala"]["value_eV"].replace("<", "").strip())
gm1, gm1e = float(kkb["gamma_minus_1_Cassini"]["value"]), float(kkb["gamma_minus_1_Cassini"]["err"])
lam_ew = val(pi_, "EW20-lambda-alpha1")
lam_max = val(pi_, "lambda_max_alpha13_assumed")
print(f"H0 = {H0:.4e} eV; rho_* = {rho_star:.3e} eV^4 = {rho_star / KGM3:.2e} kg/m^3; "
      f"rho_mean,0 = {rho_mean / KGM3:.2e} kg/m^3; rho_dust = {rho_dust / rho_star:.1e} rho_*")


def yuk_point(M, r, lam):
    """Yukawa-smoothed density at distance r from a compact mass M (r, lam in m; M in eV): M e^{-r/lam}/(4 pi r lam^2)."""
    return M * math.exp(-r / lam) / (4 * math.pi * (r / hbarc) * (lam / hbarc)**2)


def rho_obs(lam, where):
    """Lower estimate of the smoothed density at the polarization analyzer (eV^4)."""
    if where == "cosmic":  # idealized observer at the cosmic mean density today
        return rho_mean
    if lam < 1.0:  # instrument scale: the analyzing film itself
        return rho_film * min(1.0, t_film / (2 * lam))
    r = 0.0
    if where == "ground":
        r += 0.5 * rho_air if lam < RE else yuk_point(ME, RE, lam)
        r += yuk_point(Msun, AU, lam)
    elif where == "L2":
        r += yuk_point(Msc, 1.0, lam) + yuk_point(ME, dL2, lam) + yuk_point(Msun, AU, lam) + rho_sw
    return r


# --- C1: fifth force. lam >> AU: Cassini; alpha >= 1/3 for lam > lam_max excluded by Eot-Wash-type data
lim_g = abs(gm1) + 2 * gm1e
print(f"C1 Cassini |gamma-1| < {lim_g:.1e} (2 sigma) vs prediction >= 0.5: {math.log10(0.5 / lim_g):.1f} orders")
check("C1 light 0+ (lam >~ AU): |gamma - 1| >= 1/2 excluded by Cassini", 0.5 > 100 * lim_g)
m_min = hbarc / lam_max
print(f"   short-range: alpha = 1 excluded for lam > {lam_ew * 1e6:.1f} um (EW20); alpha = 1/3 taken as lam < "
      f"{lam_max * 1e6:.0f} um -> m_0+ >= {m_min:.2e} eV")

# 025's GR-compatible solution (delta(z*) = 1%, mu = 2, a1' = 1): its mass
delta, mu, a1p = 0.01, 2.0, 1.0
s_star = 4 * delta / (9 * mu)
Yc_025 = (rho_star / MP**2) / s_star
m_025 = math.sqrt(Yc_025 / (3 * a1p))
lam_025 = hbarc / m_025
print(f"   025 GR-compatible example (delta = 1%, mu = 2): m = {m_025:.1e} eV, lam = {lam_025 / 3.0857e19:.2f} kpc")
check("C1b 025's xi ~ 1 heavy solution has lam >> AU, so it is in the Cassini-excluded range", lam_025 > 1e3 * AU)

# --- C2: endpoints vs lam
print("C2 lam [m] | rho_obs ground/rho_* | L2/rho_* | cosmic/rho_* | beta_ACT/beta_MK predicted | "
      "beta/xi ground, L2, cosmic (a1' = 1)")
lam_cross = None
for lam in [1e-9, 1e-6, lam_ew, lam_max, 1e-2, 1.0, 1e3, RE, 1e9, AU, 1e14, 1e16, 3.0857e16, 3.0857e17,
            3.0857e18, 3.0857e19, 3.0857e20]:
    m = hbarc / lam
    Yc = 3 * a1p * m**2
    rg, rl, rc = rho_obs(lam, "ground"), rho_obs(lam, "L2"), rho_obs(lam, "cosmic")
    ratio = (rg - rho_star) / (rho_dust - rho_star)
    bx = [-(0.75) * (r - rho_star) / MP**2 / Yc for r in (rg, rl, rc)]
    print(f"   {lam:9.2e} | {rg / rho_star:9.2e} | {rl / rho_star:9.2e} | {rc / rho_star:8.1e} | {ratio:+10.2e} | "
          + ", ".join(f"{v:+.2e}" for v in bx))
grid = [10**(k / 20) for k in range(-180, 440)]
for lam in grid:
    if rho_obs(lam, "ground") < rho_star:
        lam_cross = lam
        break
print(f"   ground observer density drops below rho_* only for lam > {lam_cross:.2e} m = {lam_cross / 3.0857e16:.1f} pc")
r_meas = beta_act / beta_mk
r_err = r_meas * math.hypot(sig_act / beta_act, sig_mk / beta_mk)
print(f"   measured beta_ACT/beta_EK22 = {r_meas:.2f} +- {r_err:.2f}")
check("C2a for every lam < lam_cross the predicted beta_ACT/beta_MK is negative (measured > 0 at > 2 sigma)",
      all((rho_obs(l, 'ground') - rho_star) / (rho_dust - rho_star) < 0 for l in grid if l < lam_cross)
      and r_meas - 2 * r_err > 0)
check("C2b endpoint kill range (lam < lam_cross) overlaps Cassini range (lam >~ AU)", lam_cross > 10 * AU)
rmin = min(abs((rho_obs(l, 'ground') - rho_star) / (rho_dust - rho_star)) for l in grid if l <= lam_max)
print(f"   for lam <= {lam_max * 1e6:.0f} um: |beta_ACT/beta_MK| >= {rmin:.1e} (measured {r_meas:.2f})")
check("C2c in the fifth-force-allowed range the predicted |ratio| exceeds 1e10", rmin > 1e10)

# --- C3: size of beta at the HB photon-coupling bound, lam <= lam_max, alpha = 1/3 (a1 -> infinity)
alpha = 1 / 3
for lam in [lam_max, lam_ew, 1e-6]:
    m = hbarc / lam
    b_obs = math.sqrt(2 * alpha) * g_hb * rho_obs(lam, "ground") / (4 * MP * m**2)
    b_uni = math.sqrt(2 * alpha) * g_hb * (rho_star - rho_dust) / (4 * MP * m**2)
    print(f"C3 lam = {lam:.1e} m: |beta| max (ACT observer piece) = {b_obs:.1e} rad, short by "
          f"{math.log10(beta_act / b_obs):.1f} orders; universal (MK) piece {b_uni:.1e} rad, short by "
          f"{math.log10(beta_mk / b_uni):.1f} orders")
m = hbarc / lam_max
b_obs = math.sqrt(2 * alpha) * g_hb * rho_obs(lam_max, "ground") / (4 * MP * m**2)
b_uni = math.sqrt(2 * alpha) * g_hb * (rho_star - rho_dust) / (4 * MP * m**2)
check("C3a ACT-type observer piece at the HB bound is > 15 orders short", math.log10(beta_act / b_obs) > 15)
check("C3b universal piece at the HB bound is > 35 orders short", math.log10(beta_mk / b_uni) > 35)
# xi needed (alpha = 1/3, mu/a1 -> 1): g = xi sqrt(3/(2 mu a1))/M_P; universal piece xi = 4 beta Y_c/(3 kappa drho)
for a1p_ in (1.0, 1e3):
    Yc = 3 * a1p_ * m**2
    xi_u = 4 * beta_mk * Yc / (3 * (rho_star - rho_dust) / MP**2)
    xi_a = 4 * beta_act * Yc / (3 * rho_obs(lam_max, "ground") / MP**2)
    xi_hb = g_hb * MP * math.sqrt(2 * (a1p_ + 1) * a1p_ / 3)
    print(f"   a1' = {a1p_:g}: xi needed universal {xi_u:.1e}, ACT-observer {xi_a:.1e}; HB allows xi < {xi_hb:.1e}")
# larger alpha (a1' -> 0): beta ~ sqrt(alpha) lam^2; alpha needed to close the gap
a_strong = 8 * math.pi * MP**2 / mp**2
for lam in (lam_max, 1e-6, 1e-9):
    mm = hbarc / lam
    a_req = 0.5 * (beta_act * 4 * MP * mm**2 / (g_hb * rho_obs(lam, "ground")))**2
    print(f"   alpha needed at lam = {lam:.0e} m for ACT beta at the HB bound: {a_req:.1e} "
          f"(nucleon strong/gravity ~ {a_strong:.1e})")
check("C3c alpha needed exceeds strong-interaction strength at all three lam",
      all(0.5 * (beta_act * 4 * MP * (hbarc / l)**2 / (g_hb * rho_obs(l, 'ground')))**2 > a_strong
          for l in (lam_max, 1e-6, 1e-9)))
summary()
