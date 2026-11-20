"""Coordinator follow-ups: (a) low-z k_AF bounds from data/kaf_lowz_bounds.csv (task 024) with our own beta(z);
(b) map published torsion-free equivalents onto xi*eps for the n = 0 shape:
   - Yu-Yang-Qiu 2508.10532 (nonmetricity): S contains (alpha/M^2) l(Q) Ftilde F with l = M^2 ln(Q/M_P^2) (case II),
     Q = 6H^2; their Eq. 53: Delta chi ~ 40.21 alpha; Eq. 54: 8.23e-5 < alpha < 1.78e-4 from Planck PR4 0.30 +- 0.11 deg.
     In our conventions alpha ln Q F Ftilde = -(g/4) phi F Ftilde with g phi = -4 alpha ln Q, so
     |beta| = 2 alpha |ln Q_0 - ln Q_*| = 4 alpha ln E(z*), i.e. p_0hat = 4 alpha d ln(6H^2)/dt.
   - Li-Xia-Li-Zhang hep-ph/0611192: Delta chi = (delta/2) ln(R(z)/R(0)) (their Eq. 42), delta = -4 alpha_em c/(3 pi)
     (Eq. 17); WMAP3+B03 fit delta = -0.011 (+0.0074 -0.0079, 1 sigma; Sec. IV). So p_0hat = delta d ln R/dt.
   With R = 3 H0^2 (Om a^-3 + 4 OL) (radiation is traceless) and Q = 6 H^2 (incl. radiation).
(c) the theorist's table work/20261120-022-theory-ec06/ec06_backgrounds.csv (eps_z_form column) where it maps onto a shape.
"""
import csv
import math
import os

import numpy as np
from scipy.integrate import quad

from ec06_common import MEAS, EC6, ec6, DATA, ROOT, ZSTAR, ZRE, DEG, BETA_ACT, E, OMM, OML, H0_EV, beta_unit, reion_visibility

checks = []


def check(name, ok, info=""):
    checks.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}: {name} -- {info}")


K0 = beta_unit(ZSTAR, 0.0)        # beta per unit xi*eps, n = 0 (rad)


def beta_custom(zs, p_of_z):
    """1/2 int_0^zs p_0hat(z) dz/((1+z) H), p in units of H0."""
    return quad(lambda z: 0.5 * p_of_z(z) / ((1 + z) * E(1 / (1 + z))), 0, zs, epsrel=1e-10, limit=400)[0]


# ------------------------------------------------------------------ (b) published couplings
lnE = lambda z: math.log(E(1 / (1 + z)))
R_over = lambda z: OMM * (1 + z) ** 3 + 4 * OML          # R / (3 H0^2)
dlnE_dz = lambda z: (lnE(z + 1e-4 * (1 + z)) - lnE(z - 1e-4 * (1 + z))) / (2e-4 * (1 + z))
dlnR_dz = lambda z: 3 * OMM * (1 + z) ** 2 / R_over(z)
# p_0hat/H0 = kappa d ln X/dt = -kappa (1+z) E d ln X/dz ; kappa = 1
p_lnQ = lambda z: -(1 + z) * E(1 / (1 + z)) * 2 * dlnE_dz(z)
p_lnR = lambda z: -(1 + z) * E(1 / (1 + z)) * dlnR_dz(z)
bQ = beta_custom(ZSTAR, p_lnQ)     # per unit kappa (= 4 alpha for YYQ)
bR = beta_custom(ZSTAR, p_lnR)     # per unit delta (LXL)
check("P1 1/2 int kappa dlnQ/dt dt = -kappa ln E(z*) (exact total derivative)", abs(bQ + lnE(ZSTAR)) < 1e-5 * lnE(ZSTAR),
      f"{bQ:.5f} vs {-lnE(ZSTAR):.5f}")
check("P2 YYQ case II: |beta| = 4 alpha ln E(z*) reproduces their Eq. 53 (40.21 alpha) to 1%",
      abs(4 * lnE(ZSTAR) / ec6("YYQ25_chi_per_alpha") - 1) < 0.01, f"4 ln E(z*) = {4 * lnE(ZSTAR):.2f}")
lnRr = math.log(R_over(ZSTAR) / R_over(0))
check("P3 LXL: |beta| = (|delta|/2) ln(R*/R0) (numeric = closed form)", abs(-bR / (0.5 * lnRr) - 1) < 1e-5,
      f"(1/2) ln(R*/R0) = {0.5 * lnRr:.4f}")
alpha_act = BETA_ACT * DEG / (4 * lnE(ZSTAR))
b_dp, s_dp = float(MEAS["DP22-beta"]["value"]), float(MEAS["DP22-beta"]["err_plus"])
alpha_pl = (b_dp - s_dp) * DEG / (4 * lnE(ZSTAR)), (b_dp + s_dp) * DEG / (4 * lnE(ZSTAR))
check("P4 YYQ Eq. 54 range 8.23e-5 < alpha < 1.78e-4 reproduced from 0.30 +- 0.11 deg with our beta(alpha)",
      abs(alpha_pl[0] / ec6("YYQ25_alpha_lo") - 1) < 0.02 and abs(alpha_pl[1] / ec6("YYQ25_alpha_hi") - 1) < 0.02, f"{alpha_pl[0]:.3e} - {alpha_pl[1]:.3e}")
delta_act = BETA_ACT * DEG / (0.5 * lnRr)
xe_yyq = 4 * lnE(ZSTAR) / K0         # xi*eps equivalent (same total beta) per unit alpha
xe_lxl = 0.5 * lnRr / K0             # per unit |delta|
# reionization share of these shapes
zv, gv, _ = reion_visibility()
zz = zv[::20]
fr = {}
for nm, pf, tot in (("lnQ", p_lnQ, bQ), ("lnR", p_lnR, bR)):
    bz = np.array([beta_custom(z, pf) for z in zz])
    fr[nm] = (beta_custom(ZRE, pf) / tot, np.trapezoid(gv[::20] * bz, zz) / np.trapezoid(gv[::20], zz) / tot)
print(f"YYQ (ln 6H^2): |beta| = {4 * lnE(ZSTAR):.2f} alpha rad; ACT 0.215 deg -> alpha = {alpha_act:.3e}; "
      f"xi*eps(n=0, same total beta) = {xe_yyq:.2f} alpha -> {xe_yyq * alpha_act:.3e}; "
      f"matter-era local p_0hat = -12 alpha H; beta_reio/beta_rec = {fr['lnQ'][0]:.3f} (z_re) {fr['lnQ'][1]:.3f} (vis)")
print(f"LXL (ln R): |beta| = {0.5 * lnRr:.3f} |delta| rad; ACT -> |delta| = {delta_act:.3e} (their WMAP3+B03 fit "
      f"{ec6('LXL07_delta')}, i.e. |beta| = {abs(ec6('LXL07_delta')) * 0.5 * lnRr / DEG:.1f} deg); xi*eps(n=0) = {xe_lxl:.3f} |delta| -> {xe_lxl * delta_act:.3e}; "
      f"matter-era local p_0hat = -3 delta H; beta_reio/beta_rec = {fr['lnR'][0]:.3f} (z_re) {fr['lnR'][1]:.3f} (vis)")
check("P5 both published maps give the same xi*eps as our n = 0 fit when fed ACT (consistency of the bookkeeping)",
      abs(xe_yyq * alpha_act / (BETA_ACT * DEG / K0) - 1) < 1e-9 and abs(xe_lxl * delta_act / (BETA_ACT * DEG / K0) - 1) < 1e-9,
      f"{BETA_ACT * DEG / K0:.4e}")

# ------------------------------------------------------------------ (a) low-z k_AF bounds (task 024 file)
with open(os.path.join(DATA, "kaf_lowz_bounds.csv"), newline="", encoding="utf-8-sig") as fh:
    L = list(csv.DictReader(fh))
print("\nlow-z bound | z_eff | angle limit [deg] | xi*eps limit n=0 (ours) | scout | n=-1 | const | margin vs ACT n=0")
best = None
for r in L:
    if r["kind"] != "meas_angle_deg" or not r["z_eff"]:
        continue
    zeff = float(r["z_eff"])
    alim = abs(float(r["central"])) + 2 * float(r["err_1sigma_or_limit"])
    lim = {n: alim * DEG / beta_unit(zeff, n) for n in (0.0, -1.0, None)}
    scout = r["ec06_xi_eps_limit"]
    margin = lim[0.0] / (BETA_ACT * DEG / K0)
    print(f"{r['id']} | {zeff} | {alim:.2f} | {lim[0.0]:.3f} | {scout} | {lim[-1.0]:.3f} | {lim[None]:.3f} | {margin:.0f}")
    if scout:
        check(f"P6 [{r['id']}] our n = 0 xi*eps limit agrees with the scout's (10%)", abs(lim[0.0] / float(scout) - 1) < 0.1,
              f"{lim[0.0]:.3f} vs {scout}")
    best = min(best or 1e9, margin)
for r in L:
    if r["kind"] in ("upper_p0_GeVh", "meas_p0_GeVh", "meas_p0_eVh") and r["p0_limit_eV"]:
        p0 = float(r["p0_limit_eV"])
        print(f"{r['id']}: p0 < {p0:.2e} eV -> xi*eps < {p0 / H0_EV:.3g} (p_0hat today = xi eps H0; z range {r['z_min']}-{r['z_max']})")
check("P7 every low-z angle bound is >= 50x weaker than the ACT n = 0 requirement", best >= 50, f"tightest margin {best:.0f}x")

# ------------------------------------------------------------------ (c) theorist table (022)
path = os.path.join(ROOT, "work", "20261120-022-theory-ec06", "ec06_backgrounds.csv")
if os.path.exists(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        T = list(csv.DictReader(fh))
    print("\n022 table (eps = p_0hat/(xi H) per component):")
    with open(os.path.join(DATA, "ec06_theory_inputs.csv"), newline="", encoding="utf-8-sig") as fh:
        TI = {r["id"]: r for r in csv.DictReader(fh)}
    for r in T:
        form = r["eps_z_form"]
        try:
            e0 = float(r["eps_today"])
        except ValueError:
            e0 = None
        if "constant" in form and e0 is not None:
            xi = BETA_ACT * DEG / (abs(e0) * K0)
            print(f"  {r['theory']} [{r['component']}]: n = 0, |eps| = {abs(e0)} -> xi(ACT) = {xi:.3e}; "
                  f"022 says: {r['needed_xi_and_notes']}")
        elif "sqrt|Omega_K|" in form:
            # eps = 2 sqrt|OK| (1+z)/E  ->  p_0hat = xi * 2 sqrt|OK| H0 (1+z)
            OK = abs(float(TI["Omega_K"]["value"])) + 1.9 * float(TI["Omega_K"]["err"])  # 95% upper |Omega_K|
            pf = lambda z: 2 * math.sqrt(OK) * (1 + z)
            b = beta_custom(ZSTAR, pf)
            fz = beta_custom(ZRE, pf) / b
            print(f"  {r['theory']} [{r['component']}]: p_0hat = 2 xi sqrt|OK| H0 (1+z), |OK| = {OK:.4f} -> "
                  f"xi(ACT) = {BETA_ACT * DEG / b:.3e}, beta(z_re)/beta = {fz:.2f}; 022 says: {r['needed_xi_and_notes']}")
        else:
            print(f"  {r['theory']} [{r['component']}]: not a parameterized shape here ({form}); 022 says: {r['needed_xi_and_notes']}")
else:
    print("INFO: 022 table not present")
print(f"\n{sum(checks)} PASS, {len(checks) - sum(checks)} FAIL")
