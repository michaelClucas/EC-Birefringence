"""Goals 1-2: beta_total and beta_reio for EC-06 background shapes, checks, accumulation plot, ALP contrast.

Writes beta_shapes.csv, alp_contrast.csv, beta_accumulation.png. Prints PASS/FAIL lines.
"""
import csv
import os

import numpy as np

from ec06_common import (HERE, SHAPES, ZSTAR, ZRE, DEG, BETA_ACT, H0_EV, E, c008, beta_unit, beta_unit_analytic,
                         beta_unit_conformal, reion_visibility, PL18)

checks = []


def check(name, ok, info=""):
    checks.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}: {name} -- {info}")


# ------------------------------------------------------------------ beta per unit xi*eps, per shape
z_vis, g_vis, tau_cum = reion_visibility()
tau_tot = tau_cum[-1]
tau_pl = float(PL18["tau"]["value"])
check("C0 tanh reionization optical depth vs Planck tau (z_re from data)", abs(tau_tot / tau_pl - 1) < 0.05,
      f"tau = {tau_tot:.4f} vs {tau_pl} (He II at 3.5 included)")

rows = []
for name, n in SHAPES:
    b_tot = beta_unit(ZSTAR, n)
    b_re = beta_unit(ZRE, n)
    b_10 = beta_unit(10.0, n)
    bz = np.array([beta_unit(z, n) for z in z_vis[::20]])
    b_vis = np.trapezoid(g_vis[::20] * bz, z_vis[::20]) / np.trapezoid(g_vis[::20], z_vis[::20])
    b_conf = beta_unit_conformal(ZSTAR, n)
    an = beta_unit_analytic(ZSTAR, n)
    if an is not None:
        check(f"C1 [{name}] numeric = analytic beta/(xi eps)", abs(b_tot / an - 1) < 1e-7,
              f"{b_tot:.8f} vs {an:.8f}")
    check(f"C2 [{name}] cosmic-time route = conformal-time route (p_tau = a p_0hat)", abs(b_conf / b_tot - 1) < 1e-4,
          f"{b_tot:.6f} vs {b_conf:.6f}")
    xieps = BETA_ACT * DEG / b_tot
    rows.append(dict(shape=name, n="" if n is None else n, beta_per_xieps_rad=b_tot,
                     beta_per_xieps_deg=b_tot / DEG, xieps_for_ACT=xieps,
                     frac_reio_zre=b_re / b_tot, frac_reio_vis=b_vis / b_tot, frac_z_lt_10=b_10 / b_tot,
                     beta_reio_zre_deg=BETA_ACT * b_re / b_tot, beta_reio_vis_deg=BETA_ACT * b_vis / b_tot,
                     dbeta_rec_minus_reio_deg=BETA_ACT * (1 - b_vis / b_tot),
                     p0_today_eV=xieps * H0_EV))  # p_0hat(z=0) = xi eps H0 for every shape

# skeptic 020 estimate and analytic fraction
r0 = rows[0]
check("C3 skeptic 020: S_0hat = eps H needs xi eps = 2 beta / ln(1+z*) = 1.07e-3", abs(r0["xieps_for_ACT"] - 1.07e-3) < 0.01e-3,
      f"{r0['xieps_for_ACT']:.4e}")
check("C4 [H] reionization fraction = ln(1+z_re)/ln(1+z*)",
      abs(r0["frac_reio_zre"] - np.log1p(ZRE) / np.log1p(ZSTAR)) < 1e-7, f"{r0['frac_reio_zre']:.5f}")
# constant p: beta = 1/2 p t_L; compare with literature/bounds.md |p_0| = 1.1e-35 eV (t_L = 13.80 Gyr, beta = 0.215)
tL = c008.T_TODAY - c008.t_of_a(1 / (1 + ZSTAR))
rc = rows[-1]
check("C5 [const] beta/(xi eps) = 1/2 H0 t_L (t_L from the 008 cosmic-time table)",
      abs(rc["beta_per_xieps_rad"] / (0.5 * tL) - 1) < 1e-4, f"H0 t_L = {tL:.5f}")
tL_gyr = tL / H0_EV * c008.HBAR / 3.15576e16
check("C6 [const] p_0 for ACT matches literature/bounds.md 1.1e-35 eV (t_L = 13.80 Gyr)",
      abs(rc["p0_today_eV"] / 1.1e-35 - 1) < 0.05, f"p_0 = {rc['p0_today_eV']:.3e} eV, t_L = {tL_gyr:.3f} Gyr")
check("C7 fraction at reionization increases as n decreases (more weight at low z)",
      rows[1]["frac_reio_vis"] > rows[2]["frac_reio_vis"] > rows[0]["frac_reio_vis"] > rows[3]["frac_reio_vis"] > rows[4]["frac_reio_vis"],
      ", ".join(f"{r['shape']}: {r['frac_reio_vis']:.3f}" for r in rows))

with open(os.path.join(HERE, "beta_shapes.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})

print("\nshape | beta/(xi eps) [deg] | xi eps for 0.215 deg | beta_reio/beta_rec (z_re pt, vis-wtd, z<10) | beta_rec - beta_reio [deg]")
for r in rows:
    print(f"{r['shape']:>13} | {r['beta_per_xieps_deg']:9.3f} | {r['xieps_for_ACT']:.3e} | "
          f"{r['frac_reio_zre']:.3f} {r['frac_reio_vis']:.3f} {r['frac_z_lt_10']:.3f} | {r['dbeta_rec_minus_reio_deg']:.3f}")

# ------------------------------------------------------------------ ALP contrast (task-008 field ODE)
from field_evolution import evolve, phi_of_a  # noqa: E402

alp_rows = []
zgrid = np.unique(np.concatenate([np.geomspace(1e-3, ZSTAR, 300), [ZRE, ZSTAR]]))
alp_curves = {}
for m in [1e-33, 1e-32, 2e-32, 3e-32, 5e-32, 1e-31, 1e-30, 1e-29]:
    out = evolve(m / H0_EV)
    phi0 = phi_of_a(out, np.array([1.0]))[0]
    ph = phi_of_a(out, 1 / (1 + zgrid))
    acc = phi0 - ph                      # beta(z_s) propto -(g/2)(phi_0 - phi(z_s))
    acc_rec = acc[-1]
    ph_vis = phi_of_a(out, 1 / (1 + z_vis[::20]))
    b_vis = np.trapezoid(g_vis[::20] * (phi0 - ph_vis), z_vis[::20]) / np.trapezoid(g_vis[::20], z_vis[::20])
    frac = b_vis / acc_rec
    alp_curves[m] = acc / acc_rec
    alp_rows.append(dict(m_eV=m, frac_reio_vis=frac, beta_reio_vis_deg=BETA_ACT * frac))
    print(f"ALP m = {m:.0e} eV: beta_reio/beta_rec (vis-weighted) = {frac:.4f}")
check("C8 ALP with m = 1e-29 eV (frozen through recombination, decayed by reionization) has beta_reio/beta_rec < 0.01",
      abs(alp_rows[-1]["frac_reio_vis"]) < 0.01, f"{alp_rows[-1]['frac_reio_vis']:.2e}")
check("C9 slow-roll ALP m = 1e-33 eV: beta_reio/beta_rec > 0.9 (rolls at z < few)", alp_rows[0]["frac_reio_vis"] > 0.9,
      f"{alp_rows[0]['frac_reio_vis']:.3f}")
with open(os.path.join(HERE, "alp_contrast.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(alp_rows[0].keys()))
    w.writeheader()
    for r in alp_rows:
        w.writerow({k: f"{v:.6g}" for k, v in r.items()})

# ------------------------------------------------------------------ plot
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, ax = plt.subplots(figsize=(7, 4.3), dpi=100)
cols = ["#1f4e79", "#2e86c1", "#5dade2", "#e67e22", "#c0392b", "#117a65"]
for (name, n), col in zip(SHAPES, cols):
    b = np.array([beta_unit(z, n) for z in zgrid])
    lab = f"CFJ p = xi eps {name}" if n is not None else "CFJ p = const (k_AF-like)"
    ax.plot(1 + zgrid, b / b[-1], color=col, lw=2, label=lab)
for m, ls in zip([1e-32, 1e-30, 1e-29], [":", "--", "-."]):
    ax.plot(1 + zgrid, alp_curves[m], color="0.4", ls=ls, lw=1.4, label=f"ALP m = {m:.0e} eV")
gn = g_vis / g_vis.max()
ax.fill_between(1 + z_vis, 0, 0.25 * gn, color="0.85", label="reionization visibility (arb.)")
ax.axvline(1 + ZRE, color="0.6", lw=0.8)
ax.set_xscale("log")
ax.set_xlim(1, 1 + ZSTAR * 1.05)
ax.set_ylim(-0.05, 1.05)
ax.set_xlabel("1 + z_s (source redshift)")
ax.set_ylabel("beta(z_s) / beta(z*)")
ax.set_title("Rotation accumulated between z_s and today")
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
png = os.path.join(HERE, "beta_accumulation.png")
fig.savefig(png, dpi=100)
sz = os.path.getsize(png) / 1024
check("C10 plot size < 300 KB", sz < 300, f"{sz:.0f} KB")

print(f"\n{sum(checks)} PASS, {len(checks) - sum(checks)} FAIL")
