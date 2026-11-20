"""Item 2: beta/(g phi_i) vs m for a homogeneous pseudoscalar with V = m^2 phi^2 / 2.

Writes scan.csv and beta_vs_m.png. Runtime ~1.5 min on one core.

Columns (phi_i = 1):
  phibar_lss   visibility-averaged phi over recombination
  env_lss      oscillation envelope at the visibility peak (= |phi| if not oscillating)
  washout      |phibar_lss| / env_lss
  phi0         homogeneous phi today (NaN if oscillating today: then only env0 is meaningful)
  env0         envelope today
  beta_static  beta/(g phi_i) = -(phi0 - phibar_lss)/2 (only if not oscillating today)
  beta_lss     LSS term alone = +phibar_lss/2 (what survives a time average of the observer term)
  beta_maxabs  (env0 + |phibar_lss|)/2 if oscillating today: largest |beta|/(g phi_i) over the phase today;
               = |beta_static| otherwise
  rho0         rho_phi today / phi_i^2 [eV^2]
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import HERE
from field_evolution import scan

M = np.concatenate([np.logspace(-35, -25, 41), np.logspace(-24.5, -23, 4)])


def main():
    rows = scan(M)
    for r in rows:
        r["washout"] = abs(r["phibar_lss"]) / r["env_lss"]
        r["beta_static"] = -(r["phi0"] - r["phibar_lss"]) / 2 if not r["switched"] else np.nan
        r["beta_lss"] = r["phibar_lss"] / 2
        r["beta_maxabs"] = (r["env0"] + abs(r["phibar_lss"])) / 2 if r["switched"] else abs(r["beta_static"])
    keys = ["m", "phibar_lss", "env_lss", "washout", "phi0", "env0", "beta_static", "beta_lss",
            "beta_maxabs", "rho0", "switched", "a_s"]
    with open(os.path.join(HERE, "scan.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([r[k] if not isinstance(r[k], float) else f"{r[k]:.6e}" for k in keys])
    print(f"{'m [eV]':>10} {'phibar_LSS':>12} {'washout':>10} {'phi0':>12} {'env0':>10} {'beta/(g phi_i)':>15}")
    for r in rows:
        b = r["beta_static"] if not r["switched"] else r["beta_lss"]
        print(f"{r['m']:10.3e} {r['phibar_lss']:12.4e} {r['washout']:10.3e} {r['phi0']:12.4e} "
              f"{r['env0']:10.3e} {b:15.4e}")

    m = np.array([r["m"] for r in rows])
    floor = 1e-6   # numerical noise floor of the LSS average (see README)
    bl = np.array([max(abs(r["beta_lss"]), floor / 2) for r in rows])
    bs = np.array([abs(r["beta_static"]) for r in rows])
    bm = np.array([r["beta_maxabs"] for r in rows])
    wo = np.array([max(r["washout"], floor) for r in rows])
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=110)
    ax.loglog(m, bs, "k-", lw=2, label=r"$|\beta|/(g\phi_i)$, field static today")
    ax.loglog(m, bm, "C1--", label=r"max over phase today $(\mathrm{env}_0+|\bar\phi_{LSS}|)/2\phi_i$")
    ax.loglog(m, bl, "C0-", label=r"LSS term $|\bar\phi_{LSS}|/2\phi_i$ (time-averaged $\beta$)")
    ax.loglog(m, wo, "C2:", label=r"washout factor $|\bar\phi_{LSS}|/\mathrm{env}_{LSS}$")
    ax.axhline(floor, color="0.6", lw=0.8)
    ax.text(2e-35, floor * 1.6, "numerical floor", color="0.4", fontsize=7)
    ax.set_xlabel("m [eV]")
    ax.set_ylabel(r"$\beta/(g\,\phi_i)$  [rad per unit $g\phi_i$]")
    ax.set_ylim(1e-7, 2)
    ax.set_title(r"$\phi''+2\mathcal{H}\phi'+a^2m^2\phi=0$ on Planck 2018 $\Lambda$CDM", fontsize=9)
    ax.legend(fontsize=7, loc="lower left")
    ax.grid(alpha=0.3, which="both", lw=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "beta_vs_m.png"))
    print("wrote scan.csv, beta_vs_m.png")


if __name__ == "__main__":
    main()
