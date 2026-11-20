"""EC-06d: cited bounds on H(z) at z* and z = 100 versus the SNY histories of task 025 (sny_background.csv, read only)."""
import csv
import math
import os

from common import ROOT, check, load, summary, val

eb = load("expansion_history_bounds.csv")
ci = load("cosmo_inputs.csv")
zstar, Om, Nsm = val(ci, "z_star"), val(ci, "Omega_m"), val(ci, "N_eff_SM")
zeq, anu = val(eb, "z_eq_P18"), val(eb, "nu_per_species")

# allowed delta H / H (95%)
r = (1 + zstar) / (1 + zeq)                        # rho_r / rho_m at z*
f_r = r / (1 + r)                                  # radiation fraction of the total at z* (Lambda negligible)
drho_Neff = val(eb, "dNeff_P18") * anu / (1 + Nsm * anu) * f_r
dH_Neff = math.sqrt(1 + drho_Neff) - 1
dH_Oe = 1 / math.sqrt(1 - val(eb, "Omega_e_P15")) - 1
dH_Oe50 = 1 / math.sqrt(1 - val(eb, "Omega_e_z50_P15")) - 1
print(f"rho_r/rho_m(z*) = {r:.3f}, radiation fraction {f_r:.3f}")
print(f"allowed dH/H at z* (95%): Delta N_eff < 0.30 -> {dH_Neff:.2e}; Omega_e < 0.0036 -> {dH_Oe:.2e}")
print(f"allowed dH/H for z > 50 (95%): Omega_e(z>50) < 0.02 -> {dH_Oe50:.2e}")
check("both z* bounds allow < 1% in H", max(dH_Neff, dH_Oe) < 0.01)

path = os.path.join(ROOT, "work", "20261120-025-theory-ec06c-costs", "sny_background.csv")
with open(path, encoding="utf-8") as fh:
    rows = list(csv.DictReader(fh))
out = []
print("\ncase   z       H/H_LCDM(025, no rad)  H/H_LCDM(with rad)  orders in H  orders in rho   bound")
for row in rows:
    if row["reached"] != "1" or float(row["z"]) < 50:
        continue
    z, H = float(row["z"]), float(row["H"])
    R0 = float(row["H_over_LCDM"])
    HL = math.sqrt(Om * (1 + z)**3 * (1 + (1 + z) / (1 + zeq)) + 1 - Om)   # SNY units: H0 = 1
    R = H / HL
    bounds = [("dNeff", dH_Neff), ("Omega_e", dH_Oe)] if abs(z - zstar) < 1 else [("Omega_e(z>50)", dH_Oe50)]
    for name, dh in bounds:
        oH = math.log10((R - 1) / dh)
        orho = math.log10((R**2 - 1) / ((1 + dh)**2 - 1))
        print(f"{row['case']:5s} {z:8.2f}  {R0:10.3g}             {R:10.3g}         {oH:6.2f}       {orho:6.2f}      {name}")
        out.append((row["case"], z, R, name, oH, orho))
star = [o for o in out if abs(o[1] - zstar) < 1]
low = [o for o in out if o[1] < 200]
check("radiation correction: 025's H/H_LCDM(z*) is high by sqrt(1 + rho_r/rho_m) = 1.15 (025 comment says ~1.5)",
      abs(math.sqrt(1 + r) - 1.149) < 0.01, f"({math.sqrt(1 + r):.3f})")
check("every SNY case excluded at z* by > 3.5 orders in H (both bounds)", min(o[4] for o in star) > 3.5,
      f"(min {min(o[4] for o in star):.2f})")
check("every SNY case excluded at z = 100 by > 2.5 orders in H", min(o[4] for o in low) > 2.5,
      f"(min {min(o[4] for o in low):.2f})")
print(f"\nSummary: H(z*)/H_LCDM = {min(o[2] for o in star):.1f}-{max(o[2] for o in star):.1f} (radiation included); "
      f"shortfall {min(o[4] for o in star):.1f}-{max(o[4] for o in star):.1f} orders in H, "
      f"{min(o[5] for o in star):.1f}-{max(o[5] for o in star):.1f} orders in energy density (95% bounds).")
summary()
