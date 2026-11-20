"""Cross-check of EC-02's xi with the neutrino-current fit of Yin-Komatsu 2610.02667 (M/sqrt(b_J xi_nu) = 44 TeV).
Minimal EC: S^mu = -(kappa/2) J5^mu; CFJ -1/2 xi S_mu A_nu Ftilde -> |b_J|/M^2 = xi kappa/2. Order of magnitude only:
axial vs vector current, neutrino vs antineutrino, and the Planck PR4 beta used by Yin-Komatsu are not tracked."""
import csv
import math
import os

from common import ROOT, check, load, summary, val

ci = load("cosmo_inputs.csv")
bi = load("bounds_inputs.csv")
Lam = val(ci, "Reig_Lambda_star_TeV_per_sqrt_xi") * 1e12        # eV
xinu = val(ci, "xi_nu") + 2 * val(ci, "xi_nu", "err")            # 95% upper
MP = val(bi, "MP_reduced")
kappa = 1 / MP**2
xi = 2 / (kappa * xinu * Lam**2)
print(f"xi needed (EC spin torsion, nu asymmetry at the 95% bound xi_nu = {xinu}): {xi:.1e}")
with open(os.path.join(ROOT, "src", "claims_expected.csv"), encoding="utf-8") as f:
    exp = [r for r in csv.DictReader(f) if r["claim"] == "EC-02" and r["key"] == "xi_needed"][0]
print(f"CLAIMS EC-02 xi_needed = {exp['expected']}; ratio {xi / float(exp['expected']):.1f}")
check("Yin-Komatsu fit implies the same 28-29 orders as EC-02 (within a factor 10)",
      0.1 < xi / float(exp["expected"]) < 10, f"(log10 xi = {math.log10(xi):.2f})")
summary()
