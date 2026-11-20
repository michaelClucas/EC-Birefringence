"""Goal 5: low-z astrophysical CFJ (k_AF) constraints on the EC-06 background.

Today p_0hat = xi eps H0 for every shape; at z <~ few it is a slowly varying isotropic timelike CFJ vector. SME matching
(literature/bounds.md Sec. 2): p_mu = -2 (k_AF)_mu, so |k_AF| = |p|/2. Compared with:
  - bounds.csv P-kAF-CFJ90 (radio galaxies, CFJ 1990), P-kAF0-pulsar (2503.20626), P-kAF-Caloni (CMB, constant k_AF);
  - radio-galaxy UV polarization at z > 2, theta = -0.8 +- 2.2 deg (1003.4823, data/ec06_cosmo_inputs.csv),
    using the rotation the background itself accumulates from z = 2 (direct, no k_AF conversion).
"""
import csv
import os

from ec06_common import DATA, EC6, SHAPES, ZSTAR, DEG, BETA_ACT, H0_EV, beta_unit, ec6

checks = []


def check(name, ok, info=""):
    checks.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}: {name} -- {info}")


with open(os.path.join(DATA, "bounds.csv"), newline="", encoding="utf-8-sig") as fh:
    B = {r["id"]: r for r in csv.DictReader(fh)}


def upper_eV(rid):
    return float(B[rid]["value_eV"].replace("<", "").strip().split()[0])


kaf_rg, kaf_psr, kaf_cmb = upper_eV("P-kAF-CFJ90"), upper_eV("P-kAF0-pulsar"), upper_eV("P-kAF-Caloni")
check("K0 bounds.csv GeV->eV conversion consistent (CFJ90: 6e-43 GeV = 6e-34 eV)", abs(kaf_rg / 6e-34 - 1) < 1e-9, f"{kaf_rg:.1e}")

th = ec6("RG_UV_theta_deg")
th_err = float(EC6["RG_UV_theta_deg"]["err"])
z_rg = ec6("RG_UV_zmin")
print(f"bounds: radio-galaxy |k_AF| < {kaf_rg:.1e} eV, pulsar k0 < {kaf_psr:.1e} eV, CMB (constant) |k_AF| < {kaf_cmb:.1e} eV; "
      f"RG UV theta = {th} +- {th_err} deg (z > {z_rg})")
worst = 0.0
for name, n in SHAPES:
    xieps = BETA_ACT * DEG / beta_unit(ZSTAR, n)
    p0 = xieps * H0_EV
    kaf0 = p0 / 2
    b_rg = BETA_ACT * beta_unit(z_rg, n) / beta_unit(ZSTAR, n)
    worst = max(worst, b_rg)
    print(f"{name:>13}: xi eps = {xieps:.3e}, p_0hat(0) = {p0:.2e} eV, |k_AF|(0) = {kaf0:.2e} eV "
          f"(RG/{kaf0 / kaf_rg:.0e}, pulsar/{kaf0 / kaf_psr:.0e}); beta(z_s = {z_rg}) = {b_rg:.4f} deg "
          f"= {b_rg / th_err:.3f} x the RG UV error")
    check(f"K1 [{name}] |k_AF| today below the radio-galaxy and pulsar bounds", kaf0 < kaf_rg and kaf0 < kaf_psr, f"{kaf0:.2e} eV")
check("K2 rotation accumulated from z = 2 is < 0.1 x the radio-galaxy UV error for every shape (low-z probes cannot bind)",
      worst < 0.1 * th_err, f"max {worst:.3f} deg vs {th_err} deg")
# the CMB 'constant k_AF' bound is the same CMB beta re-expressed; for the const shape it must be of the same order
xc = BETA_ACT * DEG / beta_unit(ZSTAR, None) * H0_EV / 2
check("K3 [const] k_AF that fits ACT is within a factor 2 of the CMB-derived |k_AF| scale (same observable)",
      0.5 < xc / kaf_cmb < 2, f"{xc:.2e} vs {kaf_cmb:.2e} eV")
print(f"\n{sum(checks)} PASS, {len(checks) - sum(checks)} FAIL")
