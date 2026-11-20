# 20261120-005-bounds

## Question
Which independent laboratory, astrophysical and cosmological bounds must any candidate mechanism (4D EC or EC-KK) respect, in our units?

## Method
- Opened the Kostelecky-Russell Data Tables (arXiv:0801.0287v19, 5 Feb 2026) and extracted Tables D6, D10, D15, D19 and Summary Table S2 with `pdftotext`. Opened Kostelecky-Russell-Tasson (0712.4393) in full for the torsion-to-b mapping, Kostelecky 2004 (hep-th/0312310) App. A for conventions, Kostelecky-Mewes 2009 (0905.0031) for the d = 3 rotation formula, Kranas et al. (1809.10064) and Deutschmann et al. (1702.00410) in full, 2512.21888 v2 in HTML, Planck 2018 I via ar5iv for the dipole. Other numbers come from arXiv abstracts; `data/bounds.csv` says which.
- Every published value is in `data/bounds.csv` with units as published and converted to eV. Non-bound inputs (unit constants, CMB dipole, H0, Omega_m) are in `data/bounds_inputs.csv`. ACT beta is read from `data/birefringence_measurements.csv`.
- `convert_bounds.py` reads those files, rechecks the typed eV conversions (33 rows, no mismatch), derives torsion bounds and the ACT-equivalent CFJ coefficient, and writes `derived_bounds.csv`.

Rerun from the repo root:
```
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-005-bounds/convert_bounds.py
```
(on Windows pass the interpreter to pin.py as `.venv\Scripts\python.exe`).

## Result
Mapping: minimal coupling gives b_mu = -(3/4) S_mu (KRT Eqs. 2-3, xi_4 = 3/4); KRT's A_mu matches our S_mu (same signature, eps_{0123} = +1, same 1/6). The trace T_mu does not couple minimally.

| Quantity | Bound | Source |
|---|---|---|
| |S_X|, |S_Y| (Sun frame) | < 8.9e-25 eV | neutron b~_perp < 6.7e-34 GeV (68%), Allmendinger 2014 |
| |S_0| (CMB frame, FRW torsion), via Sun-CMB boost v_perp = 1.225e-3 | **< 7.3e-22 eV** | same, / v_perp |
| |S_0| from b~_T directly | < ~1e-18 eV | Data Tables S2 (order of magnitude) |
| ACT beta as constant CFJ: |(k_AF)_T| | (5.7 +- 2.0)e-36 eV | beta = |k_AF| t_L, t_L = 6.61e32 eV^-1 |
| same as |p_0| in our L = -(1/2) p A Ftilde | 1.1e-35 eV | p = -2 k_AF |
| CMB |k_AF| upper (Caloni 2023) | < 7.4e-36 eV | Data Tables D19 |
| g_agg, m < 1e-12 eV | < 8e-13 GeV^-1 = 8e-22 eV^-1 (99.7%) | Chandra NGC 1275 |
| g_agg, m < 3e-7 eV | < 5.4e-12 GeV^-1 (95%) | magnetic white dwarfs |
| g_agg, m < 0.02 eV | < 5.8e-11 GeV^-1 (95%) | CAST 2024 |
| ULA DM washout | g <~ 9.6e-16 GeV^-1 (m/1e-24 eV), ~2x below what ACT needs | 2512.21888 Eqs. 6, 12 |
| ULA = all DM | m > 2e-20 eV (95%) | Lyman-alpha, 2007.12705 |
| alpha_dot/alpha today | (1.0 +- 1.1)e-18 /yr = (1.5 +- 1.6)e-8 H0 | Yb+ clocks, 2010.06620 |
| Delta alpha/alpha | Oklo < 1.1e-8; z = 1.15: (1.3 +- 1.4)e-6; recombination (3.6 +- 3.7)e-3 | 1503.06011, 2112.05819, 1406.7482 |
| light radion, gravitational strength | m > 5.1e-3 eV | Eot-Wash ISL, lambda < 38.6 um |
| WEP | eta = (-1.5 +- 2.7)e-15 | MICROSCOPE |
| R if SM in bulk | 1/R > 1.4 TeV | LHC mUED recast |
| FRW trace torsion at BBN | |T_0| <~ 0.12 H ~ 5e-17 eV | Kranas 2019 Eq. 54 (trace only, not axial) |

Comparison, as an orientation only: the ACT beta needs |p_0| ~ 1e-35 eV today if p is constant, while the lab bound on a cosmological axial S_0 coupled minimally to neutrons is 7e-22 eV. A photon coupling p_mu = c S_mu with constant S_0 would need only c >~ 1e-14, so lab spin bounds do not exclude a CFJ-type torsion route by themselves. The tension is with theory: in minimal EC, S_0 today is set by the cosmic axial spin density, and the photon coupling c must be named (KRT: dimension-5 if gauge invariant).

## Verdict
Bound table complete for the six categories, with conversions checked by script. Strongest numbers for the program: |S_0| < 7.3e-22 eV (lab, via boost), g_agg < 8e-13 GeV^-1 for ultralight fields, Oklo |Delta alpha/alpha| < 1.1e-8, and the ALP-DM washout limit, which already disfavours the single-field ALP-DM interpretation of beta.

## Open issues
1. beta -> k_AF: the Data Tables pairing (beta = 0.35 deg <-> k^(3)_(V)00 = 7.32e-45 GeV, refs 2011.11254 + 2311.16368) is 4.47x smaller than KM 2009 Eq. (158) applied with t_L. Needs the theorist's own beta(p) formula; also the 2 between polarization angle and Stokes-space angle in KM's psi.
2. Sign of b vs S and of p vs k_AF not fixed (depends on SME gamma^5/eps conventions). Magnitudes are safe.
3. The boost bound assumes S_mu purely timelike in the CMB frame, constant and unscreened over the solar system, and minimal coupling to nucleons; nuclear-model dependence of neutron b~ is not propagated.
4. No numerical BBN bound found on the axial spin-torsion term kappa^2 s^2 in the EC Friedmann equation. Kranas 2019 bounds trace torsion only.
5. Not verified: Planck 2018 X axion isocurvature (Sec. 9.4.1), SPT-3G oscillation bound, Fedderke et al. washout original, Filzinger 2023 numbers, Cassini PPN gamma, CL of Dolan 2022 / Payez 2015 / Rogers 2023 / Kranas 2019, arXiv id of Data Tables ref [91].
6. Lammerzahl 1997 normalisation of K vs S_mu not checked (superseded, low priority).
7. KK: the relation Delta alpha/alpha <-> radion and the radion matter coupling strength are for the KK theorist; no KK-photon-specific collider bound collected.
