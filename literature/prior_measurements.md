# Prior birefringence measurements cited by the four target papers

Source column: "abstract" = read from the original paper's arXiv abstract (fetched via export.arxiv.org API on 2026-11-20); "cited in X" = number as quoted in the citing paper's text. Angles in degrees, 68% CL unless stated. Machine-readable copies are in `data/birefringence_measurements.csv`.

## Isotropic beta

| Measurement | Data | Value | CL | arXiv | Source |
|---|---|---|---|---|---|
| Minami & Komatsu 2020 | Planck 2018 (PR3) HFI, alpha from Galactic foreground EB | 0.35 +- 0.14 (2.4 sigma) | 68% | 2011.11254 | abstract |
| Diego-Palazuelos et al. 2022 | Planck PR4, nearly full sky | 0.30 +- 0.11 | 68% | 2201.07682 | abstract (authors decline cosmological significance; beta falls with larger Galactic mask) |
| Eskilt 2022 | Planck PR4 LFI+HFI, all bands, no frequency dependence | 0.33 +- 0.10 | 68% | 2201.13347 | abstract |
| Eskilt 2022, power law beta_nu = beta_0 (nu/150 GHz)^n, f_sky 0.93, no FG EB | Planck PR4 HFI + 70 GHz | beta_0 = 0.26 +- 0.11, n = -0.45 +0.61 -0.82 | 68% | 2201.13347 | abstract |
| Eskilt 2022, same, filamentary dust EB model | | beta_0 = 0.33 +- 0.12, n = -0.37 +0.49 -0.64 | 68% | 2201.13347 | abstract |
| Eskilt 2022, all bands incl. 30, 44 GHz | | beta_0 = 0.29 +0.10 -0.11, n = -0.35 +0.48 -0.47 | 68% | 2201.13347 | abstract |
| Eskilt & Komatsu 2022 | WMAP 9yr + Planck PR4, 23-353 GHz, nearly full sky | 0.342 +0.094 -0.091 (3.6 sigma); no evidence of frequency dependence | 68% | 2205.13962 | abstract; also cited in 2509.13654 Sec. I/VI and 2601.13624 Sec. 1 |
| Eskilt et al. 2023 (Cosmoglobe DR1) | LFI + WMAP only (reprocessed) | 0.35 +- 0.70 | 68% | 2305.02268 | abstract |
| Eskilt et al. 2023 (Cosmoglobe DR1) | + HFI PR4, synchrotron-dominated channels <= 70 GHz | beta_{<=70GHz} = 0.53 +- 0.28 (1.9 sigma) | 68% | 2305.02268 | abstract |
| Eskilt et al. 2023 (Cosmoglobe reanalysis) | Cosmoglobe WMAP + Planck | 0.26 +- 0.10 | 68% | 2305.02268 | cited in 2509.13654 Sec. VI (not in the abstract; body not checked) |
| Sullivan et al. 2025 | Planck PR4 map-space, SEVEM | 0.46 +- 0.04 (stat) +- 0.28 (syst) | 68% | 2502.07654 | abstract; no miscalibration correction, consistent with zero given syst |
| Sullivan et al. 2025 | Planck PR4 map-space, Commander | 0.48 +- 0.04 (stat) +- 0.28 (syst) | 68% | 2502.07654 | abstract |
| ACT DR6 team (Louis et al. 2025) | ACT DR6 EB, mean array angle <psi_i> incl. optics systematic | 0.20 +- 0.08 (0.20 +- 0.03 stat only) | 68% | 2503.14452 | cited in 2509.13654 Sec. I, IV, VI |
| Diego-Palazuelos & Komatsu 2025 | ACT DR6 | 0.215 +- 0.074 | 68% | 2509.13654 | target paper |

Not found in these four papers: a WMAP-only beta (the closest is Cosmoglobe LFI+WMAP above), ACT DR4 isotropic beta, POLARBEAR or SPT isotropic beta. None of the four cites them with numbers. SPT BB analyses remove isotropic rotation by EB deprojection (2510.07928 Sec. 3.1).

## Forecast / calibrator-based (not a measurement)
- BICEP3 polarization-angle calibrator (BICEP/Keck XVIII, arXiv:2410.12089): no beta derived (absolute angle limited by a systematic); forecast on-sky sensitivity sigma(beta) = 0.078 deg (2017-2018 data), 0.055 deg with data through 2023 (abstract). Cited in 2509.13654 Sec. VI as a needed independent test.

## Time-dependent (oscillating) rotation
| Measurement | Value | CL | arXiv | Source |
|---|---|---|---|---|
| SPT-3G ALP oscillation, periods 1-100 days | amplitude < 0.071 deg | 95% median | 2203.16567 | abstract; cited in 2512.21888 footnote 1 |
| POLARBEAR Crab Nebula, 150 GHz | amplitude < 0.065 deg | 95% median | 2403.02096 | abstract; cited in 2512.21888 Sec. 1 |

## Anisotropic birefringence (A_CB = L(L+1) C_L^{alpha alpha}/2pi)
| Measurement | Method | Value | CL | arXiv | Source |
|---|---|---|---|---|---|
| Balkenhol et al. 2025, SPT-3G D1 BB | BB power spectrum | < 1.2e-4 rad^2; < 0.53e-4 with lensing prior | 95% | 2510.07928 | target paper |
| Lonappan et al. 2025 | BB from SPTpol, ACT, POLARBEAR, BICEP | best fit 0.42 +0.40 -0.34 e-4; < 1e-4 rad^2 | 68% / 95% | 2504.13154 | abstract (v3) |
| Lonappan et al. 2025 | BK + POLARBEAR + ACT, Planck-based A_LT prior | < 0.85e-4 rad^2 | 95% | 2504.13154 | cited in 2510.07928 Sec. 3.2. Differs from the v3 abstract (< 1e-4); probably a different version or data combination, unresolved |
| Namikawa 2024 | SPTpol BB, exact (non-thin LSS) | 1e4 A_CB = 1.03 +0.91 -0.97 | 2 sigma | 2404.13771 | abstract |
| Bianchini et al. 2020 | SPTpol 150 GHz, quadratic estimator (4-pt) | < 0.10e-4 rad^2 (0.033 deg^2) | 95% | 2006.08061 | abstract; also cited in 2510.07928 Sec. 3.2 |
| Gruppuso et al. 2020 | Planck 2018, scales > ~15 deg | A^{alpha alpha} < 0.104 deg^2 (marginalised); < 0.085 deg^2 at A^{alpha T}=0 | 95% | 2008.10334 | abstract |
| Contreras et al. 2017 | Planck | < [2.2 (stat) +- 0.7 (syst)] e-5 rad^2 = 0.07 deg^2; dipole sqrt(3C_1/4pi) = 0.32 +- 0.10 +- 0.08 deg excess | 95% | 1705.06387 | abstract |
| BICEP2/Keck IX 2017 (BK14) | 150 GHz, 4-pt | amplitude bound not in abstract; g_{a gamma} <= 7.2e-2/H_I | 95% | 1705.02523 | abstract |
| Namikawa et al. 2020 (ACT) | 2 seasons, 456 deg2 | amplitude not in abstract (2-3x improvement); g_{a gamma} < 4.0e-2/H_I | 2 sigma | 2001.10465 | abstract |
| Gluscevic et al. 2012 | WMAP-7 | sqrt of quadrupole of scale-invariant spectrum < 1 deg | 68% | 1206.5546 | abstract |
| POLARBEAR 2015 | 4-pt, first season | consistent with zero; quoted as PMF < 93 nG | 95% | 1509.02461 | abstract |
| Zagatti et al. 2024 | Planck PR3/PR4 EB estimator | consistent with zero at ~2 sigma up to L = 1500 (no amplitude in abstract) | | 2401.11973 | abstract |
| Bortolami et al. 2022 | Planck Legacy/NPIPE cross-spectra | A^{alpha E} = -7.8 +- 5.6 nK deg; A^{alpha B} = 0.3 +- 4.0 nK deg | 68% | 2206.01635 | abstract |

Unit note: 1 deg^2 = (pi/180)^2 rad^2 = 3.046e-4 rad^2. A_CB in 2510.07928 is in rad^2 (Eq. 5).
