# 20261120-010 data-likelihood

## Question
Build a minimal, reusable likelihood layer to fit candidate mechanisms against: isotropic beta (ACT DR6, ACT+WMAP+Planck as published, Planck PR4), beta(nu) = beta0 (nu/nu0)^n, and SPT anisotropic A_CB. Can it reproduce the published significances and limits, and which public likelihoods run on 2 CPUs?

## Method
Files:
- `src/data.py`: loads `data/*.csv` (no measured number hardcoded in src/).
- `src/likelihood.py`:
  - (a) `IsotropicBeta(dataset)`: Gaussian on a published beta (deg, Minami-Komatsu sign, beta > 0 clockwise, as in CONVENTIONS.md; optional split-normal for asymmetric errors). Datasets: `act_dr6` (DPK25-beta-baseline), `act_wmap_planck` (DPK25-comb-EK22), `act_cosmoglobe` (DPK25-comb-CG), `planck_pr4` (ESK22-beta, Planck PR4 LFI+HFI alone), `wmap_planck_pr4` (EK22), `cosmoglobe_wmap_planck` (CG23-beta-cited). Each carries the set of underlying data; `joint()` refuses to multiply likelihoods sharing data, so ACT x (ACT+WMAP+Planck) is impossible; only the published combination is used.
  - (b) `beta_powerlaw`, `PA5Difference` (PA5 f150 - f090 psi difference 0.21 +- 0.09, optionally plus the Pi^PA alpha-difference prior variance), `ACTPsiBands` (four published EB-only per-array-band psi_i = alpha_i + beta(nu_i), stat errors, Pi^PA alpha prior marginalised analytically: psi ~ N(beta(nu), C_stat + C_alpha)), `PlanckPowerLaw` (Eskilt 2022 beta0, n as independent split normals). Per-frequency data available: the four ACT psi_i (stat only, EB only, no PA4 f220), their mean, the PA5 difference, and Planck PR4 power-law fits (3 variants). No per-frequency beta with full systematics is published for ACT. Frequencies are the nominal band labels (`data/act_dr6_bands.csv`); effective band centres are not recorded.
  - (c) `SPTLite`: SPT-3G D1 BB lite (candl_data, offset log-normal) with D_l = r D^r + A_LT D^lens + A_CB D^CB, vectorised numpy re-implementation checked against candl to 1e-14 in ln L. Uniform priors r in [0,1], A_LT in [0,3], grid marginalisation (`acb_posterior`). `SPTApprox`: Gaussian on the published negative-allowed A_CB = (-1.13 +- 0.94)e-4 (fallback, no download needed).
- `src/templates.py`: D^r (r=1, n_t=0) and D^lens from CAMB at the Planck 2018 best fit (`data/planck2018_bestfit.csv`); D^CB from the full-sky thin-LSS formula evaluated exactly with Wigner-d correlation functions on a Gauss-Legendre grid (about 1 min, cached in `artifacts/templates/`).
- `src/compare.py`: Delta chi2 vs null, Wilks p-value and Gaussian sigma, AIC/BIC with explicit k and n, Savage-Dickey B_10 for a 1-parameter nested model (numeric, plus analytic Gaussian check), `summary_1param_vs_null`.
- `src/reproduce.py`: `python -m src.reproduce`, runs `validate.py` now.
- `download_candl_data.py`: clones github.com/Lbalkenhol/candl_data at commit c4eb8fb into `artifacts/` (gitignored).
- `check_templates.py`: independent checks of the CB template (Wigner-d orthogonality, isotropic limit C^BB = 4<alpha^2>C^EE, full-sky vs flat-sky convolution at l = 200-2000).

CB template normalisation: the published SPT analysis uses the exact (non-thin-LSS) template of Lonappan et al. / Namikawa 2024, which needs a modified CLASS build (github.com/toshiyan/biref-aniso-bb, no precomputed spectra). Namikawa 2024 (2404.13771, abstract) states the exact spectrum is about 2x below thin-LSS at l >~ 100. Our raw thin-LSS template gives a limit 0.485x the published one, consistent with that. `SPTLite` by default multiplies the template by one scale s = 0.485 fixed on the A_CB >= 0, no-A_LT-prior limit; the A_LT-prior limit and the negative-allowed posterior are then independent checks.

Run: `.venv/Scripts/python tools/pin.py .venv/Scripts/python -m src.reproduce` (1 thread; about 3-4 min on the shared 2 CPUs, first run adds about 1.5 min for templates). Setup: `download_candl_data.py` once.

## Result (run.log)
Isotropic beta:
- ACT DR6: Delta chi2 = 8.44, 2.905 sigma (published 2.9) PASS. Delta AIC = 6.4 (k=1 vs 0), B_10 = 6.3 for uniform beta in [-1, 1] deg.
- ACT+WMAP+Planck (EK22): 4.552 sigma (published 4.6) PASS; B_10 = 2.3e3. With Cosmoglobe: 3.915 (3.9) PASS.
- Recomputing both combinations from ACT and EK22 / Cosmoglobe by inverse variance gives 0.2646 +- 0.0578 and 0.2309 +- 0.0595 PASS.
- Planck PR4 alone (ESK22): 3.3 sigma (INFO).
Frequency dependence (INFO, approximate):
- psi bands stat-only weighted mean 0.1917 +- 0.0318 vs published 0.192 +- 0.032 PASS (loading check).
- psi bands + Pi^PA prior, achromatic: 0.193 +- 0.075 vs full EB+TB analysis 0.212 +- 0.072: the compressed reconstruction gets the error right and the mean 0.26 sigma low (no TB, no PA4 f220, no cross-band stat correlation).
- Power law fit to psi bands: beta0 = 0.26 deg, n = 1.6, Delta chi2 = 4.1 vs achromatic (2.0 sigma, 1 extra parameter). This is the PA5 f090/f150 split; the paper attributes it to unresolved instrument systematics, so it is not evidence for chromatic beta.
SPT A_CB:
- candl own test logL reproduced; numpy vs candl agree to 5e-15 PASS.
- r < 0.259 (published 0.26) PASS; A_LT = 0.67 +- 0.20 (0.69 +- 0.18) PASS. These test the lite data, windows, CAMB templates, priors and grid.
- A_CB < 1.2e-4 reproduced by construction (calibration); raw thin-LSS template gives 5.8e-5.
- With A_LT APS prior: 4.53e-5 vs published 5.3e-5 (ratio 0.85, tolerance 20%) PASS.
- Negative allowed: (-1.04 +- 0.90)e-4 vs (-1.13 +- 0.94)e-4 PASS; P(A_CB < 0) = 0.87 (0.89).
- Gaussian fallback truncated at 0: A_CB < 1.25e-4 (published 1.2e-4) PASS.
- check_templates.py: Wigner-d orthogonality 8e-15, isotropic limit exact to 3e-10, full-sky vs independent flat-sky convolution within 1.1% at l = 200-2000 (C_L^aa restricted to L >= 30 in both). The first version compared with L >= 1 and failed by 8-15%: at L ~ 1 the flat-sky integral is not a valid approximation of the sum (difference gamma - 1/2 + ln 2 in the log-divergent IR part), not an error in the template.
- validate.py: 18 PASS, 0 FAIL.

## Verdict
The likelihood layer reproduces the published isotropic significances (2.9, 4.6, 3.9 sigma) and combinations exactly, and the SPT-3G results to 15% after a single template normalisation that matches the documented thin-LSS overestimate. Usable now for fitting candidate mechanisms: isotropic beta (published values), beta(nu) (approximate), A_CB (SPT-3G lite, or Gaussian fallback).

## Public likelihoods on 2 CPUs
- Runnable: SPT-3G D1 BB lite via candl (candl-like 2.2.0, pure numpy without JAX, <1 ms per evaluation; 3D grid marginalisation about 15-60 s depending on load). SPTpol BB lite is in the same candl_data clone (not wrapped yet). CAMB single runs (about 10-20 s) for templates.
- Not runnable here this iteration: the exact CB template (modified CLASS, C build, not attempted); ACT DR6 EB/TB beta likelihood (multi-GB spectra plus covariances from LAMBDA and the pdp79/act_dr6_analysis pipeline; not attempted, would exceed the CPU and disk budget for an iteration); Planck PR4 / WMAP EB analyses (map-level, not feasible). These are replaced by published summary values.

## Open issues
- A_CB normalisation is anchored to one published number. An exact template (biref-aniso-bb) would remove the calibration; the residual 15% gap on the A_LT-prior limit is likely the l-dependence of the thin-LSS error (Namikawa: larger at low l).
- beta(nu): no ACT per-frequency beta with systematics exists; `ACTPsiBands` is an approximation (no cross-band correlation, EB only). Nominal 90/150 GHz frequencies, not effective band centres.
- PlanckPowerLaw ignores the beta0-n correlation (not published).
- EK22 asymmetric errors are symmetrised by default (split=True available); the published 4.6 sigma combination uses symmetric errors.
- `data/planck2018_bestfit.csv` values were not re-opened from 1807.06209 in this task.
- CPU pinning: the coordinator's repin.ps1 flagged 2 unpinned python.exe from this task. The pinned jobs were launched as `.venv/Scripts/python tools/pin.py ...`; on Windows the `.venv/Scripts/python.exe` launcher stub that starts pin.py is itself not pinned (it only waits, PID observed with full affinity while the child chain was on CPUs 12-13). Also unpinned: `download_candl_data.py` (git only) and short one-off helper invocations (`python -c` import checks, BOM stripping, file edits), each a few seconds.
