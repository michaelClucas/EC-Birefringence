# data/

## birefringence_measurements.csv

One row per reported number. Written by task `work/20261120-001-lit-target-papers` from the full text of arXiv 2509.13654, 2510.07928, 2601.13624, 2512.21888, and the arXiv abstracts of the papers they cite (see `literature/prior_measurements.md`).

Columns:
- `id`: unique row id, `<first-author-tag><yy>-<short>`.
- `quantity`: what the number is, with units in the name. `beta_iso_deg` isotropic rotation (deg, Minami-Komatsu sign convention, positive = clockwise; see `literature/conventions_beta.md`); `psi_alpha_plus_beta_deg` total angle alpha_i + beta for one array-band; `alpha_miscal_deg` instrument miscalibration angle; `alpha_prior_sigma_deg` prior width; `A_CB_rad2` = L(L+1)C_L^{alpha alpha}/(2 pi) in rad^2; `A_alphaalpha_deg2` same in deg^2; `beta0_150GHz_deg`, `beta_powerlaw_index_n` for beta_nu = beta0 (nu/150 GHz)^n; `beta_osc_amplitude_deg` amplitude of a time-oscillating rotation; `g_MPl`, `f_EDE`, `H0_kms_Mpc`, `A_LT`, `r` model parameters; `g_phigamma_GeV-1_at_1e-24eV` coupling coefficient (scaling in `notes`); `washout_1_minus_J0` polarization washout fraction; `chi2_ACT_EB_38pts` goodness of fit.
- `value`: central value, or the limit when `upper_limit_flag` = 1.
- `err_minus`, `err_plus`: lower and upper uncertainty (positive numbers), at the CL given. Empty for limits.
- `CL`: confidence level in percent (68, 95). Empty when not applicable.
- `upper_limit_flag`: 1 if `value` is an upper limit.
- `dataset`, `arxiv`: data and source paper.
- `location_in_paper`: table/equation/section; "abstract" means the arXiv abstract; "cited in X" means the number was taken from paper X's text, not the original.
- `notes`: priors, caveats, conversions.

Caveats:
- Rows with `location_in_paper` = abstract were not checked against the body of the original paper.
- Planck PR4 map-space rows (SUL25) carry the systematic error only in the error columns; stat error is in `notes`.
- YIN26 CL is not stated in their Table 1 caption; 68% is assumed.
- Numbers are in the authors' own conventions. Sign conversion of beta to the repo Lagrangian convention is pending (see `literature/conventions_beta.md` Sec. 4).

Scripts must load numbers from this file rather than hardcoding them.

## bounds.csv, bounds_inputs.csv

Written by task `work/20261120-005-bounds`. Narrative and mapping conventions: `literature/bounds.md`.

`bounds.csv` columns: `id` (prefix F fermion SME, T torsion, P photon CFJ, A axion-photon, U ultralight axion, K KK-specific, C EC cosmology), `quantity`, `symbol`, `value_published` (as printed; `a +- b`, `< x`, `> x`, or an interval), `units_published`, `value_eV` (converted; dimensionless quantities unchanged; couplings in eV^-1), `CL`, `sign` (one-sided upper/lower or two-sided), `experiment`, `arxiv`, `location_in_paper` ("abstract" = only the arXiv abstract was opened; "0801.0287v19 Table Dxx" = read from the Kostelecky-Russell Data Tables), `conversion_note`.

`bounds_inputs.csv`: unit constants (exact SI/IAU definitions), Planck 2018 CMB dipole and H0/Omega_m used for derived bounds, and two theory inputs (g*, T_BBN) marked as such.

## planck2018_bestfit.csv, act_dr6_bands.csv

Written by task `work/20261120-010-data-likelihood`.

- `planck2018_bestfit.csv`: Planck 2018 base-LCDM parameters (arXiv:1807.06209 Table 2, TT,TE,EE+lowE+lensing column), used only to compute the BB templates (r, lensing, anisotropic-birefringence) for the SPT likelihood in `src/templates.py`. Values not re-opened in task 010.
- `act_dr6_bands.csv`: ACT DR6 polarized array-bands with nominal frequency (the band label) and band edges quoted in arXiv:2509.13654 Sec. II. Used by the beta(nu) likelihoods in `src/likelihood.py`. Effective band centres are not recorded; the nominal labels are an approximation for chromatic fits.
- Row `DPK25-alphaprior-corr` (0.9) in `birefringence_measurements.csv`: correlation of the Pi^PA alpha priors within a PA (2509.13654 Sec. II), added in task 010.

## ultralight_fraction_bounds.csv

Written by task `work/20261120-013-lit-nogo-check`. Per-mass bounds on the energy density today of an ultralight (DE-like or DM-like) scalar, for the EC-03' energy cap. Columns: `id`, `quantity`, `mass_eV` (single mass or range), `value` (as printed), `units`, `CL` (%), `bound_type`, `data_combination`, `arxiv`, `location_in_paper`, `note`. Fractions of ρ_DM in `note` use the Ω_c h² of the same table row or 0.12 as stated. Two rows (`ULF-P18-zc`, `ULF-NNQ25-consistency`) are not upper bounds: they record that no bound exists or that a region is consistent with data. Table values were extracted from the arXiv HTML. Narrative: `literature/nogo_literature_check.md`.
