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

## kk05_edm_inputs.csv

Written by task `work/20261120-016-kk05-edm-cs`. Frequency band of the HfF+ oscillating-EDM search (arXiv:2006.15787 abstract only), used to show that a KK-05 field rolling on a Hubble time is a static EDM for lab purposes. The static eEDM bound itself is the JILA row in `shift_uholst_inputs.csv`.

## ec06_cosmo_inputs.csv

Written by task `work/20261120-023-cosmo-ec06`. Planck 2018 reionization mid-point z_re = 7.67 +- 0.73 and sigma8 = 0.8111 +- 0.0060 (arXiv:1807.06209 Table 2, TT,TE,EE+lowE+lensing; not re-opened), the fixed tanh-model settings (Delta z = 0.5, He II at z = 3.5), the radio-galaxy UV polarization rotation theta = -0.8 +- 2.2 deg at z > 2 (arXiv:1003.4823 abstract; CL not stated there), and the LiteBIRD-class forecast ~0.05 deg on beta_reio - beta_rec (arXiv:2108.09287 abstract; a forecast, not a measurement). Also the published torsion-free CFJ-equivalent numbers used as check targets: Yu-Yang-Qiu 2508.10532 Eqs. 53-54 (case II) and Li-Xia-Li-Zhang hep-ph/0611192 delta = -0.011 (WMAP3+B03, Sec. IV).

## kaf_lowz_bounds.csv

Written by task `work/20261120-024-lit-ec06-novelty`. Low-redshift (z ≲ 4) and local bounds on a constant or slowly varying CFJ coefficient, plus one CMB reference row. Narrative: `literature/ec06_prior_art.md` Sec. 6.

Columns: `id`; `quantity`; `symbol`; `value_published` (as printed); `units_published`; `central`, `err_1sigma_or_limit` (numbers used by the script; for upper limits `central` = 0 and the second column is the limit); `kind` (`meas_angle_deg`, `upper_p0_GeVh`, `meas_p0_GeVh`, `meas_p0_eVh`, `upper_kAF_GeV`, `upper_spatial_GeV`, `reference_noconv`); `z_min`, `z_max`, `z_eff` (sample redshift range and the effective redshift used for conversion); `CL`; derived columns `kAF_T_limit_eV`, `p0_limit_eV`, `ec06_xi_eps_limit`, `conversion_note` (written by `work/20261120-024-lit-ec06-novelty/kaf_lowz.py --write`, checked without `--write`); `experiment`, `arxiv`, `location_in_paper`, `note`.

Conversions (shown per row in `conversion_note`):
- Measurements are turned into limits as |central| + 2σ (Gaussian 1σ assumed where the CL is not stated).
- Angle to constant coefficient: |(k_AF)_T| = α / t_L(z_eff), with t_L the lookback time for flat ΛCDM, H0 and Ω_m from `cosmo_inputs.csv` (Kostelecky-Mewes 0905.0031 Eq. 158). p_0 = 2|(k_AF)_T| (`literature/bounds.md` Sec. 2).
- h0-scaled values use h = H0/100 from `cosmo_inputs.csv`. GeV to eV ×1e9.
- EC-06: p_0̂ = ξεH gives β(z) = (ξε/2) ln(1+z), so ξε < 2α/ln(1+z_eff).

Caveats:
- L-CFJ90 is taken as quoted in astro-ph/9704263 Eq. (5); the 1990 paper was not opened.
- Several angle rows come from compilations (1411.6287 Table 3, 1409.8242 Table 1), not the originals.
- L-WBB17-radio is Faraday-dominated per its authors and is not a birefringence bound.
- L-CMB-LCG26 is not converted, because its k-to-angle normalization differs from the one used here by a factor of about 6.
- GRB polarization gives no d = 3 bound (CFJ rotation is energy independent), so it has no row.
- The 1003.4823 θ value is also in `ec06_cosmo_inputs.csv` (task 023, abstract); here it is from the full text (Table 2, labelled 1σ).

## ec06_theory_inputs.csv

Written by task `work/20261120-022-theory-ec06`. Same columns as `cosmo_inputs.csv`: Planck 2018 Omega_K (TT,TE,EE+lowE+lensing+BAO; 68% error inferred from the quoted 95% interval) for the closed-FRW teleparallel axial branch, and the GW170817 tensor-speed interval (1710.05834 abstract) for the Magueijo-Zlosnik tensor-speed bound.

## ec06c_costs_inputs.csv

Written by task `work/20261120-025-theory-ec06c-costs`. Columns: `id`, `quantity`, `value`, `units`, `CL`, `source_arxiv`, `location`, `notes`. Contents: photon-mass limits (PDG 2024 adopted 1e-18 eV, Ryutov 07; Proca galactic-field 1e-26 eV, Adelberger-Dvali-Gruzinov hep-ph/0306245, and Chibisov 3e-27 eV, both listed but not used by PDG; listing page 1 opened); GW amplitude birefringence kappa < 0.74 Gpc^-1 and dCS l0 < 1e3 km (2101.11153 abstract); nondynamical CS gravity LAGEOS 1/m_cs < ~1000 km and the GP-B projection 0.01 km^-1 (0708.0001 Sec. VI, Fig. 2); IAU 2015 nominal solar/terrestrial constants and the au (not re-opened); Milky Way v_c and R0 (abstract values); and theory inputs: Shie-Nester-Yo 0805.3834 Table I cases III-V and the case II parameters (Sec. V.B), used as model parameters, not measurements.
