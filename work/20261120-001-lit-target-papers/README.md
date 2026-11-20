# 20261120-001-lit-target-papers

## Question
What exactly do the four target papers (arXiv 2509.13654, 2510.07928, 2601.13624, 2512.21888) report: data, multipoles, frequencies, headline numbers with CL and location, frequency dependence, priors, alpha (miscalibration) handling, systematics, Chern-Simons conventions, anisotropy? Which earlier measurements do they cite, and what is the sign convention for beta?

## Method
- Full text read for all four papers from the arXiv HTML versions (converted to text with TeX math kept, `h2t.py`):
  - https://arxiv.org/html/2509.13654 (v2, 14 Apr 2026): full text.
  - https://arxiv.org/html/2510.07928 (v2, 14 Nov 2025): full text.
  - https://arxiv.org/html/2601.13624 (v2, 05 Aug 2026): full text.
  - https://arxiv.org/html/2512.21888 (v2, 13 Oct 2026): full text.
- Convention sources, full text: https://arxiv.org/html/2202.13919 (Komatsu 2022 review, cited by 2601.13624 for its angle convention), https://arxiv.org/html/2205.13962 (Eskilt & Komatsu), https://arxiv.org/html/2011.11254 (Minami & Komatsu; searched for a sign statement, none found).
- Abstracts only (export.arxiv.org API, `abs.py`) for the cited prior measurements listed in `literature/prior_measurements.md`.
- Figures were not read (HTML text only); all numbers come from text and tables.
- Rerun: `sh fetch.sh <scratchdir>`, then `.venv/Scripts/python tools/pin.py <python> -I h2t.py in.html out.txt`, `abs.py abs.xml abs.txt`, and `chk.py data/birefringence_measurements.csv` (CSV structure check).

Outputs: `literature/2509.13654.md`, `literature/2510.07928.md`, `literature/2601.13624.md`, `literature/2512.21888.md`, `literature/prior_measurements.md`, `literature/conventions_beta.md`, `data/birefringence_measurements.csv`, `data/README.md`.

## Result
Isotropic beta: every independent analysis gives the same sign and a value of roughly 0.2-0.35 deg. ACT DR6 gives beta = 0.215 +- 0.074 deg (68%, 2.9 sigma; 2509.13654 Table 3). Here alpha comes from an optics-model prior, so the error is set by the prior widths (0.09-0.11 deg), not by statistics (0.032 deg). Planck/WMAP give 0.342 +0.094 -0.091 deg (2205.13962) using Galactic-foreground alpha calibration. The authors' naive combinations are 0.264 +- 0.058 (4.6 sigma) and 0.231 +- 0.059 (3.9 sigma) with the Cosmoglobe reanalysis. Frequency dependence: ACT assumes achromatic beta. Its per-array total angles show 90 GHz below 150 GHz in both PA5 and PA6 (PA5 difference 0.21 +- 0.09 deg). This is treated as an unexplained instrumental issue, not a beta(nu) signal. Planck PR4 power-law fits give n consistent with 0 (e.g. n = -0.37 +0.49 -0.64 with dust EB model; 2201.13347 abstract), and WMAP+Planck shows no frequency dependence across 23-353 GHz. Anisotropy: no detection. The tightest BB-spectrum limit is SPT-3G A_CB < 0.53e-4 rad^2 (95%, with lensing prior; 1.2e-4 without; 2510.07928 Eq. 7-8). The tightest overall cited is SPTpol four-point A_CB < 0.10e-4 rad^2 (2006.08061). SPT BB data are blind to isotropic beta by construction (EB deprojection). The two theory papers take beta as input: 2601.13624 fits an EDE + Chern-Simons model (g M_Pl = 0.129 +0.031 -0.075 with ACT EB; f_EDE ~ 0.23, H0 ~ 77, tau ~ 0.135, S8 ~ 0.96, the last two in strong tension with other data). 2512.21888 finds single-field ALP DM in 1e-25-1e-22 eV excluded by washout and a fine-tuned two-field case marginally viable.

## Verdict
Literature record complete for the four target papers (full text) and their cited headline measurements (mostly abstracts). Isotropic beta ~ 0.2-0.35 deg is a candidate signal at 2.9 sigma (ACT alone) to ~4-4.6 sigma (combinations that assume independence), still systematics-limited by the polarization-angle calibration. Anisotropy is consistent with zero. Constraint for any mechanism: a nearly achromatic, isotropic rotation of order 0.2-0.3 deg with A_CB well below 1e-5-1e-4 rad^2.

## Open issues
1. Sign convention vs the repo: Komatsu 2022 uses tensor epsilon^{0123} = +1/sqrt(-g); `CONVENTIONS.md` uses epsilon^{0123} = -1/sqrt(-g). In that case our baseline L ⊃ -(g/4) phi F Ftilde corresponds to beta = -(g/2) Delta phi, not +(g/2) Delta phi. This needs a derivation in `work/` (`literature/conventions_beta.md` Sec. 4). 2601.13624 and 2512.21888 do not state their epsilon sign.
2. 2509.13654 itself does not name IAU vs COSMO/HEALPix. The COSMO identification is sourced to the Komatsu 2022 review.
3. 2601.13624 does not say how alpha was handled in the EB band powers it fits (ACT DR6 EB includes alpha + beta). Its Table 1 CL is not stated, and its M_Pl definition is not stated. It gives l <= 3400 in one place and 3500 in another.
4. 2512.21888 Eq. 16 upper bound (6.64e-14) is a decade above the single-field Eq. 6 value (6.64e-15). Possibly a typo, not verified.
5. Lonappan et al. A_CB limit: < 0.85e-4 as cited by 2510.07928, < 1e-4 in the 2504.13154 v3 abstract. Not resolved (body not read).
6. Cosmoglobe beta = 0.26 +- 0.10 (cited by 2509.13654) is not in the 2305.02268 abstract. Body not checked.
7. Prior-measurement numbers marked "abstract" were not checked against the original paper bodies. No WMAP-only beta and no ACT DR4/POLARBEAR/SPT isotropic beta are cited with numbers by these four papers.
8. Fedderke et al. 2019 (washout bound) arXiv ID not given in the citing reference list. Not looked up.
9. Public likelihoods (existence from the papers' footnotes; not downloaded or tested):
   - ACT DR6: spectra and covariances at act.princeton.edu/act-dr6-data-products (2509.13654 fn. 1). Analysis code at github.com/pdp79/act_dr6_analysis (fn. 2).
   - SPT: github.com/SouthPoleTelescope/spt_candl_data (SPT-3G D1 BB multi-frequency, 2510.07928 fn. 1). github.com/lbalkenhol/candl_data (SPT-3G D1 BB lite and SPTpol BB candl reimplementation, fn. 2).
   - Planck PR4: not linked in any of the four papers. Availability not verified in this task.
