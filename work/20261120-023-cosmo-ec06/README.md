# 20261120-023 cosmo-ec06 (EC-06, iteration 6)

The ACT DR6 β = 0.215° ± 0.074° (arXiv:2509.13654) is 2.9σ, has systematics caveats, and may not be cosmological. EC-06 is a candidate mechanism; there is no `RESULT.md`.

## Question
Candidate class EC-06: a homogeneous background torsion of order H, not sourced by spin, plus an added CFJ photon coupling -½ξ p_μA_νF̃^{μν} with p = S (axial) or T (trace). For parameterized backgrounds:
1. β_total (last scattering to today) and β_reio (the part seen by reionization-era low-ℓ EB) as functions of ξε; check against skeptic 020 and analytic integrals.
2. How β accumulates in z, contrasted with ALPs; is this an observational distinguisher, and can current data test it?
3. Fit ξε to the β likelihoods; Δχ², AIC/BIC with honest parameter counts, comparison with a free-β ALP.
4. Anisotropy A_CB from torsion perturbations vs SPT-3G.
5. Low-z astrophysical k_AF bounds.

## Method
- **Bookkeeping** (CONVENTIONS.md: β = +½∫p_τ dτ). With the FRW tetrad e^0̂ = a dτ, p_τ = a p_0̂ (lower indices; p^0̂ = -p_0̂), and dt = a dτ, so β = ½∫p_0̂ dt = ½∫₀^{z_s} p_0̂ dz/((1+z)H). Symbolic check in `check_symbolic.py` (7 PASS), including dp = 0 for a homogeneous timelike p (CFJ term gauge invariant up to a total derivative on FRW only; off FRW this needs dp = 0, a theory question for 022).
- **Shapes.** p_0̂ = ξε H(z)(1+z)^n with n = 0 (tracks H), ±0.5, ±1, and p_0̂ = ξε H₀ constant (k_AF-like). ξε is one number: ξ times the background amplitude ε in S_0̂ (or T_0̂) = εH(1+z)^n. Analytic: n = 0 gives β = (ξε/2) ln(1+z_s); n ≠ 0 gives ξε((1+z_s)^n - 1)/(2n); constant gives (ξε/2)H₀t_L.
- **Background.** Flat ΛCDM with radiation from task 008 (`common.py`, Planck 2018 inputs from `data/`). Two independent routes: cosmic-time quadrature over z, and ½∫a p_0̂ dτ on 008's conformal-time table.
- **Reionization.** tanh x_e(z) (Planck/CAMB convention: Δz = 0.5, He I with H, He II at 3.5), z_re = 7.67 from `data/ec06_cosmo_inputs.csv` (1807.06209 Table 2). Check: it gives τ = 0.0543 vs Planck 0.0544. β_reio is quoted three ways: point value at z_re, average over the reionization visibility g = τ′e^{-τ} (z < 40, includes late He II scattering), and β(z = 10).
- **ALP contrast.** 008's field ODE (`field_evolution.py`), β(z_s) ∝ φ₀ - φ(z_s), for m = 10⁻³³ to 10⁻²⁹ eV.
- **Fits** (`ec06_fit.py`). `src.likelihood.IsotropicBeta` (`act_dr6`, `act_wmap_planck` = ACT + WMAP/Planck PR4 EK22, `act_cosmoglobe`), `src.compare.summary_1param_vs_null`. Savage-Dickey prior: uniform β ∈ [-1°, 1°] mapped to ξε (as for the ALP in 014).
- **Anisotropy** (`ec06_aniso.py`). Limber projection with the linear P(k) (Eisenstein-Hu 1998 no-wiggle transfer, A_s and n_s from `data/`; check σ₈ = 0.832 vs Planck 0.811). Two assumptions for how the torsion follows perturbations (below).
- **k_AF** (`ec06_kaf.py`). |k_AF| = |p|/2 (literature/bounds.md Sec. 2) vs `data/bounds.csv` rows P-kAF-CFJ90, P-kAF0-pulsar, P-kAF-Caloni; and the direct rotation from z = 2 vs the radio-galaxy UV polarization θ = -0.8° ± 2.2° at z > 2 (1003.4823 abstract, added to `data/ec06_cosmo_inputs.csv`).
- **Published maps, low-z file, theory table** (`ec06_published_maps.py`). Maps Yu-Yang-Qiu 2508.10532 (case II) and Li-Xia-Li-Zhang hep-ph/0611192 onto ξε; recomputes the ξε limits of task 024's `data/kaf_lowz_bounds.csv` with this task's β(z); evaluates 022's `ec06_backgrounds.csv` where a row maps onto a shape. Published numbers used as check targets are in `data/ec06_cosmo_inputs.csv`.

Rerun (about 65 s, 1 core):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-023-cosmo-ec06/run_checks.py > work/20261120-023-cosmo-ec06/run.log
```
Last run (`run.log`): **104 PASS, 0 FAIL**.

## Result

### 1. β per unit ξε and the reionization share (`beta_shapes.csv`)
| shape p_0̂ | β/(ξε) | ξε for 0.215° | β_reio/β_rec at z_re / visibility-avg / z=10 | β_rec − β_reio (vis) | p_0̂ today |
|---|---|---|---|---|---|
| ξεH (n = 0) | 200.4° (½ln 1091 rad) | **1.073×10⁻³** | 0.31 / 0.23 / 0.34 | 0.166° | 1.5×10⁻³⁶ eV |
| ξεH(1+z)⁻¹ | 28.6° | 7.5×10⁻³ | 0.89 / 0.77 / 0.91 | 0.050° | 1.1×10⁻³⁵ eV |
| ξεH(1+z)^-½ | 55.6° | 3.9×10⁻³ | 0.68 / 0.55 / 0.72 | 0.097° | 5.6×10⁻³⁶ eV |
| ξεH(1+z)^+½ | 1835° | 1.2×10⁻⁴ | 0.06 / 0.04 / 0.07 | 0.206° | 1.7×10⁻³⁷ eV |
| ξεH(1+z)^+1 | 31224° | 6.9×10⁻⁶ | 0.007 / 0.004 / 0.009 | 0.214° | 9.9×10⁻³⁹ eV |
| ξεH₀ const | 27.3° | 7.9×10⁻³ | 0.95 / 0.85 / 0.97 | 0.033° | 1.1×10⁻³⁵ eV |

- Skeptic 020's ξε ≈ 1.07×10⁻³ for S_0̂ = εH is reproduced (1.0729×10⁻³). Analytic and numeric agree to 10⁻⁸; cosmic-time and conformal-time routes agree to < 6×10⁻⁵. The constant case gives p₀ = 1.13×10⁻³⁵ eV (t_L = 13.81 Gyr), matching `literature/bounds.md` (1.1×10⁻³⁵ eV).
- Sign: β > 0 needs ξε > 0 with ε defined on the lower hat component (S_0̂ = -S^0̂).

### 2. Accumulation in z and the ALP contrast (`beta_accumulation.png`, `alp_contrast.csv`)
**The premise in the task brief is wrong for n = 0.** A CFJ vector that tracks H gives dβ = (ξε/2) d ln(1+z): the rotation is spread uniformly in ln(1+z), so 66-77% of it is accumulated between reionization and recombination. The reionization-era EB sees only **0.23-0.31 of β** (0.05-0.07° for β_rec = 0.215°). Only shapes that decay faster than H into the past (n ≲ -1) or a constant p put most of the rotation at z < 10 (0.77-0.95).

ALPs (008 ODE, visibility-averaged): m = 10⁻³³ eV gives β_reio/β_rec = 0.95; 10⁻³² eV 0.80; 10⁻³¹ eV −0.06; ≥ 10⁻³⁰ eV ≈ 0 (frozen through recombination, decayed before reionization). In between, 2×10⁻³² eV gives 0.57, 3×10⁻³² eV 0.38, 5×10⁻³² eV 0.13: an ALP with m ≈ 4×10⁻³² eV reproduces the n = 0 ratio 0.23.

So the distinguisher is the **ratio R = β_reio/β_rec**, and for a fixed EC-06 shape it is a prediction: n = 0 predicts R ≈ 0.23-0.31 (β_rec − β_reio ≈ 0.17° if β_rec = 0.215°). It does not separate EC-06 from ALPs as a class, because an ALP mass can be tuned to any R in [0, 1]. What separates them is the full β(z_s) curve (smooth log accumulation vs the ALP's step at z ~ m/H), which two numbers (β_reio, β_rec) do not resolve.

**Current data.** I found no published measurement of the reionization-era isotropic angle from Planck low-ℓ EB. Sherwin & Namikawa (2108.09287, abstract) forecast ~0.05° on β_reio − β_rec for a LiteBIRD-class satellite and no Planck number; Nakatsuka-Namikawa-Komatsu (2203.08560, abstract) say this test needs a full-sky satellite. Namikawa-Murai-Naokawa (2506.20824, abstract) constrain ALP masses from Planck high-frequency EB, not a reionization angle. So current data cannot tell n = 0 from a late-rolling ALP. A LiteBIRD-class 0.05° error would separate n = 0 (Δβ ≈ 0.17°) from the constant-p shape or a late ALP (Δβ ≈ 0.01-0.03°) at about 3σ, if β_rec ≈ 0.2°.

### 3. Fits (`fit_results.csv`)
β is linear in ξε, so each fit is the free-β fit relabelled: same χ²_min, same Δχ², same B₁₀ (checked for every shape and data set).

| data | free-β ALP best [68%] | EC-06 n = 0: ξε best [68%] | const: ξε | Δχ² vs null | σ | ΔAIC k=1 | ΔAIC k_full=2 | ΔBIC (n=1) | B₁₀ |
|---|---|---|---|---|---|---|---|---|---|
| ACT DR6 | 0.215° [0.141, 0.289] | 1.07×10⁻³ [0.70, 1.44]×10⁻³ | 7.9×10⁻³ [5.2, 10.6]×10⁻³ | 8.44 | 2.91 | 6.44 | 4.44 | 8.44 | 6.3 |
| ACT + WMAP/Planck (EK22) | 0.264° [0.206, 0.322] | 1.32×10⁻³ [1.03, 1.61]×10⁻³ | 9.7×10⁻³ [7.6, 11.8]×10⁻³ | 20.72 | 4.55 | 18.72 | 16.72 | 20.72 | 2293 |
| ACT + Cosmoglobe | 0.231° [0.172, 0.290] | 1.15×10⁻³ [0.86, 1.45]×10⁻³ | 8.5×10⁻³ [6.3, 10.6]×10⁻³ | 15.33 | 3.92 | 13.33 | 11.33 | 15.33 | 158 |

Other shapes: rescale by the β/(ξε) column above (all in `fit_results.csv`).
- Parameter counts: k_eff = 1 (only ξε enters β). k_full = 2 for EC-06 with the shape fixed by theory (ξ and ε are separate; data fix the product), 3 if n is also free (ΔAIC drops by 2 more). Free-β ALP k_full = 2 (task 021). With equal k_full, EC-06 and the ALP tie exactly; EC-06 with free n loses by 2.
- BIC on one compressed β has ln n = 0, so ΔBIC = Δχ² with no penalty (as in 014).
- The combined Δχ² assumes the observed rotation is the same angle in ACT and WMAP/Planck. For n = 0, Planck/WMAP's low-ℓ part would see a smaller angle; that information is not in the compressed published combination, so this is a small inconsistency in the combined fit, not a gain.

### 4. Anisotropy (`acb_estimate.csv`)
A homogeneous p gives A_CB = 0 at linear order. If the torsion tracks the local expansion:
- **(A) kinematic**, p_0̂ = ξε θ/3 with δθ/3 = −fHδ_m/3 (linear continuity): max L(L+1)C_L/2π = **4.9×10⁻¹¹ rad²** (at L ≈ 770), 2.4×10⁻¹¹ at L = 100. Cross-checked against a coherent-cell estimate of the variance (agree within a factor 5).
- **(B) potential**, δp/p = Φ: max 3.4×10⁻¹⁴ rad² (L = 2); a Sachs-Wolfe-type estimate (βΦ)² ~ 10⁻¹⁴ agrees.
Both are ≥ 2×10⁶ below SPT-3G A_CB < 1.2×10⁻⁴ (95%). Non-linear P(k) at L ≳ 300 could raise (A) by O(10-100), still ≥ 10⁴ below. Assumptions: n = 0 shape at the ACT best fit, torsion follows linear matter perturbations, no independent torsion fluctuations (e.g. an inflationary spectrum of a torsion degree of freedom), no spatial p_i components. If the theory has its own torsion perturbations (a dynamical torsion mode), this bound has to be redone with that spectrum.

### 5. Low-z k_AF
At ξε for ACT, |k_AF| today = |p_0̂|/2 = 5×10⁻³⁹ to 5.7×10⁻³⁶ eV across shapes. Radio-galaxy CFJ bound (bounds.csv P-kAF-CFJ90) < 6×10⁻³⁴ eV: 100-10⁵ times weaker. Pulsar k⁰_AF < 10⁻²⁷ eV: irrelevant. The rotation accumulated from z = 2 is 0.0004-0.16° vs the radio-galaxy UV error 2.2° (1003.4823): ≤ 0.075 of it. The CMB constant-k_AF bound (P-kAF-Caloni, < 7.4×10⁻³⁶ eV) is the same CMB observable; the constant shape's 5.7×10⁻³⁶ eV sits just below it, consistent with it. **No low-z astrophysical probe constrains EC-06.** The literature scout (task 024) added `data/kaf_lowz_bounds.csv` while this task ran (not used in the scripts; it was still being written). Its tightest low-z row (L-AFG10-UV, ξε < 0.136 in its own conversion) is ~130 times above the n = 0 requirement, so the conclusion stands. Its CMB reference row L-CMB-LCG26 (2611.23082, |k_AF| < 3.82×10⁻³⁶ eV, 95%, abstract says "driven by the SPT-3G D1 B-mode data") is below the constant-shape value 5.7×10⁻³⁶ eV in our p = 2k_AF conversion. It is a bound on the CMB rotation itself, so it applies to every isotropic-β mechanism (ALPs included), not to EC-06 specifically; whether it holds with the miscalibration angle free, and the ×4.5 β↔k_AF conversion discrepancy in `literature/bounds.md`, are open.

### 6. Torsion-free published equivalents (`ec06_published_maps.py`)
Both are p ∝ ∂ ln(curvature) CFJ couplings, so β is an exact total derivative, and both have a matter-era p_0̂ ∝ H (the n = 0 shape locally) that switches off in the Λ era.
- **Yu-Yang-Qiu 2508.10532**, case II: (α/M²)·M² ln(Q/M_P²) F̃F with Q = 6H². In our conventions gφ = −4α ln Q, so p_0̂ = 4α d ln(6H²)/dt (matter era −12αH) and |β| = 4α ln E(z*) = 40.21α rad. This reproduces their Eq. 53 (40.21α) and their Eq. 54 range 8.23×10⁻⁵ < α < 1.78×10⁻⁴ from Planck PR4 0.30° ± 0.11° (we get 8.25×10⁻⁵ to 1.78×10⁻⁴). ACT DR6 gives **α = 9.3×10⁻⁵**, inside their range at its low edge. Same total β as our n = 0 shape for ξε = 11.5α = 1.073×10⁻³. Their fit agrees with ours. Sign not compared (their ε and F̃ conventions not checked).
- **Li-Xia-Li-Zhang hep-ph/0611192**: Δχ = (δ/2) ln(R(z)/R(0)) (their Eq. 42), so p_0̂ = δ d ln R/dt (matter era −3δH) and |β| = 9.35|δ| rad (R includes only matter and Λ). ACT gives **|δ| = 4.0×10⁻⁴**, equivalent to ξε = 2.67|δ| = 1.073×10⁻³. Their WMAP3+B03 fit δ = −0.011 (+0.0074 −0.0079) means |β| ≈ 5.9°. That fits the old WMAP3+B03 data and is excluded by the current 0.2-0.3° measurements, so it does not agree with our value; our ACT-based |δ| is 27 times smaller and sits inside their 2σ interval (they quote |δ| < 10⁻² there).
- Reionization share: ln 6H² gives 0.27 (z_re) / 0.18 (visibility); ln R gives 0.23 / 0.14. These are slightly below n = 0 because p switches off in the Λ era.
- So EC-06 with n = 0 has the same β history as these torsion-free couplings to within ~20% in the reionization share. A torsion background adds nothing observable to the CMB β over them. It only changes where the coupling comes from and what else bounds it.

### 7. Low-z bounds from task 024 (`data/kaf_lowz_bounds.csv`) with this task's β(z)
| row | z_eff | angle limit | ξε limit n = 0 (ours / scout) | n = −1 | const | margin vs ACT n = 0 |
|---|---|---|---|---|---|---|
| L-AFG10-UV (1003.4823) | 2.80 | 5.2° | 0.136 / 0.136 | 0.246 | 0.229 | 127× |
| L-WPC97-UV | 0.81 | 3.6° | 0.212 / 0.212 | 0.281 | 0.257 | 197× |
| L-Leahy-radio | 0.47 | 5.2° | 0.471 / 0.471 | 0.568 | 0.528 | 439× |
p_0 rows (applied as p_0̂ today = ξεH₀): CFJ90 ξε < 0.8, CF97 < 2.2. L-GGP15-comb (1411.6287) gives ξε < 0.032, but it includes the WMAP9 CMB, so it is not a low-z bound. Confirmed: **low-z data do not test EC-06** (tightest margin 127×, scout's numbers reproduced to < 1%).

### 8. Theorist table (022 `ec06_backgrounds.csv`)
- TEGR / flat f(T), trace T_t = 3H (n = 0, ε = 3): ξ = 3.58×10⁻⁴, matching 022 and skeptic 020.
- Closed f(T), axial branch, p_0̂ = 2ξ√|Ω_K| H₀(1+z): ξ = 1.83×10⁻² at |Ω_K| = 0.0043 (my 95% upper value from `data/ec06_theory_inputs.csv`), = 2.19×10⁻² at 022's 0.0030 (scales as |Ω_K|^-½; consistent). β_reio/β_rec = 0.65 at z_re.
- Magueijo-Złośnik tracker: with 022's GW-allowed ε_S ≤ 1.5×10⁻⁷, ξ ≥ 7×10³ (022: dead). The γ → ∞ branch, PGT 0+ (oscillating Φ) and PGT 0− (ALP) are not single shapes here; see 022.

## Verdict
- EC-06 is a viable candidate mechanism on the cosmology side: ξε = 1.07×10⁻³ (n = 0) to 7.9×10⁻³ (constant p), or as small as 7×10⁻⁶ for p growing as H(1+z), fits ACT with Δχ² = 8.4 (ACT) / 20.7 (ACT+WMAP/Planck), identical to a free-β ALP. It does not beat a free-β ALP: equal at k_full = 2, worse by ΔAIC = 2 if the shape is free.
- A_CB and low-z k_AF bounds do not constrain it (≥ 10⁶ and ≥ 10² margins). The constraints that matter are on ξ itself and on the background ε (theory 022, bounds on the torsion's other couplings).
- The n = 0 cosmology is the same as the published torsion-free couplings ∝ ∂ ln H² (2508.10532, α = 9.3×10⁻⁵ for ACT, agreeing with their Planck-based range) and ∂ ln R (hep-ph/0611192, |δ| = 4.0×10⁻⁴; their 2007 WMAP3 value is 27 times larger and fit old data). Prior art for the photon-sector phenomenology, not new.
- Testable prediction: β_reio/β_rec is fixed by the shape (0.23-0.31 for n = 0). Not testable with current data; a LiteBIRD-class satellite reaches ~3σ for n = 0 vs a late ALP. It cannot exclude ALPs as a class (an ALP mass can match any single ratio).

## Open issues
- 022's table arrived late; only rows that map onto a parameterized shape were evaluated (Sec. 8). Not used: PGT 0+ oscillating Φ history.
- Signs of the 2508.10532 and hep-ph/0611192 maps were not compared (only magnitudes).
- β_reio is a scalar proxy. The real low-ℓ EB from a z-dependent rotation needs a Boltzmann code with a time-dependent rotation (not run, per CPU rules).
- The combined ACT+WMAP/Planck fit treats the rotation as one angle; for strongly z-dependent shapes the Planck low-ℓ share differs. Effect not quantified (Planck β is dominated by ℓ > 30, recombination).
- Anisotropy assumes the torsion follows matter perturbations; a theory with an independent torsion mode needs its own spectrum. Spatial components p_i (from velocity flows) were not included; their effect was not estimated.
- 1003.4823: CL of the ±2.2° not stated in the abstract; individual source redshifts not opened.
