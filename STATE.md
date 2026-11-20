# State

Iteration: 2 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0 added in iter 2).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 2: `repin.ps1` re-pinned 2 python.exe mid-iteration and 1 after return. Cause: the `.venv` launcher stub that starts `pin.py` and short unpinned helper scripts (data task 010). Not a sustained load.
- `artifacts/` (gitignored) holds SPT candl_data (commit c4eb8fb), regenerate with `work/20261120-010-data-likelihood/download_candl_data.py`.
- `python -m src.reproduce`: currently reproduces the published β significances and SPT bounds (18 PASS). Does not yet regenerate claim numbers.

## Iteration 2 tasks
- 006 skeptic on iter 1: no sign or factor errors; all iter-1 formulas survive. EC-02 framing weakened (fixed in 007). β-viewpoint citation moved to 2202.13919 (CONVENTIONS.md updated).
- 007 theory no-go: EC-02, EC-05 dead with cited inputs; EC-04 is an axion with g_φγ ≈ 0. 29/29 PASS (rerun by overseer). Inputs in `data/cosmo_inputs.csv`.
- 008 cosmology of torsion axions: QCD-mass variants dead; generic ALP region mapped; washout reproduced from 2512.21888 to 1%. 25/25 PASS (rerun).
- 009 KK bounds: P' parity argument makes KK-01 an all-orders no-go; 2608.09299 analysed; KK-02, KK-03 live but ALP-in-disguise; radion × FF̃ dead. 33/33 PASS (rerun).
- 010 data likelihood: `src/` likelihood layer, candl SPT likelihood runs on 2 CPUs. 18/18 PASS (rerun).

## Settled (checked; skeptic-audited where noted)
- Iter-1 items: all survived skeptic 006 (EC-00, KK-01 tree level, CFJ sign, β = -(g/2)(φ_0 - φ_LSS), contortion, FRW torsion, S = -(κ/2)J₅, L_eff = -(3κ/16)J₅², CS reduction, δα/α formulas).
- Data: ACT DR6 2.905σ, ACT+WMAP+Planck 4.55σ, Cosmoglobe combination 3.9σ reproduced. Power-law β(ν) fit gives n ≈ 1.6 at 2.0σ, driven by the PA5 f150-f090 split the ACT paper attributes to instrument systematics. SPT: r, A_LT reproduce from candl; A_CB template scaled by 0.485 (fitted to the published 1.2×10⁻⁴ limit, so that limit is reproduced by construction; lensing-prior limit 4.5×10⁻⁵ vs published 5.3×10⁻⁵).
- Generic ALP baseline (008): washout 50% at m ≈ 1.1×10⁻²⁸ eV; allowed 10⁻³⁵ ≲ m ≲ 10⁻²⁷ eV, 7×10⁻²¹ ≲ g ≲ 8×10⁻¹³ GeV⁻¹. A_CB = 3.6×10⁻⁷(H_inf/φ_i)².

## Live candidate mechanisms (none has torsion doing the work)
- EC-03' torsion-scale pseudoscalar with the gluon coupling removed by assumption: g = 8α/(πf) = 4.7×10⁻²¹ GeV⁻¹, β_max = 0.136° at ρ_φ ≤ 5% ρ_DM; 0.215° needs ~12%. The 5% cap is an uncited assumption. Marginal.
- EC-04 PGT 0⁻ mode: shift-symmetric axion, g_φγ ≈ 0 unless a CFJ-type term is added by hand; needed mass takes ~10¹²⁰ tuning. Axion class.
- KK-02 Wilson-line axion: β = kαΔϑ/2π, k = 1 gives at most 0.209°. ALP in disguise.
- KK-03 KR 2-form + Green-Schwarz: g window 1.4-8×10⁻¹³ GeV⁻¹; natural λ = κ₄ excluded. ALP in disguise.

## Dead mechanisms
- EC-00 minimal EC + minimal Maxwell: β = 0 (theorem; skeptic-survived).
- EC-02 minimal EC spin source + CFJ: 28.5 orders short at ξ = 1 with cited ξ_ν ≤ 0.048 (2409.08280); ~35 orders for polarized fermion DM. Dead as a torsion mechanism; the torsion-free J₅·AF̃/M² operator at M ~ 3×10⁴ GeV (0711.4617, 2302.11140) is a separate published mechanism.
- EC-05 fermion-loop CFJ: ≥ 29 orders short; scheme ambiguity resolves to (H/m_e)² suppression on FRW.
- EC-03 torsion axion with the universal torsion-quark coupling (KR/torsion axion 1502.03694; Holst/Nieh-Yan BI fields with fermions): QCD mass 1.4-8.6×10⁻¹² eV, full washout, local oscillation ≤ 1.1×10⁻²⁰ deg, 19.3 orders short. Holst/Nieh-Yan BI without fermions: g = 0.
- KK-01 pure EC-KK: no parity-odd photon coupling at all orders (P' = diag(1,-1,-1,-1,-1) argument, 009; skeptic-proposed).
- KK-04 radion × FF̃ (2608.09299 teleparallel; EC analogue via u-Holst term, c_odd = 16c³/(1+12c²)): needs an oriented P'-odd background, and the radion must roll after recombination, which Cassini (γ-1 = -1/2) and the inverse-square-law bound (m > 5 meV) exclude.

Pending skeptic for 007, 008, 009 results.

## Open questions
- 008: the torsion-fermion coupling is shift-symmetric ∂φ·J₅. With massive fermions this may give no QCD potential and a photon coupling suppressed ~10⁻⁷⁰. If so, EC-03 dies by a different route and EC-03' has no torsion-generated g at all. Needs a theorist check.
- The u-Holst term (009) is the only place found where torsion itself generates a parity-odd photon coupling. Does it work with a field other than the radion (avoiding Cassini), and is the oriented background natural?
- EC-03' energy-density cap (5% ρ_DM) needs a cited bound for an early-dark-energy-like ultralight field.
- CS level normalization (k ∈ Z vs 6Z) quoted, not derived.
- Unresolved: Data Tables β→k_AF factor 4.47; Planck isocurvature bound unverified; no BBN bound on axial torsion.

## Direction
Both tracks are converging on a Result B / C structure: minimal EC and pure EC-KK give nothing; every live route is an axion-like field where torsion sets at most the coupling scale. Before a stop: skeptic on 007-009, literature check of the no-go, the u-Holst loose end, and `src.reproduce` regenerating every claim number.

## Next planned tasks (iter 3)
1. Skeptic on 007, 008, 009 (EC-02/05 numbers, QCD-mass argument, washout, P' argument, Cassini kill of KK-04, KK-02/03 windows).
2. Theorist: shift-symmetric ∂φ·J₅ question (does a QCD potential and a photon coupling arise; anomaly matching); u-Holst term with a non-radion field.
3. Literature scout: has the no-go (minimal EC and EC-KK give no birefringence; torsion axions reduce to ALPs) been published? Search specifically; also cite an EDE-like energy-density bound for EC-03'.
4. Data analyst: extend `src.reproduce` to regenerate the claim numbers in CLAIMS.md from data/ (calling the work/ scripts), and fit the live axion-class candidates (EC-03', KK-02) to the β likelihood vs null with Δχ²/AIC/BIC.
