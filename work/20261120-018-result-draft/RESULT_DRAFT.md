# RESULT draft, revision 2 (not RESULT.md)

Revised by the overseer in iteration 7 after skeptic 020 (blocker + majors + minors), tasks 019-024. Becomes `RESULT.md` only after a fresh skeptic fails to break it. Iteration-7 tasks 025-027 folded in.

## Question
Can Einstein-Cartan (ECSK) torsion, in 4D or in an Einstein-Cartan extension of Kaluza-Klein theory, explain the cosmic birefringence reported in ACT DR6 (β = 0.215° ± 0.074°, arXiv:2509.13654) and constrained by SPT-3G (A_CB < 1.2×10⁻⁴, arXiv:2510.07928), in the context of arXiv:2601.13624 and 2512.21888?

Caveat on the target: the ACT β is 2.9σ (4.6σ combined with WMAP and Planck; 3.9σ with the Cosmoglobe reanalysis), with the miscalibration angle set by an optics-model prior and dust EB neglected. ACT per-band angles prefer a power law n ≈ 1.6 over achromatic at 2.0σ (driven by the PA5 f150-f090 split, attributed to instrument systematics), while Planck PR4 gives n = -0.35 ± 0.48. The signal may not be cosmological.

## Answer
- **4D track: Result B for Einstein-Cartan proper; Result C only outside EC.** Minimal EC gives β = 0. Spin-sourced torsion (any coupling: CFJ, fermion loop, gauge-invariant higher-dimension operators) falls ≥ 28 orders of magnitude short. Torsion pseudoscalars coupled through ∂φ·J₅ are exactly shift symmetric and have no low-energy photon coupling. Unsourced O(H) torsion: the Magueijo-Złośnik EC extension needs ξ ≥ 1.4×10⁴ (GW speed, on its tracker), disfavoured but not excluded; TEGR torsion is not defined by TEGR's dynamics; the PGT 0⁺ (SNY) torsion gives a closed CFJ vector but only in histories excluded by H(z*), otherwise it reduces to the published ∂ ln R coupling. The one surviving torsion-based candidate mechanism is f(T) teleparallel trace torsion with a CFJ coupling ξ = 3.6×10⁻⁴ (EC-06c). Its costs: f(T) gravity with f_TT ≠ 0 (not EC; local Lorentz invariance broken, the frame fixed by f(T)'s Lorentz equation); a P-odd photon coupling added by hand; U(1) gauge invariance broken everywhere on solutions, repaired only by a Stückelberg photon mass in a narrow window 7×10⁻²⁸ ≲ m_γ ≲ 10⁻²⁶ eV; a naive EFT cutoff of ~20-60 GeV; and phenomenology identical to published torsion-free models (2508.10532). Other survivors are ALPs with a photon coupling added beyond EC, where torsion at most sets f.
- **Einstein-Cartan-Kaluza-Klein track: Result B.** Pure EC-KK has no parity-odd photon coupling at any order in a P'-invariant state (5D parity P' = diag(1,-1,-1,-1,-1)); P'-odd matter backgrounds are Planck suppressed like their 4D analogues. Radion × FF̃ (including 2608.09299) is excluded by Cassini and the inverse-square law. The only route where torsion itself generates the coupling (KK-05, u-Holst term) has no consistent charged matter: the CMB photon must be the KK photon, and no light state carries KK charge e. Surviving KK candidate mechanisms (5D Chern-Simons Wilson-line axion KK-02, Kalb-Ramond + Green-Schwarz KK-03) are ALPs in which torsion plays no role. Warped compactifications: on S¹/Z₂ (RS included) g_μ5 is projected out and the warp factor is P'-even, so any photon is a bulk/brane U(1) in the KK-02/03 class (argument, not separately checked).

## Derivation summary (details in work/)
| step | result | where |
|---|---|---|
| Minimal EC + Maxwell | F = dA is connection independent; β = 0 | 003; skeptic 006 |
| EC + Dirac | S^μ = -(κ/2)J₅^μ, L_eff = -(3κ/16)J₅², no tree-level photon coupling | 003; skeptic 006 |
| CFJ rotation | β = +½∫p_τ dτ; -(g/4)φFF̃ gives β = -(g/2)(φ_0 - φ_LSS) in our conventions | 003; skeptic 006 |
| Spin source (EC-02) | ≤ 7.5×10⁻³⁰ deg at ξ = 1; 28.5 orders short | 007; skeptic 011 |
| Fermion loop (EC-05) | ≤ 1.3×10⁻³⁰ deg; ≥ 29 orders short | 007; skeptic 011 |
| Higher-dim. gauge-invariant torsion-photon operators with EC torsion | (∇·S)FF̃/Λ² needs Λ ≲ 2×10⁻⁴³ eV; DMP F∂F̃T needs Λ ≲ 4×10⁻¹⁵ eV | skeptic 020 S3 |
| Torsion pseudoscalars (EC-03) | ∂φ·J₅ shift symmetric: no QCD mass; photon coupling (m_φ/m_e)²/96 | 012; skeptics 011, 015 |
| Green-Schwarz completion | restores φFF̃ and φGG̃, so QCD mass and washout return | 012; skeptic 015 |
| PGT 0⁻ mode (EC-04) | shift-symmetric axion, g_φγ ≈ 0 | 007; skeptic 011 |
| Unsourced O(H) torsion (EC-06) | β = (ξε/2) ln(1+z*) for p ∝ H; ξε = 1.07×10⁻³; A_CB ≤ 5×10⁻¹¹ | 022, 023; skeptic 026 |
| MZ torsion (EC-06a) | GW170817 → ξ ≥ 1.4×10⁴ (tracker only); disfavoured | 022; skeptic 026 |
| f(T) trace torsion (EC-06c) | ξ = 3.6×10⁻⁴; c_T = 1; U(1) breaks on solutions; Stückelberg window 7×10⁻²⁸-10⁻²⁶ eV | 022, 023, 025; skeptic 026 |
| PGT 0⁺ SNY (EC-06d) | closed p, but needed histories have H(z*) = 33-66 × ΛCDM | 025 |
| Pure EC-KK (KK-01) | no parity-odd photon term at any order (P') | 004, 009; skeptics 006, 011 |
| Radion × FF̃ (KK-04) | Cassini γ-1 = -1/2, ISL m > 5 meV exclude | 009; skeptic 011 |
| u-Holst χFF̃ (KK-05) | Z = 1+12c², Y = 16c³, dβ = -c dlnα; no light KK-charged electron | 012, 016, 019; skeptics 015, 020 |

## Numbers
All numbers in CLAIMS.md are regenerated by `python -m src.reproduce` (expected values in `src/claims_expected.csv`). `--checks` also reruns the symbolic checks behind the no-number verdicts (7 scripts, ~16 min). EC-06c Stückelberg window, EFT cutoff and EC-06d H(z*) ratio are not yet in the table.

## Comparison with the axion-like baseline
Free-β ALP: Δχ² = 20.7 vs β = 0 (ACT + WMAP + Planck, EK22 combination). Penalties vs the free ALP:
| candidate | β_max | penalty (EK22) | penalty (Cosmoglobe) |
|---|---|---|---|
| EC-03' (cap 2.7% ρ_DM, the nearest tabulated mass 10⁻³² eV; best-fit mass 5.6×10⁻³³ eV has no published cap) | 0.100° | 8.0 | 4.95 |
| KK-02, k = 1 | 0.209° | 0.9 | 0.14 |
| EC-06c (free ξ) | unbounded | 0 (same Δχ² as free ALP) | 0 |
No candidate fits better than a generic ALP. SPT A_CB does not discriminate.

## Distinguishability
- ALP-class survivors (KK-02/03, EC-03', EC-04): not distinguishable from the axion-like explanation by β alone. KK-02 has β ≤ kα/2 ≈ 0.21° at k = 1, so a confirmed β well above that excludes it at k = 1.
- EC-06c (and the published torsion-free 2508.10532): β accrues evenly in ln(1+z), so 23-31% is seen by reionization-era EB. A single ALP mass can reproduce any one reionization/recombination ratio, so this tests the shape, not ALPs as a class; a LiteBIRD-class 0.05° measurement separates it from a late-rolling ALP at ~3σ (023).
- All survivors are achromatic; a confirmed chromatic β excludes them.

## Literature check (literature/novelty_audit.md)
- Already published: KR torsion parity conservation (hep-th/0207165); EFT no-go without gravity (2310.09152), which minimal EC falls under after torsion elimination; derivative-coupling photon decoupling for ALPs (1708.00443, 1903.12559, 2011.10025, 2112.00553) and shift symmetry for torsion scalars with gluons (2506.11836); the KK charge-to-mass problem (gr-qc/9805018); CFJ ∝ ∂ ln(geometric scalar) phenomenology (hep-ph/0611192, 2508.10532); torsion-photon fit with β ∝ ν² (0908.0629, disfavoured by Planck PR4 n at ~4.9σ).
- Applications not found in the literature: the P' symmetry argument applied to KK photon birefringence; the EC shortfall numbers (inside known no-gos); the photon decoupling applied to torsion pseudoscalars; the Cassini exclusion applied to 2608.09299; unsourced torsion as the origin of the CFJ vector (EC-06).

## Known weaknesses
- EC-02/EC-05 source estimates use generous order-of-magnitude assumptions; shortfalls are robust by many orders.
- KK: S¹ and S¹/Z₂ only; warped case covered by an unchecked argument.
- SPT A_CB template scaled by a fitted factor; ACT EB/TB and Planck map-level likelihoods not run locally (published summaries used).
- EC-03' cap taken at the nearest tabulated mass.
- EC-06c: Stückelberg lower edge assumes vacuum mixing (plasma would lower it); Proca + timelike CFJ microcausality not checked; EFT cutoff is power counting only; photon back-reaction on the f(T) frame not computed.

## Reproduce
`python -m src.reproduce` (about 3-4 minutes on 1 CPU; run `work/20261120-010-data-likelihood/download_candl_data.py` once first). Verified on a clean clone (021).
