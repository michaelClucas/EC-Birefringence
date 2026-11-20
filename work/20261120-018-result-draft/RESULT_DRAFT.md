# RESULT draft, revision 3 (not RESULT.md)

Revised by the overseer in iteration 9 after final skeptic 029 and tasks 028, 030-033. Becomes `RESULT.md` only after a fresh skeptic fails to break it. Pending: skeptic review of 031 (EC-06a).

## Question
Can Einstein-Cartan (ECSK) torsion, in 4D or in an Einstein-Cartan extension of Kaluza-Klein theory, explain the cosmic birefringence reported in ACT DR6 (β = 0.215° ± 0.074°, arXiv:2509.13654) and constrained by SPT-3G (A_CB < 1.2×10⁻⁴, arXiv:2510.07928), in the context of arXiv:2601.13624 and 2512.21888?

Caveat on the target: the ACT β is 2.9σ (4.6σ combined with WMAP and Planck; 3.9σ with the Cosmoglobe reanalysis), with the miscalibration angle set by an optics-model prior and dust EB neglected. ACT per-band angles prefer a power law n ≈ 1.6 over achromatic at 2.0σ (driven by the PA5 f150-f090 split, attributed to instrument systematics); Planck PR4 gives n = -0.35 ± 0.48. The signal may not be cosmological.

## Answer
- **4D track: Result B (no-go) for torsion.** No torsion mechanism we examined, in Einstein-Cartan or in the named extensions, produces the observed β within independent bounds and with a consistent photon sector. Birefringence of the observed size is possible only by adding an axion-like field with a photon coupling put in by hand beyond EC (EC-03', EC-04); such a field reproduces β as any ALP does, and torsion at most sets its decay constant.
- **Einstein-Cartan-Kaluza-Klein track: Result B (no-go) for torsion.** Same structure: the surviving routes (5D Chern-Simons Wilson-line axion KK-02; Kalb-Ramond + Green-Schwarz KK-03) are ALPs with couplings added beyond EC-KK, and torsion plays no role in them.

## Kills with shortfalls (4D)
| route | reason | shortfall | where |
|---|---|---|---|
| Minimal EC + Maxwell (EC-00) | F = dA connection independent; β = 0 | exact zero | 003; skeptic 006 |
| Spin-sourced torsion + CFJ (EC-02) | relic ν asymmetry (ξ_ν ≤ 0.048, 2409.08280) | 28.5 orders at ξ = 1 | 007; skeptic 011 |
| Any fermionic torsion source (Pauli bound, condensates) | n_5 ≤ fermion number density | ≥ 25 orders (25.49) | skeptic 029 |
| Vector-current-sourced torsion (non-minimal) | baryon density | c ≈ 10³⁶ needed (36 orders) | skeptic 029 |
| Fermion-loop CFJ (EC-05) | loop factor × EC torsion | ≥ 29 orders | 007; skeptic 011 |
| Higher-dim. gauge-invariant torsion-photon operators with EC torsion | (∇·S)FF̃/Λ², DMP F∂F̃T | Λ ≲ 2×10⁻⁴³ eV / 4×10⁻¹⁵ eV needed | skeptic 020 |
| Torsion pseudoscalars via ∂φ·J₅ (EC-03: KR torsion axion, BI × Nieh-Yan, Holst-BI with fermions) | exact shift symmetry: no QCD mass; photon coupling (m_φ/m_e)²/96 | 10⁻⁷⁸-10⁻⁷⁰ suppression | 012; skeptics 011, 015 |
| Green-Schwarz completion of EC-03 | restores φGG̃ and the QCD mass, so washout | as superseded EC-03 numbers | 012 |
| MZ torsion, EC + varying Λ (EC-06a) | CFJ vector not closed in matter, Stückelberg photon needed; strong coupling at 3.6 keV | 9.4 orders in energy below 10 TeV | 031 (skeptic pending) |
| TEGR torsion (EC-06b) | not defined by TEGR's dynamics (frame choice) | n/a | 022; skeptic 026 |
| f(T) trace torsion (EC-06c, not EC) | U(1) broken on solutions; Stückelberg repair strongly coupled at 0.2-0.7 keV (7 MeV at the PDG limit) | ≥ 6 orders in energy below where photon physics is tested | 025; skeptic 028 |
| PGT 0⁺ torsion, SNY (EC-06d) | needed histories have H(z*) = 29-57 × ΛCDM vs δH/H < 0.18-0.49% (95%) | 3.8-4.5 orders | 025, 032; skeptic 028 |
| PGT 0⁻ mode (EC-04) | shift-symmetric axion, g_φγ ≈ 0 | ALP only with a coupling added by hand; m ~ H₀ needs ~4×10¹²⁰ tuning | 007; skeptic 011 |

## Kills with shortfalls (KK)
| route | reason | shortfall | where |
|---|---|---|---|
| Pure EC-KK (KK-01) | no parity-odd photon term at any order in a P'-invariant state (P' = diag(1,-1,-1,-1,-1)); P'-odd matter backgrounds are Planck suppressed as in 4D | exact zero | 004, 009; skeptics 006, 011 |
| Radion × FF̃ (KK-04, incl. 2608.09299) | radion must roll after recombination; matter couples to the Jordan metric | ~4.3 orders vs Cassini γ (4.34; 3.6-4.1 brane/2σ) | 009, 032; skeptic 011 |
| u-Holst χFF̃ (KK-05) | CMB photon must be the KK photon; no light state carries KK charge e (m ≳ eM_P/√2 ~10²¹ m_e) | ~21 orders in the electron mass | 012, 016, 019; skeptics 015, 020 |
| Warped (S¹/Z₂, RS) | g_μ5 projected out, warp factor P'-even, any photon is a bulk/brane U(1) in the KK-02/03 class | argument, not separately checked | 004, 009; skeptic 020 |

## Comparison with the axion-like baseline
Free-β ALP: Δχ² = 20.7 vs β = 0 (ACT + WMAP + Planck, EK22 combination). Penalties of the surviving ALP-class routes vs the free ALP:
| candidate | β_max | penalty (EK22) | penalty (Cosmoglobe) |
|---|---|---|---|
| EC-03' (cap 2.7% ρ_DM at the nearest tabulated mass 10⁻³² eV; no published cap at its best-fit mass 5.6×10⁻³³ eV) | 0.100° | 8.0 | 4.95 |
| KK-02, k = 1 | 0.209° | 0.9 | 0.14 |
No candidate fits better than a generic ALP. SPT A_CB does not discriminate.

## Distinguishability
- ALP-class survivors are not distinguishable from the axion-like explanation by β alone. KK-02 has β ≤ kα/2 ≈ 0.21° at k = 1, so a confirmed β well above that excludes it at k = 1.
- All survivors are achromatic; a confirmed chromatic β excludes them.
- The dead O(H)-torsion routes would have accrued β evenly in ln(1+z) (23-31% after reionization); the same phenomenology exists without torsion (2508.10532) and is separable from a late-rolling ALP at ~3σ by a LiteBIRD-class 0.05° measurement (023). This is a statement about the published torsion-free model, not a torsion prediction.

## Literature check (literature/novelty_audit.md)
- Already published: KR torsion parity conservation (hep-th/0207165); EFT no-go without gravity (2310.09152), which minimal EC falls under after torsion elimination; derivative-coupling photon decoupling for ALPs (1708.00443, 1903.12559, 2011.10025, 2112.00553) and shift symmetry for torsion scalars with gluons (2506.11836); KK charge-to-mass problem (gr-qc/9805018); CFJ ∝ ∂ ln(geometric scalar) phenomenology (hep-ph/0611192, 2508.10532); torsion-photon fit with β ∝ ν² (0908.0629, disfavoured by Planck PR4 at ~4.9σ); neutrino-current EB shapes (2610.02667), consistent with EC-02's shortfall.
- Applications not found in the literature: the P' argument for KK photon birefringence; the EC shortfall numbers (inside known no-gos); photon decoupling applied to torsion pseudoscalars; the Cassini exclusion applied to 2608.09299; strong-coupling kills of unsourced-torsion CFJ routes (f(T), MZ).

## Known weaknesses
- Order-of-magnitude source estimates for EC-02/EC-05 (shortfalls robust by many orders).
- Strong-coupling scales (EC-06a, EC-06c) are power counting on the scalar sector; the true scales may be lower.
- KK: S¹ and S¹/Z₂ only; warped case by argument.
- SPT A_CB template scaled by a fitted factor; ACT EB/TB and Planck map-level likelihoods not run locally (published summaries used).
- EC-03' cap taken at the nearest tabulated mass.
- Not regenerated by reproduce: 020's higher-dim. Λ scales, the 0908.0629 4.9σ, the LiteBIRD 3σ, KK-05's 10²¹ m_e.

## Reproduce
`python -m src.reproduce` (about 4-7 minutes on 1 CPU; run `work/20261120-010-data-likelihood/download_candl_data.py` once first) regenerates every number in CLAIMS.md (66) and the draft-only numbers in `src/draft_expected.csv`. `--checks` also reruns the symbolic checks (8 scripts). Verified on a clean clone (021, 030).
