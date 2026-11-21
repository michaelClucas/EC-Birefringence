# RESULT draft, revision 4 (not RESULT.md)

Revised by the overseer in iteration 11 after final skeptic 035 and skeptic 034. Becomes `RESULT.md` only after a fresh skeptic fails to break it. Iteration-11 tasks 036-038 folded in; 036 and 037 await skeptic review.

## Question
Can Einstein-Cartan (ECSK) torsion, in 4D or in an Einstein-Cartan extension of Kaluza-Klein theory, explain the cosmic birefringence reported in ACT DR6 (β = 0.215° ± 0.074°, arXiv:2509.13654) and constrained by SPT-3G (A_CB < 1.2×10⁻⁴, arXiv:2510.07928), in the context of arXiv:2601.13624 and 2512.21888?

Caveat on the target: the ACT β is 2.9σ (4.6σ combined with WMAP and Planck; 3.9σ with the Cosmoglobe reanalysis), with the miscalibration angle set by an optics-model prior and dust EB neglected. ACT per-band angles prefer a power law n ≈ 1.6 over achromatic at 2.0σ (driven by the PA5 f150-f090 split, attributed to instrument systematics); Planck PR4 gives n = -0.35 ± 0.48. The signal may not be cosmological.

## Answer
- **4D track: Result C (mixed).** Core finding (no-go): no photon coupling that Einstein-Cartan torsion generates by itself, minimally or through fermions, produces the observed β; every such route falls short by many orders of magnitude or loses its photon coupling to an exact shift symmetry, and every unsourced-torsion route examined (EC with varying Λ, TEGR, f(T)) either is undefined or breaks U(1) gauge invariance and becomes strongly coupled far below tested energies. β of the observed size is possible only by adding, by hand, a photon coupling to a torsion-associated field that then behaves as an axion-like particle (ALP); the costs are listed below. The one gauge-invariant unsourced-torsion route (PGT 0⁺, closed trace torsion) is excluded at every mass of the 0⁺ mode.
- **Einstein-Cartan-Kaluza-Klein track: Result C (mixed).** Core finding (no-go): pure EC-KK has no parity-odd photon coupling at any order in a P'-invariant state; radion × FF̃ is excluded; the one route where torsion generates the coupling (u-Holst) has no light charged matter. β is possible only through ALPs with couplings added beyond EC-KK (5D Chern-Simons Wilson line; Kalb-Ramond + Green-Schwarz), in which torsion plays no role or is a relabelling of the 2-form.

## Routes that can produce β, with their costs (Result C content)
| id | route | added beyond EC | coupling needed | other costs | fit vs free ALP |
|---|---|---|---|---|---|
| EC-04 | PGT dynamical 0⁻ torsion mode + CFJ term | PGT (propagating torsion), photon coupling by hand | 10⁻¹⁴ ≲ ξ < 2×10⁸ | m ~ H₀ needs ~4×10¹²⁰ tuning of the mass term | free amplitude: penalty 0 |
| EC-03' | pseudoscalar with torsion-scale f, coupling by hand | photon coupling by hand (torsion generates none) | g = 8α/(πf) = 4.7×10⁻²¹ GeV⁻¹ | energy fraction 11-19% ρ_DM needed; cited cap 2.7-6.9% at 10⁻³² eV (no published cap at its best-fit mass 5.6×10⁻³³ eV) | penalty 8.0 (EK22), 4.95 (Cosmoglobe) at 2.7% |
| EC-06a remnant | f(Λ)FF̃ or g = ξ/(√3 M_P) coupling to the MZ Λ field | MZ extension + ALP-type coupling (not the torsion CFJ) | ξ 21× the CFJ value | ALP-class; viability not computed | not fit |
| KK-02 | bulk U(1) + 5D CS, Wilson-line axion | bulk U(1), CS term, charged bulk fermion for the potential | k = 1 (allowed with the mixed CS term) | 1/L between 5.5×10⁹ GeV and ~10¹⁸ GeV; 2πLM ≈ 220-280 | β ≤ 0.209°; penalty 0.9 (EK22), 0.14 (Cosmoglobe) |
| KK-03 | Kalb-Ramond 2-form (= totally antisymmetric torsion) + Green-Schwarz term | GS term; natural λ = κ₄ excluded | g between 1.4×10⁻¹³ and 8×10⁻¹³ GeV⁻¹ | narrow window between BBN and Chandra | free amplitude in window: penalty 0 |

## Kills with shortfalls (4D)
| route | reason | shortfall | where |
|---|---|---|---|
| Minimal EC + Maxwell (EC-00) | F = dA connection independent | exact zero | 003; skeptic 006 |
| Spin-sourced torsion + CFJ (EC-02) | relic ν asymmetry (ξ_ν ≤ 0.048, 2409.08280) | 28.5 orders at ξ = 1; at the largest allowed ξ: 10.5-13.3 orders (dual-axion completion: CAST / Chandra × red giant, ξ < 9.4×10¹⁷ / 1.5×10¹⁵) or 4.2-9.1 orders (Stückelberg completion: weakest / PDG photon mass) | 007, 037; skeptic 011 |
| Any fermionic torsion source (Pauli bound) | n_5 ≤ fermion number density | ≥ 25 orders (25.49) at ξ = 1; at the bound 7.5-10.3 (dual axion) or 1.2-6.1 orders (Stückelberg) | skeptic 029, 037 |
| Vector-current-sourced torsion (non-minimal) | baryon density; the signal fixes Λ*² ≈ 7 GeV², and the induced contact term exceeds the CMS dijet limit for any Immirzi parameter | ≥ 11.1 (dual axion) or 12.0-16.9 orders (Stückelberg) | skeptic 029, 037 |

Coupling-bound caveat (037): the operator ξS·AF̃ with EC torsion is C J₅·AF̃ (dimension 6), not gauge invariant off FRW, so it needs a completion; every killing bound at large ξ depends on that completion (dual massless axion, Reig 2611.07578, or a Stückelberg photon). The weakest corner, a Stückelberg photon at the weakest photon-mass bound with a 1 MeV cutoff (a choice), leaves the Pauli-bound route 1.2 orders short; that corner ignores the Proca galactic limit (+8 orders) and a dimensional-analysis ΔN_eff thermalisation estimate (+15 orders), and the Stückelberg mode's own strong coupling (as for EC-06a/EC-06c) is not computed for it.
| Fermion-loop CFJ (EC-05) | loop factor × EC torsion | ≥ 29 orders | 007; skeptic 011 |
| Higher-dim. gauge-invariant torsion-photon operators with EC torsion | (∇·S)FF̃/Λ², F∂F̃T | Λ ≲ 2×10⁻⁴³ eV / 4×10⁻¹⁵ eV needed | skeptic 020 (estimate, not in reproduce) |
| Torsion pseudoscalars via ∂φ·J₅ (EC-03) | exact shift symmetry; photon coupling (m_φ/m_e)²/96 | 10⁻⁷⁸-10⁻⁷⁰ suppression | 012; skeptics 011, 015 |
| Green-Schwarz completion of EC-03 | restores φGG̃ and the QCD mass, so washout | superseded EC-03 numbers | 012 |
| MZ torsion, EC + varying Λ (EC-06a) | CFJ vector not closed in matter; Stückelberg photon; strong coupling at 3.6 keV | 9.4 orders in energy below 10 TeV; off-tracker best case ≤ 9 MeV | 031; skeptic 034 |
| TEGR torsion (EC-06b) | not defined by TEGR's dynamics | n/a | 022; skeptic 026 |
| f(T) trace torsion (EC-06c, not EC) | U(1) broken on solutions; Stückelberg repair strongly coupled at 0.2-0.7 keV (7 MeV at the PDG limit) | ≥ 6 orders in energy | 025; skeptic 028 |
| PGT 0⁺ SNY, O(H) branch (EC-06d) | H(z*) = 29-57 × ΛCDM vs δH/H < 0.18-0.49% (95%) | 3.8-4.5 orders | 025, 032; skeptic 028 |
| PGT 0⁺, heavy GR-compatible branch (EC-06d) | 0⁺ mode couples to the matter trace with α ≥ 1/3, unscreened: Cassini excludes Compton lengths ≳ AU; shorter lengths make the observer endpoint local, so ACT (optics calibration) and Planck/WMAP (dust-EB calibration) would see opposite signs | 3.9 orders vs Cassini; ACT/Planck sign; 23.7 orders short in size at the HB bound | 036 (skeptic pending) |

The strong-coupling comparisons use 10 TeV as the scale up to which photon physics is observed to be standard (`E_collider_safe` in data/ec06c_costs_inputs.csv, a choice representing LHC and TeV-blazar photon physics, not a single cited measurement). The verdicts do not depend on its exact value: the strong-coupling scales are keV-MeV.

## Kills with shortfalls (KK)
| route | reason | shortfall | where |
|---|---|---|---|
| Pure EC-KK (KK-01) | no parity-odd photon term at any order in a P'-invariant state (P' = diag(1,-1,-1,-1,-1)); P'-odd matter backgrounds not computed in 5D (in 4D their analogues are Planck suppressed) | exact zero | 004, 009; skeptics 006, 011 |
| Radion × FF̃ (KK-04, incl. 2608.09299) | radion must roll after recombination; matter couples to the Jordan metric | ~4.3 orders vs Cassini γ (4.34; 3.6-4.1 for brane matter or a 2σ bound) | 009, 032; skeptic 011 |
| u-Holst χFF̃ (KK-05) | CMB photon must be the KK photon; no light state carries KK charge e | ~21 orders in the electron mass (m ≳ eM_P/√2) | 012, 016, 019; skeptics 015, 020 (not in reproduce) |
| Warped (S¹/Z₂, RS) | g_μ5 projected out; warp factor P'-even; any photon is a bulk/brane U(1) in the KK-02/03 class | argument, not separately checked | 004, 009 |

## Not examined
- Metric-affine / non-metricity photon sectors (e.g. 2607.26141 in literature/).
- Weyl-invariant Einstein-Cartan with an ALP (2406.11956, noted, not checked).
- EC with a non-minimally coupled scalar (torsion ∝ ∂ ln Ω(φ), the Higgs-in-EC family); same closed-CFJ class as EC-06d heavy.
- The Hojman-Rosenbaum-Ryan gradient-torsion scheme (in literature/, not mapped onto EC-04/EC-06).
- Bosonic spin sources of torsion.
- More than five dimensions.
- KK-05 loophole: a brane photon mixing with the KK photon.

## Comparison with the axion-like baseline
Free-β ALP: Δχ² = 20.7 vs β = 0 (ACT + WMAP + Planck, EK22). Penalties per route are in the cost table above. No route fits better than a generic ALP; routes with a free amplitude (EC-04, KK-03) tie it. SPT A_CB does not discriminate.

## Distinguishability
- All surviving routes are ALP-like: not distinguishable from the axion-like explanation by β alone. KK-02 has β ≤ kα/2 ≈ 0.21° at k = 1.
- All are achromatic; a confirmed chromatic β excludes them.
- The O(H)-torsion phenomenology (β even in ln(1+z), 23-31% after reionization) belongs to dead torsion routes and to the published torsion-free model 2508.10532; a LiteBIRD-class 0.05° measurement separates that shape from a late-rolling ALP at ~3σ (023, not in reproduce).

## Literature check (literature/novelty_audit.md)
- Already published: KR torsion parity conservation (hep-th/0207165); EFT no-go without gravity (2310.09152), which minimal EC falls under after torsion elimination; derivative-coupling photon decoupling for ALPs (1708.00443, 1903.12559, 2011.10025, 2112.00553) and shift symmetry for torsion scalars with gluons (2506.11836); KK charge-to-mass problem (gr-qc/9805018); CFJ ∝ ∂ ln(geometric scalar) phenomenology (hep-ph/0611192, 2508.10532); torsion-photon fit with β ∝ ν² (0908.0629, disfavoured by Planck PR4 at ~4.9σ, not in reproduce); neutrino-current EB shapes (2610.02667).
- Applications not found in the literature: the P' argument for KK photon birefringence; the EC shortfall numbers (inside known no-gos); photon decoupling applied to torsion pseudoscalars; the Cassini exclusion applied to 2608.09299; strong-coupling kills of unsourced-torsion CFJ routes (f(T), MZ).

## Known weaknesses
- Order-of-magnitude source estimates for EC-02/EC-05.
- Strong-coupling scales (EC-06a, EC-06c) are power counting on the scalar sector.
- KK: S¹ and S¹/Z₂ only; warped case by argument; P'-odd 5D matter backgrounds not computed.
- SPT A_CB template scaled by a fitted factor; ACT EB/TB and Planck map-level likelihoods not run locally (published summaries used).
- Unresolved factor 4.47 in the Data Tables β → k_AF map (literature/bounds.md); no verdict depends on it.
- Correction log: iteration 6 called EC-06a dead on a mis-sided GW bound; iteration 7 softened it; iteration 9 killed it by strong coupling (031, skeptic 034).

## Reproduce
`python -m src.reproduce` (about 4-7 minutes on 1 CPU; run `work/20261120-010-data-likelihood/download_candl_data.py` once first) regenerates every number in CLAIMS.md and the draft-only numbers in `src/draft_expected.csv`. `--checks` reruns the symbolic checks. Verified on a clean clone (021, 030). Numbers marked "not in reproduce" above are cited or estimated.
