# 039: skeptic on 036 (EC-06d heavy 0⁺) and 037 (shortfalls at the bound)

Iteration 12, 2026-11-21. Everything here is a candidate mechanism. There is no `RESULT.md`. The ACT β = 0.215° ± 0.074° is 2.9σ with systematics caveats.

## Question

- A (036): is the heavy-0⁺ branch of SNY PGT dead at every mass? Check the mass formula, the universal α ≥ 1/3 trace coupling, the Cassini kill, the ACT-vs-Planck/WMAP endpoint argument, the measured ratio, and the HB size shortfall.
- B (037): check the operator C J₅·AF̃, the dual-axion bound ξ < 1.5×10¹⁵, the CMS dijet kill of the vector route, and the Stückelberg weak corner. Compute that corner's strong-coupling scale.

## Method

Three scripts, all inputs from `data/` plus `inputs_039.csv` (grain size and density, rock density, collider energies; sources given, collider energies recalled and not opened). `run.log` has 22 PASS, 0 FAIL, about 2 min in total on 1 pinned CPU.

- `rc_trace_indep.py` (A1, sympy). The connection is solved from the 64 linear equations of metric compatibility plus torsion. It does not use 036's contortion formula. The curvature comes from an explicit [∇_m, ∇_n]V commutator, including the torsion term, so 036's Riemann index formula is not used either. The script then redoes the Einstein-frame mass and the Brans-Dicke α for general a0.
- `endpoint_grain.py` (A2-A4):
  - measured ratio, with a Monte Carlo;
  - the 0⁺ field at the emitting dust grain, using the exact interior Yukawa solution of a uniform sphere;
  - HB-bound size versus λ_C for a ground observer, with the rock half-space included (036 counted air only).
- `stueckelberg_cutoff.py` (B):
  - sympy gauge-variation identity;
  - explicit massless helicity spinors for f f̄ → γπ from L_π = -(C/m)π(∂_νJ_μ)F̃^{μν}, which equals -(1/4m)π(dp)·F̃ with p = -(ξκ/2)J₅;
  - J = 1 partial wave, unitarity |a_1| ≤ 1;
  - a 028/031-style exchange estimate, γγ → ππ through the one-loop ⟨J₅J₅⟩ correlator;
  - an NDA red-giant emission check.

## Results

| claim | status | reason |
|---|---|---|
| 036: m² = 2M_P²(a1+M_P²)/(a1 b) | SURVIVES | Independent connection solve plus commutator gives R = R_LC + 2∇·T - (2/3)T² in our conventions. The Einstein-frame m² = -2a0μ/(a1b) follows, and a1 → ∞ gives 2M_P²/b. |
| 036: universal trace coupling α = μ/(3a1) ≥ 1/3, unscreened | SURVIVES | ω = -3a0/(2a0-2a1), so α = μ/(3a1), the f(R) value 1/3 as a1 → ∞. K = -(b/4)/Y is constant for κρ ≪ Y_c, so there is no chameleon. |
| 036: Cassini kill for λ_C ≳ AU, 3.9 orders | SURVIVES | log₁₀(0.5/6.7×10⁻⁵) = 3.87. This is the standard massless-scalaron γ = 1/2 result. |
| 036: measured ratio β_ACT/β_EK22 = +0.63 ± 0.28 | SURVIVES | Reproduced from data/ (0.629 ± 0.275). MC gives P(ratio < 0) = 0.002. ACT/MK20 = 0.61. |
| 036: ACT includes the observer endpoint, MK/EK22 cancels it | SURVIVES | ACT's Π^PA α prior comes from the optics model (data row DPK25-alphaprior, 2403.00763). In MK/EK22 the observer rotation is common to the CMB and the foreground, so the fitted α_ν absorbs it. |
| 036: MK source endpoint is the diffuse-ISM mean, so the predicted ratio is negative for all λ_C < 20 pc | **BROKEN** (partly) | The calibrating photons are emitted inside dust grains, and the field at the emitter includes the grain's own coherent field, ρ̄ = ρ_g a²/(2λ²) for λ ≫ a. That exceeds ρ_* for λ_C < 2.1 km (a = 0.1 μm) to 5.2 km (a = 0.25 μm). There the predicted ratio is positive, ~ρ̄_obs/ρ̄_grain, and model dependent. The sign kill holds only for ~5 km ≲ λ_C < 20 pc. |
| 036: "every mass excluded" | SURVIVES (repaired) | The HB size kill, recomputed with the rock half-space, holds up to λ_C ≈ 2-3×10³ km: 6.7 orders short at 1 km, 4.7 at 10 km. That overlaps the corrected sign-kill range (≥ 5 km), which overlaps Cassini. Yukawa searches with α ≥ 1/3 also cover 100 μm to AU (not in data/). |
| 036: size 23.7 orders short at 100 μm | SURVIVES | 7.13×10⁻²⁷ rad reproduced. The formula \|β\| = √(2α) g Δρ̄/(4M_P m²) re-derived: β_c = √(α/2), φ = β_c ρ̄/(M_P m²), β = (g/2)Δφ. |
| 036: universal (MK) piece 41.4 orders short | **BROKEN** | The MK piece is set by the grain, not ρ_*. It is 5.7×10⁻²⁹ rad, **25.8 orders** short. The conclusion is unchanged. |
| 037: C = ξκ/4, dimension 6, δO = t.d. - α(∂_νJ_μ)F̃^{μν} | SURVIVES | sympy identity. The variation vanishes for curl-free J, which includes homogeneous FRW. |
| 037: dual-axion ξ < 1.5×10¹⁵ (Chandra × RG) | SURVIVES | 8×10⁻¹³ GeV⁻¹ × 1.6×10⁻¹³ / (κ m_e) = 1.49×10¹⁵. The operator fixes only the product g_aγ g_ae, and every split must obey both independent bounds, so the product bound holds. The axion is exactly massless, so Chandra (m_a < 10⁻¹² eV) applies. Heavy-fermion loops give no induced g_aγ at keV, because derivative couplings decouple. Caveats: Chandra's B-field model, and the 4m_eC mapping, which is checked only via Reig Eq. 11. |
| 037: CMS dijet kill of the vector route | SURVIVES | Arithmetic checked: coefficient of J_B² < 18π/Λ². In class D the J² contact term is an exact local term of the dual formulation, so a TeV contact bound applies. LL versus VV normalisation is O(1) against an 11-order margin. |
| 037: weak corner (Stückelberg, weakest m_γ, Λ_s ≥ 1 MeV) is 1.2 orders short | **BROKEN** | 037 used Λ_s = (m/C)^{1/3} = 0.39 MeV. The explicit J = 1 partial wave gives \|a_1\| = 4.69×10⁻³ (C/m) s^{3/2}, so Λ = 2.34 MeV, a factor 6.0 higher, and ξ_max at 1 MeV rises by 214. At 037's own criterion the Pauli envelope is therefore **reachable, 1.1 orders inside** (ξ_max = 4.0×10²⁶ vs needed 3.1×10²⁵). So 037's statement "no route reaches 0.215° in either class" is false for that criterion. |
| 037 weak corner closes once the cutoff is compared with tested energies | **closes** | At the needed ξ = 3.1×10²⁵ and m_γ = 7.9×10⁻¹⁴ eV: tree unitarity of e⁺e⁻ → γπ fails at **2.3 MeV** (1.9 MeV for \|a\| = ½). The 028/031-style γγ → ππ loop exchange (N = 4) gives 1.1 MeV. That is 3.7 orders below BaBar's e⁺e⁻ → γ + invisible √s = 10.58 GeV and 5.0 below LEP2. Requiring Λ ≥ 10.58 GeV gives ξ < 3.4×10¹⁴: **Pauli 11.0 orders short, EC-02 ν 13.9**. NDA red-giant emission (γe → πe vs the g_ae-allowed rate) independently needs ξ lower by 12.6 orders. |

Weak-corner table (tree, |a_1| = 1; full table in run.log):

| route | m_γ (eV) | ξ needed | Λ (MeV) |
|---|---|---|---|
| Pauli | 7.9×10⁻¹⁴ | 3.1×10²⁵ | 2.3 |
| Pauli | 10⁻¹⁸ (PDG) | 3.1×10²⁵ | 0.055 |
| EC-02 ν | 7.9×10⁻¹⁴ | 2.9×10²⁸ | 0.24 |
| EC-02 ν | 10⁻¹⁸ | 2.9×10²⁸ | 0.0056 |

## Verdict

- **036 survives** with two corrections. The EC-06d heavy-0⁺ branch stays dead at every mass. The endpoint-sign kill is valid only for λ_C ≳ 5 km, because dust grains dominate the MK source endpoint below that. The HB size kill covers λ_C ≲ 2×10³ km, and Cassini covers λ_C ≳ AU. The "41.4 orders" universal-piece number should read about 26 orders.
- **037 survives in its verdict, but the weak-corner number is wrong.** With a correctly normalised partial wave, the Pauli corner passes 037's 1 MeV criterion by 1.1 orders. The corner closes on the strong-coupling scale: 2.3 MeV at the needed ξ, 3.7 orders in energy below e⁺e⁻ → γ + invisible data at 10.58 GeV. Requiring Λ above that energy leaves the Pauli route 11 orders short. Proposed wording: "Stückelberg completion: strongly coupled at ≤ 2.3 MeV at the needed ξ even with the weakest m_γ; ≥ 11 orders short if the cutoff must exceed 10 GeV."
- The CLAIMS/STATE sentence "weakest corner leaves the Pauli route 1.2 orders short" should be replaced.

## Open issues

1. The grain endpoint assumes the photon polarization is set at the emitting dipole inside the grain. A near-field estimate (induced magnetic moment ~ gφ_grain d) supports this, but it is not a full radiative computation. It only moves the sign-kill boundary, which the size kill covers.
2. The point-particle self-field (Thomson electron at LSS, analyzer atoms) is argued negligible (≲ 10⁻³⁷ rad) and not tracked.
3. The size kill for λ_C between 10³ km and an AU relies on Yukawa/planetary data (not in data/) and on the sign kill. Both exclude it. The ground model is a flat rock half-space.
4. The partial-wave scale is single-species. Coupled channels over N light fermions lower it by N^{1/3}, which strengthens the conclusion.
5. The collider energies (BaBar 1702.03327, LEP2) are recalled and not opened. Positronium is not a clean test: S-wave e⁺e⁻ at rest gives zero for this operator (J₅^i vanishes in S-wave and F̃^{00} = 0).
6. The S5 star-emission estimate is NDA and ignores plasma and degeneracy factors. Its margin is 12.6 orders.
