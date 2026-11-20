# 20261120-009-kk-bounds: KK candidate mechanisms against bounds, and arXiv:2608.09299

## Question
1. How does 2608.09299 (5D teleparallel KK) get a radion × FF̃ coupling when the EC-KK reduction of task 004 finds none? Is there an Einstein-Cartan analogue? Does the mechanism change α, and which bounds do they apply? This includes the skeptic's input that a proper 5D rotation P' forbids radion × F_KK F̃_KK at all orders.
2. KK-02 (Wilson-line axion from a bulk U(1) with a 5D CS term): which (L, κ_CS) give β ≈ 0.2°, and what do the bounds say?
3. KK-03 (KR 2-form = axial torsion, with a Green-Schwarz term): is λ free, and what is the viable window?
4. S¹/Z₂ with a brane photon and a bulk axion.
5. Verdict per KK claim.

## Method
- Read 2608.09299 v3 in full (arXiv HTML). Notes in `literature/2608.09299.md`.
- Derivations are in `derivation.md`. Every new step has a scripted check:
  - `check_parity.py`: P' as a proper rotation, transformation of the zero modes, the 2608 EM term = n∧dn∧dn, and the orientation of u as the P'-odd datum. Exact rational arithmetic on random polynomial fields.
  - `check_ec_uholst.py`: EC analogues. Pointwise 5D frame algebra, with exact integration over the 50 contortion components.
  - `check_orbifold.py`: CS densities under y → -y, and brane gauge invariance.
  - `check_rotation_hosotani.py`:
    - helicity evolution for L = -¼Z(F² + c_odd FF̃);
    - Hosotani potential, closed form against a direct KK-mode proper-time sum.
  - `kk_bounds.py`: all numbers, written to `results.csv`. Measured inputs come from `data/bounds.csv`, `data/birefringence_measurements.csv`, `data/cosmo_torsion_axion_inputs.csv` and the new `data/kk_bounds_inputs.csv`, which holds Cassini γ (1403.7377 Sec. 4.1.2 and Table 4), BBN N_ν (1912.01132 abstract) and Planck+BAO N_eff (1807.06209 abstract).
- Rerun from the repo root: `.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-009-kk-bounds/run_checks.py`. It takes about 5 min in total; each script runs under 3 min on one core. Output goes to `run.log`.

## Result
**Checks: 33/33 PASS** (`run.log`).

1. **The skeptic is right.** P' = diag(1,-1,-1,-1,-1) is a proper rotation, so it is a symmetry of every diffeomorphism-invariant 5D theory. Radion × F_KK F̃_KK is P'-odd, so it is forbidden at all loop orders in any P'-invariant vacuum. This extends the tree-level KK-01 no-go of task 004. A parity-odd photon term needs a P'-odd field or a P'-odd background.
2. **What 2608 does.** Their EM term is exactly the 5D abelian CS form of the KK 1-form: n∧dn∧dn = 2φ³dy∧dA∧dA, with n = e⁵ = φ(dy + √2A). The Nieh-Yan part carries no photon.
   - The term is built from the metric and the circle direction, not from torsion.
   - It is nonzero only because u (four Stückelberg scalars with no kinetic term) is set to +n "with a consistent orientation". That orientation, the sign of det ∂_μψ^a, is the P'-odd background. It is neither a boundary term nor a non-Lagrangian term.
   - The radion moves α: α ∝ φ⁻³, so |Δα/α|_{LSS→0} = |β|/|c|.
   - They apply no α, fifth-force, solar-system or stabilization bound. Stabilization is explicitly out of scope; only the Planck PR4 β = 0.30 ± 0.11° fit appears.
3. **EC analogue.**
   - 2608's term transplanted to EC (u∧T^A∧T_A with spin-sourced torsion) gives nothing for the photon.
   - u∧Nieh-Yan gives only an even F² term.
   - The "u-Holst" term c u∧e^A∧e^B∧R_AB(ω) vanishes identically in Riemannian geometry (2608's P̂₁), but in EC it gives, after the contortion is integrated out exactly, L = -(1/8κ₅)(1 + 12c²)[W² + c_odd WW̃] with c_odd = 16c³/(1 + 12c²) and W = dn. So a torsion-mediated radion × FF̃ exists. It needs the same oriented u.
4. **The whole radion × FF̃ class is excluded** (2608, and the EC u-Holst variant: candidate mechanism, not in CLAIMS yet; suggested id KK-04).
   - For L = -¼Z(F² + c_odd FF̃): |β| = (|c_odd|/2)|Δln Z| and |Δα/α| = 2|β|/|c_odd|.
   - α-variation alone requires c_odd ≳ 0.7 (CMB α), and c_odd ≳ 4 (clocks) if the radion is still rolling today. These are not fatal.
   - The fatal bound: the radion must roll after recombination (m ≲ H_* ≈ 3×10⁻²⁹ eV), but a light KK radion is Brans-Dicke with ω = 0, so γ - 1 = -1/2 against Cassini (2.1 ± 2.3)×10⁻⁵. ISL requires m > 5×10⁻³ eV, about 10²⁶ H_*, for gravitational-strength couplings.
   - Torsion does not change this.
5. **KK-02.**
   - CS level quantization (κ_CS = k/(24π²)) gives g = kα/(πf_a) and **β = kαΔϑ/(2π), independent of L and f_a**, with |β| ≤ kα/2 = 0.209° k for Δϑ ≤ π.
   - ACT 0.215 ± 0.074° needs Δϑ = 1.03π for k = 1 (a hilltop start; 0.67π to 1.38π at 1σ), 0.51π for k = 2, and 0.34π for k = 3.
   - L window: 5.5×10⁹ k GeV < 1/L < ~10¹⁸ GeV (Chandra g bound; 5D EFT). Collider KK bounds are irrelevant.
   - Mass: a charged bulk fermion with Ml ≈ 220-280 (l = 2πL) puts the Hosotani mass between H₀ and H_*, so the field rolls after recombination. With no charged bulk matter the field is massless and needs a pre-BBN kick, with f_a < 1.7×10¹⁰ GeV for k = 1 (from the stiff-energy bound).
   - α and radion: β does not care about the radion, which can be stabilized (m > 5 meV) with no α variation.
   - No QCD mass. Anisotropy is below the SPT bound for k ≤ 3.
   - Washout: the mass range is given above; task 008 does the calculation.
6. **KK-03.**
   - λ is free. Torsion has no flux quantization, so B and its dual b are non-compact. A compact B would quantize λ and turn KK-03 into KK-02.
   - b is massless and moves only after a kick, as a stiff field. BBN ΔN_ν < 0.16 then requires g > 1.4×10⁻¹³ GeV⁻¹; Chandra gives g < 8×10⁻¹³ GeV⁻¹.
   - Window: λ/κ₄ ∈ [1.5, 8.7]×10⁶, i.e. a GS scale M_s ≈ (0.8-2)×10¹⁵ GeV.
   - The natural value λ = κ₄ is excluded (ΔN_ν ~ 10¹¹).
   - The only torsion-specific signal, S⁰ ~ 10⁻⁴⁵ eV today, is about 20 orders below lab bounds.
7. **S¹/Z₂.**
   - A brane-localized photon cannot couple to the A'_5 Wilson line of a bulk U(1)' at tree level. A'_5(0)FF̃ is not gauge invariant, brane fields cannot carry U(1)' charge, and derivative couplings decouple.
   - It works only with bulk hypercharge, through the mixed 𝒜'∧𝓕_Y∧𝓕_Y term, which is orbifold-allowed (that is KK-02 again).
   - A bulk B_{μν} (even, consistent with T = H) plus a brane Green-Schwarz term (Hořava-Witten type) gives a KK-03-type coupling with a free brane coefficient.

## Verdict (all are candidate mechanisms; none is an explanation)
| claim | verdict | cost | does torsion matter? |
|---|---|---|---|
| KK-01 pure EC-KK | **dead**, now at all orders for P'-even scalars × F_KK F̃_KK | n/a | n/a |
| (suggested KK-04) radion × FF̃ via oriented u: 2608 teleparallel or EC u-Holst | **dead**: Cassini/ISL forbid the light radion needed to roll after recombination | P'-odd Stückelberg orientation, plus an unscreened ultralight radion | 2608: no (metric-only n∧dn∧dn). EC u-Holst: yes (zero without torsion), but it dies anyway |
| KK-02 bulk U(1) + CS, Wilson-line axion | **live** | bulk hypercharge, CS level k; a charged bulk fermion with Ml ≈ 220-280 (or a pre-BBN kick); k = 1 needs a hilltop start, k ≥ 2 is generic | **no: ALP in disguise** |
| KK-03 KR = torsion + GS | **live, narrow** | λ ~ 10⁶κ₄ (M_s ~ 10¹⁵ GeV), massless b with a pre-BBN kick, g within a factor 6 of the Chandra bound | **no**: torsion only fixes the kinetic normalization, which is degenerate with λ. ALP in disguise |

For the program's question, which requires torsion to matter: **neither live KK claim counts as an Einstein-Cartan-KK explanation.** The only torsion-essential KK coupling found (EC u-Holst) belongs to the excluded radion class. KK-02's β ≈ kα/2 × (Δϑ/π) prediction (0.209° for k = 1) is a property of any quantized axion-photon coupling, not of EC.

## Open issues
1. CS level normalization: k ∈ Z on spin manifolds with the mixed gauge-gravitational term, or 6Z without it. This is quoted, not derived here. If only k ∈ 6Z is allowed, β_max = 1.25° per unit and k = 6 needs Δϑ ≈ 0.17π.
2. The KK-03 and massless KK-02 bounds assume the kick happens before BBN and that the stiff energy scales as a⁻⁶ back to T = 1 MeV. A kick after BBN (no known mechanism) would relax them to the CMB ΔN_eff bound (Δb < 4×10¹⁷ GeV).
3. The sign of β relative to c_odd and k is convention-level and unchecked. Magnitudes only.
4. EC u-Holst: the GW sector and the fermion cross term (an EDM-like W·J) were not computed. Radion screening (chameleon type), which could rescue the radion class, was not explored.
5. The S¹/Z₂ brane GS coupling for the KR axion is argued, not symbolically checked.
6. Washout and the precise mass window (H₀ to H_* is order-of-magnitude) are left to task 008.
7. `data/kk_bounds_inputs.csv`: the Cassini number is taken from Will's review (1403.7377), not the original Nature paper (no arXiv). BBN N_ν comes from the abstract only.
8. KK-02's 1/L upper limit uses M_5³(2πL) = M_P² with the reduced Planck mass. O(1) factors are unchecked.
