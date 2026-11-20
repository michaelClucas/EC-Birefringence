# 20261120-003 theory-ec-photon (4D track, iteration 1)

The ACT DR6 β = 0.215° ± 0.074° (arXiv:2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Every mechanism here is a candidate mechanism, not an explanation.

## Question

1. Is the contortion formula in CONVENTIONS.md right?
2. Does minimal Einstein-Cartan with minimally coupled Maxwell give zero birefringence (including the fermion route at tree level)?
3. Which torsion components are allowed by homogeneity and isotropy?
4. Which couplings or extensions give nonzero β, with what formula and sign in our conventions, and what parameter size does β = 0.215° need?

## Method

Symbolic checks in sympy (exact rationals at random points, explicit Dirac matrices, Lie-derivative linear algebra) plus one numerical field evolution for the rotation sign. Full write-up: `derivation_ec_photon.md`. Run all checks (single core, about 2.5 min total, each script under 1 min):

```
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-003-theory-ec-photon/run_checks.py
```

Last run: `run.log`, 47 PASS, 0 FAIL.

| script | checks |
|---|---|
| check_contortion.py | contortion formula on random metric/torsion |
| check_maxwell_zero.py | F = dA connection-free; ∇-curl = F - T·A; gauge breaking; FRW conformal invariance; naive coupling on FRW torsion |
| check_dirac_torsion.py | spinor connection; c_A; torsion-squared part of R; elimination of algebraic S |
| check_frw_torsion.py | homogeneous isotropic torsion components |
| check_cfj_rotation.py | φFF̃ ↔ CFJ sign; EOM; dispersion; rotation sense (symbolic + numerical); sky mapping to β; gauge variation |
| check_nieh_yan.py | Nieh-Yan identity and θN = 6√(-g)S·∂θ |
| estimates.py | order-of-magnitude needed parameters (inputs: `data/birefringence_measurements.csv`, `inputs.json`) |

## Result

**Errors in CONVENTIONS.md (for the overseer to fix):**

1. Contortion. Correct for our index placement: **K_{λμν} = ½(T_{λμν} + T_{νλμ} - T_{μνλ})**. The printed ½(T_{λμν} - T_{μλν} - T_{νλμ}) equals K_{μνλ} of the correct one: it reproduces the torsion but is not metric compatible (it is the formula for a derivative-index-last Γ).
2. CFJ ↔ φFF̃. -(g/4)φFF̃ equals -½p_μA_νF̃^{μν} up to a total derivative for **p_μ = -g∂_μφ**, not +g∂_μφ.
3. β sign (to be cited in CONVENTIONS.md). With our ε_{0123} = +√(-g) and β > 0 clockwise on the sky (arXiv:2509.13654 as recorded in `literature/conventions_beta.md`): **β = +½∫p_τ dτ**, and L ⊃ -(g/4)φFF̃ gives **β = -(g/2)[φ_0 - φ_LSS]**. The paper's +(g/2)Δφ corresponds to Komatsu's ε^{0123} = +1/√(-g). Derived from the CFJ dispersion relation (ω_±² = k² ∓ p_τk for e_x ± ie_y), the rotation χ = -½∫p_τdτ right-handed about k (confirmed by a numerical evolution of the field equations: -0.25001 vs -0.25000 rad), and the observer-view mapping χ > 0 = counterclockwise seen. It agrees with the ε-sign conversion flagged by the scout.

**Zero for minimal EC (task 2), all PASS:**
- (a) F = dA and the Maxwell action contain no connection; on FRW the comoving Maxwell equations are the flat ones (conformal invariance), so β = 0.
- (b) Using the torsionful ∇ in F gives F - T^λ_{μν}A_λ, which breaks U(1) (δF' = -T^λ_{μν}∂_λχ). On FRW torsion it contains a CFJ term with p_μ = 2S_μ, a parity-even ⅓T_μA_νF^{μν} term, and photon mass terms ½(h² - S0²)A_i².
- (c) Dirac fermions: torsion enters only as c_A S_μJ_5^μ with c_A = 3/4; R(Γ) ⊃ -⅔T_μT^μ + (3/2)S_μS^μ; S is algebraic, S^μ = -(κ/2)J_5^μ, L_eff = -(3κ/16)J_5·J_5. No photon at tree level.

**FRW torsion (task 3), PASS:** exactly two components: T^i_{0j} = h(τ)δ^i_j (trace, T_τ = 3h) and T^i_{jk} = f(τ)[ijk] (axial, S^τ = -f/a², f pseudoscalar).

**Candidate mechanisms (task 4).** β in paper convention, S^0̂ = aS^τ the comoving-observer component, t cosmic time from last scattering. Details and actions: derivation Sec. 6.

| candidate mechanism | gauge invariant | effective term | β | needed (PRELIMINARY) |
|---|---|---|---|---|
| baseline φFF̃ | yes | -(g/4)φFF̃ | -(g/2)Δφ | \|gΔφ\| = 7.5×10⁻³ |
| (i) -½ξS_μA_νF̃^{μν} | only if dS = 0 (true on FRW background; anisotropy needs S = ∂ϑ) | CFJ, p = ξS | -(ξ/2)∫S^0̂dt | \|ξS^0̂\| = 1.1×10⁻³⁵ eV ≈ 8×10⁻³ H_0 (constant), or 2.8×10⁻⁴⁰ eV today if ∝ a⁻³ |
| (0') naive ∇-coupling | no | CFJ p = 2S + gauge-breaking terms + photon mass | -∫S^0̂dt | as (i) with ξ = 2 |
| (ii) Holst + BI field, (iii) Nieh-Yan × θ | yes | φFF̃ only if added by hand or loop-induced | -(g/2)Δφ | \|Δφ\| ≈ 1.8×10¹⁶ GeV if g ~ 1/M_P |
| (iv) fermion loop on S background | yes for closed S | CFJ, \|p\| = 2Ce²(3/4)\|S\| | ∓Ce²(3/4)∫S^0̂dt, C scheme dependent | \|S^0̂\| ≈ 4×10⁻³³ eV for one charged fermion, C = 3/16π² |
| (v) PGT axial torsion, Kalb-Ramond H = dB + CS, Itin-Hehl F²T² | see derivation | reduce to (i) or φFF̃ | as (i) / φFF̃ | not computed |

Minimal EC as the source of S: S^0̂ = -(κ/2)n_5. A relic-neutrino asymmetry |ξ_ν| = 0.1 (placeholder bound) gives |S^0̂| ≈ 2×10⁻⁶⁸ eV today; with a⁻³ scaling from last scattering (generous, ignores non-relativistic suppression) (i) needs |ξ| ≈ 10²⁸: **about 28 orders of magnitude short at ξ = 1**. Chiral fermionic dark matter is 29-30 orders short. Only a net chirality density sources the isotropic S^0̂; ordinary spin polarization averages to zero.

## Verdict

- Minimal EC with minimally coupled photon: β = 0 (settled, symbolic). The tree-level fermion route gives only the axial contact term.
- (i) and (iv) are kinematically able to produce β if a cosmological S^0̂ of order 10⁻³⁵ to 10⁻³³ eV exists. Minimal EC supplies about 10⁻⁶⁸ eV, so they need a torsion background not tied to spin density (e.g. parity-odd FRW torsion without spinors, Magueijo-Złośnik 1908.05184, or dynamical axial torsion). The fermion-loop coefficient C is regularization dependent (Jackiw-Kostelecky hep-ph/9901358), so (iv) is not predictive as it stands.
- (ii), (iii), (v-b) reduce to an axion-like φFF̃ model; torsion is a mediator, and the photon coupling is either put in by hand or comes from a loop whose size is disputed. If the coupling is universal to quarks, a QCD anomaly mass is generated (1502.03694), which points to washout. Candidate no-go for the minimal "torsion axion" route, not yet checked.
- None of the formulas here is new: CFJ rotation, the 3κ/16 contact term, the Nieh-Yan identity, and KR/BI axions are all in the prior work listed in `literature/prior_work_torsion_birefringence.md`. This task fixes them in our conventions and signs.

## Open issues

1. CONVENTIONS.md: fix the contortion formula and p_μ = -g∂_μφ, and record β = +½∫p_τdτ = -(g/2)Δφ for L ⊃ -(g/4)φFF̃.
2. Loop-induced photon coupling from ∂φJ_5/f or S_μJ_5^μ: decoupling argument ((ω/m_f)² suppression, e.g. "axions are blind to anomalies", arXiv ID to be confirmed by the literature scout) vs. anomaly coefficient used in 1502.03694. This decides (ii)-(iv).
3. A cosmological S^0̂ ~ 10⁻³⁵ eV source: check parity-odd FRW torsion without spin matter (1908.05184), and PGT 0⁻ mode dynamics (Shapiro hep-th/0103093 argues propagating axial torsion must be heavy).
4. A constant time-like CFJ vector has known stability/causality issues in the photon sector; a dynamical S = ∂ϑ avoids them. Not analysed here.
5. With algebraic S, adding (i) also gives contact terms κξJ_5·K and κξ²K·K with the gauge-variant Chern-Simons current K^μ = A_νF̃^{μν}. Not analysed.
6. Itin-Hehl F²T² route (θ ∝ ℓ²hS^0̂) needs both h and S^0̂; minimal EC with Dirac matter has h = 0. Normalization not derived.
7. inputs.json placeholders (neutrino degeneracy bound, light fermion DM mass) need cited values from the bounds analyst; cosmology/constants should move to `data/`.
8. Frequency dependence: all candidates here are achromatic at leading order (rotation independent of k for k ≫ |p|). Das-Mohanty-Prasanna-type derivative couplings give ν² scaling instead; not treated.
