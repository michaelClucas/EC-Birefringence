# EC-06 derivation: unsourced O(H) torsion backgrounds + a CFJ photon coupling

Conventions: `CONVENTIONS.md`. Every candidate here is a candidate mechanism. The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with systematics caveats and may not be cosmological.

## 0. Common formulas

Photon term (added by hand, not part of any theory below): L ⊃ -½ ξ X_μ A_ν F̃^{μν}, X = S (axial) or X = T (trace), so p_μ = ξX_μ.

- β = +½∫p_τ dτ = (ξ/2)∫X_t dt, with X_t the **lower cosmic-time** component (X_τdτ = X_tdt).
- Parameterize X_t = ε_X(z) H(z). Then **β = (ξ/2)∫₀^{z*} ε_X(z) dz/(1+z)**. For constant ε: β = (ξε/2) ln(1+z*) = 3.50 ξε rad, so 0.215° needs ξε = 1.07×10⁻³ (skeptic 020 value).
- For the axial vector S_t = -S^0̂ (S^0̂ is the comoving-observer upper component used in CLAIMS EC-01).
- Achromatic: the CFJ dispersion ω² = k² ∓ p_τk/a gives a k-independent rotation for k ≫ |p| (003).
- After reionization: β_reion/β = ∫₀^{z_re}ε dz/(1+z) / ∫₀^{z*}ε dz/(1+z). Constant ε: ln(1+z_re)/ln(1+z*) = 0.31 (z_re = 7.67) or 0.34 (z = 10).
- If p is closed, p = ξ∂θ, then β = (ξ/2)[θ₀ - θ_*] depends only on the endpoints (`ec06_numbers.py` N0). This is the ALP formula with g_φγ = -ξ.

Parity: S is axial, T is polar. ξS·AF̃ is P-even. ξT·AF̃ is P-odd (explicit parity violation in the photon sector). Both are CPT-odd in the background.

## 1. Gauge invariance (`check_gauge.py`, 8 PASS)

- G1: δ_χ(-½p_μA_νF̃^{μν}) = ∂(…) + ¼χ(dp)_{νμ}F̃^{μν}. Gauge invariant iff dp = 0.
- G2: the Euler-Lagrange expressions obey ∂_νE^ν = -¼(dp)_{νμ}F̃^{μν} identically. With a conserved current the field equations then demand (dp)·F̃ = 0 on every solution. That is not a small breaking: the coupled Maxwell system is inconsistent unless dp = 0.
- G3: on FRW, p = p_τ(τ)dτ is closed. Perturbed p stays closed only if it is the gradient of one scalar (δp_τ = ∂_τs, δp_i = ∂_is). An independent δp_τ or any transverse vector part breaks it.
- G4: a Stückelberg completion A → A + ∂σ turns the breaking into a mixing ¼σ(dp)·F̃. σ has no kinetic term unless the photon has a mass, so the completion is a massive photon with a longitudinal mode mixed in by |dp|. That mixing grows as m_γ → 0. Not quantified here (open issue 2).

So a CFJ coupling to torsion is consistent only if the torsion's photon-coupled vector is the gradient of a scalar field on and off FRW. Then it is an ALP coupling θFF̃ in disguise.

Gauge-invariant operators of higher dimension built from O(H) torsion, for example (∇·S)FF̃/M² or the Itin-Hehl F²T² family, are suppressed by (H/M)². β = 0.215° needs M ≈ 5×10⁻²⁸ eV, which is above H* = 3.3×10⁻²⁹ eV by only a factor of 14 (N6). That is not an EFT, so they are dead. Non-polynomial dimension-4 scalars such as ln|T_scalar|·FF̃ are gauge invariant, but they are a P-odd scalar-photon coupling put in by hand.

## 2. (b) Teleparallel: TEGR and f(T) (`check_teleparallel.py`, 22 PASS)

Weitzenböck connection Γ^λ_{μν} = e_A^λ∂_μe^A_ν (ω = 0).

- T1, flat FRW with the diagonal Cartesian tetrad: T_τ = 3a′/a (so T_t = 3H, ε_T = 3), S = 0, T closed (T = 3d ln a).
- T2, TEGR frame dependence: a local boost η(τ,x,y) of the same tetrad leaves the metric unchanged, but shifts T_μ by (-∂_xη, -∂_τη, 0, 0). The result is not closed, and it adds an axial S^z = 2∂_yη/(3a²). TEGR's field equations do not see the local Lorentz frame (the action changes by a boundary term), so the torsion is pure gauge there. **In TEGR a photon coupling to T or S is not defined by the theory**: β takes any value depending on a frame the dynamics does not fix. Adding the coupling makes the frame physical, so TEGR+CFJ is no longer TEGR. It is a Lorentz-violating teleparallel theory whose frame equation is then set by the ξ-term itself.
- f(T) with f_TT ≠ 0: the antisymmetric field equations f_TT S_{[μν]}{}^ρ∂_ρT = 0 select "good" tetrads. T4: the diagonal flat tetrad and both closed-FRW tetrads below satisfy them. The boosted tetrad of T2 does not, so the check has teeth. In f(T) the background torsion is physical. Torsion scalars: 6H² (flat), 6(H² - k/a²) (k = +1).
- T3, closed FRW (k = +1). Tetrads from the left- or right-invariant 1-forms of S³ (unit comoving radius, 3-Ricci = 6 checked) give the same metric, T_τ = 3a′/a, and **S^0̂ = ±2√k/a**, closed on FRW. The sign is the branch (the left/right choice) and is not fixed by the dynamics. These are the two k = +1 branches of Hohmann-Järv-Krššák-Pfeifer 1901.05472 (Eqs. 119-127; their axial content is computed here, not read from the paper). With √k = H₀√|Ω_K|: ε_S(z) = ±2√|Ω_K|(1+z)/E(z), and β = ξ√|Ω_K| H₀χ(z*) with H₀χ* = 3.12 (N2), weighted by conformal distance, so it is late-weighted.
- T5, scalar perturbations with the diagonal perturbed tetrad: T_τ = ∂_τ(3 ln a - 3Ψ), T_i = ∂_i(Φ - 2Ψ), hence (dT)_{τi} = ∂_i∂_τ(Φ + Ψ). This is closed only while Φ + Ψ is static (matter domination). It is broken in the radiation era and at late times, the ISW eras. f(T) perturbations also carry an extra Lorentz scalar mode, which T2 shows is not closed either. So the trace-CFJ coupling breaks U(1) beyond the background, at order ξ k Φ̇.
- Fermions: in TEGR the standard (equivalence) prescription couples Dirac fields like GR (no torsion coupling). If they are coupled through the Weitzenböck connection, the k = +1 axial torsion is 2√|Ω_K|H₀ ≲ 1.6×10⁻³⁴ eV ≪ 7.3×10⁻²² eV (lab). No bound either way.
- Solar system: with closed p the rotation along a ray is (ξ/2)Δθ between endpoints. Near masses θ ~ GM/r ~ 10⁻⁸, so the local signal is ~10⁻⁸ξ: negligible. The CFJ term has no A·B density in static Coulomb fields, so it gives no equivalence-principle violation at this order.

## 3. (a) Magueijo-Złośnik EC with varying Λ (`check_mz.py`, 10 PASS)

Theory (1908.05184 Eq. 7): Palatini + quasi-Euler -(3/2Λ)ε RR + Nieh-Yan/γ + quasi-Pontryagin -(3/γΛ)R∧R, with Λ a field and matter not coupled to the connection. **This is an extension of EC** (quasi-topological terms with 1/Λ). Minimal EC gives no torsion without spin. Torsion forms T^i = -T e⁰∧e^i + P ε^i_{jk} e^j∧e^k; metric Hubble rate H = g - T (their Eq. 19).

- M1: T_t = -3T, S^0̂ = -2P, S_t = +2P.
- M2, γ → ∞ (Euler term only), parity-odd branch (Eq. 37): T = (g² - P²)/2g and H = (g² + P²)/2g. Both components are O(H). On the scaling tracker x² = P²/g² = 1/3, and H²/(κρ/3) = 2/3 (N3a), a 33% deficit in the expansion rate. That is excluded by N_eff (data: N_eff ≥ 2.65 at 95% allows 5.3%). MZ themselves call the γ → ∞ scaling solutions unviable. Their wording is "an additional radiation component with half the photon+neutrino density": the sign differs from my H = g - T reading, but both are excluded. 2001.06373 finds tensor perturbations undetermined with the Euler term alone. **Dead background.**
- M3, finite γ (Pontryagin term on), from Eqs. 26, 28, 29: T/g = γx(1 - x² - 2x/γ)/(1 - x² + 2γx), x = P/g. Algebraic once P ≠ 0.
- M4/M5: the finite-γ tracker (Eqs. 57-59, 64-66 with R_Λ = 0) has R_c ≈ γ²/9 (→ 1/2 as γ → ∞), x² = R_c/(1+R_c), and an expansion-rate deficit 1 - γ²/9. N_eff alone would allow γ ≤ 0.81, i.e. ε_S ≤ 0.49.
- **GW speed kills it.** 2001.06373 finds the tracker P² = γ²κρ/27 (= our x ≈ γ/3, N3c) and tensor modes with ω² = [1 + (1+3w)γ²/9]k² (Eq. 72); in the Λ era ω² = [1 - 8√3γε_P + 84ε_P²]k², ε_P = P/√Λ (Eq. 65). GW170817 (|Δv/v| ≤ 3×10⁻¹⁵, 1710.05834) gives γ ≤ 2.3×10⁻⁷, so ε_S ≈ 2γ/3 ≤ 1.5×10⁻⁷ and ξ ≥ 7×10³ (N3d). With the paper's own γ² ≲ 10⁻¹⁵ the result is ε_S ~ 2×10⁻⁸ and ξ ≳ 5×10⁴, matching the scout (024). The Λ-era term alone gives ε_S(0) ≤ 1.2×10⁻⁸. Their validity condition |γ| ≪ 10⁻²¹ (Eq. 70) would push ξ beyond 10¹⁸. The parity-odd torsion that would rotate photons also changes the GW speed, and the GW bound wins by 7-22 orders in ξ. These equations were used as published, not rederived.
- Off FRW: P is a genuine new degree of freedom in this branch (MZ Sec. V, "c" has its own dynamics), and the torsion equation (MZ Eq. 15) ties T to dΛ ∧ R, not to a gradient. Nothing makes δS a gradient of one scalar, so by G3 the CFJ term breaks U(1) on perturbations. This is not proven for a specific mode (open issue 1).
- Fermions: MZ couple spinors to the torsion-free connection, so there is no b_μ. With the full connection, |S_0| = 2|P| ≲ 0.5H₀ ~ 7×10⁻³⁴ eV ≪ 7.3×10⁻²² eV.
- Λ variation: on the finite-γ branch Eq. 29 fixes Λ̇/Λ² = 2γP/[3(g² - P² + 2γgP)], so Λ̇ ∝ γP. A dark-energy w bound limits γx today, but the late-time solution is not known. Not computed.

## 4. (c) PGT torsion cosmologies

- 0⁻ mode (BHN 1009.5112, Chen et al. 0908.3323): S = ∂φ/√A is closed (007). CFJ is then exactly the ALP coupling g = -ξ/√A. This is CLAIMS EC-04: gauge invariant, but an ALP with the coupling put in by hand, and S ~ H needs m ~ H₀ (~10¹²⁰ tuning). Nothing new.
- 0⁺ mode (Shie-Nester-Yo 0805.3834, L = -(a₀/2)R + (b/24)R² + (a₁/8)(torsion²), no-ghost a₁ > 0, b > 0): only the trace T_t = Φ is nonzero, and it does not couple to fermions (SNY). SNY Eq. 46 gives **Φ = -(3/2)d ln(R + 6μ/b)/dt** (N0b): on FRW the trace is the gradient of the 0⁺ scalar ln(R + 6μ/b). I expect this to hold covariantly when only the 0⁺ mode propagates, as for f(R) with torsion (T_μ ∝ ∂ ln f′, Capozziello et al. 0708.3038), but it is not checked off FRW. If it holds, ξT·AF̃ = -(3ξ/4)·ln(R + 6μ/b)·FF̃-type up to a total derivative: gauge invariant, P-odd, an ALP-like coupling of the 0⁺ mode. Size: Φ/H ≈ 0.8-1.4 today (SNY cases III-V), oscillating. Early times are not reliable in SNY, so ε(z*) is unknown. The SNY background differs from ΛCDM (Minkevich et al. 0811.1430). Yo-Nester found only the 0± modes free of nonlinear constraint pathologies. Any 1± vector-torsion mode gives non-gradient T and S, so G3 kills the CFJ coupling there.
- Minkevich-type PGT cosmologies (torsion functions ~ H) were not examined.

## 4b. Relation to published scalar-curvature couplings (input from task 024)

On flat FRW the f(T) trace is T_t = 3H, i.e. p = 3ξ∂ ln a. Yu-Yang-Qiu 2508.10532 couple ln Q·FF̃ (Q = 6H², which equals the f(T) torsion scalar on flat FRW, T4). Li-Xia-Li hep-ph/0611192 use p ∝ ∂ ln R. In matter domination all three are ∂ ln a up to constants (c = -ξ). Over the full history β per unit coefficient is 10.49 (trace), 10.05 (ln Q) and 9.36 (ln R) (N7). The ln T_scalar coupling is the gauge-invariant completion of the trace CFJ term, so the flat-f(T) EC-06 row is the published phenomenology in torsion language. It differs only off FRW, where the CFJ form breaks U(1) (T5) and the ln-scalar form does not. A non-closed p needs either a gravitational CS counterterm (Li-Cai-Wang 0907.5159) or a photon mass (Hook-Huang-Zhou 2609.10671); those two are cited via 024 and not opened here.

## 5. Bounds common to all three (`ec06_numbers.py`)

- Lab b_μ: |S_0| < 7.3×10⁻²² eV vs H₀ = 1.4×10⁻³³ eV: binds only for ε > 10¹¹. Never.
- Photon-sector LV today: p₀ = ξεH₀ ~ 1.5-3.5×10⁻³⁶ eV, |k_AF| ~ p/2, vs the radio-galaxy CFJ bound 6×10⁻³⁴ eV (bounds.csv P-kAF-CFJ90): 2.5 orders below. UV radio-galaxy rotation at z > 2 (-0.8° ± 2.2°, 1003.4823): trace model gives 0.034° at z = 2. Pulsars (10⁻²⁷ eV) are far weaker. Low-z data cannot see ξε ~ 10⁻³. The CMB β is the k_AF measurement.
- No theory bounds ξ, because none of them contains ξ. The coupling is added. Bounds on ξ come from ACB (the torsion perturbations rotate too) and from the gauge-breaking cost above.
