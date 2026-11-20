# Derivations: Einstein-Cartan torsion and the photon (4D track)

Task 20261120-003. Conventions: `CONVENTIONS.md` (signature -+++, ε_{0123} = +√(-g), ∇_μV^λ = ∂_μV^λ + Γ^λ_{μν}V^ν, T^λ_{μν} = Γ^λ_{μν} - Γ^λ_{νμ}, Γ = Christoffel + K), except where this file reports an error in it. Each section names the script that checks it. All checks: `run_checks.py`.

The ACT DR6 β = 0.215° ± 0.074° (arXiv:2509.13654) is 2.9σ with systematics caveats; the signal may not be cosmological. Everything below is a candidate mechanism, not an explanation.

## 1. Contortion (check_contortion.py)

Metric compatibility ∇_μ g_{αβ} = 0 with Γ = {} + K requires K_{βμα} + K_{αμβ} = 0, i.e. K_{λμν} ≡ g_{λρ}K^ρ_{μν} is antisymmetric in its first and third indices. Torsion: T_{λμν} = K_{λμν} - K_{λνμ}. Using the antisymmetry,

T_{λμν} = K_{λμν} + K_{μνλ},  T_{μνλ} = K_{μνλ} + K_{νλμ},  T_{νλμ} = K_{νλμ} + K_{λμν},

so

**K_{λμν} = ½ (T_{λμν} + T_{νλμ} - T_{μνλ})**   (equivalently ½(T_{λμν} - T_{μνλ} - T_{νμλ})).

CONVENTIONS.md writes K_{λμν} = ½(T_{λμν} - T_{μλν} - T_{νλμ}). The script finds that this expression equals K_{μνλ} = -K_{λνμ} of the correct contortion. It reproduces the torsion but is **not metric compatible** in our index placement (it is the correct formula for the convention where the derivative index of Γ is the last one). Checked on two random metrics with random derivatives and random torsion (exact rationals).

## 2. Minimal EC + Maxwell gives zero rotation (check_maxwell_zero.py, check_dirac_torsion.py)

(a) F = dA, F_{μν} = ∂_μA_ν - ∂_νA_μ, has no connection. S_Maxwell = -¼∫√(-g) g^{μα}g^{νβ}F_{μν}F_{αβ} involves only the metric. On FRW, g = a²η, √(-g)g^{μα}g^{νβ} = η^{μα}η^{νβ} (checked), so the comoving equations are flat-space Maxwell: no rotation of the polarization plane. In minimal EC the metric itself is changed by torsion only through the spin-spin contact energy, which keeps FRW form; the photon still propagates with conformally invariant Maxwell equations.

(b) If one instead "minimally couples" with the full connection, F'_{μν} = ∇_μA_ν - ∇_νA_μ = F_{μν} - T^λ_{μν}A_λ (checked for a generic Γ; the Levi-Civita part cancels because it is symmetric). Under A → A + dχ, δF' = -T^λ_{μν}∂_λχ ≠ 0: U(1) gauge invariance is broken. On FRW torsion (S^0, trace h) the script finds

-¼F'F' = -¼FF - ½ p_μ A_ν F̃^{μν} + ⅓ T_μ A_ν F^{μν} + ½(h² - s²) A_iA_i   (flat, a = 1, s = S^0),

with **p_μ = 2 S_μ**. So this gauge-breaking prescription contains a CFJ term with p = 2S (the ξ = 2 case of candidate (i)), a parity-even gauge-breaking term with the trace vector, and photon mass terms m² ~ s² - h². It is not a consistent U(1) theory; it is recorded only because it is the "naive minimal coupling" often used in the older literature (e.g. Poplawski 1108.6100, Garcia de Andrade hep-th/0110150, both listed in `literature/prior_work_torsion_birefringence.md`).

(c) Fermion route at tree level. Dirac conventions (Srednicki): {γ^a,γ^b} = -2η^{ab}, L_D = (i/2)(ψ̄γ^μ∇_μψ - ∇_μψ̄ γ^μψ) - mψ̄ψ, γ_5 = iγ^0γ^1γ^2γ^3, J_5^μ = ψ̄γ^μγ_5ψ (J_5^0 = n_R - n_L). The spinor connection is Ω_μ = -¼ ω_{μab}γ^aγ^b (fixed by requiring ψ̄γ^aψ to be a vector), with torsion part ω_μ^a_b = K^a_{μb}. Results (checked with explicit 4×4 matrices and random torsion):

- L_D ⊃ c_A S_μ J_5^μ, **c_A = 3/4**; trace and pure-tensor torsion drop out of the Hermitian Dirac action.
- R(Γ) = R̊ - ⅔ T_μT^μ + (3/2) S_μS^μ + (pure tensor part) + total derivative, no cross terms (constant torsion in flat space fixes the coefficients).
- With L = (1/2κ)R + L_D, S is algebraic: **S^μ = -(κ/2) J_5^μ**, and **L_eff = -(3κ/16) J_5·J_5** (the standard axial-axial contact term; magnitude 3κ/16 matches the literature, and c_A = 3/4 matches the b_μ = -(3/4)A_μ map of Kostelecky-Russell-Tasson 0712.4393 as recorded in `data/bounds.csv`, up to their normalization of A_μ).

No A_μ appears. At tree level minimal EC with Dirac matter has no photon-torsion coupling, so β = 0. A photon coupling can only come from a loop (Sec. 4 (iv)) or an added term.

## 3. FRW-compatible torsion (check_frw_torsion.py)

Spatially flat FRW in conformal coordinates. Homogeneity (Lie derivative along ∂_i) makes components τ-dependent only; isotropy (Lie derivative along x^j∂_k - x^k∂_j) is then linear algebra on the 24 components. Solving it gives exactly two free functions:

T^i_{0j} = -T^i_{j0} = h(τ) δ^i_j,  T^i_{jk} = f(τ) [ijk],  all other components zero.

On g = a²η: **S^τ = -f/a²**, S^i = 0; **T_τ = 3h**, T_i = 0. f is parity odd (pseudoscalar), h parity even. The comoving-observer component of the axial vector is S^0̂ = a S^τ; its covariant conformal component is S_τ = -a S^0̂. This confirms the FRW statement in CONVENTIONS.md (with T_{i0j} = a² h δ_{ij}).

## 4. Rotation from a CFJ term, sign, and the φFF̃ baseline (check_cfj_rotation.py)

**φFF̃ ↔ CFJ.** F_{μν}F̃^{μν} = 2∂_μ(A_νF̃^{μν}) (checked), so -(g/4)φFF̃ ≅ +(g/2)∂_μφ A_νF̃^{μν} = -½ p_μA_νF̃^{μν} with **p_μ = -g ∂_μφ**. The Euler-Lagrange equations agree for p = -g∂φ and disagree for p = +g∂φ. CONVENTIONS.md states p_μ = +g∂_μφ: sign error.

**Field equation.** ∂_μF^{μν} = p_μF̃^{μν} (EL equations, constant p). Both terms are metric independent after multiplying by √(-g) on FRW (Sec. 2a and √(-g)ε^{μνρσ} = -[μνρσ]), so in comoving coordinates p_μ is the conformal-coordinate component p_τ. With A_0 = 0, ∇·A = 0, p = (p_τ,0,0,0), coordinates (x,y,z) right-handed:

A'' - ∇²A = p_τ ∇×A.

For p_μ = ∂_μθ (any FRW p_τ(τ)) the extra term ∝ ∂_{[ρ}p_{μ]} vanishes, so the same equation holds with p_τ(τ).

**Dispersion.** For e_± = e_x ± i e_y, wave along +z: ω_±² = k² ∓ p_τ k, ω_± ≈ k ∓ p_τ/2.

**Rotation.** e_+e^{iθ} + e_-e^{-iθ} = 2(cos θ e_x - sin θ e_y) with θ = p_τΔτ/2, so the polarization angle measured right-handed about the propagation direction k changes by

χ = -½ ∫ p_τ dτ.

Independent check: direct numerical evolution of the EL equations for a linearly polarized packet (p_τ = 0.05, Δτ = 10) gives χ = -0.25001 rad vs predicted -0.25000.

**Map to the sky and to the paper's β.** The photon moves along +k toward the observer, who looks along -k. A right-handed rotation about +k is counterclockwise as seen by the observer (checked: facing -z with up +y, right is +x, so x→y is right→up). `literature/conventions_beta.md` records that in arXiv:2509.13654 (and Komatsu 2202.13919, Eskilt-Komatsu 2205.13962) β > 0 is clockwise on the sky. Hence β = -χ:

**β = +½ ∫_{τ_LSS}^{τ_0} p_τ dτ.**

For the baseline L ⊃ -(g/4)φFF̃ in our conventions: **β = -(g/2)[φ(τ_0) - φ(τ_LSS)]**. The paper formula β = +(g/2)Δφ holds in Komatsu's convention ε^{0123} = +1/√(-g), which is the opposite orientation sign to ours (ε^{0123} = -1/√(-g)); flipping ε flips F̃ and hence the sign of the term. The two routes (explicit dispersion + sky mapping here, and ε-sign conversion of the review's formula) agree. Equivalently, in our conventions L ⊃ +(g/4)φFF̃ gives β = +(g/2)Δφ.

Assumptions in this mapping: (x,y,z) is right-handed with [123] = +1 (this is what ε_{0123} = +√(-g) means physically), and "clockwise on the sky" means as seen by the observer, as stated in conventions_beta.md for 2601.13624 and 2202.13919.

**Gauge invariance of a CFJ term with a field S_μ.** δ_χ(-½ξS_μA_νF̃^{μν}) = total derivative + (ξ/2)χ ∂_νS_μ F̃^{μν} (checked). It is gauge invariant iff dS = 0, i.e. S_μ = ∂_μϑ locally. The FRW background S_μ = (S_τ(τ),0,0,0) is closed, so the isotropic background is fine; perturbations with a curl part break gauge invariance, so any consistent anisotropic version needs S = ∂ϑ, which makes it an axion-type φFF̃ model with φ = ϑ.

## 5. Nieh-Yan (check_nieh_yan.py)

In coordinates, with the permutation symbol [μνρσ]:

∂_μ([μνρσ] T_{νρσ}) = ½ [μνρσ] T^λ_{μν}T_{λρσ} - [μνρσ] R_{μνρσ}(Γ),

checked on random first-order jets of metric and torsion (coefficients fitted on two configurations and confirmed on a third). Also [μνρσ]T_{νρσ} = -6√(-g) S^μ. So the Nieh-Yan density is N = -6∂_μ(√(-g)S^μ), and a pseudoscalar coupling ∫ θ N d⁴x = +6∫√(-g) S^μ∂_μθ d⁴x after integration by parts.

## 6. Candidate catalogue

Rotation formulas use β = +½∫p_τ dτ (paper convention, positive clockwise). S^0̂(t) is the comoving-observer component, ∫...dt is cosmic time from last scattering to today. Numbers: `estimates.py` (log in run.log). All PRELIMINARY, order of magnitude.

| id | action (our conventions) | U(1) gauge invariant | effective 4D term | β |
|---|---|---|---|---|
| (0) minimal EC + Maxwell + Dirac | (1/2κ)R(Γ) - ¼F² + L_D | yes | none at tree level | 0 |
| (0') naive ∇-coupling | -¼F'F', F' = F - T·A | **no** | CFJ with p = 2S, plus ⅓T_μA_νF^{μν} and mass terms | -∫S^0̂ dt |
| (i) non-minimal CFJ | -½ξ S_μA_νF̃^{μν} | only if dS = 0 (true on FRW background) | CFJ, p_μ = ξS_μ | -(ξ/2)∫S^0̂ dt |
| (ii) Holst with BI field β_BI(x), photon coupling via loop or by hand | (1/2κ)R + (1/2κ)β_BI^{-1}·Holst - (g/4)φFF̃ | yes | φFF̃ with φ = canonical BI field | -(g/2)Δφ |
| (iii) Nieh-Yan pseudoscalar | θ N = 6√(-g)S^μ∂_μθ, plus -(g/4)θFF̃ | yes | eliminating algebraic S gives kinetic term for θ and ∂θ·J_5; photon term must be added or loop-induced | -(g/2)Δθ_can |
| (iv) fermion loop on a torsion background | c_A S_μJ_5^μ (c_A = 3/4), one charged fermion loop | the induced term is, for closed S | CFJ, \|p\| = 2C e² (3/4)\|S\|, C regularization dependent (0 or 3/16π², ...) | ∓C e²(3/4)∫S^0̂ dt (sign scheme dependent) |
| (v-a) PGT with propagating axial torsion | quadratic torsion/curvature terms making S dynamical; still needs (i) or (iv) for photons | as (i)/(iv) | as (i)/(iv) with S = ∂ϑ for the 0⁻ mode | as (i) |
| (v-b) Kalb-Ramond torsion, H = dB + (CS 3-form) | T_{λμν} ∝ H_{λμν}, H dual to ∂a | yes (CS augmentation) | aFF̃ | -(g_a/2)Δa |
| (v-c) Itin-Hehl F²T² constitutive couplings | ℓ²(T∧F)∧⋆(T∧F)-type | yes | induced θ ∝ ℓ² u v from trace (h) × axial (S) | -(c/2)Δθ, θ ∝ ℓ² h S^0̂ (normalization not derived here) |

Size estimates (from `estimates.py`):

- Baseline: |gΔφ| = 2β = 7.5×10⁻³ rad; β > 0 needs gΔφ < 0 in our conventions.
- (i) constant S^0̂: |ξS^0̂| = 1.1×10⁻³⁵ eV = 7.9×10⁻³ H_0. If S^0̂ ∝ a⁻³: |ξS^0̂(today)| = 2.8×10⁻⁴⁰ eV.
- Minimal EC source S^0̂ = -(κ/2)n_5: relic-neutrino asymmetry |ξ_ν| = 0.1 (placeholder bound), 3 flavors, relativistic, n_5 ≈ 31 cm⁻³ → |S^0̂| ≈ 2×10⁻⁶⁸ eV today. With a⁻³ scaling from last scattering this needs |ξ| ≈ 10²⁸: **28 orders of magnitude short at ξ = 1**. Chiral fermionic dark matter (1 keV, fully chiral, generous) is 29-30 orders short. Real neutrinos are non-relativistic today, so n_5 is further suppressed by v. Spin densities of unpolarized or randomly polarized matter average to zero by isotropy; only a net chirality (pseudoscalar) density sources S^0̂.
- (0') same as (i) with ξ = 2, plus photon mass ~|S^0̂|, |h|.
- (iv): |p|/|S| = 2.6×10⁻³ for one charged fermion with C = 3/16π², so |S^0̂| ≈ 4×10⁻³³ eV is needed. That is far below the lab bounds on axial torsion from b_μ (~10⁻²² eV spatial, ~10⁻¹⁸ eV time-like, `data/bounds.csv` T-KRT rows), so lab bounds do not exclude it, but minimal EC supplies 10⁻⁶⁸ eV and C is not fixed.
- (ii)/(iii): if the photon coupling is only loop-induced from ∂φ J_5/f with heavy charged fermions, it is suppressed by ~(ω/m_f)² (≈ (H_LSS/m_e)² ~ 10⁻⁶⁹) according to the decoupling argument; Castillo-Felisola et al. 1502.03694 instead keep the anomaly term αQ̄²/(π f_Φ) with f_Φ = √(8/3)M_P. Which is right for a slowly varying background is an open issue. If g ~ 1/M_P is assumed, |Δφ| ≈ 2β M_P ≈ 1.8×10¹⁶ GeV; the same universal coupling to quarks generates a QCD anomaly mass for the field (1502.03694), which then washes out (2512.21888).
