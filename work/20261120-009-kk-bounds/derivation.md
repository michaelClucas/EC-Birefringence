# KK bounds: derivations (task 20261120-009)

Conventions: `CONVENTIONS.md`. KK ansatz and the KK-02 and KK-03 couplings are as in `work/20261120-004-kk-reduction/derivation.md`. Check labels [P*], [E*], [O*], [R*], [H*], [K*] refer to `run.log`. Numbers are in `results.csv`. All measured inputs come from `data/`.

## 1. The proper-rotation parity P' (skeptic input)

P' = diag(1,-1,-1,-1,-1) on (x⁰,x¹,x²,x³,y) equals exp[π(J₁₂ + J₃y)], a rotation by π in two planes [P1]. It lies in the identity component of SO(4,1), so every diffeomorphism-invariant 5D action is invariant under it, with or without 5D parity. y → -y maps S¹ to itself. On the zero modes it acts as 4D parity P combined with:
- KK photon: A'_0(x) = -A_0(Px), A'_i(x) = A_i(Px) (from g_{μ5} → -g_{μ5}).
- radion: φ' = φ(Px) (even).
- bulk U(1): a_μ is an ordinary polar vector, and the Wilson line θ = 𝒜_5 is odd.
- torsion pseudoscalars: the FRW axial component S⁰ ∝ ε^{0ijk}T_{ijk} is odd.

Result: φ^p F_KK F̃_KK is P'-odd [P4], and φ³F_KK² is P'-even [P5]. **If the vacuum is P'-invariant, no function of the radion multiplies F_KK F̃_KK at any loop order.** This holds for any field content that is P'-covariant: 5D Dirac fermions, including their masses, transform under the proper rotation, so a 5D mass is not odd under it. A parity-odd photon term therefore needs either a P'-odd field (θ in KK-02 [P6], the KR dual axion in KK-03, a torsion pseudoscalar bilinear) or a P'-odd background (an oriented vector u, a homogeneous S⁰, a net chiral asymmetry). This extends the tree-level KK-01 no-go of task 004 to all orders for P'-even scalars.

## 2. What 2608.09299 actually adds

With e⁵ = n = φ(dy + √2A) (their Eq. 30) and T⁵ = dn in Weitzenböck gauge, the A = B = 5 part of their Eq. 39 is

  (c/8) ε^{MNPQL} (dn)_{MN}(dn)_{PQ} n_L ∝ n∧dn∧dn = 2φ³ dy∧dA∧dA   [P2, the dφ∧dy terms cancel identically],

and the A = B = a part is φ dy∧T^a∧T_a with no photon [P3]. The EM term is therefore the 5D abelian CS form of the 1-form n. It is built from the metric and the Killing direction alone, so it does not need teleparallel torsion. It is P'-odd with u = +n [P7] because P'*n = -n[φ', A'] [P8]. The P'-odd datum is the orientation of u relative to ∂_y, which in their construction is the sign of det(∂_μψ^a) of the Stückelberg vacuum. That is a parity-odd background: their 5D action is a scalar under diffeomorphisms, and the vacuum ψ^a = x^a breaks P'. It is not a boundary term or a non-Lagrangian term.

## 3. Einstein-Cartan analogue

Pointwise algebra in an orthonormal frame. ω = ω̊ + K, T^A = K^A{}_B∧e^B (T̊ = 0), R = R̊ + D̊K + K∧K. The EC action has the K-quadratic part ½Q(K) (κ = 1) [E2: same normalization as 004; E4: non-degenerate]. With an oriented unit 1-form u and W = du:

- Identity: T^A∧T_A = e^A∧e^B∧(K∧K)_{AB} [E1]. So u∧NY = u∧(T∧T - e∧e∧R) ≅ -du∧e^A∧e^B∧K_{AB} (using e∧e∧R̊ = 0 [E3] and e∧e∧D̊K = d(e∧e∧K)). This is linear in K and has no odd K² piece.
- u-Holst, c u∧e^A∧e^B∧R_{AB}(ω) ≅ c[du∧e^A∧e^B∧K_{AB} + u∧T^A∧T_A]. In Riemannian geometry it vanishes identically (2608's P̂₁).

Eliminating K for each 5D term:

| 5D term | L_eff(W) | photon parity-odd? |
|---|---|---|
| c u∧T^A∧T_A (2608 transplanted to EC) | 0 [E6] | no |
| c u∧NY | -(3c²/2)W·W, even only [E5] | no |
| c u∧e∧e∧R (u-Holst) | -(3c²/2) W·W - 4c³ ⋆(u∧W∧W), exact for all c [E7, E8] | **yes** |

Adding the Maxwell term from R̊ (-W·W/8) and writing ⋆(u∧W∧W) = ½W W̃:

  L = -(1/8κ₅)(1 + 12c²)[W·W + c_odd W W̃],  c_odd = 16c³/(1 + 12c²).

For u = e⁵ of the KK vielbein, W_{ab} = e^{βφ}kF_{ab}, so in 4D L = -¼Z(φ)[F² + c_odd FF̃] with Z ∝ e^{-√3φ} ∝ 1/α. This is an EC radion-FF̃ coupling that exists only because of torsion. It needs the same oriented u (P'-odd background) as 2608.

## 4. Rotation and the α link

For L = -¼Z(F² + c_odd FF̃), integrating by parts (FF̃ = ∂K, K^μ = 2A_νF̃^{μν}) gives CFJ p_μ = -c_odd ∂_μ ln Z for the canonically normalized photon. With β = ½∫p_τ dτ (CONVENTIONS), β = -(c_odd/2)Δln Z [R1, numeric helicity evolution]. Since α ∝ 1/Z for both the KK photon and a bulk photon:

  **|Δα/α|_{LSS→0} = 2|β|/|c_odd|.**

For 2608, Z = L₅φ³ and c_odd = 2c, so β = 3cΔln φ (their Eq. 52) and |Δα/α| = |β|/|c| [R2].

## 5. KK-02 (bulk U(1) + 5D CS, Wilson-line axion)

- From 004: g = 24πκ_CSα/f_a, f_a = 1/(2πLg₄).
- **Quantization.** The CS term changes under large gauge transformations, so κ_CS = k/(24π²). With k ∈ Z the term is well-defined on spin 5-manifolds once the mixed gauge-gravitational CS term is included; otherwise k ∈ 6Z (standard result, not re-derived here). Integrating out a charged bulk fermion shifts k by ±½. Then g = kα/(πf_a), and the 4D term is the standard θ-term (k/8π²)ϑ f∧f, with period 2π for k = 1.
- **β is independent of L and f_a:** β = (g/2)Δφ = kαΔϑ/(2π). With Δϑ ≤ π (a compact field rolling to its minimum), |β| ≤ kα/2 = 0.209° k [K1].
  - ACT 0.215 ± 0.074° needs Δϑ = 1.03π (1σ: 0.67π to 1.38π) for k = 1. So k = 1 must start near the hilltop.
  - k = 2 needs Δϑ = 0.51π (0.34π to 0.69π), and k = 3 needs 0.34π: generic initial conditions.
- **L window.** Chandra g < 8×10⁻¹³ GeV⁻¹ (m < 10⁻¹² eV) gives f_a > 2.9×10⁹ k GeV, i.e. 1/L > 5.5×10⁹ k GeV. The 5D EFT needs 1/L < M_5, i.e. 1/L < M_P/√(2π) ≈ 10¹⁸ GeV. Collider KK bounds (TeV) are irrelevant.
- **Mass.** With no charged bulk matter, ϑ is massless perturbatively (the KK modes of an abelian field are neutral, and brane fields do not wind). A charged bulk Dirac fermion (charge q, mass M, circumference l = 2πL) gives V = (3/π²l⁴)Σ(1 + nx + n²x²/3)e^{-nx}cos(nqϑ)/n⁵ with x = Ml [H1: closed form matches a direct KK proper-time sum]. Then m ≈ qg₄x e^{-x/2}/(πl) [H2]. Rolling between recombination and today needs m between about H₀ and H_* (1.4×10⁻³³ to 3.3×10⁻²⁹ eV). That gives x = Ml ≈ 220-280 over the whole L window: the fermion must sit at about 35-45 × 1/L, below the 5D strong-coupling scale (Λl ≈ 24π³/g₄² ~ 10⁴). This is a single tuned-but-technically-natural number (exponential suppression), the extra-natural mechanism of hep-th/0301218.
- **Massless alternative.** ϑ is a free field kicked before BBN. Its excursion is not bounded by π, and the stiff energy needs Δ(f_aϑ) < 5.5×10¹⁰ GeV between recombination and today (BBN ΔN_ν bound, Sec. 6). With k = 1 that means f_a < 1.7×10¹⁰ GeV, inside the Chandra window. This is viable, but the initial velocity is an unexplained input.
- **Anisotropy.** Even O(π) superhorizon fluctuations of ϑ give at most (kα/2)² = 1.3×10⁻⁵ k² rad² of rotation variance, below SPT A_CB < 1.2×10⁻⁴ for k ≤ 3.
- **Radion.** β does not depend on L, so a moving radion changes only α: δα/α = √(2/3)δσ/M_P. Standard stabilization with m_σ > 5×10⁻³ eV (ISL, gravitational strength) removes α variation. No tension.
- **QCD mass.** None, as long as only hypercharge is in the bulk with the CS term (no 𝒜∧tr G∧G).

## 6. KK-03 (KR 2-form = axial torsion + Green-Schwarz)

- g = √6|λ|g₄²M_P (004). **λ is free** in EC-KK. Torsion carries no flux quantum (T is not the curvature of a bundle), so B is a non-compact 2-form, its dual b is non-compact, and nothing quantizes the GS coefficient. If one instead postulates charged strings (a compact B), then b becomes periodic, λ is quantized, and KK-03 becomes KK-02 with a compact b (β = k'αΔϑ/2π). That uses an input EC does not supply.
- **Dynamics.** b has an exact shift symmetry, and U(1) has no instantons. So b is massless, and a homogeneous b moves only by an initial velocity, ḃ ∝ a⁻³ (stiff, ρ ∝ a⁻⁶). The displacement after recombination is Δb = ḃ_* a_*³ ∫da/(a⁴H). Requiring ΔN_ν(BBN) < 0.16 (Fields et al. +2σ) gives Δb < 5.5×10¹⁰ GeV. The CMB ΔN_eff bound gives only 4×10¹⁷ GeV. With β = (g/2)Δb this needs g > 1.4×10⁻¹³ GeV⁻¹. Chandra caps g < 8×10⁻¹³ GeV⁻¹.
- **Window:** λ/κ₄ ∈ [1.5×10⁶, 8.7×10⁶]. Written as λ = 1/M_s², this is M_s ≈ (0.8-2)×10¹⁵ GeV.
- **Natural value excluded.** λ = κ₄ gives g = 9×10⁻²⁰ GeV⁻¹, and then ΔN_ν(BBN) ~ 3×10¹¹.
- **The only torsion-specific signal** is the axial torsion carried by the moving b, S⁰ ~ g_Bḃ₀ ~ 7×10⁻⁴⁶ eV, about 20 orders below the neutron b_⊥ bound after the Sun-CMB boost. So it is unobservable.

## 7. S¹/Z₂ with a brane photon (goal 4)

- **Bulk U(1)' Wilson line.** 𝒜'∧𝓕'∧𝓕' is kept for A'_μ odd and A'_5 even [O1]. The opposite assignment needs an ε(y) coefficient and kills the zero mode [O2]. A brane-localized photon cannot couple to ϑ' at tree level:
  - A'_5(0)F F̃ is not gauge invariant, because odd gauge parameters have ∂_yΛ(0) ≠ 0 [O5].
  - Brane fields cannot carry U(1)' charge, since A'_μ(0) = 0.
  - Derivative couplings ∂ϑ'·J₅ to brane fermions decouple below the fermion mass.
  
  A photon coupling arises only if hypercharge also lives in the bulk with A_Yμ even, through the mixed term 𝒜'∧𝓕_Y∧𝓕_Y. That term is orbifold-allowed [O3] and reduces to ϑ' f_Y f̃_Y [O4]. This is KK-02 with a second bulk field, and it is not a brane photon.
- **Bulk KR 2-form.** B_{μν} even (the axial torsion T_{abc} is even, so T = H is consistent with the orbifold), plus a brane-localized Green-Schwarz term dH ⊃ λ_b δ(y) dy∧F∧F (the Hořava-Witten type). This gives the 4D coupling b F F̃ with g ∝ λ_b g₄²/g_B: the KK-03 structure with a free brane coefficient. It is not symbolically checked here (only the parity counting is standard).
