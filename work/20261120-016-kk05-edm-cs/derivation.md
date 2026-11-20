# Derivations: KK-05 electron EDM, P'-odd background, CS level (task 20261120-016)

Conventions: `CONVENTIONS.md`, 5D frame η = diag(-1,1,1,1,1), index 4 ≡ 5̂. Check labels [D*], [E*], [G*], [C*] refer to `run.log`. Measured inputs come from `data/` (JILA eEDM and MICROSCOPE in `shift_uholst_inputs.csv`, clocks and Eot-Wash in `bounds.csv`, HfF⁺ oscillation band in the new `kk05_edm_inputs.csv`).

## 1. Setup and normalisation (re-derivation of 12c² and 16c³)

KK vielbein e⁴ = dy + kA (radion stabilised), W = de⁴ = kF, k = √(2κ₄). Action (κ₅ = 1):

  S = ∫ [ ½R(ω) + c(χ) ⋆(u∧e^A∧e^B∧R_AB(ω)) ],  u = e⁴ with a fixed orientation (P'-odd).

With ω = ω̊ + K and D̊e = 0, the u-Holst D̊K term integrates by parts to X∧e^A∧e^B∧K_AB with X = d(cu) = cW + c′dχ∧u. The K-dependent Lagrangian is

  L_K = ½Q(K) + c·⋆(u∧e^A∧e^B∧(K∧K)_AB) + ⋆(X∧e^A∧e^B∧K_AB).

This is implemented from scratch (sympy, exact rationals, Riemann tensor of a constant connection, 50 contortion components), not by reusing 009/012 code:

- [D1] The Levi-Civita connection of the KK vielbein solves de = -ω̊∧e and has ω̊_[ab5] = W_ab/6 (004's value).
- [D2] ½R̊ = -⅛W·W: the Maxwell term with the healthy sign, so W = kF is canonically normalised at c = 0.
- [D3] Exact elimination for c = 1/3, 2, 7 gives L_eff = -(3c²/2)W·W - 4c³⋆(u∧W∧W). With Maxwell: **Z = 1 + 12c², Y = 16c³**, confirming 009 E8 and 012 U1. So dβ = -c dlnα (012 U5) rests on verified coefficients.
- [D4] The contortion response is K_[ab5] = -(4/3)c²W_ab ± cW̃_ab exactly, and K_[abc] = 0 for 4D indices (012 U3).

## 2. Fermion coupling

5D Dirac fermion, L = iΨ̄Γ^C D_CΨ - MΨ̄Ψ with {Γ^A,Γ^B} = -2η^{AB}, Γ⁴ = iγ₅. After hermitian projection only the totally antisymmetric connection survives [D7a]:

  L ⊃ (3i/4)(ω̊ + K)_[ab5] Ψ̄Γ^{ab5}Ψ + (4D-index terms).

**Why the fermion sees the same K_W.** The full K-Lagrangian is quadratic, ½K·M·K + K·(S_W + J), where J is the fermion bilinear. The cross term after elimination is -S_W M⁻¹J = K_W·J, with K_W the contortion response to W alone. So the dipole is the minimal coupling evaluated on ω̊ + K_W. The JJ terms are four-fermion contact terms.

**EDM or MDM.** Γ^{ab5} = σ^{ab}γ₅ up to a factor and sign. A real-mass electron at rest gives L = φ†(D_E σ·E + D_B σ·B)φ per unit W [D5, D6]:

  **D_E = (3/2)|1/6 - (4/3)c²|,  D_B = (3/2)c.**

- The P'-even pieces (Levi-Civita W/6 and the c² W term) are the **electric** dipole.
- The P'-odd piece cW̃ is the **magnetic** one, since σ^{ab}γ₅W̃_ab ∝ iσ^{ab}W_ab.

This is the reverse of 012's assignment. The reason is that the KK photon is axial under P' (009 Sec. 1), so P'-parity and 4D parity of a dipole operator differ by a sign. At c = 0 the pure KK Pauli term is already an EDM, d = k/4. This is the 5D CP-violation noted by Thirring (1972; not opened).

**Physical units.** The canonical photon is F_c = √Z F, so

  **d(c) = k(3/2)|1/6 - (4/3)c²| / √(1+12c²)  →  (2/√12) c k = 0.577 c k  (c ≫ 1).**

The EDM grows linearly in c. The c² growth of K beats the √Z ~ c suppression.

- **012** used the c-odd piece and omitted Z. Its estimate (3/4)ck had the wrong piece but nearly the right slope (0.75 vs 0.577).
- **Skeptic 015** used (3/4)ck/√Z, which is the canonically normalised **MDM** piece. It saturates at 8.2×10⁻³³ e·cm [E1]. 015 drops the (4/3)c²W term, so its conclusion "no EDM bound on c" does not hold.

## 3. Conditions for the bound

1. **Vector-like electron.** The electron is a 5D bulk Dirac zero mode. On S¹ with the U(1)_y isometry intact (needed for a massless KK photon), the Dirac index is zero and zero modes are vector-like (Witten, Nucl. Phys. B186 (1981) 412; not opened). There are two ways to get chiral zero modes, and both remove the massless KK photon:
   - the S¹/Z₂ orbifold projects out g_μ5;
   - y-dependent backgrounds break the isometry.

   A single-chirality zero mode would have no tree-level Pauli term [D7b].
2. **Mass phase.** The only P'-even 5D mass is MΨ̄Ψ, which gives θ = 0. The oriented u also allows m′Ψ̄u̸Ψ = m′Ψ̄Γ⁴Ψ, a γ₅-mass. With phase θ the EDM is d = k|D_E cos θ ± D_B sin θ|/√Z.
3. **Electron charge.** The electron must couple to the KK photon with charge e. Zero modes are KK-neutral, and since A_KK is P'-axial, a vector-current coupling of a zero mode is itself P'-odd. This is the KK charge problem again (012 open issue 6). The bound is computed relative to that charge coupling.

## 4. Numbers [E2-E5]

| bound | requirement on c |
|---|---|
| JILA eEDM < 4.1×10⁻³⁰ e·cm (2212.11841 Eq. 4), θ = 0 | **c ≤ 188** |
| MICROSCOPE (012 formula), Δχ = 0.39 / 1 / 2 M_P | c ≥ 340 / 133 / 66 |
| clocks (field still rolling today; m = 10⁻³³ to 10⁻²⁹ eV) | c ≥ 5.1×10⁴ |

- EDM + clocks: excluded by a factor of 273.
- EDM + MICROSCOPE: closed for Δχ < 0.71 M_P (factor 1.8 at 0.39 M_P).
- To evade the EDM at c = c_clock, the mass phase needs |cos θ| < 3.6×10⁻³. That means the electron mass comes ≥ 99.6% from the P'-odd m′Ψ̄u̸Ψ term [E5].

## 5. Rolling χ as a CP-odd background today [D8, E6, E7]

χ is P'-even. The static odd datum is the orientation of u times c(χ₀). A rolling χ adds the source c′dχ∧u.

- **Axial background.** This source gives K_[123] = ∓c′χ̇, independent of c. Its effective term is 3c′²χ̇², which is 012's -3c′²(∂χ)². It produces an axial background for every bulk fermion of either chirality: b₀ = 1.5c′χ̇. Today ċ₀ ≈ (β/2)H₀R with R = O(1), so b₀ ≈ 8×10⁻³⁶ eV. The Eot-Wash b̃_T limit is 1.1×10⁻¹⁷ eV (bounds.csv F-e-bT), 18 orders above.
- **EDM drift.** The EDM changes by (1/d)ḋ ≈ ċ₀/c ≲ 1.3×10⁻¹⁵ per year. The HfF⁺ oscillating-EDM search covers 27 nHz to 400 mHz (2006.15787), and a field with m ≲ H_* has frequency ≲ 10⁻¹⁴ Hz. For the lab this is a static EDM, which is what the JILA bound constrains.
- **No oscillation-specific bound applies.** The CFJ/χFF̃ background is P-odd and T-even, so at linear order it gives a CPT-odd b_μ, not an EDM.

## 6. A non-dynamical P'-odd background that keeps the KK photon massless and 4D Lorentz symmetry [G1-G4]

Take a fixed vector field V = ∂_y, i.e. a fixed orientation of the S¹ fibre, as prior geometry. Then u = V♭/|V| = φ(dy + kA) = e⁵ [G1].

- **KK gauge.** Under y → y + λ(x), V keeps its components [G2]. So every term built from V and the dynamical fields is KK-gauge invariant, and there is no Stückelberg mass at any order.
- **Contrast with a winding scalar.** A dynamical winding scalar Y = wy/R shifts under the same transformation, which gives the 012 U6 mass.
- **4D Lorentz.** 4D diffeomorphisms preserve V [G3].
- **P'.** P' maps V → -V [G3].
- **Broken generators.** Exactly the y-dependent diffeomorphisms are broken [G4]. These act only on massive KK modes.

The cost is explicit breaking of 5D diffeomorphism invariance. Kostelecký (hep-th/0312310) shows explicit breaking is generically in conflict with the Bianchi identities in Riemann-Cartan geometry. Here the conflict could only appear in y-dependent sectors, which is not checked.

In the zero-mode EFT the background reduces to a fixed sign of the c(χ)³ χ-FF̃ coupling. The 4D theory is then just a parity-violating scalar-photon coupling, and torsion supplies its functional form. I found no reference that uses ∂_y as an explicit background for this purpose. 2608.09299 gets the orientation from kinetic-less Stückelberg scalars instead.

## 7. CS level [C1-C6]

- With charge-normalised A (∫F/2π ∈ Z), κ_CS ∫_6 F³ = 2π(k/6)∫c₁³ [C1].
- On spin 6-manifolds the index theorem gives ∫(c₁³/6 - c₁p₁/24) ∈ Z. CP³ with O(n) gives index (n³ - n)/6. This is checked against Hirzebruch-Riemann-Roch, χ(O(n-2)) = (n+1 choose 3) [C2].
- CP³ is spin and has ∫c₁³ = 1 for O(1). So:
  - **without** the mixed term the minimal level is **k ∈ 6Z**;
  - **with** the mixed term -(k/192π²)A∧tr(R∧R) [C4], **k ∈ Z** [C3]. Products of S² never see the obstruction [C3b].
- After reduction, ϑ has period 2π for |k| = 1, because ∫c₁² is even on spin 4-manifolds [C5].
- A massive charge-q bulk fermion shifts the level by ½q³ sign(M). k_eff must be an integer, so k_bare ∈ Z + ½Σq³ (Intriligator-Morrison-Seiberg hep-th/9702198, Witten hep-th/9603150; standard, not opened).
- The charged bulk fermion that KK-02 needs for a Hosotani mass also induces the gravitational CS term at one loop. Once it is present, the k ∈ Z (spin) case is the consistent one.

β_max = |k|α/2 for Δϑ ≤ π:
- **minimal |k| = 1: 0.209°**, with Δϑ = 1.03π for ACT's central value;
- if only 6Z: 1.25°, with Δϑ = 0.17π [C6].
