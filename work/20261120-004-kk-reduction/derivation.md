# Einstein-Cartan-Kaluza-Klein on S¹ and S¹/Z₂: 4D fields and couplings to F

Task 20261120-004-kk-reduction. Conventions: `CONVENTIONS.md` (5D signature -,+,+,+,+; Γ^λ_{μν} with the derivative index μ first; ε_{0123} = +√-g; F̃^{μν} = ½ε^{μνρσ}F_{ρσ}; CS baseline L ⊃ -(g/4) φ F F̃). Tangent indices A,B = 0..3, 5̂; 4D tangent a,b. Check numbers [Cn] refer to lines printed by `run_checks.py`.

## 0. Contortion formula (convention issue found)

With ∇_μV^λ = ∂_μV^λ + Γ^λ_{μν}V^ν and Γ = {} + K, metric compatibility ∇_μ g_{λν} = 0 gives K_{λμν} = -K_{νμλ} (antisymmetric in 1st and 3rd index), and T_{λμν} = K_{λμν} - K_{λνμ}. The solution is

  K_{λμν} = ½ (T_{λμν} + T_{νλμ} - T_{μνλ}).

The formula in `CONVENTIONS.md`, ½(T_{λμν} - T_{μλν} - T_{νλμ}), reproduces T but is antisymmetric in (λ,μ), i.e. it belongs to the convention where the derivative index is last. With the stated index order it is not metric-compatible [C: "INFO CONVENTIONS.md formula", violation O(1)]. Both formulas agree on totally antisymmetric torsion (K = ½T), so the axial sector used below is unaffected. Reported for the overseer to fix.

## 1. Action and compactification

5D Einstein-Cartan plus a Dirac fermion as spin source (first-order, vielbein e^A_M, Lorentz connection ω^{AB}_M independent):

  S = (1/2κ₅) ∫ d⁵x e R(e, ω) + ∫ d⁵x e [ -½(ψ̄Γ^M D_Mψ - \overline{D_Mψ} Γ^M ψ) - m ψ̄ψ ],  D_M = ∂_M + ¼ ω_{MAB}Γ^{AB}.

Split ω = ω̊(e) + K (Levi-Civita spin connection plus contortion). Then, for any metric and any K field [C: "R(LC+K) = R_LC + div V + Q(K)", D = 5 and 4]:

  R(ω) = R̊ + ∇̊_μ V^μ + Q(K),
  V^μ = g^{σν}K^μ_{νσ} - g^{μσ}K^ν_{νσ},
  Q(K) = g^{σν}(K^μ_{μλ}K^λ_{νσ} - K^μ_{νλ}K^λ_{μσ}).

Q contains no derivatives and no ω̊. For totally antisymmetric K, Q = -K_{ABC}K^{ABC} [C].

The Hermitian Dirac action couples only to the totally antisymmetric part of the connection, because {Γ^C, Γ^{AB}} = 2Γ^{CAB} [C]:

  L_ψ ⊃ (i/4) ω_{[CAB]} ψ̄Γ^{CAB}ψ ≡ ¼ (ω̊_{[CAB]} + K_{[CAB]}) J^{CAB},

with J^{ABC} the real bilinear ∝ ψ̄Γ^{ABC}ψ. In 5D, Γ^{ABC} = c ε^{ABCDE}Γ_{DE}/2 with c = -i in our representation [C], so the source is the dual of the 2-form bilinear ψ̄Γ_{DE}ψ. There is no independent "axial current" in 5D; Γ^5̂ is part of the Clifford algebra.

A Weyssenhoff spin fluid gives the same structure (algebraic torsion ∝ spin density); nothing below depends on the source being a Dirac field except the counting of which components are sourced.

Compactification: y ∈ [0, 2πR), zero modes only (∂_y = 0). Metric ansatz, Einstein frame, radion kept:

  ds₅² = e^{2αφ} g_{μν}dx^μdx^ν + e^{2βφ}(dy + k A_μ dx^μ)²,  α = 1/(2√3), β = -2α.

Vielbein e^a = e^{αφ}ē^a, e^5̂ = e^{βφ}(dy + kA). Physical radius L = R e^{βφ}. With κ₄ = κ₅/(2πR) and k = √(2κ₄), the gravity sector reduces to [C: "KK Ricci reduction", residual 3e-15, total-derivative coefficient fitted = -2α]

  √-ĝ R̂ = √-g [ R - ½(∂φ)² - ¼ e^{-√3 φ} k²F² - 2α □φ ],

so S₄ = ∫√-g [ (1/2κ₄)(R - ½(∂φ)²) - ¼ e^{-√3φ} F_{μν}F^{μν} ] + (torsion and fermion terms). φ here is dimensionless; the canonical radion is σ = φ/√(2κ₄) = M_P φ/√2.

## 2. 4D fields from torsion, elimination of torsion, couplings to F

### 2.1 Counting

T^A_{BC} in 5D: 5 × 10 = 50 = 5 (trace) + 10 (totally antisymmetric 3-form) + 35 (pure tensor) [C: ranks 5, 10, joint 15]. Zero modes with tangent indices (using e^5̂ avoids mixing with A_μ):

| block | # | 4D content |
|---|---|---|
| T^a_{bc} | 24 | trace vector v_a (4), axial vector S^a (4), tensor q^a_{bc} (16) |
| T^a_{b5̂} | 16 | scalar s = T^a_{a5̂} (1), 2-form t_{[ab]} (6), symmetric traceless (9) |
| T^5̂_{ab} | 6 | 2-form u_{ab} |
| T^5̂_{a5̂} | 4 | vector w_a |

The 5D 3-form part (10) splits as T_{abc} ↔ axial vector S^a (4) and T_{ab5̂} (6, a 2-form). No zero-mode component carries four 4D indices, so there is no 4D pseudoscalar linear in torsion. Parity-odd scalars from torsion are at least bilinear (ε^{abcd}t_{ab}u_{cd}, S·v, S·w, ...).

### 2.2 (a) No source: no torsion

Q is non-degenerate on the full 50-dimensional space of metric-compatible K [C: rank 50/50], and the ∇̊V term is a total derivative (it stays one for y-independent zero modes). So δS/δK = 0 with no spin source gives K = 0. Pure 5D EC = 5D GR; its S¹ reduction is ordinary KK (graviton, KK photon, radion), parity-even.

### 2.3 (b) With a 5D Dirac source

Only K_{[ABC]} (10 components) is sourced; the other 40 vanish by non-degeneracy. Varying

  e [ -(1/2κ₅) K_{ABC}K^{ABC} + ¼ K_{ABC}J^{ABC} ]

gives K_{ABC} = (κ₅/4) J_{ABC}, i.e. T_{ABC} = (κ₅/2) J_{ABC}, and

  L_eff = (κ₅/32) J_{ABC}J^{ABC}.

Zero mode ψ = χ(x)/√(2πR) (radion frozen at L = R; radion factors otherwise multiply each term):

  L_4f = (κ₄/32) [ J_{abc}J^{abc} + 3 J_{ab5̂}J^{ab5̂} ] = -(3κ₄/16) A_aA^a + (3κ₄/32) J_{ab5̂}J^{ab5̂},

using J_{abc} = ε_{abcd}A^d and ε_{abcd}ε^{abce} = -6δ_d^e [C]. A^a is the 4D axial current. The first term has the magnitude 3κ/16 of the standard 4D Einstein-Cartan (Hehl-Datta) contact term, which fixes the normalization of J. The second is a KK-specific tensor-tensor contact term.

Couplings of F in the 4D effective theory:
1. Maxwell: -¼ e^{-√3φ} F² (radion-dependent gauge coupling).
2. Pauli term from the fermion's Levi-Civita connection. For the KK vielbein the totally antisymmetric part of ω̊ is ω̊_{[ab5̂]} = (k/6) F_{ab} (times a radion factor), and ω̊_{[abc]} = 0 [C: "KK Levi-Civita spin connection", residual 6e-17]. So L ⊃ (ik/8) F_{ab} χ̄Γ^{ab5̂}χ: a Planck-suppressed (k = √2/M_P) dipole term, linear in F.
3. Minimal coupling with charge n k/R for KK modes n ≠ 0. The zero mode is neutral.

There is no torsion-F term. F sits inside ω̊, and K enters e R(ω) only through Q(K), which has no ω̊, and the total divergence ∇̊V. Eliminating K therefore produces J² terms only. The fermion axial current couples to F only through the Pauli term, which is fixed by Levi-Civita geometry and does not involve torsion.

### 2.4 (c) Parity-odd photon terms at tree level: none

Birefringence needs a term quadratic in A_μ with a parity-odd background, either φ F F̃ or p_μ A_ν F̃^{μν}. At tree level:
- Gravity sector: R̊ reduces to the parity-even KK action above. Q(K) has no F.
- Fermion sector: F enters linearly (Pauli). A term linear in F with a background ⟨χ̄Γ^{ab5̂}χ⟩ is a magnetization source, not a change in photon propagation. For unpolarized cosmological matter that background also averages to zero.
- No pseudoscalar is linear in the torsion zero modes (2.1).

So pure EC-KK with a minimally coupled spin source gives no parity-odd photon coupling at tree level. The coefficient is zero for any κ₅, R and radion.

Possible one-loop loophole (open issue): integrating out charged fermions in an axial-torsion background S_μ induces a CFJ term whose coefficient is regularization-dependent (Jackiw-Kostelecký 1999, to verify). In KK the zero modes are KK-neutral, so only the charged tower (mass ~ n/L) contributes, and ⟨S_μ⟩ ∝ κ₄⟨A_μ⟩ vanishes for unpolarized matter.

A note on literature that says "torsion contains F" in KK: that statement refers to the teleparallel (Weitzenböck) connection, whose torsion T^5̂ = de^5̂ = e^{βφ}kF. In Einstein-Cartan the torsion is the spin-sourced K, and F lives in ω̊. The two should not be mixed.

### 2.5 Orbifold S¹/Z₂ (y → -y)

A component with n₅ indices along y picks up (-1)^{n₅} (times an intrinsic sign).
- Metric: g_{μν} even, g_{μ5} odd, g_{55} even. The KK photon zero mode is projected out and the radion survives.
- Torsion: even (survive) T^a_{bc} (v, S, q: 24) and T^5̂_{a5̂} (w: 4). Odd: T^a_{b5̂} (16) and T^5̂_{ab} (6).
- Fermion: one chirality survives, so J_{ab5̂} (an L-R bilinear) vanishes and A^a survives. The 4D EC result -(3κ₄/16)A² is recovered.

So EC on S¹/Z₂ is 4D EC plus a radion, with no geometric photon. The photon has to be a bulk or brane gauge field, and the 4D no-go (task 003) applies to it unchanged.

## 3. Minimal added terms that give a parity-odd photon coupling

### 3.0 No geometric parity-odd 5D term exists

The Lorentz-scalar 5-forms built from e (1-form), T (2-form) and R (2-form) without ε are only e_A∧R^{AB}∧T_B. A 3-e term would need a 2-form with three indices, and there is none. Using the Bianchi identity DT^A = R^A{}_B∧e^B,

  d(T_A∧T^A) = 2 R_{AB}∧e^B∧T^A  [C: exterior-algebra check, exact zero at random points],

so that term is exact. In 5D there is no Holst term (Holst e∧e∧R is a 4-form). The Nieh-Yan 4-form d(e^A∧T_A) needs a 1-form multiplier, i.e. a new vector. The gravitational CS term vanishes for so(4,1), since the symmetric trace of an odd power vanishes. So any parity-odd photon coupling needs a new 5D field. This is a no-go for pure EC-KK.

### 3.1 Candidate mechanism A: bulk U(1) with a 5D Chern-Simons term (Wilson-line axion)

  S ⊃ ∫ √-ĝ [ -(1/4g₅²) 𝓕_{MN}𝓕^{MN} ] + κ_CS ∫ 𝒜∧𝓕∧𝓕   (orientation d⁴x∧dy > 0).

The zero-mode ansatz, gauge-covariant under the KK U(1), is 𝒜 = a_μdx^μ + θ(x)(dy + kA_μdx^μ). Reduction, modulo total derivatives [C: Euler-Lagrange equivalence for all 9 fields; a control with coefficient 2 fails]:

  κ_CS ∫𝒜∧𝓕∧𝓕 = 2πR κ_CS ∫ [ 3θ f∧f + 3kθ² f∧F + k²θ³ F∧F ],

where f = da is the bulk photon and F = dA is the KK photon. With f∧f = -½√-g f_{μν}f̃^{μν} d⁴x [C]:

If the observed photon is the bulk U(1) (the realistic case, since SM fermions are KK-neutral zero modes):
- g₄² = g₅²/(2πL), so 1/g₄² ∝ e^{-φ/√3}.
- Axion: ϑ = 2πRθ = ∮𝒜, period 2π. Canonical φ_a = f_a ϑ with f_a = 1/(2πL g₄) (radion frozen at L).
- L ⊃ -(g_{φγ}/4) φ_a F_γF̃_γ with

  **g_{φγ} = 6 κ_CS g₄² / f_a = 12π κ_CS g₄³ L = 24π κ_CS α_em / f_a.**

  For a quantized level κ_CS = k_CS/(24π²): g_{φγ} = k_CS α_em/(π f_a). The sign follows the sign of κ_CS.
- The 4D coupling is radiatively stable, since it comes from a topological 5D term. ϑ gets a Hosotani potential from charged bulk fields, which is the extra-natural axion setup (Arkani-Hamed et al. 2003, hep-th/0301218, to verify). Mass and f_a are left to the cosmology task.
- Orbifold: CS invariance needs 𝒜_μ odd and 𝒜_5 even. Then ϑ survives but the bulk photon zero mode does not, and the photon must be on a brane, with brane-localized ϑ F F̃ couplings (model-dependent). The other choice, 𝒜_μ even with κ_CS → κ_CS ε(y), projects ϑ out. On S¹/Z₂ mechanism A needs extra brane input.

If the observed photon is the KK photon (only for completeness, since the fermion zero modes would then be neutral): with background θ₀, the k²θ³ term gives g_{φγ} = 6κ_CS k² θ₀²/f_a = 12κ_CS θ₀²/(M_P² f_a). That is suppressed by (θ₀/M_P)², plus mixing through 3kθ₀² f∧F and the kinetic mixing ∝ kθ₀ from 𝓕². Not pursued.

### 3.2 Candidate mechanism B: Kalb-Ramond 2-form identified with axial torsion

Impose T_{ABC} = H_{ABC} = (dB)_{ABC} (totally antisymmetric). Then K = ½H and the EC term gives -(1/8κ₅)H² = -(1/12g_{B5}²)H², so EC fixes g_{B5}² = 2κ₅/3 and g_B² = 2κ₄/3 in 4D. The fermion gets the derivative coupling (1/8)H·J. Reduction: B̂ = B₂ + B₁∧(dy + kA), H₃ = dB₂ - kB₁∧F, G = dB₁. The 4D dual of B₂ is an axion b (it survives on S¹/Z₂ if B_{μν} is even).
- Without extra terms: dH₃ = -kG∧F, so the only coupling is the mixed b G F̃ term with g_mix = k = √(2κ₄) (normalization L ⊃ -(g_mix/2) b G_{μν}F̃^{μν}). A mixed coupling with an unpopulated second vector gives no rotation, only an amplitude cos Φ [C: "mixed coupling"]. So β = 0.
- With a 5D Green-Schwarz term H = dB + λ 𝒜∧𝓕 for the photon's bulk U(1): dualizing (S = ∫ -(1/2g_B²)H∧*H - σ(dH - λ f∧f), canonical b = g_Bσ) gives

  **g_{bγ} = 2|λ| g₄²/g_B = √6 |λ| g₄² M_P**  (EC-fixed g_B),

  which is Planck-suppressed for natural λ ~ κ₄: g_{bγ} ~ 4πα_em/M_P.
- In 5D a 2-form is dual to a vector (Ĥ = *d𝒞), and B_{μν} ↔ 𝒞_5. A GS term λ𝒜∧𝓕 in H then dualizes to a 𝒞∧𝓕∧𝓕 coupling, so B is mechanism A with ϑ = ∮𝒞 and f_a set by the EC-fixed g_B. This equivalence is expected from 5D duality and has not been checked symbolically (open issue).

### 3.3 Holst/Nieh-Yan-like terms

No 5D Holst term exists, and the unique ε-free term e∧R∧T is exact (3.0). Coupling the 4D Nieh-Yan form to a Wilson line, ∫𝒞∧d(e^A∧T_A) ≅ ∫𝒞_5 N_4 dy, needs the new vector 𝒞. N_4 contains T_A∧T^A and e∧e∧R, not F, so it is not a photon coupling at tree level. Dead as a direct birefringence source.

## 4. Fine-structure constant and the radion

With σ the canonical radion (φ = √2 σ/M_P) and L = R e^{-φ/√3}:
- Bulk-U(1) photon (mechanism A/B): α_em ∝ 1/L ∝ e^{φ/√3}, so **δα/α = -δL/L = δφ/√3 = √(2/3) δσ/M_P**.
- KK photon (pure KK identification): α_em ∝ e^{√3φ} ∝ L^{-3} in 4D Planck units, so **δα/α = -3δL/L = √6 δσ/M_P**.

f_a also depends on the radion. In Einstein frame the ϑ kinetic term carries √-ĝ ĝ^{55}ĝ^{μν} = e^{4αφ}, so f_a ∝ e^{2αφ} at fixed g₅ and R, which makes g_{φγ} radion-dependent too. The formulas in 3.1 are evaluated at a stabilized L = R. The numerical bound is left to the bounds task.
