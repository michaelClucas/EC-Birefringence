# 20261120-004-kk-reduction: Einstein-Cartan-Kaluza-Klein on S¹ (and S¹/Z₂)

## Question
In a 5D Einstein-Cartan theory with a spin source, reduced on S¹ (KK photon = g_{μ5}), does the photon pick up a parity-odd coupling (φ F F̃ or CFJ p_μ A_ν F̃^{μν}) that could produce cosmic birefringence? If not, what is the minimal added 5D term that does, and what is the 4D g_{φγ}?

## Method
- Action: (1/2κ₅)∫e R(e,ω) plus a 5D Dirac fermion (Hermitian, minimal). Split ω = ω̊ + K. KK ansatz in Einstein frame with the radion kept: ds² = e^{2αφ}g + e^{2βφ}(dy + kA)², α = 1/(2√3), β = -2α. Zero modes only, then the orbifold projection.
- The full derivation is in `derivation.md`. Checks are numeric/symbolic, at random points or with exact algebra:
  - `check_kk_ricci.py`: 5D Ricci reduction, using 2-jets at random points.
  - `check_torsion.py`: contortion formula, the 50 = 5 + 10 + 35 counting, R(ω̊+K) = R̊ + ∇̊V + Q(K) for generic K, non-degeneracy of Q, 5D Clifford identities, and the Pauli term from ω̊.
  - `check_cs.py`: CS reduction by Euler-Lagrange equivalence in sympy, f∧f vs F F̃, and the no-rotation result for a mixed coupling.
  - `check_forms.py`: an exterior-algebra check that the only ε-free 5D geometric 5-form is exact.
- Run from the repo root: `.venv/Scripts/python tools/pin.py .venv/Scripts/python.exe work/20261120-004-kk-reduction/run_checks.py`. About 35 s on one core. Output is in `run_checks.log`.

## Result
- **Check results: 13/13 PASS** (`run_checks.log`).
- (a) With no spin source, torsion is zero. The quadratic form Q(K) is non-degenerate on all 50 components, so pure EC-KK is ordinary KK: √-ĝR̂ = √-g[R - ½(∂φ)² - ¼e^{-√3φ}F² - 2α□φ] (checked).
- (b) A Dirac source sees only K_{[ABC]} (10 of 50 components; in 5D this is the dual of the 2-form bilinear ψ̄Γ_{DE}ψ). Solving gives K = (κ₅/4)J and L_eff = (κ₅/32)J_{ABC}J^{ABC}, which in 4D is -(3κ₄/16)A_aA^a + (3κ₄/32)J_{ab5̂}J^{ab5̂}. The first term reproduces the 4D Hehl-Datta magnitude. There is **no torsion-F coupling**: F sits in ω̊, and K enters eR(ω) only through Q(K), which has no ω̊, and a total divergence. F couples to the fermion only through the Levi-Civita Pauli term (ik/8)F_{ab}χ̄Γ^{ab5̂}χ with k = √(2κ₄). That term is linear in F, so it does not change photon propagation.
- (c) **No parity-odd photon term at tree level.** The coefficient is zero for any κ₅, R and radion. No torsion zero mode is a 4D pseudoscalar linear in T. Also, the only ε-free 5D 5-form built from (e, T, R) is e∧R∧T = ½d(T∧T), which is exact (checked). There is no 5D Holst term, and gravitational CS vanishes for so(4,1). So a parity-odd photon coupling needs a new 5D field.
- Orbifold S¹/Z₂: the KK photon is projected out. The surviving torsion is T^a_{bc} (24) and T^5̂_{a5̂} (4), and with a chiral zero mode this is 4D EC plus a radion. The 4D no-go applies to whatever photon is added.
- Minimal additions:
  - **Candidate mechanism A**: a bulk U(1) with a 5D CS term κ_CS𝒜∧𝓕∧𝓕. On S¹ it reduces to 2πRκ_CS[3θ f∧f + 3kθ² f∧F_KK + k²θ³F_KK∧F_KK] (checked). For the bulk photon: **g_{φγ} = 6κ_CS g₄²/f_a = 24πκ_CS α_em/f_a**, with f_a = 1/(2πL g₄), g₄² = g₅²/(2πL), and ϑ = ∮𝒜 = φ_a/f_a. On S¹/Z₂ the CS parity rules keep ϑ but remove the bulk photon zero mode, so a brane photon is needed.
  - **Candidate mechanism B**: a KR 2-form with T = dB. EC fixes g_B² = 2κ₄/3. Without a Green-Schwarz term, the axion couples only as b G F̃ with g_mix = √(2κ₄), to the photon and a second vector, and that gives no rotation (checked). With H = dB + λ𝒜∧𝓕: **g_{bγ} = 2|λ|g₄²/g_B = √6|λ|g₄²M_P**, which is Planck-suppressed for λ ~ κ₄. By 5D duality (B ↔ vector) this should be mechanism A (not checked).
  - Holst/Nieh-Yan-like terms: there is no 5D Holst term and the 5D NY-type term is exact. Coupling the 4D NY form to a Wilson line needs a new vector and still does not involve F. Dead as a direct photon coupling.
- Radion and α: bulk photon **δα/α = -δL/L = √(2/3) δσ/M_P**; KK photon **δα/α = -3δL/L = √6 δσ/M_P** (σ is the canonical radion).

## Verdict
Pure EC-KK (S¹ or S¹/Z₂) with a minimally coupled spin source gives **no** parity-odd photon coupling at tree level, and in 5D there is no purely geometric parity-odd term to add. Birefringence needs a new 5D field with a 𝓕∧𝓕 coupling. The minimal choice is a bulk U(1) with a 5D CS term, whose Wilson line is the axion. This is candidate mechanism A, which behaves as an axion-like particle with g = 24πκ_CS α/f_a. Torsion plays no role in it. Mechanism B (KR = axial torsion) gives the same structure once a GS term is added, with an EC-fixed decay constant near M_P. As torsion-specific mechanisms, both are an ALP in disguise.

## Open issues
1. **CONVENTIONS.md error**: with Γ^λ_{μν} and the derivative index first, the contortion is K_{λμν} = ½(T_{λμν} + T_{νλμ} - T_{μνλ}). The listed ½(T_{λμν} - T_{μλν} - T_{νλμ}) is antisymmetric in (λ,μ) and not metric-compatible in this index order (check_torsion INFO line). The two agree for totally antisymmetric torsion. The overseer should fix it.
2. One-loop CFJ from charged KK modes in an axial-torsion background: the coefficient is regularization-dependent (Jackiw-Kostelecký 1999, to verify) and ⟨S_μ⟩ = 0 for unpolarized matter. Not computed.
3. Mechanism A on S¹/Z₂ needs brane-localized ϑ F F̃ (anomaly inflow). This is model-dependent and not derived.
4. KR ↔ Wilson-line duality with a GS term (mechanism B = A) was not symbolically checked. Signs of g_mix and g_{bγ} are convention-level and unchecked; magnitudes were derived by hand.
5. Radion-frame bookkeeping for f_a(φ) and g_{φγ}(φ) during radion evolution. Formulas are given at a stabilized radius.
6. The normalization of J (the phase of ψ̄Γ^{ABC}ψ in mostly-plus) is fixed by matching the 4D 3κ/16, not derived independently. The relative factor 3 between the J_{abc} and J_{ab5̂} terms follows from the index sum.
7. Literature for the scout to verify: extra-natural inflation/quintessence (hep-th/0301218); KR torsion in Randall-Sundrum (Mukhopadhyaya, Sen, SenGupta, PRL 89 (2002) 121101); KR-induced optical activity (SenGupta and collaborators); Jackiw-Kostelecký PRL 82 (1999) 3572.
