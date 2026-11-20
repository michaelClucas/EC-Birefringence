# 20261120-006 skeptic-iter1: audit of 003 (4D EC photon) and 004 (EC-KK reduction)

The ACT DR6 β = 0.215° ± 0.074° (arXiv:2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything audited here is a candidate mechanism or a no-go for one.

## Question

Do the iteration-1 results of `work/20261120-003-theory-ec-photon` (EC-00, EC-02, CFJ ↔ φFF̃ sign, β formulas, FRW torsion, contortion, S = -(κ/2)J₅, L_eff, the 28-orders estimate) and `work/20261120-004-kk-reduction` (KK-01, CS reduction, g_φγ, δα/α) survive an independent attempt to break them? Is "β > 0 clockwise as seen by the observer" and the map χ → β right?

## Method

1. Read both READMEs, both derivations and the check scripts, looking at whether each check could fail.
2. Re-derived the following independently, without importing any 003/004 code:
   - `indep_cfj_beta.py`: L = -¼F² - (g/4)φ(t)FF̃ with my own ε (ε^{0123} = -1). The EL equations come from sympy. A plane wave goes along a **random** direction k, not the z axis. The script checks the helicity eigenvalues, integrates the ODE numerically and reads the polarization at the observer in the **HEALPix/COSMO basis** (e_θ, e_φ) of the line of sight n = -k̂, with U = 2Re(a_θ a_φ*) = -U_IAU. It then takes β = Δψ_COSMO, which is the convention of 2202.13919 (Q ± iU = Pe^{±2iβ}) and of 2509.13654 Eq. 4. It checks separately that increasing ψ_COSMO is clockwise for an observer facing n with north up. As a control it flips ε to Komatsu's ε^{0123} = +1/√-g; that must reproduce his published β = +(g/2)Δφ.
   - `indep_torsion_estimates.py`: contortion (metric compatibility, with a control on the old formula), the Riemann convention on S² (R = +2), and R(K) for constant pure-axial and pure-trace torsion. It also covers the S = -(κ/2)J₅ and L_eff algebra, FRW S^τ = -f/a², and the size estimates with inputs from `data/` plus `inputs_skeptic.csv`, including a stress test at |ξ_ν| = 1. Last, it rederives the KK-02 f_a and g_φγ normalization, which 004 did by hand.
3. Checked references on the web: 2509.13654 (HTML) and 2202.13919 (HTML).

Run from the repo root (single core, about 1 min total):
```
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-006-skeptic-iter1/indep_cfj_beta.py
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-006-skeptic-iter1/indep_torsion_estimates.py
```
Result: 9/9 and 18/18 PASS (`run.log`). My first draft had a numpy transpose bug in the contortion. It failed the metric-compatibility test, which shows that test discriminates; "reproduces T" alone does not.

## Result

| claim | verdict | reason |
|---|---|---|
| Contortion K_{λμν} = ½(T_{λμν}+T_{νλμ}-T_{μνλ}) | SURVIVES | Independent numeric check: it reproduces T and is metric compatible. The old formula fails compatibility but agrees on totally antisymmetric T. |
| CFJ ↔ φFF̃: p_μ = -g∂_μφ | SURVIVES | 003 has a control (+g fails). My EL equations with my own ε give ω_±² = k² ∓ p_t k with p_t = -gφ̇ for helicity ± about k, along a random k. |
| CFJ dispersion ω_±² = k² ∓ p_τk (e_x ± ie_y, wave along +z) | SURVIVES | Reproduced as above, with helicity defined right-handed about k̂. |
| β = +½∫p_τdτ and β = -(g/2)(φ_0 - φ_LSS) (our ε) | SURVIVES | Independent route via the HEALPix e_θ/e_φ basis, not 003's "facing" cross product: numerical β = -0.2502 vs -0.2500 rad. The control with Komatsu's ε gives +0.2502, i.e. his published +(g/2)Δφ, so our sky mapping matches his end to end. |
| "β > 0 clockwise as seen by the observer"; χ right-handed about k → β = -χ | SURVIVES, with a sourcing caveat | 2509.13654 only says "clockwise (β>0)" and names no viewpoint or IAU/COSMO. The viewpoint comes from 2202.13919 (same group): β>0 clockwise, β = -ψ_IAU, z along the line of sight. A right-handed rotation about k (pointing at the observer) appears counterclockwise to an observer looking along -k, so β = -χ is right. CONVENTIONS.md attributes the convention to 2509.13654; it should cite 2202.13919 for the viewpoint. |
| FRW torsion: only T_τ = 3h and S^τ = -f/a² | SURVIVES | S^τ = -f/a² rederived. The two-function result matches the known homogeneous-isotropic torsion literature. |
| R(Γ) ⊃ -⅔T_μT^μ + (3/2)S_μS^μ | SURVIVES | Independent constant-torsion computation with the Riemann sign fixed on S². |
| S^μ = -(κ/2)J₅^μ, L_eff = -(3κ/16)J₅² | SURVIVES | Algebra rechecked. c_A = 3/4 is consistent with Shapiro's η = 1/8 coupling for S_Shap = εT = 6S_ours (literature cross-check; I did not recompute the Dirac traces). The sign of S depends on the c_A sign, which is irrelevant to β here. |
| EC-00 (minimal EC + minimal Maxwell: β = 0) | SURVIVES | F = dA has no connection, and Maxwell is conformally invariant on FRW. The result is trivial and standard. It holds only for F = dA (the gauge-breaking ∇-curl is (0'), handled separately; p = 2S rechecked by hand). |
| EC-02 (28 orders short) | WEAKENED (number fine, framing incomplete) | Rederived: 1.1×10⁻³⁵ eV constant, 2.7×10⁻⁴⁰ eV today for a⁻³ scaling, 2.0×10⁻⁶⁸ eV at ξ_ν = 0.1, giving 10^28.1. Robust: still 10^27.1 at an extreme |ξ_ν| = 1, and a⁻³ ignores the v suppression of J₅ for massive neutrinos, so it is generous. But "28 orders short at ξ = 1" is only a statement about ξ = 1. ξ is a free Wilson coefficient. The needed ξ ≈ 10²⁸ is the neutrino-current Chern-Simons operator J₅·AF̃/M² with M ≈ 3×10⁴ GeV, which is a published live mechanism (Geng-Ho-Ng 0711.4617; Zhou-Huang-Geng 2302.11140). So EC-02 is dead as a torsion mechanism (torsion only supplies M = M_P). The four-fermion phenomenology is not dead. Also, ξ_ν = 0.1 is still an uncited placeholder, and 003 keeps cosmology inputs in a task-local json instead of `data/`. |
| KK-01 (no parity-odd photon coupling, pure EC-KK, S¹ and S¹/Z₂, tree level) | SURVIVES (tree level, zero modes, minimal action) | The structure R(ω̊+K) = R̊ + div + Q(K), with no ω̊ in Q, is correct, and the Dirac source only sees K_{[ABC]}. On S¹, KK-momentum conservation forbids tree-level exchange of n ≠ 0 modes between zero modes. The 5-form counting is right (no ε-free 5-form other than e∧R∧T, which is exact; no degree-3 so(4,1) invariant). Limits: minimal 5D action only (no 5D Gauss-Bonnet with torsion, no quadratic-torsion PGT terms), and orbifold brane terms are not covered. |
| CS reduction coefficients 3, 3k, k² | SURVIVES | EL-equivalence check with a working control (coefficient 2 fails). |
| g_φγ = 6κ_CS g₄²/f_a = 24πκ_CS α/f_a, f_a = 1/(2πLg₄) | SURVIVES | 004 did this by hand; scripted here with f∧f = -½√-g ff̃. The quantization κ_CS = k/(24π²) → g = kα/(πf_a) is plausible but unchecked. |
| δα/α = √(2/3)δσ/M_P (bulk), √6 δσ/M_P (KK photon) | SURVIVES | Rederived by hand: bulk F² carries e^{βφ} → α ∝ 1/L; KK photon α ∝ e^{√3φ}; σ = M_Pφ/√2. |

## Verdict

No sign or factor error found in 003 or 004. All formulas survive an independent rederivation along a different route (arbitrary k, HEALPix/COSMO Stokes basis, Komatsu-ε control). Nothing must be fixed in the formulas. Two framing points need fixing:
1. EC-02 should be recorded as "minimal-EC torsion gives the neutrino-current CS operator with M = M_P, 28 orders too weak". It should not read as a no-go for neutrino-asymmetry birefringence: the same operator with M ~ 30 TeV is a published live mechanism with no torsion in it.
2. The "β > 0 clockwise as seen by the observer" viewpoint should be cited to 2202.13919, not to 2509.13654, which does not state a viewpoint.

## Open issues

1. Symmetry argument (mine, not scripted; should be checked by the KK theorist): in 5D, reflecting x¹, x², x³ and y together is a proper rotation (det +1), so it is an exact symmetry of any 5D Lorentz-invariant theory on a y-symmetric S¹. In 4D it acts as the effective parity. A radion × F_KK F̃_KK term is odd under it, so it is forbidden **to all loop orders** unless the action contains an orientation-odd (CS-type) term or there is a Y-odd background (Wilson line, chirality density). This strengthens KK-01 beyond tree level. It also predicts that arXiv:2608.09299 (teleparallel KK, radion × FF̃) must contain a parity-odd teleparallel term or such a background. Check which.
2. EC-02: cite a neutrino-degeneracy bound (e.g. EMPRESS arXiv:2203.09617, or a BBN+CMB fit) and move H0, z*, T_CMB and N_eff into `data/`. `inputs_skeptic.csv` here has the same problem; it is temporary.
3. ACT DR6 maps vs the Planck β sign: 2509.13654 states agreement in sign with WMAP/Planck, which implies a common convention, but the paper does not state its Stokes convention. Low risk, unverified.
4. Not audited: EC-05 loop normalization |p| = 2Ce²(3/4)|S| (the factor 2 depends on the CFJ normalization, ½p vs k_AF), KK-03 g_bγ, and the 4.47× k_AF discrepancy already in STATE.md (it is likely the same p = 2k_AF and √(4π) normalization question).
5. KK-01 assumes the minimal 5D action. 5D Gauss-Bonnet with torsion, or quadratic-torsion terms, are not covered; issue 1 suggests they still cannot give a parity-odd term without an orientation-odd piece.
6. I did not rerun 003/004's scripts (the overseer already did). I read them: 003 check 7 (the sky mapping) is sign bookkeeping and could not fail on physics. The numerical check 5 and my HEALPix route cover it.
