# 20261120-026 skeptic: EC-06 (tasks 022, 023, 024)

ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything below is about candidate mechanisms. No `RESULT.md`.

## Question
Do the EC-06 verdicts and numbers of 022 (theory), 023 (cosmology) and 024 (literature, `literature/ec06_prior_art.md`, `literature/novelty_audit.md`, `data/kaf_lowz_bounds.csv`) survive an independent check?

## Method
- `skep_beta.py` (own code, loads only `data/`; no import from 022/023): β for p ∝ H by quadrature in z and on an a-grid; reionization share (point, z = 10, visibility-weighted with my own tanh x_e, H only and with He as a sensitivity); f(T) ξ; trace vs ln Q per-unit β; MZ GW bound with both GW170817 sides and the Λ-era coefficient of 2001.06373 Eq. 65; 1003.4823 k_AF conversion; A_CB under 023's kinematic assumption with a BBKS P(k) normalised to σ₈ (023 used Eisenstein-Hu).
- `skep_tegr.py` (sympy, own formulation in the orthonormal frame): Weitzenböck torsion of the flat FRW tetrad, of a locally boosted tetrad, of the perturbed diagonal tetrad, and of the perturbed tetrad plus the extra Lorentz scalar ζ. 10 PASS, 0 FAIL.
- Web: abs/HTML pages of 2001.06373 (ar5iv), 2508.10532 (v2 HTML), hep-ph/0611192, 1003.4823, 1710.05834, 0907.5159, 2609.10671, 2507.02362.

Rerun (repo root, each < 1 min, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-026-skeptic-ec06/skep_beta.py
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-026-skeptic-ec06/skep_tegr.py
```
Output: `run.log`.

## Result

| claim | status | reason |
|---|---|---|
| β = (ξε/2) ln(1+z*) for p_0̂ = ξεH; ξε = 1.073×10⁻³ | SURVIVES | ∫H dt = ln(1+z*) = 6.99478 by two routes; ξε = 1.0729×10⁻³. |
| Reionization share 0.23-0.31 (n = 0) | SURVIVES | point 0.309, z = 10 0.343, visibility 0.230 (0.228 with He). CLAIMS EC-06c quotes only "31%": should read 23-31%. |
| f(T): ε_T = 3, ξ = 3.6×10⁻⁴ | SURVIVES | T_τ = 3a′/a, S = 0 (sympy); ξ = 3.576×10⁻⁴. T_t = 3H holds exactly, so the formula needs no ΛCDM assumption. |
| Sign: ξ > 0 gives β > 0 (paper convention) | SURVIVES | T_τ > 0 and β = +½∫p_τdτ (CONVENTIONS) give β = (3ξ/2)ln(1+z*) > 0. The sign carries no information: ξ is a free sign. The comparison with 2508.10532 cannot be made: its F̃ uses ε_{μνρσ} with no stated sign of ε^{0123}. |
| f(T) trace matches 2508.10532 within 4% | SURVIVES | 1.5 ln(1+z*) = 10.492 vs ln E(z*) = 10.054: 4.4%. Post-reionization share 0.31 vs 0.27. Eqs. 14, 36, 53 (40.21α), 54 (8.23×10⁻⁵ < α < 1.78×10⁻⁴) confirmed on the v2 HTML. |
| (dT)_{τi} = ∂_i∂_τ(Φ+Ψ) | WEAKENED | Reproduced for the diagonal perturbed tetrad (δT_τ = -3Ψ′, δT_i = ∂_i(Φ-2Ψ)). In f(T) the scalar sector also carries the Lorentz scalar ζ (022 itself notes it). With ζ: (dT)_{τi} = ∂_i[∂_τ(Φ+Ψ) + ∇²ζ - ζ″]. The formula is incomplete. U(1) breaking survives only if the f(T) antisymmetric equation does not set ∇²ζ - ζ″ = -(Φ+Ψ)′; that was not checked. |
| TEGR: a local boost shifts T_μ by a non-closed piece and creates axial S | SURVIVES | sympy: δT_μ = (-∂_xη, -∂_τη, 0, 0), d(δT)_{τx} = ∂_x²η - ∂_τ²η ≠ 0, δS^z = 2∂_yη/(3a²); metric unchanged. A boost with η ∝ x·f(τ) also shifts the homogeneous T_τ, so the isotropic β is frame-dependent. Wording fix: with the covariant-TEGR "proper" (zero inertial) spin connection prescription the diagonal tetrad is selected and the TEGR row becomes the f(T) row numerically. "Dead" should read "not defined by TEGR dynamics; with an extra frame prescription it reduces to EC-06c". |
| MZ: c_T - 1 = γ²/18 on the matter tracker | SURVIVES | Eq. 72 of 2001.06373 (ar5iv): ω² = [1 + (1+3p/ρ)γ²/9]k², so γ²/18 at w = 0. At GW170817's z ≈ 0.01 the Λ-era leading term (Eq. 65, with κρ/Λ = Ω_m(1+z)³/Ω_Λ assumed) gives a coefficient of 0.057 instead of 0.056. The tracker relation applies there to within 4%. |
| MZ: GW170817 bound value, ε_S ≤ 1.5×10⁻⁷, ξ ≥ 7×10³ | WEAKENED (number) | c_T - 1 > 0 (GW faster than light), so the applicable side is +7×10⁻¹⁶, not 3×10⁻¹⁵ (1710.05834 abstract: "between -3e-15 and +7e-16"). This gives γ ≤ 1.12×10⁻⁷, ε_S ≤ 7.5×10⁻⁸ and **ξ ≥ 1.4×10⁴**. The kill gets stronger by 2×. |
| MZ verdict "dead" | WEAKENED | Nothing excludes ξ ≥ 10⁴. The observable is p = ξεH, and it is the same as in the live f(T) row. No theory contains ξ (022 says so). The torsion-photon vertex has strength ~ξ × (torsion fluctuation), and that fluctuation is Planck-suppressed, so it is not strongly coupled. EC-02 was dead because a bound (ξ < 9×10¹⁷) sat below its need. MZ has no such bound, only naturalness: ξ is 4×10⁷ times the f(T) value. Correct wording: "disfavoured: needs ξ ≥ 1.4×10⁴, unbounded but unnatural". Loophole not checked: the bound ties ε to γ only on the tracker. In the Λ era P ∝ a⁻³ (Eqs. 29-30), and if such a decaying mode also exists in the matter era, ε(z*) is unbounded by GW170817. |
| 023: A_CB ≤ 5×10⁻¹¹ | SURVIVES (number), WEAKENED (scope) | BBKS Limber gives max D_L = 4.6×10⁻¹¹ (L ≈ 770) vs 023's 4.9×10⁻¹¹. Assumption (A), torsion tracking the local expansion, is not what f(T) gives. On the diagonal tetrad δT_τ = -3Ψ′ (ISW-like, (3ξΔΨ/2)² ~ 10⁻¹⁷) and δT_i is a gradient (endpoint term only). The ζ contribution is unknown. Label the number "generic kinematic estimate", not an EC-06c prediction. Conclusion (no A_CB constraint) unchanged. |
| 1003.4823 conversion: \|k_AF\| < 1.65×10⁻³⁴ eV, ξε < 0.136 | SURVIVES | t_L(2.80) = 5.5075×10³² eV⁻¹, k = 5.2°/t_L = 1.648×10⁻³⁴ eV, ξε < 0.136. The 5.2° (\|θ\| + 2σ) is ~95% if ±2.2° is 1σ. The abstract states no CL. Data inconsistency: `kaf_lowz_bounds.csv` says "68 (Table 2 labelled 1σ)" while `ec06_cosmo_inputs.csv` says "not stated, 1σ assumed". Harmonize. |
| 024 spot checks (abs pages) | SURVIVES (5/5) | 2508.10532: nonmetricity parity violation, l(Q)F̃F, matches. hep-ph/0611192: gravitational-leptogenesis variant, "weak evidence" from WMAP3+B03, matches. 1003.4823: θ = -0.8 ± 2.2, 8 RGs z > 2, matches. 0907.5159: fixed p_μ, CS term in the gravity sector for consistency, matches. 2609.10671: axial vector needs a photon mass, matches. 2507.02362 is a black-hole/Riemann-Cartan paper and its abstract has no birefringence; the [FT] claim of a birefringence remark in the conclusions was not refuted (full text not opened). |

## Verdict
- No BROKEN items. The β formula, reionization share, f(T) ξ and sign, ln Q comparison, TEGR boost algebra, k_AF conversion and A_CB number all reproduce independently.
- Required fixes:
  1. MZ: use the +7×10⁻¹⁶ side, so ε_S ≤ 7.5×10⁻⁸ and ξ ≥ 1.4×10⁴. Change "dead" to "disfavoured (unnatural ξ, not excluded)", and list the tracker-only assumption.
  2. f(T) gauge breaking: add ζ, giving (dT)_{τi} = ∂_i[∂_τ(Φ+Ψ) + ∇²ζ - ζ″]. U(1) breaking is "expected", not shown, until ζ is solved from the antisymmetric f(T) equation.
  3. TEGR: reword as above (frame prescription, then = EC-06c).
  4. CLAIMS EC-06c: write 23-31%, not 31%.
  5. Label A_CB ≤ 5×10⁻¹¹ as generic (kinematic tracking).
  6. Harmonize the 1003.4823 CL between the two data files.

## Open issues
- MZ matter-era decaying mode (P ∝ a⁻³ branch) and whether it evades the γ-ε link: not checked.
- ζ in f(T) scalar perturbations (Golovnev-Koivisto 1808.05565): solve the antisymmetric equation to settle the U(1) breaking.
- 2001.06373 Eq. 65 normalisation of Λ assumed (Λ ↔ 3H₀²Ω_Λ); only the γ² term was used.
- 2507.02362 full text not opened.
