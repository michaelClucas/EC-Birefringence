# 20261120-022 theory EC-06 (4D track, iteration 6)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism. There is no `RESULT.md`.

## Question
Skeptic 020 blocker: can cosmological torsion of order H that is not sourced by spin, plus a non-minimal CFJ coupling -½ξ X_μA_νF̃^{μν} (X = axial S or trace T), give β? Three theories: (a) Magueijo-Złośnik EC with varying Λ (1908.05184), (b) teleparallel TEGR / f(T), (c) PGT torsion cosmologies. For each: background (minimal EC or extension, components, size vs H), gauge invariance on and off FRW, β formula and the post-reionization share, independent bounds, verdict.

## Method
sympy derivations in our conventions, plus numbers from `data/`. Write-up: `derivation_ec06.md`. Rerun (about 1.5 min, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-022-theory-ec06/run_checks.py
```
Last run: `run.log`, **54 PASS, 0 FAIL**.

| script | content |
|---|---|
| check_gauge.py | CFJ gauge variation ∝ χ(dp)·F̃; field-equation divergence = -¼(dp)·F̃ (inconsistent unless dp = 0); closure of FRW and perturbed p; Stückelberg mixing |
| check_teleparallel.py | flat FRW Weitzenböck torsion; TEGR local-boost dependence; k = +1 S³ tetrads (axial ±2√k/a); f(T) good-tetrad equation; scalar-perturbation closure |
| check_mz.py | MZ torsion forms → our S, T; γ → ∞ and finite-γ algebra; tracker R_c(γ); expansion-rate deficit |
| ec06_numbers.py | β, ξ needed, reionization share, N_eff and GW170817 limits on the MZ tracker, Ω_K limit on k = +1, lab / k_AF / radio-galaxy comparisons, comparison with the published ln Q / ln R forms; writes `ec06_backgrounds.csv` |

New data: `data/ec06_theory_inputs.csv` (Ω_K, 1807.06209; GW170817 speed, 1710.05834). Inputs from 023: `data/ec06_cosmo_inputs.csv` (z_re, radio-galaxy UV rotation).

## Result
Define ε_X = X_t/H (lower cosmic-time component). Then **β = (ξ/2)∫ε_X dz/(1+z)**. It is achromatic. With constant ε, ξε = 1.07×10⁻³ and 31% of β accrues after z_re = 7.67.

| theory | background | components (ε) | p closed? | ξ needed | post-reion share |
|---|---|---|---|---|---|
| (a) MZ, γ → ∞ | extension of EC (quasi-Euler 1/Λ) | ε_S = 1.73, ε_T = -1.5 on tracker | no (P is a new d.o.f.) | — | — |
| (a) MZ, finite γ (Pontryagin) | extension of EC (+ quasi-Pontryagin) | tracker P/H = γ/3, ε_S ≈ 2γ/3. N_eff alone: γ ≤ 0.81. **GW170817 tensor speed** (2001.06373 Eq. 72, c_T - 1 = γ²/18): γ ≤ 2.3×10⁻⁷, ε_S ≤ 1.5×10⁻⁷; Λ era ε_S(0) ≤ 1.2×10⁻⁸ | FRW yes, off-FRW no | **≥ 7×10³** (5×10⁴ with the paper's γ² ≲ 10⁻¹⁵; ≫ 10¹⁸ with its Eq. 70) | — |
| (b) TEGR | GR-equivalent; torsion is frame gauge | anything (local boost shifts T, creates S) | — | undefined | — |
| (b) f(T), flat | f(T), good tetrad | ε_T = 3, S = 0 | FRW yes; off-FRW (dT)_{τi} = ∂_i∂_τ(Φ+Ψ) ≠ 0 | 3.6×10⁻⁴ (P-odd coupling) | 0.31 |
| (b) f(T)/teleparallel, k = +1 | closed FRW axial branch | ε_S = ±2√\|Ω_K\|(1+z)/E(z): ≤ 0.11 today, ≤ 5×10⁻³ at z* | FRW yes, off-FRW no | ≥ 2.2×10⁻² (\|Ω_K\| ≤ 0.0030) | 0.65 |
| (c) PGT 0⁻ | PGT | S = ∂φ/√A | yes | = EC-04 (ALP, ~10¹²⁰ tuning) | — |
| (c) PGT 0⁺ (SNY) | PGT, trace only | T_t = Φ = -(3/2)d ln(R+6μ/b)/dt, Φ/H ≈ 0.8-1.4 today, z* unknown | FRW yes (gradient); off-FRW expected yes, not checked | ~10⁻³ if ε ~ 1 | unknown |

- Independent bounds do not bind any row. Lab b_μ (|S_0| < 7.3×10⁻²² eV) needs ε > 10¹¹. The radio-galaxy k_AF bound is 2.5 orders above p₀ ~ 10⁻³⁶ eV. The UV radio-galaxy rotation at z > 2 allows ±2.2°, and the models give 0.03°. Solar-system torsion effects are ~10⁻⁸ξ. ξ is bounded by nothing in these theories, because none of them contains it.
- Gauge-invariant higher-dimension torsion-photon operators ((∇·S)FF̃/M², F²T²) need M ≈ 5×10⁻²⁸ eV, only 14 H*. Dead.
- Non-closed p makes Maxwell + CFJ inconsistent ((dp)·F̃ = 0 forced on solutions) unless a Stückelberg/massive photon is added.

## Verdict
- **(a) MZ: dead.** The γ → ∞ (Euler-only) branch: 33% expansion-rate deficit on the tracker, MZ's own "unviable", undetermined tensor modes (2001.06373). The finite-γ (Pontryagin) branch: the same parity-odd torsion makes gravitational waves travel at c_T - 1 = γ²/18 on the tracker. Our tracker x = P/g ≈ γ/3 agrees with 2001.06373's P² = γ²κρ/27. GW170817 then gives ε_S ≤ 1.5×10⁻⁷ and ξ ≥ 7×10³. That is 2×10⁷ times the f(T) value, and far outside a perturbative coupling. The scout's estimate (ε ≲ 2×10⁻⁸, ξ ≳ 5×10⁴) is reproduced with the paper's own γ² ≲ 10⁻¹⁵. The torsion is also not closed off FRW.
- **(b) TEGR: dead (ill-defined).** Its torsion is pure local-Lorentz gauge, and β depends on an unfixed frame. **f(T): candidate mechanism** with the cleanest numbers (ξ = 3.6×10⁻⁴ on the trace, or ξ ≥ 0.022 on the k = +1 axial branch, which needs Ω_K < 0 near the current limit). Costs: P-odd coupling for the trace; U(1) broken at the perturbation level wherever Φ+Ψ evolves; f(T)'s own perturbation problems.
- **(c) PGT:** 0⁻ = EC-04 (ALP in disguise, tuned). **0⁺ (SNY) is the only case with a closed p** (on FRW, and expected off FRW). The coupling is then -(3ξ/4)ln(R+6μ/b)FF̃, an ALP-like coupling of the scalar-torsion mode: gauge invariant, P-odd, background history disputed and unknown at z*.
- **Prior art (task 024):** the flat f(T) trace coupling p = ξT = 3ξ∂ln a has the same phenomenology as Yu-Yang-Qiu 2508.10532 (ln Q FF̃, Q = 6H²) and Li-Xia-Li hep-ph/0611192 (∂ ln R). In matter domination ∂ln Q = ∂ln R = -3∂ln a, so they coincide with c = -ξ. Over the full history β differs by 4% from ln Q (N7) and 12% from ln R, with a post-reionization share of 0.31 vs 0.27 for ln Q. On flat FRW, Q = 6H² is also the f(T) torsion scalar (T4), so the gauge-invariant completion of the trace-CFJ coupling is **exactly** the published c ln(T_scalar)FF̃ form. EC-06 (b, flat) is therefore not new phenomenology. New here: the torsion reading, the TEGR frame-dependence no-go, and the off-FRW failure of the CFJ form. The k = +1 axial branch (conformal-distance weighting, 65% post-reionization) is a different β(z) shape, and the scout found no prior art for it. Off-FRW costs of a non-closed p are a gravitational CS counterterm (0907.5159) or a photon mass (2609.10671); these were not opened here, the citations come via 024. Teleparallel Nieh-Yan (2007.08038) is GW-only.
- **For Result C** the cost list is common to all live rows: (1) the photon coupling ξ ~ 10⁻³-10⁻² is added by hand, and no theory predicts it; (2) gauge invariance holds only if the torsion vector is a gradient, which holds on FRW but fails off FRW for f(T) (it holds only for PGT 0±); (3) ξε ≈ 10⁻³ is a free choice, not tuned against anything but not explained either; (4) theory choice: the background is not minimal EC in any case (f(T), or PGT with R² terms). Where p is closed, the mechanism is an ALP whose field is a torsion scalar.

## Open issues
1. MZ: the GW bound uses 2001.06373 Eqs. 65-72 as published (not rederived). It applies the matter-tracker relation c_T - 1 = γ²/18 at GW170817's z ≈ 0.01 (Λ era); the separate Λ-era term 42ε_P² gives ε_S(0) ≤ 1.2×10⁻⁸ independently of γ. Off-FRW closure of S is not derived (moot once dead).
2. Quantify the Stückelberg cost: longitudinal-photon mixing ∝ ξ k Φ̇/m_γ vs photon-mass bounds (not in data/).
3. SNY 0⁺: verify T_μ ∝ ∂_μ ln(R + 6μ/b) covariantly (needs the full connection equation), confirm the sign of Φ vs our T_t, integrate the SNY background back to z*.
4. MZ tracker deficit sign: my H = g - T gives a deficit, while MZ describe extra radiation. Both are excluded at γ → ∞; at finite γ the N_eff bound uses |deficit| ≤ 5.3%. The matter-era change of effective G at recombination (same 1 - γ²/9) may be bounded more tightly by the CMB than by N_eff.
5. f(T) strong-coupling/extra-mode literature (e.g. Golovnev-Koivisto 1808.05565, Blagojević-Nester 2006.15303) not opened here.
6. Anisotropic β from torsion perturbations vs SPT A_CB, and the β(z) shapes in `ec06_backgrounds.csv`, are for task 023.
7. GW-speed test for the live rows: the MZ kill came from parity-odd torsion modifying tensor propagation. f(T) on flat FRW has c_T = 1 (literature, not rechecked). The k = +1 axial branch and the PGT 0⁺/0⁻ backgrounds were not checked for c_T or GW birefringence.
8. Minkevich PGT cosmologies and the k = -1 teleparallel branches not computed.
