# 20261120-011 skeptic on iteration 2 (tasks 007, 008, 009)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats. Everything here is a candidate mechanism.

## Question
Do the iteration-2 results survive an adversarial audit? In particular:
- the EC-02 numbers;
- whether a shift-symmetric ∂φ·J₅ torsion coupling gives a QCD mass and a photon coupling (008 readings A vs B);
- the 008 washout and ODE;
- the P' all-orders argument;
- the 2608.09299 description;
- the Cassini/ISL kill of KK-04;
- the KK-02 and KK-03 windows.

## Method
- Read the READMEs, derivations, scripts and logs of 007, 008 and 009. Read 2608.09299 (arXiv HTML) and compared it with `literature/2608.09299.md`.
- Wrote three scripts of my own. None imports any work/ code, and every measured input is loaded from `data/`:
  - `s1_ec02_independent.py`: own ΛCDM integral of n₅ ∝ a⁻³ from LSS to today; massless-ν envelope; Reig ξ bound.
  - `s2_shift_symmetry.py`:
    - two-flavour chiral toy with the GG̃/η′ term. It compares reading A (keep only the anomaly piece of ∂·J₅) with reading B (the full field-dependent chiral rotation, which moves φ into both θ_G and arg m_q).
    - φ-π⁰ kinetic mixing.
    - B₁ computed from the Feynman-parameter integral directly, not from the closed form.
    - Each test has a control that fails when the physics is wrong: reading A must reproduce the textbook axion mass and does.
  - `s3_pprime.py`:
    - P' as exp[π(J₁₂+J₃y)] with its spinor lift.
    - Pull-back of the KK metric under the diffeomorphism.
    - Parity of φᵖFF̃ and φᵖF².
    - A bulk-U(1) control: θff̃ must come out even.
    - A schematic KK-fermion-tower probe with and without a Wilson line.
- Rerun: `.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-011-skeptic-iter2/<script>` (each script under 1 min, 1 thread). Log: `run.log`, **25 PASS, 0 FAIL**.

## Result

| claim | verdict | reason |
|---|---|---|
| EC-02 numbers (7.5×10⁻³⁰ deg, 28.5 orders, ξ needed 2.9×10²⁸; generous 27.5) | **SURVIVES** | Own integral gives 7.47×10⁻³⁰ deg, ξ = 2.88×10²⁸ (generous case: 6.1×10⁻²⁹ deg, 27.55 orders). n₅ = helicity × speed is correct; the axial current is C-even, so ν and ν̄ add with their helicities. |
| EC-02 Reig bound ξ < 9×10¹⁷ | SURVIVES (completion-dependent) | Reproduced as 9.40×10¹⁷. It holds only for Reig's KR/dual-axion completion. The death of EC-02 rests on the 28-order size, not on this bound. |
| EC-05 (r, ≥ 29 orders, (H/m_e)² decoupling) | **SURVIVES** | r_JK = 9αΣQ²/8π re-derived by hand (|b| = ¾|S|, |p| = 2|k|). B₁ → -r/12 reproduced from the Feynman-parameter integral. |
| EC-04 (0⁻ = shift-symmetric axion, g_φγ ≈ 0, ~10¹²⁰ tuning) | **SURVIVES** | w = A/m² = 1.5M_P²/H₀² ≈ 4×10¹²⁰. g ≈ 0 agrees with my result on the QCD question below. |
| EC-03 QCD mass 1.4-8.6×10⁻¹² eV, full washout, local amplitude 1.1×10⁻²⁰ deg (008 reading A) | **BROKEN** | A coupling that enters only as ∂φ·J₅ has an exact shift symmetry, so V(φ) is flat. Reading A keeps the GG̃ piece of ∂·J₅ and drops 2imq̄γ₅q, which is not allowed for massive quarks. In the toy, reading B gives V_min(φ) flat to 10⁻¹⁵ and a massless eigenstate under π⁰ mixing; reading A gives the textbook mass. The 1502.03694 and 0911.2698 "QCD axion" claims rest on the massless-fermion ABJ step. |
| EC-03 status "dead" | **SURVIVES, by a different route** | Reading B: no mass, but the photon coupling decouples as B₁ ≈ -(ω/m_f)²/12 ~ 10⁻⁷⁰. Neutrinos are light but neutral. So β ≈ 0. |
| EC-03' g = 8α/(πf) "torsion-scale coupling", β_max = 0.136° | **BROKEN as a torsion mechanism** | Under the correct reading B, torsion generates no low-energy φFF̃. "Removing the gluon coupling" while keeping the photon anomaly has no consistent derivative-coupling realization. g = 8α/(πf) is then a value put in by hand, so EC-03' is a generic ALP (an "ALP in disguise", like KK-02). The 0.136° number itself is arithmetically right. A torsion origin needs a non-derivative φFF̃ term, e.g. a Green-Schwarz-modified Bianchi identity (the KK-03 structure). |
| 008 washout vs 2512.21888, ODE ICs and units | **SURVIVES** | φ″ + 2𝓗φ′ + a²m²φ = 0, with τ in 1/H₀ and μ = m/H₀, is consistent. The frozen start at a = 10⁻⁹ is fine. The Bessel-J₁/₄, WKB-vs-full-ODE and adiabatic checks can fail. 2512.21888 Eqs. 4/6/12 hand-checked: φ₀ = 2.48×10¹² GeV, g_wash = 9.6×10⁻¹⁶. Half washout at 1.1×10⁻²⁸ eV ≈ 3H_* is plausible. Cosmetic: the `washout` column of `scan.csv` exceeds 1 at a node (1.78×10⁻²⁸ eV). |
| A_CB = 3.6×10⁻⁷(H_inf/φ_i)² | SURVIVES | (g/2)²(H/2π)²·S_L with g = 2β/(Fφ_i). Units are rad², matching `data/` (A_CB_rad2). |
| KK-01 P' all-orders argument | **SURVIVES** | Checks: P' ∈ SO⁺(4,1), det +1, time orientation kept, spinor lift exists (S² = +1), so 5D masses and spin structures are covariant. Pull-back gives A₀ → -A₀, A_i → +A_i, φ even. φᵖFF̃ is odd and φᵖF² even, while the bulk control θff̃ is even. The KK-tower probe gives 0 without a Wilson line, nonzero with w = 0.13, and 0 at w = ½. Loophole (stated in 009): any P'-odd background or condensate, including a Wilson line at a non-symmetric minimum. |
| 2608.09299 description (`literature/2608.09299.md`) | **SURVIVES** | Eqs. 3, 4, 17-18 ("with a consistent orientation"), 30, 39, 41 (2φ³ *F F term), 38 (-¼φ³F²), 14 (Jordan φR, no (∂φ)²), 52 and the PR4 fit all match. Stabilization is "beyond the scope". The paper does not say where SM matter lives. |
| KK-04 Cassini / ISL kill | **SURVIVES** | Brane matter on M₄×S¹ sees the induced metric ĝ_μν = g_μν + O(A²), which is the Jordan metric. Bulk zero modes couple to it as well. With φR and no kinetic term this is Brans-Dicke ω = 0, γ = ½, a gravitational-strength coupling to all matter wherever it sits. Evading it needs matter coupled to the φ-dressed (Einstein-frame) metric, which no brane position gives, or screening. The photon's own φ³F² also makes α vary, adding equivalence-principle bounds. |
| u-Holst c_odd = 16c³/(1+12c²) | not re-derived | Belongs to the excluded radion class. 009's E7/E8 are exact-integration checks that can fail. |
| KK-02 β = kαΔϑ/2π ≤ 0.209° | **SURVIVES, minor caveat** | Hand reduction: (k/24π²)𝒜∧𝓕∧𝓕 → (k/8π²)ϑ f∧f = (kα/4π)ϑFF̃, so β = kαΔϑ/2π. Caveat: if the potential comes from a bulk fermion of charge q > 1, the minima are spaced 2π/q, so Δϑ ≤ π/q and β ≤ kα/2q. The k ∈ Z vs 6Z question is still quoted, not derived. |
| KK-03 window (Δb < 5.5×10¹⁰ GeV, g 1.4-8×10⁻¹³ GeV⁻¹) | **SURVIVES** | Hand estimate of the stiff-fluid BBN → LSS scaling gives Δb ≈ 6×10¹⁰ GeV. |

## Verdict
No sign or number errors in 007 or 009. In 008, the reading-A quantities are wrong: the QCD mass, the full washout and the 1.1×10⁻²⁰ deg local amplitude. Reading B is correct, because a purely derivative torsion-fermion coupling gives neither a QCD potential nor a low-energy photon coupling.

**Required fixes (overseer):**
1. EC-03 stays dead, with the reason replaced by "shift symmetry: m = 0 and g_eff ~ (H/m_f)² ~ 10⁻⁷⁰". Drop the QCD-mass numbers from CLAIMS.
2. EC-03' loses its torsion-generated coupling. Reclassify it as an ALP with g set by hand. The generic-ALP region applies, and the 0.136° figure only describes one arbitrary choice of g.
3. Close the STATE open question on ∂φ·J₅ accordingly.

## Open issues
1. A torsion route to a non-derivative φFF̃ in 4D: a Green-Schwarz-modified torsion Bianchi identity, dH ∝ FF̃. If it includes GG̃, the QCD mass returns. Theorist.
2. KK-02 charge-q caveat and the k ∈ Z vs 6Z normalization.
3. KK-04: radion screening (chameleon) was not explored. It is the only escape left.
4. Process note: one file-edit helper in this task ran under system python without pin.py (no computation, under 1 s).
