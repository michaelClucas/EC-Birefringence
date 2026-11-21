# 20261121-043 theory: Result C costs (KK-02, KK-03, EC-06a remnant, EC-04)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Every route below is a candidate mechanism with a hand-added photon coupling (ALP class). No `RESULT.md`.

## Question
Skeptic 040 MAJOR 2: quantify the costs left unquantified in the draft's cost table.
1. KK-02: the initial-angle tuning of the Wilson-line axion (β = kαΔϑ/2π), for ACT DR6 (central and 1σ) and for the joint ACT + Planck value β = 0.277° ± 0.057° (Eskilt, arXiv:2608.06480), and the k needed.
2. KK-03: λ/κ₄ for the g window and the pre-BBN kick.
3. EC-06a remnant (MZ Λ-field with f(Λ)FF̃): the allowed g range, consistency with Chandra and the β fit, and the ghost question.
4. EC-04: confirm the ~4×10¹²⁰ tuning.

## Method
- New data: `data/alp_costs_inputs.csv`. It holds the 2608.06480 value (abstract opened in this task), 2001.06373 Eq. 70, the 0.057 GW-speed coefficient from 026/034, and two stated choices for the conversion estimate. Every other input is loaded from existing `data/` files.
- Scripts. Each prints PASS/FAIL, uses one thread and runs in under 3 min.
  - `kk02_hilltop.py`:
    - closed-form full-relaxation tuning;
    - the field equation ϑ̈ + 3Hϑ̇ + m² sin ϑ = 0, integrated for 460 initial angles (fine grid 790) and 40 masses (fine grid 60) between H₀ and 200 H₀;
    - a grid-convergence check.
  - `kk03_kick.py`: an independent recomputation of 009's window, plus the kick temperature.
  - `ec06a_remnant.py`:
    - recomputes 034's tracker integrals and g limits;
    - naive Chandra transfer;
    - an estimate of photon-φ conversion suppression from the field's c_s dispersion.
  - `ec04_tuning.py`.
- Rerun from the repo root (about 2.5 min):
  `.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261121-043-theory-alp-costs/run_checks.py`
- `run.log`: **37 PASS, 0 FAIL**.

## Result

| route | cost | value | source |
|---|---|---|---|
| KK-02, k = 1 | ceiling under full relaxation (start from rest, |Δϑ| ≤ |ϑ_i| ≤ π) | 0.209°: 0.08σ below ACT, 1.19σ below 2608.06480. ACT central needs ϑ_req = 1.028π, so it is unreachable in this regime | `kk02_hilltop.py` H1, H3 |
| KK-02, k = 1 | ACT 1σ lower edge (0.141°) | start within 0.326π of the hilltop: P = 0.33 (θ_i uniform, full relaxation). Joint P = 0.28 (θ_i uniform, ln m uniform on [H₀, H_*]) | H2, H6, joint |
| KK-02, k = 1 | β ≥ ACT central 0.215° | only by **overshoot**: the field starts near the hilltop, m ≈ 6-200 H₀, and is caught today past the minimum (|Δϑ| up to 1.59π at m ≈ 27 H₀). Joint P ≈ **0.7%** (6.7-7.3×10⁻³). At the best mass, 12% of initial angles | H7, H8, H10 |
| KK-02, k = 1 | 2608.06480 central 0.277° (ϑ_req = 1.325π) | overshoot only, m ≈ 9-44 H₀. Joint P ≈ **1.3-1.4×10⁻⁴**. 1σ lower edge 0.220°: P ≈ 0.55% | H9 |
| KK-02, k = 2 | the generic level needed for 2608.06480 | ACT central ϑ_req = 0.514π (joint P = 0.42); 2608.06480 central 0.663π (joint P = 0.29). With k ∈ 6Z (no mixed CS term), k = 6: 0.171π / 0.221π | H4, joint |
| KK-03 | Green-Schwarz coefficient | λ/κ₄ = (1.47-8.67)×10⁶, i.e. M_s = (0.83-2.0)×10¹⁵ GeV (ACT central). Range over the ACT 1σ interval: (0.96-1.98)×10⁶ to 8.67×10⁶. For 2608.06480: (1.89-8.67)×10⁶. Natural λ = κ₄ gives ΔN_ν = 3.5×10¹¹ | C1-C3, C8 |
| KK-03 | pre-BBN kick | b is exactly massless and needs a velocity. Its stiff energy at 1 MeV must be ρ_b/ρ_rad = 7.5×10⁻⁴ (g = Chandra edge) to 2.6×10⁻² (BBN edge), i.e. ΔN_ν = 0.0046-0.16, a factor 35 in energy. The velocity must originate at T ≲ 6 MeV (g_min) to 37 MeV (g_max), or the Universe was kination-dominated above that temperature. No source mechanism is known. Timing needs no tuning: 68% of Δb accrues between z* and z = 500 | C4-C7 |
| EC-06a remnant | vacuum coupling g = ξ/(√3M_P) (lower limits; every γ bound is an upper bound) | ACT: ≥ 1.2×10⁻²⁰ (N_eff only), ≥ 7.4×10⁻¹⁴ (GW-speed formula outside its regime), ≥ 0.24 (k* consistency at 24 Hz), ≥ 8.2 GeV⁻¹ (2001.06373 Eq. 70). ×1.29 for 2608.06480 | M1, M2, M4 |
| EC-06a remnant | Chandra and β fit | β is fitted on every branch (ξγ = 3.46×10⁻² for ACT, 4.46×10⁻² for 2608.06480). Naive Chandra transfer (g < 8×10⁻¹³) allows only the N_eff and GW-speed branches, which are exactly the branches outside the regime where MZ's own GW analysis holds. The self-consistent branches exceed the Chandra bound by 2.9×10¹¹ (k*) and 10¹³ (Eq. 70). With φ's own dispersion ω = c_s k (c_s ≤ 0.69 today), the photon-φ mismatch is qL ≈ 7×10²⁸ (1 keV, 1 kpc), so conversion is suppressed by ~10⁻⁵⁷ and Chandra does not bound g at all (estimate) | M3, M5, M6 |
| EC-04 | 0⁻ mass tuning | w = A/m² = 1.5(M_P/H₀)² = **4.3×10¹²⁰** at m = H₀ (1.7×10¹¹⁹ to 4.8×10¹²¹ for m = 0.3-5 H₀). The coefficient of (∇·S)² must be 10¹²⁰ instead of O(1); with O(1) coefficients the mode sits at m ≈ 3×10¹⁸ GeV | T1, T2 |

Cost strings for the RESULT table:
- KK-02: "k = 1: β ≤ 0.209° unless the field overshoots (start near the hilltop, m ≈ 6-200 H₀); P(β ≥ 0.215°) ≈ 0.7%, P(β ≥ 0.277°) ≈ 1.4×10⁻⁴, P(β ≥ ACT −1σ) = 0.28 (θ_i uniform, ln m uniform on [H₀, H_*]); k = 2 generic (P = 0.42 / 0.29)."
- KK-03: "λ/κ₄ = (1.5-8.7)×10⁶ (M_s ≈ 0.8-2×10¹⁵ GeV); massless b needs a velocity created at T ≲ 6-37 MeV carrying ΔN_ν = 0.005-0.16 at BBN, no known source."
- EC-06a remnant: "g ≥ 1.2×10⁻²⁰ to ≥ 8 GeV⁻¹ depending on the MZ γ bound; the GW-consistent branches (g ≥ 0.24 GeV⁻¹) fail naive Chandra by 10¹¹-10¹³ and survive only through the c_s dispersion (estimate); θ-sector ghost question open."
- EC-04: "(∇·S)² coefficient tuned to 4×10¹²⁰ so that m ~ H₀."

## Verdict
- **KK-02.** 009's "k = 1 needs a hilltop start" understated the problem. Under full relaxation k = 1 cannot reach 0.215° at all. It does so only through the overshoot channel, with about 0.7% joint probability, and reaches 0.277° with about 10⁻⁴. k = 2 reaches both values generically. Still a live candidate mechanism, ALP in disguise. Its cost at k = 1 is a percent-level (ACT) or 10⁻⁴-level (2608.06480) coincidence of initial angle and mass, or else a CS level k ≥ 2.
- **KK-03.** 009's numbers are confirmed: λ ≈ 10⁶-10⁷ κ₄, and a kinetic-energy injection at T ≲ 6-37 MeV with no known source. Live and narrow.
- **EC-06a remnant.** The g range is 10⁻²⁰-10 GeV⁻¹ (lower limits). Only branches outside MZ's own regime of validity pass naive Chandra. The self-consistent branches survive only through the in-medium dispersion estimate. Viability is not established, and a ghost would kill it.
- **EC-04.** 4.3×10¹²⁰ is confirmed. It means the dimensionless coefficient of (∇·S)² in the PGT action must be ~10¹²⁰, instead of O(1), so that the torsion mode is as light as H₀.

## EC-06a ghost question
MZ's Λ field has no kinetic term in the action (1908.05184 Eq. 7). Its kinetic term, -(N/2)Z^{μν}∂θ∂θ with N = 4M_P²/(3γ²), appears only after torsion is eliminated, and its sign follows the overall sign of MZ's action, which 031 took as given (031 T4, open issue 5). Skeptic 034's K3 found that the θ sector adds -γ²H²/9 to H² on the tracker: a negative energy density. It could not split this into a kinetic part (sign of N) and a θ·GB coupling part, so a ghost is neither shown nor excluded. Two calculations would settle it, and both are well defined:
- Torsion eliminated, θ·GB with a kinetic term is a shift-symmetric Horndeski theory (a linear scalar-GB coupling corresponds to G₅ ∝ ln X). Compute the quadratic action for scalar and tensor perturbations on the tracker FRW background and evaluate the Kobayashi-Yamaguchi-Yokoyama no-ghost conditions 𝒢_T > 0 and 𝒢_S > 0 (arXiv:1105.5723).
- Equivalently, do a Hamiltonian count of MZ with torsion kept, to check that the θ momentum has a positive-definite kinetic matrix.

If 𝒢_S < 0 on the tracker, EC-06a and its f(Λ)FF̃ remnant both die independently of every bound above.

## Open issues
1. KK-02 probabilities use a prior: θ_i uniform, and ln m uniform on [H₀, H_*] (M l uniform in the Hosotani mass). A different mass prior changes the 0.7% and 1.4×10⁻⁴ figures but not the full-relaxation numbers. Masses above 200 H₀ are treated as fully relaxed; the overshoot tail there is ≤ 3×10⁻⁶. The grid-convergence check gives about 8% (ACT) and 12% (2608.06480) differences.
2. In the KK-02 overshoot channel the field is still oscillating today with a period of a few Gyr. β is effectively static, but the late-time oscillation and its energy density (negligible for f ≪ M_P) were not checked against other observables.
3. The KK-03 kick temperature holds g* = 10.75 fixed above 1 MeV; muons change this by O(10%) near 37 MeV. A kination era above that temperature was not checked against other early-universe bounds.
4. The EC-06a conversion estimate assumes φ keeps the dispersion ω = c_s k at keV momenta, i.e. that the classical torsion-eliminated EFT holds there. This is unverified. Lab bounds (light shining through walls, vacuum birefringence) would be suppressed the same way and were not computed. No IGM conversion bound is in `data/`.
5. 2608.06480's dust-mitigated value has no central value in its abstract. Only the face-value 0.277° ± 0.057° is used.
