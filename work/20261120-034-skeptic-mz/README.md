# 20261120-034 skeptic: audit of 031 (EC-06a, Magueijo-Złośnik), iteration 10

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism. There is no `RESULT.md`.

## Question
Does 031's kill of EC-06a hold? EC-06a is MZ torsion (1908.05184, 2001.06373) with a hand-added CFJ coupling -½ξ S_μA_νF̃^{μν}. The audit covers four items:
1. The S solution and non-closure in matter.
2. The 3.6 keV strong-coupling scale and alternative repairs.
3. The off-tracker analysis.
4. The k* regime argument and the 24 Hz input.

## Method
- Opened the ar5iv HTML of 1908.05184 and compared Eqs. 7 and 14-16 and the FRW Eqs. 26-30 with 031's transcription.
- Opened the ar5iv HTML of 2001.06373 for Eqs. 27-28, 44-50, 65, 67-70 and 72.
- Opened the PDF of 1710.05832 (PRL 119, 161101), extracting its text with a scratch script in `artifacts/034/` (not committed).
- Own script `skeptic_mz.py`, written independently of 031's code: no imports from 031, own loader, numbers from `data/`. `run.log` has **15 PASS, 0 FAIL**.
  - K1: Stückelberg gauge variation, computed from scratch.
  - K2: vertex and power counting.
  - K3: normalisation by an independent route (MZ Eqs. 26-29 → Friedmann deficit).
  - K4: separate-universe closure test.
  - K5: transport law integrated numerically.
  - K6: frequency dependence of the k* condition.
  - K7: alternative gauge-invariant repair.
- Rerun from the repo root (about 20 s, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-034-skeptic-mz/skeptic_mz.py
```

## Results

| claim (031) | status | reason |
|---|---|---|
| MZ Eq. 15 is T^[a∧e^b] = -(3/2Λ²)dΛ∧R^ab + (3/4γΛ²)ε^abcd dΛ∧R_cd; Eqs. 14 and 16 and FRW 26-30 transcribed correctly | SURVIVES | Matches the ar5iv text verbatim. Matter couples only to the tetrad, and Λ has no kinetic term in Eq. 7. |
| S_μ = (2/3γ)(κτ - Λδ)_μ^ν∂_νθ + O(∂θ·DT) | SURVIVES | Not rederived from scratch. It is consistent with K3 (the Eqs. 26-29 energy) and with the K4 tracker (K4c reproduces P = γH/3). |
| S is a gradient in vacuum | SURVIVES, with a qualification | Only at leading order. The -(2/3)(⋆DT)·∂θ term is O((∂θ)²) and is not closed in general, even in vacuum. This does not matter for the kill. |
| S is not closed in matter; (dS)_{τi} = ∂_iY ≠ 0; closure "on a measure-zero set" | SURVIVES, strengthened | K4 uses dust with Λ = 0 and the transport law ∇·(𝒢S) = (γ/4)GB. GB = C² + f(ρ) is local in ρ at linear order (K4a), so s = -u·S obeys ṡ = γκρ/3 along each worldline. A growing mode gives δs = 3s̄δ ∝ t^{-1/3}, and no choice of θ data cancels it. So closure fails for **every** solution with a growing density mode, not only generic ones. |
| Λ_SC = (16π)^{1/4}(2√3M_Pm/(ξ(1-c_s)))^{1/2} = 3.6 keV (ξ = 1.5×10⁴, PDG m) | SURVIVES | K1 gives δ(-½p·AF̃) = -¼σ(dp)·F̃ + total derivative. That yields c = ξ(1-c_s)/(2√3M_Pm) (K2a), and power counting gives M ~ c²E⁴ (K2b). Reproduced: 3.58 keV. With the other unitarity convention (\|a₀\| ≤ ½) it is 3.0 keV. In any dense environment 1-c_s → 1, giving 2.0 keV. The c_s momentum enhancement (unused by 031) only lowers it. Shortfall: 9.4 orders in energy, 18.9 in m. |
| Canonical normalisation S = -(1/√3M_P)[∂φ - (1-c_s)(∂φ)_⊥], γ-free | SURVIVES in magnitude | K3, an independent route: Eqs. 26-29 at O(γ²) give H² = (Λ+κρ)/3 + 5P² - 2γHP, which is -γ²H²/9 on the tracker (matching 022/031). So the torsion sector carries energy O(M_P²P²), meaning S ~ ∂φ/M_P with γ cancelling. The 1/√3 itself is not independently checked. An O(1) change moves Λ_SC by O(1), against a 9-order shortfall. |
| Stückelberg is the only repair | SURVIVES, with one alternative 031 states but does not quantify | (a) Higgsing the photon does not help. The canonical Goldstone h₂ couples through ξ∂φ∂h₂F̃/(√3M_Pm) whatever the radial-mode mass, so the same dim-6 operator and the same 3.6 keV appear. (b) Coupling the photon only to the exact-gradient vacuum piece S_vac = (1/γ)∂lnΛ, i.e. f(Λ)FF̃, is gauge invariant with no Stückelberg field. See the next row. |
| "Gauge-invariant remnant is an ALP coupling g = ξ/(√3M_P)" | WEAKENED (wording and numbers) | g = ξ/(√3M_P) is a perturbative statement about ∂φ around a fixed background. The exact gauge-invariant coupling is f(Λ)FF̃. For f = (ξ/γ)lnΛ, β comes only from c_s S_t (K7): ∫c_sε dln a = 0.216γ vs 4.58γ, so ξ must be **21× larger** than 031's. The vacuum coupling g = ξ/(√3M_P) then depends on which γ bound holds: 1.2×10⁻²⁰ GeV⁻¹ (N_eff only), 7×10⁻¹⁴ (GW speed), 0.24 (k* at 24 Hz), 8 GeV⁻¹ (2001.06373 Eq. 70). In matter g_eff = c_s g → 0, so stellar and helioscope bounds do not apply directly. |
| The remnant "is not torsion" | SURVIVES, with a nuance | The photon couples to the scalar Λ. Torsion enters only as the auxiliary field whose elimination gives Λ its kinetic term (N = 4M_P²/3γ²) and its GB source. By the repo rule (a mechanism names its coupling), this is ALP class with a hand-added coupling. The draft's "torsion plays no role" is slightly too strong for MZ; "torsion supplies only the scalar's dynamics" is accurate. |
| Homogeneous mode grows relative to H, ε_h(z*)/ε_h(0) = 10⁻⁴; no decaying branch | SURVIVES | K5: forward integration of the transport law reproduces the closed-form tracker (0.6667, 0.4155). Integrating backward from today with ε(0) shifted by 1 moves ε(z*) by 1.04×10⁻⁴. |
| Best case ≤ 8 MeV (tuned cancellation branch, N_eff only) | WEAKENED (number only) | The plain tracker with no GW bound at all (γ at the N_eff limit, ξ = 2.4×10⁻³) gives **8.9 MeV**, slightly above 8 MeV. The best case is 9 MeV, 6.0 orders short, not 6.1. No verdict change. |
| 2001.06373 Eq. 65 P-free γ² term (031 open issue 3) | resolved in part | ar5iv shows Eq. 65 (and Eq. 72) do contain a γ² term with no P. On the cancellation branch, though, k* is still far below k_LIGO unless γ < 10⁻¹⁷, so the speed formula's applicability is again set by k*. Irrelevant to the verdict. |
| k* = H/ε_S; LIGO 3×10¹² above k* at the 026 bound; ξ ≥ 4.7×10¹⁶ | SURVIVES, strengthened | Agrees with 2001.06373 Eq. 50, k* = (a/6P)(Λ+κρ+6P²). The paper's own Eq. 69 states \|k*\| ≫ k_LIGO is "necessary for the consistency" of the linearised equations. Its Eq. 70 gives \|γ\| ≪ 10⁻²¹, i.e. ξ ≫ 1.6×10¹⁸, stronger than 031's number. 031 cites Eqs. 69-70 but not this bound. |
| 24 Hz input (search snippet) | SURVIVES, now cited | 1710.05832 (PRL 119, 161101), p. 161101-1, main text: the inspiral lasted about 100 s "calculated starting from 24 Hz" in the detectors' sensitive band (text read from the PDF). Fig. 1's data are band-passed 30 Hz-2 kHz. The `mz_inputs.csv` note "seen via search snippet, not opened in full" can be updated by the overseer. **Dependence:** ξ_k ∝ f. For 10 Hz-1 kHz, ξ ≥ 2×10¹⁶-2×10¹⁸ (K6). The verdict does not depend on f at all: the kill uses ξ ≥ 1.5×10⁴, or even no GW bound. |

## Verdict on EC-06a
**Dead; 031's kill survives.** The kill is robust to the GW-sector question. If 2001.06373's linear tensor analysis is trusted, ξ ≥ 10¹⁶-10¹⁸ and Λ_SC ≤ meV. If it is not, ξ ≥ 2.4×10⁻³ (N_eff) and Λ_SC ≤ 9 MeV at the PDG photon mass. Either way the result is at least 6 orders in energy below 10 TeV. Non-closure in matter is stronger than 031 states (K4). The only gauge-invariant escape is f(Λ)FF̃, a hand-added ALP coupling to the MZ Λ-field, which is ALP class under the program's rules.

Requested changes (overseer):
1. In the CLAIMS EC-06a row, change "≤ 8 MeV" to "≤ 9 MeV (tracker or tuned branch, no GW bound)".
2. Replace the remnant sentence with "gauge-invariant remnant: f(Λ)FF̃, ALP class; with f ∝ lnΛ it needs ξ 21× larger, and g in vacuum ranges 10⁻²⁰-10 GeV⁻¹ depending on the γ bound".
3. Cite 2001.06373 Eq. 70 (|γ| ≪ 10⁻²¹).
4. Update the 24 Hz source note.

## Open issues
1. **f(Λ)FF̃ remnant viability.** Under the paper's own k* consistency, the vacuum coupling is g ≳ 0.2 GeV⁻¹. That would conflict with photon-ALP conversion in intergalactic fields if c_s ~ 0.7 there, but no IGM bound is in `data/`, and the field's c_s < 1 dispersion changes the conversion. This was not checked. It is ALP class in any case.
2. Sign of the θ sector energy. K3 gives a negative O(γ²) Friedmann contribution. I could not close the split between the kinetic and GB-coupling energies, so the ghost question (031 open issue 5) remains open. A ghost would kill EC-06a and the f(Λ) remnant independently.
3. K4 uses Λ = 0 dust and drops C² (second order). With Λ ≠ 0 the gradient term c_s enters, but non-closure already holds in the Λ → 0 limit relevant to galaxies and z*.
4. Not rederived: the exact algebraic solution of Eq. 15 (031 B1-B3); 031 T4's value N = 4M_P²/3γ² (only its magnitude was checked).
5. Local nonlinear physics near masses (031 open issue 4) is untouched. With c_s ~ Λ/κρ ≪ 1 in any dense matter, φ is extremely non-relativistic there, which likely makes things worse, not better.
