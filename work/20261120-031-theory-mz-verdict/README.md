# 20261120-031 theory: EC-06a (Magueijo-Złośnik) verdict, iteration 9

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism. There is no `RESULT.md`.

## Question
029 BLOCKER 1. Does EC-06a survive? EC-06a is MZ: an EC extension with Λ as a field plus quasi-Euler and Pontryagin terms (1908.05184, 2001.06373), and a CFJ coupling -½ξ S_μA_νF̃^{μν}. Four sub-questions:
1. Is S closed off FRW?
2. Where is the strong-coupling scale (025/028 machinery)?
3. Is there an off-tracker branch that evades GW170817?
4. Verdict.

## Method
- MZ covariant field equations (Eqs. 14-16) and the FRW equations (Eqs. 26-30) taken from the ar5iv HTML of 1908.05184.
- Tensor-mode equations (Eqs. 44-48, (III)) and the Eq. 65/72 structure taken from the ar5iv HTML of 2001.06373.
- GW170817 band (24 Hz) and distance (40 Mpc) from 1710.05832. The 24 Hz figure was read in a search snippet of the PDF; the abstract gives 40 Mpc.
- sympy:
  - Eq. 15 on FRW reproduces MZ Eqs. 28-29.
  - Exact algebraic solution of Eq. 15 for the axial vector, on random curvature.
  - Linear-perturbation closure.
  - Leading-order transport law vs MZ Eq. 30.
  - θ effective action, Stückelberg vertex, tensor singular scale.
- Numbers loaded from `data/`. New file: `data/mz_inputs.csv`.
- Write-up: `derivation_mz.md`.
- Rerun from the repo root (about 30 s, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-031-theory-mz-verdict/run_checks.py
```
`run.log`: **36 PASS, 0 FAIL** (check_offfrw 12, check_dynamics 14, mz_numbers 10).

## Result
1. **Not closed.** Eq. 15 solves exactly to S_μ = (2/3γ)(κτ_μ^ν - Λδ_μ^ν)∂_νθ + O(∂θ·DT), with θ = 3/(2Λ).
   - In vacuum it is a gradient.
   - In matter it is not. At linear order, (dS)_{τi} = ∂_iY with Y = (2κ/3γ)[(ρ+p)δθ_c′ + p′δθ_c + θ′δρ_c], which is gauge invariant and vanishes only for pure time shifts. θ carries its own initial data.
   - In canonical form, S = -(1/√3M_P)[∂φ - (1-c_s)(∂φ)_⊥] with c_s = Λ/(Λ+κρ). So U(1) breaks wherever ρ ≠ 0.
   - The breaking factor 1 - c_s is Ω_m today and about 1 at z* and in galaxies. A Stückelberg photon is required.
2. **Strong coupling.** The θ (Λ-field) exchange with the Stückelberg vertex c = ξ(1-c_s)/(2√3M_Pm) gives Λ_SC = (16π)^{1/4}(2√3M_Pm/(ξ(1-c_s)))^{1/2}.
   - At ξ = 1.5×10⁴: **3.6 keV at the PDG limit m = 10⁻¹⁸ eV**, and 0.4 eV at 10⁻²⁶ eV.
   - A 10 TeV-safe theory needs m ≈ 8 eV, 10^18.9 times the PDG limit.
   - 029's ~0.1 GeV is right as dim-7 arithmetic (0.16 GeV), but it is not the lowest scale.
3. **No loophole.**
   - The general solution is P = tracker + C/(a³H²). The homogeneous mode grows relative to H (ε_h(z*)/ε_h(0) = 10⁻⁴), so no branch has ε large at z* and small today, except a tuned cancellation of C near z = 0.
   - That cancellation evades GW170817's speed only if 2001.06373 Eq. 65's P-free γ² term does not apply. With N_eff alone it would need ξ ≈ 3×10⁻³, but Λ_SC is still ≤ 8 MeV.
   - New gravity-sector problem: the tensor equation degenerates at k* = H/ε_S, γ-independent. 2001.06373 requires k* ≫ k_LIGO. At the 026 bound, LIGO sits 3×10¹² above k*, so the c_T = 1 + γ²/18 bound (ξ ≥ 1.4×10⁴) was used outside its regime. Respecting k* needs ξ ≥ 4.7×10¹⁶ (tracker) or ξ ≥ 1.5×10¹⁴ (cancellation branch).

## Verdict
**EC-06a: dead.** Named reason: photon-sector strong coupling. Breaking U(1) in matter forces a Stückelberg photon, and Λ-field exchange then becomes strongly coupled at Λ_SC ≤ 3.6 keV (tracker, ξ ≥ 1.5×10⁴, m ≤ PDG). That is **9.4 orders in energy** below the 10 TeV where γγ physics is standard, or **18.9 orders in m_γ**.

The best case is the tuned off-tracker branch with no GW bound: 8 MeV, still 6.1 orders short in energy and 12.2 in mass.

Independently, the gravity sector's own consistency condition (k* ≫ k_LIGO) raises the required ξ to ≥ 10¹⁴-10¹⁷.

Not EC proper in any case: MZ adds quasi-topological terms with 1/Λ and an Λ field. The only gauge-invariant remnant is an ALP coupling of the Λ-field (g = ξ/√3M_P), which is ALP class with a hand-added coupling and not torsion.

Consequence for STATE: 4D is Result B for torsion, with the ALP statement. This assumes the overseer and skeptic accept this.

## Open issues
1. The tensor-torsion part q of Eq. 15 is not tracked in the linearised Eq. 16. It could add O(1) terms to the transport law, but it cannot restore closure, because S is fixed algebraically by Eq. 15.
2. The ghost reading of the k* sign flip is inferred from the EOM coefficient A/(A²+B²), not from a reduced action. The paper calls it a singularity.
3. Whether 2001.06373 Eq. 65's γ² term is P-independent at P = 0 was not rederived. It matters only for the cancellation branch's GW-speed bound, not for the verdict.
4. Local physics near masses (Eq. 16 with GB ≫ Λ²): either Λ_local ~ √GB, an O(1) change of local gravity, or large local torsion. Not solved. This could be an independent kill.
5. Power counting only, with O(1) factors untracked. The c_s < 1 momentum enhancement, which lowers Λ_SC further, is not used. The θ kinetic sign follows MZ's action sign.
6. The 025/028 Stückelberg window edges were not recomputed for MZ. They are not needed for the verdict.
7. `data/mz_inputs.csv`: the 24 Hz figure is from a search snippet of 1710.05832, with the full text not opened. The 10 TeV threshold is a choice, as in 028.
