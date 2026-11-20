# Open gaps in the torsion-birefringence literature

Task 20261120-002. These are the scout's judgment, drawn from `prior_work_torsion_birefringence.md`. Each gap is written as a task candidate for the overseer. Order-of-magnitude numbers marked "scout estimate" are unchecked and must be redone in `CONVENTIONS.md` units before use.

## Highest value (likely to settle a route fast)

**G1. EC fermion-current route as a single chain: torsion → CFJ → β, with the Reig no-go applied.** (4D track, likely no-go)
No paper combines three published pieces:
- Minimal EC: S_μ ∝ κ J⁵_μ.
- Dobado-Maroto: a fermion loop gives a CFJ term ∝ (α/6π) ΣQ² S_μ (regularization dependent per Jackiw-Kostelecky).
- Reig 2611.07578: the resulting j·A F̃ operator needs 1/Λ*² ≈ (10 TeV)⁻² and is excluded anyway.

Task:
1. Derive the coefficient of J⁵_μ A_ν F̃^{μν} in minimal EC. Include the scheme dependence and the self-consistency (the same fermions source S and run in the loop).
2. Compute β for the largest allowed relic axial asymmetry (neutrinos, ADM).
3. Compare with the KRT lab bound and with Reig's chain.

Scout estimate: about 35 orders of magnitude short on S, or about 28 orders on 1/Λ*². That would be a clean Result-B entry with a cited mechanism of death.

**G2. Torsion-descended axion (KR or BI/Nieh-Yan) with its forced QCD coupling, confronted with β.** (4D track, likely no-go unless extended)
Castillo-Felisola et al. 1502.03694 show the torsion axion has f ≈ 4×10¹⁸ GeV, a photon coupling via the anomaly, and a gluon coupling that is generic because torsion couples universally to quarks. That gives a QCD mass of ~10⁻¹² eV. Lattanzi-Mercuri identify it with the QCD axion. No paper asks whether it can make β. Task:
1. Compute β = (g/2)Δφ for this axion with a realistic misalignment amplitude. Apply the washout argument of 2512.21888 (an ALP with m ≫ 10⁻²³ eV oscillates before recombination).
2. State what extension would remove or cancel the QCD potential (non-universal torsion couplings à la Lucat-Prokopec or Chakraborty-Lahiri, a second PQ field), and what each costs.

**G3. KK radion-driven Chern-Simons birefringence (Rao et al. 2608.09299) redone in Einstein-Cartan rather than teleparallel, with KK constraints applied.** (KK track)
Rao et al. get Δα = 3c ln(ϕ₀/ϕ_LSS), with c·Δlnϕ ≈ 1.75×10⁻³ from Planck PR4. In the sections read, none of the following is done:
- (a) The Einstein-Cartan (curvature + torsion) version of the 5D parity-odd term.
- (b) Confronting the same radion excursion with fine-structure-constant variation. In KK, α_em depends on ϕ, and CMB, quasar and atomic-clock limits on Δα_em/α are ~10⁻³ at z~1100 and much tighter at low z (bounds analyst to cite). This forces c ≳ O(1) and may conflict with the GW Nieh-Yan bound (Wu et al. 2110.13870).
- (c) Radion stabilization and fifth-force limits.
- (d) An ACT DR6 fit and the anisotropy (SPT A_CB).

This is the most direct "prior art plus missing constraints" target for the KK track.

## Medium value

**G4. KK-reduced massless Kalb-Ramond (torsion) axion against CMB β and early-universe energy density.**
Maity-SenGupta-Sur hep-th/0507210 (S¹ scheme) give a closed-form β ≈ h t_L and need ρ_KR ~ ρ_DM today for 0.2°, using radio galaxies at z = 0.4. For a massless dual axion φ̇ ∝ a⁻³, so ρ_KR ∝ a⁻⁶ (stiff), and the rotation integral is dominated by early times. No paper checks the CMB-β value against N_eff/BBN and CMB limits on a stiff component, or computes the resulting EB spectral shape. Scout expectation: strongly constrained. It needs a number.

**G5. Itin-Hehl F²T² induced axion θ = c₁ u v in FRW.**
This is a gauge-invariant, parity-even coupling that gives an isotropic rotation ∝ d(uv)/dt, using exactly the two FRW-allowed torsion components. No paper evaluates it: what sources u (trace) and v (axial) in EC and its extensions, how big ℓ² can be (Preuss white-dwarf bound k² ≲ (19 m)²), and what β results. In minimal EC with Dirac matter u = 0, so θ = 0. Task: show that, then find the minimal extension giving u ≠ 0 (Magueijo-Złośnik parity-violating FRW torsion is one candidate) and compute β.

**G6. Das-Mohanty-Prasanna-type dimension-5 coupling: chromaticity test.**
The only published torsion-coupling fit to CMB rotation (0908.0629) predicts β ∝ ν². ACT DR6 and Planck report β per frequency. Task: take the per-frequency β consistency from `literature/2509.13654.md` and the Planck notes and bound or kill the ν² scaling. Also fix the inconsistent T₁ definition and check which FRW components it corresponds to.

**G7. No joint confrontation of any torsion mechanism with ACT DR6 β + SPT A_CB + lab torsion bounds.**
The SAME-AIM papers use pre-Planck data (Das et al.), radio galaxies (Maity et al.) or Planck PR4 alone (Rao et al.). None uses the ACT DR6 EB/TB β, the SPT-3G A_CB bound, or the KRT torsion bounds jointly. Any surviving candidate from G1-G6 should be fit this way. This is the program's Result-A step 4.

## Lower value / bookkeeping

**G8. Scheme dependence of the torsion → CFJ coefficient in the cosmological setting.**
Jackiw-Kostelecky show k_AF ∝ b is undetermined perturbatively. No paper fixes it for torsion in a finite-temperature, expanding background with a real asymmetric fermion population (a medium effect, as opposed to vacuum b_μ). If G1 is borderline (it is not expected to be), this matters.

**G9. EC-KK reduction with a 5D Chern-Simons A∧F∧F or 5D Nieh-Yan term.**
Kohler gr-qc/9808004 and Shankar-Wali 0904.4038 do EC-KK without parity-odd terms. Rao et al. do parity-odd KK without EC. Nobody does both in one model (the overlap of these is the KK track's core task). The A₅ Wilson-line pseudoscalar (extranatural, hep-th/0301218) coupled through a 5D CS term is a standard non-torsion baseline the KK track should compare with.

**G10. Homogeneous torsion without spin sources.**
Boehmer-Bronowski show Weyssenhoff fluids are incompatible with FRW. Magueijo-Złośnik show parity-violating FRW torsion can exist without spinors when Λ varies. Nobody connects such S₀ backgrounds to photons. Lower priority, because the photon still needs one of the couplings above.
