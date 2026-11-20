# 20261120-007 theory-nogo-minimal (4D track, iteration 2)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism.

## Question

- A. EC-02: with minimal EC (S = -(κ/2)J5) and the CFJ term -½ξS·AF̃, what is the largest β from every cited source of isotropic axial charge n5, and what bounds ξ?
- B. EC-05: the same S fed into the fermion-loop CFJ term. What are the coefficient and its ambiguity (Dobado-Maroto, Jackiw-Kostelecky), and how does Reig 2611.07578 apply?
- C. EC-04: which PGT sector gives a propagating parity-odd torsion mode on FRW? What are its FRW equation, its g_φγ, and its bounds?

## Method

Symbolic checks in sympy and numerical cosmology in scipy. Write-up: `derivation_nogo.md`. Measured inputs come from the new `data/cosmo_inputs.csv` (cited rows), `data/birefringence_measurements.csv` and `work/20261120-005-bounds/derived_bounds.csv`. Nothing measured is hardcoded. Rerun (single core, about 5 min total, each script under 1.5 min):

```
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-007-theory-nogo-minimal/run_checks.py
```

| script | content |
|---|---|
| check_helicity.py | J5^0 per particle = h·v (Dirac spinors), axial current C-even, helicity conserved in FRW and by the axial-torsion coupling |
| ec02_nogo.py | n5 for each source, ∫n5 dt from LSS, max β at ξ = 1, ξ needed, Reig cross-check, helioscope bound on ξ |
| ec05_loop.py | DM and JK coefficients converted to our normalisation, B1 decoupling series, max β |
| check_pgt_0minus.py | pseudoscalar curvature = -6∇·S for axial torsion; auxiliary-field reduction to one scalar; ghost/tachyon conditions; propagator pole; ∂K = ½FF̃ |
| ec04_estimates.py | 0⁻ axion parameters, thawing and kinetic FRW solutions, bounds |

Last run: `run.log`, 29 PASS, 0 FAIL.

## Result

**A. EC-02** (β = (ξκ/4)∫n5 dt; n5 = helicity × speed, so non-relativistic relic neutrinos are velocity-suppressed today, but the integral is dominated near z ~ 1000 where they are still relativistic):

| source (cited input) | max β at ξ = 1 | shortfall | ξ needed |
|---|---|---|---|
| ν asymmetry, \|ξ_ν\| ≤ 0.048 (95%, 2409.08280 T4, flavours equilibrated) | 7.5×10⁻³⁰ deg | 28.5 orders | 2.9×10²⁸ |
| ν, ξ_μ,τ ≤ 0.57 from N_eff alone (1807.06209), generous | 6.3×10⁻²⁹ deg | **27.5 orders** | 3.4×10²⁷ |
| fully helicity-polarised fermion DM, Lyman-α speed bound (1702.01764) | 1.6×10⁻³⁶ deg | 35 orders | 10³⁵ |
| electrons/baryons | 0 (isotropic); envelope 39 orders short | | |
| B-field spin polarisation (1502.01594) | 0 isotropic; spatial S ≲ 3×10⁻⁷⁷ eV | | |

What bounds ξ:
- The ξ-term breaks U(1) off the FRW background, because S ∝ J5 is not closed (same variation as Reig Eq. 4).
- With S = -(κ/2)J5 it is the dimension-6 operator J5·AF̃/Λ² with Λ = M_P√(2/ξ). The needed ξ ~ 10²⁸ means Λ ≈ 2-6×10⁴ GeV.
- Reig's Kalb-Ramond/dual-axion completion with a universal axial current gives g_aγ g_ae = ξκm_e. The helioscope bound then requires **ξ < 9×10¹⁷**, which excludes the needed value by ≳ 9.6 orders.

Our time integral reproduces Reig's Λ* to 20% (12 vs 9.8 TeV).

**B. EC-05**: |p| = r|S| with:
- r = (α/π)ΣQ² (Dobado-Maroto, dim. reg.)
- r = (9α/8π)ΣQ² (JK nonperturbative, = 9/8 DM)
- r = 0 (Coleman-Glashow, Pauli-Villars)

With ΣQ² = 8, r ≤ 0.021, so max β = 1.3×10⁻³⁰ deg (generous ν source), **29.2 orders short** (30.1 for the cited ξ_ν bound). Reig does not bind here: EC-05's scale Λ = M_P√(2/r) = 2.4×10¹⁹ GeV is far above Reig's limit Λ > 4×10⁹ GeV. The no-go comes from the size. On FRW, S = ∂θ (homogeneous ⇒ closed), so the ambiguity is resolved by decoupling (B1 ≈ -q²/12m_f²), with a further (H/m_e)² ~ 10⁻⁷⁰ suppression.

**C. EC-04**:
- The only healthy light parity-odd torsion mode is the PGT 0⁻ mode: axial torsion S_0(t) with kinetic term from the pseudoscalar curvature squared (Chen et al. 0908.3323 Eqs. 7, 11, 64; Baekler-Hehl-Nester 1009.5112).
- For axial torsion E ∝ ∇·S (checked: ε^{μνρσ}R_{μνρσ} = -6∇·S). L = (A/2)S² + (w/2)(∇·S)² + ¾S·J5 reduces to one pseudoscalar φ with m² = A/w, healthy iff A > 0 and w > 0, and S = ∂φ/√A.
- It is an axion with universal derivative coupling to fermion axial currents, f = (2/3)√A ≈ 2×10¹⁸ GeV (EC value), g_ae = 2.6×10⁻²². Its shift symmetry is exact, so it gets no QCD mass. It has **no photon coupling** (the loop decouples, g ≈ 10⁻⁹⁹ GeV⁻¹).
- Adding the EC-01 term gives g_φγ = -ξ/√A = 3.4×10⁻¹⁹ ξ GeV⁻¹.
- FRW: φ̈ + 3Hφ̇ + m²φ ≈ 0. β = 0.215° needs a thawing field with m ~ H₀ and θ_i ≈ 0.06/ξ, so S_0 ≈ 2×10⁻³⁵/ξ eV today and Ω_φ ~ 10⁻³/ξ², or an equivalent kinetic mode.
- Bounds: lab (ξ ≳ 10⁻¹⁴), CAST (ξ < 2×10⁸), red giants (no constraint at f ~ M_P).
- Requires w ≳ 4×10¹²⁰ (or ≳ 10¹⁰³ if A is lowered to the red-giant limit). With O(1) coefficients, m ~ M_P.

## Verdict

- **EC-02: dead as a torsion mechanism.** Minimal EC only supplies the scale M = M_P for the neutrino-axial-current operator J5·AF̃/M². The signal needs M ≈ 2-6×10⁴ GeV, 27.5-28.5 orders in ξ beyond ξ = 1 with cited ν asymmetries. Every other source is 35-39 orders short or gives no isotropic n5.
  - This is not a no-go for neutrino-asymmetry birefringence in general. The vector/left-handed-current operator with a new TeV-scale M is a published candidate (Geng-Ho-Ng 0711.4617; Zhou-Huang-Geng 2302.11140). Torsion plays no role in it, and Reig 2611.07578 argues that it too is excluded.
- **EC-05: dead** (preliminary, pending skeptic). It is 29-30 orders short even with the largest regularisation (JK, ΣQ² = 8 with no decoupling), and the physical (gradient-background) value is smaller by another ~70 orders.
- **EC-04: not excluded by bounds, but it is not a torsion mechanism for β.** The PGT 0⁻ mode is a Planck-coupled ALP with no photon coupling of its own. β needs the EC-01 coupling put in by hand (ξ ~ 1) and a 0⁻ mass ~H₀, which takes a ~10¹²⁰ tuning of the E² coefficient. It merges into the EC-03 / axion class. Recommend recording it as "ALP in disguise, unnatural", like KK-02.

## Open issues

1. EC-04 signs: the ghost/tachyon conditions are derived here in our normalisation (A > 0, w > 0). Mapping them onto the PGT coefficients (b₋, m₋ = A₀ + 2A₃) of 0908.3323 and Sezgin-van Nieuwenhuizen needs their sign conventions. The text extraction of Eq. 64 lost the sign.
2. EC-04 full theory: does a huge b₋E² affect other sectors (the 0⁺ mode via BHN's RX term, gravitational waves)? Is the 0⁻ mode's E² term radiatively stable? Full cosmology is a later task.
3. Exact shift symmetry ⇒ no QCD potential for any purely derivatively coupled torsion pseudoscalar. This contradicts the "QCD mass" threat listed for EC-03 (1502.03694). The EC-03 task should check whether the Nieh-Yan/Holst field there couples only through ∂φ·J5.
4. Helioscope bound (8.1×10⁻²³ GeV⁻¹) and red-giant bound are taken as quoted in Reig; the original papers were not opened.
5. The thermal-relic Lyman-α bound is used as a proxy for the DM speed bound; a non-thermal velocity distribution could shift the DM numbers by O(1). They stay 35 orders short.
6. `work/20261120-003-theory-ec-photon/inputs.json` placeholders (ξ_ν = 0.1, m_DM = 1 keV) are superseded by `data/cosmo_inputs.csv`. Its "28 orders" is confirmed (28.5).
7. CONVENTIONS.md: no errors found in the parts used.
