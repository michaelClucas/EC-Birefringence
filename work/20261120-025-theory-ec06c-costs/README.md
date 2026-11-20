# 20261120-025 theory: EC-06c repair costs, f(T) frame and GW speed, EC-06d at recombination (iteration 7)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism. There is no `RESULT.md`.

## Question
1. EC-06c (f(T) trace torsion, CFJ -½ξT_μA_νF̃^{μν}, ξ = 3.6×10⁻⁴) breaks U(1) off FRW. What does each repair cost against bounds: (a) a Stückelberg/Proca photon, (b) a gravitational Chern-Simons counterterm (0907.5159), (c) a closed replacement of T_μ? Skeptic 026 asked that the f(T) Lorentz scalar ζ be included.
2. Is c_T = 1 in the f(T) background? Do f(T) cosmology bounds bite? What selects the frame?
3. EC-06d (SNY 0⁺ PGT): T_t at z*, the ξ needed, and closure off FRW.

## Method
sympy (second-order f(T) action with a local-boost scalar, rotated tetrads, tensor modes; Stückelberg identity; Proca+CFJ dispersion), a numeric two-mode conversion check, a backward integration of the SNY background equations, and an order-of-magnitude bound script. Derivation: `derivation_ec06c_costs.md`. Inputs: `data/ec06c_costs_inputs.csv` (new, cited) plus existing `data/` files. Rerun from the repo root (about 4 min, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-025-theory-ec06c-costs/run_checks.py
```
Last run: `run.log`, **70 PASS, 0 FAIL** (check_stueckelberg.py takes about 3 min, the others under 25 s).

| script | content |
|---|---|
| check_ft_lorentz.py | f(T) boost scalar ζ solved from the second-order action, dT on solutions, rotations fixed by the antisymmetric equation, tensor c_T and friction, weak-field T_sc |
| check_stueckelberg.py | Stückelberg identity, Proca+timelike-CFJ dispersion, mixing element and conversion amplitude |
| ec06c_numbers.py | size of dp (local and cosmological), photon-mass window, EFT cutoff, gravitational CS sizes, ln T endpoint dependence, torsion-dragged scalar bounds |
| sny_background.py | SNY Eqs. 43-47 integrated back to z* for Table I cases; writes `sny_background.csv` |
| check_sny_limit.py | SNY covariant closure (Eq. 20); GR-like regime: ε, β/ξ and the Friedmann shift |

## Result
**f(T) frame and dT (answers skeptic 026).** With the extra Lorentz scalar ζ, (dT)_{τi} = ∂_i[(Φ+Ψ)′ + ∇²ζ - ζ″]. The f(T) antisymmetric equation (∝ f_TT, empty in TEGR) fixes ζ uniquely and non-dynamically: ∇²ζ = -3Ψ′ + 3(a″/a′)Ψ - 3(a′/a)(Φ+2Ψ). ζ is not zero, and **dT ≠ 0 on solutions**, even in matter domination with static potentials ((dT)_{τx} = 15kφ₀sin kx/τ, where 022's ζ = 0 formula gave 0). Sizes: |dT| = |C|H|∇Φ| (C = -7.5 in matter domination, -4.4 today) for sources at rest in the cosmic frame, and ≈ v|∇∇Φ| for moving sources (dominant locally). So U(1) is broken on linear perturbations everywhere. Local rotations of the frame are also fixed (∇ω = 0) by f_TT ≠ 0, but the CFJ term sources the same equation.

**What selects the frame:** the f(T) antisymmetric equations, for any f_TT ≠ 0 (f_TT drops out of the solution). EC-06c is effectively TEGR/ΛCDM-like dynamics plus a photon coupling in the frame that f(T) selects. β needs only H(z), so f(T) background bounds (1308.6142: viable f(T) ≈ ΛCDM) do not constrain it. In the f_TT → 0 limit the frame becomes strongly coupled to the photon term itself. Skeptic 026's proper-frame TEGR agrees with EC-06c on FRW only. It leaves ζ unspecified (ζ = 0 gives 022's T5).

**GW speed:** c_T = 1 exactly for any f(T) on flat FRW. Only the friction changes, by f_TT T̄′/f_T (agrees with 1801.05827). GW170817 does not bind.

| repair | what it needs | bounds | status |
|---|---|---|---|
| (a) Stückelberg/Proca photon | photon mass; π couples as -(1/4m)π(dp)·F̃, conversion amplitude ≈ ∫\|dp\|dl/(4m) | lower: m ≳ 7×10⁻²⁸ eV (solar-surface photons, moving-Sun term, A < 0.1, vacuum mixing); upper: 10⁻²⁶ eV (Proca galactic field, hep-ph/0306245) or 10⁻¹⁸ eV (PDG 2024) | **consistent in a narrow window**, 7×10⁻²⁸-10⁻²⁶ eV (wider with PDG only). The mass also removes the timelike-CFJ tachyonic band (m > p/2 ≈ 2×10⁻³² eV at z*). Naive EFT cutoff of the π-matter-photon operator ~20-60 GeV in the window |
| (b) gravitational CS (0907.5159) | their fix is for a fixed external p (diffeomorphism/stress-tensor problem), with the field constraint qRR̃ = 16πG pFF̃ | O(κ²p) coefficient is 10¹⁰³ below LAGEOS (0708.0001) and far below GW birefringence (2101.11153) | **not a repair**: the U(1) identity ∂E = -¼(dp)F̃ is purely photon-sector, and p = ξT is dynamical, so there is no diffeomorphism problem to fix |
| (c1) ln T_sc FF̃ (= 2508.10532 ln Q on FRW) | gauge invariant; β from the endpoint values of ln T_sc | T_sc ≈ 6H² + 2\|∇Φ\|²: ln(T_obs/T_*) = -19.8 (cosmological observer), +25.9 (ground), +11.1 (L2) | **environment dependent**: sign flip and a factor 2.3 between ground and L2 observers. Not the 023 fit model. Applies to ln Q and ln R couplings too |
| (c2) torsion-dragged scalar, -(M²/2)(∂θ - T)² + θFF̃ | massless scalar θ tracking ∫T (θ̇ = 3H on FRW with zero stress energy) | M ≥ 6×10⁶ GeV (CAST), M ≥ 2.6×10⁻³H_inf (SPT-3G A_CB), gravity change (M/M_P)² negligible | **consistent**, but it is an ALP-class model whose rolling rate is set by torsion |

**EC-06d (SNY 0⁺).** SNY Eq. 20 is covariant: T_μ = -(3/2)∂_μ ln|R + 6μ/b| holds exactly off FRW (spin-free matter). So p is closed and the coupling is gauge invariant. The published cases III-V, integrated back, reach z* with Φ/H → 3.00 and need **ξ = (4.2-4.5)×10⁻⁴**. Those histories are torsion-dominated, with **H(z*) = 33-66 × H_ΛCDM**, so they are excluded. In the GR-compatible regime (heavy 0⁺ mode), ε = (9/2)R/Y and δH²/H² = (9μ/4)(R/Y), so **ξ ≈ 3μβ/δ ≳ 1 for a ≤ 1% change of H(z*)**. Torsion is then ~10⁻²H, and the coupling is the Li-Xia-Li ∂ln R type. c_T = 1 is argued, not derived.

## Verdict
- **EC-06c is consistent only with a repair, at a cost.** (b) is not a repair. (c1) fixes U(1) but makes β depend on the local gravitational field at the observer, so it is not the cosmological model fitted in 023. Two survivors:
  - **(a) Stückelberg photon**, costs: (1) f(T) teleparallel gravity, not EC, with f_TT ≠ 0 and the frame set by f(T)'s Lorentz equations (which must dominate the photon's back-reaction); (2) the P-odd non-minimal coupling ξ = 3.6×10⁻⁴, added by hand; (3) a photon mass 7×10⁻²⁸ ≲ m_γ ≲ 10⁻²⁶ eV (about one decade), with a longitudinal mode mixed in by dp; (4) a low naive EFT cutoff, ~20-60 GeV; (5) the β phenomenology is the published one (2508.10532, hep-ph/0611192).
  - **(c2) torsion-dragged scalar**, costs: a new massless scalar with M ≳ max(6×10⁶ GeV, 3×10⁻³H_inf); costs (1), (2) and (5) above. It is an ALP-class model in which torsion only fixes the rolling rate.
- **c_T = 1 holds** in f(T) for any f. The GW-speed test cannot kill EC-06c.
- **EC-06d: closed off FRW (gauge invariant), but the needed O(H) torsion at z* occurs only on SNY histories excluded by H(z*). A GR-compatible SNY history needs ξ ≳ 1 and is the ∂ln R coupling in torsion variables.** Candidate mechanism, effectively dead as "O(H) torsion". It survives only as an ALP-like heavy-scalar coupling with ξ ~ 1.

## Open issues
1. Photon back-reaction on the f(T) frame (L1d and its boost analogue) against f_TT T′ M_P²: the size of f_TT needed for the frame to be gravity-dominated is not computed. The f(T) strong-coupling literature (1808.05565 body, Blagojević-Nester 2006.15303) was not opened.
2. The Stückelberg lower limit uses vacuum degenerate mixing and A < 0.1 (a choice). Plasma (ω_p ≫ m) detunes the mixing and would lower it. Microcausality of Proca + timelike CFJ was not checked. The EFT cutoff is power counting only.
3. The ln T endpoint calculation assumes the weak-field diagonal-tetrad T_sc near Earth. The same issue for 2508.10532 (ln Q) and hep-ph/0611192 (ln R) is asserted from the structure of the coupling, and their texts were not checked.
4. The Proca galactic limit (10⁻²⁶ eV) is not adopted by PDG. If it is dropped, the window extends to 10⁻¹⁸ eV.
5. SNY: matter-only equations (no radiation, BBN not evaluated). The GR-regime δ ≤ 1% is a stated choice, not a fit. The cases are SNY's illustrative parameters, not a fit. c_T = 1 for SNY is argued from Eq. 33, not derived. Spin sources would add axial torsion (P ≠ 0).
6. In (c2), stability of the θ-tracking with the M²T² term in the f(T) Lorentz equation is not checked.
