# State

Iteration: 1 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt`.
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. `tools/repin.ps1` ran at start and after subagents returned in iter 1; it had nothing to fix. `tools/pin.py` fixed in iter 1 to resolve relative executables (shutil.which).

## Iteration 1 tasks (all in `work/`)
- 001 lit-target-papers: the four papers read in full (HTML v2). `literature/<id>.md`, `prior_measurements.md`, `conventions_beta.md`, `data/birefringence_measurements.csv`.
- 002 lit-prior-torsion: prior work catalogue (`literature/prior_work_torsion_birefringence.md`, `open_gaps.md`). About 10 papers read in full, about 70 abstract-only (tagged [ABS]).
- 003 theory-ec-photon: zero-birefringence for minimal EC, FRW torsion, β per candidate. 47/47 checks PASS (rerun by overseer).
- 004 kk-reduction: 5D EC on S¹ and S¹/Z₂. 13/13 checks PASS (rerun by overseer).
- 005 bounds: `data/bounds.csv`, `literature/bounds.md`, conversion script reproduces the table.

## Settled (iter 1, checked symbolically; not yet skeptic-audited)
- Data target: ACT DR6 β = 0.215 ± 0.074° (2509.13654 Table 3); ACT+WMAP+Planck 0.264 ± 0.058° (0.231 ± 0.059° with Cosmoglobe). β assumed achromatic; PA5 f150 − f090 difference 0.21 ± 0.09°. SPT A_CB < 1.2×10⁻⁴ (95%), < 0.53×10⁻⁴ with lensing prior. SPT is blind to isotropic β.
- Conventions: β = +½∫p_τ dτ for CFJ vector p; L ⊃ -(g/4)φFF̃ gives β = -(g/2)(φ_0 - φ_LSS) in our conventions. CONVENTIONS.md corrected (contortion, CFJ sign).
- Minimal 4D EC + minimally coupled Maxwell: β = 0 exactly (F = dA, conformal invariance on FRW). Building F with the torsionful ∇ breaks U(1). EC + Dirac: S^μ = -(κ/2)J₅^μ, L_eff = -(3κ/16)J₅², no photon coupling at tree level.
- FRW torsion: only trace h (T_τ = 3h) and axial f (S^τ = -f/a²).
- Pure 5D EC-KK on S¹ or S¹/Z₂: no parity-odd photon coupling at tree level, for any κ₅, R, radion. Torsion source gives only four-fermion terms and a Planck-suppressed Pauli term (linear in F, no propagation effect). No 5D Holst term exists (e∧R∧T is exact). On S¹/Z₂ the KK photon is projected out.
- Needed size: a CFJ term needs |p_0| ≈ 1.1×10⁻³⁵ eV (≈ 8×10⁻³ H₀) if constant.
- Lab bound on cosmological axial torsion coupled minimally to fermions: |S_0| < 7.3×10⁻²² eV (neutron b_⊥, arXiv:1312.3225, via Sun-CMB boost). Lab bounds alone do not exclude a CFJ-type torsion route with coupling ξ ≳ 10⁻¹⁴.

## Live candidate mechanisms
- EC-01 non-minimal CFJ torsion-photon term -½ξS·AF̃: lives only with a non-minimal torsion source; in minimal EC the source is ~28 orders short (see Dead).
- EC-03 torsion pseudoscalar (Barbero-Immirzi field / Nieh-Yan / Kalb-Ramond torsion axion) × FF̃: needs Δφ ≈ 1.8×10¹⁶ GeV for g ~ 1/M_P. Threats: loop-generated photon coupling disputed; torsion couples to quarks, so a QCD mass (~10⁻¹² eV, 1502.03694) and washout (2512.21888) are likely.
- EC-04 dynamical axial torsion (Poincaré gauge theory): not yet examined.
- KK-02 bulk U(1) + 5D Chern-Simons term, Wilson-line axion: g_φγ = 24πκ_CS α/f_a. An ALP in disguise; torsion plays no role. Needs α-variation, KK-scale, washout bounds.
- KK-03 Kalb-Ramond 2-form + Green-Schwarz term: g = √6|λ|g₄²M_P. Not yet confronted with bounds.

## Dead mechanisms (preliminary, pending skeptic)
- EC-00 minimal EC + minimal Maxwell: zero birefringence (theorem, 003).
- EC-02 minimal EC with fermion spin as torsion source + CFJ term: relic-neutrino asymmetry gives |S_0| ~ 2×10⁻⁶⁸ eV today, about 28 orders short at ξ = 1 (003, order-of-magnitude, neutrino asymmetry bound is a placeholder, needs a cited value).
- KK-01 pure EC-KK: no parity-odd photon coupling at tree level (004).

## Open questions
- Fermion-loop CFJ from axial torsion (EC-05): coefficient is regularization dependent (hep-ph/9901358); does Reig arXiv:2611.07578 close it? Order estimate: ~35 orders short (002, unchecked).
- Closest KK prior art: arXiv:2608.09299 (5D teleparallel KK, radion × FF̃, fit to Planck PR4). Our EC-KK finds no tree-level radion FF̃: why does teleparallel differ?
- Torsion axion: is the QCD mass unavoidable, and which extension avoids it?
- Data-side: unresolved 4.47x discrepancy between the Data Tables β→k_AF entry and Kostelecky-Mewes Eq. 158. Planck axion isocurvature bound not verified. No BBN bound on axial torsion found.
- Which public likelihoods can be run locally under the 2-CPU cap (ACT DR6, SPT candl_data, Planck PR4).

## Next planned tasks (iter 2)
1. Skeptic on 003 + 004 (zero-coupling theorems, CFJ sign, KK no-go).
2. Theorist: turn EC-02 and EC-05 into a clean no-go with cited inputs (neutrino asymmetry bound, chiral DM densities, Reig 2611.07578 applied to the loop coefficient), shortfall in orders of magnitude. Examine EC-04 (Poincaré gauge theory dynamical axial torsion).
3. Cosmologist/theorist: torsion-axion route EC-03: QCD mass from the torsion-quark coupling, washout, anisotropy, g_aγγ bounds; what is needed to evade.
4. KK theorist: compare with 2608.09299; KK-02/KK-03 vs α-variation (Oklo, clocks), radion stabilization, KK-scale bounds, washout; S¹/Z₂ brane-photon variant.
5. Data analyst (light): Gaussian β likelihood from published numbers plus SPT candl A_CB if feasible on 2 CPUs.
