# 20261120-015 skeptic on iter 3 (tasks 012 and 014)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything below is about candidate mechanisms.

## Question
Try to break:
- **012**: (a) exact shift symmetry of ∂φ·J₅ couplings (m = 0, anomaly shift vs quark-mass phases in θ̄); (b) photon decoupling g_eff/g_anom ≈ (m_φ/m_e)²/96; (c) the Green-Schwarz argument; (d) KK-05: dβ = -c dlnα, MICROSCOPE c ≥ 340, clocks c ≳ 5×10⁴, EDM c ≲ 140; (e) 4D Holst/Nieh-Yan + Maxwell without fermions gives g = 0.
- **014**: AIC/BIC parameter counts, Savage-Dickey prior dependence, whether overlapping data are ever multiplied, and whether `src/claims_expected.csv` values come from scripts (3+ rows rerun).

## Method
Scripts (all run through `tools/pin.py`, 1 thread, each < 2 min; output in `run.log`, 18 PASS, 0 FAIL):

| script | content |
|---|---|
| `formfactor_dispersion.py` | (b) re-derived by a different route from 012: the absorptive part of the pseudoscalar triangle from δ(1 - xyr) on the Feynman simplex, Im P(s) = (2π/s) ln((1+v)/(1-v)), then P(r) rebuilt by an unsubtracted dispersion integral. Gives P(0) = 1 (heavy-fermion cancellation of the anomaly), slope 1/12, curvature 1/90, C = 12·8 = 96. |
| `kk05_check.py` | (d) dβ/dlnα = -c and β(c) in sympy from Z = 1+12c², Y = 16c³; MICROSCOPE and clock bounds recomputed from `data/`; EDM redone with the canonical photon normalisation. |
| `spotcheck_fit.py` | Imports 014 `fit_candidates.py`, runs it with its output redirected to `artifacts_fit/` here, checks three `claims_expected.csv` rows, and recomputes B₁₀ for a free-β model at several prior widths. |
| `common15.py` | data/ loaders. |

Also rerun unchanged: 012 `bounds_uholst.py` and `check_formfactor.py` (key lines in `run.log`). Literature: 1903.12559 abstract (photon couplings come entirely from non-anomalous loops).

## Result

| claim | verdict | reason |
|---|---|---|
| (a) m = 0 for ∂φ·J₅; θ̄ φ-independent | SURVIVES | Action depends on φ only through ∂φ, so the constant shift is an exact symmetry regardless of the anomaly. Under the chiral rotation θ shifts by -2Σc_qφ/f and arg det M by +2Σc_qφ/f; they cancel. Matches 2506.11836 for gluons. |
| (b) g_eff/g_anom = (m_φ/m_e)²/96 | SURVIVES (minor caveat) | Dispersion route gives B₁ = -r/12 - r²/90 and C = 96.000 exactly with electron only; 1.44×10⁻⁷⁸ (6×10⁻³³ eV), 3.99×10⁻⁷⁰ (10⁻²⁸ eV) reproduced. Caveat: 96 assumes the light quarks enter hadronically (π⁰/η mixing, suppressed by (m_φ/m_π)²). A perturbative u,d treatment with current masses would shift C by a few %. Irrelevant to the verdict (70+ orders). |
| (c) GS term restores φFF̃ and φGG̃ + QCD mass; not EC | SURVIVES | dH = c_GS(F∧F + trG∧G) gives a non-derivative -c_GS φX after dualisation; ∫trG∧G is integer, so continuous shift symmetry is broken to discrete and θ̄ depends on φ. Nothing in EC supplies the CS pieces. Hypercharge-only GS = ALP with free g (KK-03-type). |
| (d1) dβ = -c dlnα, β(c) = 2c - arctan(2√3c)/√3 | SURVIVES (conditional) | Re-derived in sympy. It depends only on dY/dZ = 2c, i.e. on the U1 coefficients 12c² and 16c³, which I did not re-derive (only their ratio is consistent with the contortion response K = cW̃ + (4/3)c²W). Sign vs β convention still unchecked (012 open issue 4). |
| (d2) MICROSCOPE c ≥ 340 (Δχ = 0.39 M_P) | SURVIVES | 340.2 recomputed. Δχ is a choice: 133 at Δχ = M_P. |
| (d3) clocks c ≳ 5×10⁴ | SURVIVES | Analytic matter-era estimate (χ̇₀ ≈ H₀Δχ) gives 8.1×10⁴; 012's ODE min over m is 5.1×10⁴ (max 8×10⁵). Assumes the field still rolls today. |
| (d4) EDM estimate c ≲ 140 | **BROKEN** | 012 takes d = (3/4)c·k against the unnormalised field strength. With Z = 1+12c², the canonical photon is √Z·A and the measured e refers to it, so d_phys = (3/4)c k/√(1+12c²) → (3/4)k/√12 = 8.2×10⁻³³ e·cm for c ≫ 1. That's 500× below JILA 4.1×10⁻³⁰ for every c: no upper bound on c. (Same answer from d/e = d_bare/q_bare with the KK charge q = k n/R.) If the electron is KK-neutral (task 004), the EDM experiment does not probe the KK photon at all. The P'-even (4/3)c²W piece gives an MDM Δa_e ≈ 3×10⁻¹⁷ at c = 5×10⁴, not binding. |
| (d5) KK-05 "excluded at estimate level" | **BROKEN** | Follows from (d4). No computed bound closes the window: c ≳ 5×10⁴ is allowed. KK-05 stays alive only as **strongly disfavoured** by structural costs: an explicit non-dynamical P'-odd background, the CMB photon = KK photon (fails for KK-neutral SM fermions), and Z ≈ 3×10¹⁰ times its Maxwell value. |
| (e) 4D Holst/NY + Maxwell, no fermions: g = 0 | SURVIVES | With F = dA the Maxwell term does not contain the connection; torsion is sourced by ∂φ only and its elimination gives terms in φ alone. Tree-level exact; without charged matter there are no loops. With fermions it reduces to (a)/(b). |
| 014 AIC parameter counts | WEAKENED | k_full = 2 for EC-03' and KK-02, but 1 for the ALP. A physical ALP also needs m and φ_i (g fixed or not), so ΔAIC_full penalises the named models by 2 relative to the ALP for bookkeeping, not physics. Use k_full = 2 for the ALP too, or drop the column. ΔAIC_eff (k = 1 everywhere) is fine. |
| 014 BIC | SURVIVES | n = 1 → ln n = 0 is stated; n = 6 joint is labelled. |
| 014 Savage-Dickey: larger B₁₀ for bounded models is a prior artifact | SURVIVES | A free-β model with prior ±0.418° gives B₁₀ = 15.1 (ACT) and 5.46×10³ (ACT+WMAP+Planck), exactly KK-02 k = 2. B₁₀ depends only on prior width and the likelihood mass kept. Minor: EC-03' uses β_max at the best mass; marginalising over m on the plateau would lower its B₁₀ further. |
| 014 shared data never multiplied | SURVIVES | `Joint` raises on overlapping components; fits use act_dr6 or act_wmap_planck separately, never both. The joint BIC adds SPT-3G BB (disjoint). |
| claims_expected.csv values come from scripts | SURVIVES | Rerun from source: penalty_cap27_AWP 8.017 (8.0), penalty_k1_AWP 0.897 (0.9), beta_max_cap27 0.0998 (0.100), KK-05 c_min_MICROSCOPE 340, c_min_clocks 5.13×10⁴, c_max_EDM 144 (140), g_ratio_coeff 96. All match. The c_max_EDM row reproduces faithfully but the physics is wrong (d4). |

## Verdict
- 012 part A survives: EC-03 dead by photon decoupling, the GS route is beyond EC.
- 012 part B: the KK-05 exclusion is **broken**. The EDM bound vanishes once the photon is canonically normalised. KK-05 is not excluded by any computed bound; it is strongly disfavoured by structural requirements (P'-odd non-dynamical background, KK photon as CMB photon, c ≳ 5×10⁴). This doesn't change the program direction: no live mechanism where torsion does the work survives without those costs, but the CLAIMS row must not say "excluded".
- 014 survives, except that ΔAIC_full counts the ALP inconsistently.

## Required fixes (overseer)
1. CLAIMS.md KK-05: drop "EDM estimate c ≲ 140" and "excluded at estimate level". Replace with "no bound from EDM after canonical normalisation (d ≈ 8×10⁻³³ e·cm, c-independent); needs c ≳ 5×10⁴, P'-odd background, KK photon = CMB photon; strongly disfavoured". STATE.md live-candidates line likewise.
2. `src/claims_expected.csv` row KK-05 c_max_EDM: remove, or replace with d_phys(c→∞) once a theorist rederives it.
3. 012 `bounds_uholst.py` B4/B5: include 1/√(Z₀(1+12c²)) (to be done by a theorist, not here).
4. 014: set k_full = 2 for the ALP or drop ΔAIC_full.
5. Optional: 012 F5 text should say C = 96 assumes hadronic light quarks.

## Open issues
- The U1 coefficients 12c² and 16c³ (50-component contortion elimination) were not independently re-derived; dβ = -c dlnα rests on them.
- The EDM argument assumes Z₀ = 1 for k = √(2κ₄) and that the electron's charge is the KK charge. For a KK-neutral electron the EDM bound does not apply at all; either way it does not bound c.
- Quevillon-Smith 1903.12559 read at abstract level only; the P(0) = 1 normalisation is checked as a property of the triangle integral, and the identification of the "1" with the anomaly is the standard low-energy theorem.
- `artifacts_fit/fit_results.csv` here is the regenerated 014 output (identical numbers), kept for traceability.

## Rerun
```
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-015-skeptic-iter3/formfactor_dispersion.py
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-015-skeptic-iter3/kk05_check.py
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-015-skeptic-iter3/spotcheck_fit.py
```
