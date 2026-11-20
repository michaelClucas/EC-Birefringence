# State

Iteration: 4 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 4: `repin.ps1` found nothing to fix at start and after return.
- `artifacts/` (gitignored): SPT candl_data and reproduce outputs.
- iter 3: commit 87caa4e included a run.log with an absolute local path (from a candl SyntaxWarning). Scrubbed in e239b04; it remains in git history (no history rewrites). Since then: overseer runs the path scan as a gate before every commit, and `src.reproduce` filters the warning and scrubs paths from its own output (017).
- `python -m src.reproduce` (~3 min, 1 CPU; run with `-W ignore::SyntaxWarning` through `tools/pin.py`): published β/A_CB validation (18 PASS) and all 36 live numbers in CLAIMS.md (expected values in `src/claims_expected.csv`), plus fit summary. Exit 0 after iter-4 audit.

## Iteration 4 tasks
- 015 skeptic on 012/014: shift symmetry, photon decoupling (C = 96 re-derived by dispersion integral), Green-Schwarz, MICROSCOPE/clock bounds, 4D Holst/NY g = 0 survive. Claimed 012's EDM bound c ≲ 140 BROKEN (photon normalisation Z omitted; canonical EDM saturates below JILA, so no bound). 014 AIC k_full column weakened (free ALP should also count 2 parameters, or drop the column). Scripts 18/18 PASS (rerun).
- 016 KK-05 EDM + CS level: re-derived Z = 1+12c², Y = 16c³ exactly; finds 012 and 015 both assigned the dipole pieces wrongly: the P'-even pieces are the EDM and grow linearly in c, so JILA gives c ≤ 188 for a bulk Dirac electron with P'-even mass, excluding the minimal form by a factor 273 vs clocks (1.8 vs MICROSCOPE at Δχ = 0.39 M_P). Escapes: P'-odd electron mass (≥ 99.6%), non-bulk electron, Δχ > 0.71 M_P. A fixed fibre vector is a P'-odd background that keeps the KK photon massless and 4D Lorentz invariance (cost: explicit breaking of y-dependent 5D diffeos). CS level k ∈ Z with mixed gauge-gravitational CS term (6Z without). 28/28 PASS (rerun).
- 017 reproduce: all CLAIMS numbers covered; path scrubbing. Overseer re-pointed the KK-05 EDM row from 012 to 016 (c ≤ 188).
- 018 overseer: RESULT draft in `work/20261120-018-result-draft/RESULT_DRAFT.md`.

## Settled
- Minimal EC and EC with minimally coupled fermions give no isotropic β (EC-00, EC-02, EC-05), skeptic-survived.
- Torsion pseudoscalars coupled via ∂φ·J₅ are exactly shift symmetric: no QCD mass, photon coupling (m_φ/m_e)²/96 (EC-03 dead), skeptic-survived twice.
- Any ALP route needs a photon coupling added beyond EC; Green-Schwarz completion brings back the QCD mass.
- Pure EC-KK has no parity-odd photon coupling at any order; radion × FF̃ excluded (KK-01, KK-04).
- KK-02: CS level |k| = 1 allowed with the mixed CS term; β ≤ 0.209°.
- Fits: no candidate beats a free-β ALP; SPT A_CB does not discriminate.

## Live candidate mechanisms
- KK-02, KK-03: ALPs in disguise, torsion plays no role.
- EC-03', EC-04: ALPs with a hand-set coupling; EC-03' disfavoured by the cited energy cap.
- KK-05: the only route where torsion generates the coupling. DISPUTED: 016 says excluded in minimal form (EDM), 015 said no EDM bound. Strongly disfavoured either way by structural costs (fixed P'-odd background breaking y-dependent diffeos, KK photon = CMB photon, c ≳ 5×10⁴).

## Dead mechanisms
EC-00, EC-02, EC-03, EC-05, KK-01, KK-04.

## Open questions
- KK-05 EDM dispute: which dipole structure is P'-even for the KK photon (axial under P') and whether the -(4/3)c²W term contributes to the EDM. Needs a third, independent adjudication.
- 014 AIC k_full: set free-ALP k_full = 2 or drop the column (minor, no verdict depends on it).
- KK-05 non-minimal escapes: are the escapes natural or tuned (P'-odd electron mass ≥ 99.6%)?
- Minor, no verdict depends on them: Data Tables β→k_AF factor 4.47; Planck isocurvature bound; no BBN bound on axial torsion; S¹/S¹/Z₂ only (warped not treated).

## Direction toward a stop
4D track: Result B (no-go) with a C footnote (only ALPs with couplings added beyond EC). KK track: Result B if KK-05 is excluded, Result C (torsion-generated coupling possible only with a fixed P'-odd background, non-minimal fermion sector, c ≳ 5×10⁴) if it survives in non-minimal form. The draft states both cases.

## Next planned tasks (iter 5)
1. Adjudicator (fresh agent, did neither 015 nor 016): settle the KK-05 EDM dispute with an independent derivation.
2. Fresh skeptic on `work/20261120-018-result-draft/RESULT_DRAFT.md` + CLAIMS.md + the repo, for one full iteration (the stop condition's skeptic).
3. Data analyst: fix the 014 AIC column; confirm reproduce on a clean checkout (fresh clone into artifacts/, run download + reproduce).
