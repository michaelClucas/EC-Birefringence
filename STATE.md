# State

Iteration: 3 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 3: `repin.ps1` re-pinned 1 python.exe mid-iteration (the `.venv` launcher stub pattern from iter 2); the skeptic reported one sub-second unpinned file-edit helper. No sustained unpinned load.
- `artifacts/` (gitignored): SPT candl_data and reproduce outputs.
- iter 3: commit 87caa4e included a run.log with an absolute local path (from a candl SyntaxWarning). Scrubbed in the next commit; it remains in git history (no history rewrites). Overseer now runs the path scan as a gate before every commit.
- `python -m src.reproduce` (~2 min, 1 CPU): published β/A_CB validation (18 PASS) and every number in CLAIMS.md's last two columns (expected values in `src/claims_expected.csv`; EC-03 QCD-mass rows kept as "superseded"). Exit 0 after iter-3 audit. Newer numbers (EC-03' 2.7% cap, fit penalties, KK-05) are not yet in the reproduce table.

## Iteration 3 tasks
- 011 skeptic on 007-009: EC-02, EC-05, EC-04, washout, P', KK-04, KK-02/03 survive. 008's QCD-mass reasoning BROKEN (∂φ·J₅ is exactly shift symmetric). EC-03' broken as a torsion mechanism. Scripts 25/25 PASS (rerun).
- 012 theory: shift symmetry confirmed (m = 0; agrees with 2506.11836 for gluons); photon coupling decouples as (m_φ/m_e)²/96 (new, not in 2506.11836); Green-Schwarz term restores φFF̃ but also φGG̃ and the QCD mass; u-Holst with a non-radion scalar = KK-05, excluded at estimate level. 37/37 PASS (rerun).
- 013 literature: no-go only partly published (hep-th/0207165 for KR torsion parity; 2310.09152 EFT no-go without gravity, which minimal EC falls under after torsion elimination). KK P' all-orders argument, shortfall numbers, radion Cassini kill not found in literature. Cited ultralight caps in `data/ultralight_fraction_bounds.csv`.
- 014 data: reproduce extended; fits of axion-class candidates. No candidate beats a free-β ALP; SPT A_CB does not discriminate.

## Settled
- All iter-1/2 items listed in CLAIMS.md as skeptic-survived.
- Minimal EC and EC with minimally coupled fermions give no isotropic β: EC-00 exact; EC-02 28.5 orders short; EC-05 ≥ 29 orders short. After torsion elimination minimal EC is GR + LEFT-type four-fermion operators, so 2310.09152's no-go applies.
- Torsion pseudoscalars coupled via ∂φ·J₅ (KR, BI×Nieh-Yan, PGT 0⁻, Holst-BI with minimal fermions) are exactly shift symmetric: no QCD mass, and the photon coupling decouples (~10⁻⁷⁰). EC-03 dead.
- Any working route needs a photon coupling added beyond EC (by hand, Green-Schwarz term, or 5D Chern-Simons term); then the field is an ALP and torsion at most sets f.
- Pure EC-KK: no parity-odd photon coupling at any order (P'). Radion × FF̃ excluded (Cassini, ISL).
- Fits (ACT+WMAP+Planck): free-β ALP Δχ² = 20.7 vs null; EC-03' at 2.7% cap β ≤ 0.100°, penalty 8.0; KK-02 k = 1 β ≤ 0.209°, penalty 0.9.

## Live candidate mechanisms
None in which torsion does the work.
- KK-02, KK-03: ALPs in disguise, viable as ALPs; torsion plays no role.
- EC-03', EC-04: ALPs with a hand-set coupling; EC-03' disfavoured by the cited energy cap.
- KK-05: the only route found where torsion generates the parity-odd photon coupling; excluded at estimate level (EDM estimate c ≲ 140 vs MICROSCOPE c ≥ 340); needs a non-dynamical P'-odd background and the KK photon as the CMB photon. Skeptic pending; EDM bound is an estimate, needs firming.

## Dead mechanisms
EC-00, EC-02, EC-03, EC-05, KK-01, KK-04 (reasons and numbers in CLAIMS.md).

## Open questions
- KK-05: firm up the EDM bound (cited, with the coupling of χ to the electron EDM operator derived); is there any P'-odd background that is non-dynamical and does not give the KK photon a mass?
- CS level normalization (k ∈ Z vs 6Z) for KK-02, quoted not derived.
- Minor: Data Tables β→k_AF factor 4.47; Planck isocurvature bound; no BBN bound on axial torsion. None affects a verdict.

## Direction toward a stop
Both tracks are at a Result B structure with a C-type footnote: Einstein-Cartan torsion (minimal or with the named extensions) does not produce the observed β within bounds; the effect is possible only by adding an axion-like field with a photon coupling put in beyond EC, where torsion plays no role (cost: new field + coupling + for EC-03' an energy fraction above the cited cap). Remaining before RESULT.md: skeptic on 012 (shift-symmetry photon decoupling, GS, KK-05); firm KK-05 EDM bound; reproduce covering iter-3 numbers; then a fresh skeptic on a RESULT draft for one full iteration.

## Next planned tasks (iter 4)
1. Skeptic on 012 and 014 (photon decoupling formula, GS argument, KK-05 coefficient c_odd and dβ = -c dlnα, EDM estimate, fit parameter counts).
2. Theorist: KK-05 EDM bound derived and cited; P'-odd background question; CS level normalization for KK-02.
3. Data analyst: add iter-3 numbers (EC-03 g_eff ratio, EC-03' 2.7%/6.9% β_max, fit penalties, KK-05 c bounds) to `src/claims_expected.csv` and reproduce.
4. Overseer: draft RESULT candidate text in `work/` (not RESULT.md) for the iter-5 skeptic.
