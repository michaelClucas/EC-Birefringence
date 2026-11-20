# State

Iteration: 9 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 9: `repin.ps1` re-pinned 1 python.exe and 2 agent bash shells mid-iteration (shells match the repo path in their command line; harmless).
- `artifacts/` (gitignored): candl_data, reproduce outputs, `clean_checkout/` (safe to delete).
- Path leak in commit 87caa4e (iter 3) remains in history; path-scan gate before every commit since. 033 fixed a reproduce crash path that printed "OK" and leaked a path.
- `python -m src.reproduce` (~4-7 min): 18 validation + 66/66 CLAIMS numbers + draft-only numbers (`src/draft_expected.csv`, checked against RESULT_DRAFT.md); exit 0 after iter-9 audit. `--checks`: symbolic check scripts (8 as of iter 8; 031 and 032 have run_checks.py too). Not regenerated: 020 higher-dim Λ scales, 0908.0629 4.9σ, LiteBIRD 3σ, KK-05 10²¹ m_e.

## Iteration 9 tasks
- 031 theory: EC-06a (MZ) dead: CFJ vector not closed in matter, Stückelberg photon needed, strong coupling at 3.6 keV (9.4 orders in energy below 10 TeV); no off-tracker escape; GW bound was outside 2001.06373's regime (respecting it needs ξ ≥ 4.7×10¹⁶). 36/36 PASS (rerun). Open: tensor torsion in the Λ equation, local physics near masses, one input (24 Hz) from a search snippet.
- 032 literature/bounds: 028's compact-star inputs verified into data/compact_star_inputs.csv (Sirius B radius citation was wrong; correct source 1806.00061; no verdict change); cited H(z*) bound (δH/H < 0.18-0.49%, 95%) gives EC-06d 3.8-4.5 orders short with radiation-corrected ratio 29-57; KK-04 4.34 orders vs Cassini; 2610.02667 added (no novelty change). PASS (rerun).
- 033 data: reproduce gaps closed; draft-number checking added.
- 018: draft revision 3 written: both tracks Result B for torsion, ALP statement explicit, shortfall tables per kill.

## Settled
- 4D: every torsion route examined is dead with a named reason and shortfall (see draft rev 3 tables); EC-06a pending skeptic.
- KK: KK-01, -04, -05 dead with shortfalls.
- Survivors are ALP-class with hand-added couplings (EC-03', EC-04, KK-02, KK-03); torsion plays no role.

## Open questions
- Skeptic review of 031 (EC-06a strong-coupling kill; non-closure in matter; k* regime argument; the 24 Hz input).
- KK-05 brane-photon mixing loophole (low priority).

## Direction toward a stop
Draft rev 3 states Result B for both tracks. Remaining: skeptic on 031, and a fresh final skeptic on draft rev 3 (one full iteration, the stop condition). If neither finds a blocker, write RESULT.md, commit, push, send a push notification, end the loop.

## Next planned tasks (iter 10)
1. Skeptic on 031.
2. Fresh final skeptic (new agent, none of the work, no earlier skeptic READMEs until its own list is formed) on RESULT_DRAFT.md rev 3 + CLAIMS.md + repo, judging the draft under both outcomes of (1).
