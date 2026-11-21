# State

Iteration: 11 (completed 2026-11-21)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 11: `repin.ps1` re-pinned 1 python.exe mid-iteration.
- Path leak in commit 87caa4e (iter 3) remains in history; path-scan gate before every commit since.
- `python -m src.reproduce`: 18 validation + 74/74 CLAIMS numbers + 5/5 draft numbers, exit 0 after iter-11 audit. `--checks` auto-discovers run_checks.py in folders > 035 plus 031, 032, 034, 035 (full --checks not rerun this iteration). 036/037 numbers not yet in the reproduce tables.
- An iteration-7 data agent (task 027, work already committed) reported a stream stall on 2026-11-21; no work lost.

## Iteration 11 tasks
- 036 theory: EC-06d heavy-0⁺ branch dead at every 0⁺ mass (unscreened α ≥ 1/3 coupling to the matter trace: Cassini 3.9 orders for λ ≳ AU; local observer endpoint for shorter λ gives opposite-sign ACT vs Planck/WMAP; 23.7 orders short in size at the HB bound). Not the same model as hep-ph/0611192 (CLAIMS corrected). 36/36 PASS (rerun).
- 037 theory: shortfalls at the coupling bound; stronger dual-axion bound ξ < 1.5×10¹⁵ (Chandra × red giant); vector-current route dead via CMS dijet; every large-ξ bound is completion dependent; weakest corner (Stückelberg, weakest m_γ, 1 MeV cutoff) leaves the Pauli route 1.2 orders short. PASS (rerun).
- 038 data: 031/032/034/035 numbers into reproduce; --checks auto-discovery.
- 018: draft revision 4 (Result C both tracks, cost table, kill tables with shortfalls at the bound, Not examined section, 029/035 leftover fixes). README.md now carries the systematics caveat.

## Settled
- Every photon coupling EC torsion generates by itself fails (4D and KK); every unsourced-torsion route examined is dead (EC-06a-d).
- β is reachable only by ALP-like fields with couplings added beyond EC (EC-03', EC-04, KK-02, KK-03, EC-06a f(Λ) remnant) at the costs in the draft's cost table.

## Open questions
- Skeptic review of 036 and 037.
- 037 weak corner: Stückelberg completion of the C J₅·AF̃ operator; its strong-coupling scale (as computed for EC-06a/c) was not computed and would likely close the 1.2-order corner.
- KK-05 brane-photon mixing loophole (listed under Not examined).

## Next planned tasks (iter 12)
1. Skeptic on 036 + 037, including the Stückelberg strong-coupling scale for the 037 weak corner.
2. Fresh final skeptic on draft rev 4 (new agent; judge under both outcomes of 1).
3. Data: add 036/037 numbers to reproduce.
If the final skeptic finds no blocker, write RESULT.md.
