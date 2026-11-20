# State

Iteration: 7 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 7: `repin.ps1` re-pinned 1 python.exe mid-iteration (launcher-stub pattern).
- `artifacts/` (gitignored): candl_data, reproduce outputs, `clean_checkout/` (safe to delete).
- Path leak in commit 87caa4e (iter 3) remains in history; path-scan gate before every commit since.
- `python -m src.reproduce` (~3 min): 18 validation + all CLAIMS numbers incl. EC-06 (exit 0 after iter-7 audit). `--checks` (~16 min) also reruns 7 symbolic run_checks scripts (003, 004, 009, 012, 016, 019, 022; 025 picked up automatically). Not yet covered: 025's numbers (Stückelberg window, EFT cutoff, SNY H(z*) ratio).

## Iteration 7 tasks
- 025 theory: f(T) Lorentz scalar ζ is fixed non-dynamically and dT ≠ 0 on solutions, so U(1) breaks everywhere; c_T = 1 in f(T); gravitational CS is not a repair; surviving repair is a Stückelberg photon with 7×10⁻²⁸ ≲ m_γ ≲ 10⁻²⁶ eV, naive cutoff ~20-60 GeV; SNY (EC-06d) closed p but needed histories excluded by H(z*). 70/70 PASS (rerun).
- 026 skeptic on 022-024: nothing broken. MZ bound corrected to ξ ≥ 1.4×10⁴ and verdict softened to "disfavoured, not excluded"; ζ missing from 022's gauge formula (addressed by 025); TEGR wording; 23-31% reionization share; A_CB labelled generic; 1003.4823 CL aligned in data/ (overseer fixed).
- 027 data: EC-06 numbers in reproduce; `--checks` option.
- 018 draft: revision 2 in `work/20261120-018-result-draft/RESULT_DRAFT.md` with all iter 5-7 findings; no PENDING items.

## Settled
- 4D Einstein-Cartan proper: no β from minimal EC, spin-sourced torsion (any coupling), or torsion pseudoscalars with ∂φ·J₅ (EC-00, -02, -03, -05).
- Unsourced O(H) torsion: MZ disfavoured (ξ ≥ 1.4×10⁴), TEGR undefined, SNY reduces to the ∂ ln R coupling.
- KK track: Result B (KK-01, -04, -05 dead; KK-02/03 ALPs where torsion plays no role).

## Live candidate mechanisms
- EC-06c f(T) trace-torsion CFJ (not EC): ξ = 3.6×10⁻⁴; costs listed in CLAIMS.md and the draft. Needs skeptic review of the photon-mass window and the 20-60 GeV cutoff (a cutoff below collider energies for the photon sector would be a serious problem).
- ALP class (KK-02, KK-03, EC-03', EC-04).

## Open questions
- EC-06c: is the naive EFT cutoff ~20-60 GeV in conflict with collider photon physics (then EC-06c is dead); plasma effect on the Stückelberg lower edge; Proca + timelike CFJ microcausality; photon back-reaction on the f(T) frame.
- KK-05 brane-photon mixing loophole (low priority).

## Direction toward a stop
Both tracks have a draft answer: KK Result B; 4D Result B for EC proper with a Result C only in f(T) teleparallel gravity (or dead if the cutoff kills EC-06c). Remaining: skeptic on 025; a fresh final skeptic on draft revision 2 for one full iteration; 025 numbers into reproduce. If the final skeptic finds no blocker, write RESULT.md.

## Next planned tasks (iter 8)
1. Skeptic on 025 (ζ solution, dT ≠ 0, Stückelberg window, EFT cutoff vs collider bounds, SNY H(z*) exclusion).
2. Fresh final skeptic (new agent, did none of the work) on RESULT_DRAFT.md rev 2 + CLAIMS.md + repo, against the OVERSEER_PROMPT stop criteria.
3. Data analyst: add 025 numbers to reproduce; run `--checks` once on the clean clone.
