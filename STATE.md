# State

Iteration: 6 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 6: `repin.ps1` found nothing to fix.
- `artifacts/` (gitignored): candl_data, reproduce outputs, `clean_checkout/` (safe to delete).
- Path leak in commit 87caa4e (iter 3) remains in history; path-scan gate before every commit since.
- Rule slip iter 6: task 023 reported briefly creating an empty file outside the repo by mistake and deleting it.
- `python -m src.reproduce`: 18 validation + 36/36 CLAIMS numbers (verified on a clean clone in iter 5). EC-06 numbers are NOT yet in the reproduce table.
- Note: `work/20261120-024-lit-ec06-novelty/kaf_lowz.py` must be run from the repo root.

## Iteration 6 tasks (EC-06: unsourced O(H) torsion + CFJ)
- 022 theory: MZ torsion dead (GW170817 → ξ ≥ 7×10³); TEGR torsion is a frame choice (dead); f(T) trace torsion live with costs (ξ = 3.6×10⁻⁴, U(1) breaks off FRW); PGT 0⁺ open. 54/54 PASS (rerun).
- 023 cosmology: ξε = 1.073×10⁻³ for p ∝ H; for that shape only 23-31% of β accrues after reionization (the brief's premise was wrong); a single ALP mass reproduces any one reionization/recombination ratio, so the ratio tests a shape, not ALPs as a class; LiteBIRD-class 0.05° would separate n = 0 from a late ALP at ~3σ. Fit = free-β fit relabelled (Δχ² 8.44 / 20.72). A_CB ≤ 5×10⁻¹¹. Reproduces 2508.10532's coupling map. 104/104 PASS (rerun).
- 024 literature: EC-06 phenomenology published without torsion (2508.10532 ln(6H²); hep-ph/0611192 ln R); torsion origin not found. Low-z k_AF bounds in `data/kaf_lowz_bounds.csv` (best 1.7×10⁻³⁴ eV, 1003.4823) don't test EC-06. Novelty audit in `literature/novelty_audit.md`: each "new" item is a known general result applied for the first time, or new numbers inside known no-gos. PASS (rerun).

## Settled
- 4D EC proper (minimal, spin-sourced, fermion-loop, torsion pseudoscalars with ∂φ·J₅): no β (EC-00, -02, -03, -05), skeptic-survived.
- KK: pure EC-KK no β at any order (KK-01); radion × FF̃ excluded (KK-04); KK-05 dead (no light KK-charged electron); KK-02/03 are ALPs, torsion plays no role.
- Unsourced O(H) torsion: MZ (EC extension) dead by GW speed; TEGR torsion unphysical.

## Live candidate mechanisms
- EC-06c f(T) trace-torsion CFJ: needs ξ = 3.6×10⁻⁴. Costs: theory is f(T) teleparallel gravity, not EC (local Lorentz broken); non-minimal photon coupling; U(1) gauge invariance broken off FRW unless a photon mass or gravitational CS counterterm is added (need to check those against photon-mass bounds and GW/CS bounds); phenomenology identical to published torsion-free models. GW-speed test for f(T) background not done. Skeptic pending.
- EC-06d PGT 0⁺ (SNY): closed p on FRW; size at z* and off-FRW closure unknown.
- ALP-class (KK-02, KK-03, EC-03', EC-04): torsion plays no role.

## Open questions
- EC-06c: does the gauge-invariance repair survive bounds (Stückelberg photon mass needed vs m_γ < ~10⁻¹⁸ eV; gravitational CS coefficient vs GW birefringence / Gravity Probe B)? GW speed in the f(T) background?
- EC-06d: evaluate ε at recombination and off-FRW closure.
- KK-05 brane-photon mixing loophole (low priority).

## Direction toward a stop
KK track: ready for Result B. 4D track: Result B for Einstein-Cartan proper and its named extensions; a Result C in a teleparallel f(T) theory (not EC) at the cost of a non-minimal coupling and broken U(1) off FRW, unless iter 7 kills it. The RESULT draft needs: KK verdict fix, corrected novelty list (literature/novelty_audit.md), EC-06 section, 020's MINOR items.

## Next planned tasks (iter 7)
1. Theorist: EC-06c gauge-invariance repair costs (photon mass, gravitational CS) vs bounds; GW speed in f(T) with the needed background; EC-06d ε(z*) and off-FRW closure.
2. Skeptic on 022 + 023 (EC-06 verdicts, the f(T) ξ, the 31% reionization share, MZ GW bound).
3. Data analyst: add EC-06 numbers to `src/claims_expected.csv` and reproduce.
4. Overseer: revise RESULT draft.
