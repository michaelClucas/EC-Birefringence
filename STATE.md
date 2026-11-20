# State

Iteration: 5 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 5: `repin.ps1` re-pinned 2 python.exe during the clean-checkout run (launcher-stub pattern); nothing at start or after return.
- `artifacts/` (gitignored, ~650 MB): SPT candl_data, reproduce outputs, `clean_checkout/` clone from task 021 (safe to delete).
- Path leak in commit 87caa4e (iter 3) remains in history; scrubbed since. Path-scan gate before every commit.
- `python -m src.reproduce` (~3 min main checkout, ~8 min clean clone; run with `-W ignore::SyntaxWarning` through `tools/pin.py`): 18 validation PASS + 36/36 live CLAIMS numbers. Verified on a clean clone (021): no dependence on uncommitted files. Not covered: the symbolic checks behind EC-00 and KK-01 (run separately via the work/ `run_checks.py` files).

## Iteration 5 tasks
- 019 adjudicator (KK-05 EDM): the P'-even pieces are the EDM (016 right, 012 and 015 wrong on assignment; 012 also off by 2 and missing Z). JILA bound c ≤ 188 holds only for a KK-neutral electron with a non-geometric charge, which breaks the KK gauge symmetry. With the consistent KK-momentum charge: no EDM bound, but the electron mass is ~M_P scale. 9/9 PASS (rerun).
- 020 final skeptic on the RESULT draft: numbers, signs, units, cited bounds (2301.08361, 2409.08280, Chandra) and reproduce verified. BLOCKER: the 4D no-go only covers spin-sourced torsion; unsourced O(H) torsion (EC with varying Λ, teleparallel FRW torsion, PGT torsion cosmologies) with a CFJ coupling needs only ξε ≈ 10⁻³. MAJOR: KK verdict logic inconsistent (now resolved: KK-05 dead by charge problem); literature novelty overstated (photon decoupling of derivative-coupled ALPs is standard: 1708.00443, 1903.12559, 2011.10025; only the torsion application is new; draft omitted 0908.0629, a torsion-photon birefringence fit). 8 MINOR items in its README. 6/6 PASS (rerun).
- 021 data: AIC k_full fixed (free ALP = 2); clean-checkout reproduce passes.

## Settled
- 4D: no birefringence from minimal EC or from spin-sourced torsion (EC-00, EC-02, EC-05); torsion pseudoscalars with ∂φ·J₅ couplings have no photon coupling (EC-03); ALP routes need a coupling added beyond EC (EC-03', EC-04).
- KK: pure EC-KK gives nothing at any order (KK-01); radion × FF̃ excluded (KK-04); KK-05 (torsion-generated χFF̃ via u-Holst) dead because no light state can carry KK charge e; KK-02/KK-03 are ALPs in which torsion plays no role.
- Fits: no candidate beats a free-β ALP.

## Live / unexamined
- EC-06 (NEW, unexamined): background torsion of order H not tied to spin (Magueijo-Złośnik EC with varying Λ 1908.05184; teleparallel; PGT torsion cosmologies) + a non-minimal CFJ coupling ξ. Needs ξε ≈ 1.1×10⁻³ (skeptic estimate). Questions: is such a background a solution of EC or only of an extension; gauge invariance of the CFJ term beyond FRW (p must be closed); anisotropy from torsion perturbations vs SPT A_CB; β(z) history (reionization vs recombination) vs ALP; what bounds ξ (photon-sector LV, the torsion's own coupling to fermions in that theory, solar-system torsion); prior art (0908.0629 Das-Mohanty-Prasanna; open gap G10 in literature/open_gaps.md).
- KK-02, KK-03, EC-03', EC-04: ALPs.
- KK-05 loophole: brane photon mixing with the KK photon (low priority).

## Direction toward a stop
KK track is ready for Result B (with KK-02/03 as ALPs where torsion plays no role) after the draft is corrected. The 4D track cannot stop until EC-06 is settled: it is either a Result C (a non-spin torsion background plus an added photon coupling, with stated costs) or a no-go with a named reason.

## Next planned tasks (iter 6)
1. Theorist: EC-06 class: derive the background torsion in each named theory (Magueijo-Złośnik, teleparallel/TEGR and f(T), PGT torsion cosmologies), which components are nonzero (axial vs trace) on FRW, the CFJ couplings allowed and their gauge invariance off-FRW, β formula, and what independently bounds ξ.
2. Cosmologist: β(z) and β for parameterized backgrounds S_0 = εH, T_τ = ε′·aH etc.; reionization vs recombination split; anisotropy from torsion perturbations vs SPT; fit to the β likelihood with Δχ² vs null and vs free ALP.
3. Literature scout: prior art on non-spin cosmological torsion + photon couplings (0908.0629 and citing papers, Magueijo-Złośnik follow-ups, teleparallel birefringence papers), and fix the novelty list per skeptic 020 MAJOR 4.
4. Overseer: update RESULT draft (KK verdict, novelty list, MINOR items from 020) after 1-3.
