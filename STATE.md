# State

Iteration: 10 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 10: `repin.ps1` re-pinned 1 python.exe mid-iteration.
- Path leak in commit 87caa4e (iter 3) remains in history; path-scan gate before every commit since.
- `python -m src.reproduce`: 18 validation + 66/66 CLAIMS numbers + 4/4 draft numbers; `--checks` 8/8 (skeptic 035 run, exit 0). Not regenerated: 031's 3.6 keV and 9.4 orders, KK-04's 3.6-4.1, EC-03′ best-fit mass 5.6×10⁻³³ eV, and the iter-9 list.

## Iteration 10 tasks
- 034 skeptic on 031: EC-06a kill survives (non-closure holds for every solution with a growing density mode; 3.58 keV re-derived; off-tracker best case is ≤ 9 MeV; k* argument strengthened, ξ ≫ 1.6×10¹⁸; 24 Hz cited to PRL 119, 161101). Missed remnant: f(Λ)FF̃, ALP-class. 15/15 PASS (rerun).
- 035 final skeptic on draft rev 3: BLOCKER 1: EC-06d heavy-0⁺ GR-compatible branch (ξ ≳ 1, closed p) was set aside as "not O(H) torsion" but never killed; β depends on ln Y at the endpoints, and the observer endpoint may follow local density (ground vs L2 ratio -7.8×10³; script PASS on rerun). BLOCKER 2: label. EC-04 and KK-03 produce β within bounds at stated costs, which is the program's Result C definition; call is Result C on both tracks with the B-type no-go for EC-generated couplings as the main content. MAJOR 3: shortfalls quoted at ξ = 1, not at the cited bound ξ < 9×10¹⁷ (then ~10.5 and ~7.5 orders; vector-current route needs c ≈ 2×10¹⁸, not excluded by any repo bound). MAJOR 4: no "not examined" scope section (metric-affine photon sectors, Weyl-invariant EC, EC + non-minimal scalar, Hojman-Rosenbaum-Ryan gradient torsion, bosonic spin sources, D > 5, KK-05 brane-photon loophole). MINOR: 10 TeV threshold is "a choice" in data; four 029 fixes unapplied; uncovered numbers. Verified at source: CAST g_aγ·g_ae, ξ_ν, Proca galactic limit.

## Decision (overseer)
Accept 035's labelling: both tracks are Result C (effect possible only at stated costs, via fields with couplings added beyond EC), with the no-go for every coupling EC generates itself as the main content. This matches OVERSEER_PROMPT's definitions.

## Live / open
- EC-06d heavy-0⁺ branch: needs the 0⁺ mass / Compton length (SNY parameters), local vs universal endpoint, matter coupling of the 0⁺ mode (fifth-force bounds), and whether it reduces to the published ∂ ln R coupling (hep-ph/0611192). Possible kill (environment-dependent β) or a Result C entry.
- Vector-current-sourced torsion at the ξ bound (c ≈ 2×10¹⁸): any bound on c (baryon contact interactions M_P/√c ~ 10⁹ GeV)?
- ALP class with costs: EC-03', EC-04, KK-02, KK-03, EC-06a f(Λ) remnant.

## Next planned tasks (iter 11)
1. Theorist: EC-06d heavy-0⁺ branch: SNY 0⁺ mass and matter coupling, endpoint dependence, fifth-force bounds, verdict.
2. Theorist (or same): shortfalls at the cited ξ bound for EC-02 / Pauli / vector-current routes; is c ≈ 2×10¹⁸ excluded?
3. Data analyst: regenerate the uncovered numbers; cite or relabel the 10 TeV threshold.
4. Overseer: draft rev 4: Result C labels with B-type core, cost table for each surviving route, "not examined" scope section, remaining 029 fixes. Then a fresh final skeptic in iter 12.
