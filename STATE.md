# State

Iteration: 12 (completed 2026-11-21)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 12: `repin.ps1` re-pinned 2 python.exe and several agent bash shells mid-iteration (shells match the repo path; harmless).
- Path leak in commit 87caa4e (iter 3) remains in history; path-scan gate before every commit since.
- `python -m src.reproduce`: 18 validation + 79/79 CLAIMS numbers + all draft numbers, exit 0 after the iter-12 audit (a "Failed to import candl_data" warning is printed and is harmless). `--checks`: 14/14 scripts PASS (040's run).

## Iteration 12 tasks
- 039 skeptic on 036/037: EC-06d heavy branch still dead at every mass (sign kill holds for ~5 km-20 pc since the Planck/WMAP source endpoint is inside dust grains; size kill with ground rock covers λ_C ≲ 2×10³ km; universal piece ~26 orders, not 41.4). 037 survives; the 1.2-order weak corner closes: Stückelberg cutoff 2.3 MeV at the needed ξ, 3.7 orders below BaBar's 10.58 GeV; Pauli route then 11.0 and EC-02 13.9 orders short (collider energies recalled, not opened). 22/22 PASS.
- 040 final skeptic on draft rev 4: verdict "both tracks survive with listed fixes, Result C with a B-type core"; no published paper states the overall result; six bounds verified at source; honesty clean. BLOCKER (procedural): reproduce must cover 036/037 numbers on a committed checkout (now done except as noted below). MAJOR: quantify remaining Result C costs (EC-06a remnant g range 10⁻²⁰-10 GeV⁻¹ from 034 plus ghost question; KK-03 λ/κ₄ ≈ 1.5-8.7×10⁶ and the pre-BBN kick; KK-02 hilltop-start tuning); KK-05 ~21 orders and 020's higher-dim cutoffs not in reproduce. MINOR: headline "only"/"many orders" → "among routes examined"; new measurement arXiv:2608.06480 (β = 0.277° ± 0.057°, 4.8σ) not in data/ (KK-02 penalty ≈1.4, EC-03' ≈9.7 against it); EC-04 upper limit should use Chandra (ξ < ~3×10⁶); Cassini convention mismatch KK-04 vs EC-06d; 10 TeV "choice" note; 2212.13871 missing from literature.
- 041 data: 036/037 numbers into reproduce.
- Overseer: draft updated with 039's corrections (weak corner closed, EC-06d endpoint wording, +13/+15 orders).

## Settled
- Both tracks: Result C (mixed) with a B-type core. Every photon coupling EC torsion generates by itself fails (4D and KK); every unsourced-torsion route examined is dead; β is reachable only by ALP-like fields with couplings added beyond EC, at the costs in the draft's cost table. Two final skeptics' verdicts (035 relabel, 040 survives-with-fixes) agree on this.

## Remaining before RESULT.md (iter 13)
1. Data: add 2608.06480 to data/ (verified at source) and refit penalties; EC-04 Chandra bound; KK-05 ~21 orders and 020 cutoffs into reproduce; cite BaBar 10.58 GeV search.
2. Theory: quantify KK-02 hilltop tuning, KK-03 λ/κ₄ and pre-BBN kick, EC-06a remnant range and ghost question (short).
3. Literature: 2212.13871; Cassini convention alignment.
4. Overseer: apply wording fixes; ask skeptic 040 to confirm its fixes are applied (it found them; confirmation only); then write RESULT.md, commit, push, notify, stop.
