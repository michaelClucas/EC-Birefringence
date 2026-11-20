# State

Iteration: 8 (completed 2026-11-20)

The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with stated systematics caveats (miscalibration angle set by an optics-model prior, dust EB neglected). The signal may not be cosmological. Nothing below is an explanation; every mechanism is a candidate mechanism. No `RESULT.md`.

## Environment
- `.venv`: Python 3.12.14 (uv), packages pinned in `requirements.txt` (candl 2.2.0).
- CPU: logical CPUs 12-13 only, via `tools/pin.py`. iter 8: `repin.ps1` re-pinned 2 python.exe mid-iteration (launcher-stub pattern).
- `artifacts/` (gitignored): candl_data, reproduce outputs, `clean_checkout/` (safe to delete).
- Path leak in commit 87caa4e (iter 3) remains in history; path-scan gate before every commit since.
- `python -m src.reproduce` (~3-6 min): 18 validation + 57/57 CLAIMS numbers, exit 0 after iter-8 audit. `--checks` (~28 min on the clean clone): 8/8 symbolic check scripts PASS (030). Not covered: KK-04 Cassini shortfall (~4.3 orders), Cosmoglobe penalties (4.95, 0.14), 028's Λ_η numbers, 029's Pauli-bound numbers.

## Iteration 8 tasks
- 028 skeptic on 025: ζ solution and dT ≠ 0 survive; gravitational CS not a repair (survives); Stückelberg window BROKEN (white dwarfs push the lower edge to ~10⁻²⁴ eV; WD/NS inputs recalled, in work/028/inputs_028.csv, not yet cited); new strong-coupling scale from the f(T) Lorentz mode η: Λ_η ≈ 0.2-0.7 keV (7 MeV at the PDG limit), far below where photon physics is standard. EC-06c dead. EC-06d survives as dead. 15/15 PASS (rerun).
- 029 final skeptic on draft rev 2: BLOCKER: 4D verdict "Result B for EC proper" contradicts EC-06a (MZ, an EC extension) which is only disfavoured by naturalness (ξ ≥ 1.4×10⁴); its off-FRW U(1) breaking and the off-tracker (decaying-mode) loophole were never examined. MAJOR: "Result B" headlines must say ALP-class routes with hand-added couplings reproduce β; missing shortfalls (KK-04 vs Cassini ~4.3 orders; EC-06d needs a cited H(z*) bound); reproduce gaps. Checked and holding: EC-proper no-go vs any fermionic source (≥ 25.5 orders short; the draft's "≥ 28" is too strong in general); vector-current-sourced torsion needs c ≈ 10³⁶; six bounds verified at source (1312.3225, 2212.11841, 2002.11761, 1710.05834, 2209.15487, Cassini via Will); novelty list consistent; 2610.02667 missing from literature/. Pauli-bound script PASS (rerun).
- 030 data: 025 numbers in reproduce; `--checks` on clean clone 8/8.

## Settled
- 4D Einstein-Cartan proper: no β from minimal EC or from any fermionic torsion source (≥ 25.5 orders short), or torsion pseudoscalars with ∂φ·J₅.
- O(H) torsion: TEGR undefined; f(T) (EC-06c) dead by strong coupling; SNY (EC-06d) dead as O(H) torsion.
- KK track: KK-01, -04, -05 dead; KK-02/03 are ALPs with hand-added couplings, torsion plays no role.

## Live / open
- EC-06a Magueijo-Złośnik (EC extension: Λ as field + Euler/Pontryagin): disfavoured only by naturalness (ξ ≥ 1.4×10⁴ on the tracker). Open: off-FRW U(1) breaking and strong-coupling scale with the 025/028 machinery (029 estimates ~0.1 GeV if the structure carries over); off-tracker solutions where GW170817 may not bound ε at z*.
- ALP class (KK-02, KK-03, EC-03', EC-04).

## Direction toward a stop
KK track: Result B for torsion, with the explicit statement that ALP-class routes with hand-added couplings (KK-02/03) reproduce β. 4D track: depends on EC-06a. If the 025/028 machinery kills MZ, 4D is Result B with the same ALP statement; if not, 4D is Result C with MZ's costs quantified (ξ ≥ 1.4×10⁴, extension terms, U(1) status).

## Next planned tasks (iter 9)
1. Theorist: EC-06a MZ: off-FRW closure of the CFJ vector, strong-coupling scale (025/028 method), and the off-tracker loophole for the GW bound. Verdict dead / Result C with costs.
2. Literature + bounds: cite the WD/NS inputs from 028 (move to data/), a CMB bound on H(z*) for EC-06d, 2610.02667 into literature/, and the KK-04 Cassini shortfall in orders.
3. Data analyst: reproduce coverage for KK-04 shortfall, Cosmoglobe penalties, Λ_η numbers, Pauli-bound numbers; commit-clean run.
4. Overseer: draft rev 3 (headline wording per 029 MAJOR 2; "≥ 25.5 orders" for general fermionic sources; shortfalls), then a fresh final skeptic in iter 10.
