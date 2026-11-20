# 20261120-030 reproduce the task 025 numbers (EC-06c, EC-06d)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. EC-06c and EC-06d are candidate mechanisms. There is no `RESULT.md`.

## Question
CLAIMS.md's EC-06c and EC-06d rows gained numbers from task 025 that `python -m src.reproduce` did not cover. Do they regenerate from 025's scripts, and does the full reproduce (default, and `--checks` on a clean clone of HEAD) pass?

## Method
- `src/claims_numbers.py`: three new sources in `EC06_SOURCES`, each a 025 script executed in-process with `exec_script()` (cwd = repo root, `__file__` pointed at `artifacts/reproduce/`, so `sny_background.csv` is written there and not into work/):
  - `025-num`: `ec06c_numbers.py`
  - `025-sny`: `sny_background.py`
  - `025-lim`: `check_sny_limit.py`

  The three share one `common.py` PASS/FAIL counter, so each `_EC-06.025-*_checks_ok` flag uses that script's own `FAIL:` lines and exit code.
- Eight new lines in `EC06` and rows in `src/claims_expected.csv` (claims_text = exact CLAIMS.md text):

| claim | CLAIMS text | getter | rule |
|---|---|---|---|
| EC-06c | 7×10⁻²⁸ ≲ m_γ ≲ 10⁻²⁶ eV (lower) | `m_low` = max(conversion A < 0.1, CFJ stability) | sig1 |
| EC-06c | same (upper) | `m_up_proca` (data/ec06c_costs_inputs.csv, hep-ph/0306245) | sig1 |
| EC-06c | ~20-60 GeV (lower) | (m_low M_P²/ξ)^(1/3) | oom |
| EC-06c | ~20-60 GeV (upper) | (m_up_proca M_P²/ξ)^(1/3) | oom |
| EC-06d | ξ ≈ 4.3×10⁻⁴ | mean of `xi_needed` at z* over SNY cases III-V | sig2 |
| EC-06d | H(z*) = 33-66 × ΛCDM (lower, upper) | min, max of `H_over_LCDM` at z* over III-V | sig2 |
| EC-06d | ξ ≳ 1 | `xi_per(1e-2)` (μ = 1, \|δH²/H²(z*)\| ≤ 10⁻²) | oom |

- `sny_background` raises if fewer than 3 cases reach z*.
- Timing label in `regenerate()` changed to "022-026 (EC-06)".

Commands (repo root; clean clone at `artifacts/clean_checkout`, updated with `git pull` to b95b350 plus the uncommitted src/ edits copied in):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --log work/20261120-030-reproduce-025/run.log
(in artifacts/clean_checkout) ../../.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py ../../.venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --checks --log ../reproduce/030_clean_checks.log
```
Both with OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS = 1.

## Result
Default run in the main checkout (`run.log`, path-free): **exit 0**, 18 validation PASS, all 16 `_checks_ok` flags PASS, **57/57 live claim numbers match**, 4/4 superseded reproduced. 390 s wall (383 s inside reproduce), EC-06 block 61 s (was 43 s).

| number | regenerated |
|---|---|
| m_γ lower | 6.64×10⁻²⁸ eV |
| m_γ upper | 1.0×10⁻²⁶ eV (data input) |
| cutoff lower | 22.3 GeV |
| cutoff upper | 54.9 GeV |
| ξ (SNY, z*) | 4.335×10⁻⁴ (cases 4.17, 4.48, 4.35) |
| H(z*)/H_ΛCDM | 33.2, 65.8 |
| ξ (GR-like) | 1.12 |

The first default run used sig1 for the cutoff and failed on the upper end: 54.9 GeV rounds to 5×10¹, not the printed 6×10¹ (025 printed 5.5×10¹⁰ eV and the README/CLAIMS rounded up). Since CLAIMS prints it with "~" (a power-counting estimate), both cutoff rows now use `oom`, as "~10¹⁸ GeV" in KK-02 does. This is flagged below.

`--checks` on the clean clone (`run_checks.log`, path-free): **exit 0**, 1696 s wall (1685 s inside reproduce). Default part 57/57 match. Checks 8/8 scripts PASS, 1326 s, with 025's `run_checks.py` picked up automatically:

| task | PASS lines | time |
|---|---|---|
| 003 | 47 | 225 s |
| 004 | 13 | 79 s |
| 009 | 33 | 322 s |
| 012 | 37 | 197 s |
| 016 | 28 | 176 s |
| 019 | 9 | 24 s |
| 022 | 54 | 75 s |
| 025 | 70 | 228 s |

The clone's `git status` after the run shows only the copied src/ files. Times are slower than 027's because the two pinned CPUs were shared with other agents.

No file in another work/ folder changed.

## Verdict
Every number in CLAIMS.md's EC-06c and EC-06d rows now regenerates from the 025 scripts inside `src.reproduce`, and both the default run and `--checks` on a clean clone of HEAD exit 0. The cutoff's upper end matches only at order of magnitude: its value is 55 GeV.

## Open issues
- CLAIMS.md "~20-60 GeV": the computed range is 22-55 GeV. Suggest "~20-55 GeV" (then the upper row can go back to sig rules with expected 55, digits 2).
- "ξ ≈ 4.3×10⁻⁴" is checked as the mean of the three SNY cases; single cases span 4.2-4.5×10⁻⁴ (025 README). Case V alone (4.35) is borderline at 2 significant figures.
- "10⁻²⁶ eV" is a data input (Proca galactic limit), so its row checks that the data file and CLAIMS agree, not a computation.
- "ξ ≳ 1" depends on 025's stated choices (μ = 1, δ ≤ 1%); the row reproduces the number, not the choice.
- The clean clone's stash holds the old task-021 patch to `fit_candidates.py`; HEAD already contains it.
