# 20261120-027 reproduce EC-06 numbers, `--checks` option

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. EC-06 is a candidate mechanism. There is no `RESULT.md`.

## Question
1. Do the EC-06, EC-06a-d numbers in CLAIMS.md (predicted β and needed-parameter columns) regenerate from the 022/023/024 scripts, and does `python -m src.reproduce` cover them?
2. Skeptic 020 MINOR 12: can one command rerun the symbolic checks behind the verdicts that carry no number?
3. Skeptic 026 corrections (EC-06a GW side, EC-06c reionization share): source them from the script that prints them.

## Method
- `src/claims_numbers.py`: new table `EC06`, one line per number (`"claim.key": (source, getter)`), and `EC06_SOURCES`, each run once per reproduce:
  - `022`: `ec06_numbers.py` executed in-process with its `__file__` pointed at `artifacts/reproduce/`, so `ec06_backgrounds.csv` is written there and not into work/.
  - `023`: `ec06_common.beta_unit(z*, n = 0)` gives ξε. `ec06_aniso.py` is executed the same way as 022 (`acb_estimate.csv` goes to artifacts/), and A_CB is the max over its models A and B.
  - `024`: `kaf_lowz.py` (read-only mode) executed with cwd = repo root. ξε < 0.14 is the L-AFG10-UV row it recomputes.
  - `026`: `skep_beta.py` (read-only, cwd = repo root). It gives the c_T > 1 side of GW170817 (ε_S = 2γ/3, ξ = ξε/ε_S) and the visibility-averaged reionization share.
  - New helper `exec_script()`. Each source adds a `_EC-06.<src>_checks_ok` flag that requires no FAIL line or nonzero exit.
- `src/claims_expected.csv`: 14 EC-06 rows. claims_text is exact CLAIMS.md text, with the 026 wording for the corrected rows.
- `src/reproduce.py --checks` (default run unchanged) runs these scripts:
  - `run_checks.py` of 003, 004, 009, 012, 016 and 022, plus `019/adjudicate_edm.py`.
  - Any `work/20261120-025-*/run_checks.py` that exists at run time.
  - Each script runs in its own `tools/pin.py` process with 1 thread and cwd = repo root. Full output goes to `artifacts/reproduce/checks_<task>.log` (path-scrubbed).
  - A script fails on a nonzero exit, a `FAIL`/`[FAIL]` line, a nonzero "N FAIL" count, or "SOME FAIL".
  - Output files (.log/.csv/.png/.json/.txt) that the scripts write into their own work folder are snapshotted and restored afterwards, so work/ stays unchanged. 009, 012 and 016 rewrite their `run.log`. Those three files were restored.
- Skipped, not verbatim in CLAIMS.md: ALP Δχ² = 20.7, "Δχ > 0.71 M_P", "≥ 99.6%", "factor 273".

Commands (repo root):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --log work/20261120-027-reproduce-ec06/run.log
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --checks --log <file>
```

## Result
`run.log` is the default run with the 026-corrected rows. `run_checks.log` is the `--checks` run made before the 026 corrections, when EC-06a/EC-06c still used the 022 values. Both logs are path-free.

| claim | CLAIMS text | regenerated | source |
|---|---|---|---|
| EC-06 | ξε ≈ 1.07×10⁻³ | 1.0729×10⁻³ | 023 |
| EC-06 | ξε < 0.14 | 0.136 | 024 |
| EC-06 | A_CB ≤ 5×10⁻¹¹ | 4.93×10⁻¹¹ | 023 aniso |
| EC-06a | ε_S ≤ 7.5×10⁻⁸ (026; CLAIMS still 1.5×10⁻⁷) | 7.48×10⁻⁸ | 026 |
| EC-06a | ξ ≥ 1.4×10⁴ (026; CLAIMS still 7×10³) | 1.43×10⁴ | 026 |
| EC-06c | ε_T = 3 | 3 | 022 |
| EC-06c | 23-31% (026; CLAIMS still 31%) | 23.0 (visibility) / 30.9 (z_re) | 026 / 022 |
| EC-06c | ξ = 3.6×10⁻⁴ | 3.576×10⁻⁴ | 022 |
| EC-06c | ξ ≥ 0.022 | 0.0219 | 022 |
| EC-06d | T_t = -(3/2) d ln(R+6μ/b)/dt | -3/2 (sympy) | 022 |
| EC-06d | Φ/H ≈ 0.8-1.4 | 0.8, 1.4 | 022 (transcribed SNY Table 1, not computed) |

- Default run: 18 validation PASS. 45/49 live numbers match on text and value. 4/4 superseded numbers are reproduced. 260 s, of which EC-06 takes 43 s.
- The 4 rows not matching are the 026-corrected texts (ε_S ≤ 7.5×10⁻⁸, ξ ≥ 1.4×10⁴, and the 2 values of "23-31%"). Their values match, but the text is not yet in CLAIMS.md, so the exit code is 1. It becomes 0 once the overseer applies 026's wording.
- Before the 026 corrections (`run_checks.log`), the run gave 48/48 match and exit 0, including the old 1.5×10⁻⁷, 7×10³ and 31%.
- `--checks`: 7/7 scripts PASS, 697 s for the checks:

  | task | time |
  |---|---|
  | 003 | 130 s |
  | 004 | 41 s |
  | 009 | 166 s |
  | 012 | 116 s |
  | 016 | 114 s |
  | 019 | 44 s |
  | 022 | 86 s |

  The whole run with `--checks` took 963 s. The default part stays under 5 min.
- No file in another work/ folder changed. `git status` shows only src/ and this folder from this task.

## Verdict
Every EC-06 number in CLAIMS.md regenerates from the work/ scripts. After the CLAIMS.md update for 026, `src.reproduce` should exit 0 with 49/49. The Φ/H ≈ 0.8-1.4 range is a literature transcription and is only checked against 022's table. The derivation checks behind the no-number verdicts are now rerunnable with `--checks`.

## Open issues
- CLAIMS.md needs 026's texts (ε_S ≤ 7.5×10⁻⁸, ξ ≥ 1.4×10⁴, 23-31%) for exit 0. If the wording differs, edit claims_text in `src/claims_expected.csv`.
- 025 had no `run_checks.py` when this ran (only `check_ft_lorentz.py`, `check_stueckelberg.py`). `--checks` will include it automatically once one exists. If 025 revises EC-06c/EC-06d numbers, change one getter line in `EC06` and one row in `claims_expected.csv`.
- `--checks` with all tasks is about 12 min on the shared CPUs, so it stays opt-in.
- The `--checks` restore of output files could overwrite a concurrent edit by another agent to a .log/.csv in those folders. This is unlikely, since those tasks are closed.
- `data/ec06_cosmo_inputs.csv` shows an uncommitted change from another agent. This task does not write data/.
