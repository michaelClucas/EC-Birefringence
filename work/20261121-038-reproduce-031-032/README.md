# 20261121-038 reproduce 031/032/034/035 numbers (iteration 11)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Every mechanism here is a candidate mechanism. There is no `RESULT.md`.

## Question
Can `python -m src.reproduce` regenerate the numbers skeptic 035 listed as uncovered (EC-06a 3.6 keV, 9.4 orders, ≤ 9 MeV, ξ ≥ 4.7×10¹⁶, ξ 21× larger; EC-06d -7.8×10³; KK-04 3.6-4.1; EC-03′ best-fit mass 5.6×10⁻³³ eV), and does `--checks` cover the check scripts of tasks 031, 032, 034 and 035?

## Method
All changes are in `src/`. Work scripts are executed read-only and not edited.

- `src/claims_numbers.py`: new sources in `EC06_SOURCES`: `031-mz` (`work/20261120-031/mz_numbers.py`), `034-mz` (`work/20261120-034/skeptic_mz.py`), `035-end` (`work/20261120-035/ec06d_endpoint.py`, reads `local_inputs.csv` next to itself). Each one's `_checks_ok` flag requires exit 0 and no FAIL line (035's script always exits 0, so both of its PASS lines are required).
- New rows in `src/claims_expected.csv` (claims_text = exact CLAIMS.md text):

| claim | CLAIMS text | getter | rule | regenerated |
|---|---|---|---|---|
| EC-06a | 3.6 keV | 031 `lam_sc(xi_gw, m_PDG, 1 - c_s = Ω_m)` | sig2 | 3.583 keV |
| EC-06a | 9.4 orders | 031 `log10(E_safe / Λ_SC)` | dec1 | 9.446 |
| EC-06a | ξ ≥ 4.7×10¹⁶ | 031 `xi_k` (k* at 24 Hz) | sig2 | 4.699×10¹⁶ |
| EC-06a | ≤ 9 MeV | 034 `lam_sc(xi_noGW, m_PDG, Ω_m)` | sig1 | 8.941 MeV |
| EC-06a | ξ 21× larger | 034 `J / Jalt` | sig2 | 21.19 |
| EC-06d | -7.8×10³ | 035 ground / L2 ratio `r` | sig2 | -7792 |
| KK-04 | 3.6-4.1 (lower) | 032 `kk04_shortfall`, min of bulk 2σ, brane 1σ, brane 2σ | dec1 | 3.630 |
| KK-04 | 3.6-4.1 (upper) | max of the same three | dec1 | 4.094 |

- `src/draft_expected.csv`: EC-03′ best-fit mass "5.6×10⁻³³ eV" (grid mass of β_max under the 5% cap, already computed in `ec03()`), sig2.
- `src/reproduce.py` `--checks`: `CHECKS` now also lists 025's `run_checks.py` (was an optional glob), 031 and 032 `run_checks.py`, 034 `skeptic_mz.py` and 035 `ec06d_endpoint.py`. Any `work/<date>-<NNN>-*/run_checks.py` with NNN > 035 is added automatically (`AUTO_CHECKS_AFTER = 35`). Skeptic folders without a `run_checks.py` are only picked up if named. 007, 008 and 023 `run_checks.py` stay out (their numbers are regenerated in-process; 008's is the slow EC-03′ scan). At the end of this task 036 (`theory-pgt0plus-heavy`) was in progress with no `run_checks.py` yet; it will be picked up automatically once it has one. No 037 folder existed.

Command (repo root, OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS = 1):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --log work/20261121-038-reproduce-031-032/run.log
```
The four new check scripts were also run once through `run_checks()` alone (CHECKS restricted to them): 4/4 PASS (031: 36 PASS; 032: 14; 034: 15; 035: 2; 63 s; no output files changed). The full `--checks` run was not repeated.

## Result
`run.log` (path-free): **exit 0**, 238 s.
- Validation 18/18 PASS; all `_checks_ok` flags PASS, including 031-mz, 034-mz, 035-end.
- Live CLAIMS numbers: **74/74 match** (was 66). Superseded 4/4.
- Draft numbers: 5/5 match. Two are "not in draft": the Cosmoglobe penalties `| 4.95 |` and `| 0.14 |` (see open issues). The new 5.6×10⁻³³ eV text is in the draft.
- Path scan: 0 hits.

## Verdict
Every number on 035's uncovered list is now regenerated and matches CLAIMS.md / the draft. `--checks` covers 031, 032, 034 and 035, and later run_checks.py files automatically.

## Open issues
1. The current (uncommitted) draft no longer contains the table-cell texts `| 4.95 |` and `| 0.14 |`; it now writes "4.95 (Cosmoglobe)". Values still match; `draft_text` in `src/draft_expected.csv` should be updated to the new wording once the draft settles.
2. Still not regenerated: 020's Λ scales, DMP 4.9σ, LiteBIRD 3σ, KK-05 10²¹ m_e, the 029 Pauli numbers' draft text (033 list).
3. The "9.4 orders" is measured against `E_collider_safe` = 10 TeV in `data/mz_inputs.csv`, which 035 calls a choice, not a cited bound.
4. EC-06d -7.8×10³ depends on 035's `local_inputs.csv` (local densities) and a hardcoded h = 0.674 and m_p in its script, not on data/.
