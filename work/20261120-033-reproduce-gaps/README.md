# 20261120-033 reproduce gaps (iteration 9)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Every mechanism here is a candidate mechanism. There is no `RESULT.md`.

## Question
Can `python -m src.reproduce` regenerate the numbers skeptic 029 listed as uncovered and the iteration-8/9 additions to CLAIMS.md (028's Λ_η and white-dwarf edge, EC-06d z = 10 ratio and ξ ≥ 0.1, KK-04 Cassini shortfall, 032's radiation-corrected H(z*) ratios and shortfall), plus the draft-only numbers (Cosmoglobe penalties, 029 Pauli bound)?

## Method
All changes are in `src/`. Work scripts are executed in-process (read-only), never edited.

- `src/claims_numbers.py`: new sources in `EC06_SOURCES`:
  - `028-win`: `work/20261120-028/window_cutoff.py` (runs with its own folder as `__file__`, since it reads `inputs_028.csv` next to itself; writes nothing).
  - `029-pauli`: `work/20261120-029/pauli_bound.py`.
  - `032-kk04`: `work/20261120-032/kk04_shortfall.py`.
  - `032-hz`: `work/20261120-032/hz_bounds.py` (reads 025's committed `sny_background.csv`).

  The source runner (`_run_table`) gives each namespace `_src(name)`, so a getter can read another source: Λ_η at the window's lower edge uses 025's `m_low`, not 028's literal 7×10⁻²⁸.
- New CLAIMS.md rows in `src/claims_expected.csv` (claims_text = exact CLAIMS.md text):

| claim | CLAIMS text | getter | rule | regenerated |
|---|---|---|---|---|
| EC-06c | 0.2-0.7 keV (lower) | `lam_eta(m_low(025), α = 1)` | dec1 | 0.192 keV |
| EC-06c | 0.2-0.7 keV (upper) | `lam_eta(m_proca, 1)` | dec1 | 0.747 keV |
| EC-06c | 7 MeV | `lam_eta(m_pdg, 1)` | sig1 | 7.47 MeV |
| EC-06c | ~10⁻²⁴ eV | W1, 0.6 M_⊙ WD, A < 0.1 | oom | 2.55×10⁻²⁴ eV |
| EC-06d | H/H_ΛCDM is already 5.9 at z = 10 | max over III-V at z = 10 (025 `sny_background`) | sig2 | 5.876 |
| EC-06d | ξ ≥ 0.1 | 025 `xi_per(0.1)` | sig1 | 0.1125 |
| EC-06d | H(z*) = 29-57 × ΛCDM (re-sourced from 032; 025's 33-66 omitted radiation) | min/max over III-V at z* | sig2 | 28.89, 57.28 |
| EC-06d | 3.8-4.5 orders | min/max orders in δH/H at z*, over cases and both 95% bounds | dec1 | 3.757, 4.494 |
| KK-04 | ~4.3 orders | log10(0.5 / 2.3×10⁻⁵), bulk matter | dec1 | 4.337 |

- Draft-only numbers: new `src/draft_expected.csv` (column `draft_text`), checked by `compare_draft()` against `work/20261120-018-result-draft/RESULT_DRAFT.md`. A wrong or missing value fails. Text not in the draft is reported as "not in draft" and does not fail.

| claim | draft text | source | regenerated | in draft |
|---|---|---|---|---|
| EC-02 | ≥ 25.5 orders | 029 Pauli bound, g = 100, ρ_f = ρ_ΛCDM | 25.49 | no |
| EC-02 | c ≈ 10³⁶ | 029 `c_need` (vector-current source, ξ = 1) | 1.69×10³⁶ | no |
| EC-03' | \| 4.95 \| | 014 fits on `act_cosmoglobe` | 4.946 | yes |
| KK-02 | \| 0.14 \| | same | 0.138 | yes |

- `src/reproduce.py`:
  - `fits()` sets `fit.DATASETS` to all three combinations, so 014's fits also run on ACT + Cosmoglobe WMAP+Planck. The fit summary prints that block.
  - New step 4b: the draft table.
  - Bug fix: an exception inside the steps used to print "reproduce: OK", and its traceback went out unscrubbed (an absolute path appeared on stdout in the first run). It is now caught, printed through the scrubbed stderr, and counted as a failure (exit 1).
- `src/claims_expected.csv`: 030's `cutoff_max_GeV` row had an unquoted comma in `source`, which spilled into `status`. Now quoted (no effect on the comparison).

Command (repo root, OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS = 1):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --log work/20261120-033-reproduce-gaps/run.log
```

## Result
`run.log` (path-free): **exit 0**, 150 s.
- Validation: 18/18 PASS.
- `_checks_ok` flags: 19/19 PASS, including the new 028-win, 029-pauli, 032-kk04 and 032-hz.
- Live CLAIMS numbers: **66/66 match** (was 57). Superseded: 4/4.
- Draft numbers: 4/4 match. The two Pauli numbers are "not in draft" yet.
- Path scan: 0 hits.

ACT + Cosmoglobe penalties vs the free ALP: EC-03' 2.7% cap 4.95, KK-02 k = 1 0.14. These agree with skeptic 020's S4 hand values.

## Verdict
Every number on the 029 list that the task named, plus the iteration-8/9 CLAIMS additions, is now regenerated and matches. Reproduce exits 0 on the current tree. This includes CLAIMS.md and data/ edits by the overseer and task 032 that are not yet committed.

## Open issues
1. **Not regenerated** (029's list; not in this task's scope):
   - higher-dimension operator scales Λ ≲ 2×10⁻⁴³ eV and 4×10⁻¹⁵ eV (020 `skeptic_checks.py` S3);
   - DMP β ∝ ν² at ~4.9σ (literature);
   - LiteBIRD 0.05° at ~3σ (023);
   - KK-05 m ≳ eM_P/√2 ~ 10²¹ m_e.
2. **WD edge inputs.** The white-dwarf edge (~10⁻²⁴ eV) uses 028's `inputs_028.csv` (recalled values). 032 added `data/compact_star_inputs.csv` with opened sources (Sirius B radius differs by 1%, typical WD unchanged), but 028's script still reads its local file. Its order of magnitude is unaffected.
3. **KK-04 brane matter.** The shortfall row uses bulk matter (γ - 1 = -1/2). 032 also prints brane matter (4.09 orders) and 2σ values (3.87, 3.63). Those are not in CLAIMS.md.
4. **Dependency on 032.** The EC-06d z* ratios and the shortfall depend on 032's uncommitted script and `data/expansion_history_bounds.csv`. They must be committed together with src/. 025's no-radiation ratios (33-66) are no longer checked.
5. **"≥ 25.5 orders".** The regenerated value is 25.49. It rounds to 25.5 at one decimal, but "≥ 25.5" is strictly 0.01 too strong. "≥ 25.4" or "≈ 25.5" would be exact.
