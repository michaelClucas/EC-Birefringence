# 20261120-017 reproduce: iter-3 claim numbers, fit summary, path-free logs

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Every mechanism here is a candidate mechanism.

## Question
1. Does `python -m src.reproduce` regenerate every number now in CLAIMS.md, including the iter-3 ones (EC-03 photon decoupling, EC-03' 2.7%/6.9% caps, fit penalties, KK-05 c bounds)?
2. Can it print a compact fit summary (model, β_max, penalty vs free ALP, Δχ² vs null) for ACT DR6 and ACT+WMAP+Planck?
3. Can its output be made immune to absolute local paths (iter 3: a candl SyntaxWarning put one into a committed run.log)?

## Method
Changes in `src/` only. No other work/ folder's scripts edited; they are imported and called.

- `src/claims_expected.csv`: 11 new rows. `claims_text` is the exact substring of the CLAIMS.md row; `source` names the script.

| claim | key | CLAIMS.md text | expected | rule | source |
|---|---|---|---|---|---|
| EC-03 | g_ratio_min, g_ratio_max | 10⁻⁷⁸-10⁻⁷⁰ | 1e-78, 1e-70 | exp | 012 `check_formfactor` F5: min/max of the printed g_eff/g_anom values |
| EC-03 | g_ratio_coeff | (m_φ/m_e)²/96 | 96 | sig 2 | 012 `check_formfactor.ratio()` |
| EC-03' | beta_max_cap27_deg | ≤ 0.100° at the cited 2.7% cap | 0.100 | dec 3 | 014 `fit_candidates` |
| EC-03' | cap_P18BOSS_pct, cap_P18_pct | 2.7% cap; caps 6.9% at 10⁻³² eV | 2.7, 6.9 | dec 1 | 014 `fit_candidates.ultralight_caps` (data/ultralight_fraction_bounds.csv) |
| EC-03' | penalty_cap27_AWP | Δχ² penalty vs free ALP 8.0 (ACT+WMAP+Planck, 2.7% cap) | 8.0 | dec 1 | 014 `fit_candidates.main` |
| KK-02 | penalty_k1_AWP | Δχ² penalty vs free ALP 0.9 at k = 1 | 0.9 | dec 1 | 014 `fit_candidates.main` |
| KK-05 | c_min_MICROSCOPE, c_min_clocks, c_max_EDM | MICROSCOPE needs c ≥ 340; clocks c ≳ 5×10⁴; EDM estimate c ≲ 140 | 340, 5e4, 140 | sig 2/1/2 | 012 `bounds_uholst` B2, B3, B4 |

- New comparison rule `exp`: same decimal exponent in scientific notation. "10⁻⁷⁰" stands for 3.99×10⁻⁷⁰, which the existing `oom` rule (rounded log10) would read as 10⁻⁶⁹.
- `src/claims_numbers.py`:
  - `ec03_photon()`, `kk05()`: import the 012 top-level scripts with stdout captured; flag `[FAIL]` in their output.
  - `KK05` dict: one line per KK-05 number, `key -> (work folder, module, getter)`. If task 016 revises KK-05, change that line (if the source moves) and the matching CSV row (text and expected value).
  - `fit_values()`: reads the fit numbers from 014's `fit_results.csv` rows.
- `src/reproduce.py`, new order: validation → regenerate → fits (014's `HERE` patched so `fit_results.csv` goes to `artifacts/reproduce/`, not work/014) → claim table → fit summary → path scan.
- Paths:
  - `warnings.filterwarnings("ignore", SyntaxWarning)` before any import, and `-W ignore::SyntaxWarning` in the command.
  - stdout/stderr are wrapped: every write is scrubbed (repo root → `<repo>`, home → `<home>`, any other `X:\...` or `/Users/...` → `<abs-path>`). `--log FILE` tees the scrubbed output.
  - Final step rescans the log and the text outputs in `artifacts/reproduce/` and rewrites them if anything is found. Tested on sample strings (backslash, forward slash, lower case, `/c/...` form, home, foreign drive).

Run (from the repo root, 1 CPU, 1 thread; ~2.5 min):
```
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/Scripts/python tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --log work/20261120-017-reproduce-iter3/run.log
```

## Result
`run.log`: exit 0, 149 s.
- Validation: 18 PASS.
- Claim numbers: **36/36 live match**, 0 mismatch, 0 not regenerated. 4 superseded EC-03 QCD-mass numbers still reproduced. All script check flags PASS (007, 008, 009, 012 form factor, 012 bounds).
- New numbers: g_eff/g_anom 1.44×10⁻⁷⁸ to 3.99×10⁻⁷⁰, coefficient 96; EC-03' β_max(2.7%) = 0.0998°; caps 2.675% and 6.942%; penalties 8.02 (EC-03' 2.7%) and 0.90 (KK-02 k = 1); KK-05 c ≥ 340.2, c ≥ 5.13×10⁴, c ≲ 144.
- Path scan: 0 absolute paths on write, 0 in artifacts, 0 in the log.

Fit summary (from the log):

| model | ACT DR6 β_max | penalty | Δχ² vs null | ACT+WMAP+Planck penalty | Δχ² vs null |
|---|---|---|---|---|---|
| ALP (free β) | free | 0 | 8.44 | 0 | 20.72 |
| EC-03' 2.7% (P18+BOSS) | 0.100 | 2.42 | 6.02 | 8.02 | 12.70 |
| EC-03' 5% (008 assumed) | 0.136 | 1.13 | 7.31 | 4.84 | 15.88 |
| EC-03' 6.9% (P18) | 0.161 | 0.54 | 7.90 | 3.17 | 17.55 |
| EC-03' 12% (excluded by P18) | 0.211 | 0.00 | 8.44 | 0.82 | 19.89 |
| KK-02 k = 1 | 0.209 | 0.01 | 8.43 | 0.90 | 19.82 |
| KK-02 k = 2 | 0.418 | 0 | 8.44 | 0 | 20.72 |

## Verdict
Every number in CLAIMS.md is regenerated and matches. No CLAIMS.md change needed. The EC-04 "~4×10¹²⁰" mismatch reported by 014 is resolved (CLAIMS.md now prints it at 1 significant figure).

## Open issues
- The EC-03 range is the min/max of the three values 012 prints (m_φ = 6×10⁻³³ eV, 10⁻²⁸ eV, q ~ H(z*)). If 012's mass list changes, the range follows it.
- The candl SyntaxWarning only fires when candl's `.pyc` is rebuilt, so this run did not trigger it. The filter and the scrubber were tested separately, not against a live candl recompile.
- KK-05 c ≲ 140 is printed as 144 by 012 (sig-2 match). If 016 firms up the EDM bound, update the `KK05` line in `src/claims_numbers.py` and the CSV row.
- The `act_wmap_planck` combination and BIC caveats of 014 still apply.
