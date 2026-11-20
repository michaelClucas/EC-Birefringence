# 20261120-021 clean-reproduce: AIC k_full fix and clean-checkout test

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Every mechanism here is a candidate mechanism.

## Question
1. Skeptic 015 found that 014's ΔAIC_full gave the free-β ALP k_full = 1 and the bounded models (EC-03', KK-02) k_full = 2. Fix the count and confirm no CLAIMS.md number changes.
2. Stop condition: does `python -m src.reproduce` regenerate every claim number from data/ on a clean checkout (no uncommitted or gitignored files from the original working tree)?
3. Does `src.reproduce` still exit 0 in the main checkout after the fix?

## Method
**AIC fix.** In `work/20261120-014-data-reproduce-fit/fit_candidates.py` the free ALP now has k_full = 2: amplitude (g φ_i) and mass/timing m. A physical ALP needs m on the plateau and an initial amplitude exactly as EC-03' (φ_i, m) and KK-02 (Δϑ, 2πLM) do, so all non-null models count 2. This was chosen over dropping the column because it is the honest count and keeps the column informative. Docstring updated. k_eff (= 1 everywhere) is unchanged.

**Clean checkout.**
- `git clone . artifacts/clean_checkout` (HEAD 2453ace), then `git apply` of this task's one-file patch (the AIC fix), so the clone is HEAD plus exactly what this task changes.
- No new venv: the main repo's `.venv/Scripts/python` was called by relative path (`../../.venv/...`). The venv has no `.pth` pointing at the repo, so `import src` resolves to the clone's `src/` (cwd), and `fit_candidates.py` inserts its own repo root.
- Inside the clone: `work/20261120-010-data-likelihood/download_candl_data.py` (clones the public candl_data at the pinned commit into the clone's `artifacts/`), then `src.reproduce`, pinned, 1 thread. The clone's `artifacts/templates/` was empty, so the SPT BB templates (CAMB) and the CB calibration scale were regenerated from scratch.
- `compare_clone.py` compares the clone's regenerated templates and `cb_scale.json` with the main checkout's cached ones, and the clone's `fit_results.csv` with the committed `work/20261120-014.../fit_results.csv`.
- Main checkout: `src.reproduce` once more.

Commands (from the repo root; for the clone, from `artifacts/clean_checkout` with `../../.venv/Scripts/python`):
```
git clone . artifacts/clean_checkout
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-010-data-likelihood/download_candl_data.py
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --log <log>
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-021-clean-reproduce/compare_clone.py
```

## Result
Full output in `run.log` (path-free).

| run | exit | runtime | validation | live claim numbers | path scan |
|---|---|---|---|---|---|
| download in clone | 0 | 77 s | | | |
| `src.reproduce` in clean clone | 0 | 454 s (460 s wall) | 18 PASS | 36/36 match, 4/4 superseded reproduced | 0 paths |
| `src.reproduce` in main checkout | 0 | 184 s (189 s wall) | 18 PASS | 36/36 match | 0 paths |

- The clone took ~270 s longer: it rebuilt the CAMB templates and the CB scale (cached in the main checkout's `artifacts/templates/`), and ran while 2 other agents shared the 2 CPUs.
- Regenerated templates equal the cached ones (max relative difference 0.0); CB scale identical (0.48508).
- Clone `fit_results.csv` vs committed 014 csv: only the two ALP rows differ, only in `k_full` (1 → 2), `dAIC_kfull` (−2.00: 6.44 → 4.44 ACT DR6, 18.72 → 16.72 ACT+WMAP+Planck) and `dBIC_n6_joint_kfull` (6.65 → 4.86, 18.93 → 17.13). Penalties, Δχ², σ, B₁₀, ΔAIC_eff unchanged. ΔAIC_full now ranks the models exactly as ΔAIC_eff does (all shifted by 2), so it carries no extra information but is no longer biased against the named models.
- No CLAIMS.md number changes: the CLAIMS.md fit numbers are β_max(2.7%) and the two penalties, which do not depend on k.
- Dependence on non-committed files: none found. `git status` in the clone after the run shows only the applied patch and the gitignored `artifacts/`. The only gitignored inputs are the candl_data download (fetched by the committed script at a pinned commit) and the template/CB-scale caches (regenerated, identical). Nothing in the clone read from the original checkout's `artifacts/` or work tree, and no committed file in the clone was rewritten by the run.

## Verdict
Stop-condition item met: on a clean clone of HEAD (+ the AIC patch), the download script followed by `python -m src.reproduce` regenerates all 36 live CLAIMS.md numbers from data/ and exits 0. No `src/` fix was needed. The ALP now counts k_full = 2 like the other models; no CLAIMS.md number changes.

## Open issues
- `work/20261120-014-data-reproduce-fit/fit_results.csv` and its README still show the old ALP ΔAIC_full (6.44 / 18.72) and say k_full = 1. Not regenerated here (only `fit_candidates.py` was allowed); `src.reproduce` writes the corrected csv to `artifacts/reproduce/`. Overseer may regenerate the committed csv or add a note to 014.
- The clean test needs network access to GitHub (candl_data) and reuses the main `.venv`; a fresh venv from `requirements.txt` was not tested.
- candl prints "Failed to import `candl_data`" in both runs (the package is not pip-installed; the repo points candl at the YAML path directly). Harmless: validation 18 PASS including the candl cross-check.
- The clone was run with the main checkout's `tools/pin.py` (`../../tools/pin.py`); the file is identical in the clone.
