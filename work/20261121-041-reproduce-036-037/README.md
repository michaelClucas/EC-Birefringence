# 20261121-041 reproduce 036/037 numbers (iteration 12)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Every mechanism here is a candidate mechanism. There is no `RESULT.md`.

## Question
Can `python -m src.reproduce` regenerate the new EC-06d numbers of task 036 in CLAIMS.md ("3.9 orders", "+0.63 ± 0.28", "23.7 orders", "α ≥ 1/3") and the 036/037 numbers in the draft's 4D kill table and caveat paragraph, and do they match the printed text?

## Method
All changes are in `src/`. The 036 and 037 scripts are executed in-process, read-only, not edited.

- `src/claims_numbers.py`: new sources in `EC06_SOURCES`:
  - `036-mode`: `check_scalar_mode.py` (sympy SNY 0⁺ reduction).
  - `036-end`: `endpoint_bounds.py` (Cassini, ACT/MK endpoint ratio, HB size).
  - `037-sf`: `shortfalls.py` (β_max at the largest allowed ξ per completion class).
  - Check flags: 036 needs exit 0, no `FAIL:` and `== N PASS, 0 FAIL ==`. 037 always exits 0, so it needs no `FAIL`, exactly 12 `PASS:` lines and `ALL PASS`.
- `src/claims_expected.csv` (EC-06d, exact CLAIMS text):

| CLAIMS text | getter | rule | regenerated |
|---|---|---|---|
| 3.9 orders | 036 `log10(0.5/lim_g)`, lim_g = \|γ-1\| + 2σ Cassini | dec1 | 3.873 |
| +0.63 ± 0.28 (value) | 036 `r_meas` = β_ACT/β_EK22 | dec2 | 0.6287 |
| +0.63 ± 0.28 (error) | 036 `r_err` | dec2 | 0.2769 |
| 23.7 orders | 036 `log10(beta_act/b_obs)` at λ = 100 μm, HB bound | dec1 | 23.72 |
| α ≥ 1/3 | 036 sympy `alpha` = (a1+1)/(3a1), limit a1 → ∞ | sig4 | 0.3333 |

- `src/draft_expected.csv` (exact draft text; 037 unless noted). D = dual-axion class (CAST product, CAST × RG, Chandra × RG). S-a = Stückelberg, Λ_s > T_BBN, Retino weak-end or PDG m_γ.

| draft text | getter | rule | regenerated |
|---|---|---|---|
| 10.5-13.3 orders | EC-02 ν, min / max over D | dec1 | 10.49 / 13.29 |
| ξ < 9.4×10¹⁷ / 1.5×10¹⁵ | `xi_D_CAST` / `xi_D_ChRG` | sig2 | 9.399×10¹⁷ / 1.485×10¹⁵ |
| 4.2-9.1 orders | EC-02 ν, S-a Retino / PDG | dec1 | 4.187 / 9.084 |
| 7.5-10.3 | Pauli, min / max over D | dec1 | 7.52 / 10.32 |
| 1.2-6.1 orders | Pauli, S-a Retino / PDG | dec1 | 1.22 / 6.118 |
| ≥ 11.1 | vector J_B, min over D contact (ALPS II) | dec1 | 11.06 |
| 12.0-16.9 orders | vector J_B, S-a Retino / PDG | dec1 | 11.95 / 16.85 |
| Λ*² ≈ 7 GeV² | 2/(κ · β_obs/β(ξc = 1)) | sig1 | 7.026 |
| +8 orders | Pauli weakest corner: log10(m_γ,Retino/m_γ,Proca gal) | floor0 | 12.90 |
| +15 orders | 3 log10(Λ_s,NDA/T_BBN) | floor0 | 15.70 |
| 3.9 orders vs Cassini | 036, same key as CLAIMS | dec1 | 3.873 |
| 23.7 orders short in size | 036, same key as CLAIMS | dec1 | 23.72 |

- `src/reproduce.py`: docstring only. 036 and 037 `run_checks.py` were already picked up by `--checks` (`AUTO_CHECKS_AFTER = 35`).

Command (repo root, OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS = 1), run once:
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning -m src.reproduce --log work/20261121-041-reproduce-036-037/run.log
```

## Result
`run.log` (path-free, 0 hits): **exit 1**, 604 s (CPU shared with two other runs). The only failure is one draft number.
- Validation 18/18 PASS. Every `_checks_ok` flag passes, including 036-mode, 036-end and 037-sf.
- Live CLAIMS numbers: **79/79 match** (was 74), including all 5 new EC-06d rows. Superseded 4/4.
- Draft numbers: 22/23 match. **MISMATCH: "+8 orders"**. The sentence says the weakest corner (Pauli, Stückelberg, Retino weak-end m_γ, 1 MeV cutoff, 1.2 orders short) ignores the Proca galactic limit (+8 orders). Swapping the Retino mass for the Proca galactic mass at that corner moves 1.2 to 14.1, so it adds **12.9 orders**. The 8.0 is log10(m_γ,PDG/m_γ,Proca gal), which is the shift from the PDG row (6.1 to 14.1), not from the corner. 037's README has the same wording.
- "+15 orders" matches only as a floor: the shift is 15.7 (3 log10(1.72×10¹¹ eV / 1 MeV)), and it rounds to 16.
- The pre-existing "c ≈ 10³⁶" is still "not in draft" (non-fatal).

## Verdict
Every 036/037 number in CLAIMS.md and in the draft is now regenerated from the work scripts. All match except the draft's "+8 orders", which understates the Proca shift at the corner it describes: 12.9 orders, not 8. The kill is stronger than the draft says, so the verdict does not change. The text needs fixing, either to "+13 orders" or to "+8 orders relative to the PDG mass". No CLAIMS or draft edits were made.

## Open issues
1. "+8 orders": overseer to fix the draft wording (see Result). After the fix, update `expected` and `digits` of `EC-02,corner_shift_Proca` in `src/draft_expected.csv`. If the PDG-relative reading is chosen, also change its getter to `log10(mg_PDG/mg_gal)`. Reproduce exits 0 once this row matches.
2. "+15 orders" is a floor of 15.7. "+16" would be the rounded value.
3. The vector "12.0" is 11.95, so it matches only after rounding to one decimal.
4. The `EC-02` draft rows sit under claim id EC-02 even for the Pauli and vector routes, following 029's existing rows. CLAIMS.md does not carry 037's numbers yet (037's proposed wording is not in the EC-02 row).
5. "α ≥ 1/3" is regenerated as the a1 → ∞ limit of α = (a1+1)/(3a1). That α > 1/3 holds for all a1 > 0 is 036's sympy check B8f, which is part of the 036-mode flag.
