# 20261121-042 data: final fixes for skeptic 040

## Question
Close skeptic 040's data items: (5) add arXiv:2608.06480 and refit penalties against it; (6) EC-04 upper edge from Chandra instead of CAST; (8) cite the energies behind "photon physics standard" (BaBar 10.58 GeV of skeptic 039; 10 TeV); (3) put KK-05 "~21 orders" and 020's higher-dimension cutoffs into `src.reproduce`.

## Method
- Opened 2608.06480 (abstract and arXiv HTML full text), BaBar 1702.03327 (PDF pp. 1-6), ATLAS 2102.13405 (HTML), ATLAS 1702.01625 (abstract). Rows: `data/birefringence_measurements.csv` `ESK26-*`, new `data/collider_photon_bounds.csv`; notes in `literature/2608.06480.md` and `data/README.md`.
- `src/likelihood.py`: datasets `act_planck_esk26` (0.277 ± 0.057) and `act_planck_esk26_dust` (0.236 ± 0.067), components {ACT_DR6, PLANCK}, so `joint()` refuses products with the DPK25 combinations (checked).
- `src/reproduce.py`: both added to `FIT_DATASETS`; 014's `fit_candidates.main` fits them unchanged (fit summary prints them).
- `src/claims_numbers.py`: new sources `020-skep` (020 `skeptic_checks.py`, in-process) and `042-fix` (this folder's `run_checks.py`); new `DRAFT_NUMBERS` and `fit_values` keys; 13 rows added to `src/draft_expected.csv`.
- `run_checks.py` (read-only, 22 PASS / 0 FAIL, `run.log`): C1 EC-04 Chandra edge, g_φγ = ξ/√A with A = 3/(2κ), i.e. 0.8165 ξ/M_P (007 `ec04_estimates`), g from `data/bounds.csv` A-Chandra; CAST version reproduces 007's 1.7×10⁸. C2 KK-05 e M_P/√2 recomputed independently. C3 020's hardcoded ξ_EC02 = 2.9×10²⁸ equals 007's regenerated value. C4 EC-06a Λ_SC against the ATLAS 2.36 TeV diphoton maximum. C5 hand Gaussian penalties vs the reproduce fit rows (all 14 ESK26 rows agree to 10⁻³).
- `src.reproduce` (no `--checks`), pinned, 1 thread: `reproduce.log`.

## Result
2608.06480 (Eskilt, joint ACT DR6 + Planck PR4 NPIPE, cross-spectra included, no WMAP; 68% CL): β = 0.277° ± 0.057° (4.8σ), dust-mitigated 0.236° ± 0.067° (3.5σ). Shares ACT DR6 and Planck PR4 data with DPK25-comb-EK22/CG; used as an alternative combination only.

| quantity | baseline | dust-mitigated | draft string (baseline / dust) |
|---|---|---|---|
| free ALP Δχ² vs null | 23.62 | 12.41 | `Δχ² = 23.6 vs null (2608.06480)` / `Δχ² = 12.4 vs null (2608.06480 dust-mitigated)` |
| EC-03' 2.7% cap penalty | 9.67 | 4.13 | `penalty 9.7 (2608.06480)` / `penalty 4.1 (2608.06480 dust-mitigated)` |
| EC-03' 6.9% cap penalty | 4.16 | 1.26 | `penalty 4.2 (2608.06480)` / `penalty 1.3 (2608.06480 dust-mitigated)` |
| KK-02 k = 1 penalty | 1.42 | 0.16 | `penalty 1.4 (2608.06480)` / `penalty 0.16 (2608.06480 dust-mitigated)` |

- EC-04: Chandra (g < 8×10⁻¹³ GeV⁻¹, m < 10⁻¹² eV) gives **ξ < 2.39×10⁶**, 1.86 orders below the CAST edge (1.7×10⁸, printed 2×10⁸). Window still open (10⁻¹⁴ ≲ ξ). Skeptic 040's "~3×10⁶" was a rescaling of the rounded 2×10⁸; the direct value is 2.4×10⁶. Draft string: `10⁻¹⁴ ≲ ξ < 2.4×10⁶`. Suggested CLAIMS text: `10⁻¹⁴ ≲ ξ < 2.4×10⁶ (Chandra)`; the `EC-04.xi_max_Chandra` key is regenerated now and can move to `claims_expected.csv` once CLAIMS.md carries the text (the old `xi_max_CAST` row then becomes superseded).
- BaBar 1702.03327: single-photon (e⁺e⁻ → γ + invisible) data at √s = 10.58 GeV (Υ(4S), 5.9 fb⁻¹), 10.36 GeV (Υ(3S)) and 10.02 GeV (Υ(2S)), 53 fb⁻¹ total (p. 4, Table I). Skeptic 039's 10.58 GeV is confirmed.
- 10 TeV: no single measurement; it is a convention. Cited instead: ATLAS 2102.13405, diphoton spectrum consistent with the SM, highest event m_γγ = 2.36 TeV; ATLAS 1702.01625, γγ → γγ at m_γγ > 6 GeV. Against 2.36 TeV, EC-06a's Λ_SC = 3.6 keV is **8.8 orders** below (9.4 against 10 TeV); every 10 TeV shortfall drops by 0.63 orders. Draft string: `8.8 orders in energy below 2.36 TeV`.
- KK-05: e M_P/√2 = 5.2×10²⁶ eV = 1.0×10²¹ m_e, 21.01 orders; matches the draft's `~21 orders in the electron mass`.
- 020 cutoffs: Λ = 1.96×10⁻⁴³ eV ((∇·S)FF̃/Λ²) and 4.1×10⁻¹⁵ eV (DMP p²T operator); match the draft's `Λ ≲ 2×10⁻⁴³ eV / 4×10⁻¹⁵ eV`.
- `src.reproduce`: exit 0 (210 s); validation 18/18; claim numbers 79/79 live match; draft numbers 36/36 regenerated and matching, 10 not yet in the draft (the 8 ESK26 strings, the EC-04 Chandra string, the ATLAS string). Path scan clean.

## Verdict
All five data items closed. No verdict changes: against 2608.06480, KK-02 k = 1 sits 1.2σ low (penalty 1.4), EC-03' at the 2.7% cap is penalised 9.7 (4.1 dust-mitigated). EC-04 stays a live ALP-with-added-coupling route with a 1.9-order narrower window.

## Open issues
- `--checks` was not rerun here; this folder's `run_checks.py` is picked up automatically (folder number > 035), as is any `run_checks.py` task 043 adds.
- Task 043 (parallel) also wrote ESK26 into `data/alp_costs_inputs.csv` from the abstract and says the dust-mitigated value has no central value; the full text gives 0.236° ± 0.067° (Sec. III). The overseer may want one source of truth (`birefringence_measurements.csv`).
- The sign convention of 2608.06480 is not restated in the paper; taken equal to EK22's (same author, positive value).
- 039's `stueckelberg_cutoff.py` reads BaBar/LEP2 energies from its own local inputs file, not from `data/collider_photon_bounds.csv` (not my script to edit). LEP2 209 GeV is not cited here.
