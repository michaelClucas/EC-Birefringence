# 20261120-013 Literature check of the no-go; ultralight energy-density bound

## Question
1. Has any paper already published (a) minimal 4D EC and pure EC-KK give no birefringence, (b) torsion-sourced CFJ routes are tens of orders short, (c) torsion axions either get a QCD mass and wash out or are generic ALPs with f ~ M_P, (d) the 2608.09299 radion × FF̃ route is excluded by Cassini / ISL? Does anything contradict these?
2. Replace task 008's uncited 5% ρ_DM cap with a cited bound for an ultralight field frozen at recombination, m ~ 10⁻³³-10⁻³² eV.
3. Do 2601.13624 or 2512.21888 use a torsion-motivated coupling, and what do they say about EC-03', KK-02, KK-03?

## Method
- `lit_search.py`: arXiv API keyword search, INSPIRE query (`refersto:recid:`, `fulltext:`), and arXiv `id_list` verification with abstracts. Network only. Logs: `arxiv_search.log`, `inspire_fulltext.log`, `citing.log`, `abstracts1-4.log`.
- Full texts read from arXiv HTML (downloaded to a scratch folder, not committed) for 2310.09152, 2506.11836, hep-th/0207165, 2608.09299v3, 2601.13624, 2512.21888, 2303.02852, 2503.18924, 1410.2896, 1708.05681, 2301.08361, 2503.14454, 2511.18917. Table values were extracted from the HTML `alttext` of the table cells (1708.05681 Table 3, 2301.08361 Table 2) and checked against the row order.
- `energy_check.py`: independent spectator-field ODE (written from scratch, inputs from `data/`) for the density EC-03' needs, compared with the new bounds. Log: `energy_check.log`.

Rerun:
```
PYTHONIOENCODING=utf-8 .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-013-lit-nogo-check/lit_search.py ids 2506.11836 2310.09152
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-013-lit-nogo-check/energy_check.py
```

## Result
Full write-up: `literature/nogo_literature_check.md`. Data: `data/ultralight_fraction_bounds.csv`.

1. No paper computes cosmic birefringence from EC torsion. Pieces exist separately:
   - hep-th/0207165: fully antisymmetric (KR-dual) torsion leaves the U(1) sector parity conserving in 4D and in ADD/RS compactifications; minimal A-torsion coupling breaks U(1). Partial confirmation of EC-00; KK-01 (geometric KK photon, P' argument) not addressed.
   - 2310.09152: no SM/LEFT operator gives the isotropic β; gravity omitted. Minimal EC = GR + four-fermion contact, so EC-00/02/05 fall under it (our mapping).
   - 2506.11836: a torsion scalar coupled only via ∂φ·J_A has an exact shift symmetry and gets no potential; earlier gravitational-axion constructions (1502.03694, Mercuri, Lattanzi-Mercuri) are called flawed.
   - 2608.09299: stabilisation of the radion "beyond the scope"; no Cassini/ISL/α_em bound applied; 0 citations on INSPIRE.
2. Per-mass bounds: m = 10⁻³³ eV: none (degenerate with Λ, 1708.05681). m = 10⁻³² eV: Ω_a h² < 0.014 (Planck 2015, 1708.05681), < 0.00833 (Planck 2018) and < 0.00321 (Planck 2018 + BOSS, 2301.08361), i.e. 11.8%, 6.9%, 2.7% of ρ_DM at 95%. 1806.10608: no bound for 1+z_c ≲ 10. EC-03' needs ρ_φ,0 = 0.11 ρ_DM at m ≈ 6×10⁻³³ eV (reproduces 008), ≤ 0.21 only for 5-8×10⁻³³ eV, 0.19 at 10⁻³² eV. 3/3 PASS in `energy_check.log`.
3. Neither 2601.13624 nor 2512.21888 mentions torsion, KK or KR. 2601.13624 g M_Pl = 0.129 (+0.031 −0.075) vs torsion g_T M_P = 0.011 (reduced) or 0.057 (non-reduced); not directly applicable since EC-03' is not an EDE field.

## Verdict
- No-go published: **partly**. Nothing contradicts any status in `CLAIMS.md`.
- **Contradiction of reasoning:** 2506.11836 says the torsion axion gets no QCD potential from a derivative coupling, which removes EC-03 "reading A" (QCD mass 1.4-8.6×10⁻¹² eV from 1502.03694). EC-03 stays dead through reading B. By the same cancellation applied to the photon anomaly (our extension, not in the paper), EC-03' has no torsion-generated g_T; its coupling must be added by hand, making it a generic ALP at f ~ M_P.
- **EC-03' energy cap:** at 10⁻³² eV the required 19% is excluded by Planck 2018 (6.9%). In the 5-8×10⁻³³ eV window there is no published per-mass bound; the field is DE-degenerate there and the relevant constraint is a late-time BAO+SN fit, not a DM fraction. Conservative choice (nearest tabulated point, 2.7%): β_max ≈ 0.10°. Candidate mechanism status unchanged: marginal at best.

## Open issues
1. Theorist: confirm that the 2506.11836 shift-symmetry cancellation (their Eqs. 9-11) also removes φFF̃ from a derivative torsion coupling with massive charged fermions. If yes, update EC-03 rationale and EC-03' "torsion origin of g".
2. Data analyst: fit a sub-dominant thawing field + Λ with m = 5-8×10⁻³³ eV to DESI DR2 BAO + SN (2503.18924 shows the region is allowed at 1σ but gives no bound) to obtain an actual cap for EC-03'.
3. 2406.11956 (Weyl-invariant EC ALP) read at abstract level. If it solves strong CP it has a QCD mass and washes out; not checked for a photon coupling.
4. Shapiro hep-th/0103093 and Kostelecky-Russell-Tasson 0712.4393 read at abstract level only.
5. 1806.10608 bound read at abstract level ("degenerate with DE if 1+z_c ≲ 10").
