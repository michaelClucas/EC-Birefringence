# 20261120-029 final skeptic on RESULT_DRAFT revision 2

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Every mechanism below is a candidate mechanism. There is no `RESULT.md`.

## Question
Does `work/20261120-018-result-draft/RESULT_DRAFT.md` (rev 2), with `CLAIMS.md` and the repo, meet the OVERSEER_PROMPT stop criteria (Result B or C per track, literature check, skeptic, reproduce)? Try to break it. The f(T) cost analysis (025) is audited separately by 028; here only: does the 4D answer hold under either 028 outcome?

## Method
- Read CLAIMS.md, CONVENTIONS.md, OVERSEER_PROMPT criteria, the draft, and the work/ READMEs behind the verdicts (007, 022, 025 numbers only). Formed my own list before opening the earlier skeptic READMEs (020, 026), then checked their fixes.
- Weakest links tested myself:
  1. EC-02/05 "spin-sourced torsion (any coupling) ≥ 28 orders short": only shown for cited sources (relic-ν asymmetry, polarised DM, baryons). `pauli_bound.py` bounds ANY fermionic axial source (a condensate or exotic dark fermion included) by the Pauli limit, and the vector-current source that a non-minimal (Holst/Freidel-Minic-Takeuchi) coupling allows.
  2. KK-04 Cassini γ - 1 = -1/2: hand check of the radion matter coupling (bulk vs brane matter).
  3. KK-02 β_max = kα/2: hand reduction of the level-k 5D CS term.
  4. The MZ (EC-06a) status against the 4D verdict wording (logic check, no new numbers).
- Bounds opened at the source (arXiv abstract or full text, 2026-11-20): 1312.3225 (neutron b_⊥), 2212.11841 (JILA eEDM, full text p. 1), 2002.11761 (Eöt-Wash ISL), 1710.05834 (GW170817 speed), 2209.15487 (MICROSCOPE), 1403.7377 (Will review, Cassini).
- Literature: independent searches for (i) torsion/teleparallel origin of a cosmological CFJ vector for CMB birefringence, (ii) the KK P′ argument; abstracts of 0908.0629, 2409.15509, 2412.05982, 2610.02667 opened.
- Honesty: README.md, CLAIMS.md, draft, `git log --format=%B`.
- Reproduce: `python -m src.reproduce --checks` once, pinned, 1 thread: `reproduce_checks.log` (scrubbed log) and `reproduce_checks.stdout`.

Rerun my script (seconds, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-029-final-skeptic-rev2/pauli_bound.py
```

## Findings

| # | item | severity | evidence | required fix |
|---|---|---|---|---|
| 1 | **4D verdict contradicts EC-06a.** The draft says "Result B for Einstein-Cartan proper; Result C only outside EC" and "The one surviving torsion-based candidate mechanism is f(T) ... (EC-06c)". But CLAIMS EC-06a (Magueijo-Złośnik, an EC extension: Λ as a field + Euler/Pontryagin terms) is "disfavoured, not excluded: no data excludes ξ ≳ 10⁴" (skeptic 026). Result B needs "the bound or theorem that kills each" named extension; MZ has none, only naturalness. So MZ is a second surviving candidate mechanism, and it is inside an EC extension. This holds under **both** 028 outcomes: if EC-06c dies, the draft would read as a clean 4D Result B while MZ is still not excluded. | **BLOCKER** | CLAIMS EC-06a row; 026 README rows 32-33 and open issue (decaying P ∝ a⁻³ mode, unchecked: off the tracker GW170817 may not bound ε(z*) at all, so even ξ ~ 10⁻³ may be allowed); 022 open issue 1: off-FRW closure of the MZ torsion "not derived (moot once dead)", and MZ is no longer dead. | Either (a) kill MZ with a bound: derive the off-FRW U(1) breaking for p = ξS_MZ and apply the 025 machinery (Stückelberg window, cutoff ~ (mM_P²/ξ)^{1/3}, which at ξ ≥ 1.4×10⁴ is ~3×10⁻³ of the f(T) value if the same structure holds, i.e. ~0.1 GeV), and check the decaying-mode branch; or (b) relabel: 4D = Result C in two theories (MZ EC extension with ξ ≥ 1.4×10⁴ on the tracker and an unquantified U(1) cost; f(T)), with each cost quantified. Remove "The one surviving torsion-based candidate mechanism is f(T)". |
| 2 | **"Result B" labels while live, non-torsion survivors produce β.** KK track: KK-02 and KK-03 produce β within all bounds in EC-KK + an added 5D CS / GS term. 4D: EC-03′, EC-04 (+ CFJ by hand). The OVERSEER definition of B is "each named extension examined cannot produce the observed β within independent bounds". The draft is honest that these are ALPs where torsion plays no role, but the one-line answers ("KK track: Result B") overstate. | MAJOR | Draft Answer, KK bullet; CLAIMS KK-02/03 "live". | One-line answers: "B for every route where torsion generates or sources the coupling; β is possible only in extensions that add an ALP-photon coupling by hand (KK-02, KK-03, EC-03′, EC-04), where torsion plays no role (formally Result C with costs: added CS/GS term, 2πLM ≈ 220-280, g window ×6)". |
| 3 | **Result B "shortfall in orders of magnitude" missing for several kills.** KK-04: γ - 1 = -1/2 vs Cassini 2.3×10⁻⁵ (1σ) is ~4.3 orders, not stated. EC-06d: "H(z*) = 33-66 × ΛCDM" is not translated into a bound with CL (which observable, what allowed δH/H). KK-05 kill is a mass gap (~10²¹ m_e), fine. | MAJOR (stop criterion B) | Draft Derivation table and Answer. | Add the orders-of-magnitude shortfall and the bound/CL for each killed item, in one column. |
| 4 | "Spin-sourced torsion (any coupling) ≥ 28 orders short" is shown only for cited sources. A homogeneous fermion condensate (Poplawski 1102.5667 is in `literature/` but never bounded) or an exotic chiral dark fermion was not covered. My Pauli bound closes it: for g helicity states carrying the **entire** ΛCDM energy density as a maximally chiral degenerate gas, β(ξ = 1) ≤ 6.9×10⁻²⁷ deg, ≥ 25.5 orders short (S/H ≤ 2×10⁻²⁸ at z*). So the no-go survives, but the stated general number is not 28. | MINOR | `pauli_bound.py`, `pauli_bound.log` (PASS). | Reword: "≥ 28.5 orders for cited sources; ≥ 25.5 orders for any fermionic source (Pauli bound)", or cite this script. |
| 5 | Non-minimal fermion-torsion couplings (Holst + Freidel-Minic-Takeuchi, hep-th/0507253, listed in `literature/` only) let the **vector** baryon current, which is homogeneous and nonzero, source S^0. Never examined. My estimate: needs coefficient c ≈ 1.7×10³⁶ at ξ = 1 (36 orders), which implies a baryon contact interaction at ≲ 2 GeV. Dead, but absent from the scope statement. | MINOR | `pauli_bound.log` part 2. | One line in the scope/derivation summary. |
| 6 | KK-04 applicability: γ - 1 = -1/2 is the Brans-Dicke ω = 0 value for y-independent bulk matter. For brane-localised matter my hand estimate gives α₀² = 1/6, γ - 1 = -2/7. Still excluded by Cassini by ~4 orders, so the verdict holds; ISL is quoted at \|α\| = 1 (2002.11761: λ < 38.6 μm, 95%) while the radion coupling is α ≈ 1/3 to 1/2 (slightly weaker bound). Irrelevant to the kill (a radion heavy enough to evade ISL cannot roll after recombination). | MINOR | hand derivation (README text only, not scripted). | Say "bulk matter; brane matter gives -2/7" and note the α ≠ 1 ISL caveat. |
| 7 | KK-01 "P′-odd matter backgrounds are Planck suppressed like their 4D analogues" is stated as fact. 020 item 5 asked for "(not computed in 5D)". | MINOR | Draft KK bullet; 020 row 5. | Add "(not computed in 5D)". |
| 8 | A_CB ≤ 5×10⁻¹¹ is a generic estimate assuming torsion follows linear matter perturbations (CLAIMS says so; 026 asked for the label). The draft derivation table quotes it without the label. | MINOR | Draft table row "Unsourced O(H) torsion". | Add "(generic estimate, not an f(T) prediction)". |
| 9 | Axion-baseline comparison omits KK-03 and EC-04 (both "live"). KK-03's g window 1.4-8×10⁻¹³ GeV⁻¹ has no penalty or β_max. | MINOR | Draft comparison table. | Add rows (or state they reach the free-ALP fit, with the reason). |
| 10 | Unresolved factor 4.47 between the β → k_AF map of KM Eq. 158 and the SME Data Tables row (005 `derived_bounds.csv`, "unresolved") is not listed. It only affects a cross-check with ≥ 13 orders of margin (lab S_0), so no verdict changes. | MINOR | `work/20261120-005-bounds/derived_bounds.csv` row kV00-MK-check. | List under known weaknesses or resolve (likely cosmic vs conformal time). |
| 11 | Literature: 2610.02667 (EB spectra for the ν-asymmetry current vs axion DE/DM) is not in `literature/`. Relevant to EC-02, which reduces to that operator, and to distinguishability. Does not change novelty labels. My searches (teleparallel/EC torsion + CFJ + CMB; KK parity + radion FF̃) found nothing that contradicts `novelty_audit.md`; 0908.0629 authors/claim and 2409.15509 (GW only) confirmed. | MINOR | WebSearch/WebFetch 2026-11-20. | Add 2610.02667 to `literature/`. |
| 12 | Stale text. Draft "Numbers": "EC-06c Stückelberg window, EFT cutoff and EC-06d H(z*) ratio are not yet in the table", but `src/claims_expected.csv` now has them. 018 README "Open issues: KK-05 EDM bound, KK-02 CS level normalization ... TBD" (settled by 016/019). Iteration-6 commit says "MZ ... dead" (superseded; history cannot change, so RESULT.md should state the correction). | MINOR | files named. | Edit. |
| 13 | Reproduce passes only with task 030's uncommitted src edits (57/57); at HEAD the 025 numbers are not covered. Cosmoglobe penalties and a few literature/skeptic numbers in the draft are not regenerated (list below). | MAJOR (stop criterion: reproduce every number in the claim) | `reproduce_checks.log`, `reproduce_rerun.log` | Commit 030; add the Cosmoglobe penalties to the fit summary; either regenerate or mark as cited the other listed numbers in RESULT.md. |

## Bounds checked at source
| bound | repo value | source says | verdict |
|---|---|---|---|
| neutron \|b̃_⊥\|, 1312.3225 | < 6.7×10⁻³⁴ GeV, 68% | "< 6.7×10⁻³⁴ GeV", 68% CL, ³He-¹²⁹Xe | correct. Used at 68% CL; 95% would be ~2× weaker, irrelevant (lab S_0 bound has ≥ 13 orders margin over the need) |
| JILA eEDM, 2212.11841 | < 4.1×10⁻³⁰ e cm, 90%; d_e = (-1.3 ± 2.0 ± 0.6)×10⁻³⁰ | same (full text p. 1) | correct. Only feeds KK-05 c ≤ 188, which no longer carries the KK-05 verdict |
| Eöt-Wash ISL, 2002.11761 | λ < 38.6 μm, \|α\| = 1, 95% | same | correct value; applicability caveat item 6 |
| GW170817 speed, 1710.05834 | -3×10⁻¹⁵ ≤ Δc/c ≤ 7×10⁻¹⁶ | same | correct; 026's use of the +7×10⁻¹⁶ side is right for c_T > 1. CL listed as "90 (approx)" in data; the abstract gives none, fine |
| MICROSCOPE, 2209.15487 | η = (-1.5 ± 2.3 ± 1.5)×10⁻¹⁵, 1σ stat | same | correct |
| Cassini γ, via 1403.7377 | (2.1 ± 2.3)×10⁻⁵, 68% | Will 2014 is the cited review | correct |

## Reproducibility
- `--checks` run (`reproduce_checks.log`, 1803 s): validation 18/18 PASS; derivation checks **8/8 scripts PASS** (003, 004, 009, 012, 016, 019, 022, 025; 281 PASS lines, 0 FAIL); path scan clean. Exit 1: 49 claim numbers match, 8 "NOT REGENERATED" (the 025 rows). Cause: a race, not a defect. Task 030 edited `src/claims_numbers.py` (14:49) after my run had imported it and `src/claims_expected.csv` (15:02) before the compare step, both uncommitted.
- Rerun of plain `src.reproduce` on the current tree (`reproduce_rerun.log`, 402 s): **57/57 live claim numbers match, exit 0.** This depends on 030's **uncommitted** src changes; at HEAD (b95b350) the 025 numbers are not covered. The overseer must commit 030 before a stop.
- Draft numbers NOT regenerated by reproduce (MINOR unless noted):
  - Cosmoglobe penalties 4.95 (EC-03′) and 0.14 (KK-02): the fit summary prints only ACT DR6 and ACT+WMAP+Planck (EK22). Validation prints the Cosmoglobe β but not the penalties.
  - Higher-dimension operator scales Λ ≲ 2×10⁻⁴³ eV and 4×10⁻¹⁵ eV (skeptic 020 S3, a skeptic's script, not in reproduce).
  - DMP β ∝ ν² disfavoured at ~4.9σ (literature/ec06_prior_art.md).
  - LiteBIRD 0.05° separation at ~3σ (023).
  - KK-05 m ≳ eM_P/√2 ~ 10²¹ m_e (stated in 019/020; not a claims_expected row).
  - Free-ALP Δχ² = 20.7 is printed and checked against the published 4.6σ in validation, fine.
  - Chromaticity n = 1.61 (2.0σ) and Planck n = -0.35 ± 0.48 are printed as INFO lines, fine.

## Verdict per track
- **4D: does not survive as stated.** Under either 028 outcome the line "Result B for EC proper; Result C only outside EC" and "the one surviving torsion-based candidate mechanism is f(T)" conflict with EC-06a (MZ, an EC extension, not excluded). Fix item 1 first. The EC-proper no-go itself (EC-00, -02, -03, -05) survives, with items 4-5 strengthening it (any fermionic source ≥ 25.5 orders short). If 028 kills EC-06c and item 1 is fixed by a real MZ kill, 4D becomes a clean Result B; if MZ is not killed, 4D is Result C (MZ and possibly f(T)) with costs to be quantified.
- **KK: survives with listed fixes** (items 2, 3, 6, 7). The torsion-generated routes are dead with cited bounds; the headline must say that β remains possible only through ALP couplings added by hand.

## Open issues
- MZ off-FRW U(1) and the decaying-mode branch (item 1) are the only open physics needed for a 4D stop.
- I did not re-derive the KK-02 CS normalisation symbolically; my hand reduction agrees with kαΔϑ/2π (level-k CS → (kα/4π)ϑFF̃ in -(g/4)φFF̃ normalisation).
- The radion brane-matter coupling (item 6) is a hand estimate.
