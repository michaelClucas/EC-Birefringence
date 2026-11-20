# Novelty audit of the program's "new" items

Task 20261120-024 (literature scout), 2026-11-20. Answers skeptic 020 MAJOR item 4 (and the literature part of MINOR 6). It builds on `literature/nogo_literature_check.md` (task 013) and `literature/ec06_prior_art.md`. Logs: `work/20261120-024-lit-ec06-novelty/` (`abstracts1-4.log` verify every arXiv ID by API; `inspire_kk.log`, `citing*.log`, `inspire_search*.log` hold the searches).

Labels:
- **new (not found)**: no paper found that states it. This is a search result, not proof of novelty.
- **known result applied**: a published general result, applied here to torsion or KK. The application may be new; the physics is not.
- **published**: already stated in the literature for this setting.

Every finding is a candidate mechanism or a no-go for one. ACT β is 2.9σ with systematics caveats.

## Summary

| # | Item (CLAIMS id) | Label | Cite as |
|---|---|---|---|
| 1 | KK P′ all-orders argument (KK-01, KK-04 background) | **known symmetry argument applied; the statement for KK photon birefringence is new (not found)** | P′ = diag(1,-1,-1,-1,-1) ∈ SO(4,1)₀; 4D parity in KK as a 5D proper transformation. Nearest literature: hep-th/0207165 (parity conservation of the U(1) sector with KR torsion in ADD/RS), Witten NPB 186 (1981) 412 (KK on manifolds with isometries gives a vector-like, parity-symmetric 4D spectrum). 2608.09299 is consistent (parity-even at c = 0) but does not state it |
| 2 | EC shortfall numbers: EC-02 28.5 orders, EC-05 ≥ 29 orders | **known result applied; the numbers are new** | qualitative: Dobado-Maroto astro-ph/9706044 (needed torsion ≫ spin-sourced torsion); general: Nakai-Namba-Obata 2310.09152 (no SM/LEFT operator, four-fermion included, gives β; minimal EC + Dirac is GR + an axial four-fermion contact term); Reig 2611.07578 (j·AF̃ route excluded) |
| 3 | Torsion-pseudoscalar photon decoupling (EC-03) | **known result applied; partly published for torsion** | general: Bauer-Neubert-Thamm 1708.00443, Quevillon-Smith 1903.12559, Bonnefoy-Di Luzio-Grojean 2011.10025, Quevillon-Smith-Vuong 2112.00553; torsion, for gluons: Karananas-Shaposhnikov-Zell 2506.11836 |
| 4 | Cassini kill of radion × FF̃ (KK-04, the 2608.09299 mechanism) | **known result applied; the application to 2608.09299 is new** | light gravitationally coupled modulus excluded by Cassini: Bertotti-Iess-Tortora Nature 425 (2003) 374, Will 1403.7377; KK scalar as Brans-Dicke field and its solar-system exclusion: Overduin-Wesson gr-qc/9805018 Secs. 3.5, 5; rolling φFF̃ field vs fifth forces: Carroll astro-ph/9806099 |
| 5 | KK-05 killed by the charge problem | **known result applied (classic)** | KK charge-to-mass relation: Overduin-Wesson gr-qc/9805018 Sec. 4.2 ("electron mass ... the Planck mass", ~22 orders), Bailin-Love Rep. Prog. Phys. 50 (1987) 1087; KK modes as the standard example of the q/m relation in the weak gravity conjecture, hep-th/0601001 |

Missing from the draft and to be added (skeptic 020 item 4): Das-Mohanty-Prasanna 0908.0629 (torsion-photon birefringence fit), and for EC-06 the published p ∝ H phenomenology (Li-Xia-Li hep-ph/0611192, Yu-Yang-Qiu 2508.10532). See `ec06_prior_art.md`.

## Item details

### 1. KK P′ all-orders argument
- **Program claim (009):** P′ = diag(1,-1,-1,-1,-1) is a proper rotation (det = +1, a rotation by π in the (x¹,x²) and (x³,y) planes), so it is a symmetry of any diffeomorphism-invariant 5D theory on S¹ or S¹/Z₂. Radion × F_KK F̃_KK is P′-odd, so it is absent at all orders in a P′-invariant state.
- **Literature:** 013 found no paper with this argument. 020 did not either. This task searched INSPIRE (titles: parity AND Kaluza-Klein; "parity violation" AND extra dimension; full text KK AND parity AND Chern-Simons AND radion AND pseudoscalar, 36 hits screened; `inspire_kk.log`) and found nothing.
  - The ingredients are textbook: 5D local Lorentz invariance contains P′, and y → -y is an isometry of S¹.
  - Defining 4D parity in 5D models as (x⃗, y) → (-x⃗, -y) is common practice in 5D model building, but **I found no citable paper that states it**.
  - Closest published statements:
    - Mukhopadhyaya, Sen, SenGupta hep-th/0207165 [FT in 013]: brane photon, KR torsion, ADD/RS. The U(1) sector stays parity conserving.
    - Witten 1981 (INSPIRE 10244, NPB 186 412, pre-arXiv, not opened): KK reduction on manifolds with isometries gives vector-like 4D fermions. That is the same "the reduction inherits a parity" logic for fermions.
- **Wording for RESULT:** "a standard symmetry argument (P′ is in the identity component of SO(4,1)); we found no paper applying it to KK photon birefringence."

### 2. EC shortfall numbers (EC-02, EC-05)
- **Literature:**
  - Dobado-Maroto astro-ph/9706044 [KEY in 002]: their fermion-loop CFJ needs |S| ~ 10⁻³⁰ eV, which they say requires torsion not sourced by ordinary spin. This is the qualitative statement.
  - Nakai-Namba-Obata 2310.09152 [FT in 013]: no SMEFT/LEFT operator, four-fermion terms included, gives the ICB angle; the CνB is too small. Minimal EC with Dirac matter reduces to exactly such an operator after torsion is eliminated, so the 4D spin-sourced no-go is a special case of their result.
  - Reig 2611.07578 [ABS, verified]: the neutrino-current CS operator needs Λ* ≈ 10 TeV and is excluded.
- **Not published:** the specific numbers (ξ ≈ 3×10²⁸, 28.5 orders; EC-05 ≥ 29 orders) and the explicit chain κJ₅ → S → CFJ → β.
- **Wording:** "an EC-specific instance of the 2310.09152 no-go, with the shortfall quantified (numbers new)."

### 3. Torsion-pseudoscalar photon decoupling (EC-03)
- **General result:** a pseudoscalar coupled to fermions only through ∂_μφ J₅^μ has an exact shift symmetry, and its induced photon coupling is g_eff ∝ B₁(τ_f) → -m_φ²/(12m_f²) for m_f ≫ m_φ, i.e. it decouples. Sources:
  - Bauer-Neubert-Thamm 1708.00443 [ABS; used for the loop function in task 012].
  - Quevillon-Smith 1903.12559 [ABS]: axion couplings come from non-anomalous loops, "blind to anomalies".
  - Bonnefoy-Di Luzio-Grojean 2011.10025 [ABS].
  - Quevillon-Smith-Vuong 2112.00553 [ABS].
- **Torsion-specific:** Karananas-Shaposhnikov-Zell 2506.11836 [FT in 013]: for a scalar of gravitational (torsion) origin with purely derivative fermion couplings, "no effect, perturbative or non-perturbative" breaks the shift symmetry. The gluon anomaly term and the mass term cancel (their Eqs. 9-11), so there is no QCD potential. The photon version is the same identity with QED in place of QCD. 013 did not find it stated for photons, and I did not re-read 2506.11836 for that here.
- **Wording:** "the standard ALP-EFT derivative-coupling decoupling (1708.00443, 1903.12559, 2011.10025), applied to torsion pseudoscalars; for gluons this was shown by 2506.11836." The C = 96 coefficient and the (m_φ/m_e)²/96 numbers are the program's own (012); only the application is new.

### 4. Cassini kill of radion × FF̃ (KK-04)
- **General result:**
  - An unscreened light scalar with gravitational-strength coupling to matter is excluded by Cassini γ - 1 (value in `data/kk_bounds_inputs.csv` row gamma_minus_1_Cassini; Bertotti-Iess-Tortora 2003, INSPIRE 639311; Will 1403.7377 [ABS]).
  - The KK scalar is a Brans-Dicke field with small ω, excluded by solar-system tests (Overduin-Wesson gr-qc/9805018 [FT], Sec. 3.5: ω must exceed ~500; Sec. 5: Viking time delay).
  - Carroll astro-ph/9806099 [ABS] is the standard paper on rolling-field couplings to matter versus the FF̃ coupling.
- **Application:** 2608.09299 [FT in 013] says moduli fixing is "beyond the scope" and applies no PPN/fifth-force/α-variation bound. INSPIRE showed 0 citing papers in 013, and none appeared in this task's searches.
- **Wording:** "known solar-system exclusion of a light KK radion, applied to the 2608.09299 mechanism (application new)."

### 5. KK-05 charge problem
- **General result:**
  - A state with KK charge has q/m fixed by the compactification: charge e needs m ~ eM_P. This is the classic reason 5D KK was abandoned as a theory of the electron (Overduin-Wesson gr-qc/9805018 [FT] Sec. 4.2, "some twenty-two orders of magnitude"; Bailin-Love 1987, INSPIRE 253064, not opened).
  - In the weak-gravity literature KK modes are the standard example of q/m set by M_P (Arkani-Hamed-Motl-Nicolis-Vafa hep-th/0601001 [ABS]).
- **Application:** KK-05 is the program's own construct. Its kill (the CMB photon must be the KK photon, so the electron must carry KK charge) is this classic result. 016/019's EDM analysis and the brane-photon mixing loophole are program-specific.
- **Wording:** "fails on the classic KK charge-to-mass problem."

## MINOR 6 (literature part): operators beyond CFJ
- Gauge-invariant torsion-photon operators are catalogued in:
  - Itin-Hehl gr-qc/0307063 (F²T², 17 terms).
  - Bahamonde-Maggiolo-Pfeifer 2507.02362 (linear F·R̃_{[μν]}, F·T·T; Eq. 17, Eq. 23 has F·S·T pieces).
  - De Felice et al. 2607.26141 (metric-affine constitutive framework).
  - Das-Mohanty-Prasanna 0908.0629 (F∂F̃T).
- No paper evaluates any of them on an FRW torsion background with numbers, except DMP: a constant, unspecified T₁, with β ∝ ν², disfavoured at ≈ 4.9σ by Planck PR4 chromaticity (`ec06_prior_art.md` Sec. 4). The S3 numbers of skeptic 020 for spin-sourced torsion are new.

## ID verification
All 78 arXiv IDs in this file, `ec06_prior_art.md` and `data/kaf_lowz_bounds.csv` were returned by one arXiv API `id_list` query on 2026-11-20 with matching titles and authors (`abstracts5.log`). Pre-arXiv items (Witten 1981, Bailin-Love 1987, Bertotti et al. 2003, CFJ 1990, Kanazawa et al. 2018) were checked on INSPIRE by DOI or title.
