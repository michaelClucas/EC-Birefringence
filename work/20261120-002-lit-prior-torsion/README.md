# 20261120-002 Literature scout: prior work on torsion and birefringence

## Question
What is already published on torsion-induced birefringence, Einstein-Cartan cosmology, Einstein-Cartan-Kaluza-Klein and Kalb-Ramond torsion axions, and the gauge-invariance problem? The goal is to keep the program from rederiving a known result and calling it new. Which parts of our question does existing work already settle?

## Method
- INSPIRE-HEP API (title and abstract queries sorted by citations) and arXiv API (`abs:` boolean queries). Scripts are in this folder:
  - `inspire_search.py "<query>"`
  - `arxiv_search.py "<search_query>"`
  - `arxiv_verify_ids.py <chars> <ids...>`
  - Run with `PYTHONIOENCODING=utf-8` through `tools/pin.py`. Note: pin.py needs the absolute path of `.venv/Scripts/python.exe` as the child command; a relative path fails with WinError 2.
- Query topics: torsion AND birefringence; torsion AND CFJ / Chern-Simons / Lorentz violation; Kalb-Ramond AND birefringence / optical activity / extra dimension; Barbero-Immirzi / Holst / Nieh-Yan AND photon / anomaly / axion; Kaluza-Klein AND torsion / birefringence; "cosmic birefringence" AND (torsion, KR, Holst, Immirzi, Nieh-Yan, KK, extra dimension, dilaton, string axion, neutrino current); Weyssenhoff; four-fermion; Poplawski; Mavromatos torsion axion; teleparallel birefringence; one extended web search for "torsion explanation of cosmic birefringence".
- Verification: every arXiv ID cited was fetched with the arXiv API `id_list` (title, authors, journal ref matched). Pre-arXiv classics (CFJ 1990, Hojman et al. 1978, de Sabbata-Gasperini 1981, Benn-Dereli-Tucker 1980, Nieh-Yan 1982, Hehl et al. RMP 1976, Duncan-Kaloper-Olive 1992) were verified by DOI via INSPIRE.
- Read depth:
  - Key sections read from the full text (arXiv HTML or ar5iv, tool-assisted extraction of the action, formulas and numbers; not a line-by-line read): 0908.0629, astro-ph/9706044, gr-qc/0307063, 0810.2355, gr-qc/0006097, hep-th/0507210, 1502.03694, 2608.09299, 2611.07578, 2604.03973 (partially; 75 pp).
  - Abstract only: everything else, about 70 papers. Tagged [ABS] in the literature file.
  - The raw PDFs fetched by the tool were not parsed: no PDF library is in `.venv`, and I did not add one.

## Result
- `literature/prior_work_torsion_birefringence.md`: about 80 papers in sections A (torsion-photon couplings; CFJ from fermion loops; SME torsion bounds; Holst/BI/Nieh-Yan; PGT/teleparallel/metric-affine), B (EC cosmology), C (EC-KK, KR torsion axion, other KK), D (gauge invariance) and E (fermion-current birefringence). Each entry gives the action or coupling, the result and formula where available, and relevance to the tracks. Flags are SAME-AIM and NO-GO.
- `literature/open_gaps.md`: 10 task candidates. G1-G3 are highest value.

Main findings (candidate statements, from the literature, not yet checked by the team):
1. No paper derives β from torsion in EC (or an EC extension) and fits it to Planck PR4 or ACT DR6 with lab torsion bounds applied.
2. Partial SAME-AIM papers:
   - Das-Mohanty-Prasanna 0908.0629: a dimension-5 non-minimal coupling, fit to WMAP+BOOMERanG; predicts β ∝ ν².
   - Maity-SenGupta-Sur hep-th/0507210: a KK-reduced massless KR axion, compared with radio galaxies; needs ρ_KR ~ ρ_DM for 0.2°.
   - Rao et al. 2608.09299: 5D teleparallel KK, where the radion multiplies an EM Chern-Simons term; fit to Planck PR4, with no α_em-variation or ACT constraint in the parts read. This is the closest prior art for the KK track.
3. The fermion route is published in pieces:
   - The fermion loop turns axial torsion into a CFJ term (Dobado-Maroto; regularization dependent per Jackiw-Kostelecky).
   - The resulting current-photon operator has a 2026 no-go (Reig 2611.07578) even at the 10 TeV scale. EC gives roughly the Planck scale.
4. The torsion-dual axion (KR/BI) has f ≈ 4×10¹⁸ GeV and a generic QCD coupling. That gives m ~ 10⁻¹² eV (1502.03694), which points to washout.
5. Weyssenhoff spin fluids are incompatible with FRW (gr-qc/0601089).

## Verdict
The literature does not already answer the program's question, so the program is not a rederivation. Several sub-results are published and must be cited, not rederived: the Dobado-Maroto coefficient, the Itin-Hehl FRW axion θ ∝ uv, the Maity et al. KR rotation, Rao et al.'s KK radion CS rotation, and Reig's operator no-go. The fastest route to a 4D-track Result B is G1 + G2 in `open_gaps.md`. For the KK track it is G3.

## Open issues
- Only abstracts were read for most papers. The [ABS] entries could hide details, especially Kohler gr-qc/9808004 and Shankar-Wali 0904.4038 (EC-KK), Kostelecky-Russell-Tasson 0712.4393 (exact torsion components and bounds; bounds analyst should read in full), and de Sabbata-Gasperini 1981 (no arXiv; DOI only).
- Rao et al. 2608.09299: about 35k characters at the end were not read. Their constraint discussion may be there.
- Das et al. 0908.0629 defines T₁ inconsistently, and the units of its bound are unclear. Do not import the number without rederiving it.
- The order-of-magnitude estimates in the literature files (S ~ 10⁻⁶⁸ eV from relic asymmetries, versus ~10⁻³² eV needed) are scout estimates, not checked derivations. They must not be quoted as results.
- Formulas are in each paper's conventions and not converted to `CONVENTIONS.md`.
- Coverage gaps: the teleparallel and metric-affine GW-birefringence literature is large and only sampled. Mavromatos's KR/RVM series is only partly listed. Pre-1990 work on photon-torsion was found only through INSPIRE title search.
