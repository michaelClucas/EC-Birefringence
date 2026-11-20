# 20261120-024 Literature: EC-06 prior art, low-z k_AF bounds, novelty audit

Every finding is about a candidate mechanism. ACT DR6 β = 0.215° ± 0.074° is 2.9σ with systematics caveats and may not be cosmological. There is no `RESULT.md`.

## Question
1. Is EC-06 already published? EC-06 is an unsourced O(H) background torsion (Magueijo-Złośnik varying-Λ EC, teleparallel, PGT) with a CFJ-type photon coupling. Prior art to check: MZ follow-ups, teleparallel/f(T)/Nieh-Yan birefringence, PGT photon couplings, Das-Mohanty-Prasanna 0908.0629 and its citers, and any CFJ vector tied to H.
2. What are the low-z astrophysical bounds on a constant or slowly varying k_AF?
3. Correct the program's novelty list (skeptic 020, MAJOR 4).

## Method
- `lit_search.py` (copied from task 013, with two modes added: `inspire_id`, `citeq`) queries the arXiv API and INSPIRE.
  - Citation lists: 1908.05184 (37 citing), 0908.0629 (19), 2001.06373, 2007.08038, 2104.05917, hep-ph/0611192 (47), 0810.0403 (128), with full-text filters.
  - INSPIRE full-text and title queries for teleparallel/Nieh-Yan/f(T)/PGT/metric-affine + birefringence, CFJ + teleparallel, CFJ + Hubble + torsion, and KK parity.
  - Logs: `citing*.log`, `inspire_search*.log`, `inspire_kk.log`, `arxiv_search*.log`. The second query in `citing2.log` and the `refersto:arxiv:` queries in `inspire_search.log` used a bad syntax and returned all of INSPIRE; ignore them.
- Read the relevant sections of the arXiv HTML/ar5iv full text (tool-assisted extraction, not line by line) for: 1908.05184, 2001.06373, 0908.0629, 2507.02362, 2609.10671, hep-ph/0611192, 2508.10532, 1003.4823, 1411.6287, astro-ph/9704263, 1409.8242, 1702.01700, 2611.23082, gr-qc/9805018.
- All 78 arXiv IDs cited in the outputs were verified in one arXiv API `id_list` call (`abstracts5.log`; per-batch abstracts in `abstracts1-4.log`). Pre-arXiv items were checked on INSPIRE.
- `kaf_lowz.py` loads the published values from `data/kaf_lowz_bounds.csv` and H0, Ω_m, z*, ACT β from `data/`. It writes and checks the converted columns.

Rerun (repo root, a few seconds, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-024-lit-ec06-novelty/kaf_lowz.py
```
Output in `run.log` (PASS). Searches: `PYTHONIOENCODING=utf-8 <same prefix> work/20261120-024-lit-ec06-novelty/lit_search.py citeq <arXiv id> "<filter>"`.

## Result
Outputs: `literature/ec06_prior_art.md`, `literature/novelty_audit.md`, `data/kaf_lowz_bounds.csv` (+ `data/README.md` section).

1. **EC-06 prior art: partly published.**
   - The p ∝ H CFJ phenomenology is published:
     - Li-Xia-Li hep-ph/0611192: p ∝ ∂_μ ln R, WMAP3+B03.
     - Yu-Yang-Qiu 2508.10532: ln Q with Q = 6H², Planck β, α ≈ 0.8-1.8×10⁻⁴.
     - Both are achromatic, have a log history, use a dimensionless coupling of 10⁻⁴-10⁻², and need no new light field.
   - In matter domination both equal the TEGR trace-torsion CFJ (p = ξT, T_0̂ = 3H).
   - Not found: any paper deriving such a CFJ vector from an unsourced torsion background.
   - MZ torsion exists only in an EC extension (Λ a field plus Euler/Pontryagin terms). No follow-up couples it to photons. Its GW-speed bound (2001.06373 Eq. 68, γ² ≲ 10⁻¹⁵) gives ε ≲ 2×10⁻⁸ (scout estimate), so ξ ≳ 5×10⁴.
   - Teleparallel and Nieh-Yan papers are GW-only in 4D.
   - PGT: 2507.02362 only suggests birefringence as future work.
   - DMP 0908.0629: unspecified constant torsion, operator F∂F̃T, β ∝ ν². Planck PR4 n = -0.35 ± 0.48 disfavours n = 2 at ≈ 4.9σ. None of its 19 citers redoes the fit.
   - Gauge invariance off FRW: a non-closed p needs a gravitational CS term (0907.5159) or a photon mass (2609.10671).
2. **Low-z k_AF bounds:**
   - Best: UV radio galaxies at z = 2.3-3.6 (1003.4823), θ = -0.8° ± 2.2°, giving |(k_AF)_T| < 1.7×10⁻³⁴ eV (p₀ < 3.3×10⁻³⁴ eV). This is ~30 times above the constant k_AF that fits ACT (5.7×10⁻³⁶ eV).
   - For EC-06, ξε = 1.07×10⁻³ gives β(z = 2.8) = 0.041°, while the best low-z bound is ξε < 0.14. Low-z astrophysics does not test EC-06.
   - Local bounds (pulsar RM, geomagnetic spatial components) are irrelevant.
   - GRBs give no d = 3 bound.
3. **Novelty list corrected (`literature/novelty_audit.md`):**
   - KK P′ argument: a known symmetry argument; its application to KK photon birefringence was not found.
   - EC shortfall numbers: new numbers inside known no-gos (2310.09152, Dobado-Maroto).
   - Torsion-pseudoscalar decoupling: the known ALP-EFT result (1708.00443, 1903.12559, 2011.10025, 2112.00553) applied; shown for gluons by 2506.11836.
   - Cassini kill: the known solar-system exclusion of a light KK scalar, newly applied to 2608.09299.
   - KK-05 charge kill: the classic KK charge-to-mass problem (gr-qc/9805018 Sec. 4.2).

## Verdict
- EC-06 cannot be claimed as a new mechanism. Its phenomenology is published in curvature and nonmetricity form.
- A torsion version can add only two things: the origin of the background, and that origin's own bounds (MZ GW speed, BBN on trace torsion, gauge invariance off FRW).
- In the only quantified unsourced axial background (MZ, 2001.06373) the needed ξ is ≳ 10⁴.
- The RESULT draft's novelty wording should follow `novelty_audit.md`.

## Open issues
- The MZ estimate ε ≲ 2×10⁻⁸ uses P ≈ γH/3 (2001.06373 Eq. 25, early tracking) and |S_0̂| = 2|P|. The theorist must redo it in CONVENTIONS.md units and check whether Eq. 68 applies to the late P ∝ a⁻³ branch.
- Equation numbers taken from tool-assisted reading should be opened by whoever uses them: 2001.06373 Eqs. 25, 68, 70; 1908.05184 Eq. 37; hep-ph/0611192 Eq. 40; 2508.10532 Eq. 36.
- CFJ 1990 and Wardle-Perley-Cohen were not opened; their values come via astro-ph/9704263 and the 1011.4865/1409.8242 compilations. The CL is not stated for most low-z rows.
- 2611.23082 uses a k-to-angle normalization different from ours by a factor of about 6 (same issue as `data/bounds.csv` P-kV00-MK). Not resolved here.
- No citable paper was found stating 4D parity = (x⃗, y) → (-x⃗, -y) in KK reductions. A targeted search of 5D orbifold/UED reviews could turn one up.
- Not re-read: whether 2506.11836 states the photon version of its anomaly/mass cancellation.
- Skeptic 020 MINOR 8 (no ULA bound at 5.6×10⁻³³ eV) was not re-searched here.
