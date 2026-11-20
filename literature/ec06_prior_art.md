# Prior art for candidate class EC-06 (unsourced O(H) background torsion + photon coupling)

Task 20261120-024 (literature scout), 2026-11-20. Logs and scripts: `work/20261120-024-lit-ec06-novelty/`.

EC-06 (CLAIMS.md): a homogeneous torsion background of order H that is not sourced by spin (Magueijo-Złośnik varying-Λ EC, teleparallel FRW torsion, PGT torsion cosmologies), plus a non-minimal CFJ-type photon coupling -½ξ p_μ A_ν F̃^{μν} with p built from torsion. With p_0̂ = ξεH, β(z) = (ξε/2) ln(1+z), so ACT β = 0.215° needs ξε = 1.07×10⁻³ (`kaf_lowz.py`, from `data/`). Every finding here is about a candidate mechanism. ACT β is 2.9σ with systematics caveats.

Every arXiv ID below was returned by the arXiv API `id_list` query on 2026-11-20 (title and authors match; `abstracts1-4.log`). Pre-arXiv items were checked on INSPIRE. Reading depth: **[FT]** relevant sections of the arXiv HTML/ar5iv text read with tool-assisted extraction (not line by line); **[ABS]** abstract only; **[LIST]** title in a search listing.

## Verdict

**EC-06 is partly published.** The phenomenology is published, the torsion origin is not.

- **Published:** a CFJ vector that is the gradient of the log of a background geometric scalar, so p_0 ∝ H with no new light field, giving an achromatic β ∝ ln[X(z)/X(0)], fitted to CMB data with a dimensionless coupling of order 10⁻⁴-10⁻². This exists in two versions:
  - Ricci scalar: Li, Xia & Li hep-ph/0611192 (∂_μ ln R, fitted to WMAP3+B03). Background: Li & Zhang 0810.0403, Xia 1201.4457.
  - Nonmetricity scalar: Yu, Yang & Qiu 2508.10532 (ln Q with Q = 6H², fitted to Planck β; α ≈ 10⁻⁴).
  
  In the matter era ∂_t ln R = ∂_t ln Q = -3H. This equals minus the TEGR trace torsion T_0̂ = 3H. So the teleparallel-trace version of EC-06 (p = ξT, β = (3ξ/2)ln(1+z), skeptic 020 S1) is the same model as these two through matter domination. It differs only in the radiation era (R: Ṙ ≈ 0; Q: ∂ ln Q = -4H; TEGR: 3H) and in the Λ era.
- **Published, photons not treated:** the Magueijo-Złośnik parity-odd FRW torsion (1908.05184) and its follow-ups. Nobody couples it to photons. The follow-up 2001.06373 bounds its amplitude through the GW speed. That bound gives ε ≲ 2×10⁻⁸ (scout estimate below), so ξ ≳ 5×10⁴ would be needed.
- **Not found:** any paper coupling photons through a CFJ term (or F∂F̃T, F²T²) to an unsourced cosmological torsion background (MZ, TEGR/f(T), PGT) and computing β. Also not found: any discussion of the gauge invariance of p = ξS off FRW for such a background.
- **Closest torsion-photon fit:** Das-Mohanty-Prasanna 0908.0629. It uses a different operator (dimension 5, β ∝ ν²) on an unspecified constant torsion. Planck PR4 chromaticity disfavours that scaling at ~4.9σ (below).

## 1. Magueijo-Złośnik and follow-ups

INSPIRE record 1749561 has 37 citing papers (`citing.log`), all titles screened. A full-text filter on birefringence, photon or electromagnetic (`citing2.log`) leaves 2001.06373 and 2409.15509 (both GW only), 2211.05097 (quantum cosmology), plus PGT/EC reviews and dark-matter spin papers. **None couples the MZ torsion to the Maxwell field.**

| Paper | Depth | Content relevant to EC-06 |
|---|---|---|
| Magueijo & Złośnik, **1908.05184**, PRD 100 084036 (2019) | [FT] | The action is not plain EC (their Eq. 7). Λ is promoted to a field, and the action contains a Palatini term, a quasi-Euler term ∝ 1/Λ, Nieh-Yan, and a Pontryagin term ∝ 3/(γΛ). FRW torsion is T^i = -T(t)e⁰∧e^i + P(t)ε^i_{jk}e^j∧e^k: T is the trace (parity-even) part, P the axial (parity-odd) part. P ≠ 0 needs Λ̇ ≠ 0 (Eq. 37). For finite γ, P = 0 is not consistent. In vacuum P = c₀/a. In tracking solutions P² scales with ρ, with Ω_P/Ω_matter → 1/2 as γ → ∞ (an extra radiation-like component, which they say conflicts with data) and ≈ γ²/9 for small γ. No Maxwell field; photons enter only as a fluid. They name CMB TB as a possible future probe of the Weyl background. **So the MZ background exists only in an extension of EC (varying Λ + quasi-topological terms), not in EC with constant Λ.** |
| Alexander, Jenks, Jiroušek, Magueijo, Złośnik, **2001.06373**, "Gravity waves in parity-violating Copernican Universes" | [FT] | The same theory. Torsion T̄^I = T e^I e^0 + P ε^{IJK}e_J e_K (Eq. 16). Early tracking P = (γ/3)√(κρ/3) (Eq. 25), i.e. P ≈ γH/3. Late P ∝ a⁻³ (Eq. 30). GW speed: deviation O(γ²), and LIGO gives **γ² ≲ 10⁻¹⁵ (Eq. 68)**. Linear consistency needs \|γ\| ≪ 10⁻²¹ (Eq. 70, not reconciled with their γ = 10⁻⁵ example). GW only: photons, Maxwell and CMB polarization are not analysed. |
| Alexander, Daniel, Magueijo **2207.08885**; Alexandrov, Speziale, Złośnik **2104.03753**; Magueijo-Złośnik-Speziale **2006.05766**; Magueijo **2104.11529**, **2110.05920**, **2306.08390** | [LIST]/[ABS] | Canonical structure, quantum cosmology, varying constants. No photons. |
| Barriga, Izaurieta, Lepe **2409.15509** | [ABS] | GW polarization rotation ("gravitational Faraday") from a dark-matter axial spin tensor. GW only, spin-sourced. |
| Cruz, Izaurieta, Lepe **2005.04550**; Izaurieta et al. 2004.13163, 2004.06058 | [ABS]/[LIST] | Phenomenological FRW torsion ansätze (thermodynamics, H₀ tension). No photons. |

**Scout estimate (unchecked, for the theorist to redo in CONVENTIONS.md units):** T^I_{JK} = 2Pε^I_{JK} gives \|S_0̂\| = 2\|P\| with S^μ = (1/6)ε^{μνρσ}T_{νρσ}. Tracking P ≈ γH/3 then gives ε ≡ \|S_0̂\|/H ≈ 2γ/3. The GW bound γ ≲ 3.2×10⁻⁸ gives ε ≲ 2×10⁻⁸, so ξ ≳ 1.07×10⁻³/2×10⁻⁸ ≈ 5×10⁴ for ACT β. This holds only in the 2001.06373 theory (Pontryagin coupling tied to γ). The late-time P ∝ a⁻³ falls faster than H, which lowers β further.

## 2. Teleparallel, f(T), Nieh-Yan modified teleparallel

INSPIRE full text "cosmic birefringence" AND teleparallel: 51 hits; AND "Nieh-Yan": 20 (`inspire_search.log`). All titles screened. Citations of 2007.08038 (recid 1807227) and 2104.05917 (recid 1858057) filtered on "photon": 17 and 10 (`citing2.log`).

| Paper | Depth | Photons? |
|---|---|---|
| Li, Rao, Zhao **2007.08038** (JCAP 11 (2020) 023) | [ABS] | No. Nieh-Yan-modified teleparallel gravity with a scalar, GW velocity birefringence only. |
| Li, Rao, Tong **2104.05917** | [ABS] | No. Same model with local Lorentz covariance restored, general FRW. |
| Wu, Zhu, Niu **2110.13870** | [ABS] | No. LIGO/Virgo bound on the Nieh-Yan coupling (GW). |
| Hohmann & Pfeifer **2203.01856** | [ABS] | No. GW birefringence in spatially curved NGR cosmology (needs spatial curvature). |
| Hohmann & Pfeifer **2012.14423**, "Teleparallel axions and cosmology" | [ABS] | No photon coupling. Notes cosmologically symmetric teleparallel geometries with nonzero axial torsion (spatially curved branch). Relevant as a source of an axial background, not as a photon paper. |
| Rao, Zheng, Hou **2608.09299** | (already [KEY]) | Yes, but 5D KK: radion × FF̃ (KK-04). The only teleparallel paper with a photon parity-odd term. |
| Xu & Zhang 2610.29182; Zhang-Feng-Gao 2404.02922, 2307.00330; Zhai et al. 2508.06984; Cai et al. 2112.04794 | [LIST]/[ABS] | No (GW / scalar-induced GW; one about chiral GW imprinting on CMB). |
| f(T) + "cosmic birefringence" full-text query | [LIST] | 105 hits, all reviews or unrelated (`citing2.log` section). No f(T) photon paper found. |
| Landry 2512.20337, 2605.25752, 2606.16982 | [LIST] | Maxwell field as a source in f(T) solutions. Minimal, parity-even. |
| Doyran, Sert, Adak **2311.11972** | [ABS] | Non-minimal EM coupling in symmetric teleparallel (nonmetricity) gravity, static solutions. No birefringence. |

The Nieh-Yan teleparallel line is GW-only in 4D. None of these papers writes a CFJ or φFF̃ term built from the teleparallel torsion vector.

## 3. PGT / Riemann-Cartan torsion with photon couplings

| Paper | Depth | Content |
|---|---|---|
| Bahamonde, Maggiolo, Pfeifer **2507.02362** | [FT] | General linear gauge-invariant F-torsion couplings (Eq. 17: k₃F^{μν}R̃_{[μν]} + k₄₋₆ F·T·T). Eq. 23 contains ε-contracted F·S·T pieces, which are dropped in the reduced model. No FRW solution, no bound. Conclusions: the mechanism "could offer a natural explanation" of cosmic birefringence, as future work. **Not a computation of β.** |
| Trukhanova & Obukhov **2410.01355** | [ABS] | Quadratic PGT. Axial torsion sourced by EM helicity. Rotation only in external E/B fields. |
| Kranas, Tsagas, Barrow **1809.10064** | (bounds.md) | FRW trace torsion φ ∝ H, BBN bound -0.0058 < φ/H < 0.0194. In our notation \|T_0̂\| = 6\|φ\| ≲ 0.12H at BBN. Bounds a trace-vector background (ε′ ≲ 0.1), not a photon coupling. Not binding for ξ′ε′ ~ 10⁻³. |
| Boudet, Bombacigno, Moretti **2209.14394** | [ABS] | Metric-affine CS gravity: FRW torsion without spin (sourced by a parity-odd term); GW birefringence only. |
| Obukhov PGT primer **2206.05205**; Ni review **0912.5057** | [LIST]/[ABS] | Reviews. Ni discusses photon pseudoscalar coupling and polarization rotation tests generically. |
| Silva **2605.29369** | [ABS] | CFJ-type coupling to "axial contortion" of a chiral optical medium (condensed-matter analogue). Not cosmology. |

## 4. Das-Mohanty-Prasanna 0908.0629 and its citers

[FT] (ar5iv), IJMPD 22 (2013) 1350011 (DOI 10.1142/S0218271813500119, from `prior_work_torsion_birefringence.md`).
- Coupling (Eq. 4): L = -¼F² + ξ₁T^{αλ}_ρ F_{αν}∂_λF̃^{ρν} + ξ₂T^{σγ}_δ F_{σν}∂_γF^{δν}. Gauge invariant (built from F), with one extra derivative relative to CFJ.
- Background: a "background torsion" with no source named (not spin, not FRW-derived). Only the combination T₁ enters. It is defined inconsistently: T⁰⁰₃ + T³³₀ after Eq. (12), T⁰³₀ - T³³₀ after Eq. (20). Neither is an FRW-allowed component (CONVENTIONS: only axial S^0 and trace T_{i0j} ∝ δ_ij). T₁ is taken constant along the path, and its time dependence is not discussed.
- Dispersion ω± = p(1 ± ξ₁pT₁). Rotation α = 2ξ₁p²T₁t (Eq. 20), so **β ∝ ν²**.
- Fit: WMAP+BOOMERanG α = -2.4° ± 1.9° at 100 GHz gives ξ₁T₁ = (-3.35 ± 2.65)×10⁻²² GeV⁻¹ (Eq. 43), with H₀ = 72, z = 1100, Ω_m = 0.3. ξ₂ does not affect propagation.
- **Against current data:** β ∝ ν² is a power-law index n = 2. Planck PR4 gives n = -0.35 (+0.47/-0.48) (`data/birefringence_measurements.csv` row ESK22-n-all, 2201.13347), so n = 2 is ≈ 4.9σ away. ACT per-band ψ prefer n = 1.61 at 2.0σ over achromatic (skeptic 020 item 10), which does not rescue it either. The DMP operator is disfavoured for an isotropic signal unless T₁ is direction dependent (then it is an anisotropic effect).

INSPIRE: 19 citing papers (`citing.log`). None redoes the fit with Planck/ACT or names a torsion source. Kanazawa, Yoshioka, Lambiase, Vilasi, IJGMMP 15 (2018) 1850176 (INSPIRE 1698169, no arXiv) [ABS]: photon rotation around a scalar-tensor "screwed" cosmic string with torsion. Wavelength- and redshift-dependent rotation, local not isotropic. The rest are leptogenesis in CPT-violating KR backgrounds (Mavromatos et al.), CMB parity reviews, and Fedderke-Graham-Rajendran 1903.02666 (axion DM).

## 5. CFJ vector tied to cosmic expansion (any origin)

| Paper | Depth | Content |
|---|---|---|
| **Li, Xia & Li, hep-ph/0611192** | [FT] | L_int = -c ∂_μ ln R J^μ_{(B-L)L} (Eq. 2). Through the anomaly this becomes (2cα/3π)∂_μ ln R A_ν F̃^{μν} (Eq. 14). Rotation Δα = (δ/2) ln[R(0)/R(z)] (Eq. 40), δ = -4αc/3π (Eq. 17), achromatic. WMAP3+B03: δ = -0.011, consistent with 0 at 2σ (\|δ\| < 10⁻²). No new light field: the vector is ∂ ln R with dark energy = Λ. The gauge-invariance question does not arise because p is an exact gradient. **Same phenomenology as EC-06 with p ∝ H.** |
| Li & Zhang **0810.0403** | [ABS] | General-relativistic rotation formula for a derivative coupling of a dark-energy scalar or f(R) to currents. Includes spatial fluctuations. |
| Xia **1201.4457**; Xia-Li-Zhang 0908.1876, 0803.2350 | [ABS]/[LIST] | CMB fits of the CS rotation angle in this framework. |
| Davoudiasl, Kitano, Kribs **hep-ph/0403019** | [ABS] | Origin of the ∂_μR J^μ operator (gravitational baryogenesis). |
| **Yu, Yang & Qiu, 2508.10532** | [FT] | (α/M²) l(Q) F̃F (Eq. 14) ≡ 2 l_Q ∂_μQ F̃^{μν}A_ν (Eq. 17). Δχ = (2α/M²)[l(Q_rec) - l(Q_0)] (Eq. 36), Q = 6H². Planck β = 0.30° ± 0.11°: l = Q needs α ~ 10¹⁰⁹ (excluded as unnatural); **l = M² ln(Q/M_P²) gives Δχ ≈ 40.2α and 8.2×10⁻⁵ < α < 1.8×10⁻⁴**. The log case is the EC-06 phenomenology in symmetric teleparallel gravity (p ∝ ∂ ln H²). Torsion is mentioned only as a possible generalization. |
| Li, Cai, Wang **0907.5159** | [ABS] | A fixed (non-gradient) CFJ p_μ makes the minimally coupled stress tensor non-conserved. Consistency needs a gravitational CS term with a coefficient fixed by p. **Directly relevant to EC-06 off FRW:** if p is not closed, the gravity side must be modified too. |
| Hook, Huang, Zhou **2609.10671**, "Ultralight axial dark matter" | [FT] | -(g/4)V_μK^μ with K^μ = ε^{μνρσ}A_νF_{ρσ} (Eqs. 1.5, 3.2) for a dynamical axial vector. "Not gauge invariant, so for consistency one must include a photon mass" (m_γ ~ 10⁻¹⁵ eV benchmark). Rotation dθ/dz = (g/2)(V⁰ - V^z) (Eq. 4.2), CFJ-like. They treat oscillating DM only, not a static V⁰ ~ H. **Relevant to EC-06's gauge-invariance question:** an axial torsion vector that is not a gradient off FRW is in this class. |
| Galaverni, Finelli, Paoletti **2301.07971** | [ABS] | β(z) histories (phenomenological and EDE/quintessence) produce ℓ-dependent spectra that are not degenerate with miscalibration. Applies to the ln(1+z) history of EC-06. |
| Carroll **astro-ph/9806099** | [ABS] | Quintessence φ̇ ~ HM_P with φFF̃: the original "rotation from a field rolling on a Hubble time" plus the fifth-force argument. |

No paper was found that names torsion as the origin of a Hubble-tracking CFJ vector. INSPIRE title queries "Chern-Simons" + "Lorentz violation" + cosmolog* and "Lorentz violation" + "cosmic birefringence" return 0. Full-text CFJ + Hubble + torsion gives 8 hits, none relevant beyond 2609.10671 (`inspire_search.log`).

## 6. Low-z bounds on a constant or slowly varying k_AF

Table in `data/kaf_lowz_bounds.csv`, derived columns from `work/20261120-024-lit-ec06-novelty/kaf_lowz.py` (PASS).

| row | z | bound | \|(k_AF)_T\| (eV) | EC-06 ξε |
|---|---|---|---|---|
| L-CFJ90 (CFJ 1990 via astro-ph/9704263 Eq. 5) | 0.4-1.5 | p₀ ≤ 1.7×10⁻⁴² h₀ GeV | < 5.7×10⁻³⁴ | - |
| L-CF97 (astro-ph/9704263 Eq. 7) | 0.3-2.0 | p₀ = (1.25 ± 1.71)×10⁻⁴² h₀ GeV | < 1.6×10⁻³³ | - |
| **L-AFG10-UV (1003.4823)** | 2.34-3.57 | θ = -0.8° ± 2.2° | **< 1.7×10⁻³⁴** | < 0.14 |
| L-Leahy (via 1411.6287 Tab. 3) | 0.43-0.97 | 1.6° ± 1.8° | < 3.8×10⁻³⁴ | < 0.47 |
| L-WPC97 (via 1409.8242 Tab. 1) | 0.811 | 1.4° ± 1.1° | < 1.8×10⁻³⁴ | < 0.21 |
| L-pulsarRM25 (2503.20626) | kpc | < 10⁻³⁶ GeV | < 10⁻²⁷ | - |
| L-geomag26 (2602.14217) | lab | spatial \|k_Z\| < 4×10⁻²⁵ GeV | spatial only | - |

- The best low-z bound on a constant isotropic CFJ coefficient is UV radio-galaxy polarization at z ≈ 2.8, \|(k_AF)_T\| < 1.7×10⁻³⁴ eV (p₀ < 3.3×10⁻³⁴ eV). It is ~30 times weaker than the value a constant k_AF needs for ACT β (5.7×10⁻³⁶ eV).
- For EC-06 (β ∝ ln(1+z)), ξε = 1.07×10⁻³ predicts 0.012° at z = 0.47, 0.018° at z = 0.81 and 0.041° at z = 2.8. Low-z data give ξε < 0.14 at best, two orders above the needed value. **Low-z astrophysics does not test EC-06.**
- GRB polarization (Kostelecky-Mewes 1301.5367, 0809.2846) bounds only operators with d ≥ 4. A d = 3 (CFJ) rotation is energy independent, so it does not depolarize a broadband signal and gives no bound. Not tabulated.
- CMB reference: Lembo, Caloni, Gerbino 2611.23082 (17 Nov 2026) give \|k^(3)_(V)00\| < 7.75×10⁻⁴⁵ GeV (95%, from BB). Their normalization differs from the KM Eq. 158 map by ~6 (see the CSV note and `data/bounds.csv` P-kV00-MK). By their own relations the bound corresponds to α ≲ 0.54°, so it does not contradict ACT.

## 7. What this means for EC-06 (scout's reading; for the theorist and cosmologist)

1. The p ∝ H CFJ phenomenology with ξε ~ 10⁻³ is published twice (∂ ln R, ∂ ln Q). An EC-06 result can be new only in (a) the torsion origin of the background and (b) the bounds specific to that origin: the background's own field equations, the GW-speed bound in the MZ theory, BBN on trace torsion, and the matter couplings of the torsion.
2. MZ torsion lives in an extension (varying Λ + Euler/Pontryagin terms), and the only published amplitude bound (2001.06373) suppresses ε to ≲ 10⁻⁸. That pushes ξ to ≳ 10⁴-10⁵.
3. For TEGR the trace vector is exactly T_μ = ∂_μ ln a³ on FRW in the diagonal-tetrad gauge. That makes p = ξT gauge-safe on the background and identical in matter domination to Li-Xia-Li/Yu-Yang-Qiu. Off FRW, T_μ is not closed in general. 0907.5159 and 2609.10671 show what a non-closed p costs: a gravitational CS counterterm, or a photon mass.
4. The ℓ-dependence of a ln(1+z) history (reionization vs recombination split) separates it from miscalibration (2301.07971). The cosmologist's task 2 should use this.
