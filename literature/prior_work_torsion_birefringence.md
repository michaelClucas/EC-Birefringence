# Prior work: torsion, photons, and birefringence

Task 20261120-002 (literature scout, theory side). Compiled 2026-11-20.

Every arXiv ID below was checked against the arXiv API (title, authors, journal ref) on 2026-11-20. Pre-arXiv papers were checked against INSPIRE-HEP by DOI. Nothing here is cited from memory without that check.

Formulas are written in each paper's own conventions. They are NOT converted to `CONVENTIONS.md` yet. Anyone who uses one must convert it and show the check.

Reading depth tags:
- **[KEY]**: the action, the rotation formula and the numbers were read from the full text (arXiv HTML or ar5iv rendering, read with tool-assisted extraction of the relevant sections). Not a line-by-line read.
- **[ABS]**: abstract and metadata only.

Flags:
- **SAME-AIM**: the paper already tries what this program does (torsion or a torsion extension, giving an isotropic β compared with CMB data).
- **NO-GO**: the paper gives a no-go, or a bound that kills a route.

Summary of the search: no paper found derives β from torsion in Einstein-Cartan (or an EC extension) and fits it to the Planck PR4 or ACT DR6 β with lab torsion bounds applied. The nearest papers are Das-Mohanty-Prasanna 2009 (fits a non-minimal torsion coupling to pre-Planck CMB rotation), Rao et al. 2026 (5D teleparallel KK, compared with Planck PR4 β), Maity-SenGupta-Sur 2005 (KK-reduced Kalb-Ramond axion, compared with radio-galaxy rotation), and the neutrino-current papers together with Reig 2026's no-go for them. The neutrino-current route is the effective operator that a torsion-mediated fermion-current route reduces to.

---

## A. Torsion-photon couplings and torsion-induced birefringence

### A1. Non-minimal, gauge-invariant torsion-Maxwell couplings (classical)

| Paper | Coupling / action | Result | Relevance |
|---|---|---|---|
| Itin & Hehl, gr-qc/0307063, PRD 68 (2003) 127701 **[KEY]** | Lists all lowest-order gauge-invariant couplings of F to torsion. The F·T² class changes the Maxwell equations. The F²·T² class (17 terms, length ℓ²) only changes the constitutive tensor. Preuss's ℓ² (T_a∧F)∧⋆(T^a∧F) is one member. | The totally antisymmetric part of the constitutive tensor is an induced axion θ = ℓ² × (quadratic in torsion). In FRW torsion with trace part T_{0αβ} = (u/ℓ)δ and axial part T_{αβγ} = (v/ℓ)ε, they get **θ = c₁ u v**, so the optical activity goes as d(uv)/dt. FRW gives no birefringence (single light cone) but a torsion-dependent speed of light. No rotation formula and no numbers. | Very relevant to the 4D track. It is a parity-even, gauge-invariant route to an isotropic rotation, and it needs **both** FRW torsion components (trace and axial). Minimal EC with Dirac spin sources produces only the axial part, so θ = 0 there unless something sources the trace part. No one has put numbers on this. |
| Rubilar, Obukhov & Hehl, gr-qc/0305049, CQG 20 (2003) L185 **[ABS]** | Preuss's ℓ² T² F² coupling | O(3)-symmetric torsion gives vacuum birefringence. Isotropic FRW torsion changes the speed of light. | Parity-even light-cone effects, not a polarization-plane rotation. Background reading. |
| Preuss, Haugan, Solanki & Jordan, gr-qc/0405068, PRD 70 (2004) 067101 **[ABS]** | Same Preuss coupling | Polarimetry of the magnetic white dwarf RE J0317-853 gives k² ≲ (19 m)² | A bound on ℓ² for the Itin-Hehl class. Bounds analyst should record it. |
| Preuss, Solanki, Haugan & Jordan, gr-qc/0507071, PRD 72 (2005) 042001 **[ABS]** | Gauge-invariant non-minimal torsion-Maxwell coupling in Poincaré gauge theory | Explicit phase shift between polarization modes. White-dwarf bound on the coupling. | The PGT version of the same. |
| **Das, Mohanty & Prasanna, 0908.0629**, IJMPD (2013), DOI 10.1142/S0218271813500119 **[KEY]** **SAME-AIM (partial)** | L = -¼F² + ξ₁ T^{αλ}_ρ F_{αν} ∂_λ F̃^{ρν} + ξ₂ T^{σγ}_δ F_{σν} ∂_γ F^{δν}. Gauge invariant because it is built from F only. Dimension 5 in the torsion-coupling sense. | Circular modes ω± = p(1 ± ξ₁ p T₁), so the rotation is **α = 2 ξ₁ p² T₁ t**, which scales as **frequency squared**. Torsion is taken constant over the light path. Using WMAP+BOOMERanG α = (-2.4 ± 1.9)° at 100 GHz gives ξ₁T₁ = (-3.35 ± 2.65)×10⁻²² GeV⁻¹. ξ₂ has no effect on propagation. | The only paper found that fits a torsion-photon coupling to CMB rotation. Problems: T₁ is defined inconsistently in two places (T^{03}_0 - T^{33}_0 versus T^{00}_3 + T^{33}_0). Neither is obviously an FRW-allowed component. The units of the bound are unclear (p is plugged in as "100 GHz"). It predates Planck. **The ν² scaling should be confronted with the multi-frequency β consistency in ACT/Planck.** |
| Hehl & Obukhov, gr-qc/0001010, Lect. Notes Phys. 562 (2001) 479 **[ABS]** | Premetric electrodynamics. Gravity (including torsion) enters only through the constitutive law. | Lists the possible non-minimal coupling schemes and the Ni "abelian axion" piece. | Framework for D, and for why minimal coupling sees no torsion. |
| De Felice, Heisenberg, Olmo & Pastor-Marcos, 2607.26141 (2026) **[ABS]** | Metric-affine EM sectors that keep U(1) and projective invariance: a scalar prefactor built from non-Riemannian invariants, or a rank-4 constitutive tensor | Geometric-optics limit. Polarization-dependent propagation, birefringence and mode mixing. No CMB application. | The most recent general framework. A candidate source for the most general gauge-invariant torsion-photon operator list. |
| Kruglov, 0804.4011, EPJC 63 (2009) 343 **[ABS]** | Propagating axial torsion (from higher-derivative gravity) mixing with photons in an external B field | Ellipticity and rotation in a lab B field | Lab or astrophysical B-field effect. Not isotropic cosmic β. |

### A2. Fermion loops and Carroll-Field-Jackiw (CFJ) from axial torsion

| Paper | Coupling / action | Result | Relevance |
|---|---|---|---|
| Carroll, Field & Jackiw, PRD 41 (1990) 1231, DOI 10.1103/PhysRevD.41.1231 (INSPIRE-verified, pre-arXiv) **[ABS]** | ½ p_μ A_ν F̃^{μν} with a constant p_μ | Defines the CFJ term and the original radio-galaxy limit | The baseline in `CONVENTIONS.md`. |
| **Dobado & Maroto, astro-ph/9706044**, MPLA 12 (1997) 3003 **[KEY]** | Charged fermions with minimal torsion coupling -(1/8)ψ̄ γ₅ S̸ ψ, where S_β = ε_{μναβ}T^{μνα}. Integrating out the fermions (dim. reg.) gives **L_eff ⊃ -(α/24π)(Σ Q_i²) ε^{μναβ} F_{μν} A_α S_β**, with ΣQ² = 8 for 3 SM generations. | A CFJ term with p_μ ∝ (α/6π)·8·S_μ. To explain the (later discredited) Nodland-Ralston anisotropy they need \|S\| ~ 10⁻³⁰ eV. They note this would need massless or condensed torsion, not Planck-mass torsion. | **The basic "torsion → fermion loop → CFJ" calculation.** Our program's indirect route is this one. The coefficient is regularization dependent (see Jackiw-Kostelecky). The authors pick dim. reg. and do not discuss the ambiguity. |
| Jackiw & Kostelecky, hep-ph/9901358, PRL 82 (1999) 3572 **[ABS]** | A constant axial vector b_μ in the fermion sector (b_μ ψ̄γ^μγ₅ψ, which is exactly how axial torsion enters) | The induced CFJ coefficient k_AF ∝ b_μ is **regularization dependent**: definite and nonzero non-perturbatively, undetermined perturbatively. | **Key caveat for any torsion → CFJ claim.** The coefficient is not fixed by the low-energy theory. |
| Alexander, Ochoa & Kosowsky, 0810.2355, PRD 79 (2009) 063524 **[KEY]** | g ε^{μναβ} A_μ T_ν F_{αβ}, with T_ν a torsion-related external vector. Motivated as a radiative correction from torsion-fermion coupling (citing de Sabbata-Gasperini). | Boltzmann transport: rotation rate ω_QU = g·T. Circular polarization is down by ω_QU/k. They only require gT ≪ H₀ and give no number. | CMB transport formalism for a torsion-CFJ term. No fit. |
| Garcia de Andrade, hep-th/0111251 (2001) **[ABS]** | Identifies the CFJ vector with the axial torsion vector | Converts CFJ limits into 10⁻³³ eV < \|S^μ\| < 10⁻³² eV | Earliest "CFJ bound means torsion bound" mapping. The O(1) proportionality is assumed, not derived. |
| Garcia de Andrade, hep-th/0110150 (2001) **[ABS]** | Proca field minimally coupled to torsion | Reproduces a Nodland-Ralston-type term, with a torsion-generated photon mass. Primordial torsion ~10⁻²⁹ eV from photon-mass limits. | Breaks gauge invariance through the photon mass. Of historical interest only. |
| Garcia de Andrade, 1105.1573, MPLA 27 (2012) 1250059 **[ABS]** | de Sabbata-Gasperini perturbative modified Maxwell equations with torsion | Torsion bounds from Faraday rotation and the CMB | Low weight. Same author series as 2404.13517, 2412.05982, 2412.16617 (Gao & Garcia de Andrade: torsion, dark photons, Barbero-Immirzi, Chern-Simons electrodynamics; all [ABS]). |
| Soo, hep-th/9805090, PRD 59 (1999) 045006 **[ABS]** | Chiral fermions on Riemann-Cartan space, Pauli-Villars | The Nieh-Yan form contributes to the ABJ anomaly with a coefficient ∝ (regulator mass)² | The anomaly route from axial torsion. The Λ² coefficient is the reason for the long-running Nieh-Yan anomaly debate. |
| Chandia & Zanelli, hep-th/9702025, PRD 55 (1997) 7580 **[ABS]** | Nieh-Yan form N = T^a∧T_a - R_{ab}∧e^a∧e^b | N is a topological invariant on spaces with torsion and appears in the chiral anomaly | Background for Nieh-Yan as a θ-like term. |
| Nieh & Yan, J. Math. Phys. 23 (1982) 373, DOI 10.1063/1.525379 (INSPIRE, pre-arXiv) **[ABS]** | The Nieh-Yan identity | | Citation for the Nieh-Yan term. |

### A3. Lorentz-violation (SME) bounds on torsion

| Paper | Content | Relevance |
|---|---|---|
| **Kostelecky, Russell & Tasson, 0712.4393**, PRL 100 (2008) 111102 **[ABS]** | Maps minimal torsion-fermion couplings onto SME b_μ-type coefficients. Constrains 19 of the 24 torsion components to ≲ 10⁻³¹ GeV. | **The main lab bound on background torsion.** Bounds analyst must record the table and the exact components. Any S₀ that makes β through fermions must sit under this. |
| Russell, 0803.2228 (2008 proceedings) **[ABS]** | Dual-maser bounds on axial and mixed-symmetry torsion | Companion to the above |
| Heckel et al., 0808.2673, PRD 78 (2008) 092006 **[ABS]** | Spin-polarized torsion pendulum. Preferred-frame and CP tests with polarized electrons. | Source of the electron b_μ bounds used for torsion |
| Kostelecky & Mewes, 0809.2846, ApJ 689 (2008) L1 **[ABS]** | SME photon sector. Dimension-3 (k_AF, i.e. CFJ) bounds from WMAP. | The photon-side counterpart: limits on an isotropic CFJ coefficient from the CMB |
| Ni, 0912.5057, Rep. Prog. Phys. 73 (2010) 056901 **[ABS]** | Review of spin in gravity: torsion, PGT, teleparallel, metric-affine, pseudoscalar-photon (Ni's φ F F̃), CMB polarization rotation tests | Broad review. Good for the cross-check list. |
| Shapiro, hep-th/0103093, Phys. Rept. 357 (2002) 113 **[ABS]** | Review of torsion: matter couplings, anomalies, bounds, propagating antisymmetric torsion | Argues that consistent propagating axial torsion must be heavier than the heaviest fermion it couples to. This works against a light torsion mode that could carry β. |
| Carroll & Field, gr-qc/9403058, PRD 50 (1994) 3867 **[ABS]** | Propagating-torsion theories. Only a massive scalar or vector torsion mode couples to matter. | Torsion decays fast and produces no long-range field. This works against a cosmological torsion background in propagating-torsion models. |

### A4. Holst term, Barbero-Immirzi (BI) field and Nieh-Yan pseudoscalar

| Paper | Coupling / action | Result | Relevance |
|---|---|---|---|
| Perez & Rovelli, gr-qc/0505081, PRD 73 (2006) 044013 **[ABS]** | Holst action with fermions | γ_BI enters the torsion-induced four-fermion coupling | A constant BI parameter changes only contact interactions. No photon effect. |
| Freidel, Minic & Takeuchi, hep-th/0507253, PRD 72 (2005) 104002 **[ABS]** | Holst plus non-minimal fermion-torsion coupling | Effective four-fermion Lagrangian with two parameters. Parity violation (axial-vector term) only with non-minimal coupling. | The parity-odd fermion structure that a fermion-loop route would start from. |
| Taveras & Yunes, 0807.2652, PRD 78 (2008) 064070 **[ABS]** | BI parameter promoted to a scalar field in the Holst action | Torsion is sourced by ∂(BI). The effect equals a scalar field with a non-canonical kinetic term (k-inflation-like). | Origin of the BI-field idea. No photon coupling at tree level. |
| Calcagni & Mercuri, 0902.0957, PRD 79 (2009) 084004 **[ABS]** | BI field coupled to the Nieh-Yan invariant | BI field becomes a canonical pseudoscalar, minimally coupled to gravity | Cleaner than Holst. The BI field is an axion candidate. |
| Mercuri, 0903.2270 (2009); Mercuri, 0902.2764, PRL 103 (2009) 081302 **[ABS]** | BI as a θ-angle of the Nieh-Yan term. Peccei-Quinn mechanism in gravity. | | Sets up "BI field = axion" |
| Lattanzi & Mercuri, 0911.2698, PRD 81 (2010) 125015 **[ABS]** | BI field identified with the QCD axion through the Nieh-Yan-modified action | PQ scale fixed by the model. Dark-matter abundance. Strong bound on r. | If the BI/torsion axion is the QCD axion, its mass is QCD-generated (see next row). |
| **Castillo-Felisola, Corral, Kovalenko, Schmidt & Lyubovitskij, 1502.03694**, PRD 91 (2015) 085017 **[KEY]** | KR axion (Lagrange multiplier for zero Nieh-Yan density) and BI axion, shown to be equivalent. Universal torsion-fermion coupling through the axial current, plus the anomaly. | Effective couplings: **Φ F∧F with coefficient α Q̄²/(π f_Φ)**, Φ G G̃ with N_q/f_Φ, a gravitational term, and a contact (J⁵)². **f_Φ = √(8/3) M_P ≈ 4×10¹⁸ GeV**, fixed with no free parameter. QCD mass **m ~ 10⁻¹² eV**. Isocurvature gives H_I ≤ 10¹⁰ GeV if it is the dark matter. | **Important for the 4D track (candidate NO-GO, scout reading, unchecked):** torsion couples universally to quarks, so the gluon anomaly term is generic. The torsion axion then gets a QCD mass ≫ H at recombination, oscillates early, and its field excursion between last scattering and today is far too small for β (see washout, 2512.21888). Escaping this needs an extra (non-QCD) potential or a cancelling PQ structure. The team should verify the chain: it would kill the minimal "torsion axion" route quickly. |
| Baekler & Hehl, 1105.3504, CQG 28 (2011) 215017 **[ABS]** | Most general quadratic torsion and curvature Lagrangian with both parities, including the Holst, Nieh-Yan and Pontryagin boundary terms | The BI parameter can only be discussed properly with torsion-squared pieces included | Reference action set for "EC extensions" |
| Bhattacharjee & Chatterjee, 1101.0118, PRD 83 (2011) 106007 **[ABS]** | Gauge-invariant coupling of form fields to torsion (string-inspired). The BI field as a torsion source. | Polarization rotation of EM and GW waves | Background |
| Adshead, Brahma & Das, 2507.10696, JCAP 03 (2026) 057 **[ABS]** | Pseudoscalar coupled to the Pontryagin and Nieh-Yan terms in Einstein-Cartan-Palatini | The Nieh-Yan coupling raises the effective decay constant. Parity-violating torsion during inflation. Chiral GWs. | Shows that a pseudoscalar times Nieh-Yan in EC changes the axion kinetic term. Relevant if β comes from such a field with an anomaly coupling to photons. |
| Mir, Shabir, Tiwari & Faizal, 2604.03973, Fortsch. Phys. (2026) **[KEY, partly; 75 pp, only sections skimmed]** **SAME-AIM (schematic only)** | BI parameter as a pseudoscalar Goldstone Φ with (∂Φ)J_{B-L}/M*, Φ R R̃, and a (c_A/f_Φ) Φ F F̃ term "fixed by anomaly matching". EC torsion is eliminated algebraically into a four-fermion term. | Claims a relation tying η_B, tensor chirality and cosmic birefringence Δα. In the sections read, the photon coupling is not carried into the cosmology, there is no explicit Δα formula, and there is no fit to Planck or ACT. | Mentions torsion and birefringence together, but is not a torsion β prediction. Low weight. Listed so nobody claims "first BI-birefringence link". |
| Pollari, 2308.03145 (2023); Sur, Dutta & Chothe, 2007.04906, EPJC 81 (2021) 306 **[ABS]** | BI-field cosmology (inflation and quintom) | | Background for BI-field dynamics. No photons. |

### A5. Poincaré gauge, teleparallel and metric-affine birefringence (mostly gravitational waves)

| Paper | Content | Relevance |
|---|---|---|
| Boudet, Bombacigno, Moretti & Olmo, 2209.14394, JCAP 01 (2023) 026 **[ABS]** | Metric-affine Chern-Simons gravity. FRW torsion and nonmetricity solutions. GW birefringence from tensor-torsion coupling. | GW birefringence, not photon. Shows FRW torsion can be sourced by a parity-odd term without spin matter. |
| Wu et al., 2110.13870, PRD 105 (2022) 024035 **[ABS]** | Nieh-Yan modified teleparallel gravity: GW parity-violation constraints from LIGO/Virgo | GW-side bound on the Nieh-Yan coupling. Pairs with Rao et al. (KK) below. |
| Magueijo & Złośnik, 1908.05184, PRD 100 (2019) 084036 **[ABS]** | EC with varying Λ: parity-violating FRW solutions with torsion even without spinors | **Homogeneous isotropic FRW can carry parity-odd torsion (S₀ ≠ 0) without spin matter** in these extensions. This is a source of cosmological S₀ that avoids the tiny spin density, so it is relevant to whether S₀ can be large. |
| Magueijo, Zlosnik & Kibble, 1212.0585, PRD 87 (2013) 063504 **[ABS]** | Einstein-Cartan-Holst with chiral spinors: torsion becomes a 4-fermion self-interaction, possibly parity-violating | Background |
| Helfer, 1408.2729, PRD 94 (2016) 124011 **[ABS]** | Without torsion: high-frequency GWs rotate in a chiral fermion density by 2πN₅ℓ_P². Leads to EB from primordial tensors. | Gravitational analogue of the axial-current route. Planck-suppressed. |

### A6. Old cosmic-anisotropy (Nodland-Ralston era)

- Dobado & Maroto astro-ph/9706044 and Garcia de Andrade hep-th/0110150 tried to explain the Nodland-Ralston radio-galaxy anisotropy with torsion. That signal was later refuted (not re-checked here). These papers are of historical interest only, but the Dobado-Maroto coefficient is useful.

---

## B. Einstein-Cartan cosmology, spin fluids, four-fermion interaction

| Paper | Content | Relevance |
|---|---|---|
| Hehl, von der Heyde, Kerlick & Nester, RMP 48 (1976) 393, DOI 10.1103/RevModPhys.48.393 (INSPIRE) **[ABS]** | Foundational review of EC (U₄ theory) | Standard reference for the algebraic torsion-spin relation and the "spin density matters only at extreme densities" argument |
| Boehmer & Bronowski, gr-qc/0601089, Ukr. J. Phys. 55 (2010) 607 **[ABS]** | **Weyssenhoff fluids (Frenkel, Papapetrou-Corinaldesi, Tulczyjew conditions) are incompatible with the cosmological principle.** FRW allows only the trace vector and axial vector torsion parts. | **NO-GO for the classical spin-fluid source of homogeneous torsion.** The averaged-spin-fluid cosmologies need re-examination. Matches the FRW component statement in `CONVENTIONS.md`. |
| Brechet, Hobson & Lasenby, 0807.2523, CQG 25 (2008) 245016 **[ABS]** | 1+3 covariant dynamics of an irrotational Weyssenhoff fluid. Bounce if spin density exceeds shear. Spin terms negligible at late times. | Spin effects scale as a⁻⁶ and are negligible today |
| Poplawski, 1007.0587, PLB 694 (2010) 181 **[ABS]** | ECKS with spin fluid of quarks and leptons: bounce, then flatness and horizon solved | Torsion matters only at ~10⁵² K-scale densities |
| Poplawski, 1102.5667, GRG 44 (2012) 491 **[ABS]** | Torsion-induced four-fermion interaction as dark energy (quark condensate), with parity-odd pseudoscalar curvature and non-minimal terms | One of the few "late-time torsion" proposals: a condensate gives a constant contact energy |
| Poplawski, 1108.6100 (2011) **[ABS]** | Maxwell field minimally coupled to the affine connection in ECKS. **Breaks U(1)** and gives an effective monopole density. Significant only at extreme energies. | Explicit example of what goes wrong when F is built with the full connection (see D) |
| Khriplovich & Rudenko, 1210.7306, JCAP 11 (2012) 040 **[ABS]** | Axial-axial gravitational four-fermion interaction conflicts with flat or closed FRW. A vector-vector piece could fix it. | Constraint on the EC four-fermion term in the early universe |
| Lucat & Prokopec, 1512.06074, JCAP 10 (2017) 047 **[ABS]** | Generalized EC with arbitrary spin-torsion strength. One-loop renormalized four-fermion term. Bounce. | Torsion couplings to fermions can be non-universal (a free parameter). This matters for the size of S₀. |
| Shaposhnikov, Shkerin, Timiryasov & Zell, 2008.11686, PRL 126 (2021) 161301 **[ABS]** | EC four-fermion and scalar-fermion contact terms as a dark-matter production portal | Modern EC-as-EFT framing. Non-minimal coefficients are free. |
| Mudrunka & Nakayama, 2610.26237 (2026) **[ABS]** | EC four-fermion dark-matter production plus graviton s-channel, with the Holst term | Up-to-date EC contact-term treatment |
| Chakraborty & Lahiri, 2405.20373 (2024) **[ABS]** | Non-universal torsion-fermion couplings give a parity-violating current-current interaction tested in e⁻e⁻ and e⁻D scattering | Collider and atomic bounds on non-universal axial couplings |
| Pasmatsiou, Tsagas & Barrow, 1611.07878, PRD 95 (2017) 104007 **[ABS]** | Kinematics of EC universes with torsion | Background |
| Palle, 1405.3435, JETP 118 (2014) 587; Palle, 2106.08136 (2021) **[ABS]** | EC cosmology with "angular momentum of the Universe" (torsion from light Majorana neutrino spin). Claims help with CMB anomalies and the Hubble tension. 2106.08136 mentions a Planck parity-violation hint but derives no rotation (per web summary). | Non-mainstream. Not a β derivation. |
| Izaurieta, Lepe & Quinzacara, 2608.11453 (2026) **[ABS]** | Weyssenhoff post-bounce acceleration gives at most ~0.46 e-folds. Torsion "fuel" dilutes as a⁻⁶. | Confirms the a⁻⁶ scaling of spin-fluid torsion energy |

Size of cosmological torsion today (scout's order-of-magnitude estimate, **unchecked, for the theorist to redo in `CONVENTIONS.md` units**): in minimal EC, S ~ κ·n₅ with n₅ the net axial number density. For a relic neutrino asymmetry ξ ~ 0.05, n₅ ~ 10⁻¹³ eV³ and κ ~ 1.7×10⁻⁵⁶ eV⁻², so S ~ 10⁻⁶⁸ eV. A CFJ-type β ~ 0.2° over a Hubble time needs p₀ ~ 10⁻³⁵ eV. With the Dobado-Maroto factor ~3×10⁻³, that means S ~ 10⁻³² eV. **Gap ~35 orders of magnitude.** The KRT lab bound (~10⁻²² eV) is not the binding constraint. The binding one is the EC field equation itself.

---

## C. Kaluza-Klein track: EC-KK reductions, Kalb-Ramond torsion axion, KK birefringence

### C1. Einstein-Cartan-Kaluza-Klein reductions

| Paper | Content | Relevance |
|---|---|---|
| Kohler, gr-qc/9808004, IJMPA 15 (2000) 1235 **[ABS]** | Modified 5D KK, partially parallelized, so torsion appears. **The EM field strength is the fifth component of the torsion 2-form.** Gives Einstein-Cartan-Maxwell. The scalar does not induce 4D torsion. | Direct EC-KK prior art. The photon is literally torsion here, so the KK theorist must compare with it. |
| Shankar & Wali, 0904.4038, MPLA 25 (2010) 2121 **[ABS]** | 5D KK in EC with torsion confined to the extra dimension, fixed by constraints. Torsion cancels the KK scalar and vector contributions to 4D Ricci. Cosmic acceleration. | A constrained EC-KK. Shows 5D torsion can cancel KK fields, so a parity-odd coupling is not automatic. |
| Chang, Lebedev, Loinaz & Takeuchi, hep-ph/0005236, PRL 85 (2000) 3765 **[ABS]** | Large extra dimensions with gravity-only bulk and non-symmetric connection: brane fermions induce torsion, giving a universal axial contact interaction | Z-pole bound: scale > 28 TeV (n=2). SN1987A: √n M_S ≥ 210 TeV with light sterile ν. Torsion-specific KK collider and astro bounds. |
| Geng & Luo, 1612.00166, CQG (2017) **[ABS]** | Teleparallel (NGR) KK reduction. Torsion vector as conformal gauge field. | Teleparallel KK, parity-even |
| Lu, 1402.3650; Lu & Huang, 1306.3028 **[ABS for 1402.3650; 1306.3028 seen in search listing only]** | KK-like models of de Sitter/Poincaré gauge gravity with spin fluid | Background |
| Shaya, 2609.02909 (2026) **[ABS]** | EC-KK on S²×CP²: eliminating torsion gives a quartic spin-current term with nonzero vacuum value. Modulus stabilization. | Recent EC-KK with higher-dimensional internal space. Shows the torsion-squared vacuum energy matters for moduli. |
| **Rao, Zheng, Hou, Ou & Zhu, 2608.09299 (2026)** **[KEY]** **SAME-AIM (KK, teleparallel)** | 5D teleparallel: 𝕋/2 + (c/8) η_AB ε T^A T^B n. KK ansatz with **dynamical radion ϕ**, cylinder condition. 4D result (their Eq. 41): (cL₅/4)∫√-g [ϕ η_ab ⋆T^{aμν}T^b_{μν} + 2ϕ³ ⋆F^{μν}F_{μν}], i.e. **a Nieh-Yan term plus an EM Chern-Simons term with the radion as the "axion" (coefficient ∝ ϕ³)**. The Riemannian version gives gravitational CS with a ghost. | **Δα ≈ 3c ln[ϕ(η₀)/ϕ(η_LSS)]**, achromatic at leading order. Fit to Planck PR4 Δα = 0.30° ± 0.11° gives c·ln(ϕ₀/ϕ_LSS) = (1.75 ± 0.64)×10⁻³. Prediction: Δω²_EM = 6 Δω²_GW. Not compared with ACT. Not EC (curvature-free teleparallel). | **Closest existing paper to the KK track.** Our KK program must cite it and state what differs: (i) EC (curvature + torsion) instead of teleparallel, (ii) the same radion that drives β also changes α_em in KK. In the sections read they do not apply fine-structure-constant variation, fifth-force or radion-stabilization bounds (not confirmed for the unread remainder). |

### C2. Kalb-Ramond field as torsion; KR axion and photon optical activity

| Paper | Coupling / action | Result | Relevance |
|---|---|---|---|
| Majumdar & SenGupta, gr-qc/9906027, CQG 16 (1999) L89 **[ABS]** | Gauge-invariant Maxwell-EC coupling. Torsion = KR field strength H, augmented by the U(1) Chern-Simons 3-form. | Parity-violating torsion | The standard "torsion = H = dB + A∧F" construction. Gauge invariance comes from the CS augmentation, not from coupling A to the connection. |
| Kar, Majumdar, SenGupta & Sinha, gr-qc/0006097, EPJC 23 (2002) 357 **[KEY]** | S ⊃ R/16πG - ¼F² - ½H̃² + G^{-1/2} T·H̃, with H̃ = ∂B + (1/3)√G A∧F and T = √G H̃. Dual axion H with an effective ½ φ F⋆F coupling (O(1) after absorbing √G). | Massless axion, □H = 0, ∂₀H = h/R². Rotation Δφ = 2 tan⁻¹(h̃/η) (radiation era), 2 tan⁻¹(h'/3η³) (matter era). Wavelength independent. No numbers. | Root of the KR-torsion birefringence line |
| Maity & SenGupta, hep-th/0311142, CQG 21 (2004) 3379 **[ABS]** | KR in the Randall-Sundrum bulk. The KR-EM CS coupling is warp-enhanced on the brane. | Huge optical rotation, so RS + bulk KR is in tension with radio-galaxy data | **NO-GO-type:** RS warping over-enhances the KR-photon coupling |
| Maity, Majumdar & SenGupta, hep-th/0401218, JCAP 06 (2004) 005 **[ABS]** | Same in RS, applied to CMB TB/EB | Non-observation of CMB parity-odd correlations forces fine-tuned KR-axion parameters | **Early CMB-parity confrontation of a KR-torsion axion (pre-Planck)** |
| **Maity, SenGupta & Sur, hep-th/0507210**, PRD 72 (2005) 066012 **[KEY]** **SAME-AIM (radio galaxies, not CMB)** | 4D action from S¹ KK (Scheme I: α=β=1) or RS (Scheme II): M_P²(R - ½(∂φ)²) - (1/12)e^{-2αφ}H² - e^{-αφ}F² with H = dB + (2β/M_P)A∧F. Axion ξ dual to H, plus a dilaton. | Analytic FRW solution ξ = (1/α)tanh[αh∫dη/a²]. Rotation Δ = 2 tan⁻¹[tanh(h t_L/2)] ≈ h t_L, achromatic. Radio-galaxy limit 6° at z=0.4 gives h ≲ 5×10⁻⁴² GeV, ρ_KR ≲ 2×10⁻⁴⁶ GeV⁴. **Δ ~ 0.2° needs ρ_KR ~ ρ_DM today.** RS needs fine tuning. | **The KK-reduced KR torsion axion has been compared with polarization rotation before, but not with the CMB β or with the a⁻⁶ (stiff) scaling of the massless KR axion energy density back to recombination and BBN.** |
| Das, Jain & Mukherji, hep-ph/0011279, IJMPA 16 (2001) 4011 **[ABS]** | Heterotic-string axion plus dilaton | Dilaton modifies the rotation. Radio-polarization limits on string parameters. | Background |
| Alexandre, Mavromatos & Tanner, 0804.2353, PRD 78 (2008) 066001 **[ABS]** | Time-dependent string backgrounds (graviton, dilaton, KR) | Axion-induced optical activity | Background |
| Ho, Kao, Bamba & Geng, 1008.0486, EPJC (2015), DOI 10.1140/epjc/s10052-015-3426-5 **[ABS]** | CPT-even dimension-6 CS-like term with dynamical KR and scalar fields, plus neutrino current | Sizable CMB rotation | The KR completion that Reig 2611.07578 shows is dual to an axion |
| Bhattacharjee & Chatterjee, 1101.0118 (see A4) | Heterotic KR with YM and gravitational CS | Rotation of EM and GW polarization | |
| Mukhopadhyaya, Sen & SenGupta, hep-th/0204242, PRL 89 (2002) 121101 **[ABS]** | Bulk KR (torsion) in RS: zero mode exponentially suppressed on the brane | Explains why torsion looks weak | KK-specific suppression mechanism. Opposite to the CS-coupling enhancement in hep-th/0311142, which shows the warp dependence is coupling specific. |
| Chakrabarty, Mukhopadhyaya & SenGupta, 1604.00885, PRD 95 (2017) 015007 **[ABS]** | KR (bulk torsion) in RS. CS couplings to gauge bosons from anomaly cancellation. The lowest KK mode is an axion with unsuppressed γγ coupling and suppressed fermion coupling. | Diphoton signal | Collider-side bound and phenomenology for a KK torsion axion |
| Mavromatos (several): 2111.05675 (Universe 7 (2021) 480), 2404.18741, 1510.02457, 1312.5230 **[ABS]** | String-model-independent KR (torsion) axion, gravitational anomaly, RVM inflation, leptogenesis from constant KR torsion backgrounds | Mostly gravitational-anomaly and leptogenesis, not photon β. Says KR torsion backgrounds are negligible today. | Background for the KR axion as a universal string axion |
| Duncan, Kaloper & Olive, NPB 387 (1992) 215, DOI 10.1016/0550-3213(92)90052-D (INSPIRE) **[ABS]** | Axion hair and dynamical torsion from anomalies | Original "torsion axion from anomalies" | Background citation |

### C3. Other KK / extra-dimension photon-rotation proposals

| Paper | Content | Relevance |
|---|---|---|
| Ganguly & Parthasarathy, astro-ph/0308129, PRD 68 (2003) 106005 **[ABS]** | Vacuum optical activity from the KK scalar in an external EM field | Needs a background field. Not isotropic β. |
| Arkani-Hamed, Cheng, Creminelli & Randall, hep-th/0301218, PRL 90 (2003) 221302 **[ABS]** | Extranatural inflation: A₅ Wilson line on S¹ as an inflaton, with non-local potential. Mentions quintessence. | If a 5D CS term A∧F∧F is present, the A₅ zero mode couples as A₅ F F̃. A natural shift-symmetric KK pseudoscalar for β. Not torsion. |

---

## D. Gauge invariance: why minimal photon-torsion coupling is dropped, and who couples anyway

| Paper | Content | Relevance |
|---|---|---|
| **Hojman, Rosenbaum & Ryan, PRD 17 (1978) 3141**, DOI 10.1103/PhysRevD.17.3141 (INSPIRE) **[ABS]** | "Gauge invariance, minimal coupling, and torsion": torsion derived from a scalar potential so that minimal coupling stays gauge invariant | Classic resolution: restrict torsion to a gradient (trace) form |
| **de Sabbata & Gasperini, PRD 23 (1981) 2116**, DOI 10.1103/PhysRevD.23.2116 (INSPIRE) **[ABS]**; also PLA 83 (1981), DOI 10.1016/0375-9601(81)90510-7 | "Gauge invariance, semiminimal coupling, and propagating torsion": couple the photon to torsion only in the field equations (semiminimal), keeping U(1) | Classic alternative. Cited by Alexander et al. as the origin of the CFJ-like term. |
| Benn, Dereli & Tucker, PLB 96 (1980) 100, DOI 10.1016/0370-2693(80)90221-X (INSPIRE) **[ABS]** | Gauge-field interactions in spaces with arbitrary torsion | Form-language statement of the minimal-coupling problem |
| Poplawski, 1108.6100 (see B) | Explicit: minimal coupling of F to the connection breaks U(1) | |
| Nieh, 1701.07132 (2017) **[ABS]** | Maxwell-Dirac on Riemann-Cartan space: restores gauge invariance with torsion by modifying couplings, motivated by scale invariance | Modern attempt at keeping minimal coupling consistent |
| Fresneda, Baldiotti & Pereira, 1404.3231, Braz. J. Phys. 45 (2015) 353 **[ABS]** | Gauge-invariant propagating-torsion model coupled to Maxwell and charges | Another gauge-invariant model |
| de Andrade & Pereira, gr-qc/9708051, IJMPD 8 (1999) 141 **[ABS]** | Teleparallel equivalent of GR: EM couples to and sources torsion while keeping gauge invariance | Teleparallel version, which differs from EC |
| Itin & Hehl gr-qc/0307063; Rubilar et al. gr-qc/0305049; Das et al. 0908.0629 (see A1) | Non-minimal couplings built from F, so gauge invariant by construction | The standard escape from the problem |
| Majumdar & SenGupta gr-qc/9906027 (see C2) | Gauge invariance from the H = dB + A∧F (CS-augmented) torsion | The string/KR escape |

The consensus across these papers: F = dA does not involve the connection. Minimally coupling A through ∇ breaks U(1) unless torsion is restricted (Hojman et al.) or the coupling is moved into the field equations (de Sabbata-Gasperini). Every photon-torsion effect therefore comes from (i) gauge-invariant non-minimal operators in F and T, (ii) fermion loops (CFJ, regularization dependent), or (iii) a torsion-dual pseudoscalar (KR/BI/Nieh-Yan axion) with an anomaly or CS coupling to F F̃.

---

## E. Fermion-current birefringence: the effective operator behind the torsion-fermion route

In minimal EC, S_μ is algebraically ∝ κ J⁵_μ. Putting this into a Dobado-Maroto-type CFJ term gives an operator of the form (coefficient) × J_μ A_ν F̃^{μν}. That is exactly the "fermion-current Chern-Simons" operator studied in this line:

| Paper | Content | Relevance |
|---|---|---|
| Geng, Ho & Ng, 0711.4617; 0706.0080; 0801.0024 (2007) **[search listing only, IDs verified via arXiv API search]** | Neutrino-current-induced cosmological birefringence; neutrino asymmetry and birefringence; CPT-conserving version | Origin of the operator j_μ A_ν F̃^{μν}/M² |
| Zhou, Huang & Geng, 2302.11140, JCAP 07 (2023) 053 **[ABS]** | Rederives the β formula (corrects an earlier mistake). Neutrino (EMPRESS ξ_νe) or asymmetric-DM currents can give the Planck β. | **SAME-AIM for the operator, not for torsion** |
| Kasai & Obata, 2610.09748 (2026) **[ABS]** | Neutrino-current CS: first anisotropic-β analysis (lepton isocurvature), Affleck-Dine lepton asymmetry explains isotropic β | Anisotropy template for the current route |
| Yin & Komatsu, 2610.02667 (2026) **[ABS]** | EB spectrum shape for a neutrino-asymmetry current vs axion DE vs DM: the neutrino case is suppressed at ℓ ≲ 20 | Discriminating test (LiteBIRD) |
| **Reig, 2611.07578 (2026)** **[KEY]** **NO-GO** | The operator -(b_J/2M²) j_μ A_ν F̃^{μν} is not gauge invariant. The gauge-invariant UV completion needs a KR 2-form dual to a massless axion. The signal fixes g_aγ/f_ℓ = 1/Λ*² with Λ* ≈ 10 TeV. SU(2)_L gives g_aγ g_ae ≈ 5×10⁻¹² GeV⁻¹, excluded by electron g-2 + ALPS II, helioscopes and stellar cooling. Anomaly-free conserved currents fail differently. | **Kills the neutrino-current route. Our torsion-mediated fermion route maps onto the same operator with 1/Λ*² ~ ακ (scout reading), i.e. ~28 orders weaker than needed.** Torsion is not mentioned in the paper, so the mapping must be done by our theorist. |

---

## Verification log

All arXiv IDs in this file were returned by the arXiv API (`id_list` query) with matching titles on 2026-11-20. Scripts: `work/20261120-002-lit-prior-torsion/arxiv_verify_ids.py`, `arxiv_search.py`, `inspire_search.py`. DOIs of pre-arXiv papers were returned by INSPIRE-HEP. Not separately verified: Lu & Huang 1306.3028 (seen only in a search listing).
