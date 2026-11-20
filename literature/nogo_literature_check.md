# Literature check of the program's no-go statements

Task 20261120-013 (literature scout), 2026-11-20. Logs and scripts: `work/20261120-013-lit-nogo-check/`.

Every arXiv ID below was returned by the arXiv API (`id_list`, title and authors match) on 2026-11-20 (`abstracts*.log`). Reading depth:
- **[FT]**: relevant sections read in the arXiv HTML full text (tool-assisted grep plus reading of the surrounding paragraphs).
- **[ABS]**: abstract only.
- **[LIST]**: title seen in a search listing only.

Verdict labels per item: **confirmation** (a paper already states our result), **partial** (states part of it, or the same result in a narrower setting), **contradiction**, **not addressed**.

## Searches run

arXiv API (`arxiv_search.log`): torsion + "cosmic birefringence"; "Einstein-Cartan" + birefringence; torsion + "Carroll-Field-Jackiw"; torsion + "CMB polarization"; "Kalb-Ramond" + birefringence; "Kaluza-Klein" + "cosmic birefringence"; "extra dimension" + "cosmic birefringence"; "Nieh-Yan" + birefringence; "Barbero-Immirzi" + birefringence; torsion + "polarization rotation"; torsion + birefringence + photon; "Kaluza-Klein" + birefringence; radion + birefringence; "Wilson line" + birefringence; teleparallel + birefringence + electromagnetic; "torsion axion"; torsion + "Lorentz violation" + photon; "Kalb-Ramond" + axion + QCD; "gravitational axion" + birefringence; Holst + photon + parity; "Einstein-Cartan" + photon + parity; ULA + DESI + dark energy.

INSPIRE full text, 2020-2026 (`inspire_fulltext.log`): "cosmic birefringence" AND "Einstein-Cartan" (19 hits); "cosmic birefringence" AND "Kaluza-Klein" AND torsion (13); "Kalb-Ramond" AND "cosmic birefringence" (18). All titles screened.

INSPIRE citing papers (`citing.log`): 2608.09299 (record 3189313) has **0 citations**. 1502.03694 (record 1344204) has 53; all titles screened, none is about birefringence; the relevant one is 2506.11836 (below).

Result of the search: no paper computes cosmic birefringence from Einstein-Cartan torsion and states it is zero or too small. The pieces of the no-go exist in separate papers, listed below. Nothing found contradicts a verdict in `CLAIMS.md`. One paper (2506.11836) contradicts the reasoning behind EC-03 "reading A" (the QCD mass), while supporting the verdict.

---

## (1a) Minimal 4D EC with a minimally coupled photon gives β = 0 (EC-00)

**Verdict: partial.** No paper states it for cosmic birefringence. Two papers state the two halves.

| Paper | Depth | What it shows | Relation to EC-00 |
|---|---|---|---|
| Mukhopadhyaya, Sen & SenGupta, **hep-th/0207165** | [FT] Sec. 2.2-2.3, Sec. 4 | Minimal coupling of A to a 2-index torsion breaks U(1) (Sec. 2.2). For fully antisymmetric (KR-dual) torsion, with the gauge-invariant CS-augmented H, the parity-violating gauge-sector Lagrangian "vanishes identically", so the U(1) sector is parity conserving; parity violation survives only in the spin-1/2 sector (Sec. 2.3, Eq. 16 and text after it; Sec. 4). | Confirms that EC torsion of the type sourced by Dirac spin gives no parity-odd photon term at tree level. Does not compute β or discuss cosmology. |
| Nakai, Namba & Obata, **2310.09152** | [FT] Sec. 2, Sec. 6 | In a homogeneous neutral background only a CS-type (E·B) operator can give isotropic β. No SM / LEFT operator, four-fermion interactions included, reproduces the reported angle (Sec. 3: CS-type operators in SMEFT/LEFT, four-fermion terms listed there); the CνB contribution is computed and too small (Sec. 4, 6). Gravity is explicitly omitted (Sec. 2: the only healthy non-minimal term R̃FF̃ is parity even). | Minimal EC with Dirac matter is GR plus an axial four-fermion contact term, which is a LEFT operator, so EC-00 and EC-02/05 fall inside their no-go. That mapping is ours; torsion is not mentioned. |
| Shapiro, **hep-th/0103093** (Phys. Rept.) | [ABS] | Review: torsion couplings to matter, anomalies, bounds; consistent propagating axial torsion must be heavier than the heaviest fermion coupled to it. | Standard reference for "gauge fields do not couple minimally to torsion". The gauge-invariance statement was not re-read in the full text in this task. |
| Kostelecky, Russell & Tasson, **0712.4393** | [ABS] | Minimal torsion-fermion couplings map onto SME b_μ; 19 of 24 components bounded at ~10⁻³¹ GeV. | Fermion-sector bound only. Does not discuss a photon-sector (k_AF) term from torsion. Not addressed. |

## (1b) Pure EC-KK gives no parity-odd photon coupling at all orders (KK-01)

**Verdict: not addressed (closest is partial).**
- hep-th/0207165 [FT] Sec. 3-4 extends the KR-torsion parity statement to ADD and RS compactifications, with the photon on the brane: gravity and U(1) sectors stay parity conserving. That is not the S¹ reduction with a geometric KK photon that KK-01 covers.
- Rao et al. **2608.09299** [FT, re-read v3]: the parity-odd EM term exists only through the explicit 5D parity-violating term with coefficient c; at c = 0 the reduction is parity even. Consistent with KK-01 but not stated as a theorem. The P' = diag(1,-1,-1,-1,-1) argument of task 009 was not found in any paper.

## (2) Torsion-sourced routes (spin density + CFJ, fermion loop) fall short by 28-35 orders (EC-02, EC-05)

**Verdict: partial.**
- Dobado & Maroto, astro-ph/9706044 (already [KEY] in `prior_work_torsion_birefringence.md`): need |S| ~ 10⁻³⁰ eV and remark that this requires massless or condensed torsion, not spin-sourced EC torsion. Qualitatively the same conclusion, no number for the EC spin density.
- Nakai, Namba & Obata, 2310.09152 [FT]: CνB with parity-violating neutrino-photon interactions is too small (Sec. 6). Same operator class as EC-02.
- Reig, 2611.07578 (already [KEY]): the j·A F̃ operator is not gauge invariant and its UV completion is excluded. Torsion not mentioned.
- No paper puts the EC spin density through the CFJ chain with numbers. The 28.5 and ≥ 29 order shortfalls are not in the literature.

## (3) Torsion axions: QCD mass and washout, or a generic ALP with f ~ M_P (EC-03, EC-03', EC-04, KK-02, KK-03)

**Verdict: partial confirmation, plus one contradiction of reasoning (not of verdict).**

| Paper | Depth | What it shows | Relation |
|---|---|---|---|
| Karananas, Shaposhnikov & Zell, **2506.11836** | [FT] Sec. 2, Outlook | A scalar of gravitational (torsion) origin coupled to fermions only through ∂_μφ J_A^μ has an exact shift symmetry: "No effect, perturbative or non-perturbative, can break this symmetry", so no potential is generated. The gluon-anomaly term and the fermion-mass term cancel exactly (their Eqs. 9-11). They conclude earlier gravitational-axion arguments (they cite 1502.03694, Mercuri 0902.2764, Lattanzi-Mercuri) are flawed. Shift-breaking couplings Z(φ)∂φJ (Eq. 13) can reach QCD at loop level but are "either too weak to be viable or suffer from a severe fine-tuning problem". Only Weyl-invariant EC (2406.11956) is left as a route. | **Contradicts the reasoning of EC-03 reading A** (torsion axion gets a QCD mass 1.4-8.6×10⁻¹² eV from 1502.03694). **Confirms reading B** (task 008 open issue 1, STATE open question 1): no QCD potential. The paper treats gluons only. Applying the same cancellation to the photon anomaly with massive charged fermions gives no low-energy φFF̃ from the derivative coupling (task 008 reading B); that extension is ours. EC-03 stays dead (no mass, and g_eff ≈ 0). **EC-03' loses its torsion origin of g**: with a pure derivative coupling, g = 8α/(πf) is not generated, so φFF̃ must be added by hand. |
| Karananas, Shaposhnikov & Zell, **2406.11956** | [ABS] | Weyl-invariant EC + SM contains one extra ALP that solves strong CP; its small mass set by tiny local-Lorentz gauge couplings. | If it solves strong CP, it couples to GG̃ and gets a QCD mass, which washes out β by the EC-03 argument. Not checked further. Not addressed for photons. |
| Karananas, **1805.08781** | [ABS] | A torsion-related ALP could protect the axion solution against gravitational effects. | Not addressed for photons. |
| Castillo-Felisola et al., 1502.03694 | (already [KEY]) | Torsion axion with universal quark coupling, f = √(8/3)M_P, QCD mass ~10⁻¹² eV. | Basis of EC-03 reading A; criticised by 2506.11836. |
| de la Cruz Dombriz, Maldonado Torralba & Mota, **2112.03957** | [ABS] | The stable 0⁻ pseudoscalar of quadratic PGT as cold dark matter. | Same mode as EC-04. Photon coupling not mentioned in the abstract. Not addressed. |
| Trukhanova & Obukhov, **2410.01355** | [ABS] | Quadratic PGT: axial torsion sourced by EM helicity and spin density; polarization rotation of a plane wave only in external E or B fields. | Not isotropic β. Consistent with "torsion needs an added coupling" (EC-01/EC-04). |
| Girmohanta, Qiu & Wang, **2303.02852** | [FT, grep] | Quintessence axion from a gauged chiral U(1) on two branes in a flat 5th dimension; f ≈ 10¹⁷ GeV; explains isotropic β via an anomaly-induced φFF̃. Brane-Higgs axion, not a Wilson line. | Extra-dimension axion explaining β as an ordinary ALP. Supports "KK-02/KK-03 are ALPs in disguise"; geometry only sets f and the potential. |
| Nakagawa, Nakai & Qiu, **2503.18924** | [FT] Sec. II-III | Subdominant axion with 3.0×10⁻³³ ≲ m ≲ 1.1×10⁻³² eV (2H₀-7H₀), f ~ 10¹⁷ GeV, c_γ ~ O(1), explains β = 0.3° and a DESI-like w₀-w_a at 1σ. | Same mass window as EC-03' (best m ≈ 6×10⁻³³ eV). Generic-ALP prior art for the EC-03' parameter region. |
| Maity, Majumdar & SenGupta hep-th/0401218; Maity, SenGupta & Sur hep-th/0507210 | (already in prior_work) | KR torsion axion vs pre-Planck CMB TB/EB and radio galaxies. | The KR-axion washout / QCD question is not discussed. |
| Mavromatos et al., 2001.03465, 2108.02152 | [LIST] | String-inspired KR (torsion) axion: anomalies, inflation, axion DM. | Not about photon β. Not addressed. |

No paper was found that applies the polarization washout (2512.21888 type) to a torsion or KR axion.

## (4) Radion × FF̃ (2608.09299 type) excluded by Cassini / inverse-square-law (KK-04)

**Verdict: not addressed. No contradiction.**
- 2608.09299 v3 [FT]: "The stabilization of φ and the associated moduli-fixing problem are important topics in their own right, though they lie beyond the scope of the present work" (paragraph after the 4D field equations, Sec. II). No Cassini, PPN, fifth-force or α_em-variation bound is applied anywhere (grep for cassini/fifth/PPN/Brans/fine-structure: no hits in the physics text).
- INSPIRE: 0 citing papers as of 2026-11-20. The KK-04 kill is not in the literature.

---

## Goal 3: do 2601.13624 or 2512.21888 treat a torsion-motivated coupling?

**No.** Full-text grep of both (torsion, Einstein-Cartan, Kaluza-Klein, Kalb-Ramond, extra dimension, string): no torsion or KK motivation.
- **2601.13624** (Yin et al.) [FT]: L = -½(∂φ)² - V - ¼F² - ¼gφFF̃ with V = Λ⁴[1 - cos(φ/f)]³, θ_i fixed to 2.7 (Sec. II-III). Fit gives g M_Pl = 0.129 (+0.031, -0.075) for the ACT-EB combination (Table 1, already in `data/birefringence_measurements.csv`); M_Pl convention not stated.
  - EC-03': g_T M_P = 4.67×10⁻²¹ GeV⁻¹ × 2.435×10¹⁸ GeV = 0.011 (reduced) or 0.057 (non-reduced 1.22×10¹⁹ GeV). The non-reduced value sits at their 68% lower edge (0.054); the reduced value is ~10× below the central value. EC-03' is not an EDE field (frozen through recombination, m ~ 6×10⁻³³ eV), so this fit does not constrain it directly. An EDE version of a torsion-fixed coupling would be marginal at best.
  - KK-02, KK-03: nothing relevant.
- **2512.21888** (Shao, Obata, Zhang) [FT]: multi-field ALP DM, 10⁻²⁵-10⁻²² eV; motivation is only that the washout makes the single-field coupling too small. Relevant statement for our candidates: an ALP that oscillates before recombination (m ≫ 10⁻²⁸ eV) gives β only through the observer term and is washout-limited. KK-02 and KK-03 escape this only for m ≲ 10⁻²⁸ eV (the task 008 generic-ALP window), which is where their windows are already set.

---

## Goal 2: energy-density bound for an ultralight field frozen at recombination, m ~ 10⁻³³-10⁻³² eV

Numbers in `data/ultralight_fraction_bounds.csv` (cited, per mass).

| mass | bound on today's density | source |
|---|---|---|
| 10⁻³³ eV | **none** (Ω_a h² = 0.203 ± 0.119, flat posterior, degenerate with Λ) | 1708.05681 Table 3 |
| 10⁻³² eV | Ω_a h² < 0.014 (95%), i.e. < 11.8% of Ω_c h² | 1708.05681 Table 3 (Planck 2015 T+P+lens) |
| 10⁻³² eV | Ω_a h² < 0.00833 (95%), < 6.9% of ρ_DM | 2301.08361 Table 2 (Planck 2018) |
| 10⁻³² eV | Ω_a h² < 0.00321 (95%), < 2.7% of ρ_DM | 2301.08361 Table 2 (Planck 2018 + BOSS) |
| 1+z_c ≲ 10 | no bound on f_zc (degenerate with DE) | 1806.10608 abstract |
| 10⁻³²-10^-25.5 eV | Ω_a/Ω_d < 0.05, mass-marginalised | 1410.2896 abstract, Table 1 |

Consequences for EC-03' (independent spectator ODE in `work/20261120-013-lit-nogo-check/energy_check.py`, reproduces task 008):
- β = 0.215° with g_T needs ρ_φ,0/ρ_DM = 0.11 at the optimum m ≈ 6×10⁻³³ eV, ≤ 0.21 only for 5×10⁻³³ ≲ m ≲ 8×10⁻³³ eV, and 0.19 at 10⁻³² eV.
- At 10⁻³² eV this is excluded by Planck 2018 alone (< 6.9%) and by Planck+BOSS (< 2.7%).
- At 6×10⁻³³ eV (3H = m at z ≈ 0.6) there is **no published per-mass CMB/LSS bound**: 1708.05681 and 1806.10608 show the field is degenerate with dark energy below about 10⁻³² eV, and 2503.18924 shows a subdominant field there is compatible with DESI-era w₀-w_a. The 5% cap of task 008 is neither confirmed nor the right kind of bound: the binding constraint in this window is a late-time expansion fit (BAO + SN), which nobody has published for a sub-dominant thawing field plus Λ with this mass. Using the nearest tabulated point (10⁻³² eV, Planck+BOSS 2.7%) as a conservative cap gives β_max ≈ 0.136° × √(2.7/5) ≈ 0.10°.
- 2506.11836 makes the cap question secondary: if the torsion coupling is purely derivative, g_T is not generated at all.

## Summary for the overseer

- Is the no-go published? **Partly.** The gauge-sector parity statement for EC/KR torsion (hep-th/0207165), the "no SM/LEFT operator gives β" result (2310.09152), and the "derivative torsion coupling gives no axion potential" result (2506.11836) exist separately. Nobody has assembled them into "EC torsion cannot explain cosmic birefringence", and the EC-KK all-orders P' argument, the 28-35 order numbers, and the Cassini/ISL kill of the 2608.09299 radion are not in the literature.
- Contradiction: only with EC-03 reading A (QCD mass from 1502.03694). 2506.11836 says that mass is not generated. The EC-03 verdict holds through reading B. EC-03' loses its torsion-generated coupling under the same argument (our extension of 2506.11836 to photons, to be checked by the theorist).
