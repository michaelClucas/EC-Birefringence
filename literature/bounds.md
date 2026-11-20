# Independent bounds for torsion-birefringence candidate mechanisms

Task `work/20261120-005-bounds`, iteration 1. Machine-readable table: `data/bounds.csv` (published values, conversions to eV). Inputs that are not bounds (unit constants, CMB dipole, H0, Omega_m): `data/bounds_inputs.csv`. Derived numbers: `work/20261120-005-bounds/derived_bounds.csv`, produced by `convert_bounds.py`.

Units follow `CONVENTIONS.md`: natural units, eV, signature (-,+,+,+), S^mu = (1/6) eps^{mu nu rho sigma} T_{nu rho sigma}, T_mu = T^lambda_{mu lambda}. Unit factors: 1 GeV = 1e9 eV; 1 GeV^-1 = 1e-9 eV^-1; 1 m^-1 = hbar c = 1.9733e-7 eV; 1 yr^-1 = 2.0858e-23 eV; H0 = 67.4 km/s/Mpc = 1.438e-33 eV.

Verification status: "abstract" in `location_in_paper` means only the arXiv abstract was opened. Rows citing the Kostelecky-Russell Data Tables (arXiv:0801.0287v19, 5 Feb 2026) were read from the tables themselves; the primary papers behind them were opened only where stated.

## 1. Torsion couplings to fermions (SME b_mu)

### Mapping torsion -> b_mu

Kostelecky, Russell, Tasson, PRL 100, 111102 (2008), arXiv:0712.4393:
- Eq. (1): T_{lambda mu nu} = (1/3)(g_{lambda mu} T_nu - g_{lambda nu} T_mu) - eps_{lambda mu nu kappa} A^kappa + M_{lambda mu nu}, with a trace vector T_mu (exact index contraction not legible in the extracted PDF text) and A_mu = (1/6) eps_{mu alpha beta gamma} T^{alpha beta gamma}.
- Eq. (2): effective constant-torsion fermion Lagrangian with couplings xi_j^(4), xi_j^(5). "The special case of minimal coupling is recovered for xi_4^(4) = 3/4 with other couplings zero."
- Eq. (3): b_mu - m g^(A)_mu = -(xi_2^(4) - m xi_8^(5)) T_mu - (xi_4^(4) - m xi_9^(5)) A_mu (the m vs 2m factor in front of xi_8, xi_9 was not legible in the extracted PDF text; irrelevant for minimal coupling). So with minimal coupling **b_mu = -(3/4) A_mu**, and the trace T_mu does not couple (xi_2 = 0).
- Eq. (6): minimal-coupling bounds |A_T| < 2.9e-27 GeV, |A_X| < 2.1e-31 GeV, |A_Y| < 2.5e-31 GeV, |A_Z| < 1.0e-29 GeV (from 2004-2006 dual-maser and torsion-pendulum data). Check: 1.6e-31 / 2.1e-31 = 0.76 = 3/4.

Conventions of KRT (via their Ref. [5], Kostelecky PRD 69 105009 (2004), hep-th/0312310, App. A, opened): local Minkowski metric diag(-1,+1,+1,+1), eps_{0123} = +1, connection with derivative index first. These match `CONVENTIONS.md`, so **A_mu(KRT) = S_mu(ours)** at the level of magnitudes. The overall sign of b relative to S also depends on the SME sign of the b-term and the gamma^5 convention; it does not matter for magnitude bounds. A theorist should confirm the sign before using it in a sign-sensitive argument.

Second mapping in the literature: Shapiro, Phys. Rep. 357, 113 (2002) writes the axial vector without the 1/6, S^mu = eps^{alpha beta mu nu} T_{alpha beta nu}, and the minimal coupling as eta = 1/8, so b = (1/8) x 6 A = (3/4) A in KRT normalisation. The two agree. (Shapiro was not opened in this task; this is a consistency note only.)

Lammerzahl, PLA 228, 223 (1997), gr-qc/9704047 (abstract): axial torsion K <= 1.5e-15 m^-1 = 3.0e-22 eV from Hughes-Drever experiments. Normalisation of K vs S_mu not checked; superseded anyway.

### Current b~ bounds (Data Tables v19)

| Species | Coefficient | Value (GeV) | Source |
|---|---|---|---|
| electron | b~_X, b~_Y | (-0.9 +- 1.4)e-31 each | Heckel et al. PRD 78 092006 (0808.2673), Table D6 |
| electron | b~_Z | (-0.3 +- 4.4)e-30 | same |
| electron | (1/2)(b~_T + d~_- - 2g~_c - 3g~_T + 4d~_+ - d~_Q) | (0.9 +- 2.2)e-27 | same |
| neutron | b~_X, b~_Y | (0.1 +- 1.6)e-33, (2.5 +- 1.6)e-33 | Brown et al. PRL 105 151604 (1006.5425), abstract and Table D15 |
| neutron | b~_perp | < 6.7e-34 (68% CL) | Allmendinger et al. PRL 112 110801 (1312.3225), abstract |
| neutron | b~_T | < 1e-27 (order of magnitude) | Data Tables S2 / D15 ref [91] |
| proton | b~_perp | < 6e-32 | Brown 2010, Table D10 |
| proton | b~_T | < 1e-26 (order of magnitude) | Data Tables S2 / D10 ref [91] |

The b~ are combinations, e.g. b~_J = b_J - m d_{J0} - (1/2) eps_{JKL}(m g_{KLT} + H_{KL}); with torsion as the only source of Lorentz violation and minimal coupling, b~ = b. Neutron and proton values from nuclear comagnetometers depend on nuclear spin models.

### Converted to torsion (derived, `derived_bounds.csv`)

Minimal coupling, |S| = (4/3)|b|. Measurements quoted as c +- sigma are turned into a bound |c| + 2 sigma.

| Quantity | Bound (eV) | From |
|---|---|---|
| |S_X|, |S_Y| (Sun-centred frame) | 8.9e-25 | neutron, Allmendinger 2014 (68% CL) |
| same | 4.9e-22 | electron, Heckel 2008 |
| |S_T| directly | 1.3e-18 | neutron b~_T (order of magnitude) |
| |S_0| in CMB frame via boost | **7.3e-22** | neutron b~_perp / v_perp |

Boost argument: an FRW-compatible axial torsion is purely timelike, S_mu = (S_0, 0) in the CMB frame. In the Sun-centred frame it acquires spatial components S_J = -gamma v_J S_0, constant in that frame, which sidereal searches measure. The Sun moves at v = 369.82 km/s = 1.234e-3 c toward (l, b) = (264.021, 48.253) deg (Planck 2018 I, arXiv:1807.06205, Table 3, opened via ar5iv). That direction has Dec = -6.9 deg, so v_perp = v cos(Dec) = 1.225e-3 and |S_0| < |S_perp| / v_perp. This gives a bound on the cosmological S_0 about 2000 times stronger than the direct b~_T bound. It assumes the torsion is minimally coupled to neutrons, constant over the solar system, and not screened locally (in minimal EC, torsion is algebraic and vanishes outside spin-polarized matter, so a cosmological S_0 exists today only if there is a net cosmic axial spin density or an extension makes torsion propagate).

## 2. Photon sector: CFJ (k_AF) and the CMB rotation

SME photon term (Data Tables Table P62): L contains (1/2)(k_AF)_kappa eps^{kappa lambda mu nu} A_lambda F_{mu nu}. Ours: -(1/2) p_mu A_nu Ftilde^{mu nu} = -(1/4) p_mu eps^{mu nu rho sigma} A_nu F_{rho sigma}. Matching gives **p_mu = -2 (k_AF)_mu** up to the eps-sign convention (SME flat-space papers use signature (+,-,-,-); sign not checked). Spherical form: k^(3)_(V)00 = -sqrt(4 pi) (k_AF)_T (Data Tables P70; Kostelecky-Mewes PRD 80 015020 (2009), arXiv:0905.0031, Table VI).

Cosmological rotation, KM 2009 Eq. (158): for d = 3, delta psi = Int_0^z dz' / ((1+z') H(z')) Sum_jm Y_jm k^(3)_(V)jm, i.e. delta psi = -(k_AF)_T t_L for an isotropic constant coefficient, with t_L the lookback time. This is consistent with the axion baseline beta = (1/2) Int p_0 dt for p = -2 k_AF.

| Row | Value | Source |
|---|---|---|
| k^(3)_(V)00 from Minami-Komatsu beta = 0.35 deg | (7.32 +- 2.94)e-45 GeV | Data Tables D19, refs 2011.11254 + 2311.16368 |
| k^(3)_(V)00 Planck 2015 | (-0.29 +- 0.05 +- 0.26)e-43 GeV | D19, 1605.08633 |
| |k_AF| (CMB) | < 7.4e-45 GeV | D19, Caloni et al. 2212.04867 |
| k^(3)_(V)10 (dipole) | (0.09 +- 0.06)e-43 GeV | D19, Contreras et al. 1705.06387 |
| CFJ radio galaxies | < 6e-43 GeV | D19, CFJ PRD 41 1231 (1990) |
| k^0_AF pulsars | < 1e-36 GeV | D19, 2503.20626 |

The CMB beta is itself the dominant k_AF measurement. ACT DR6 beta = 0.215 +- 0.074 deg (`data/birefringence_measurements.csv`, row DPK25-beta-baseline) corresponds, for a constant isotropic coefficient and t_L = 13.80 Gyr = 6.61e32 eV^-1, to **|(k_AF)_T| = (5.7 +- 2.0)e-36 eV, |p_0| = 1.1e-35 eV, |k^(3)_(V)00| = 2.0e-35 eV**.

Open discrepancy: applying KM Eq. (158) to beta = 0.35 deg gives k^(3)_(V)00 = 3.3e-35 eV, while the Data Tables row (from Nilsson & Le Poncin-Lafitte) gives 7.3e-36 eV, a factor 4.47 smaller. Not resolved here; possibly a different time variable (conformal time, distance) or a 2 between polarization angle and Stokes angle. Use our own formula from `work/` (theorist), not the tabulated conversion.

Structural note (KRT 2008, closing paragraphs): F_{mu nu} defined with the torsionful covariant derivative breaks U(1); gauge-invariant nonminimal torsion-photon couplings in a constant background are dimension 5 in the SME. A CFJ term with p_mu proportional to S_mu is gauge invariant only up to a total derivative if dS = 0 (S_mu a gradient), which holds for FRW S_mu = (S_0(t), 0).

## 3. Axion-photon coupling g_agg

| Bound | g (GeV^-1) | g (eV^-1) | Mass range | CL | Source |
|---|---|---|---|---|---|
| CAST 2024 | < 5.8e-11 | < 5.8e-20 | m <~ 0.02 eV | 95 | 2406.16840 abstract |
| HB, globular clusters | < 0.66e-10 | < 6.6e-20 | wide (not in abstract) | 95 | 1406.6053 abstract |
| R2 parameter | < 0.47e-10 | < 4.7e-20 | not in abstract | not stated | 2207.03102 abstract |
| SN1987A gamma rays | < 5.3e-12 | < 5.3e-21 | m < 4.4e-10 eV | not stated | 1410.3747 abstract |
| Magnetic white dwarfs | < 5.4e-12 | < 5.4e-21 | m <~ 3e-7 eV | 95 | 2203.04319 abstract |
| Chandra NGC 1275 | < 6-8e-13 | < 8e-22 | m < 1e-12 eV | 99.7 | 1907.05475 abstract |

For ultralight pseudoscalars (m < 1e-12 eV), which is where a birefringence field must sit, the best astrophysical bound is **g < 8e-13 GeV^-1** (Chandra). IAXO projections were not collected. Hoof & Schulz (2212.09764) improve SN1987A by a factor 1.4 (no number in abstract). These apply to the Holst/Barbero-Immirzi field, the Nieh-Yan pseudoscalar, a KK Wilson-line axion and the Kalb-Ramond axion, each once its coupling is normalised to L = -(g/4) phi F Ftilde.

## 4. Ultralight axion dark matter

- Lyman-alpha: m > 2e-20 eV (95%), Rogers & Peiris 2007.12705 (abstract); m > 2.0e-21 eV (2 sigma, conservative), Irsic et al. 1703.04683 (abstract). Both for ULA = all DM.
- Fraction: Omega_a/Omega_d < 0.05 (95%) for 1e-32 < m < 1e-25.5 eV, Hlozek et al. 1410.2896; f < 0.1 for m <= 1e-26 eV and f < 0.01 for 1e-30 to 1e-28 eV, Rogers et al. 2301.08361 (CL not in abstract).
- Washout (2512.21888, opened, Sec. 2): 1 - J0(g phi_*) <~ 5.8e-3 at 95% (Eq. 11, from Fedderke et al. 2019); hence g <~ 9.6e-16 GeV^-1 (m/1e-24 eV)(f Omega_DM h^2/0.11933)^(-1/2) (Eq. 12). Required for ACT beta: 2.16e-15 to 6.64e-15 GeV^-1 (m/1e-24 eV) kappa^(-1/2) (beta/0.215 deg) (Eq. 6). Washout sits about 2x below the required coupling, so single-field ALP DM in 1e-25 to 1e-22 eV is disfavoured. The paper uses f = 0.01 from Winch et al. 2024 for multi-field fits and says ULA can be all DM only for m > 1e-24 eV.
- Isocurvature: r f_dm < 0.08 (m_a/1e-27 eV)^(-1/2), Liu, Hu & Grin 2608.07803 (abstract). The Planck 2018 X axion-I beta_iso limit (Sec. 9.4.1 of 1807.06211) could not be reached (PDF too large, ar5iv truncated): **not verified**.
- SPT-3G oscillating-rotation bound < 0.071 deg (95%) is quoted in 2512.21888 footnote 1 (Ferguson et al. 2022); original not opened, not in the CSV.

## 5. Kaluza-Klein-specific

- alpha variation. Today: alpha_dot/alpha = (1.0 +- 1.1)e-18 yr^-1 = (2.1 +- 2.3)e-41 eV = (1.5 +- 1.6)e-8 H0 (Lange et al. 2010.06620, abstract; Filzinger et al. 2301.03433 improve this, no number in abstract). Oklo: |Delta alpha/alpha| < 1.1e-8 (95%) (1503.06011). Quasar z = 1.15: (1.3 +- 1.3 +- 0.4)e-6 (ESPRESSO, 2112.05819). Recombination: (3.6 +- 3.7)e-3 (68%) (Planck int. XXIV, 1406.7482). No BBN-only alpha bound collected.
  In KK on S^1 the 4D gauge coupling depends on the radius, so a rolling radion shifts alpha with Delta alpha/alpha = O(1) x Delta ln R; these numbers bound that. The exact coefficient depends on frame and ansatz (KK theorist, task 004).
- Equivalence principle: MICROSCOPE eta(Ti,Pt) = (-1.5 +- 2.3 stat +- 1.5 syst)e-15 (2209.15487).
- Inverse-square law: gravitational-strength (|alpha| = 1) Yukawa range < 38.6 um (95%), Lee et al. 2002.11761, so a light radion/dilaton with gravitational-strength matter couplings needs m > 5.1e-3 eV. Radion couplings are O(1) times gravitational, so this applies up to O(1) factors.
- Size of the extra dimension. If SM fields propagate in 5D (UED on S^1/Z2): 1/R > 1.4 TeV (95%), R < 1.4e-19 m (Deutschmann, Flacke, Kim 1702.00410, Conclusions). In classical KK the photon is g_{mu 5} and charged-particle quantization ties e to 1/(R M_P), which puts R near the Planck length; collider and ISL bounds are then irrelevant and the radion constraints above are the binding ones. No KK-photon-specific collider bound was collected.

## 6. EC cosmology: spin/torsion in the Friedmann equation

Kranas, Tsagas, Barrow, Iosifidis, PRD 99 (2019), arXiv:1809.10064 (opened): FRW torsion S_abc = 2 phi h_a[b u_c] (their Eq. 13), with S_abc = Gamma_a[bc] = (1/2) T_abc in our notation. This is pure trace torsion: our T_mu = -6 phi u_mu, axial part zero. Primordial He-4 (0.2409 <= Y <= 0.2489) gives -0.005813 < phi/H < 0.019370 (Eq. 54). At T ~ 1 MeV, H ~ 4.5e-16 eV, so |T_0| <~ 5e-17 eV.

This does **not** bound the axial S_0 that a CFJ-type mechanism needs. Minimal EC with Dirac fermions sources only axial torsion (S proportional to kappa times the axial current), whose effect on the Friedmann equation is a kappa^2 (spin density)^2 term. No paper giving a numerical BBN bound on that term was found in this task; it is listed as an open item.

## Unverified or missing

- Data Tables rows: primary papers not opened except Brown 2010 and Allmendinger 2014 (abstracts) and Heckel 2008 (abstract, which carries no b numbers).
- CL missing for Dolan 2022, Payez 2015, Rogers 2023, Kranas 2019.
- Planck 2018 X isocurvature, SPT-3G oscillation bound, Fedderke washout original, Cassini PPN gamma, BBN on alpha and on axial spin density: not collected or not verified.
- Galactic-pole constants in `data/bounds_inputs.csv` are standard values not opened here; they only set v_perp at the 1% level.
