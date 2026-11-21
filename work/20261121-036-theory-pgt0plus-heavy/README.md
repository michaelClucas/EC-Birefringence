# 036: EC-06d heavy-0⁺ (Shie-Nester-Yo PGT) branch: mass, matter coupling, endpoints, verdict

Iteration 11, theorist. Answers BLOCKER 1 of skeptic 035 and the STATE "Live / open" EC-06d item.

## Question
In SNY PGT (0805.3834) the 0⁺ torsion trace is T_μ = -(3/2)∂ln Y, Y = R + 6μ/b, and a CFJ term p = ξT is gauge invariant. In the heavy (GR-compatible) branch:
1. What is m_0⁺ in terms of the couplings, does the mode couple to matter, and what mass range do the GR background, fifth-force data and published bounds allow?
2. Is the observer endpoint local or cosmological, what are β for ground, L2 and cosmic observers, and do ACT and Planck/WMAP exclude it?
3. Is it the Li-Xia-Li ∂ln R coupling (hep-ph/0611192)?
4. Verdict.

## Method
- Reduced SNY to a scalar-tensor theory by computing the Riemann-Cartan scalar for pure-trace torsion in our conventions and integrating out the algebraic torsion (sympy, `check_rc_trace.py`, `check_scalar_mode.py`). The mass is computed three ways: Einstein frame, linearized Jordan trace equation, and SNY's own FRW equations plus their Eq. 53 period. Matter coupling comes from the Brans-Dicke ω and from the Einstein-frame conformal coupling. The photon coupling g follows from the canonical normalization.
- Numerics in `endpoint_bounds.py`: Yukawa-smoothed observer densities for ground, L2 and cosmic-mean observers versus λ_C, the calibration endpoints of ACT (optics model) and Minami-Komatsu (Galactic-foreground EB), the Cassini and Eöt-Wash fifth-force bounds, and the HB photon-coupling bound.
- All numbers are loaded from data/. New file `data/pgt0plus_inputs.csv` holds the Eöt-Wash range (2002.11761, abstract), constants, local densities (copied from 035's scratch inputs with sources) and two flagged task assumptions.
- `run_checks.py` runs everything (36 PASS, 0 FAIL, about 30 s; `run.log`).

Derivation: `derivation_pgt0plus.md`.

## Result
1. **Mass:** m²_0⁺ = -2a0μ/(a1 b) = 2M_P²(a1 + M_P²)/(a1 b). Three independent routes agree, the formula reproduces SNY's Eq. 53 period, and the a1 → ∞ limit is the Starobinsky scalaron 2M_P²/b. Y_c = 6μ/b = 3(a1/M_P²)m². No-ghost a1 b > 0 is recovered.
2. **Matter coupling:** none directly. The trace does not couple to minimally coupled Dirac, scalar or gauge fields. Through metric mixing (the F(χ)R term) the mode couples universally to the matter energy-momentum trace with Yukawa strength **α = μ/(3a1) ≥ 1/3**. If light in the solar system it gives **|γ - 1| > 1/2**. There is no screening: s = κρ/Y_c ≤ 10⁻¹⁵ wherever the mass is allowed.
3. **Allowed mass:**
   - (a) The GR-compatible background needs only m ≳ 15H(z*). 025's δ = 1% example has m = 6×10⁻²⁸ eV, λ_C ≈ 10 kpc.
   - (b) Cassini excludes λ_C ≳ AU by 3.9 orders, which includes that example. Eöt-Wash (|α| = 1 excluded above 38.6 μm; α = 1/3 taken as λ_C < 100 μm) leaves **m_0⁺ ≳ 2×10⁻³ eV**.
   - (c) No other published SNY/PGT bound was found. SNY's own cases are the dead O(H) branch.
4. **Endpoint:** local for every allowed mass. ACT DR6 calibrates with an optics-model prior, so its β includes the observer endpoint. MK-type Planck/WMAP β is calibrated on Galactic-foreground EB, so the observer endpoint cancels.
   - Predicted β_ACT/β_MK = (ρ̄_obs - ρ_*)/(ρ̄_dust - ρ_*) is negative for every λ_C < 20 pc, and ≤ -3×10¹⁷ for λ_C ≤ 100 μm.
   - Measured: +0.63 ± 0.28 (0.215/0.342).
   - The endpoint kill (λ_C < 20 pc) and the Cassini kill (λ_C ≳ AU) overlap, so every mass is excluded.
   - The β/ξ table for ground, L2 and cosmic observers is in run.log (section C2). Ground and L2 differ by up to 10⁴-10⁵ for λ_C between 1 m and R_E. Both are dominated by the analyzer film for λ_C < 1 m.
5. **Size:** at the HB bound g < 6.6×10⁻¹¹ GeV⁻¹ and λ_C ≤ 100 μm:
   - ACT-type observer piece: |β| ≤ 7×10⁻²⁷ rad, 23.7 orders short.
   - Universal (MK-measured) piece: |β| ≤ 2×10⁻⁴⁴ rad, 41.4 orders short.
   - ξ: the universal piece needs ξ ≈ 4×10⁴⁹; HB allows ξ < 2×10⁸.
   - Raising α instead needs α ≥ 10⁴⁷, above strong-interaction strength.
6. **Li-Xia-Li:** not the same model. In the heavy branch ln Y ≈ const + R/Y_c is linear in R and suppressed by s ≪ 1. It equals ln R only for R ≫ Y_c, the torsion-dominated branch already dead (025, 032). The CLAIMS EC-06d / 025 wording "heavy branch = published ∂ln R coupling" should be corrected.

## Verdict
**EC-06d heavy-0⁺ branch: dead** (candidate mechanism excluded). There are three independent kills:
- Cassini for λ_C ≳ AU: 3.9 orders.
- The ACT-vs-MK endpoint sign and ratio for λ_C < 20 pc: predicted ratio ≤ -3×10¹⁷ in the fifth-force-allowed range, measured +0.63.
- Magnitude at the HB photon bound in the allowed range: 23.7 orders short for the ACT observer piece, 41 orders for the universal piece.

With the O(H) branch (025/028/032), EC-06d is dead in full. No Result C cost list is needed.

## Open issues
- EK22 (2205.13962) is taken to use the Minami-Komatsu foreground-EB calibration, as MK20 does (repo row: abstract only). MK20 alone gives the same conclusion (ratio 0.61).
- The α = 1/3 Eöt-Wash range is not in the abstract. 100 μm is an assumption, a factor 2.6 above the |α| = 1 value, worth 0.8 orders in β. Yukawa bounds between 100 μm and an AU are not loaded into data/. They are not needed, because the endpoint kill covers that range.
- The analyzer-film model (0.1 μm silicon, ρ t/(2λ)) is illustrative. The kills have ≥ 15 orders of margin against it.
- Linear theory only around the GR background. Nonlinear SNY solutions inside dense matter are not considered; s ≤ 10⁻¹⁵ suggests they are irrelevant.
- The overall sign of β differs between 025 (β = -(3ξ/4)ln(Y_obs/Y_src)) and 035's script. Ratios and magnitudes are unaffected.
- Li-Xia-Li's own ∂ln R model (no torsion) is still environment-dependent at the observer (025 Sec. 2c1) and has not been killed. It is out of scope for the torsion question.
