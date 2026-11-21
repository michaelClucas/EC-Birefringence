# EC-06d heavy-0⁺ branch: derivation

Conventions: CONVENTIONS.md (signature -+++, derivative-index-first connection, T_μ = T^λ_{μλ}, κ = 1/M_P²). SNY = Shie-Nester-Yo arXiv:0805.3834. Where units M_P = 1 are used, a0 = -M_P² = -1.

## 1. Reduction of SNY to a scalar-tensor theory [check_rc_trace.py A1-A5, check_scalar_mode.py B1-B3]

SNY Lagrangian: L = -(a0/2)R + (b/24)R² + (a1/8)(T_{νσμ}T^{νσμ} + 2T_{νσμ}T^{μσν} - 4T_μT^μ), with R the Riemann-Cartan scalar. With only the 0⁺ mode excited (SNY Eqs. 21-23: axial and tensor torsion vanish without spin sources) the torsion is pure trace, T^λ_{μν} = (1/3)(δ^λ_ν T_μ - δ^λ_μ T_ν). Computed in our conventions:

- R = R_LC + 2∇·T - (2/3)T_μT^μ (A3 flat with generic T_μ(x); A4 FRW spot check of covariance; A1 sign of R_LC).
- SNY torsion-square combination on pure trace = -(8/3)T², so the a1 term is -(a1/3)T² (A5; SNY index order T_{ij}{}^k mapped to ours).

Auxiliary field: (b/24)R² = (b/12)χR - (b/24)χ², χ = R on shell. With Ω = -a0/2 + bχ/12, the torsion enters algebraically: after integrating 2Ω∇·T by parts, L_T = -2∂Ω·T - ((2Ω + a1)/3)T². Varying T:

  **T_μ = -3∂_μΩ/(2Ω + a1) = -(3/2)∂_μ ln Y,  Y = χ + 6μ/b,  μ = a1 - a0** (B1).

This is SNY Eq. 20 (Y = R + 6μ/b), recovered independently, including SNY's μ = a1 - a0. On shell L_T = (b/8)(∂χ)²/Y (B2), so the theory is

  L = ½[F R_LC - K(∂χ)²] - V + L_m(g),  F = -a0 + bχ/6,  K = -(b/4)/Y,  V = (b/24)χ².

Einstein-frame kinetic coefficient Z = K/F + (3/2)(F′/F)² = a1 b/(4F²Y): positive iff a1 b > 0, which is SNY's no-ghost condition (B3).

## 2. Mass of the 0⁺ mode [B4-B7]

Three routes agree:
- Einstein frame (a0 = -1): m² = U″(0)/Z(0) with U = V/F² gives m² = 2μ/(a1 b) (B4).
- Linearized Jordan trace equation, general a0: □χ - m²χ = (2μ/(a1 b))T_m with **m²_{0⁺} = -2a0μ/(a1 b) = 2M_P²(a1 + M_P²)/(a1 b)** (B5). Static heavy limit: χ = -κT_m = κρ, the GR value of R (B5c).
- SNY's own FRW equations (their Eqs. 43-47) linearized at H = 0 give R̈ = -m²R - (2μ/(a1 b))T_m with the same m², and SNY Eq. 53's oscillation period 2π√(-a1 b/(2a0μ)) is exactly 2π/m (B6).
- a1 → ∞ (torsion frozen): m² → 2M_P²/b, the Starobinsky scalaron of (M_P²/2)(R + R²/(6m²)) (B7).

Useful relation: Y_c ≡ 6μ/b = 3(a1/M_P²)m² (B10a). The "heavy branch" s = R/Y_c ≪ 1 is m² ≫ R/(3a1′), a1′ = a1/M_P².

## 3. Matter coupling [B8]

Minimally coupled Dirac fields couple only to the axial torsion S^μ; scalars and gauge fields (F = dA) do not couple to torsion. So the trace mode has **no direct coupling to ordinary matter** (as SNY state). It couples **indirectly through the metric**: matter sits in the Jordan frame, and the 0⁺ field mixes with the metric through F(χ)R_LC. That is a universal coupling to the matter energy-momentum trace T_m, as for the f(R) scalaron:

- Brans-Dicke form: ω = KF/F′² = -3/(2μ) (a0 = -1).
- Yukawa strength **α = 1/(3 + 2ω) = μ/(3a1) = (1 + 1/a1′)/3 ≥ 1/3** for every allowed a1 > 0. Einstein-frame cross-check α = 2Q², Q = F′/(2F√Z) (B8c).
- If the mode is light on solar-system scales: γ - 1 = -2μ/(4μ - 3), **|γ - 1| > 1/2** for every a1 > 0 (B8d-e).

There is no chameleon screening. In any environment reachable in the allowed mass range, s = κρ/Y_c ≤ 10⁻¹⁵, so the mass and the coupling are their vacuum values.

## 4. Photon coupling and β [B9-B10]

CFJ p = ξT = -(3ξ/2)∂ln Y. With p = -g∂φ (CONVENTIONS) this is g φ = (3ξ/2)ln Y, an ALP-type coupling of the canonical 0⁺ field φ = √Z δχ with **g = ξ√(3/(2μ′a1′))/M_P** (B9). Rotation (025 sign): β = -(3ξ/4)ln(Y_obs/Y_src). Skeptic 035's script uses the opposite overall sign; that does not affect any ratio or magnitude. In the heavy limit:

  β ≈ -(3ξ/4)κ(ρ̄_obs - ρ̄_src)/Y_c,  |β| = √(2α) g |Δρ̄| /(4 M_P m²) (B10b),

where ρ̄ is the density Yukawa-smoothed over λ_C = 1/m. The second form shows what the branch is: a gravitational-strength scalar with an ALP photon coupling, sourced by matter.

## 5. Allowed mass range [endpoint_bounds.py C1]

- λ_C ≳ AU: the mode is effectively massless in the solar system, |γ - 1| ≥ 1/2 against Cassini |γ - 1| < 6.7×10⁻⁵ (2σ, data/kk_bounds_inputs.csv). Excluded by 3.9 orders. This includes 025's GR-compatible ξ ~ 1 example (δ(z*) = 1%, μ = 2): m = 6×10⁻²⁸ eV, λ_C ≈ 10 kpc.
- Short range: α ≥ 1/3 unscreened Yukawa. Eöt-Wash (2002.11761) excludes |α| = 1 for λ > 38.6 μm. α = 1/3 is not quoted, so λ_C < 100 μm (m_0⁺ ≥ 2×10⁻³ eV) is used as a generous stand-in. The range between 100 μm and an AU is covered by other Yukawa searches (not loaded into data/), and independently by Sec. 6 below.
- Published SNY/PGT bounds: none found beyond SNY itself. SNY's Table I cases have m ~ H₀, the O(H) branch already killed (025, 032).

## 6. Endpoints [C2]

β depends only on ρ̄ at the two endpoints. Two calibration schemes matter:
- **ACT DR6** (2509.13654) fixes the instrument angle with an optics-model prior. Its β therefore contains the observer endpoint ρ̄_obs, the smoothed density at the polarization analyzer.
- **Planck/WMAP Minami-Komatsu** (MK20, EK22) fix the angle with Galactic-foreground EB. The foreground photons share the observer endpoint, so that term cancels, and β_MK ∝ ρ̄_dust - ρ_* ≈ -ρ_* (ρ_dust ~ 5×10⁻⁴ ρ_*).

Predicted β_ACT/β_MK = (ρ̄_obs - ρ_*)/(ρ̄_dust - ρ_*). ρ_* = ρ_m(z*) = 3.5×10⁻¹⁸ kg/m³. For a ground observer, ρ̄_obs is set by the analyzing film (λ < 1 m), the air (≥ ρ_air/2), Earth, or the Sun. It stays above ρ_* for every λ_C < 20 pc. The predicted ratio is negative there, and |ratio| ≥ 3×10¹⁷ for λ_C ≤ 100 μm. The measured ratio is +0.63 ± 0.28. **The endpoint kill covers λ_C < 20 pc and Cassini covers λ_C ≳ AU. Together they exclude every mass.** An L2 observer (Planck's own instrument endpoint) is equally local; the table in run.log gives ground, L2 and cosmic-mean β/ξ at a1′ = 1.

## 7. Size at the photon-coupling bound [C3]

At λ_C ≤ 100 μm with g < 6.6×10⁻¹¹ GeV⁻¹ (HB, 1406.6053, valid to ~10 keV) and α = 1/3:
- ACT-type observer piece (analyzing film, 0.1 μm silicon assumed): |β| ≤ 7×10⁻²⁷ rad, **23.7 orders** short of 3.75×10⁻³ rad.
- Universal piece (what MK measures): |β| ≤ 2×10⁻⁴⁴ rad, **41.4 orders** short.
- In ξ: the universal piece needs ξ ≈ 4×10⁴⁹ (a1′ = 1), while HB allows ξ < 2×10⁸ (√(a1′μ′) more for large a1′, with β/ξ falling by the same 1/a1′).
- Raising α (a1′ → 0) gains √α, but closing the gap needs α ≥ 10⁴⁷ at 100 μm, and more at shorter range. That exceeds the nucleon strong/gravity ratio (1.7×10³⁸).

Since β ∝ 1/m², heavier masses only lower β.

## 8. Relation to Li-Xia-Li hep-ph/0611192

LXL couple ∂_μ ln R (metric Ricci scalar, no new field) through the B-L current and the anomaly. Here the photon couples to ln Y = ln(Y_c + R). In the heavy (GR-compatible) branch, ln Y ≈ const + R/Y_c: the coupling is **linear in R/Y_c**, suppressed by s ≪ 1, not ln R. The two coincide only for R ≫ Y_c, the torsion-dominated SNY histories killed by H(z*) (025, 032). So the identification "heavy-0⁺ = published ∂ln R coupling" in CLAIMS EC-06d / 025 is wrong in the heavy limit. The heavy branch is a gravitationally coupled scalar with a hand-added ALP coupling. LXL's own model also has an environment-dependent observer endpoint (025 Sec. 2c1); it is not a torsion route and was not re-examined here.

## 9. Lab k_AF

For λ_C ≤ 100 μm, p = -(3ξ/2)∇ln Y is non-zero only within λ_C of density changes, and its line integral depends only on endpoints. The homogeneous k_AF bounds in data/kaf_lowz_bounds.csv (radio galaxies, the geomagnetic spatial k_AF) assume a uniform background and do not apply directly. They are not needed, since Secs. 5-7 already exclude the branch.
