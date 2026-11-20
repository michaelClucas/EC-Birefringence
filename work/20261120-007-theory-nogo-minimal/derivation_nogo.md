# Derivations: EC-02, EC-05, EC-04 (task 20261120-007)

Conventions: `CONVENTIONS.md`. Every finding is a candidate mechanism. The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Inputs: `data/cosmo_inputs.csv` (new), `data/birefringence_measurements.csv`, `work/20261120-005-bounds/derived_bounds.csv`.

Settled in 003 and used here: CFJ term L = -½ p_μ A_ν F̃^{μν} gives β = +½∫p_τ dτ; with p = ξS (EC-01), β = -(ξ/2)∫S^0̂ dt (cosmic time, last scattering to today); EC + Dirac gives L ⊃ (3/4) S_μ J5^μ and S^μ = -(κ/2) J5^μ.

## A. EC-02: minimal EC spin source + CFJ coefficient ξ

### A0. Formula
S^0̂ = -(κ/2) n5 with n5 ≡ J5^0̂, so

  β = +(ξκ/4) ∫ n5 dt.

Only the isotropic, time-like n5 enters. A spatial spin density gives a spatial S^i and a direction-dependent rotation, not an isotropic β.

### A1. n5 is helicity times speed (check_helicity.py)
For a plane wave of momentum p and helicity h: u†γ5u / u†u = h|p|/E = h v (exact; = h for m = 0). The axial current is C-even (C(γ^μγ5)^T C⁻¹ = +γ^μγ5), so an antiparticle of helicity h also contributes +hv. Helicity commutes with the free FRW Dirac Hamiltonian α·p + a m β for every a, and with γ5, so free streaming and the axial-torsion coupling both preserve it. Hence

  n5 = Σ_species Σ_h h ∫ d³p/(2π)³ f_h(p) p/E.

Relic neutrinos were produced relativistic as ν (h = -1) and ν̄ (h = +1). For Dirac neutrinos the sterile helicities are unpopulated; for Majorana neutrinos "ν" and "ν̄" are the two helicities of one field. Both give n5 = -Σ_i ∫ [f(ξ) - f(-ξ)] (p/E_i), i.e. the lepton asymmetry times the mean speed of each mass state. Today relic neutrinos with m ≥ 0.0087 eV are non-relativistic, so n5 today is suppressed by ⟨p⟩/m. At last scattering T_ν ≈ 0.18 eV, so all states are still relativistic there, and the time integral is dominated by z ≳ 100 (97% of ∫a⁻³dt), so the suppression barely matters for β (cases nu-A and nu-B differ by 2%).

Difference from the published neutrino-current operators (Geng-Ho-Ng 0711.4617; Zhou-Huang-Geng 2302.11140; Reig 2611.07578): those couple a vector or left-handed current j^μ, whose j^0 contains the number density Δn and is not velocity-suppressed. Torsion couples only to the axial part.

### A2. Inputs (data/cosmo_inputs.csv)
- ξ_ν = 0.024 ± 0.012 (68%), BBN+CMB+BAO, flavours equilibrated: 2409.08280 Table 4. Used |ξ_ν| ≤ 0.048 (mean + 2σ).
- Generous case without flavour equilibration: ξ_μ = ξ_τ saturating N_eff ≤ 2.99 + 2×0.17 (1807.06209 abstract) with ΔN_eff = (30/7)(ξ/π)² + (15/7)(ξ/π)⁴ per species: ξ_μ,τ ≤ 0.57.
- Masses (case nu-B): normal ordering, m1 = 0, Δm² from NuFit-6.0 (2410.05380 Table 1).
- Cosmology: Planck 2018 Table 2 (1807.06209); T_CMB (0911.1955); N_eff^SM = 3.044 (2012.02726).

### A3. Other sources of isotropic n5
- Electrons and baryons: no process gives them a net helicity, and thermal equilibrium has ⟨h⟩ = 0. Envelope (every electron the same helicity, thermal speed, n_e ≤ n_b): 39 orders short.
- Magnetic-field spin polarisation: P = μ_B B/T ≈ 1×10⁻¹⁰ at last scattering for the Planck bound B_1Mpc < 4.4 nG (1502.01594). This polarises spin, not helicity (⟨σ·p̂⟩ = 0 over an isotropic momentum distribution), so its isotropic n5 is 0. The spatial S^i ≲ 3×10⁻⁷⁷ eV, random in direction between coherence patches.
- Chiral dark matter. A fully spin-aligned population is not isotropic. The maximal isotropic case is a fully helicity-polarised fermion with n5 = n_DM⟨v⟩. The speed is bounded by free streaming: thermal-relic mass ≥ 5.3 keV (1702.01764, 2σ) gives ⟨v⟩₀ = 1.3×10⁻⁸ (1.4×10⁻⁵ at LSS). A degenerate, Pauli-limited one-helicity gas with the same ⟨v⟩ bound needs m ≥ 3.3 keV. n5 ∝ a⁻⁴. Both cases come out about 35 orders short.

### A4. Result (ec02_nogo.py, ξ = 1)

| source | n5(LSS) [eV³] | \|S^0̂\|(LSS) [eV] | max β [deg] | shortfall [orders] | ξ needed | Λ_eff = M_P√(2/ξ) |
|---|---|---|---|---|---|---|
| ν asymmetry, \|ξ_ν\| ≤ 0.048, massless (envelope) | 1.5×10⁻⁴ | 1.2×10⁻⁵⁹ | 7.5×10⁻³⁰ | 28.5 | 2.9×10²⁸ | 2.0×10⁴ GeV |
| same, NO masses m1 = 0 | 1.5×10⁻⁴ | 1.2×10⁻⁵⁹ | 7.3×10⁻³⁰ | 28.5 | 2.9×10²⁸ | 2.0×10⁴ GeV |
| ν, ξ_μ,τ ≤ 0.57 from N_eff only (generous) | 1.3×10⁻³ | 1.1×10⁻⁵⁸ | 6.3×10⁻²⁹ | **27.5** | 3.4×10²⁷ | 5.9×10⁴ GeV |
| helicity-polarised thermal WDM, 5.3 keV | 3.3×10⁻¹¹ | 2.8×10⁻⁶⁶ | 1.0×10⁻³⁶ | 35.3 | 2×10³⁵ | 7 GeV |
| helicity-polarised degenerate DM | 5.4×10⁻¹¹ | 4.5×10⁻⁶⁶ | 1.6×10⁻³⁶ | 35.1 | 1×10³⁵ | 9 GeV |
| all electrons one helicity (envelope) | 2.8×10⁻¹⁵ | 2.4×10⁻⁷⁰ | 1.1×10⁻⁴⁰ | 39.3 | 2×10³⁹ | 0.08 GeV |

Today |S^0̂| from relic neutrinos is ≲ 10⁻⁶⁸ eV, so the lab bound 7.3×10⁻²² eV is irrelevant. Cross-check: our time integral, applied to Reig's operator with ξ_νe = 0.05, gives Λ* = 12 TeV against Reig's 9.8 TeV (Eq. 3; their fit is to the Yin-Komatsu β, not exactly 0.215°).

### A5. What bounds ξ
1. **Gauge invariance.** δ(-½ξS_μA_νF̃^{μν}) = ¼ξχ(∂_νS_μ - ∂_μS_ν)F̃^{μν} + total derivative (same structure as Reig Eq. 4). In EC, S = -(κ/2)J5 is algebraic and not closed off the FRW background, so the term breaks U(1) as an operator. ξ is not a free parameter of a consistent EC theory unless the term is completed.
2. **Dimension.** With S = -(κ/2)J5, the ξ-term is the dimension-6 operator (ξκ/4) J5_μ A_ν F̃^{μν}, scale Λ_eff = M_P√(2/ξ). EC fixes Λ_eff ≈ M_P at ξ ~ 1. The needed ξ ≈ 3×10²⁷ to 3×10²⁸ means Λ_eff ≈ 2-6×10⁴ GeV, i.e. the neutrino-current operator of 0711.4617 / 2302.11140 with a new scale. At that point torsion contributes nothing; it only fixes M = M_P. KRT (0712.4393, photon-sector paragraph) also note that the leading gauge-invariant torsion-photon couplings are dimension 5 and need their own scale.
3. **UV completion (Reig 2611.07578).** The gauge-invariant completion is a Kalb-Ramond field dual to a massless axion with g_aγ/f = 1/Λ² (Eq. 9). For a universal axial current, the electron coupling is g_ae = 2m_e/f, so g_aγ g_ae = ξκm_e. The helioscope bound g_aγg_ae < 8.1×10⁻²³ GeV⁻¹ (quoted by Reig) gives **ξ < 9.4×10¹⁷**. The most generous needed ξ, 3.4×10²⁷, is excluded by 9.6 orders. The axial current is not conserved (mass term, anomaly), so Reig's "conserved anomaly-free current" escape does not exist here.

## B. EC-05: fermion-loop CFJ from axial torsion

Coefficients converted to |p| = r|S| in our normalisation (ec05_loop.py, symbolic):
- Dobado-Maroto astro-ph/9706044: coupling -(1/8)ψ̄γ5 S̸_DM ψ with S_DM = εT = 6S (Eq. 4). L_eff = -(α/24π)ΣQ² ε^{μναβ}F_{μν}A_αS_β (Eq. 11), from dimensional regularisation at leading order in K²/m² (Eq. 10). This gives **r_DM = (α/π)ΣQ²**.
- Jackiw-Kostelecky hep-ph/9901358: coupling -b̸γ5 (Eq. 1), k = (3/16π²)e²b (Eq. 16, "nonperturbative"), |p| = 2|k|, |b| = (3/4)|S|. This gives **r_JK = (9α/8π)ΣQ²** = (9/8) r_DM, the same as 003's C = 3/16π².
- Coleman-Glashow routing and Pauli-Villars give r = 0 (JK Eq. 17 and text). JK conclude that perturbation theory leaves r undetermined.

With ΣQ² = 8 (all SM charged fermions, no decoupling; DM text under Eq. 12), r_JK = 0.021 and r_DM = 0.019. EC-05 is EC-02 with ξ → r:

| S source | r | max β [deg] | shortfall |
|---|---|---|---|
| ν asymmetry, \|ξ_ν\| ≤ 0.048 | JK 0.021 | 1.6×10⁻³¹ | 30.1 |
| ν, N_eff-only (generous) | JK 0.021 | 1.3×10⁻³⁰ | **29.2** |

A constant background would need |S^0̂| ≈ 5×10⁻³⁴ eV through the loop (ΣQ² = 8).

**Reig applied properly.** The loop operator is (rκ/4) J5·A F̃ with Λ_eff = M_P√(2/r) ≈ 2.4×10¹⁹ GeV. Reig's dual-axion plus helioscope argument, applied to a universal axial current, requires Λ > 3.6×10⁹ GeV. EC-05 satisfies it easily. Reig's bound does not kill EC-05; its size does, at 29-30 orders. Reig would bind only if something lowered Λ to the 10⁴ GeV the signal needs, and then it excludes that by about 10 orders (A5.3).

**The regularisation ambiguity does not survive on FRW.** A homogeneous S_μ = (S_0(t), 0) satisfies dS = 0, so S = ∂θ. A gradient insertion carries momentum q ~ H, and its amplitude is fixed by the anomalous Ward identity: the derivative coupling decouples as B1(τ) ≈ -q²/(12m_f²) (series checked symbolically; Bauer-Neubert-Thamm 1708.00443 loop function). This gives an extra suppression (H_LSS/m_e)²/12 ≈ 4×10⁻⁷⁰ relative to the JK value. The JK/DM nonzero values are the strictly constant-b, zero-momentum definition. The "maximum" quoted above uses them anyway.

## C. EC-04: dynamical axial torsion (Poincaré gauge theory)

**Which sectors propagate a parity-odd scalar.** PGT has two healthy spin-0 connection modes, 0⁺ and 0⁻ (Chen, Ho, Nester, Wang, Yo 0908.3323, with earlier analyses cited there). On FRW, the 0⁻ mode is the axial torsion S_0(t), and it couples to the fermion axial spin density (their Eq. 11). The 0⁻ dynamics come from the pseudoscalar-curvature-squared term (b₋/24)E² (their Eq. 7). Its mass parameter is m₋ ≡ A₀ + 2A₃, the EC term plus the axial-torsion-squared coefficient. At late times it is a harmonic oscillator with ω₋² ∝ m₋/b₋ (their Eq. 64; the sign was not legible in our text extraction). Baekler-Hehl-Nester 1009.5112 generalise this with parity-odd couplings (RX, VA terms). A propagating 1⁻/1⁺ axial vector would need an (∂_{[μ}S_{ν]})² kinetic term. A vector cannot propagate its spin-1 and spin-0 parts both healthily, and Shapiro (hep-th/0103093) argues propagating axial-vector torsion must be heavier than the fermions it couples to. So the only light candidate is 0⁻.

**Reduction (check_pgt_0minus.py, all symbolic).** For pure axial torsion on flat space, the pseudoscalar curvature is X = ε^{μνρσ}R_{μνρσ}(Γ) = -6 ∂_μS^μ exactly, with no S² term. The axial sector is therefore

  L = (A/2) S_μS^μ + (w/2)(∇·S)² + (3/4) S_μ J5^μ,

where A = 3/(2κ) in pure EC (shifted by PGT torsion-squared terms) and w ∝ b₋. An auxiliary field Φ eliminates S, with S = (∂Φ - (3/4)J5)/A, and leaves

  L = -(1/2A)(∂Φ - ¾J5)² - Φ²/(2w).

This is one scalar, ghost-free iff A > 0 and non-tachyonic iff w > 0, with m² = A/w. Consistently, the S propagator (Aη + wkk)⁻¹ has a single pole at k² = -A/w. Setting w = 0 recovers S = -(κ/2)J5.

**The 0⁻ mode is an axion.** With φ = Φ/√A: S_μ = ∂_μφ/√A, and the fermion coupling is (1/2f)∂φ·J5 with universal f = (2/3)√A = 2.0×10¹⁸ GeV (EC value of A). That gives g_ae = m_e/f = 2.6×10⁻²².
- The shift symmetry φ → φ + c is exact (only ∂φ appears), so QCD generates no potential.
- There is no tree-level photon coupling. The loop through charged fermions decouples as (m_φ/m_f)²/12, giving g_eff ≈ 10⁻⁹⁹ GeV⁻¹ for m_φ ~ H₀. **g_φγ ≈ 0 for the pure PGT 0⁻ mode.**
- Adding the EC-01 term -½ξS·K (gauge invariant for S = ∂φ/√A, by ∂_μK^μ = ½FF̃, checked) gives **g_φγ = -ξ/√A ≈ -0.82 ξ/M_P = 3.4×10⁻¹⁹ ξ GeV⁻¹**, and β = -(ξ/2)Δθ with θ = φ/√A.

**FRW.** φ̈ + 3Hφ̇ + m²φ = (3/(4√A)) ∇_μJ5^μ ≈ 0, since ∇·J5 = a⁻³d(a³n5)/dt is nonzero only through the slow relativistic-to-NR transition of relic neutrinos. S^0̂ = ∓φ̇/√A. The signal needs |ξΔθ| = 2β = 7.5×10⁻³. Two ways to get it (ec04_estimates.py):
- Thawing field with m ≈ H₀ to 2H₀: θ_i ≈ 0.06/ξ to 0.016/ξ, Ω_φ ≈ 10⁻³/ξ², today |S_0| ≈ 2×10⁻³⁵/ξ eV. This is the standard "Planck-coupled ALP with Δφ ≈ 2×10¹⁶ GeV/ξ" picture.
- Massless kinetic mode: θ̇ ∝ a⁻³, energy fraction 3×10⁻⁵/ξ² at last scattering.

**Bounds.**
- Lab: |S_0| < 7.3×10⁻²² eV gives only ξ ≳ 10⁻¹⁴.
- CAST g < 5.8×10⁻¹¹ GeV⁻¹ gives ξ < 1.7×10⁸.
- g_ae ~ 10⁻²² is far below the red-giant bound.
- Naturalness: m ≲ H₀ needs w = A/m² ≳ 4×10¹²⁰. With O(1) PGT coefficients, the 0⁻ mass is ~M_P and S averages to zero. Making A small instead lowers f. The red-giant bound g_ae < 1.6×10⁻¹³ (quoted in 2611.07578, ref. [38]) needs f > 3×10⁹ GeV, so √A ≳ 5×10⁹ GeV and still w ≳ 10¹⁰³ for m ≲ H₀.

