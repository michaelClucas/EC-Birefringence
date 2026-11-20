# EC-06a (Magueijo-Złośnik) verdict: derivation

Conventions: `CONVENTIONS.md`. Everything here is a candidate mechanism. The ACT DR6 β = 0.215° ± 0.074° is 2.9σ with systematics caveats and may not be cosmological. Checks: `check_offfrw.py` (A-E), `check_dynamics.py` (T1-T6), `mz_numbers.py` (N1-N5); `run.log` 36 PASS, 0 FAIL.

## 0. Theory and photon coupling

MZ (1908.05184, Eq. 7, ar5iv HTML opened): Palatini action plus -(3/2Λ)(ε_abcd + (2/γ)η_ac η_bd)(R - Λee/3)(R - Λee/3) - (2/γ)T^aT_a, with Λ a field without kinetic term and matter not coupled to the connection. Field equations (opened):

- Eq. 14 (tetrad): ε_abcd(e^b R^cd - ⅓Λ e^b e^c e^d) = -2κτ_a, the Einstein equation with the full curvature.
- Eq. 15 (connection): T^[a∧e^b] = -(3/2Λ²) dΛ∧R^ab + (3/4γΛ²) ε^abcd dΛ∧R_cd.
- Eq. 16 (Λ): ε_abcd(R^ab R^cd - (Λ²/9) e^a e^b e^c e^d) + (2/γ) R^ab R_ab = 0.

Write θ ≡ 3/(2Λ). Then dθ = -(3/2Λ²)dΛ, and Eq. 15 reads T^[a∧e^b] = dθ∧[R^ab - (1/2γ) ε^ab_cd R^cd]. The photon term is added by hand: -½ξ S_μ A_ν F̃^{μν}.

- A1-A2: our forms code rebuilds MZ's torsion from ω^i_0 = g e^i, ω^ij = -Pε_ijk e^k. On FRW, Eq. 15 reproduces MZ Eqs. 28 and 29 for exactly one ε sign, which is ε^0123 = -1, the same as ours. So the transcription and the conventions agree.

## 1. Off-FRW closure (item 1)

**Exact algebraic solution for S.** Eq. 15 is linear in the torsion once the curvature is given. B1-B2 use random rational curvature, with no symmetry beyond antisymmetry in each index pair. For any v = dθ, the axial part of the solution is

S_A = (1/3γ) E_A^F v_F - (2/3) W_A^F v_F.

Here E_A^F is the dual of ε_ABCD e^B∧R^CD, which is exactly the Eq. 14 combination. W is the dual of R^A_B∧e^B, which equals the dual of DT^A by the Bianchi identity (C1). Eq. 14 fixes E = 2𝒢 with 𝒢_μ^ν ≡ κτ_μ^ν - Λδ_μ^ν. The normalisation is checked against MZ Eq. 29 on FRW (A3). Hence

**S_μ = (2/3γ)(κτ_μ^ν - Λδ_μ^ν)∂_νθ - (2/3)(⋆DT)_μ^ν∂_νθ.**

The second term is quadratic in ∂θ, since T itself is linear in ∂θ. On FRW it is the 6gPΛ̇/Λ² piece of Eq. 29, which is O(γ²) relative to the first term on the tracker.

**Closure.**

- The Λ piece is -(2Λ/3γ)∂θ = (1/γ)∂ ln Λ, a gradient. **In vacuum S is closed at leading order.**
- The matter piece (2κ/3γ)τ·∂θ is not closed in general. For a perfect fluid at linear order in Newtonian gauge (D1): (dS)_{τi} = ∂_i Y, with
  Y = (2κ/3γ)[(ρ̄+p̄)δθ_c′ + p̄′δθ_c + θ̄′δρ_c], where δθ_c = δθ + θ̄′v and δρ_c = δρ + ρ̄′v.
  Φ and Ψ drop out. Y is gauge invariant and vanishes for a pure time-shift mode (D2). For dust, closure needs δθ_c′ = -θ̄′δ_c at every point and time (D3), so θ's comoving perturbation would have to track the growing density contrast.
- Nothing in the theory enforces that relation. θ (that is, Λ) carries its own initial data: MZ count one extra degree of freedom, and the transport law in Sec. 2 is first order along u, so δθ_c is free on the initial slice. Closure holds only on a measure-zero set of solutions.
- **U(1) breaks wherever there is matter.** In effective-field form (Sec. 2, E1): S ∝ (φ̇, c_s∇φ) with c_s = Λ/(Λ+κρ). Then dS = (c_s - 1)∂_t∂_iφ, which vanishes only where κρ = 0. The breaking factor is 1 - c_s = κρ/(Λ+κρ). It equals Ω_m ≈ 0.32 for the cosmic mean today, about 1 at z*, and about 1 in galaxies and any matter. By 022 G2, Maxwell + CFJ is then inconsistent without a Stückelberg completion. Gravitational CS is not a repair (028).

## 2. Leading-order dynamics and the θ normalisation

- T1: at leading order in γ (P = O(γ), T = O(γ²), Λ̇ = O(γ²)), MZ Eq. 30 is the FRW form of
  **∇_μ(𝒢^{μν}S_ν) = (γ/4)(GB - (8/3)Λ²)**, with GB = 24H²(Ḣ+H²).
  The covariant form follows from Eq. 16. The Pontryagin difference with torsion is exact, 2d(K∧R̊) + O(K²), and for axial K it is ∝ ∇(G̊·S). The axial K drops out of the Euler density at linear order. The tensor-torsion part q is not tracked (open issue 1).
- T2: the tracker is P = γH/3 - (γΛ²/9)∫a³dt/(a³H²), exact for any flat history. So ε_S = 2P/H = 2γ/3 in the radiation and matter eras, and 0.42γ today (N1). The general solution adds P_h = C/(a³H²). T3 reproduces 2001.06373 Eq. 27 (radiation tracker), P² = γ²κρ/27 (matter), and their Eq. 28 (δ_P ∝ a³).
- T4: the θ-dependent action is -(M_P²/2)θ[εRR + (2/γ)RR] - 3M_P²/θ. C2 checks εRR = +GB·vol. Matching the transport equation to the Euler-Lagrange equation of L_θ = -(N/2)Z^{μν}∂_μθ∂_νθ - (M_P²/2)θ GB - 3M_P²/θ, with Z = 𝒢·𝒢, gives **N = 4M_P²/(3γ²)**. For dust plus Λ, Z = Λ²h - (Λ+κρ)²uu. The canonical field φ = √N(Λ+κρ)θ has L = ½φ̇² - ½c_s²(∇φ)², and
  **S = -(1/√3 M_P)[∂φ - (1-c_s)(∂φ)_⊥]**, so γ cancels.
  The gradient piece is an ALP coupling of φ with g = ξ/(√3 M_P), which is 3.5×10⁻¹⁵ GeV⁻¹ at ξ = 1.5×10⁴. The non-gradient piece carries the U(1) breaking.

## 3. Strong coupling (item 2)

- T5: a Stückelberg photon (A → A + ∂π/m) turns the non-gradient piece into a π-φ-γ vertex c·π ∂_t∂_iφ F̃^{i0}, with **c = ξ(1-c_s)/(2√3 M_P m)**. This is the same structure as 028's η vertex, with the Λ-field in place of the f(T) Lorentz mode.
- γγ → φφ by π exchange, or γπ → γπ by φ exchange, has amplitude ~c²E⁴. Partial-wave unitarity, |A|/16π ≤ 1, gives
  **Λ_SC = (16π)^{1/4}(2√3 M_P m/(ξ(1-c_s)))^{1/2}**.
  This is conservative: for c_s < 1 the φ momenta at fixed energy are E/c_s, which lowers Λ_SC further, and that is not used.
- Numbers (N5), with 1 - c_s = 0.315 (the cosmic mean, the smallest value used):

| ξ | branch | m = 10⁻²⁶ eV | 10⁻²⁴ eV | 10⁻¹⁸ eV (PDG) | m for Λ_SC = 10 TeV |
|---|---|---|---|---|---|
| 1.5×10⁴ | tracker, GW speed (026) | 0.36 eV | 3.6 eV | **3.6 keV** | 7.8 eV (10^18.9 × PDG) |
| 2.9×10⁻³ | cancellation branch, N_eff only (best case) | 0.8 keV | 8 keV | **8 MeV** | 1.5×10⁻⁶ eV (10^12.2 × PDG) |
| 4.7×10¹⁶ | tracker, k* condition | 2×10⁻⁷ eV | 2×10⁻⁶ eV | 2 meV | — |

- Above Λ_SC, γγ → (Λ-field quanta) is strongly coupled. TeV-blazar pair production (√s ~ 0.1-2 MeV) and LEP/LHC γγ physics see standard photons. So the theory needs a UV completion of the photon-Λ sector below ~keV-MeV, and none is known.
- 029's ~0.1 GeV is reproduced as the dim-7 arithmetic (m M_P²/ξ)^{1/3} = 0.16 GeV at 10⁻²⁶ eV (73.8 GeV at the PDG limit). It is not the lowest scale. The θ-exchange scale is 0.4 eV - 4 keV over the same range.
- The 025/028 Stückelberg window edges come from dS near compact objects. For MZ those edges depend on the local, nonlinear θ solution, which was not derived. The strong-coupling verdict does not use the window: it fails at the PDG upper limit itself.

## 4. Gravity-sector consistency: the tensor singular scale

- T6, from 2001.06373 Eqs. 44-48 and (III): the coefficient of ℋ̈ in the tensor equation is A/(A²+B²). With Λ̇ from Eq. 29 (exact), A = 0 at **k_P* = (Λ+κρ)/(6P) = H/ε_S**, independent of γ. Only the helicity with sign(k_P) = sign(P) reaches it.
- 2001.06373 (Eqs. 69-70) itself requires k* ≫ k_LIGO.
- At the 026 bound (γ = 1.1×10⁻⁷), GW170817 at 24 Hz sits at k/k* = 3×10¹² (N3). The relation c_T - 1 = γ²/18 (Eq. 72) is an expansion for k ≪ k*, so 026's ξ ≥ 1.4×10⁴ was obtained outside its regime.
- Imposing k* > k_GW today gives ε_S(0) < 1.45×10⁻²⁰, so **ξ ≥ 4.7×10¹⁶**.
- Our reading, not derived from an action: the kinetic coefficient changes sign across k*, so one graviton helicity is ghost-like above k*.

## 5. Off-tracker branches (item 3)

- The homogeneous mode has ε_h ∝ (aH)⁻³. It grows through the radiation and matter eras: ε_h(z*)/ε_h(0) = 1.0×10⁻⁴ (N4a). **There is no decaying branch with ε(z*) ≫ ε(0).** 026's "P ∝ a⁻³" is the Λ-era behaviour only. In the matter era P_h is constant and ε_h ∝ a^{3/2}.
- The only way to make ε small today while it is large at z* is to tune C so that the tracker and homogeneous modes cancel near z = 0. That gives ε(0) = 0 and |ε| ≤ 1.0×10⁻³γ on 0 ≤ z ≤ z_GW = 0.009.
  - With only the N_eff bound on the 1 - γ²/9 Friedmann deficit (γ ≤ 0.69), ξ ≥ 2.9×10⁻³ (N4b).
  - 2001.06373 Eq. 65 contains a γ² term with no P in it. If that term holds at P = 0, γ ≤ 1.1×10⁻⁷ again, whatever the branch (not rederived).
  - With the k* condition along the GW path, ξ ≥ 1.5×10¹⁴.
  - In every case the strong-coupling kill holds: ≤ 8 MeV at the PDG mass.
- Λ evolution: d ln Λ/d ln a ~ γε. That is negligible for SN/BAO on the tracker and ≲ 0.3 in the cancellation branch at γ ~ 0.7, which was not fitted. BBN: ε_h is suppressed as a³ relative to the tracker in the radiation era. The tracker deficit γ²/9 is what N_eff bounds.

## 6. Caveats

Power counting with O(1) factors untracked. Leading order in γ (the γ ~ 0.7 cancellation branch has O(γ²) ≈ 50% corrections, still ≪ the 6 orders of shortfall). Tensor-torsion q not tracked in the θ equation. Local nonlinear θ near masses not solved. The sign of the θ kinetic term follows MZ's action sign (magnitude only is used).
