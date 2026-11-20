# Derivations: shift symmetry of torsion pseudoscalars, and the u-Holst term with a hidden scalar (task 20261120-012)

Conventions: `CONVENTIONS.md`. Check labels [S*], [H*], [F*], [U*], [B*] refer to `run.log`. Measured inputs come from `data/` (new file `data/shift_uholst_inputs.csv`).

## A. Shift symmetry

### A1. The structure of every torsion-pseudoscalar coupling to fermions

All variants reach the fermions through the minimal Dirac coupling to the totally antisymmetric contortion, i.e. S_μ J₅^μ. After the torsion is eliminated, the pseudoscalar φ enters as follows.

- **KR torsion axion** (1502.03694). T = H = dB, dualised with a Lagrange multiplier for dH = 0. Result [S2]: L = -½(∂φ)² + a ∂φ·J₅ - (a²/2)J₅². φ appears only through ∂φ.
- **BI field × Nieh-Yan** (0902.0957, 0911.2698). β·d(e^a∧T_a) = -dβ∧e^a∧T_a + total derivative, and e^a∧T_a is the axial torsion 3-form. So β also appears only through ∂β.
- **PGT 0⁻** (task 007): S = ∂φ/√A. Derivative only.
- **BI field in the Holst term, minimal fermions.** Here 1/γ(x) multiplies e∧e∧R(ω), which is not a derivative coupling. Eliminating the 24 contortion components with γ̄ = 1/γ(x) dynamical gives [H2-H4c]

  L_eff = -(3/4κ) (∂γ̄ + J₅/2)² / (1 + γ̄²).

  This is a perfect square. It reproduces Perez-Rovelli (gr-qc/0505081 Eq. 15) for the J₅² coefficient, -(3κ/16)γ²/(1+γ²) [H2], and the Taveras-Yunes kinetic structure 1/(1+γ̄²) (0807.2652) [H4]. Two fields can be defined:
  - Canonical field: γ̄ = sinh(φ/F) with F = √(3/2) M_P.
  - G(φ), defined by ∂G = ∂γ̄/(1+γ̄²).

  In terms of these, L = -½(∂φ)² - (3/4κ)∂G·J₅ - (3κ/16) J₅²/cosh²(φ/F). The ∂G·J₅ term is derivative-type. The last term breaks the shift symmetry [H5], but it is even in γ̄ [H6] and parity-even, so it cannot shift θ̄. This is the "non-derivative J² coupling" route of Karananas-Shaposhnikov-Zell (2506.11836, App. B, Sec. on quantum corrections).

### A2. No QCD potential from ∂φ·J₅ (question 1)

A field-dependent chiral rotation q → exp(i c_q φγ₅/f) q removes the derivative coupling. It produces φ in two places:
- the colour anomaly: θ → θ - 2Σc_qφ/f;
- the quark-mass phases: ρ_q → ρ_q + 2c_qφ/f.

So θ̄ = θ + arg det M is φ-independent [S1a]. The argument does not depend on the rotation: the action depends on φ only through ∂φ, so it is invariant under φ → φ + const. This holds exactly, perturbatively and non-perturbatively, because the variation ∫∂(const)·J₅ vanishes identically whatever the anomaly. The anomaly appears only when the operator is rewritten, never in the variation.

In the LO two-flavour chiral Lagrangian, V(φ) is flat to machine precision [S3a]. Dropping the mass phases, as in the reading-A bookkeeping of 008 (1502.03694, 0911.2698), instead gives V = -K√(m_u² + m_d² + 2m_u m_d cos θ̄(φ)) [S3b], with the textbook curvature. That reading double-counts: it keeps half of an exact field redefinition.

**The QCD-mass kill of 008 does not hold for the KR, NY-BI and PGT-0⁻ variants: m_QCD = 0 exactly.**

This agrees with Karananas-Shaposhnikov-Zell 2506.11836 (their Eqs. 7-11, "Proof to the contrary"), who call the earlier constructions [refs 33-39 there, including 1502.03694, Mercuri, Lattanzi-Mercuri] unsuccessful for this reason. They treat gluons only.

For the Holst-BI variant (A1), the shift symmetry is broken, but only by the parity-even (κ/cosh²)J₅² term. This produces no θ̄ coupling at tree level. It does produce a potential V ~ (3κ/16)⟨J₅²⟩/cosh²(φ/F). That potential is UV-sensitive: the perturbative ⟨J₅²⟩ is power divergent. It is not a QCD-axion mass. Its natural size is not computed here (open issue).

### A3. The Green-Schwarz loophole (skeptic 011 open issue)

If the KR field strength is H = dB + c_GS(ω_CS(A) + ω_CS(G)), then dH = c_GS(F∧F + tr G∧G). The dual Lagrange multiplier now multiplies dH - c_GS X. After dualisation, L contains -c_GS φX, a non-derivative coupling [S2, S2c]. Under φ → φ + const the action shifts by c_GS·const·∫X, and ∫tr G∧G is an integer instanton number. So:
- the gluon part breaks the shift symmetry to a discrete one;
- θ̄ = θ̄₀ + λφ/f [S1c];
- the QCD potential and mass return, with curvature K m_u m_d/(m_u+m_d)λ²/f² [S3c].

If c_GS is set to the anomaly coefficient, this reproduces the 1502.03694 / 0911.2698 numbers exactly: m = 1.4-8.6×10⁻¹² eV, full washout (008). EC-03 with GS is therefore dead by 008's original route.

**The GS term is an addition beyond Einstein-Cartan.** In EC the axial torsion is algebraic and has no Bianchi identity of its own. Even in the KR identification T = dB, nothing in EC supplies the Chern-Simons pieces. In 4D the SM axial current anomaly is a global anomaly, not a gauge anomaly, so consistency does not require GS. With a GS term for hypercharge only (no gluons), the model is the 4D analogue of KK-03: an ALP whose g = c_GS-dependent coupling is a free parameter, with torsion fixing only the kinetic normalisation.

### A4. Photon coupling (question 2)

With only ∂φ·J₅, the φγγ amplitude is (α/πf)Σ_f N_cQ_f²c_f B₁(τ_f), with B₁ = 1 - τf(τ)² and τ = 4m_f²/q² (1708.00443; Quevillon-Smith 1903.12559, "axions are blind to anomalies"). The two terms have separate origins:
- the "1" is the anomaly;
- τf² is the 2imψ̄γ₅ψ (pseudoscalar) triangle. Its closed form was checked against an independent Feynman-parameter integral [F1].

Limits:
- Heavy fermions: B₁ = -q²/12m² - q⁴/90m⁴ [F2].
- q ≫ m: B₁ → 1 [F3].
- On shell, q² → m_φ², because □φ = m_φ²φ.
- The electron dominates (Σ N_cQ² = 8 [F4]).

This gives g_eff/g_anom ≈ (m_φ/m_e)²/96:
- 1.4×10⁻⁷⁸ at m_φ = 6×10⁻³³ eV;
- 4×10⁻⁷⁰ at 10⁻²⁸ eV;
- 4×10⁻⁷¹ for q ~ H(z*) [F5].

So the EC-03' torsion value g_T = 8α/(πf) = 4.67×10⁻²¹ GeV⁻¹ becomes g_eff ≲ 10⁻⁹⁰ GeV⁻¹.

Plasma terms are not computed. Their natural scale, n_e/m_e³ ≈ 1.4×10⁻²⁹ at z*, is itself fatal [F6]. The shift symmetry also forces any medium term to be a ∂φ-type CFJ term, which is the EC-05 class of task 007.

**2506.11836 does not treat photons; the extension is ours.** The same exact-shift-symmetry argument kills a potential, but not a coupling. A shift-symmetric theory still allows ∂φ-couplings, so for photons one needs the extra decoupling statement, B₁ → 0 for vector-like fermions heavier than q. That is the B₁ computation here. Chiral gauge bosons are different: W, Z do not decouple (Quevillon-Smith; Bonnefoy et al. 2011.10025). This is irrelevant for CMB photons below m_e.

Nakai et al. 2310.09152 show that no SMEFT/LEFT operator on a cosmological background gives isotropic β. They omit gravity. Our 4D EC results are partly an application of it: after torsion is integrated out, minimal EC + SM is SM + Planck-suppressed J² and ∂φ·J operators, i.e. LEFT operators. EC-00, EC-02, EC-05 and the derivative-coupled pseudoscalars then fall under their statement. The additions here are:
- the explicit torsion elimination, which shows that the EC operators are of LEFT type;
- the extra light field φ, which is outside LEFT. Its photon coupling is computed here, not covered by them.

## B. u-Holst term with a non-radion scalar

### B1. 4D (re-check of iter 1/2)

In 4D EC + Maxwell with F = dA, the photon does not touch the connection. With a pseudoscalar coupled to the Holst or Nieh-Yan density and no fermions, torsion is sourced by ∂φ alone, and L_eff contains no A. So g = 0 exactly. The claim holds.

With fermions, the result is the A1/A4 structure: ∂φ·J₅, plus the Holst-only J₅² term. g_eff is suppressed as above.

### B2. 5D, photon = KK photon, c → c(χ)

Take L ⊃ c(χ) u∧e^A∧e^B∧R_AB(ω) with u = e⁵ oriented (the P'-odd background of 009) and χ a P'-even bulk scalar zero mode. After parts, the linear source is X = d(c u) = cW + c′dχ∧u, and the odd quadratic coefficient is c pointwise. Eliminating the 50 contortion components [U1] gives

  L_eff = -(3c²/2) W·W - 4c³ ⋆(u∧W∧W) - 3c′² (∂χ)²,

with no W×dχ cross term. The (∂χ)² piece has a healthy sign [U1b]. Adding Maxwell from R̊ gives L = -¼Z₀(1+12c(χ)²)F² - ¼Z₀·16c(χ)³FF̃. The radion is assumed stabilised.

If instead u = dχ (exact), then X = 0 and K = 0, so there is no photon term [U2]. This is consistent with the P' theorem: a P'-even field alone cannot give FF̃.

For L = -¼ZF² - ¼YFF̃ with both Z and Y time-dependent, the helicity phase difference is ∫Y′/Z dτ [U4, numerical]. So β = ½|∫Y′/Z dτ|. The earlier formula of 009, -(c_odd/2)Δln Z, is the special case Y = c_odd Z with constant c_odd.

With Y = 16c³ and Z = 1+12c² [U5]:

  **dβ = -c d ln α** (exact),  β(c) = 2c - arctan(2√3 c)/√3 + const.

The parity-odd photon term generated by torsion therefore always comes with an α variation of fixed ratio, set by c.

### B3. Matter couplings of χ and bounds

χ has no direct matter coupling. It couples through α(χ) (a dilaton-type d_e coupling) and, if fermions live in the bulk, through the torsion-induced dipole term. In the large-c regime Δc = β/2 ≪ c:
- **CMB α** (1406.7482): c ≥ 0.34 [B1].
- **MICROSCOPE** (|d_e| < 4×10⁻⁵, 2212.06032 Sec. III.3): c ≥ √2βM_P/(d_eΔχ). That is 133, 340 or 1330 for Δχ = 1, 0.39 or 0.1 M_P [B2]. Δχ = 0.39 M_P is the 008 energy-cap excursion.
- **Clocks** (2010.06620): in the matter era χ̇₀ ≈ H₀Δχ for every m that rolls after recombination, because the oscillation amplitude dilutes as a^{-3/2}. This gives c ≳ 5-8×10⁴ across 10⁻³³ to 10⁻²⁹ eV [B3]. The quasar bound (2112.05819) is weaker: c ≳ 10³ for thawing fields, down to ≲ 1 for m ≥ 10⁻²⁹ eV.
- **Electron EDM** (order of magnitude). The contortion response is K_[ab5] = c W̃_ab + (4/3)c² W_ab, exact [U3]. The c-odd piece is P'-odd, so it is the EDM-type one. With the Pauli normalisation of task 004, (3/4)K_[ab5]χ̄Γ^{ab5}χ with k = √(2κ₄), this gives d ~ 2.8×10⁻³² c e·cm. JILA (2212.11841, 4.1×10⁻³⁰ e·cm) then gives c ≲ 140. This applies only if the electron couples to the KK photon [B4].
- **Orientation.** A dynamical field that orients u, for example a compact bulk scalar with winding w, gives the KK photon a Stückelberg mass: (∂π - (w/R)A)² [U6]. So the orientation must be a non-dynamical explicit P'-breaking background, as in 2608.09299's kinetic-less Stückelberg fields.
- **Which photon.** SM fermions are KK-neutral zero modes (task 004). The KK photon is therefore not the photon that Thomson-scatters the CMB, unless the KK charge problem is solved by an unstated mechanism.

**Verdict B.** A torsion-essential, non-gravitationally-coupled version exists: c(χ)·u-Holst, candidate mechanism KK-05 (suggested). It costs:
- an explicit P'-breaking background;
- c ≳ 5×10⁴ (clocks; α̇ ≈ (β/c)H₀ today is unavoidable);
- the program's assumption that the CMB photon is the KK photon.

The EDM estimate (c ≲ 10²) would close the window [B5]. Without it, the mechanism survives only at c ≳ 5×10⁴, where the photon kinetic term is 3×10¹⁰ times its Maxwell value.
