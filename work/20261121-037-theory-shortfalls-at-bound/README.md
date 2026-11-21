# 20261121-037 theory: shortfalls at the bound (iteration 11)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism. There is no `RESULT.md`.

## Question

Skeptic 035 MAJOR 3: the EC-02, Pauli (any fermionic source) and vector-current shortfalls were quoted at ξ = 1.

1. Is ξ < 9×10¹⁷ the right bound? Find every independent bound on the coefficient of J₅·AF̃ (and of J·AF̃ for the vector route), and convert to ξ (or ξc).
2. Recompute β_max and the shortfall at the strongest allowed ξ / ξc. Does any route reach 0.215° within all bounds?
3. Which bounds depend on the completion?

## Method

- Operator: -½ξS·AF̃ with S = -(κ/2)J gives C J·AF̃ with C = ξκ/4. This is dimension 6 (coefficient Λ⁻²). In Reig's normalisation, 1/Λ*² = 2C. EC sources S universally, so J₅ includes electrons and quarks with the same coefficient as relic neutrinos.
- Gauge structure (`check_gauge.py`, sympy):
  - δO = total derivative − α(∂_νJ_μ)F̃^{μν}.
  - The variation is zero for a closed J and for the homogeneous FRW background. It is nonzero for any inhomogeneous J.
  - So the operator needs a completion. Reig 2611.07578 (HTML opened) argues that a local, Lorentz-invariant completion with a massless photon must contain a light pseudoscalar: a Kalb-Ramond field dual to a massless axion. The other option is a massive (Stückelberg) photon. I use these two classes:
    - **D (dual massless axion).** Axial current: g_aγ g_ae = 4m_eC = ξκm_e. Conserved vector current: the axion coupling to the current is a total derivative, but a contact term 1/(2f²) J·J is left, with f = g_aγ Λ*². I reproduced Reig Eqs. 11, 12 and S14 to within 10%.
    - **S (Stückelberg).** A → A + ∂π/m_γ gives (C/m_γ) π (∂J)F̃, which becomes strongly coupled at Λ_s = (m_γ/C)^{1/3}. Two criteria:
      - S-a: Λ_s ≥ T_BBN (Reig's criterion; T_BBN = 1 MeV is a choice).
      - S-b (NDA): π has to decouple above T = m_μ, otherwise ΔN_eff ≥ 0.39, which exceeds the N_eff 95% room of 0.29. Decoupling uses Γ ~ T⁷/Λ_s⁶ < H.
- Vector route with no completion (`check_fmt.py`): the contact terms that the torsion itself generates, from Freidel-Minic-Takeuchi (hep-th/0507253 Eqs. 18 and 23, opened via ar5iv).
- `shortfalls.py` recomputes β(ξ = 1) for each route. It reproduces 007 (7.5×10⁻³⁰ and 6.3×10⁻²⁹ deg) and 029 (6.91×10⁻²⁷ deg, ξc = 1.7×10³⁶), then scales β linearly in ξ. All inputs are loaded from `data/`. New cited rows are in `data/shortfall_bound_inputs.csv`.

Rerun (about 25 s, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261121-037-theory-shortfalls-at-bound/run_checks.py
```
Last run: `run.log`, 25 PASS, 0 FAIL.

## Result

Shortfall is log₁₀(0.215°/β_max).

| route | bound (class) | max ξ or ξc | source | max β (deg) | shortfall (orders) |
|---|---|---|---|---|---|
| EC-02 ν, cited ξ_ν | CAST g_aγg_ae (D) | ξ < 9.4×10¹⁷ | 2611.07578 + CAST (Reig ref. 45) | 7.0×10⁻¹² | 10.5 |
| EC-02 ν, cited ξ_ν | Chandra g_aγ × RG g_ae (D) | ξ < 1.5×10¹⁵ | 1907.05475 × 2007.03694 | 1.1×10⁻¹⁴ | **13.3** |
| EC-02 ν, cited ξ_ν | Λ_s > T_BBN, PDG m_γ (S-a) | ξ < 2.4×10¹⁹ | PDG 2024 m_γ < 10⁻¹⁸ eV | 1.8×10⁻¹⁰ | 9.1 |
| EC-02 ν, cited ξ_ν | Λ_s > T_BBN, Proca galactic m_γ (S-a) | ξ < 2.4×10¹¹ | hep-ph/0306245 via PDG | 1.8×10⁻¹⁸ | 17.1 |
| EC-02 ν, cited ξ_ν | Λ_s > T_BBN, weakest m_γ (Retino 7.9×10⁻¹⁴ eV) (S-a) | ξ < 1.9×10²⁴ | 1302.6168 | 1.4×10⁻⁵ | 4.2 |
| EC-02 ν, generous N_eff | same rows | | | | 9.6 / 12.4 / 8.2 / 16.2 / 3.3 |
| Pauli, any fermion | CAST product (D) | 9.4×10¹⁷ | as above | 6.5×10⁻⁹ | 7.5 |
| Pauli, any fermion | Chandra × RG (D) | 1.5×10¹⁵ | as above | 1.0×10⁻¹¹ | **10.3** |
| Pauli, any fermion | S-a: PDG / galactic / weakest m_γ | 2.4×10¹⁹ / 2.4×10¹¹ / 1.9×10²⁴ | as above | 1.6×10⁻⁷ / 1.6×10⁻¹⁵ / 1.3×10⁻² | 6.1 / 14.1 / **1.2** |
| vector J_B (c ≈ 2×10¹⁸ at ξ = 9×10¹⁷) | contact 1/(2f²)J_B² < CMS dijet, g_aγ = ALPS II (D) | ξc < 1.5×10²⁵ | 1803.08030, ALPS II via Reig | 1.9×10⁻¹² | **11.1** (12.5 with CAST g, 14.3 with Chandra g) |
| vector J_B | S-a: PDG / weakest m_γ | ξc < 2.4×10¹⁹ / 1.9×10²⁴ | as above | 3.0×10⁻¹⁸ / 2.4×10⁻¹³ | 16.9 / 12.0 |
| vector J_B | FMT torsion contact only (no completion) | c < 3.3×10¹⁵/√(γ²+1) | 1803.08030 | not a bound on ξc | evaded for Immirzi γ < 0.04 (ξ = 9.4×10¹⁷) or γ < 0.22 (ξ = 2.4×10¹⁹) |
| all J₅ routes | S-b NDA decoupling (PDG / weakest m_γ) | ξ < 4.7×10³ / 3.7×10⁸ | NDA, not a published bound | | EC-02 24.8 / 19.9; Pauli 21.8 / 16.9 |

Answers:

1. **ξ < 9×10¹⁷ is not the strongest bound inside its own completion.**
   - In class D the CAST product bound is weaker than combining the independent bounds. The axion is exactly massless, so the Chandra g_aγ < 8×10⁻¹³ GeV⁻¹ applies, and the red-giant g_ae < 1.6×10⁻¹³ applies too. Reig shows that the trapping window and electron g-2 close larger g_ae. Together these give **ξ < 1.5×10¹⁵**.
   - Completion-free bounds on the operator alone (LEP e⁺e⁻→γγ contact scale, lepton g-2 loops, γγ→νν̄ cooling in stars and SN1987A) are all TeV-class by dimensional analysis, i.e. ξ ≳ 10²⁹. They do not reach the needed ξ ≈ 3×10²⁸. I did not compute them precisely. The operator is CP-even, so it gives no EDM.
   - The Dror-Lasenby-Pospelov energy enhancement for non-conserved currents (1705.06726, 1707.01503) is for a new massive vector. Here its analogue is the Stückelberg longitudinal photon, which is class S.
2. **No route reaches 0.215° within the bounds of either class.** For the vector route, the dual-axion contact term kills it independently of the Immirzi parameter, because the signal fixes Λ*² ≈ 7 GeV² for the baryon current. The only reason the skeptic's c ≈ 2×10¹⁸ was not excluded is that no completion had been imposed. The torsion contact by itself has a small-γ loophole (`check_fmt.py`). There is no Result C entry.
3. **Completion dependence:**
   - Every bound that kills these routes depends on the completion: class D uses the axion bounds, class S uses the photon-mass bounds plus a cutoff criterion. The gauge-variation identity is the one completion-independent statement. Reig's claim that D and S exhaust the options was not re-derived here.
   - The weakest corner is class S with the weak end of the weakest photon-mass bound and a 1 MeV cutoff: the Pauli envelope is 1.2 orders short and cited ν sources are 4.2 orders short. That corner ignores two things: the Proca galactic limit, which applies specifically to a Stückelberg mass and adds 8 orders, and π thermalisation (NDA), which adds 15 orders.

## Verdict

- EC-02 stays dead as a torsion mechanism at the largest allowed ξ: **≥ 10.5 orders short in class D (13.3 with Chandra × RG)**, and 4.2-17 orders short in class S depending on the m_γ bound used.
- The Pauli envelope covers any fermionic source: ≥ 7.5 orders short in class D (10.3 with Chandra × RG). In class S it is 1.2-14 orders short.
- The vector-current candidate mechanism is dead in both classes, ≥ 11 orders short.
- Proposed CLAIMS wording: "short by ≥ 10.5 / 7.5 orders at the largest ξ allowed in the dual-axion completion (≥ 13.3 / 10.3 combining Chandra and red-giant bounds); Stückelberg completion 4-17 orders depending on the m_γ bound".

## Open issues

1. Class S is the soft spot. The 1.2-order Pauli margin depends on the T_BBN cutoff choice (shortfall changes by 3·log₁₀(E/1 MeV)) and on which photon-mass bound is used. A real computation of π emission or thermalisation would settle it. The NDA estimate says it is overwhelmingly excluded.
2. Reig's exhaustiveness claim (massless photon ⇒ light pseudoscalar) is taken from the paper, not re-derived. Skeptic 011 calls it completion dependent.
3. The CMS contact normalisation (Eichten-Lane-Peskin, flavour-universal LL) was assumed from the abstract, not checked in the body. Using the LL limit for a VV operator is conservative. This is an O(1) effect, irrelevant to an 11-order margin.
4. ALPS II g_aγ and CAST's g_aγg_ae are as quoted in Reig. The red-giant g_ae was taken from the abstract via a search snippet.
5. The completion-free TeV-class estimates (LEP γγ, g-2, SN1987A) are dimensional analysis, not computed. A careful SN1987A νν̄ emission rate could approach Λ ~ 10 TeV, still below the needed 20-60 TeV.
6. With a universal axial current, the dual axion also couples to nucleons (g_aN = (m_N/m_e) g_ae). SN1987A and neutron-star bounds on g_aN were not added, and could only strengthen class D.
