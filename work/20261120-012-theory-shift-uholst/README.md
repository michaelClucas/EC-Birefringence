# 20261120-012 theory: shift symmetry of torsion pseudoscalars; u-Holst with a hidden scalar

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism.

## Question
- **A.** A torsion pseudoscalar couples to fermions as (∂φ/f)J₅. Two sub-questions:
  1. Does this give φGG̃ and a QCD mass? Does 008's QCD-mass kill of EC-03 hold?
  2. Is g_φγ the anomaly value, or is it decoupled?
  
  Added on request: the Green-Schwarz loophole (skeptic 011), 2506.11836 (Karananas-Shaposhnikov-Zell) and 2310.09152 (Nakai et al.).
- **B.** Does the 5D u-Holst term give a parity-odd photon coupling when the rolling field is a hidden scalar rather than the radion? Re-check the 4D Holst/Nieh-Yan "g = 0" claim. Derive the matter couplings and test them against bounds.

## Method
Derivations are in `derivation.md`. Each step has a script with PASS/FAIL checks:

| script | content |
|---|---|
| `check_chiral_rotation.py` | θ̄ bookkeeping in sympy for four cases: derivative only, reading A, GS term, KSVZ control. KR dualisation with and without a GS-modified Bianchi identity. LO 2-flavour chiral Lagrangian V(φ), minimised numerically. |
| `check_holst_bi.py` | 4D, all 24 contortion components eliminated, with a dynamical Barbero-Immirzi field and minimally coupled fermions. Compared with Perez-Rovelli (gr-qc/0505081 Eq. 15) and the Taveras-Yunes (0807.2652) kinetic structure. |
| `check_formfactor.py` | Derivative-coupling φγγ form factor B₁. Closed form checked against an independent Feynman-parameter integral, plus series and limits. Suppression numbers. |
| `check_uholst_chi.py` | 5D: c(χ)·u-Holst with the 50 contortion components eliminated; exact-u case; contortion dipole response; helicity ODE with Z(τ) and Y(τ) both varying; dβ = -c dlnα; winding scalar Higgsing the KK photon. |
| `bounds_uholst.py` | CMB α, MICROSCOPE d_e, clocks, quasar (from a field ODE on ΛCDM), and an EDM estimate. |

Inputs come from `data/`. New file `data/shift_uholst_inputs.csv` holds:
- MICROSCOPE d_e < 4×10⁻⁵ (2212.06032 Sec. III.3);
- eEDM < 4.1×10⁻³⁰ e·cm (2212.11841 Eq. 4);
- PDG m_γ, not opened and not used in a check.

Literature read: 2506.11836 (HTML), gr-qc/0505081, 0807.2652 and 2212.11841 (PDF text), and the abstracts of 1903.12559, 2011.10025, 2212.06032 (HTML) and 2310.09152.

Rerun from the repo root (about 1.5 min, one core):
```
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-012-theory-shift-uholst/run_checks.py
```
Last run (`run.log`): **37 PASS, 0 FAIL**.

## Result

### A. Shift symmetry
1. **No QCD mass for KR, BI×Nieh-Yan or PGT 0⁻.**
   - φ enters only through ∂φ [S2], so the action is exactly invariant under φ → φ + const.
   - A chiral rotation moves φ into θ_GG and into the quark-mass phases with opposite signs, so θ̄ is φ-independent [S1a].
   - V(φ) is flat in the chiral Lagrangian [S3a].
   - 008's reading A keeps the anomaly shift but drops the mass phases, which half-applies a field redefinition. **m_QCD = 0 exactly.** The 008 kill does not hold.
   - This agrees with 2506.11836 (Eqs. 7-11), which shows the same for gluons and calls the 1502.03694 and Mercuri constructions unsuccessful.
2. **Holst-BI with minimal fermions** is new here. The exact torsion elimination gives L_eff = -(3/4κ)(∂γ̄ + J₅/2)²/(1+γ̄²).
   - It reproduces Perez-Rovelli's γ²/(1+γ²) and the Taveras-Yunes kinetic structure.
   - In the canonical field it splits into ∂G(φ)·J₅ plus -(3κ/16)J₅²/cosh²(φ/F), with F = √(3/2)M_P.
   - The J₅² term breaks the shift symmetry, but it is parity-even and even in γ̄, so there is no θ̄ coupling and no QCD-axion mass. It does give a UV-sensitive non-topological potential, not computed here.
3. **Photon (question 2).** For a derivative coupling, g_eff = g_anom·Σ N_cQ²B₁(τ_f)/8, with B₁ → -m_φ²/12m_f².
   - g_eff/g_anom = 1.4×10⁻⁷⁸ at m_φ = 6×10⁻³³ eV, 4×10⁻⁷⁰ at 10⁻²⁸ eV, and 4×10⁻⁷¹ for q ~ H(z*).
   - EC-03''s g_T = 4.67×10⁻²¹ GeV⁻¹ becomes ≲ 10⁻⁹⁰ GeV⁻¹.
   - Plasma terms are bounded in scale by n_e/m_e³ ~ 10⁻²⁹ (estimate, not computed).
   - 2506.11836 does not treat photons. The decoupling result is what adds to it here: shift symmetry alone forbids a potential, not a ∂φ-coupling.
4. **Green-Schwarz loophole.** With H = dB + c_GS(ω_CS(A) + ω_CS(G)), dualisation gives an explicit -c_GS φ(FF̃ + GG̃) [S2c]. This restores the photon coupling, and also θ̄ = θ̄₀ + λφ/f and the textbook QCD mass [S1c, S3c]. With c_GS at the anomaly value this is exactly 008's reading A: m = 1.4-8.6×10⁻¹² eV, full washout, dead. **The GS/CS term is an addition beyond EC.** EC axial torsion is algebraic, and the SM axial anomaly is global, so nothing requires the term. With GS for hypercharge only, the model is a 4D KK-03-type ALP with a free coupling.
5. **2310.09152.** After torsion is integrated out, minimal EC + SM gives SM plus Planck-suppressed J² and ∂φ·J (LEFT-type) operators. EC-00, EC-02 and EC-05 are therefore applications of their no-go; gravity is omitted there, and the torsion elimination done here supplies it. The extra light pseudoscalar is outside LEFT and is covered by A3-A4 here.

### B. u-Holst with a hidden scalar
- **4D**: with F = dA, torsion never sees the photon. φ×Holst or φ×Nieh-Yan without fermions gives g = 0 exactly; the iter 1/2 claim holds. With fermions the case reduces to A.
- **5D, c → c(χ)** [U1]: L_eff = -(3c²/2)W² - 4c³⋆(u∧W∧W) - 3c′²(∂χ)². There is no cross term, and the χ kinetic correction is healthy. So -¼Z₀(1+12c²)F² - ¼Z₀16c³FF̃ with χ-dependent c is a torsion-generated χFF̃. It is zero without torsion, and zero for exact u [U2].
- For any Z(τ), Y(τ): β = ½|∫Y′/Z dτ| [U4]. Then **dβ = -c dlnα exactly** [U5]: the odd term always comes with α variation.
- χ is not gravitationally coupled. Its matter coupling is α(χ), plus a fermion dipole if fermions live in the bulk. Bounds:
  - CMB α: c ≥ 0.34.
  - MICROSCOPE: c ≥ 340 at Δχ = 0.39 M_P (133 at Δχ = M_P).
  - Clocks: c ≳ 5×10⁴, nearly m-independent for 10⁻³³ to 10⁻²⁹ eV, because χ̇₀ ≈ H₀Δχ in the matter era.
  - EDM estimate: c ≲ 140.
- **Further costs:**
  - The orientation must be non-dynamical. A winding scalar that orients u gives the KK photon a mass [U6].
  - SM fermions are KK-neutral (004), so the KK photon is not the CMB photon without an extra mechanism.

## Verdict
- **EC-03 (torsion axion with universal torsion-fermion coupling):** still dead, by a different route. The QCD-mass kill of 008 is wrong for derivative couplings (m = 0). The cause of death is g_eff ≤ 4×10⁻⁷⁰ g_anom, so β ≈ 0. The GS variant restores 008's mass and kill, but the GS term is not EC.
- **EC-03':** there is no torsion-generated photon coupling. The torsion-scale g = 8α/(πf) exists only as a non-derivative φFF̃ put in by hand (or by a GS term with no gluon part). EC-03' should be reclassified as a **generic ALP with f set at M_P by hand; torsion plays no role.**
- **Holst-BI (variant a) with fermions:** no θ̄ coupling and no photon coupling. It has a UV-sensitive J₅²-induced potential. Not a β source.
- **u-Holst × hidden scalar (suggested KK-05):** this is a candidate mechanism in which torsion does the work and the rolling field is not gravitationally coupled. It is **excluded at estimate level, otherwise extremely disfavoured.** It needs:
  - an explicit P'-breaking background;
  - the CMB photon to be the KK photon;
  - c ≳ 5×10⁴ (clocks, from the exact dβ = -c dlnα).
  
  The EDM estimate (c ≲ 140) would close the window.

## Open issues
1. The Holst-BI J₅² potential size (UV-sensitive ⟨J₅²⟩, cf. 2506.11836 Sec. "Quantum corrections") is not computed. It cannot give a θ̄ coupling at tree level. Loop-generated P-odd terms need a γ̄-odd interaction, and none exists with minimal coupling. Non-minimal couplings (ζ̃ of 2506.11836) are not analysed.
2. Plasma (finite-density) corrections to the derivative-coupling φγγ vertex are bounded only by a dimensional estimate (n_e/m_e³ ~ 10⁻²⁹).
3. EDM estimate:
   - Pauli normalisation from 004 with O(1) factors unchecked.
   - EDM-vs-MDM assignment by c-parity (P'-odd ⇒ EDM), not by explicit 4D spinor reduction.
   - Meaningful only if the electron couples to the KK photon.
4. The sign of dβ = -c dlnα relative to the β convention is unchecked; magnitudes only.
5. The clock bound assumes a free quadratic potential and a field that is still moving today. A field that stops rolling after z ~ few would need a non-standard potential.
6. The KK charge problem (SM fermions KK-neutral) is not addressed. It affects every KK-photon mechanism, including 009's KK-04.
7. 1903.12559 and 2011.10025 were read at abstract level. The B₁ formula follows 1708.00443 as used in 008 and was independently re-derived here [F1, F2].
8. The PDG photon-mass row in `data/shift_uholst_inputs.csv` was not opened. It is not used numerically.
9. Errors found in other folders (not edited):
   - 008 Sec. 1 table and Verdict: the QCD mass for KR/NY-BI is an artefact of reading A.
   - STATE.md/CLAIMS.md EC-03 "QCD mass 1.4-8.6×10⁻¹² eV" should become "m = 0 (shift symmetry); dead by photon decoupling", and EC-03' "g torsion-fixed" should become "g by hand".
   - 009 derivation Sec. 4: the general formula is β = ½∫Y′/Z. Their constant-c_odd version is a special case and is correct as stated.

## Files
`common.py`, `check_chiral_rotation.py`, `check_holst_bi.py`, `check_formfactor.py`, `check_uholst_chi.py`, `bounds_uholst.py`, `run_checks.py`, `run.log`, `derivation.md`. New data: `data/shift_uholst_inputs.csv`.
