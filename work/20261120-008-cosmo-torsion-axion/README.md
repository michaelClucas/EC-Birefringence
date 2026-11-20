# 20261120-008 cosmo-torsion-axion (EC-03, iteration 2)

The ACT DR6 β = 0.215° ± 0.074° (arXiv:2509.13654) is 2.9σ, has systematics caveats, and may not be cosmological. Everything here is a candidate mechanism.

## Question

Candidate mechanism EC-03: a torsion-associated pseudoscalar φ with L ⊃ -(g/4)φFF̃, so β = -(g/2)(φ_0 - φ̄_LSS) (CONVENTIONS.md). Variants: (a) Barbero-Immirzi (BI) field from the Holst term, (b) BI/pseudoscalar coupled to the Nieh-Yan term, (c) Kalb-Ramond (KR) torsion axion.

1. Per variant: g_φγ (origin, size, fixed or free), f, forced gluon coupling, QCD mass.
2. β/(gφ_i) vs m from the field ODE on Planck 2018 ΛCDM, with recombination washout.
3. QCD-mass variants: can β reach 0.2° with bounded g and bounded field excursion?
4. A variant without QCD mass: which (m, g) region gives β, and what the torsion origin adds over a generic ALP.
5. Anisotropic birefringence A_CB from inflationary δφ = H_inf/2π vs SPT-3G.

## Method

- **Literature.** Full text read (pdftotext) for 1502.03694 (Castillo-Felisola et al., KR/BI torsion axion), 0911.2698 (Lattanzi-Mercuri, BI = QCD axion), 0902.0957 (Calcagni-Mercuri, BI + Nieh-Yan), 0807.2652 (Taveras-Yunes, BI in Holst), 1511.02867 (Grilli di Cortona et al., QCD axion mass and photon coupling). Abstract only: 1903.02666 (Fedderke-Graham-Rajendran washout), 1903.12559, 1708.00443.
- **Inputs.** New file `data/cosmo_torsion_axion_inputs.csv`: Planck 2018 ΛCDM (1807.06209 Table 2), T_CMB, recombination constants (RECFAST values), QCD axion m_a f_a and c_γ (1511.02867 Eq. 22, Eq. 2/17), torsion-axion scales (1502.03694 Eq. 24, 0911.2698), ρ_loc and the Rayleigh factor (2512.21888), BK18 r. Bounds from `data/bounds.csv` (A-Chandra, U-frac-Hlozek, U-frac-Rogers23, U-washout-J0); β and A_CB from `data/birefringence_measurements.csv`.
- **Background and recombination** (`common.py`, `recombination.py`). Flat ΛCDM with photons and massless neutrinos. Hydrogen-only Peebles three-level atom with the RECFAST fudge factor 1.14. No helium, no reionization. Visibility g(z) normalized on 200 < z < 2000.
- **Field** (`field_evolution.py`, `beta_scan.py`). φ'' + 2(a'/a)φ' + a²m²φ = 0 in conformal time (DOP853, rtol 1e-10), spectator field, starting frozen at a = 1e-9. Once m/H = 300, the solution switches to WKB, φ = A(a_s/a)^{3/2}cos(m(t - t_s) + δ), with phase-accurate cosmic time. φ̄_LSS = ∫g(z)φ(z)dz. The washout factor is |φ̄_LSS|/envelope(z_peak).
- **Variants and regions** (`variants.py`). QCD numbers for three f_a choices. Generic-ALP minimum coupling g_min(m) = 2β/(F φ_i,max), with F = |φ_0 - φ̄_LSS|/φ_i (best phase today) and φ_i,max set by ρ_φ(today) ≤ cap × ρ_DM. The 2512.21888 local-DM regime is reproduced at m = 1e-24 eV.
- **Anisotropy** (`anisotropy.py`). Superhorizon δφ follows the background, so δφ_LSS/δφ_i = φ̄_LSS/φ_i. The massless transfer function T(k) is integrated to τ*. Projection: S_L = 2L(L+1)∫dk/k T² j_L²(kD).

Rerun (about 8 min total, one core, each script < 2.5 min):
```
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-008-cosmo-torsion-axion/run_checks.py > work/20261120-008-cosmo-torsion-axion/run.log
```
Last run (`run.log`): 25 PASS, 0 FAIL.

## Result

### 1. Variants: coupling, f, gluon coupling, QCD mass

| variant | photon coupling | f | gluon GG̃? | m |
|---|---|---|---|---|
| (a) BI field in Holst (0807.2652) | none. Pure gravity: massless scalar with non-canonical kinetic term, no potential, no matter coupling. With fermions it reduces to (b)/(c) | not defined in 0807.2652 | only through fermions, as (c) | 0 without fermions |
| (b) BI × Nieh-Yan (0902.0957; 0911.2698) | 0902.0957: canonical massless pseudoscalar, no gauge coupling. 0911.2698: anomaly term, NY "charges" = 1 for every fermion (their footnote 27) | f = 1.98×10¹⁸ GeV, fixed (text after their Eq. 4) | yes, universal | 0911.2698 Eq. 5: 3.04×10⁻¹² eV; 1511.02867 gives 2.88×10⁻¹² eV |
| (c) KR torsion axion (1502.03694) | anomaly via ∂φ·J₅ → φ∂·J₅ with the massless-fermion ABJ equation (their Eqs. 17, 22, 23): coefficient αQ̄²/(πf), Q̄² = 8 | f = √(8/3)M_P = 3.98×10¹⁸ GeV, fixed (Eq. 24); θ_eff = θ + N_qφ/f, so f_a = f/N_q = 6.6×10¹⁷ GeV | yes. Torsion couples to every quark's axial current | 1.43×10⁻¹² eV (f_a = f, as in their Eq. 27) to 8.6×10⁻¹² eV (f_a = f/6) |

Both 1502.03694 and 0911.2698 compute the photon and gluon couplings at one loop, from the anomaly, with no free parameter. Universal unit charges give E/N = 8/3 (checked: = 2Q̄²/N_q). With QCD mixing, g_φγ = (α/2πf_a)|8/3 - 1.92| = 0.75 α/(2πf_a), which is **2.2×10⁻²² to 1.3×10⁻²¹ GeV⁻¹** across the three f_a choices. The papers differ by up to a factor 6 in f_a (KR √(8/3)M_P vs BI 1.98×10¹⁸ GeV, N_q in or out), but the β conclusion below does not depend on f_a.

**Caveat on the coupling itself (unresolved, for the theorist and skeptic).** In all three variants the field enters only through ∂_μφ J₅^μ, which is shift symmetric. 1502.03694 trades this for φ × (anomaly) using ∂·J₅ with massless fermions. With massive quarks and leptons, ∂·J₅ also contains 2im ψ̄γ₅ψ. A field-dependent chiral rotation then moves φ into arg det M as well as into θ_GG. θ̄ = θ + arg det M is φ-independent, so no QCD potential arises. The photon amplitude from a derivative coupling to a heavy fermion goes as B₁ → -ω²/(12m_f²) (heavy-fermion limit checked symbolically; loop function from 1708.00443, not opened here). For the cosmological field, ω ~ H_LSS, which gives a suppression of about 3.5×10⁻⁷⁰. Two readings follow:
- **(A) anomaly terms physical** (papers as written): QCD mass plus g as in the table.
- **(B) shift symmetry exact**: no QCD mass, and the photon coupling is suppressed by about 10⁻⁷⁰, so β ≈ 0.

Both readings fail (Sec. 3). This was open issue 2 of task 003.

### 2. β/(gφ_i) vs m (`beta_vs_m.png`, `scan.csv`)

- m ≲ 3×10⁻³⁴ eV: frozen field, slow roll, β/(gφ_i) ≈ μ²[I(τ₀) - I(τ*)]/2. This is the analytic slow-roll limit (check: ODE vs analytic agree to 2e-6).
- 6×10⁻³³ ≲ m ≲ 3×10⁻²⁹ eV: plateau β/(gφ_i) ≈ +0.50. The field is frozen through recombination and has decayed by today.
- Washout: φ̄_LSS/φ_i = 0.99 (10⁻²⁹), 0.92 (3×10⁻²⁹), 0.43 (10⁻²⁸), -0.009 (10⁻²⁷), -3×10⁻⁶ (3×10⁻²⁷). Below 10⁻⁶ for m ≳ 5×10⁻²⁷ eV, which is the numerical noise floor (the true value is smaller). Half washout at m ≈ 1.1×10⁻²⁸ eV.
- For m ≳ 3×10⁻²⁸ eV only the observer term -(g/2)φ_0 survives. It is set by the local DM field (2512.21888 regime). Reproduced here: φ₀ = 2.48×10¹² GeV, g_req = 2.18×10⁻¹⁵, g_washout = 9.60×10⁻¹⁶ GeV⁻¹ at 10⁻²⁴ eV. Their Eqs. 4, 6, 12 match to 1%, and the ratio 2.27 rules out a single field there.

### 3. QCD-mass variants: β with bounded g and excursion

With m ≈ 1.4 to 8.6×10⁻¹² eV ≫ H(z*) ≈ 3×10⁻²⁹ eV, the field oscillates long before recombination. Across the visibility function there are about mσ_t ~ 10¹⁵ oscillations, so φ̄_LSS = 0 at any precision (washout complete). The only remaining β is the local observer term. Its amplitude, assuming all local DM (0.4 GeV/cm³) is this field, is

β_amp = α|E/N - 1.92|√(2ρ_loc)/(4π m_a f_a) = **1.9×10⁻²² rad = 1.1×10⁻²⁰ deg**

This is independent of f_a, because both g and m scale as 1/f_a (sympy PASS). It oscillates with a 0.5 to 3 ms period, so it averages to zero in any CMB map. **Shortfall: 10^19.3 even for the instantaneous amplitude. The static isotropic β is exactly zero.** The Chandra bound (8×10⁻¹³ GeV⁻¹) does not bind: g is 9 orders below it. The excursion bound φ_i ≤ πf_a is weaker than the DM bound (0911.2698 needs θ_i ≤ 1.2×10⁻⁴ for Ω ≤ Ω_DM). For comparison, a field frozen at πf_a would give β = α|E/N - 1.92|/4 = 0.078° (0.28° without QCD mixing). The QCD mass is what kills it, not the coupling size.

Reading (B) of the caveat gives g_eff ~ 10⁻⁷⁰ × the anomaly value, so β ≈ 0 there too. **Verdict: EC-03 (a), (b), (c) as published (universal torsion-fermion coupling) is dead as a source of β. This is a candidate no-go, pending the skeptic.**

### 4. Without QCD mass (by assumption) and the generic-ALP baseline (`g_required.png`)

**Generic ALP** (g, m, φ_i free). Energy cap: ρ_φ,0 ≤ 5% ρ_DM (Hlozek, 10⁻³² to 10^-25.5 eV), ≤ 1% for 10⁻³⁰ to 10⁻²⁸ (Rogers), and an assumed 5% for m ≤ 10⁻³² eV.
- Minimum coupling: g_min = **7.4×10⁻²¹ GeV⁻¹ at m ≈ 6×10⁻³³ eV** (φ_i = 9.4×10¹⁷ GeV = 0.39 M_P).
- g_min ≈ 9×10⁻²¹ across 10⁻³² to 10⁻³⁰ eV, 2×10⁻²⁰ across 10⁻³⁰ to 6×10⁻²⁹, then rising steeply above 10⁻²⁸ eV.
- Allowed region: **10⁻³⁵ ≲ m ≲ 10⁻²⁷ eV and g_min(m) ≤ g ≤ 8×10⁻¹³ GeV⁻¹** (Chandra), a band about 8 orders wide. Above that, the single-field 2512.21888 result applies.
- In (g, f) terms, φ_i ≤ πf requires g f ≥ 2β/(πF) = 2.2×10⁻³. For an anomaly coupling g = (α/2πf)·c, that means c ≳ 1.9 at the plateau.

**Torsion-fixed coupling without a QCD mass.** This assumes reading (A) for photons only: no gluon coupling, which universal torsion couplings do not allow, so it needs a non-universal or cancelling quark sector. Then g_T = 8α/(πf_CF) = 4.67×10⁻²¹ GeV⁻¹, fixed.
- With the energy cap, the best case is **β_max = 0.136° at m ≈ 6×10⁻³³ eV**, below the ACT 1σ lower edge of 0.141°.
- The excursion cap π f does not bind; the energy cap does. Since β ∝ √(cap), β = 0.215° needs ρ_φ,0 ≈ 12% of ρ_DM (about 3% of ρ_crit) at m ~ H₀. That is the DE-like regime, where no ULA bound exists in `data/bounds.csv`.
- So this case is **marginal**: about 1.0σ short at a 5% cap, and it fits only with a 12% early or late dark-energy-like component.

**What the torsion origin adds over a generic ALP.**
1. It fixes f ~ M_P and therefore g ~ α/M_P. That sits at the bottom edge of the allowed band, which is non-trivial: g_T is within a factor 1.6 of g_min.
2. It forces a universal coupling to quark axial currents. That brings in the QCD mass, which kills β (reading A), or kills both mass and coupling (reading B).

Torsion does not explain why the field would be light (m ~ H₀) and still photon-coupled. That requires removing the gluon coupling by hand, at which point the model is an ALP with f set at M_P. The torsion origin adds a coupling-size prediction and a fatal QCD constraint. It adds nothing else.

### 5. Anisotropic birefringence

δβ(n̂) = (g/2)δφ_LSS(n̂) with P_δφ = (H_inf/2π)². Using g = 2β/(Fφ_i):

A_CB = (βR/F)² (H_inf/(2πφ_i))² S_L,   with R = φ̄_LSS/φ_i.

S_L = 1.00 (L=2), 0.96 (10), 0.78 (30), 0.55 (50), 0.13 (100), 0.015 (200). Modes that entered the horizon before τ* have decayed. On the plateau (R = F = 1), A_CB = 3.6×10⁻⁷ (H_inf/φ_i)². SPT-3G A_CB < 1.2×10⁻⁴ therefore requires only **H_inf/φ_i < 18**. With H_inf ≤ 4.7×10¹³ GeV (BK18 r < 0.036) and φ_i at the energy cap (3 to 9×10¹⁷ GeV), the prediction on the plateau is A_CB ≲ 10⁻¹⁴ (10⁻¹⁵ to 9×10⁻¹⁵). That is about 10 orders below SPT. The worst case in the scan is m = 10⁻³⁵ eV, where F is tiny and g must be large: A_CB ≤ 2×10⁻⁹, still 5 orders below. The isotropic-β candidates in Sec. 4 have no anisotropy constraint and no detectable anisotropy signal. (For the QCD variants, the relevant isocurvature constraint is on DM, H_I ≲ 10¹⁰ GeV per 0911.2698, and has nothing to do with β.)

## Verdict

- **EC-03 (a)/(b)/(c) as published, with universal torsion-quark coupling: dead (candidate no-go, pending skeptic).** Reading A: QCD mass 1.4 to 8.6×10⁻¹² eV, full washout, local-oscillation amplitude 1.1×10⁻²⁰ deg (f_a independent), 19 orders short, and the static β is 0. Reading B (shift symmetry exact): no mass but g_eff suppressed by about 10⁻⁷⁰, β ≈ 0.
- **(a) Holst BI without fermions and (b) Calcagni-Mercuri NY without fermions:** g = 0. A photon coupling has to be added by hand, which makes them a generic ALP.
- **Torsion-scale coupling g = 8α/(πf), f = √(8/3)M_P, with the gluon coupling removed by assumption:** marginal candidate mechanism. β_max = 0.136° at m ≈ 6×10⁻³³ eV with ρ_φ ≤ 5% ρ_DM; 0.215° needs about 12%. Requires an extension that removes the QCD coupling, and it is in reading-B tension (electron loops decouple), so the photon coupling would itself need to be a non-derivative φFF̃ term.
- **Generic ALP baseline:** β = 0.215° allowed for 10⁻³⁵ ≲ m ≲ 10⁻²⁷ eV, 7×10⁻²¹ ≲ g ≲ 8×10⁻¹³ GeV⁻¹. Anisotropy is unconstrained by SPT (H_inf/φ_i < 18).

## Open issues

1. Reading (A) vs (B): does a shift-symmetric ∂φ·J₅ torsion coupling generate φGG̃ and φFF̃ at low energy with massive fermions? The argument here (θ̄ invariance, B₁ decoupling) says no. If the NY divergence absorbed by the BI field (Mercuri) breaks the shift symmetry non-perturbatively, that changes it. Theorist and skeptic task.
2. Energy cap for m ≤ 10⁻³² eV is an assumption (5% of ρ_DM). The torsion-fixed case is decided by this number. A quintessence-type bound on a thawing m ~ H₀ field (w(z) constraints, DESI) is needed. Bounds analyst.
3. Recombination is hydrogen-only Peebles (κ = 1 at z = 1090.8 vs Planck 1089.92). Reionization is ignored: its low-ℓ EB bump sees φ(z ~ 8), which differs from φ̄_LSS for m ≳ 10⁻³² eV. This matters for the β(ℓ) shape (EDE/2601.13624-type analyses), not for the order of magnitude.
4. The spectator approximation: the energy density of φ is not fed back into H. At the 5 to 12% level near m ~ H₀ this shifts F by O(10%).
5. Anisotropy: S_L uses a sharp LSS and a massless transfer function. The observer-term dipole and the SPT L-range weighting were not modeled. All are irrelevant at A_CB ≲ 10⁻¹⁴.
6. f normalization differs between papers (1502.03694 √(8/3)M_P with N_q = 6 vs 0911.2698 1.98×10¹⁸ GeV, N = 1). Only the torsion-fixed g_T depends on it (factor up to 6). The QCD verdict does not.
7. 1903.02666, 1903.12559 and 1708.00443 were read at abstract level only. The washout number is taken via 2512.21888 Eq. 11 (`data/bounds.csv` U-washout-J0).
8. CONVENTIONS.md: no errors found. 2512.21888 writes β = +(g/2)Δφ (Komatsu ε sign). Only magnitudes are compared with it here.

## Files
`common.py`, `recombination.py`, `field_evolution.py`, `beta_scan.py` (→ `scan.csv`, `beta_vs_m.png`), `variants.py` (→ `g_required.png`, `variants_summary.csv`), `anisotropy.py`, `check_ode.py`, `check_numbers.py`, `check_anisotropy.py`, `run_checks.py`, `run.log`. New data: `data/cosmo_torsion_axion_inputs.csv`.
