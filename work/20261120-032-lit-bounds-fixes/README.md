# 20261120-032 literature and bounds fixes (iteration 9)

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here concerns candidate mechanisms. There is no `RESULT.md`.

## Question
1. Are 028's recalled compact-star and corona inputs (`work/20261120-028-skeptic-ec06c-costs/inputs_028.csv`) correct at their sources?
2. Which cited CMB bound on H(z) at high z kills EC-06d's SNY histories (H(z*) = 33-66 × ΛCDM per 025), and by how many orders?
3. What is KK-04's shortfall in orders against Cassini?
4. What is arXiv:2610.02667, and does it change any verdict or novelty statement?

## Method
- Opened at the source (arXiv abstract pages; full text through arXiv/ar5iv HTML or PDF text):
  - 1703.10625 (abs);
  - 1801.00876 (abs: wrong paper);
  - 1806.00061 (HTML, Tables 3-4);
  - 1809.01240 (abs, HTML Sec. 3.1.7);
  - 1611.00629 (abs, PDF text Sec. 3);
  - 1912.05702 (abs);
  - 1112.3850 (abs, PDF text Eq. 5);
  - 1302.6168 (abs);
  - 1807.06209 (PDF text: Sec. 7.5.2, Eq. 67b, Table 2);
  - 1502.01590 (ar5iv HTML: Sec. 5.1.5, Table 3);
  - 2003.07355 (abs, no number used);
  - 2610.02667 (abs, HTML Secs. I-IV).
- New data: `data/compact_star_inputs.csv`, `data/expansion_history_bounds.csv` (both documented in `data/README.md`). Literature: `literature/2610.02667.md`.
- Scripts (they load `data/`, plus 025's `sny_background.csv` and 028's `inputs_028.csv` read-only): `verify_inputs.py`, `hz_bounds.py`, `kk04_shortfall.py`, `ec02_yinkomatsu.py`. Rerun from the repo root (a few seconds, 1 thread):
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-032-lit-bounds-fixes/run_checks.py > work/20261120-032-lit-bounds-fixes/run.log
```
`run.log`: 14 PASS, 0 FAIL.

## Result

**1. Compact-star inputs (`data/compact_star_inputs.csv`).**

| row | 028 value | at source | status |
|---|---|---|---|
| M Sirius B | 1.018 ± 0.011 M_⊙ | same (1703.10625 abstract) | correct |
| R Sirius B | 0.008098 ± 0.000046 R_⊙, "1801.00876" | 0.00802 ± 0.00011 R_⊙ (1806.00061 Table 4) | **citation wrong**: 1801.00876 is a maths paper (random lifts). The value is within 1σ (+1%), but the error bar was not found. Effect on 028 W1: Sirius B edge ×1.02 |
| Teff Sirius B | 25000 K | 25922 ± 296 K (1806.00061 Table 3) | rounded, 4% low; it enters only the scale height |
| typical WD mass | 0.6 M_⊙ | "masses around 0.6 M_⊙" (1611.00629 Sec. 3) | correct (representative) |
| R(0.6 M_⊙) | 0.0125 R_⊙ | not printed in 1611.00629; the Nauenberg zero-T relation gives 0.0126 (8750 km) | consistent, partially verified |
| J0030 M, R | 1.34 (+0.15 -0.16) M_⊙, 12.71 (+1.14 -1.19) km | same (1912.05702 abstract) | correct |
| Baumbach-Allen n_e | 10⁸(2.99x⁻¹⁶ + 1.55x⁻⁶ + 0.036x⁻¹·⁵) cm⁻³ | same, 1112.3850 Eq. (5) | coefficients correct, but **validity is 1.05-3 R_⊙**. 028 W2a uses it at ~19 R_⊙ (extrapolated). This affects only the solar plasma edge, and the Sun is not the binding object |
| Retinò et al. (028 text) | ~8×10⁻¹⁴ eV | 1.4×10⁻⁴⁹ to 3.4×10⁻⁵¹ kg = 7.9×10⁻¹⁴ to 1.9×10⁻¹⁵ eV (1302.6168) | 028 quoted the weak end of a range |

None of these changes 028's verdict on EC-06c.

**2. EC-06d vs H(z*) (`data/expansion_history_bounds.csv`, `hz_bounds.py`).**
- Allowed δH/H at z* (95%):
  - **0.18%** from Planck 2015 XIV early dark energy, Ω_e < 0.0036 (TT,TE,EE+lowP+BSH, Table 3);
  - **0.49%** from Planck 2018 ΔN_eff < 0.30 (TT,TE,EE+lowE+lensing+BAO, Sec. 7.5.2). Extra radiation is 24% of the total at z* (z_eq = 3402).
- For z > 50: Ω_e < 0.02 gives ≈ 1%.
- Correction to 025: its H_ΛCDM omits radiation. Including it lowers SNY's ratio by √(1 + ρ_r/ρ_m) = 1.15 (025's script comment says ~1.5). The ratios become **H(z*)/H_ΛCDM = 28.9 (case V), 53.8 (IV), 57.3 (III)**.
- Shortfall at z* = log₁₀[(H/H_ΛCDM - 1)/(δH/H)_max]:
  - **3.8-4.1 orders** with ΔN_eff;
  - **4.2-4.5 orders** with Ω_e;
  - 4.9-6.0 orders in energy density.
- At z = 100 (ratio 9.9-19.5) the shortfall is 2.9-3.3 orders (Ω_e(z > 50) ≲ 0.02).
- At z = 10 (ratio 3.0-5.9) I found no cited direct H bound. The z = 10 EDE limit is only in 1502.01590 Fig. 11.

**3. KK-04 vs Cassini (`kk04_shortfall.py`).**
- Bulk matter (ω = 0): |γ - 1| = 1/2 is **4.34 orders** beyond 2.3×10⁻⁵, and 3.87 orders beyond the 2σ value 6.7×10⁻⁵ (|2.1| + 2×2.3).
- Brane matter (029: 2/7): 4.09 and 3.63 orders.
- Cassini is from `data/kk_bounds_inputs.csv` (Will 1403.7377, verified at source by 029, not re-opened here).

**4. 2610.02667 (Yin & Komatsu, `literature/2610.02667.md`).**
- What it does: EB spectrum shapes for a neutrino-asymmetry current J_μA_νF̃^{μν} (a⁻³ source) vs axion DE/DM, with a Planck PR4 fit M/√(b_Jξ_νe) = 44.3 TeV.
- EC-02 is this operator with |b_J|/M² = ξκ/2. Their fit implies ξ ≈ 1.3×10²⁹ (ξ_νe at 95% bound), against CLAIMS's 3×10²⁸: the same 28-29 orders (`ec02_yinkomatsu.py`).
- No torsion/KK content. Novelty statements unchanged, so nothing was appended to `novelty_audit.md`.

## Verdict
- 028 inputs: now cited. One wrong citation (Sirius B radius) and one model used outside its validity range (Baumbach-Allen at 19 R_⊙). No verdict changes.
- EC-06d: **dead, excluded by 3.8-4.5 orders in H(z*)** (95% CMB bounds), using the radiation-corrected ratio 29-57 rather than 33-66.
- KK-04: **4.3 orders** short of Cassini (3.6-4.3 over bulk/brane matter and 1σ/2σ).
- 2610.02667: consistent with EC-02 dead. No novelty change.

## Open issues
1. The H(z*) bounds are for specific parameterisations (extra radiation, constant early-DE fraction). The SNY histories are a different expansion law, but with δH/H ~ 30-57 no parameterisation choice matters at the order-of-magnitude level. No model-independent H(z*) reconstruction or z ≈ 10 bound was found.
2. 1502.01590 is Planck 2015. I found no Planck 2018 Ω_e number in 1807.06209.
3. The 0.6 M_⊙ WD radius was checked by formula, not by a printed source value. The origin of the 028 Sirius B error ±0.000046 is unknown.
4. 028's solar plasma edge (W2a) should use a model valid beyond 3 R_⊙ (e.g. a Leblanc-type solar-wind density). It does not affect the WD-driven verdict.
5. The 025 script comment "including [radiation] raises H_ΛCDM(z*) by ~1.5" is wrong (1.15). The overseer may want 025's 33-66 updated to 29-57 in CLAIMS/draft.
