# 20261120-028 skeptic: EC-06c repair costs and EC-06d (audit of task 025), iteration 8

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism. There is no `RESULT.md`. I did not do task 025.

## Question
Does 025 hold up on: (1) the f(T) Lorentz-scalar ζ solution and dT ≠ 0 on solutions; (2) gravitational CS (0907.5159) as a non-repair; (3) the Stückelberg window 7×10⁻²⁸ ≲ m_γ ≲ 10⁻²⁶ eV (the A < 0.1 criterion, plasma, the Proca galactic limit); (4) the naive EFT cutoff of 20-60 GeV and whether it conflicts with photon physics; (5) EC-06d: the H(z*) exclusion and ξ ≳ 1 in the GR-compatible regime?

## Method
- `eta_sector.py` (sympy): compares 025's ζ solution with Golovnev-Koivisto 1808.05565 Eq. (6) (paper HTML opened, Sec. 7, Eqs. 6 and 9). Rebuilds the boost-sector second-order action and extracts the η-η part (the stiffness of the Lorentz mode η = ∂_xζ) and the CFJ source in the η equation. It uses only the generic torsion helpers from 025's `common.py`. The action assembly and the extraction are new.
- `window_cutoff.py`: lower edge of the photon-mass window from photons leaving the Sun, white dwarfs and a neutron star (025's moving-source mechanism, unchanged); plasma detuning (Baumbach-Allen corona, WD exosphere scale height); the four window choices; 025's dim-7 cutoff; a new cutoff from η exchange.
- Inputs: `data/` plus `inputs_028.csv` in this folder (compact-object parameters, coronal model, each cited; values recalled, not re-opened; the overseer should move them to `data/` if they are kept).
- Literature opened: 1808.05565 (HTML, Sec. 7-8), hep-ph/0306245 (abstract).
- Rerun from the repo root, about 30 s each, 1 thread:
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-028-skeptic-ec06c-costs/eta_sector.py
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-028-skeptic-ec06c-costs/window_cutoff.py
```
`run.log`: 15 PASS, 0 FAIL.

## Results

| claim (025) | verdict | reason |
|---|---|---|
| 1. ζ fixed non-dynamically by the f(T) antisymmetric equation, ∇²ζ = -3Ψ′ + 3(a″/a′)Ψ - 3(a′/a)(Φ+2Ψ) | **SURVIVES** | Identical to GK 1808.05565 Eq. (6) (E1). GK also find ζ constrained, not propagating. One correction: GK Eq. (9) gives slip Φ - Ψ = -12 f_TT H(H′-H²)ζ/f_T, while 025 set Φ = Ψ. That is an O(f_TT T̄) correction and does not change dT ≠ 0. |
| 1. dT ≠ 0 on solutions, U(1) broken everywhere | **SURVIVES** | It follows from ζ ≠ 0 and 025's δT_μ (reproduced in E3). |
| 1/open issue 1. photon back-reaction on the frame | **resolved classically, not in quantum processes** | η-η action (E2c): (M_P²/2)·8 f_TT H² (∂_xη)² = (2/3)α M_P²(∂_xη)², with α ≡ f_TT T̄. It has only spatial gradients, so η is a constraint field. TEGR has zero stiffness (E2b). A 1 W/mm² circularly polarised beam (K ~ 10⁻⁴ eV³) induces δT ~ ξK/(αM_P²) ~ 10⁻⁶²/α eV, against ~10⁻²⁶ eV from Earth's moving-source term. So classical back-reaction is negligible. The same small stiffness gives the low cutoff in item 4. Rotations have zero quadratic stiffness (025 L1a shows the FRW action is exactly ω-independent), so the photon source in the ω equation is balanced only through mixing. Not quantified. |
| 2. gravitational CS is not a repair | **SURVIVES** | (θ₂/4)RR̃ contains no A_μ, so it cannot enter the Maxwell equation or its divergence identity ∂E = -¼(dp)·F̃. |
| 3. lower edge 7×10⁻²⁸ eV (Sun, A < 0.1) | **BROKEN** | The same mechanism applied to white dwarfs gives m ≳ 2.5×10⁻²⁴ eV (typical 0.6 M_⊙ WD) and 1×10⁻²³ eV (Sirius B) at A < 0.1, or 2.5×10⁻²⁵ to 1×10⁻²⁴ eV at A < 1 (W1). Plasma does not help for WDs: the detuned layer is ~4 km of an 8700 km radius (W2b). For the Sun, plasma moves the vacuum-mixing region out to ~19 R_⊙ and lowers the solar edge to 2×10⁻³⁰ eV (W2a). So plasma lowers the solar edge but does not remove the bound, and the Sun is not the binding object anyway. A < 0.1 is not tied to any observation. A test that could justify it is WD fluxes against Gaia-parallax radii (1611.00629, few-% agreement, not opened), which supports P = A² ≲ 0.1. The WD conclusion holds for both A < 0.1 and A < 1. A neutron star (J0030, vacuum mixing, magnetosphere ignored) would give 2.7×10⁻¹⁸ eV, above the PDG limit. QED and plasma detuning in the magnetosphere were not evaluated. |
| 3. upper edge 10⁻²⁶ eV (Proca galactic, hep-ph/0306245) | **WEAKENED** (model dependent) | ADG's abstract: the bound holds only if the galactic field is in the Proca regime. For a Higgs mass, vortices invalidate it. A pure Stückelberg mass is in the Proca regime, but a UV completion below the cutoff (item 4) could be Higgs-like. PDG lists the bound but does not adopt it. The PDG limit 10⁻¹⁸ eV (Ryutov) is itself disputed: Retinò et al. 2016 (Cluster) argue for a weaker ~8×10⁻¹⁴ eV (recalled, not opened). |
| 3. window "open by about one decade" | **BROKEN** | With the Proca galactic limit the window is **closed** under both A criteria (W3). With PDG only it is open from ~10⁻²⁴ to 10⁻²³ eV up to 10⁻¹⁸ eV. A neutron-star bound could close that too. The window's existence depends on choices 025 did not make: which upper limit, which compact objects, and the magnetosphere. |
| 4. naive cutoff (mM_P²/ξ)^{1/3} ≈ 20-60 GeV | **SURVIVES as arithmetic, but it is not the lowest scale** | Reproduced as a dim-7 operator: π + γ + stress-energy via gravitational exchange, amplitude ~ ξE³/(mM_P²). It is 23-55 GeV, or 53-128 GeV with (4π)^{1/3} (C1). At the window shifted to 10⁻²³ eV it is ~0.6-1.3 TeV. |
| 4. new: η (Lorentz-mode) exchange | **BREAKS EC-06c** | The vertex -(ξ/2m)π(η_xx - η_ττ)F̃^{τx} is exchanged through η's stiffness (2/3)αM_P² and gives a ππγγ contact term with amplitude ~ 3ξ²E⁴/(8αM_P²m²). Unitarity gives Λ_η = (64π αM_P²m²/3ξ²)^{1/4} (C2). That is **0.2-0.7 keV** in 025's window (α = 1; lower for α < 1, Λ ∝ α^{1/4}), **24 keV** at 10⁻²³ eV and **7 MeV at the PDG limit**. Λ_η ≥ 10 TeV would need m ≈ 2×10⁻⁶ eV, which is 10¹² above PDG. Above Λ_η, γγ → ππ (invisible) is strongly coupled. The EFT then predicts unitarity-saturated photon disappearance in keV-TeV photon collisions. TeV blazar opacity on the EBL and CMB (√s ~ 0.1-2 MeV) and LEP/LHC γγ physics see none of this. So the theory needs a UV completion of the photon/f(T)-Lorentz sector below ~keV-MeV, and none is known. This is the f(T) strong-coupling problem (GK Sec. 8) made photon-observable by the CFJ coupling. |
| 5. EC-06d: SNY histories have H(z*) = 33-66 × ΛCDM | **SURVIVES** | `sny_background.csv`: case III has H/H_ΛCDM = 65.8 at z* and already 5.9 at z = 10, so the exclusion does not rest on z* alone. |
| 5. GR-compatible SNY needs ξ ≳ 1 | **SURVIVES** (algebra checked; δ = 3bρ/8 not re-checked against SNY) | β = (3ξ/4)Δln Y ≈ (3ξ/4)s*, s* = 4δ/(9μ), so ξ = 3μβ/δ ≥ 1.1 (δ ≤ 1%) or ≥ 0.11 (δ ≤ 10%). Any plausible δ gives ξ ≫ 3.6×10⁻⁴. BBN would likely tighten δ (raising ξ). Not evaluated. |

## Verdict
- **EC-06c with the Stückelberg repair: dead.** Two independent problems. (i) The photon-mass window is closed by white-dwarf photon conversion if the Proca galactic limit is used. Without it, only ~10⁻²⁴-10⁻¹⁸ eV remains, possibly closed by neutron stars. (ii) The overseer's kill test: the Lorentz-mode (η) exchange gives a cutoff Λ_η ≲ 7 MeV for every allowed m (≲ 1 keV in 025's window), far below photon energies where photon physics is observed to be standard. The 20-60 GeV estimate in 025 is correct as arithmetic, but it misses this lower scale. Caveat: the cutoff estimate is power counting from the linearised boost sector on FRW. O(1) factors and the rotation sector (zero quadratic stiffness, likely worse) are not tracked. Nothing found lifts it by the ≥ 10¹⁰ needed.
- (c2) the torsion-dragged scalar was not audited here. Its -(M²/2)(∂θ - T)² term adds an η stiffness ~M²(∂η)², small next to (2/3)αM_P² for M ≪ M_P. It also couples θ fluctuations to η, so it may share the η strong coupling. Not evaluated. It is an ALP-class model, as 025 says.
- EC-06d: dead as O(H) torsion, as 025 says.
- Items 1 and 2 survive. Item 5 survives.

## Open issues
1. Λ_η uses the 1D scalar boost sector. The vector boosts and the rotations were not done. Rotations have no quadratic stiffness on flat FRW, so their contribution may be stronger coupling still.
2. Neutron-star conversion (2.7×10⁻¹⁸ eV at A < 0.1) needs magnetospheric QED and plasma detuning before it can be used.
3. WD and NS parameters and the Retinò et al. limit were recalled, not opened. `inputs_028.csv` should move to `data/` with the sources checked.
4. The A criterion for WDs (P = A² ≲ 0.1 from Gaia radii) is argued, not fitted. The conclusion holds even at A < 1.
5. SNY Friedmann correction δ = 3bρ/8 not re-derived from 0805.3834; BBN not evaluated.
