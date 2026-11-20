# Conventions

Fixed by the overseer in iteration 1. Every derivation and every imported result is converted to these. Changes are made only by the overseer and logged at the bottom.

## Units and signature
- Natural units ħ = c = k_B = 1. Masses and energies in eV unless stated.
- Metric signature (-,+,+,+). ε_{0123} = +√(-g) (so ε^{0123} = -1/√(-g)).
- κ = 8πG = 1/M_P², with M_P = 2.435×10¹⁸ GeV the reduced Planck mass.
- FRW: ds² = a²(τ)(-dτ² + dx²), conformal time τ, a(today) = 1. H = ȧ/a in cosmic time.

## Connection, torsion, contortion
- Covariant derivative on vectors: ∇_μ V^λ = ∂_μ V^λ + Γ^λ_{μν} V^ν (derivative index first, lower-left).
- Torsion: T^λ_{μν} = Γ^λ_{μν} - Γ^λ_{νμ} (antisymmetric in the last two indices).
- Contortion K^λ_{μν} defined by Γ^λ_{μν} = {λ μν} + K^λ_{μν} (Christoffel plus contortion), with metric compatibility (K_{λμν} = -K_{νμλ}). Then K_{λμν} = ½(T_{λμν} + T_{νλμ} - T_{μνλ}) with index order (λ, μ, ν) matching Γ^λ_{μν}. Checked symbolically in work/20261120-003-theory-ec-photon/check_contortion.py and work/20261120-004-kk-reduction (both PASS).
- Irreducible decomposition:
  - trace vector T_μ = T^λ_{μλ}
  - axial vector S^μ = (1/6) ε^{μνρσ} T_{νρσ}
  - pure tensor part q^λ_{μν}: the remainder (traceless, ε-contraction zero).
- FRW-compatible torsion (homogeneous, isotropic): only T_{ijk}-type axial part S^0 ≡ S_0(τ) and the trace part with T_{i0j} ∝ h(τ) δ_{ij}. Any other component must be justified.

## Electromagnetism and birefringence
- F_{μν} = ∂_μ A_ν - ∂_ν A_μ. Maxwell Lagrangian L = -¼ F_{μν}F^{μν}.
- Dual: F̃^{μν} = ½ ε^{μνρσ} F_{ρσ}.
- Chern-Simons baseline: L ⊃ -(g_{φγ}/4) φ F_{μν} F̃^{μν}. Carroll-Field-Jackiw: L ⊃ -½ p_μ A_ν F̃^{μν}. The two are equivalent up to a total derivative with **p_μ = -g_{φγ} ∂_μ φ** (corrected iter 1; checked in work/20261120-003-theory-ec-photon/check_cfj_rotation.py).
- Rotation angle β: the sign convention of arXiv:2509.13654 (Minami-Komatsu): β > 0 is clockwise on the sky as seen by the observer (IAU). Wording and polarization convention in `literature/conventions_beta.md`. Note: Komatsu's review uses ε^{0123} = +1/√(-g), opposite to ours, so formulas imported from it flip sign.
- **In our conventions:** a photon with CFJ vector p gives β = +½∫p_τ dτ (conformal time, last scattering to today), and L ⊃ -(g/4)φFF̃ gives **β = -(g_{φγ}/2)(φ_0 - φ_LSS)**. Derived and checked numerically (field evolution) in work/20261120-003-theory-ec-photon (derivation Sec. on CFJ rotation).
- Angles: β is reported in degrees in `data/`; all code works in radians internally and converts on input/output.

## Kaluza-Klein track
- 5D indices M, N = 0..3, 5; 4D indices μ, ν. Fifth coordinate y ∈ [0, 2πR) for S¹.
- 5D signature (-,+,+,+,+). Same torsion/contortion definitions with 5D indices.
- KK ansatz is stated explicitly in each derivation (including whether the radion is kept and in which frame, Einstein or Jordan). Conversion to 4D quantities above must be shown.

## Changelog
- iter 1: initial version.
- iter 1 (audit): contortion formula corrected (old form failed metric compatibility; only the axial part had been used, so nothing downstream changed). CFJ ↔ φFF̃ sign corrected to p = -g∂φ. β↔Δφ relation fixed with sign. FRW torsion: trace T_τ = 3h, axial S^τ = -f/a² (work/20261120-003).
