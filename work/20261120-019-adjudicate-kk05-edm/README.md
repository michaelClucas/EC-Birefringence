# 20261120-019: adjudicate the KK-05 electron EDM (012 vs 015 vs 016)

## Question
In the KK-05 candidate mechanism (5D EC-KK on S¹, c(χ)·u∧e∧e∧R(ω), Z = 1+12c², Y = 16c³), what electron EDM does the torsion background induce for a 5D bulk Dirac electron coupled to the KK photon, and does JILA bound c? Three answers disagree:
- 012: c-odd piece is the EDM, d ≈ (3/4)ck, c ≲ 140.
- 015: same piece but canonically normalised, d → (3/4)k/√12 ≈ 8×10⁻³³ e·cm, no bound.
- 016: assignment reversed; P'-even pieces give the EDM, d = k(3/2)|1/6 - (4/3)c²|/√(1+12c²), c ≤ 188.

## Method
`adjudicate_edm.py` (sympy, exact, symbolic W_ab, c, mass phase θ; ~1 min). It shares no code with 012/015/016:
- Levi-Civita spin connection obtained by solving de^A + ω^A_B∧e^B = 0 (016 wrote it by hand).
- Explicit 5D gammas: Bjorken-Drell γ^μ, Γ⁴ = s·iγ₅. Both inequivalent reps s = ±1 are run.
- Fermion spin term from the symmetrised Lagrangian, Ψ̄(i/2){Γ^C, Ω_C}Ψ (016 used a numerical hermitian projection).
- Nonrelativistic limit: upper 2×2 block of γ⁰H, giving L ⊃ φ†(D_E σ·E + D_B σ·B)φ per unit W.
- General mass phase m e^{iθγ₅}, removed by the chiral rotation ψ → e^{-iθγ₅/2}ψ.
- Contortion response K_[ab5] = -(4/3)c²W ± cW̃ is taken as input. It is not in dispute: 009, 012 U3 and 016 D4 derived it three times.
- Literature sanity check: massive KK modes of a massless bulk fermion must have g = 1 (Ohkuwa 1984 Osaka thesis; Hosoya, Ishikawa, Ohkuwa, Yamagishi, PLB 134 (1984) 44; abstract/secondary level only).

Output is in `run.log`: 9 PASS, 0 FAIL.

## Result

### 4D reduction
- The only terms linear in F_KK that reach a zero mode come from the spin connection. The charge term -kA_a∂_y vanishes for n = 0.
- After symmetrisation only ω_[ABC] survives [A3]. The relevant pieces are ω_[ab5] = W_ab/6 (Levi-Civita, [A1]) plus K_[ab5].
- They enter as ∝ ω_[ab5] Ψ̄Γ^{ab}Γ⁵Ψ with Γ⁵ = ±iγ₅. So every such term is a σ^{ab}γ₅ structure:
  - a P'-even coefficient (W) gives iψ̄σγ₅ψF, which is electric;
  - the P'-odd coefficient (W̃) gives iψ̄σγ₅ψF̃ ∝ ψ̄σψF, which is magnetic.
- Whether this is an EDM in the lab depends on the phase of the fermion mass relative to it.

For a real mass (θ = 0), both reps and both K signs [A4]:
- **D_E = (3/2)(1/6 - (4/3)c²)**, **D_B = (3/2)c**;
- for general θ, D_E(θ) = D_E cos θ ∓ D_B sin θ.

This is exactly 016 D5/D6.

### P' bookkeeping
- P' = diag(1,-1,-1,-1,-1) has determinant +1. It is a rotation in the 4 spatial dimensions, so every 5D-covariant term, including MΨ̄Ψ, is P'-even.
- The KK photon g_{μ5} is P'-axial (A₀ → -A₀, A_i → A_i), so E → E and B → -B.
- σ·E is therefore P'-even and σ·B is P'-odd. The P'-even contortion pieces give the EDM and the P'-odd cW̃ gives the MDM.
- The 4D parity with y fixed is three reflections. It flips MΨ̄Ψ, so the 5D mass is the CP violation that makes the c = 0 term an EDM (Thirring 1972, not opened).

### Sanity check at c = 0 [A5]
For n ≠ 0 and M = 0, the mass is the KK γ₅-mass (n/R)Γ⁴, so θ = ±π/2.
- The rotated W/6 term is then purely magnetic, D_B' = ∓1/4.
- Together with the Dirac moment q/2m (q = kn/R, m = |n|/R) this gives **g = 1** for both reps and both signs of n. That is the known KK result.
- This fixes the normalisation of the W/6 term and shows it is not "intrinsically an EDM". It is an EDM for a real-mass zero mode (d = k/4 ≈ 9.5×10⁻³³ e·cm, far below JILA) and magnetic for a KK-charged mode.

### Who is right, and the erroneous steps
- **012 is wrong.** derivation.md B3 says: "The c-odd piece is P'-odd, so it is the EDM-type one."
  - That treats P' as ordinary parity with a polar photon. The KK photon is P'-axial, so the assignment is reversed.
  - It also uses (3/4)K_[ab5] counting each unordered (a,b) once, a factor-2 undercount (3/2 is correct).
  - It omits Z.
- **015 is wrong on the physics and the coefficient.** (d4) keeps 012's c-odd piece as "the EDM" and its 3/4, adds 1/√Z, and assigns the (4/3)c²W piece to an MDM. The P'-even piece is the EDM for a real-mass electron.
- **016 is correct for the computation it states.** It has the right assignment, the 3/2 normalisation, canonical normalisation, d(c) = k(3/2)|1/6 - (4/3)c²|/√(1+12c²) → 0.577ck, and **c ≤ 188** for θ = 0 [A6, reproduced: 187.6].

### Which θ applies, and the limit on 016's conclusion
The θ = 0 bound needs a KK-neutral zero mode (real mass) that nevertheless carries charge e under the KK photon. That charge coupling is not KK-gauge invariant:
- under y → y + λ(x), a y-independent ψ does not transform, so eψ̄γψA_KK is not invariant;
- KK charge is KK momentum, q = kn/R;
- 016 Sec. 6 itself keeps KK gauge invariance exact (V = ∂_y), so that background allows no such term.

With the only consistent charge (n ≠ 0), θ = arctan(n/MR) → π/2 and the roles swap: the c-odd piece becomes electric.
- d = (3/2)ck/√(1+12c²) saturates at 1.6×10⁻³² e·cm, 250 times below JILA, for every c [A7].
- So the electric piece is 015's form, with 3/2 in place of 3/4.
- Such an electron is impossible anyway: m ≥ e√Z/k ≈ 2×10²⁷ eV at c = 1 [A8].

### Final d(c, θ) in physical units
**d = (k/e)|(3/2)(1/6 - (4/3)c²) cos θ ∓ (3/2)c sin θ| / √(1+12c²)**, with k/e = 3.79×10⁻³² cm.
- θ = 0 (real-mass zero mode with an ad hoc charge): c ≤ 188 (JILA 2212.11841 Eq. 4, 90% CL).
- θ = π/2 (geometric KK charge): no bound.

## Verdict for KK-05
- 016's coefficients and EDM/MDM assignment stand. 012's assignment and 015's "the (4/3)c²W piece is magnetic" are wrong.
- "Excluded by EDM (c ≤ 188 vs clocks c ≥ 5×10⁴)" holds only if the electron is a real-mass zero mode with a non-geometric, KK-gauge-breaking charge coupling. That is the unsolved KK charge problem (012 open issue 6, 016 condition 3), which 016 lists but does not count against the bound.
- For a gauge-consistent KK charge the EDM gives no bound, but no light charged electron exists.
- Recommended CLAIMS wording: KK-05 is **excluded conditionally** ("EDM excludes the minimal form if the electron is a real-mass bulk zero mode with charge e under the KK photon; that coupling breaks KK gauge invariance"). In either branch KK-05 is not viable as a description in which SM electrons scatter CMB photons through the KK photon. It is a candidate mechanism without a consistent matter sector.

## Open issues
- K_[ab5] is taken from 009/012/016, not re-derived here.
- The g = 1 literature result was read at abstract/secondary level only (Ohkuwa thesis abstract, later citing papers). The original PLB was not opened. Thirring 1972 was not opened.
- Leading nonrelativistic order only. One-loop KK-tower contributions to d are not computed.
- A brane-localised electron (not bulk) coupled to the KK photon via a 4D charge is a different model. Its dipole would come from the brane-induced contortion coupling, which is not computed.
- No θ between 0 and π/2 is motivated by any consistent model. 016's P'-odd mass m'Ψ̄u̸Ψ is the only other source and gives the same formula.

## Rerun
```
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-019-adjudicate-kk05-edm/adjudicate_edm.py
```
