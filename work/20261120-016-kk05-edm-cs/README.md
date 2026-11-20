# 20261120-016 KK-05 electron EDM, P'-odd background, CS level

The ACT DR6 β = 0.215° ± 0.074° (2509.13654) is 2.9σ with systematics caveats and may not be cosmological. Everything here is a candidate mechanism.

## Question
1. **KK-05.** χ is coupled through c(χ)·u∧e∧e∧R(ω) in 5D EC-KK. How does it reach the electron EDM? Give a definite JILA bound on c. Does a rolling χ act as a CP-odd background today (an oscillating EDM)? Verdict.
2. Is there a non-dynamical P'-odd background that leaves the KK photon massless and 4D Lorentz invariant?
3. KK-02: is the CS level quantised in Z or 6Z? Give β_max for the minimal level.

The coordinator relayed skeptic 015's objection: with Z = 1+12c² included, d = (3/4)ck/√Z ≈ 8×10⁻³³ e·cm, so there is no EDM bound. Checked independently below.

## Method
Derivations are in `derivation.md`. Rerun from the repo root (about 1.5 min, one core):
```
.venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-016-kk05-edm-cs/run_checks.py
```

| script | content |
|---|---|
| `check_dipole.py` | Exact sympy elimination of the 50 contortion components, written independently of 009/012. Reproduces Z = 1+12c², Y = 16c³ and K_[ab5]. Builds the Levi-Civita KK connection, explicit 5D gamma matrices, the hermitian spin-connection term and rest-frame EDM/MDM coefficients. Also covers the chirality check and the rolling-χ axial background. |
| `edm_bound.py` | c_EDM from JILA. MICROSCOPE and clock bounds (same construction as 012). θ-phase escape. b₀ from the rolling χ. EDM drift vs the oscillating-EDM band. |
| `check_background.py` | Fixed fibre vector V = ∂_y: KK-gauge invariance, 4D diffeos, P'-oddness, which generators are broken. |
| `check_cs_level.py` | Index theorem on CP³ (Â vs HRR), k ∈ Z vs 6Z, gravitational CS partner, ϑ period, β_max. |

New data: `data/kk05_edm_inputs.csv` (HfF⁺ oscillating-EDM band, 2006.15787 abstract).

Last run (`run.log`): **28 PASS, 0 FAIL**.

## Result
1. **The 012 coefficients hold** [D3, D4], from an independent exact computation:
   - Z = 1+12c², Y = 16c³, so dβ = -c dlnα stands.
   - K_[ab5] = -(4/3)c²W ± cW̃.
2. **The EDM is the P'-even part, not the c-odd part** [D5, D6].
   - Per unit W = kF: D_E = (3/2)|1/6 - (4/3)c²| and D_B = (3/2)c.
   - The c-odd cW̃ term is a **magnetic** moment, because the KK photon is P'-axial.
   - Canonically normalised: d = k(3/2)|1/6 - (4/3)c²|/√(1+12c²) → 0.577ck = 2.2×10⁻³² c e·cm. This **grows linearly in c**.
   - Skeptic 015's formula is the MDM piece and saturates at 8.2×10⁻³³ e·cm [E1]. 015 dropped the c² term, so its "no bound" conclusion does not hold.
   - 012 had the wrong piece and no Z. Its 140 was close to the right number by accident.
3. **JILA bound:** **c ≤ 188** [E2]. This holds for an electron that is:
   - a 5D bulk Dirac zero mode, which is necessarily vector-like while the KK photon is massless (D7b; Witten 1981);
   - given a P'-even bulk mass (θ = 0);
   - coupled to the KK photon with charge e.
   
   For comparison, plain KK (c = 0) gives 9.5×10⁻³³ e·cm.
4. **Other bounds on c** [E3, E4]:
   - MICROSCOPE: c ≥ 340 at Δχ = 0.39 M_P, c ≥ 133 at 1 M_P.
   - Clocks: c ≥ 5.1×10⁴.
   - **Window closed:** the clocks exceed the EDM bound by a factor of 273. EDM + MICROSCOPE alone close it for Δχ < 0.71 M_P.
5. **Escape requires a tuned phase** [E5]. The electron mass must come ≥ 99.6% from the P'-odd term m′Ψ̄u̸Ψ (|cos θ| < 3.6×10⁻³). That term is allowed only by the background and is non-minimal.
6. **Rolling χ** [D8, E6, E7]:
   - χ itself is P'-even.
   - It induces an axial background b₀ = 1.5c′χ̇ ≈ 8×10⁻³⁶ eV for all bulk fermions, 18 orders below Eot-Wash.
   - The EDM drifts by ≲ 10⁻¹⁵ per year, far below the 27 nHz bottom of the HfF⁺ oscillation band, so it acts as a static EDM.
   - No oscillating-EDM bound applies. The χFF̃ background alone is P-odd and T-even, so it gives no EDM.
7. **P'-odd background: yes** [G1-G4].
   - A fixed fibre vector V = ∂_y (prior geometry) gives u = V♭/|V| = e⁵.
   - It is invariant under KK U(1) gauge transformations, so the KK photon stays massless at all orders.
   - It is invariant under 4D diffeos, so 4D Lorentz symmetry is unbroken.
   - It is odd under P'.
   - Cost: explicit breaking of the y-dependent 5D diffeos (Kostelecký hep-th/0312310 consistency, unchecked for massive KK modes). In the zero-mode EFT it is just a fixed-sign parity-violating coupling.
   - I found no literature using this construction.
8. **CS level** [C1-C6]:
   - On spin manifolds, k ∈ Z with the mixed term -(k/192π²)A∧tr R∧R, and k ∈ 6Z without it. The obstruction example is CP³ with O(1).
   - A charged bulk fermion shifts k_bare by ½q³ and induces the mixed term, so k ∈ Z is the relevant case.
   - **Minimal |k| = 1: β_max = 0.209°**, needing Δϑ = 1.03π. If only 6Z were allowed: 1.25°, needing Δϑ = 0.17π.

## Verdict
- **KK-05 (candidate mechanism):** **excluded in its minimal form.** With a bulk Dirac electron and a P'-even mass, the EDM gives c ≤ 188 while the clocks need c ≥ 5.1×10⁴, a factor of 273. MICROSCOPE alone gives a factor of 1.8 at Δχ = 0.39 M_P. Outside the minimal form it is **strongly disfavoured**. Surviving requires **all** of the following:
  - a non-dynamical P'-odd background (V = ∂_y prior geometry);
  - an unsolved KK charge mechanism making the KK photon the CMB photon;
  - either an electron mass tuned to ≥ 99.6% P'-odd, or an electron that is not a bulk fermion (then its dipole is not fixed);
  - c ≳ 5×10⁴ (Z ~ 3×10¹⁰), or a potential that stops χ before today together with Δχ > 0.71 M_P (window 133 ≤ c ≤ 188 at Δχ = M_P).
- **KK-02:** β_max = 0.209° at |k| = 1 is confirmed as the minimal level (spin manifolds, gravitational CS term present).

## Open issues
1. The overall signs of D_E and D_B and of β relative to c are not tracked. Only magnitudes are used. The relative sign of ω̊ and K matters only near c ≈ 0.35, where D_E dips; the bound is at c ≫ 1.
2. The bound assumes the electron's charge coupling to the KK photon. That coupling is P'-odd for a zero mode and has no known origin (the KK charge problem). If the CMB photon is not the KK photon, KK-05 gives no β at all.
3. The θ escape is parametrised, not modelled. A radiatively generated m′ (gravitational, P'-odd) has not been estimated.
4. Explicit breaking of y-dependent diffeos (Sec. 6): Bianchi-identity consistency in the massive KK sector is not checked.
5. References not opened: Thirring 1972 (5D CP violation), Witten 1981 (no chiral KK fermions), hep-th/9702198 and hep-th/9603150 (half-integer CS shifts), hep-th/0312310 (explicit diffeo breaking). The CS quantisation is derived here from the index theorem [C2, C3], not quoted.
6. Errors in other folders (not edited):
   - 012 Sec. B3 and `bounds_uholst.py` B4 assign the EDM to the c-odd piece and omit 1/√Z. The correct coefficient is 0.577ck (c_EDM = 188, not 144). The verdict "window closed" survives.
   - Skeptic 015's "no EDM bound" uses only the MDM piece.
   - CLAIMS/STATE KK-05 "EDM estimate c ≲ 140" should become "c ≤ 188 (JILA, bulk Dirac electron, θ = 0)".
   - 009 open issue 1 is resolved here.

## Files
`common.py`, `check_dipole.py`, `edm_bound.py`, `check_background.py`, `check_cs_level.py`, `run_checks.py`, `run.log`, `derivation.md`. New data: `data/kk05_edm_inputs.csv` (plus a section in `data/README.md`).
