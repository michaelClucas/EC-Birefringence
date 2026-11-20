# Sign convention for the birefringence angle beta

Sources read in full text: arXiv:2509.13654 (primary), arXiv:2601.13624, arXiv:2512.21888, and the convention references they point to: Komatsu 2022 review arXiv:2202.13919 (cited by 2601.13624 for its angle convention) and Eskilt & Komatsu arXiv:2205.13962. Minami & Komatsu arXiv:2011.11254 was searched and contains no explicit sign statement.

## 1. Sign of beta on the sky

- arXiv:2509.13654, Sec. I: "the rotation can only be clockwise ( beta>0 ) or counterclockwise ( beta<0 )". So positive beta = clockwise on the sky.
- The paper does not say "as seen by the observer" and does not name IAU or COSMO/HEALPix. Its mixing matrix (Eq. 4) gives E^o = cos(2psi) E - sin(2psi) B, B^o = sin(2psi) E + cos(2psi) B with psi = alpha + beta, i.e. (E + iB)^o = e^{2i psi}(E + iB) at zero leakage, and C^EB,o = (1/2) sin(4 beta)(C^EE - C^BB), C^TB,o = sin(2 beta) C^TE.
- Cross-checks:
  - arXiv:2205.13962 (Eskilt & Komatsu), Sec. I after Eq. 3: same convention, positive beta is a clockwise rotation of the plane of linear polarization on the sky, with E^o = E cos 2beta - B sin 2beta, B^o = E sin 2beta + B cos 2beta (same rotation as 2509.13654 Eq. 4).
  - arXiv:2601.13624 (Yin et al.), text after Eq. 4: positive beta is clockwise on the sky, with the z-axis along the observer's line of sight; follows Komatsu 2022. Their Eq. 10, Delta_E +- i Delta_B = e^{+-2i beta}(...), matches 2509.13654 Eq. 4.
  - arXiv:2512.21888: no clockwise/counterclockwise statement. Its beta = (1/2) g (phi_obs - phi_LSS) has the same form as 2601.13624, so it is presumably the same convention (not stated by the authors).
  - arXiv:2510.07928: only anisotropic power C_L^{alpha alpha}, sign-independent.

## 2. IAU vs CMB (COSMO/HEALPix) polarization convention

Taken from the Komatsu 2022 review (arXiv:2202.13919, Sec. "Rotation of the plane of linear polarisation" / Stokes parameter discussion), which 2601.13624 cites as its convention source and which is by the same group as 2509.13654:
- IAU: position angle psi with Q +- iU = P e^{+-2i psi}; psi > 0 is counterclockwise on the sky.
- CMB community: right-handed coordinates with z along the observer's line of sight (pointing toward the sky), which flips the sign of U and of the position angle relative to IAU.
- The review adopts the CMB convention and writes beta = -psi_IAU. Then beta > 0 is clockwise as seen by the observer, and Q +- iU = P e^{+-2i beta} = Q_IAU -+ i U_IAU.
- Our reading (not a phrase used in any of these papers): this "CMB convention" is the COSMO convention that HEALPix uses by default (U_COSMO = -U_IAU). None of the four target papers names COSMO or HEALPix explicitly. The overseer should treat "Minami-Komatsu convention = COSMO/HEALPix, beta = -psi_IAU, clockwise positive" as sourced to 2202.13919, not to 2509.13654.

## 3. beta vs change of the pseudoscalar field

| Source | CS term | Stated beta formula | Notes |
|---|---|---|---|
| 2202.13919 Eq. 7-8 | I_CS = int d^4x sqrt(-g) ( -(alpha/4f) chi F Ftilde ), F Ftilde = sum F_{mu nu} epsilon^{mu nu mu' nu'} F_{mu' nu'} / (2 sqrt(-g)), symbol epsilon^{0123} = +1 | beta = (alpha/2f) [chi(eta_0) - chi(eta_LS)] (Sec. text after Eq. 12) | metric a^2 diag(-1, 1); F Ftilde = -4 B.E. Also: if L_CS = +(1/4) g_a chi F Ftilde, then beta = -(1/2) g_a int chi' d eta |
| 2601.13624 Eq. 1, 4 | L ⊃ -(1/4) g phi F_{mu nu} Ftilde^{mu nu} | beta(t) = (g/2)[phi(t0) - phi(t)] | signature and epsilon sign not stated |
| 2512.21888 Eq. 1, 2 | S ⊃ int sqrt(-g) ( -(1/4) g_{phi gamma} phi F_{mu nu} Ftilde^{mu nu} ) | beta = (1/2) g_{phi gamma} (phi(t0, 0) - phi(t_LSS, d_LSS nhat)) | FLRW ds^2 = -dt^2 + a^2 dx^2; epsilon sign not stated |

Summary in their convention: L ⊃ -(g/4) phi F Ftilde gives beta = +(g/2)(phi_today - phi_emission), factor 1/2, positive = clockwise (CMB/COSMO convention).

## 4. Warning for the repo conventions (for the theorist; not settled here)

- Komatsu 2022 uses epsilon^{0123} = +1 as a symbol divided by sqrt(-g), i.e. tensor epsilon^{0123} = +1/sqrt(-g).
- `CONVENTIONS.md` fixes epsilon_{0123} = +sqrt(-g), so epsilon^{0123} = -1/sqrt(-g). This is the opposite sign for the upper-index tensor.
- So Ftilde^{mu nu} in our convention is minus Komatsu's, and our baseline L ⊃ -(g/4) phi F Ftilde equals +(g/4) phi F Ftilde in Komatsu's convention. By the review's own alternative formula this would give beta = -(g/2) Delta phi in our convention. This needs an explicit derivation in `work/`; it is flagged, not concluded. 2601.13624 and 2512.21888 do not state their epsilon sign, so their formulas can't be converted without assuming one.
