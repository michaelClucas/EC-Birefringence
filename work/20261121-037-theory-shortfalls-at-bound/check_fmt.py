"""037: vector-current-sourced axial torsion in EC + Holst + non-minimal fermion coupling (Freidel-Minic-Takeuchi,
hep-th/0507253, opened 2026-11-21 via ar5iv: Eq. 18 contortion, Eq. 23 four-fermion term).

FMT Eq. 18: axial (epsilon) part of the contortion proportional to gamma^2/(gamma^2+1) (A + (alpha/gamma) V).
FMT Eq. 23: S_int = -(3/2) pi G gamma^2/(gamma^2+1) (A^2 + (2 alpha/gamma) A.V - alpha^2 V^2), pi G = kappa/8.
Normalise to minimal EC (gamma -> infinity, alpha = 0): S = -(kappa/2) [ gamma^2/(gamma^2+1) ] (J5 + (alpha/gamma) J).
Define c (as in skeptic 029): S^0 = -(c kappa/2) n_V  =>  c = gamma alpha/(gamma^2+1).

Checks:
 1. gamma -> inf, alpha = 0 reproduces minimal EC contact 3 kappa/16 J5^2
 2. V^2 contact coefficient = (3 kappa/16) c^2 (gamma^2 + 1) >= (3 kappa/16) c^2 for real gamma
 3. the A.V cross term / (2 x A^2 term) equals the S response ratio alpha/gamma (consistency of Eqs. 18 and 23)
 4. loophole: at fixed photon-coupling product zeta r (zeta = xi gamma^2/(gamma^2+1), r = alpha/gamma),
    the V^2 contact = (3 kappa/16) gamma^4 r^2/(gamma^2+1) -> 0 as gamma -> 0. So the contact bounds c, not xi*c,
    unless the xi bound is placed on zeta (it is: the J5 part of the same operator carries zeta).
"""
import sympy as sp

g, al, k, xi = sp.symbols("gamma alpha kappa xi", positive=True)
r = sp.symbols("r", positive=True)


def report(name, ok):
    print(("PASS" if ok else "FAIL") + ": " + name)
    return bool(ok)


ok = True
pref = sp.Rational(3, 2) * (k / 8) * g**2 / (g**2 + 1)  # |FMT Eq. 23 prefactor|, pi G = kappa/8
cAA, cAV, cVV = pref, pref * 2 * al / g, -pref * al**2
ok &= report("gamma -> inf, alpha = 0: A^2 coefficient -> 3 kappa/16 (minimal EC)",
             sp.simplify(sp.limit(cAA, g, sp.oo) - 3 * k / 16) == 0)
c = g * al / (g**2 + 1)
ok &= report("|V^2 coefficient| = (3 kappa/16) c^2 (gamma^2 + 1)", sp.simplify(-cVV - 3 * k / 16 * c**2 * (g**2 + 1)) == 0)
ok &= report("hence |V^2| >= (3 kappa/16) c^2 for every real gamma (minimum as gamma -> 0 at fixed c)",
             sp.simplify((-cVV) / (3 * k / 16 * c**2) - 1 - g**2) == 0)
ok &= report("cross term / (2 A^2 term) = alpha/gamma = S-response ratio of FMT Eq. 18", sp.simplify(cAV / (2 * cAA) - al / g) == 0)
zeta = xi * g**2 / (g**2 + 1)
VV_r = sp.simplify((-cVV).subs(al, r * g))
ok &= report(f"V^2 contact at fixed r = alpha/gamma: {VV_r} -> 0 as gamma -> 0",
             sp.limit(VV_r, g, 0) == 0 and sp.simplify(VV_r - 3 * k / 16 * g**4 * r**2 / (g**2 + 1)) == 0)
ok &= report("vector-current photon coupling is (zeta r kappa/4) J.A Ft with zeta = xi gamma^2/(gamma^2+1); xi*c = zeta*r",
             sp.simplify((xi * c).subs(al, r * g) - zeta * r) == 0)
print("Conclusion: torsion-induced quark/lepton contact terms bound c = gamma alpha/(gamma^2+1), but small Immirzi gamma")
print("lets r = alpha/gamma (what beta needs at fixed zeta) grow with no contact cost: the FMT contact alone is not a bound on xi*c.")
print("ALL PASS" if ok else "SOME FAIL")
