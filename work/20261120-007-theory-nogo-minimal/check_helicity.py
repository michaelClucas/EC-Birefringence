"""Axial charge density of fermions with definite helicity (derivation Sec. A1).

Claims checked (Dirac representation, signature (-,+,+,+) so gamma^0 squares to -1 is avoided by
working with the Hermitian matrices beta = gamma^0_(+---) and alpha_i; J5^0 = psi^dag gamma5 psi):
 1. For a positive-energy plane wave of momentum p and helicity h = +-1:
       u^dag gamma5 u / u^dag u = h |p| / E          (axial charge per particle = h v)
 2. Massless limit: gamma5 eigenvalue = h (chirality = helicity).
 3. Charge conjugation: C (gamma^mu gamma5)^T C^-1 = + gamma^mu gamma5 while C gamma^mu^T C^-1 = -gamma^mu,
    so the axial current is C-even: an antiparticle of helicity h contributes +h v like a particle.
 4. Free Dirac Hamiltonian in FRW (conformal time, rescaled spinor): H = alpha.p + a m beta.
    The helicity operator Sigma.p commutes with H for any a(tau): helicity is conserved for
    free-streaming relic fermions, while v = |p|/E decays as the momentum redshifts.
 5. gamma5 commutes with Sigma.p, so the axial-torsion coupling S_0 psi^dag gamma5 psi does not flip helicity.
"""
import sympy as sp

from common import report

I2 = sp.eye(2)
Z2 = sp.zeros(2)
sx = sp.Matrix([[0, 1], [1, 0]])
sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sz = sp.Matrix([[1, 0], [0, -1]])
sig = [sx, sy, sz]


def blk(a, b, c, d):
    return sp.Matrix(sp.BlockMatrix([[a, b], [c, d]]))


beta = blk(I2, Z2, Z2, -I2)
alpha = [blk(Z2, s, s, Z2) for s in sig]
g5 = blk(Z2, I2, I2, Z2)
Sigma = [blk(s, Z2, Z2, s) for s in sig]
# gamma matrices in the (+,-,-,-) Bjorken-Drell form; gamma5 = i g0 g1 g2 g3
g0 = beta
gi = [beta * a for a in alpha]
gam = [g0] + gi


def small(M, tol=1e-25):
    """Exact expression evaluated to 40 digits; all entries below tol."""
    return all(abs(sp.N(x, 40)) < tol for x in M)


def main():
    ok = True
    ok &= report("gamma5 = i g0 g1 g2 g3", sp.simplify(sp.I * g0 * gi[0] * gi[1] * gi[2] - g5) == sp.zeros(4))

    good1 = good2 = True
    # exact (m, p) pairs with rational E, including the massless case
    # exact test directions (rotation invariance makes three generic ones sufficient)
    dirs = [(sp.pi / 3, sp.pi / 4), (sp.Rational(2, 7) * sp.pi, sp.Rational(5, 3) * sp.pi), (sp.pi / 2, 0)]
    for (th, ph), h, (m, p) in [(d, hh, mp) for d in dirs for hh in (1, -1)
                                for mp in ((3, 4), (5, 12), (0, 1))]:
        m, p = sp.Integer(m), sp.Integer(p)
        E = sp.sqrt(p**2 + m**2)
        n = [sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)]
        sn = sum((n[i] * sig[i] for i in range(3)), sp.zeros(2))
        # 2-spinor helicity eigenvector chi_h of sigma.n
        if h == 1:
            chi = sp.Matrix([sp.cos(th / 2), sp.exp(sp.I * ph) * sp.sin(th / 2)])
        else:
            chi = sp.Matrix([-sp.exp(-sp.I * ph) * sp.sin(th / 2), sp.cos(th / 2)])
        good1 &= small(sn * chi - h * chi)
        u = sp.Matrix.vstack(sp.sqrt(E + m) * chi, sp.sqrt(E - m) * h * chi)
        ud = u.conjugate().T
        Hd = sum((alpha[i] * p * n[i] for i in range(3)), sp.zeros(4)) + m * beta
        good1 &= small(Hd * u - E * u)
        ratio = (ud * g5 * u)[0] / (ud * u)[0]
        good2 &= small(sp.Matrix([ratio - h * p / E]))
    ok &= report("u(p,h) solves the free Dirac equation", good1)
    ok &= report("u^dag g5 u / u^dag u = h|p|/E, and = h for m = 0", good2)

    C = sp.I * gam[2] * gam[0]
    Ci = C.inv()
    v_odd = all(sp.simplify(C * g.T * Ci + g) == sp.zeros(4) for g in gam)
    a_even = all(sp.simplify(C * (g * g5).T * Ci - g * g5) == sp.zeros(4) for g in gam)
    ok &= report("C gamma^T C^-1 = -gamma (vector current C-odd)", v_odd)
    ok &= report("C (gamma gamma5)^T C^-1 = +gamma gamma5 (axial current C-even)", a_even)

    a, m = sp.symbols("a m", positive=True)
    px, py, pz = sp.symbols("p_x p_y p_z", real=True)
    pv = [px, py, pz]
    H = sum((alpha[i] * pv[i] for i in range(3)), sp.zeros(4)) + a * m * beta
    hel = sum((Sigma[i] * pv[i] for i in range(3)), sp.zeros(4))
    ok &= report("[Sigma.p, alpha.p + a m beta] = 0 (helicity conserved in FRW)",
                 sp.simplify(hel * H - H * hel) == sp.zeros(4))
    ok &= report("[gamma5, Sigma.p] = 0 (axial torsion does not flip helicity)",
                 sp.simplify(g5 * hel - hel * g5) == sp.zeros(4))
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
