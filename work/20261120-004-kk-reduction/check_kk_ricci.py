"""Check 1: KK reduction of the 5D Ricci scalar (y-independent zero modes).

Ansatz (Einstein frame, radion kept):
  ds5^2 = e^{2 a phi} g_mn dx^m dx^n + e^{2 b phi} (dy + A_m dx^m)^2,
  a = 1/(2 sqrt 3), b = -2a.
Claim:
  sqrt(-g5) R5 = sqrt(-g) [ R - 1/2 (dphi)^2 - 1/4 e^{-6 a phi} F^2 ] + c sqrt(-g) Box phi
with c a constant (the Box term is a total derivative).

Method: the identity is pointwise in the 2-jets of (g_mn, phi, A_m). At each test
point we draw random jets (value, gradient, symmetric Hessian), build the 5D metric
jet with the product/chain rule, and compute both sides with geom.py. c is fitted
by least squares over several points; PASS needs a residual at round-off level.
"""
import numpy as np

from geom import christoffel, ricci_scalar

rng = np.random.default_rng(1)
ETA4 = np.diag([-1.0, 1, 1, 1])
AL = 1 / (2 * np.sqrt(3))
BE = -2 * AL


class Jet:
    """Second-order jet of a scalar in 4 variables (y-derivatives are zero)."""

    def __init__(self, v, d, dd):
        self.v, self.d, self.dd = float(v), np.asarray(d, float), np.asarray(dd, float)

    @staticmethod
    def rand(scale, v0=0.0):
        h = rng.normal(size=(4, 4)) * scale
        return Jet(v0 + rng.normal() * scale, rng.normal(size=4) * scale, h + h.T)

    def __add__(self, o):
        o = o if isinstance(o, Jet) else Jet(o, np.zeros(4), np.zeros((4, 4)))
        return Jet(self.v + o.v, self.d + o.d, self.dd + o.dd)

    __radd__ = __add__

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return Jet(self.v * o, self.d * o, self.dd * o)
        return Jet(self.v * o.v, self.d * o.v + self.v * o.d,
                   self.dd * o.v + np.outer(self.d, o.d) + np.outer(o.d, self.d) + self.v * o.dd)

    __rmul__ = __mul__


def jexp(f, k):
    """exp(k f)"""
    e = np.exp(k * f.v)
    return Jet(e, e * k * f.d, e * (k * f.dd + k * k * np.outer(f.d, f.d)))


def metric_arrays(M, n):
    g = np.zeros((n, n)); dg = np.zeros((n, n, n)); ddg = np.zeros((n, n, n, n))
    for i in range(n):
        for j in range(n):
            J = M[i][j]
            g[i, j] = J.v
            dg[:4, i, j] = J.d
            ddg[:4, :4, i, j] = J.dd
    return g, dg, ddg


def one_point():
    g4 = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            g4[i][j] = g4[j][i] = Jet.rand(0.1, ETA4[i, j])
    phi = Jet.rand(0.4)
    A = [Jet.rand(0.4) for _ in range(4)]
    ea, eb = jexp(phi, 2 * AL), jexp(phi, 2 * BE)
    M5 = [[None] * 5 for _ in range(5)]
    for m in range(4):
        for n in range(4):
            M5[m][n] = ea * g4[m][n] + eb * A[m] * A[n]
        M5[m][4] = M5[4][m] = eb * A[m]
    M5[4][4] = eb
    g5, dg5, ddg5 = metric_arrays(M5, 5)
    gi5, G, dG = christoffel(g5, dg5, ddg5)
    lhs = np.sqrt(-np.linalg.det(g5)) * ricci_scalar(gi5, G, dG)
    g, dg, ddg = metric_arrays(g4, 4)
    gi, G4, dG4 = christoffel(g, dg, ddg)
    R4 = ricci_scalar(gi, G4, dG4)
    sg = np.sqrt(-np.linalg.det(g))
    dA = np.array([a.d for a in A])          # dA[n, m] = d_m A_n
    F = dA.T - dA                            # F[m, n] = d_m A_n - d_n A_m
    F2 = np.einsum("ma,nb,mn,ab->", gi, gi, F, F)
    kin = gi @ phi.d @ phi.d
    box = np.einsum("mn,mn->", gi, phi.dd) - np.einsum("mn,lmn,l->", gi, G4, phi.d)
    known = sg * (R4 - 0.5 * kin - 0.25 * np.exp(-6 * AL * phi.v) * F2)
    return lhs - known, sg * box, abs(lhs)


def run():
    rows = np.array([one_point() for _ in range(8)])
    c = np.linalg.lstsq(rows[:, 1:2], rows[:, 0], rcond=None)[0][0]
    resid = np.max(np.abs(rows[:, 0] - c * rows[:, 1]))
    ok = resid < 1e-10 * (1 + rows[:, 2].max()) and abs(c + 2 * AL) < 1e-10
    return [("KK Ricci reduction: sqrt(-g5)R5 = sqrt(-g)[R - (dphi)^2/2 - e^{-sqrt3 phi}F^2/4 - 2a Box(phi)]",
             bool(ok), f"max residual {resid:.2e}, fitted Box coefficient {c:.9f} vs -2a = {-2*AL:.9f}")]


if __name__ == "__main__":
    for name, ok, det in run():
        print(("PASS" if ok else "FAIL"), name, "|", det)
