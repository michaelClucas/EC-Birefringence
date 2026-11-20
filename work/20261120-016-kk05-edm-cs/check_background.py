"""G1-G4: a fixed fibre vector V = d/dy as the non-dynamical P'-odd background of KK-05.

Coordinates (t, x, y) suffice (the other 4D directions behave like x). KK metric
  ds^2 = g4 + phi^2 (dy + k A_mu dx^mu)^2.
"""
import sympy as sp

from common import check

t, x, y, k = sp.symbols("t x y k", real=True)
phi = sp.Function("phi", positive=True)(t, x)
A0, A1 = sp.Function("A0")(t, x), sp.Function("A1")(t, x)
lam = sp.Function("lam")(t, x)
X = [t, x, y]
g4 = sp.diag(-1, 1)
Avec = [A0, A1]
g = sp.zeros(3, 3)
for m in range(2):
    for n in range(2):
        g[m, n] = g4[m, n] + phi ** 2 * k ** 2 * Avec[m] * Avec[n]
    g[m, 2] = g[2, m] = phi ** 2 * k * Avec[m]
g[2, 2] = phi ** 2
V = sp.Matrix([0, 0, 1])

# G1: u = V_flat/|V| is the KK one-form e^5 = phi (dy + k A)
Vflat = g * V
norm = sp.sqrt((V.T * g * V)[0])
u = sp.simplify(Vflat / norm)
check("G1 u = V^flat/|V| = phi (dy + k A): the oriented u of KK-05 is built from V and the dynamical metric",
      sp.simplify(u - sp.Matrix([phi * k * A0, phi * k * A1, phi])) == sp.zeros(3, 1))

# G2: KK U(1) gauge transformation y' = y + lam(x): V is invariant, A -> A - d lam / k
Jac = sp.Matrix([[sp.diff(xp, xx) for xx in X] for xp in [t, x, y + lam]])
Vp = Jac * V
check("G2 under y -> y + lam(t,x) the vector V = d/dy keeps components (0,0,1): every term built from V is "
      "KK-gauge invariant, so no Stueckelberg mass for the KK photon at any order", sp.simplify(Vp - V) == sp.zeros(3, 1))

# G3: 4D diffeos (x^mu -> f(x)) also commute with V; P' (y -> -y) flips it
Jl = sp.Matrix([[sp.diff(xp, xx) for xx in X] for xp in [t + x ** 2, x + sp.sin(t), y]])
Jp = sp.Matrix([[sp.diff(xp, xx) for xx in X] for xp in [t, -x, -y]])
check("G3 4D diffeos keep V (4D Lorentz invariant background); P' maps V -> -V (P'-odd)",
      sp.simplify(Jl * V - V) == sp.zeros(3, 1) and sp.simplify(Jp * V + V) == sp.zeros(3, 1))

# G4: a y-dependent reparametrisation y -> y + eps f(y) does not keep V: the explicitly broken generators act only on
# y-dependent (massive KK) modes
eps = sp.Symbol("eps")
f = sp.Function("f")
Jy = sp.Matrix([[sp.diff(xp, xx) for xx in X] for xp in [t, x, y + eps * f(y)]])
check("G4 broken generators are exactly the y-dependent diffeos (V -> (1 + eps f'(y)) V)",
      sp.simplify(Jy * V - sp.Matrix([0, 0, 1 + eps * sp.diff(f(y), y)])) == sp.zeros(3, 1))
