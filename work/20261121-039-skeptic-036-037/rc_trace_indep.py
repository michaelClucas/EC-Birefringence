"""Skeptic 039, item A1: independent check of 036's Riemann-Cartan scalar for pure-trace torsion and the 0+ mass.

Independent route (does not use 036's contortion formula or its Riemann index formula):
  1. Flat metric. Solve the 64 linear equations  Gamma_{b m a} + Gamma_{a m b} = 0  (metric compatibility,
     Gamma_{b m a} = g_{br} Gamma^r_{ma}) and  Gamma^l_{mn} - Gamma^l_{nm} = T^l_{mn}  for Gamma.
  2. Get the curvature from the commutator on a generic vector field V:
        [nabla_m, nabla_n] V^r = R^r_{s m n} V^s - T^l_{mn} nabla_l V^r
     computed by explicit second covariant derivatives (nabla_m acting on the tensor nabla_n V^r).
  3. R = g^{sn} R^m_{s m n} and compare with 036: R = R_LC + 2 d.T - (2/3) T.T.
  4. Mass: Einstein-frame algebra for L = 1/2 [F R - K (d chi)^2] - V with 036's F, K, V and general a0,
     compared with m^2 = -2 a0 mu/(a1 b), mu = a1 - a0; and the Brans-Dicke alpha = mu/(3 a1).
"""
import itertools

import sympy as sp

npass = nfail = 0


def check(name, ok):
    global npass, nfail
    npass += bool(ok)
    nfail += (not ok)
    print(("PASS" if ok else "FAIL") + ": " + name)


X = sp.symbols("t x y z", real=True)
R4 = range(4)
g = sp.diag(-1, 1, 1, 1)
gi = g.inv()
Tv = [sp.Function(f"T{i}")(*X) for i in R4]
d = sp.eye(4)
Tup = [[[sp.Rational(1, 3) * (d[l, n] * Tv[m] - d[l, m] * Tv[n]) for n in R4] for m in R4] for l in R4]

# 1. solve for the connection
G = {(l, m, n): sp.Symbol(f"G{l}{m}{n}") for l, m, n in itertools.product(R4, R4, R4)}
eqs = []
for b, m, a in itertools.product(R4, R4, R4):
    eqs.append(sum(g[b, r] * G[(r, m, a)] for r in R4) + sum(g[a, r] * G[(r, m, b)] for r in R4))
for l, m, n in itertools.product(R4, R4, R4):
    if m < n:
        eqs.append(G[(l, m, n)] - G[(l, n, m)] - Tup[l][m][n])
sol = sp.solve(eqs, list(G.values()), dict=True)
check("connection uniquely determined by metric compatibility + torsion", len(sol) == 1 and len(sol[0]) == 64)
Gam = [[[sol[0][G[(l, m, n)]] for n in R4] for m in R4] for l in R4]

# 2. commutator on a generic vector
V = [sp.Function(f"V{i}")(*X) for i in R4]
DV = [[sp.diff(V[r], X[n]) + sum(Gam[r][n][s] * V[s] for s in R4) for r in R4] for n in R4]  # DV[n][r] = nabla_n V^r


def DDV(m, n, r):  # nabla_m (nabla_n V^r), the tensor W_n^r = DV[n][r]
    return (sp.diff(DV[n][r], X[m]) + sum(Gam[r][m][s] * DV[n][s] for s in R4)
            - sum(Gam[s][m][n] * DV[s][r] for s in R4))


Rsc = 0
for s, n in itertools.product(R4, R4):
    if gi[s, n] == 0:
        continue
    m = None
    # Ricci R_{sn} = R^m_{s m n}: extract R^r_{s m n} as the coefficient of V^s in comm + T nabla V
    for mm in R4:
        comm = sp.expand(DDV(mm, n, mm) - DDV(n, mm, mm) + sum(Tup[l][mm][n] * DV[l][mm] for l in R4))
        Rsc += gi[s, n] * comm.coeff(V[s])
Rsc = sp.simplify(Rsc)
divT = sum(gi[m, m] * sp.diff(Tv[m], X[m]) for m in R4)
TT = sum(gi[m, m] * Tv[m] ** 2 for m in R4)
print("   R (flat, pure trace) =", sp.simplify(Rsc))
check("commutator check: no derivative-of-V terms left (torsion term in the commutator is right)",
      all(sp.expand(DDV(0, 1, 2) - DDV(1, 0, 2) + sum(Tup[l][0][1] * DV[l][2] for l in R4)).coeff(sp.diff(V[k], X[j])) == 0
          for k in R4 for j in R4))
check("R = 2 d.T - (2/3) T.T  (036 A3)", sp.simplify(Rsc - (2 * divT - sp.Rational(2, 3) * TT)) == 0)
check("Shapiro-normalisation alternative -(4/3) T.T is NOT what our conventions give",
      sp.simplify(Rsc - (2 * divT - sp.Rational(4, 3) * TT)) != 0)

# 4. mass from the reduced scalar-tensor theory with general a0 (a0 = -M_P^2)
a0, a1, b, chi = sp.symbols("a0 a1 b chi", real=True)
mu = a1 - a0
Y = chi + 6 * mu / b
F = -a0 + b * chi / 6
K = -(b / 4) / Y
Vp = b / 24 * chi**2
# Einstein frame: g_E = (F/F0) g; canonical kinetic coefficient Z (units F0 = M_P^2) and potential U
F0 = -a0
Z = K / F + sp.Rational(3, 2) * (sp.diff(F, chi) / F) ** 2
U = Vp * F0**2 / F**2
m2 = sp.simplify((sp.diff(U, chi, 2) / (Z * F0)).subs(chi, 0))
check(f"Einstein-frame mass m^2 = {sp.factor(m2)} equals -2 a0 mu/(a1 b)", sp.simplify(m2 - (-2 * a0 * mu / (a1 * b))) == 0)
omega = sp.simplify((K * F / sp.diff(F, chi) ** 2).subs(chi, 0))
alpha = sp.simplify(1 / (3 + 2 * omega))
check(f"Brans-Dicke omega = {omega}, alpha = 1/(3+2 omega) = mu/(3 a1)", sp.simplify(alpha - mu / (3 * a1)) == 0)
MP2 = sp.Symbol("M2", positive=True)
check("a1 -> infinity gives the Starobinsky scalaron 2 M_P^2/b",
      sp.simplify(sp.limit(m2.subs(a0, -MP2), a1, sp.oo) - 2 * MP2 / b) == 0)
print(f"== {npass} PASS, {nfail} FAIL ==")
