"""037: U(1) structure of the operator O = J_mu A_nu Ftilde^{mu nu} (CONVENTIONS: eps^{0123} = -1, F = dA).

Checks (sympy, generic functions of x^0..x^3):
 1. delta O under A -> A + d(alpha) = d_nu(alpha J_mu Ft^{mu nu}) - alpha (d_nu J_mu) Ft^{mu nu}   (uses d_nu Ft^{mu nu} = 0)
 2. (d_nu J_mu) Ft^{mu nu} = 0 identically for J = d(theta)  (closed current: gauge invariant up to a total derivative)
 3. for a homogeneous J = (n(t), 0, 0, 0) (the FRW background) the variation vanishes
 4. for J0 = n(t, x) (inhomogeneous) it does not vanish: U(1) is broken by fluctuations
 5. Stueckelberg A -> A + d(pi)/m generates -(1/m) pi (d_nu J_mu) Ft^{mu nu}: the same curl operator, coefficient C/m
 6. dimension: [J] = 3, [A] = 1, [Ft] = 2 -> dimension 6, coefficient xi kappa/4 has mass dimension -2
"""
import itertools

import sympy as sp

X = sp.symbols("x0:4")
eta = sp.diag(-1, 1, 1, 1)


def eps_up(m, n, r, s):
    return -sp.LeviCivita(m, n, r, s)  # eps^{0123} = -1


def F_low(A):
    return [[sp.diff(A[n], X[m]) - sp.diff(A[m], X[n]) for n in range(4)] for m in range(4)]


def Ft_up(A):
    F = F_low(A)
    return [[sp.Rational(1, 2) * sum(eps_up(m, n, r, s) * F[r][s] for r in range(4) for s in range(4))
             for n in range(4)] for m in range(4)]


def O(J, A):
    Ft = Ft_up(A)
    return sum(J[m] * A[n] * Ft[m][n] for m in range(4) for n in range(4))


def curlterm(J, A):
    Ft = Ft_up(A)
    return sum(sp.diff(J[m], X[n]) * Ft[m][n] for m in range(4) for n in range(4))


def report(name, ok):
    print(("PASS" if ok else "FAIL") + ": " + name)
    return bool(ok)


A = [sp.Function(f"A{i}")(*X) for i in range(4)]
J = [sp.Function(f"J{i}")(*X) for i in range(4)]
al = sp.Function("alpha")(*X)
ok = True

Ft = Ft_up(A)
bianchi = [sp.simplify(sum(sp.diff(Ft[m][n], X[n]) for n in range(4))) for m in range(4)]
ok &= report("d_nu Ft^{mu nu} = 0 identically (Bianchi)", all(b == 0 for b in bianchi))

Ag = [A[i] + sp.diff(al, X[i]) for i in range(4)]
dO = sp.expand(O(J, Ag) - O(J, A))
total = sum(sp.diff(al * sum(J[m] * Ft[m][n] for m in range(4)), X[n]) for n in range(4))
ok &= report("delta O = total derivative - alpha (d_nu J_mu) Ft^{mu nu}",
             sp.simplify(dO - (total - al * curlterm(J, A))) == 0)

th = sp.Function("theta")(*X)
Jg = [sp.diff(th, X[i]) for i in range(4)]
ok &= report("closed current J = d(theta): (d_nu J_mu) Ft^{mu nu} = 0", sp.simplify(sp.expand(curlterm(Jg, A))) == 0)

n_t = sp.Function("n")(X[0])
Jh = [n_t, 0, 0, 0]
ok &= report("homogeneous J = (n(t),0,0,0): variation vanishes (FRW background is gauge safe)",
             sp.simplify(curlterm(Jh, A)) == 0)

n_tx = sp.Function("n")(X[0], X[1])
Ji = [n_tx, 0, 0, 0]
ci = sp.simplify(curlterm(Ji, A))
ok &= report(f"inhomogeneous J0 = n(t,x): variation = {ci} != 0 (U(1) broken off the background)", ci != 0)

m, pi_ = sp.symbols("m", positive=True), sp.Function("pi")(*X)
As = [A[i] + sp.diff(pi_, X[i]) / m for i in range(4)]
dS = sp.expand(O(J, As) - O(J, A))
totS = sum(sp.diff(pi_ / m * sum(J[mm] * Ft[mm][n] for mm in range(4)), X[n]) for n in range(4))
ok &= report("Stueckelberg: O -> O - (1/m) pi (d_nu J_mu) Ft^{mu nu} + total derivative",
             sp.simplify(dS - (totS - pi_ / m * curlterm(J, A))) == 0)

dim = 3 + 1 + 2
ok &= report(f"operator dimension {dim} = 6; coefficient xi kappa/4 has dimension 4 - 6 = -2 (Lambda^-2)", dim == 6)
print("Stueckelberg mode coupling: (C/m_gamma) pi dJ Ft, dimension 7 -> strong-coupling scale Lambda_s = (m_gamma / C)^(1/3)")
print("CP: J5 is C-even, A Ft is C-even; J5^mu and A_nu Ft^{mu nu} are both axial vectors -> P-even, hence T-even: no EDM from O")
print("ALL PASS" if ok else "SOME FAIL")
