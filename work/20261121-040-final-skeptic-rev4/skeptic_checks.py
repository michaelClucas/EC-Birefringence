"""Final skeptic 040: independent checks of three weak kills and two cost numbers.

T1  (div S) F Ftilde / Lambda^2 on EC spin-sourced torsion (skeptic 020 estimate): in FRW with S^mu = -(kappa/2) J5^mu
    and a conserved axial current, div S = 0 identically, so the operator gives beta = 0 (stronger than 020's cutoff).
T2  EC-03 photon decoupling: B1 = 1 - tau f(tau)^2 expanded from the Feynman-parameter form (not the arcsin form used in
    012), and the /96 normalisation (per-electron 1/12 over sum N_c Q^2 = 8).
T3  Orbifold / warped KK: g_{mu 5} and A_5 are Z2-odd under y -> -y, so their y-independent zero modes vanish; the RS warp
    factor exp(-2k|y|) is even. Also: with a torsionful connection, nabla_[M A_N] = d A - T.A, so F = dA (gauge invariant)
    is torsion independent in any dimension.
T4  Cassini shortfall convention: KK-04 (4.34, 1 sigma) vs EC-06d heavy (3.87, |central| + 2 sigma), same data row.
T5  EC-04 tuning: (M_P/H0)^2.
Inputs from data/. Run from the repo root.
"""
import csv
import math
import os

import sympy as sp

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
NP = NF = 0


def check(name, ok, msg=""):
    global NP, NF
    NP += bool(ok)
    NF += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {msg}")


def rows(name):
    with open(os.path.join(ROOT, "data", name), newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


# T1 ------------------------------------------------------------------------------------------
t = sp.symbols("t", positive=True)
a = sp.Function("a")(t)
Q5 = sp.symbols("Q5")                      # comoving axial charge, conserved for massless fermions
kap = sp.symbols("kappa", positive=True)
J0 = Q5 / a**3                             # J5^0 in cosmic time, homogeneous
S0 = -kap / 2 * J0
divS = sp.diff(a**3 * S0, t) / a**3        # nabla_mu S^mu = a^-3 d(a^3 S^0)/dt on flat FRW
check("T1 div S = 0 for conserved J5 on FRW", sp.simplify(divS) == 0)
# with explicit non-conservation (mass term / anomaly), div S = -(kappa/2) d_t J^0 + 3H J^0 = -(kappa/2) * source
J0b = sp.Function("J")(t)
divS_b = sp.diff(a**3 * (-kap / 2 * J0b), t) / a**3
nonc = sp.diff(J0b, t) + 3 * sp.diff(a, t) / a * J0b     # nabla_mu J5^mu on FRW
check("T1b div S = -(kappa/2) x nabla.J5 (axial non-conservation)", sp.simplify(divS_b + kap / 2 * nonc) == 0)
print("  T1: the operator is sourced only by axial non-conservation (2 i m psibar g5 psi, anomaly), each further"
      " suppressed; 020's (H/Lambda)^2 estimate is an upper bound, so the kill holds a fortiori.")

# T2 ------------------------------------------------------------------------------------------
x, y, z = sp.symbols("x y z", positive=True)   # z = q^2/m^2
integrand_series = sp.series(1 / (1 - x * y * z), z, 0, 3).removeO()
tauf2 = 2 * sp.integrate(sp.integrate(integrand_series, (y, 0, 1 - x)), (x, 0, 1))
B1 = sp.expand(1 - tauf2)
check("T2 B1 = -z/12 - z^2/90 (Feynman-parameter expansion)", sp.simplify(B1 + z / 12 + z**2 / 90) == 0, str(B1))
qsum = 3 * 1 + 3 * 3 * (sp.Rational(4, 9) + sp.Rational(1, 9))
check("T2 sum N_c Q^2 = 8 and electron share (1/12)/8 = 1/96", qsum == 8 and sp.Rational(1, 12) / qsum == sp.Rational(1, 96))

# T3 ------------------------------------------------------------------------------------------
X = sp.symbols("x0:4")
yy = sp.symbols("y", real=True)
kk = sp.symbols("k", positive=True)
# tensor transformation of g_{mu 5} under y -> -y: Jacobian diag(1,1,1,1,-1)
Jac = sp.diag(1, 1, 1, 1, -1)
G = sp.Matrix(5, 5, lambda i, j: sp.Symbol(f"g{min(i, j)}{max(i, j)}"))
Gp = Jac.T * G * Jac
check("T3 g_{mu5} -> -g_{mu5}, g_{mu nu}, g_55 even under y -> -y",
      all(sp.simplify(Gp[m, 4] + G[m, 4]) == 0 for m in range(4)) and Gp[4, 4] == G[4, 4] and Gp[0, 1] == G[0, 1])
# odd field with y-independent zero mode: h(-y) = -h(y) and h const => h = 0
h0 = sp.symbols("h0")
check("T3 Z2-odd zero mode vanishes", sp.solve(sp.Eq(h0, -h0), h0) == [0])
warp = sp.exp(-2 * kk * sp.Abs(yy))
check("T3 RS warp factor even under y -> -y", sp.simplify(warp.subs(yy, -yy) - warp) == 0)
# F with torsionful connection
n = 5
A = [sp.Function(f"A{i}")(*X, yy) for i in range(n)]
coords = list(X) + [yy]
Gam = [[[sp.Symbol(f"G{l}_{m}{nn}") for nn in range(n)] for m in range(n)] for l in range(n)]
ok = True
for m in range(n):
    for nn in range(n):
        cov = (sp.diff(A[nn], coords[m]) - sum(Gam[l][m][nn] * A[l] for l in range(n))) - \
              (sp.diff(A[m], coords[nn]) - sum(Gam[l][nn][m] * A[l] for l in range(n)))
        dA = sp.diff(A[nn], coords[m]) - sp.diff(A[m], coords[nn])
        Tor = [Gam[l][m][nn] - Gam[l][nn][m] for l in range(n)]
        ok &= sp.simplify(cov - (dA - sum(Tor[l] * A[l] for l in range(n)))) == 0
check("T3 nabla_M A_N - nabla_N A_M = (dA)_MN - T^L_MN A_L in 5D (so gauge-invariant F = dA is torsion-free)", ok)

# T4 ------------------------------------------------------------------------------------------
kb = rows("kk_bounds_inputs.csv")["gamma_minus_1_Cassini"]
c, s = float(kb["value"]), float(kb["err"])
o1 = math.log10(0.5 / s)
o2 = math.log10(0.5 / (abs(c) + 2 * s))
print(f"  T4 Cassini: 0.5/sigma -> {o1:.2f} orders (KK-04 headline); 0.5/(|c|+2 sigma) -> {o2:.2f} (EC-06d heavy)")
check("T4 KK-04 4.34 and EC-06d 3.87 reproduced (different conventions for the same bound)",
      abs(o1 - 4.34) < 0.01 and abs(o2 - 3.87) < 0.01)

# T5 ------------------------------------------------------------------------------------------
H0 = float(rows("cosmo_inputs.csv")["H0"]["value"])
hbar = 6.582119569e-16
H0_eV = H0 / 3.0856775814913673e19 * hbar      # km/s/Mpc -> s^-1 -> eV
MP = 2.435e27                                  # reduced Planck mass in eV (CONVENTIONS.md)
r = (MP / H0_eV) ** 2
print(f"  T5 H0 = {H0_eV:.3e} eV, (M_P/H0)^2 = {r:.2e}")
check("T5 EC-04 tuning ~ (M_P/H0)^2 within a factor 2 of the quoted 4e120", 0.5 < 4e120 / r < 2.0, f"ratio {4e120 / r:.2f}")

print(f"== {NP} PASS, {NF} FAIL ==")
