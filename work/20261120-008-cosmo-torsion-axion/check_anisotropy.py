"""Checks for the anisotropy estimate (PASS/FAIL)."""
import numpy as np
import sympy as sp

from anisotropy import S_L, transfer, TAU_STAR, D


def report(name, ok, detail):
    print(f"{'PASS' if ok else 'FAIL'}: {name} -- {detail}")


# 1. sympy: A_CB = (g/2)^2 (H/2pi)^2 R^2 with g = 2 beta/(F phi_i) -> (beta R/F)^2 (H/(2 pi phi_i))^2
g, H, R, b, F, p = sp.symbols("g H R beta F phi_i", positive=True)
acb = (g / 2) ** 2 * (H / (2 * sp.pi)) ** 2 * R**2
report("A_CB in terms of beta", sp.simplify(acb.subs(g, 2 * b / (F * p)) - (b * R / F) ** 2 * (H / (2 * sp.pi * p)) ** 2) == 0, "sympy")

# 2. scale-invariant projection: 2 L(L+1) int dx/x j_L(x)^2 = 1 (T = 1)
s = S_L([2, 10, 50])
report("S_L = 1 for T(k) = 1 (L = 2, 10, 50)", np.all(np.abs(s - 1) < 1e-2), f"{s}")

# 3. superhorizon modes frozen, subhorizon decay
t_small = transfer(1e-3 / TAU_STAR)
t_big = transfer(200 / TAU_STAR)
report("T(k tau* = 1e-3) = 1", abs(t_small - 1) < 1e-4, f"{t_small:.6f}")
report("|T(k tau* = 200)| < 0.02 (decays inside horizon)", abs(t_big) < 0.02, f"{t_big:.2e}")
