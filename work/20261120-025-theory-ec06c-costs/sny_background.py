"""EC-06d: Shie-Nester-Yo (0805.3834) PGT scalar-torsion (0+) background, integrated back toward z* (task 025).

SNY equations (their Eqs. 43-47, matter only, scaled units 8 pi G = 1, a0 = -1, H(t0 = 1) = 1 = 1/T0):
    a' = aH,  H' = mu R/(6 a1) + Tr/(6 a1) - 2H^2,  Phi' = -a0 R/(2 a1) + Tr/(2 a1) - 3H Phi + Phi^2/3,
    R' = -(2/3)(R + 6mu/b) Phi,   constraint (b/18) Y (3H - Phi)^2 - (b/24) R^2 - 3 a1 H^2 = rho,
    Y = R + 6mu/b,  mu = a1 - a0,  Tr = -rho (dust), rho ~ a^-3.
Covariant connection equation (their Eq. 20): grad R = -(2/3) Y T_mu, so T_mu = -(3/2) d ln|Y| holds off FRW
(given P = Q = 0, their Eqs. 21-23): p = xi T is closed exactly, and beta = (xi/2) int Phi dt = -(3 xi/4) ln(Y0/Y_src).
Sign: SNY T_i = T_{ij}^j with T_ij^mu = 2 d_[i e_j]^mu + ..., derivative index first: same as our T_mu = T^l_{ml},
so our lower T_t = Phi.

Checks: period formula of their case II, rho(t0) = 0.83 of case III, constraint preserved along the
integration, int Phi dt = -(3/2) ln Y. Output: Phi/H and Y at z = 1, 10, 100, 1090 (or where the solution
breaks down), xi needed for beta = 0.215 deg.
"""
import csv
import math
import os

import numpy as np
from scipy.integrate import solve_ivp

from common import ROOT, check, load, summary

HERE = os.path.dirname(os.path.abspath(__file__))
inp = load("ec06c_costs_inputs.csv")
beta_rad = math.radians(float(load("birefringence_measurements.csv")["DPK25-beta-baseline"]["value"]))
zstar = float(load("cosmo_inputs.csv")["z_star"]["value"])
a0 = -1.0

# case II period check (their Sec. V.B: mu = 1.2, b = 4, period 3.63)
mu2, b2 = float(inp["SNY_caseII_mu"]["value"]), float(inp["SNY_caseII_b"]["value"])
a1_2 = mu2 + a0
per = 2 * math.pi * math.sqrt(-a1_2 * b2 / (2 * a0 * mu2))
check("Y1 case II oscillation period 2 pi sqrt(-a1 b/(2 a0 mu)) = 3.63 with a1 = mu + a0 (SNY Eq. 53)",
      abs(per - 3.63) < 0.01, f"({per:.3f})")


def run(case):
    mu = float(inp[f"SNY_case{case}_mu"]["value"])
    b = float(inp[f"SNY_case{case}_b"]["value"])
    Phi0 = float(inp[f"SNY_case{case}_Phi"]["value"])
    Y0 = float(inp[f"SNY_case{case}_Y"]["value"])
    a_now = float(inp["SNY_a_now"]["value"])
    a1 = mu + a0
    R0 = Y0 - 6 * mu / b
    H0 = 1.0
    rho0 = b / 18 * Y0 * (3 * H0 - Phi0)**2 - b / 24 * R0**2 - 3 * a1 * H0**2

    def rhs(t, s):
        a, H, Phi, R, I = s
        rho = rho0 * (a_now / a)**3
        Tr = -rho
        return [a * H, mu * R / (6 * a1) + Tr / (6 * a1) - 2 * H**2,
                -a0 * R / (2 * a1) + Tr / (2 * a1) - 3 * H * Phi + Phi**2 / 3,
                -2.0 / 3.0 * (R + 6 * mu / b) * Phi, Phi]

    def constraint(s):
        a, H, Phi, R, _ = s
        Y = R + 6 * mu / b
        return (b / 18 * Y * (3 * H - Phi)**2 - b / 24 * R**2 - 3 * a1 * H**2) / (rho0 * (a_now / a)**3) - 1

    def hit(z):
        f = lambda t, s: s[0] - a_now / (1 + z)
        f.terminal = False
        return f

    zs = [1.0, 10.0, 100.0, zstar]
    events = [hit(zz) for zz in zs]

    def blow(t, s):
        return 1e8 - abs(s[1]) - abs(s[2])
    blow.terminal = True
    sol = solve_ivp(rhs, (1.0, -5.0), [a_now, H0, Phi0, R0, 0.0], method="LSODA", rtol=1e-10, atol=1e-12,
                    events=events + [blow], dense_output=True)
    rows = []
    for zz, te, ye in zip(zs, sol.t_events[:4], sol.y_events[:4]):
        if len(te) == 0:
            rows.append(dict(case=case, z=zz, reached=0))
            continue
        s = ye[0]
        Y = s[3] + 6 * mu / b
        rows.append(dict(case=case, z=zz, reached=1, t=te[0], H=s[1], Phi=s[2], eps=s[2] / s[1], Y=Y,
                         intPhi=-s[4], cons=constraint(s), lnY=math.log(abs(Y / Y0)) if Y * Y0 > 0 else float("nan")))
    # earliest point reached
    sl = sol.y[:, -1]
    zmax = a_now / sl[0] - 1
    return rows, zmax, rho0, sol


out = []
for case in ("III", "IV", "V"):
    rows, zmax, rho0, sol = run(case)
    if case == "III":
        check("Y2 case III: rho(t0) from the constraint = 0.83 (SNY Sec. V.C)", abs(rho0 - 0.83) < 0.01, f"({rho0:.3f})")
    print(f"   case {case}: rho0 = {rho0:.3f}; backward integration reached z = {zmax:.3g} (stop: |H|+|Phi| > 1e8 or t = -5)")
    for r in rows:
        if r["reached"]:
            print(f"     z = {r['z']:7.1f}: t = {r['t']:.4f}, H = {r['H']:.4g}, Phi = {r['Phi']:.4g}, Phi/H = {r['eps']:.3f}, "
                  f"Y = {r['Y']:.4g}, constraint residual = {r['cons']:.1e}")
            check(f"Y3 case {case} z = {r['z']:.0f}: Friedmann constraint preserved (< 1e-5)", abs(r["cons"]) < 1e-5)
            check(f"Y4 case {case} z = {r['z']:.0f}: int_t^t0 Phi dt = -(3/2) ln(Y0/Y) (closed-form beta)",
                  abs(r["intPhi"] - (-1.5) * (-r["lnY"])) < 1e-6 * max(1, abs(r["intPhi"])),
                  f"({r['intPhi']:.6g} vs {1.5 * r['lnY']:.6g})")
        else:
            print(f"     z = {r['z']:7.1f}: not reached")
        out.append(r)
    out.append(dict(case=case, z=zmax, reached=-1))

# beta per unit xi: beta = (xi/2) int_{t_src}^{t0} Phi dt
print("\n   beta/xi = (1/2) int Phi dt from z to today, and xi for beta = 0.215 deg:")
res = []
for r in out:
    if r.get("reached") == 1:
        bx = 0.5 * r["intPhi"]
        r["beta_per_xi"] = bx
        r["xi_needed"] = beta_rad / bx if bx != 0 else float("inf")
        print(f"     case {r['case']} z = {r['z']:7.1f}: beta/xi = {bx:+.4f} rad, xi = {r['xi_needed']:+.3g}")
        res.append(r)
reached_star = [r for r in res if abs(r["z"] - zstar) < 1]
check("Y5 at least one SNY case reaches z* in its own (matter-only) equations", len(reached_star) > 0,
      f"({len(reached_star)} of 3)")
# expansion history vs flat LCDM (matter + Lambda, radiation omitted: including it raises H_LCDM(z*) by ~1.5)
Om = float(load("cosmo_inputs.csv")["Omega_m"]["value"])
for r in res:
    r["H_over_LCDM"] = r["H"] / math.sqrt(Om * (1 + r["z"])**3 + 1 - Om)
    print(f"     case {r['case']} z = {r['z']:7.1f}: H_SNY/H_LCDM = {r['H_over_LCDM']:.3g}")
star = [r["H_over_LCDM"] for r in reached_star]
check("Y6 all three SNY cases have H(z*) > 10 H_LCDM(z*) (torsion-dominated, radiation-like early universe: excluded by the CMB)",
      len(star) == 3 and min(star) > 10, f"(min {min(star):.3g})")
check("Y7 on these torsion-dominated histories Phi/H -> 3 at z* (same eps_T = 3 as flat f(T))",
      all(abs(r["eps"] - 3) < 0.01 for r in reached_star))
with open(os.path.join(HERE, "sny_background.csv"), "w", newline="", encoding="utf-8") as f:
    keys = ["case", "z", "reached", "t", "H", "Phi", "eps", "Y", "intPhi", "beta_per_xi", "xi_needed", "H_over_LCDM", "cons"]
    w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    for r in out:
        w.writerow(r)
print("   wrote sny_background.csv (reached = -1 rows give the earliest z the integration reached)")
summary()
