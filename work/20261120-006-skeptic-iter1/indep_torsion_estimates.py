"""Skeptic, independent checks of the torsion algebra (003) and the size estimates (003, 004).

No code imported from work/20261120-003 or -004. numpy + sympy only.

 A. Contortion K_{lmn} = 1/2 (T_{lmn} + T_{nlm} - T_{mnl}) (Gamma^l_{mn}, derivative index m):
    reproduces T and is metric compatible (K_{lmn} = -K_{nml}) for a random constant metric;
    the iter-1 CONVENTIONS.md formula fails metric compatibility unless T is totally antisymmetric.
 B. Riemann convention sanity on the unit 2-sphere (R = +2), then R(Gamma = K) for constant
    torsion in flat space: coefficient of T_m T^m and S_m S^m (claim: -2/3 and +3/2).
 C. Algebra: (1/2k)(3/2) S^2 + (3/4) S.J5 -> S = -(k/2) J5, L_eff = -(3k/16) J5^2.
 D. FRW: T^i_{jk} = f [ijk] on g = a^2 eta gives S^tau = -f/a^2 (eps_{0123} = +sqrt(-g)).
 E. Size estimates: needed xi S^0hat, minimal-EC S from a relic neutrino asymmetry, shortfall.
 F. KK-02 normalization: g = 6 kappa_CS g4^2 / f_a = 24 pi kappa_CS alpha / f_a from the
    reduced CS term and the theta kinetic term (by hand in 004, scripted here).
"""
import csv
import itertools
import math
import os
import sys

import numpy as np
import sympy as sp
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
res = []


def rep(name, ok, info=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + ": " + name + (("  | " + info) if info else ""))


def levi(n):
    e = np.zeros((n,) * n)
    for p in itertools.permutations(range(n)):
        e[p] = np.linalg.det(np.eye(n)[list(p)])
    return e


E4 = levi(4)


def part_A(rng):
    h = rng.normal(scale=0.2, size=(4, 4)); g = np.diag([-1., 1, 1, 1]) + (h + h.T) / 2
    Tu = rng.normal(size=(4, 4, 4)); Tu = Tu - Tu.transpose(0, 2, 1)       # T^l_{mn}
    Tl = np.einsum("lr,rmn->lmn", g, Tu)
    K = 0.5 * (Tl + Tl.transpose(1, 2, 0) - Tl.transpose(2, 0, 1))
    # numpy: X.transpose(1,2,0)[l,m,n] = X[n,l,m];  X.transpose(2,0,1)[l,m,n] = X[m,n,l]
    rep("A. K reproduces torsion: K_{lmn} - K_{lnm} = T_{lmn}", np.allclose(K - K.transpose(0, 2, 1), Tl))
    rep("A. K metric compatible: K_{lmn} = -K_{nml}", np.allclose(K, -K.transpose(2, 1, 0)))
    Kold = 0.5 * (Tl - Tl.transpose(1, 0, 2) - Tl.transpose(1, 2, 0))       # 1/2(T_lmn - T_mln - T_nlm)
    rep("A. old CONVENTIONS.md formula NOT metric compatible for generic T (control)",
        not np.allclose(Kold, -Kold.transpose(2, 1, 0)))
    Ta = Tl - Tl.transpose(1, 0, 2) + Tl.transpose(1, 2, 0) - Tl.transpose(0, 2, 1) + Tl.transpose(2, 0, 1) - Tl.transpose(2, 1, 0)
    Ta /= 6  # totally antisymmetric part
    Ko = 0.5 * (Ta - Ta.transpose(1, 0, 2) - Ta.transpose(1, 2, 0))
    Kn = 0.5 * (Ta + Ta.transpose(1, 2, 0) - Ta.transpose(2, 0, 1))
    rep("A. old and new agree for totally antisymmetric torsion (K = T/2)", np.allclose(Ko, Kn) and np.allclose(Kn, Ta / 2))


def ricci_scalar_const(g, Gam):
    """Constant Gamma^r_{m s} (derivative index m) in a constant metric.
    R^r_{s m n} = G^r_{m l} G^l_{n s} - G^r_{n l} G^l_{m s};  R_{sn} = R^r_{s r n}."""
    Rm = np.einsum("rml,lns->rsmn", Gam, Gam) - np.einsum("rnl,lms->rsmn", Gam, Gam)
    Ric = np.einsum("rsrn->sn", Rm)
    return np.einsum("sn,sn->", np.linalg.inv(g), Ric)


def part_B(rng):
    th, ph = sp.symbols("theta phi")
    x = [th, ph]; g = sp.diag(1, sp.sin(th)**2); gi = g.inv()
    G = [[[sum(gi[r, l] * (sp.diff(g[l, m], x[s]) + sp.diff(g[l, s], x[m]) - sp.diff(g[m, s], x[l])) / 2
               for l in range(2)) for s in range(2)] for m in range(2)] for r in range(2)]
    Rm = lambda r, s, m, n: (sp.diff(G[r][n][s], x[m]) - sp.diff(G[r][m][s], x[n])
                             + sum(G[r][m][l] * G[l][n][s] - G[r][n][l] * G[l][m][s] for l in range(2)))
    R = sp.simplify(sum(gi[s, n] * Rm(r, s, r, n) for r in range(2) for s in range(2) for n in range(2)))
    rep("B. Riemann convention gives R = +2 on the unit 2-sphere", R == 2)
    eta = np.diag([-1., 1, 1, 1]); ei = eta
    # pure axial: T_{lmn} = eps_{lmnk} S^k, eps_{0123} = +1 flat; check S = (1/6) eps^{mnrs} T_{nrs}
    S = rng.normal(size=4)
    Tl = np.einsum("lmnk,k->lmn", E4, S)
    eps_up = -E4
    S_back = np.einsum("mnrs,nrs->m", eps_up, Tl) / 6  # eps^{mnrs} T_{nrs}: all indices on T lower
    rep("B. T_{lmn} = eps_{lmnk} S^k inverts S^m = (1/6) eps^{mnrs} T_{nrs}", np.allclose(S_back, S))
    K = 0.5 * (Tl + Tl.transpose(1, 2, 0) - Tl.transpose(2, 0, 1))
    Gam = np.einsum("rl,lmn->rmn", ei, K)
    Rs = ricci_scalar_const(eta, Gam)
    S2 = S @ eta @ S
    rep("B. R(K) = +(3/2) S_m S^m for pure axial constant torsion", np.isclose(Rs, 1.5 * S2), f"R/S^2 = {Rs / S2:.6f}")
    # pure trace: T^l_{mn} = (1/3)(d^l_n T_m - d^l_m T_n)
    Tv = rng.normal(size=4)
    Tu = (np.einsum("ln,m->lmn", np.eye(4), Tv) - np.einsum("lm,n->lmn", np.eye(4), Tv)) / 3
    rep("B. trace convention T_m = T^l_{ml} recovered", np.allclose(np.einsum("lml->m", Tu), Tv))
    Tl = np.einsum("lr,rmn->lmn", eta, Tu)
    K = 0.5 * (Tl + Tl.transpose(1, 2, 0) - Tl.transpose(2, 0, 1))
    Rt = ricci_scalar_const(eta, np.einsum("rl,lmn->rmn", ei, K))
    T2 = Tv @ eta @ Tv
    rep("B. R(K) = -(2/3) T_m T^m for pure trace constant torsion", np.isclose(Rt, -2 / 3 * T2), f"R/T^2 = {Rt / T2:.6f}")


def part_C():
    k, S, J = sp.symbols("kappa S J")
    L = 1 / (2 * k) * sp.Rational(3, 2) * S**2 + sp.Rational(3, 4) * S * J
    sol = sp.solve(sp.diff(L, S), S)[0]
    rep("C. S = -(kappa/2) J5", sp.simplify(sol + k * J / 2) == 0)
    rep("C. L_eff = -(3 kappa/16) J5^2", sp.simplify(L.subs(S, sol) + 3 * k * J**2 / 16) == 0)


def part_D():
    a, f = sp.symbols("a f", positive=True)
    g = sp.diag(-a**2, a**2, a**2, a**2)
    sq = sp.sqrt(-g.det())
    T_l = lambda l, m, n: (a**2 * f * sp.LeviCivita(l, m, n) if min(l, m, n) >= 1 else 0)  # T_{ijk} = a^2 f [ijk]
    S0 = sp.Rational(1, 6) * sum((-sp.LeviCivita(0, n, r, s) / sq) * T_l(n, r, s)
                                 for n in range(1, 4) for r in range(1, 4) for s in range(1, 4))
    rep("D. FRW axial torsion: S^tau = -f/a^2", sp.simplify(S0 + f / a**2) == 0)


def load_inputs():
    d = {}
    with open(os.path.join(ROOT, "data", "bounds_inputs.csv"), newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            d[r["id"]] = float(r["value"])
    with open(os.path.join(ROOT, "data", "birefringence_measurements.csv"), newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r["id"] == "DPK25-beta-baseline":
                d["beta_deg"] = float(r["value"])
    with open(os.path.join(HERE, "inputs_skeptic.csv"), newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            d[r["id"]] = float(r["value"])
    return d


def part_E():
    c = load_inputs()
    beta = math.radians(c["beta_deg"])
    H0 = c["H0"] / (c["Mpc"] / 1e3) * c["hbar"]           # km/s/Mpc -> s^-1 -> eV
    MP = c["MP_reduced"]; kap = 1 / MP**2
    Tg = c["kB"] * c["T_cmb"]
    Om = c["Omega_m"]
    Or = (math.pi**2 / 15 * Tg**4 * (1 + 7 / 8 * (4 / 11)**(4 / 3) * c["N_eff"])) / (3 * H0**2 * MP**2)
    Ol = 1 - Om - Or
    H = lambda x: H0 * math.sqrt(Or / x**4 + Om / x**3 + Ol)
    aL = 1 / (1 + c["z_star"])
    t = quad(lambda x: 1 / (x * H(x)), aL, 1, limit=200)[0]
    I3 = quad(lambda x: x**-3 / (x * H(x)), aL, 1, limit=200)[0]
    need_const = 2 * beta / t
    need_a3 = 2 * beta / I3
    rep("E. needed |xi S^0hat| (constant) ~ 1.1e-35 eV [003]", 1.0e-35 < need_const < 1.25e-35, f"{need_const:.3e} eV = {need_const / H0:.2e} H0")
    rep("E. needed |xi S^0hat(today)| for a^-3 scaling ~ 2.8e-40 eV [003]", 2.5e-40 < need_a3 < 3.1e-40, f"{need_a3:.3e} eV")
    Tnu = (4 / 11)**(1 / 3) * Tg
    out = {}
    for xi in (c["xi_nu_generous"], c["xi_nu_Neff_extreme"]):
        n5 = 3 * Tnu**3 / (6 * math.pi**2) * (math.pi**2 * xi + xi**3)
        S = kap / 2 * n5
        out[xi] = (n5, S, math.log10(need_a3 / S))
        cm3 = (c["hbarc"] * 1e2)**3
        print(f"INFO: xi_nu = {xi}: n5 = {n5 / cm3:.1f} cm^-3, |S^0hat| = {S:.2e} eV, shortfall 10^{math.log10(need_a3 / S):.1f}")
    rep("E. minimal-EC neutrino source ~2e-68 eV, shortfall ~28 orders at xi_nu = 0.1 [003]",
        1.5e-68 < out[c["xi_nu_generous"]][1] < 2.5e-68 and 27.5 < out[c["xi_nu_generous"]][2] < 28.5)
    rep("E. shortfall stays > 26 orders even for an extreme |xi_nu| = 1", out[c["xi_nu_Neff_extreme"]][2] > 26)
    # effective operator scale: (xi kappa/2) J5 . A Ft  ==  J5 . A Ft / M^2 with M = MP sqrt(2/xi)
    xi_need = need_a3 / out[c["xi_nu_generous"]][1]
    print(f"INFO: equivalent dimension-6 scale M = MP*sqrt(2/xi_need) = {MP * math.sqrt(2 / xi_need) / 1e9:.2e} GeV "
          "(compare neutrino-current CS literature, e.g. 2302.11140)")


def part_F():
    kcs, g4, L, R, alpha = sp.symbols("kappa_CS g_4 L R alpha", positive=True)
    # reduced term 2 pi R kappa_CS * 3 theta f^f with theta = vartheta/(2 pi R), f = g4 f_c,
    # f^f = -(1/2) sqrt(-g) f f~  ->  coefficient of vartheta f_c f~_c:
    coef = 2 * sp.pi * R * kcs * 3 / (2 * sp.pi * R) * g4**2 * (-sp.Rational(1, 2))
    # theta kinetic: -(1/(2 g5^2)) * (2 pi R) * (L/R) * (R/L)^2 (d theta)^2, g5^2 = 2 pi L g4^2
    g5sq = 2 * sp.pi * L * g4**2
    kin = sp.Rational(1, 2) / g5sq * 2 * sp.pi * R * (L / R) * (R / L)**2 / (2 * sp.pi * R)**2  # coefficient of (d vartheta)^2 /(-1)
    fa = sp.sqrt(2 * kin)
    rep("F. f_a = 1/(2 pi L g4)", sp.simplify(fa - 1 / (2 * sp.pi * L * g4)) == 0)
    g = sp.simplify(-4 * coef / fa)  # -(g/4) phi F F~ with vartheta = phi/f_a
    rep("F. g_phigamma = 6 kappa_CS g4^2 / f_a = 24 pi kappa_CS alpha / f_a",
        sp.simplify(g - 6 * kcs * g4**2 * (2 * sp.pi * L * g4)) == 0
        and sp.simplify((6 * kcs * g4**2).subs(g4, sp.sqrt(4 * sp.pi * alpha)) - 24 * sp.pi * kcs * alpha) == 0)


def main():
    rng = np.random.default_rng(11)
    part_A(rng); part_B(rng); part_C(); part_D(); part_E(); part_F()
    return all(res)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
