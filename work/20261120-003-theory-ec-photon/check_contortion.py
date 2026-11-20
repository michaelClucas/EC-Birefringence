"""Task 1: check the contortion formula of CONVENTIONS.md on random metric/torsion data at a point.

Gamma^l_{mn} = Christoffel + K^l_{mn}. Requirements:
  (a) Gamma^l_{mn} - Gamma^l_{nm} = T^l_{mn}
  (b) nabla_m g_{ab} = d_m g_{ab} - Gamma^r_{ma} g_{rb} - Gamma^r_{mb} g_{ar} = 0
"""
import itertools
import sys

import sympy as sp

from common import (R4, christoffel, contortion_conventions_md, contortion_correct, lower1,
                    raise1, random_metric_point, random_torsion, report, rng)


def test_formula(kfun, seed):
    r = rng(seed)
    g, dg = random_metric_point(r)
    T = random_torsion(r)
    Tl = lower1(g, T)
    K = raise1(g, kfun(Tl))
    C = christoffel(g, dg)
    G = [[[C[l][m][n] + K[l][m][n] for n in R4] for m in R4] for l in R4]
    tors_ok = all(sp.expand(G[l][m][n] - G[l][n][m] - T[l][m][n]) == 0
                  for l, m, n in itertools.product(R4, R4, R4))
    comp_ok = all(sp.expand(dg[m][a][b] - sum(G[rr][m][a] * g[rr, b] + G[rr][m][b] * g[a, rr]
                                                 for rr in R4)) == 0
                  for m, a, b in itertools.product(R4, R4, R4))
    return tors_ok, comp_ok, Tl


def which_permutation(seed):
    """Identify what the CONVENTIONS.md expression equals in terms of the correct K."""
    r = rng(seed)
    g, _ = random_metric_point(r)
    Tl = lower1(g, random_torsion(r))
    Kc = contortion_correct(Tl)
    Kg = contortion_conventions_md(Tl)
    found = []
    for p in itertools.permutations(range(3)):
        for sgn in (1, -1):
            ok = True
            for idx in itertools.product(R4, R4, R4):
                j = tuple(idx[p[i]] for i in range(3))
                if sp.expand(Kg[idx[0]][idx[1]][idx[2]] - sgn * Kc[j[0]][j[1]][j[2]]) != 0:
                    ok = False
                    break
            if ok:
                found.append((p, sgn))
    return found


def main():
    allok = True
    for seed in (1, 2):
        t, c, _ = test_formula(contortion_correct, seed)
        allok &= report(f"corrected K_lmn = 1/2(T_lmn + T_nlm - T_mnl): torsion reproduced (seed {seed})", t)
        allok &= report(f"corrected K_lmn: metric compatible (seed {seed})", c)
    t, c, _ = test_formula(contortion_conventions_md, 1)
    print(f"INFO: CONVENTIONS.md K_lmn = 1/2(T_lmn - T_mln - T_nlm): torsion reproduced={t}, "
          f"metric compatible={c}")
    allok &= report("CONVENTIONS.md contortion formula is inconsistent with its own definitions "
                    "(expected: fails)", not (t and c))
    found = which_permutation(3)
    names = "lmn"
    for p, sgn in found:
        print(f"INFO: CONVENTIONS.md expression == {'+' if sgn > 0 else '-'}"
              f"K_{{{''.join(names[i] for i in p)}}} (correct K, index order l,m,n of Gamma^l_mn)")
    return allok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
