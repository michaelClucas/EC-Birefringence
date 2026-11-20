"""Task 3: homogeneous isotropic torsion on spatially flat FRW (conformal coordinates).

Impose Lie_xi T^l_{mn} = 0 for the three translations d_i and the three rotations
x^j d_k - x^k d_j. Translations make all components functions of tau only. For tau-only
components the rotation condition is linear algebra:
 Lie_xi T^l_mn = xi^r d_r T^l_mn - T^r_mn d_r xi^l + T^l_rn d_m xi^r + T^l_mr d_n xi^r
(first term zero because d_i T = 0). Solve for all 24 independent components.
Expected: 2 free functions, T^i_{0j} = -T^i_{j0} = h delta^i_j and T^i_{jk} = f [ijk].
Then compute S^m and T_m for g = a^2 eta. Parity: also report what an improper
reflection x -> -x does (f is a pseudoscalar).
"""
import itertools
import sys

import sympy as sp

from common import ETA, R4, axial_vector, report, trace_vector, zeros3, perm

X = sp.symbols("x0:4")


def lie(T, xi):
    out = zeros3()
    J = [[sp.diff(xi[a], X[b]) for b in R4] for a in R4]  # J[a][b] = d_b xi^a
    for l, m, n in itertools.product(R4, R4, R4):
        v = -sum(T[r][m][n] * J[l][r] for r in R4)
        v += sum(T[l][r][n] * J[r][m] for r in R4)
        v += sum(T[l][m][r] * J[r][n] for r in R4)
        out[l][m][n] = sp.expand(v)
    return out


def main():
    ok = True
    syms = {}
    T = zeros3()
    for l in R4:
        for m in R4:
            for n in range(m + 1, 4):
                s = sp.Symbol(f"t{l}{m}{n}")
                syms[s] = (l, m, n)
                T[l][m][n] = s
                T[l][n][m] = -s
    eqs = []
    for j, k in ((1, 2), (2, 3), (3, 1)):
        xi = [0, 0, 0, 0]
        xi[j], xi[k] = X[k], -X[j]
        L = lie(T, xi)
        eqs += [L[l][m][n] for l, m, n in itertools.product(R4, R4, R4)]
    eqs = [e for e in set(eqs) if e != 0]
    sol = sp.solve(eqs, list(syms), dict=True)[0]
    free = sorted({s for v in sol.values() for s in v.free_symbols} | (set(syms) - set(sol)),
                  key=str)
    ok &= report("rotation invariance leaves exactly 2 free components", len(free) == 2)
    Tsol = [[[sp.sympify(T[l][m][n]).subs(sol) for n in R4] for m in R4] for l in R4]
    print("INFO: free components:", free)
    nz = {(l, m, n): Tsol[l][m][n] for l, m, n in itertools.product(R4, R4, R4)
          if m < n and Tsol[l][m][n] != 0}
    print("INFO: nonzero T^l_{mn} (m<n):", nz)
    h, f = sp.symbols("h f")
    ref = zeros3()
    for i in (1, 2, 3):
        ref[i][0][i], ref[i][i][0] = h, -h
    for i, j, k in itertools.product((1, 2, 3), repeat=3):
        ref[i][j][k] = f * perm(i, j, k)
    # match: solution space must equal span{h-structure, f-structure}
    mp = {}
    for fs in free:
        Tone = [[[Tsol[l][m][n].subs({x: (1 if x == fs else 0) for x in free}) for n in R4]
                 for m in R4] for l in R4]
        mp[fs] = Tone
    def is_span(Tone):
        cs = sp.symbols("ch cf")
        eqs2 = [Tone[l][m][n] - (cs[0] * ref[l][m][n].subs({h: 1, f: 0}) + cs[1] * ref[l][m][n].subs({h: 0, f: 1}))
                for l, m, n in itertools.product(R4, R4, R4)]
        return bool(sp.solve(eqs2, cs, dict=True))
    ok &= report("solution space = span{T^i_0j = h delta_ij, T^i_jk = f [ijk]}",
                 all(is_span(mp[s]) for s in free))
    # check reference torsion is invariant
    inv = True
    for j, k in ((1, 2), (2, 3), (3, 1)):
        xi = [0, 0, 0, 0]
        xi[j], xi[k] = X[k], -X[j]
        L = lie(ref, xi)
        inv &= all(L[l][m][n] == 0 for l, m, n in itertools.product(R4, R4, R4))
    ok &= report("reference torsion is rotation invariant", inv)
    a = sp.Symbol("a", positive=True)
    g = a**2 * ETA
    S = [sp.simplify(v) for v in axial_vector(g, ref)]
    Tv = [sp.simplify(v) for v in trace_vector(ref)]
    print("INFO: on g = a^2 eta: S^m =", S, ", T_m =", Tv)
    ok &= report("S^tau = -f/a^2 (so f = -a^2 S^tau), S^i = 0", S == [-f / a**2, 0, 0, 0])
    ok &= report("T_tau = 3h, T_i = 0", Tv == [3 * h, 0, 0, 0])
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
