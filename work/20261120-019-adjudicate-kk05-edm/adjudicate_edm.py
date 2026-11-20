"""Task 20261120-019: independent adjudication of the KK-05 electron EDM (012 vs 015 vs 016).

5D frame, eta = diag(-1,1,1,1,1), index 4 = fifth direction. Clifford {G^A, G^B} = -2 eta^AB
(Bjorken-Drell gammas for A = 0..3), G^4 = s*i*gamma5 with s = +-1 (the two inequivalent 5D reps).
KK vielbein (radion fixed): e^a = dx^a, e^4 = dy + k A, W = de^4 = k F (constant frame components).
Dual frame: E_a = d_a - k A_a d_y, E_4 = d_y.

Route (independent of 016's code):
  * Levi-Civita spin connection obtained by SOLVING de^A + w^A_B ^ e^B = 0 with sympy (016 wrote it by hand).
  * Fermion spin term taken from the symmetrised Lagrangian (i/2)[Psibar G^C D_C Psi - (D_C Psi)bar G^C Psi]
    = Psibar (i/2){G^C, Omega_C} Psi  (016 used a numerical hermitian projection).
  * Exact sympy algebra with symbolic W_ab, c, theta; nonrelativistic limit read off the upper 2x2 block
    of gamma0*H (positive-energy spinors satisfy gamma0 u = u for L = psibar(i gamma.d - m)psi, m > 0).
  * Sanity check against literature: massive KK modes of a massless 5D Dirac fermion must have g = 1.
Contortion response K_[ab4] = -(4/3)c^2 W_ab +- c Wt_ab is taken from 016 D4 / 012 U3 / 009 E8 (three
derivations agree; not in dispute).
"""
import csv
import itertools
import os

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data")
NFAIL = 0


def check(name, ok, detail=""):
    global NFAIL
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    NFAIL += (not ok)


def info(msg):
    print(f"[INFO] {msg}")


def load(name):
    with open(os.path.join(DATA, name), newline="", encoding="utf-8-sig") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


D = 5
ETA = [-1, 1, 1, 1, 1]
I = sp.I


def perm_sign(p):
    p, s = list(p), 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


# ------------------------------------------------------------------ gamma matrices
s0, Z2 = sp.eye(2), sp.zeros(2)
PAULI = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.Matrix([[1, 0], [0, -1]])]


def blk(a, b, c_, d):
    return sp.Matrix(sp.BlockMatrix([[a, b], [c_, d]]))


g0 = blk(s0, Z2, Z2, -s0)
gi = [blk(Z2, s, -s, Z2) for s in PAULI]
g5 = I * g0 * gi[0] * gi[1] * gi[2]
E4 = sp.eye(4)


def gammas(s):
    return [g0] + gi + [s * I * g5]


# ------------------------------------------------------------------ symbolic field strength
Wsym = sp.symbols("W01 W02 W03 W12 W13 W23", real=True)
W = sp.zeros(D, D)
for (a, b), v in zip(itertools.combinations(range(4), 2), Wsym):
    W[a, b], W[b, a] = v, -v
EPS4 = {p: perm_sign(p) for p in itertools.permutations(range(4))}
Wt = sp.zeros(D, D)          # Wt_ab = (1/2) eps_abcd W^cd, eps_0123 = +1
for a in range(4):
    for b in range(4):
        Wt[a, b] = sp.Rational(1, 2) * sum(EPS4.get((a, b, cc, d), 0) * ETA[cc] * ETA[d] * W[cc, d]
                                           for cc in range(4) for d in range(4))


def levi_civita_connection():
    """Solve de^A = -w^A_B ^ e^B for w_ABC (antisym in AB), with de^4 = (1/2) W_ab e^a e^b, de^a = 0."""
    syms, w = [], {}
    for A, B in itertools.combinations(range(D), 2):
        for C in range(D):
            x = sp.Symbol(f"w_{A}{B}{C}")
            syms.append(x)
            w[A, B, C], w[B, A, C] = x, -x
    for A in range(D):
        for C in range(D):
            w[A, A, C] = 0
    eqs = []
    for A in range(D):
        for C, Dd in itertools.combinations(range(D), 2):
            de = W[C, Dd] if (A == 4) else 0
            # (w^A_B ^ e^B)_{CD} = w^A_{B C} delta^B_D - w^A_{B D} delta^B_C = w^A_{D C} - w^A_{C D}
            wedge = ETA[A] * (w[A, Dd, C] - w[A, C, Dd])
            eqs.append(sp.expand(de + wedge))
    sol = sp.solve(eqs, syms, dict=True)[0]
    return {key: sp.sympify(val).subs(sol) for key, val in w.items()}


def alt(w, A, B, C):
    idx = (A, B, C)
    return sum(perm_sign(p) * w[idx[p[0]], idx[p[1]], idx[p[2]]] for p in itertools.permutations(range(3))) / 6


def spin_generators(G, w_test):
    """Pick the sign of Sigma^AB = sgn (1/4)[G^A, G^B] for which G^A is covariantly constant."""
    for sgn in (1, -1):
        Sig = {(A, B): sgn * sp.Rational(1, 4) * (G[A] * G[B] - G[B] * G[A]) for A in range(D) for B in range(D)}
        Om = sum((sp.Rational(1, 2) * w_test[A][B] * Sig[A, B] for A in range(D) for B in range(D)), sp.zeros(4))
        if all(sp.simplify(Om * G[A] - G[A] * Om + sum((ETA[A] * w_test[A][B] * G[B] for B in range(D)), sp.zeros(4)))
               == sp.zeros(4) for A in range(D)):
            return Sig, sgn
    raise RuntimeError("no consistent spinor generator sign")


def spin_term(w, G, Sig):
    """H with L contains Psibar H Psi from (i/2){G^C, Omega_C}, Omega_C = (1/2) w_ABC Sigma^AB."""
    H = sp.zeros(4)
    for C in range(D):
        Om = sum((sp.Rational(1, 2) * w[A, B, C] * Sig[A, B] for A in range(D) for B in range(D)), sp.zeros(4))
        H += I / 2 * (G[C] * Om + Om * G[C])
    return H.applyfunc(sp.expand)


def nr_coeffs(H):
    """Upper 2x2 block of gamma0 H, expanded as sum_i (aE_i E_i + aB_i B_i) sigma_i.
    W_0i = -E_i, W_23 = B_x, W_31 = B_y, W_12 = B_z. Returns (D_E, D_B) defined by L contains phi^dag (D_E s.E + D_B s.B) phi."""
    blk2 = (g0 * H)[:2, :2]
    comps = [sp.expand(sp.Rational(1, 2) * (s * blk2).trace()) for s in PAULI]
    Wsub_E = {Wsym[0]: -1, Wsym[1]: 0, Wsym[2]: 0, Wsym[3]: 0, Wsym[4]: 0, Wsym[5]: 0}   # E_x = 1
    Wsub_B = {Wsym[0]: 0, Wsym[1]: 0, Wsym[2]: 0, Wsym[3]: 0, Wsym[4]: 0, Wsym[5]: 1}    # B_x = 1
    unit = sp.expand(sp.Rational(1, 2) * blk2.trace())
    return (sp.simplify(comps[0].subs(Wsub_E)), sp.simplify(comps[0].subs(Wsub_B)),
            [sp.simplify(x.subs(Wsub_E)) for x in comps[1:]] + [sp.simplify(x.subs(Wsub_B)) for x in comps[1:]],
            sp.simplify(unit))


def main():
    c, th = sp.symbols("c theta", real=True)

    # A1: Levi-Civita connection by solving the structure equation
    w0 = levi_civita_connection()
    a0_ok = all(sp.simplify(alt(w0, a, b, 4) - W[a, b] / 6) == 0 for a in range(4) for b in range(4))
    a0_4d = all(sp.simplify(alt(w0, a, b, cc)) == 0 for a, b, cc in itertools.combinations(range(4), 3))
    check("A1 solved KK Levi-Civita connection: w0_[ab4] = W_ab/6, w0_[abc] = 0 (agrees with 004/016)", a0_ok and a0_4d)

    # contortion response, totally antisymmetric (only that part reaches the fermion; checked in A3)
    def full_connection(sgK):
        w = dict(w0)
        Ka = -sp.Rational(4, 3) * c ** 2 * W + sgK * c * Wt
        for a in range(4):
            for b in range(4):
                if a == b:
                    continue
                for p in itertools.permutations(range(3)):
                    idx = (a, b, 4)
                    key = (idx[p[0]], idx[p[1]], idx[p[2]])
                    w[key] = w[key] + perm_sign(p) * Ka[a, b] / 2   # each unordered {a,b} visited twice
        return w

    rows = {}
    for s in (1, -1):
        G = gammas(s)
        cliff = all(sp.simplify(G[A] * G[B] + G[B] * G[A] + 2 * ETA[A] * (A == B) * E4) == sp.zeros(4)
                    for A in range(D) for B in range(D))
        wt = [[0] * D for _ in range(D)]
        vals = [3, -1, 2, 5, -4, 1, 7, -2, 3, 1]
        for n_, (A, B) in enumerate(itertools.combinations(range(D), 2)):
            wt[A][B], wt[B][A] = vals[n_], -vals[n_]
        Sig, sgn = spin_generators(G, wt)
        herm_ok = all(sp.simplify(g0 * Sig[A, B].H * g0 + Sig[A, B]) == sp.zeros(4) for A in range(D) for B in range(D))
        check(f"A2 (s={s:+d}) Clifford algebra; Sigma sign {sgn:+d}; g0 Sigma^dag g0 = -Sigma (symmetrised form valid)",
              cliff and herm_ok)

        # A3: only w_[ABC] survives in (i/2){G, Omega}
        wr = {}
        for A in range(D):
            for B in range(D):
                for C in range(D):
                    wr[A, B, C] = sp.Rational((3 * A + 5 * B * B + 7 * C + A * B * C) % 11 - 5) if A < B else 0
        for A in range(D):
            for B in range(A):
                for C in range(D):
                    wr[A, B, C] = -wr[B, A, C]
        wa = {(A, B, C): alt(wr, A, B, C) for A in range(D) for B in range(D) for C in range(D)}
        check(f"A3 (s={s:+d}) only the totally antisymmetric connection enters the fermion term",
              sp.simplify(spin_term(wr, G, Sig) - spin_term(wa, G, Sig)) == sp.zeros(4))

        for sgK in (1, -1):
            w = full_connection(sgK)
            H = spin_term(w, G, Sig)
            # general mass phase: L = psibar(i g.d - m e^{i theta g5})psi; rotate psi = e^{-i theta g5/2} psi'
            U = sp.cos(th / 2) * E4 - I * sp.sin(th / 2) * g5
            Hr = (U * H * U).applyfunc(sp.expand)
            DE, DB, off, unit = nr_coeffs(Hr)
            rows[s, sgK] = (sp.simplify(DE), sp.simplify(DB), off, unit)

    # A4: zero mode, real mass (theta = 0)
    DEcands = set()
    ok4 = True
    for (s, sgK), (DE, DB, off, unit) in rows.items():
        DE0, DB0 = sp.simplify(DE.subs(th, 0)), sp.simplify(DB.subs(th, 0))
        info(f"s={s:+d}, sign(Wt in K)={sgK:+d}: theta=0  D_E = {DE0},  D_B = {DB0};  "
             f"general theta: D_E = {sp.simplify(sp.trigsimp(DE))}, D_B = {sp.simplify(sp.trigsimp(DB))}")
        ok4 &= sp.simplify(DE0 ** 2 - (sp.Rational(3, 2) * (sp.Rational(1, 6) - sp.Rational(4, 3) * c ** 2)) ** 2) == 0
        ok4 &= sp.simplify(DB0 ** 2 - (sp.Rational(3, 2) * c) ** 2) == 0
        ok4 &= all(sp.simplify(x) == 0 for x in off) and unit == 0
        # general theta formula of 016: |D_E cos th +- D_B sin th|, and D_B' = -+D_E sin th + D_B cos th
        dE = sp.Rational(3, 2) * (sp.Rational(1, 6) - sp.Rational(4, 3) * c ** 2)
        dB = sp.Rational(3, 2) * c
        ok4 &= any(sp.simplify(sp.trigsimp(DE - (sa * dE * sp.cos(th) + sb * dB * sp.sin(th)))) == 0
                   for sa in (1, -1) for sb in (1, -1))
    check("A4 real-mass zero mode: D_E = (3/2)(1/6 - (4/3)c^2) (P'-even pieces -> EDM), D_B = (3/2)c "
          "(P'-odd piece -> MDM), both reps, both K signs; mass phase mixes them as D_E cos th +- D_B sin th "
          "(= 016 D5/D6)", ok4)
    info("012 used (3/4)c for the c-odd piece: the 3/4 counts each unordered pair (a,b) once in "
         "(3/4)K_[ab5] Psibar G^{ab5} Psi; summing over both orderings gives 3/2, as found here.")

    # A5: literature sanity check, massive KK modes of a massless bulk fermion must have g = 1
    #   Psi ~ e^{i n y/R}: i G^4 d_y -> -(n/R) G^4 = -(n/R) s i g5, i.e. mass operator m e^{i th g5} with
    #   m = |n|/R, th = s*sign(n)*pi/2. Charge: i psibar g^a (-k A_a)(i n/R) psi = (k n/R) psibar g^a psi A_a with
    #   mostly-plus A_a = (-Phi, A): L contains -q Phi psi^dag psi, q = k n/R. Dirac part gives L contains (q/2m) s.B (g = 2).
    #   g = 2 + 4 m k D_B'/q with D_B' the sigma.B coefficient per unit W of the rotated Pauli term at c = 0.
    ok5 = True
    for s in (1, -1):
        for nsgn in (1, -1):
            DE, DB, _, _ = rows[s, 1]
            thv = s * nsgn * sp.pi / 2
            DBp = sp.simplify(DB.subs({th: thv, c: 0}))
            DEp = sp.simplify(DE.subs({th: thv, c: 0}))
            g = sp.simplify(2 + 4 * nsgn * DBp)       # m/q = 1/(k*sign(n))
            info(f"A5 s={s:+d}, sign(n)={nsgn:+d}: KK mode at c=0 has D_E' = {DEp}, D_B' = {DBp} -> g = {g}")
            ok5 &= (g == 1 and DEp == 0)
    check("A5 massless 5D fermion, KK modes n != 0: the Levi-Civita W/6 Pauli term is purely MAGNETIC and gives "
          "g = 1 (literature: universal g = 1 of massive KK states, Ohkuwa 1984 thesis; Hosoya et al. PLB 134 (1984) 44)",
          ok5)
    info("So the same W/6 term is an EDM only because a real 5D mass M (P_4-odd, Thirring 1972) sets theta = 0. "
         "A KK mode with both M and n has tan(theta) = n/(M R) and EDM ~ M/sqrt(M^2+n^2/R^2).")

    # A6: canonical normalisation, numbers
    shu = load("shift_uholst_inputs.csv")
    binp = load("bounds_inputs.csv")
    cosmo = load("cosmo_torsion_axion_inputs.csv")
    import math
    MP = float(binp["MP_reduced"]["value"])
    HBARC = float(binp["hbarc"]["value"])
    alpha = float(cosmo["alpha_em"]["value"])
    dJ = float(shu["eEDM_JILA"]["value"])
    k = math.sqrt(2) / MP
    e = math.sqrt(4 * math.pi * alpha)
    cm = HBARC * 100

    def d_ecm(cv, thv):
        dE = 1.5 * (1 / 6 - 4 / 3 * cv ** 2)
        dB = 1.5 * cv
        return k * max(abs(dE * math.cos(thv) + dB * math.sin(thv)), abs(dE * math.cos(thv) - dB * math.sin(thv))) \
            / math.sqrt(1 + 12 * cv ** 2) / e * cm

    lo, hi = 1.0, 1e6
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if d_ecm(mid, 0.0) < dJ else (lo, mid)
    cmax = math.sqrt(lo * hi)
    info(f"A6 k/e = {k / e * cm:.3e} cm; d(c=0) = {d_ecm(0, 0):.2e} e cm; large-c slope (theta=0) "
         f"{d_ecm(1e5, 0) / 1e5:.3e} e cm per unit c")
    check("A6 theta = 0 (real-mass zero mode + non-geometric charge): JILA gives c <= 188 (016 E2 reproduced)",
          185 < cmax < 191, f"c_max = {cmax:.1f}")
    sat = 1.5 * k / math.sqrt(12) / e * cm
    worst = max(d_ecm(cv, math.pi / 2) for cv in [10 ** (x / 10) for x in range(-30, 80)])
    info(f"A7 theta = pi/2 (KK-charged mode, M = 0): d = (3/2)c k/sqrt(1+12c^2) -> {sat:.2e} e cm; "
         f"max over c in [1e-3,1e8] = {worst:.2e} e cm; JILA {dJ:.1e}")
    check("A7 for the only gauge-consistent charge (KK momentum, theta -> pi/2) the EDM saturates "
          "below JILA: no bound on c (015's form, with 3/2 in place of 3/4)", worst < dJ, f"ratio {worst / dJ:.1e}")
    for cv in (1.0, 188.0, 5.1e4):
        mmin = math.sqrt(4 * math.pi * alpha) * math.sqrt(1 + 12 * cv ** 2) / k
        info(f"A8 KK-charged electron with physical charge e at c = {cv:g}: m >= e sqrt(Z)/k = {mmin:.2e} eV")
    print(f"SUMMARY: {NFAIL} FAIL")


if __name__ == "__main__":
    main()
