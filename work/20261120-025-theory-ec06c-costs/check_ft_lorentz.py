"""f(T) on flat FRW: what fixes the frame, is T_mu closed on solutions, tensor speed (task 025).

L1  Rotations: the diagonal FRW tetrad rotated by omega(tau,x,y,z) about z gives the same metric, T_tau and
    torsion scalar, but T_i gains a curl (-d_y omega, d_x omega, 0), (dT)_{xy} = +(d_x^2 + d_y^2) omega.
    The f(T) antisymmetric equation f_TT S_[mu nu]^tau T' = 0 is violated unless grad omega = 0, so f(T)
    (f_TT != 0) removes this curl at linear order; in TEGR (f_TT = 0) it is free. The CFJ term adds a
    photon source (xi/2) curl(eK) to the same equation (photon back-reaction on the frame).
L2  Scalar boost sector: tetrad = boost(eps*eta(tau,x)) x diag perturbed tetrad. Second-order action
    e f(T): TEGR part (f = T) is eta-independent up to a boundary term; the eta equation is
    f_TT * [ ... ] = 0 and fixes eta algebraically. On that solution (dT)_{tau x} != 0 even in matter
    domination with static Phi = Psi.
L3  Tensor modes: second-order action gives h'' + (2a'/a + f_T'/f_T) h' - d_z^2 h = 0, so c_T = 1 exactly
    for any f(T) on flat FRW (only the friction changes).
L4  Static weak field: torsion scalar of the diagonal tetrad with static Phi = Psi is +2 (grad Phi)^2
    (same sign as the FRW value +6H^2): T_sc ~ 6H^2 + 2|grad Phi|^2 has no zero, but near any mass
    |grad Phi|^2 >> H^2, so ln T_sc is set by the local field (relevant to the ln T completion, repair c).
"""
import sympy as sp

from common import EPS, R4, check, is_closed, series_inv, summary, torsion_scalar, torsion_series, trace_vec, \
    trunc, weitzenbock_torsion

tau, x, y, z = sp.symbols("tau x y z", real=True)
X = [tau, x, y, z]
a = sp.Function("a", positive=True)(tau)
ap = sp.diff(a, tau)
eta_m = sp.diag(-1, 1, 1, 1)


def superpot(Tt, gm):
    """Same superpotential as task 022 check_teleparallel.py (S_r^{mn})."""
    gi = gm.inv()
    tv = trace_vec(Tt)
    tvu = [sum(gi[n, b] * tv[b] for b in R4) for n in R4]
    A = [[[sum(gi[n, b] * Tt[m][b][r] for b in R4) for r in R4] for n in R4] for m in R4]
    B = [[[sum(gm[r, al] * gi[m, be] * gi[n, ga] * Tt[al][be][ga] for al in R4 for be in R4 for ga in R4)
           for n in R4] for m in R4] for r in R4]
    Sp = [[[None] * 4 for _ in R4] for _ in R4]
    for r in R4:
        for m in R4:
            for n in R4:
                K = -sp.Rational(1, 2) * (A[m][n][r] - A[n][m][r] - B[r][m][n])
                Sp[r][m][n] = sp.Rational(1, 2) * (K + int(r == m) * tvu[n] - int(r == n) * tvu[m])
    Tsc = sp.simplify(sum(Sp[r][m][n] * Tt[r][m][n] for r in R4 for m in R4 for n in R4))
    return Sp, Tsc


# ---------------------------------------------------------------- L1 remnant rotations
e0 = sp.diag(a, a, a, a)
g0 = sp.diag(-a**2, a**2, a**2, a**2)
T0 = weitzenbock_torsion(e0, X)
Sp0, Ts0 = superpot(T0, g0)
for label, om in (("omega(x,y)", sp.Function("omega")(x, y)), ("omega(tau,x,y,z)", sp.Function("omega")(tau, x, y, z))):
    Rm = sp.eye(4)
    Rm[1, 1] = Rm[2, 2] = sp.cos(om)
    Rm[1, 2], Rm[2, 1] = -sp.sin(om), sp.sin(om)
    e1 = Rm * e0
    g1 = sp.simplify(e1.T * eta_m * e1)
    T1 = weitzenbock_torsion(e1, X)
    tv1 = trace_vec(T1)
    Sp1, Ts1 = superpot(T1, g1)
    same_metric = sp.simplify(g1 - g0) == sp.zeros(4)
    print(f"   rotation {label}: T_mu = {[sp.simplify(v) for v in tv1]}, T_sc = {Ts1}")
    check(f"L1a {label}: same metric, same T_tau = 3a'/a, same torsion scalar",
          same_metric and sp.simplify(tv1[0] - 3 * ap / a) == 0 and sp.simplify(Ts1 - Ts0) == 0)
    Y0 = sp.Matrix(4, 4, lambda m, n: sum(g0[n, k] * Sp0[m][k][0] for k in R4))
    Y1 = sp.Matrix(4, 4, lambda m, n: sp.simplify(sum(g1[n, k] * Sp1[m][k][0] for k in R4)))
    dY = sp.simplify(Y1 - Y0)
    print(f"   rotation {label}: Y = g_nk S_m^(k tau) change =", dY.tolist())
    A1 = sp.simplify(Y1 - Y1.T)
    check(f"L1b {label}: antisymmetric S_[tau i]^tau = (d_y w, -d_x w)/2: f(T) Lorentz eq forces grad w = 0",
          sp.simplify(A1[0, 1] - sp.diff(om, y) / 2) == 0 and sp.simplify(A1[0, 2] + sp.diff(om, x) / 2) == 0
          and sp.simplify(A1[0, 3]) == 0 and sp.simplify(A1[1, 2]) == 0)
    curl = sp.simplify(sp.diff(tv1[2], x) - sp.diff(tv1[1], y))
    check(f"L1c {label}: T_i = (-d_y w, d_x w, 0), (dT)_xy = +(d_x^2+d_y^2) w: not closed",
          sp.simplify(tv1[1] + sp.diff(om, y)) == 0 and sp.simplify(tv1[2] - sp.diff(om, x)) == 0
          and sp.simplify(curl - sp.diff(om, x, 2) - sp.diff(om, y, 2)) == 0 and not is_closed(tv1, X))
# CFJ variation under an infinitesimal remnant rotation: delta S = -1/2 xi int e dT_i K^i
xi = sp.Symbol("xi")
w = sp.Function("w")(tau, x, y, z)
Kx, Ky = sp.Function("Kx")(tau, x, y, z), sp.Function("Ky")(tau, x, y, z)  # densities e K^i
dL = -xi / 2 * (-sp.diff(w, y) * Kx + sp.diff(w, x) * Ky)
EL = sp.euler_equations(dL, [w], X)[0].lhs
check("L1d CFJ contribution to the omega (rotation) equation: (xi/2)(d_y K^x - d_x K^y) (photon source for the frame)",
      sp.simplify(EL - xi / 2 * (sp.diff(Kx, y) - sp.diff(Ky, x))) == 0 or
      sp.simplify(EL + xi / 2 * (sp.diff(Kx, y) - sp.diff(Ky, x))) == 0, f"({sp.simplify(EL)})")

# ---------------------------------------------------------------- L2 boost sector, second-order action
Phi = sp.Function("Phi")(tau, x)
Psi = sp.Function("Psi")(tau, x)
eta = sp.Function("eta")(tau, x)
base = sp.diag(a * (1 + EPS * Phi), a * (1 - EPS * Psi), a * (1 - EPS * Psi), a * (1 - EPS * Psi))
Lb = sp.eye(4)
Lb[0, 0] = Lb[1, 1] = 1 + EPS**2 * eta**2 / 2
Lb[0, 1] = Lb[1, 0] = EPS * eta
eb = (Lb * base).applyfunc(lambda v: trunc(v))
gb = (eb.T * eta_m * eb).applyfunc(lambda v: trunc(v))
gbi = series_inv(gb)
Tb = torsion_series(eb, X)
Tsc, vb = torsion_scalar(Tb, gb, gbi)
edet = trunc(a**4 * (1 + EPS * Phi) * (1 - EPS * Psi)**3)
Tbar = 6 * ap**2 / a**4
check("L2a background torsion scalar = 6 a'^2/a^4 = 6H^2", sp.simplify(Tsc.coeff(EPS, 0) - Tbar) == 0)
T1s, T2s = sp.expand(Tsc).coeff(EPS, 1), sp.expand(Tsc).coeff(EPS, 2)
print("   delta_1 T_sc =", sp.simplify(T1s))

# first-order trace vector including the boost
v1 = [sp.expand(c).coeff(EPS, 1) for c in vb]
print("   delta_1 T_mu =", [sp.simplify(c) for c in v1])
exp_v = [-3 * sp.diff(Psi, tau) - sp.diff(eta, x), sp.diff(Phi - 2 * Psi, x) - sp.diff(eta, tau), 0, 0]
check("L2b boosted trace: dT_tau = -3Psi' - eta_x, dT_x = (Phi-2Psi)_x - eta'",
      all(sp.simplify(v1[m] - exp_v[m]) == 0 for m in R4))
dTtx = sp.simplify(sp.diff(v1[1], tau) - sp.diff(v1[0], x))
print("   (dT)_{tau x} =", dTtx)

# action: e f(T) to second order; F1 = f_T(Tbar(tau)), F2 = f_TT(Tbar(tau)), f0 = f(Tbar)
F1, F2, f0 = sp.Function("F1")(tau), sp.Function("F2")(tau), sp.Function("f0")(tau)
e1d, e2d = sp.expand(edet).coeff(EPS, 1), sp.expand(edet).coeff(EPS, 2)
L2 = sp.expand(a**4 * (F1 * T2s + F2 * T1s**2 / 2) + e1d * F1 * T1s + e2d * f0)
ELeta = sp.euler_equations(L2, [eta], [tau, x])[0].lhs
rules = {sp.Derivative(F1, tau): F2 * sp.diff(Tbar, tau), sp.Derivative(f0, tau): F1 * sp.diff(Tbar, tau)}
ELeta = sp.simplify(sp.expand(ELeta.subs(rules)))
check("L2c TEGR (F2 = 0, F1 const) gives no eta equation (frame is pure gauge)",
      sp.simplify(ELeta.subs(F2, 0)) == 0)
red = sp.factor(sp.simplify(ELeta / F2))
print("   eta equation / f_TT =", red)
check("L2d f(T): eta equation is proportional to f_TT", sp.simplify(sp.diff(red, F2)) == 0)
# general solution: eta_xx fixed by Phi, Psi; then (dT)_{tau x}
exx = sp.Symbol("exx")
sol = sp.solve(red.subs(sp.Derivative(eta, (x, 2)), exx), exx)
check("L2e eta_xx is fixed uniquely by Phi, Psi (no free function once f_TT != 0)", len(sol) == 1)
target = -3 * sp.diff(Psi, tau, x) + 3 * sp.diff(a, tau, 2) / ap * sp.diff(Psi, x) - 3 * ap / a * sp.diff(Phi + 2 * Psi, x)
check("L2f eta_xx = d_x[-3Psi' + 3(a''/a')Psi - 3(a'/a)(Phi + 2Psi)]", sp.simplify(sol[0] - target) == 0)
# matter era, a = c tau^2, Phi = Psi = phi0 cos(kx) static; eta = E(tau) sin(kx)
c0, k, phi0 = sp.symbols("c0 k phi0", positive=True)
Et = sp.Function("E")(tau)
md = {Phi: phi0 * sp.cos(k * x), Psi: phi0 * sp.cos(k * x), eta: Et * sp.sin(k * x)}
eq = sp.simplify((red.subs(md).doit().subs(a, c0 * tau**2).doit()) / sp.sin(k * x))
Esol = sp.solve(eq, Et)
print("   matter era: E(tau) =", Esol)
dT_md = sp.simplify(dTtx.subs(md).subs(Et, Esol[0]).doit().subs(a, c0 * tau**2).doit())
print("   matter era, on the f(T) solution: (dT)_{tau x} =", sp.factor(dT_md))
check("L2g matter era, static Phi = Psi: (dT)_{tau x} = 15 k phi0 sin(kx)/tau (1 + 2/(k tau)^2) != 0",
      sp.simplify(dT_md - 15 * k * phi0 * sp.sin(k * x) / tau * (1 + 2 / (k * tau)**2)) == 0)
check("L2h (cross-check with 022 T5) at eta = 0, static Phi = Psi, dT = 0",
      sp.simplify(dTtx.subs({Phi: phi0 * sp.cos(k * x), Psi: phi0 * sp.cos(k * x), eta: 0}).doit()) == 0)

# ---------------------------------------------------------------- L3 tensor modes
h = sp.Function("h")(tau, z)
et = sp.diag(a, a * (1 + EPS * h / 2), a * (1 - EPS * h / 2), a)
gt = (et.T * eta_m * et).applyfunc(lambda v: trunc(v))
gti = series_inv(gt)
Tt = torsion_series(et, X)
Tsc_t, vt = torsion_scalar(Tt, gt, gti)
check("L3a tensor mode: delta_1 T_sc = 0 and delta_1 T_mu = 0", sp.expand(Tsc_t).coeff(EPS, 1) == 0
      and all(sp.expand(c).coeff(EPS, 1) == 0 for c in vt))
edt = trunc(et.det())
L2t = sp.expand(a**4 * F1 * sp.expand(Tsc_t).coeff(EPS, 2) + sp.expand(edt).coeff(EPS, 2) * f0)
ELh = sp.expand(sp.euler_equations(L2t, [h], [tau, z])[0].lhs.subs(rules))
c_tt = ELh.coeff(sp.diff(h, tau, 2))
c_zz = ELh.coeff(sp.diff(h, z, 2))
c_t = ELh.coeff(sp.diff(h, tau))
print("   tensor EL:", sp.simplify(ELh))
check("L3b c_T^2 = -coef(h_zz)/coef(h_tautau) = 1 for any f(T)", sp.simplify(-c_zz / c_tt - 1) == 0)
fric = sp.simplify(c_t / c_tt)
check("L3c friction = 2a'/a + f_TT Tbar'/f_T (only the damping changes)",
      sp.simplify(fric - 2 * ap / a - F2 * sp.diff(Tbar, tau) / F1) == 0, f"({fric})")

# ---------------------------------------------------------------- L4 static weak field
P = sp.Function("P")(x, y, z)
es = sp.diag(1 + EPS * P, 1 - EPS * P, 1 - EPS * P, 1 - EPS * P)
gs = (es.T * eta_m * es).applyfunc(lambda v: trunc(v))
Ts, _ = torsion_scalar(torsion_series(es, X), gs, series_inv(gs))
grad2 = sum(sp.diff(P, v)**2 for v in (x, y, z))
check("L4 static Phi = Psi = P: torsion scalar = +2 eps^2 |grad P|^2 (same sign as +6H^2)",
      sp.simplify(Ts - 2 * EPS**2 * grad2) == 0, f"({sp.factor(Ts)})")
summary()
