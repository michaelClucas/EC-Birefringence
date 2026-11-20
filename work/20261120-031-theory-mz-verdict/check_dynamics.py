"""Task 031, items 2-3: leading-order-in-gamma dynamics of the MZ axial torsion, theta normalisation, tensor singular scale.

T1  MZ Eq. 30 (FRW, from 1908.05184) at leading order in gamma (P = gamma p, T = O(gamma^2), Ldot = O(gamma^2))
    equals (3/(2 gamma)) a^-3 d_t[a^3 (Lambda + kappa rho) 2P] = (3/8) GB - Lambda^2,   GB = 24 H^2 (Hdot + H^2):
    the FRW form of  nabla_mu (G^{mu nu} S_nu) = (gamma/4) (GB - (8/3) Lambda^2).
T2  P = gamma H/3 - (gamma Lambda^2/9) int(a^3 dt)/(a^3 H^2) solves it exactly for any flat history (the tracker);
    the general solution adds C/(a^3 H^2).
T3  radiation tracker = 1908.05184/2001.06373 Eq. 27 (P = gamma H0 sqrt(Omega_r)/(3 a^2)); matter tracker
    P^2 = gamma^2 kappa rho/27 (2001.06373 Sec. VIII); homogeneous mode / tracker ~ a^3 in radiation (their Eq. 28).
T4  theta = 3/(2 Lambda) effective action. The theta-dependent action is -(M_P^2/2) theta [eps RR + (2/gamma) RR]
    (from 1908.05184 Eq. 7). With S = (2/(3 gamma)) G.dtheta, L_theta = -(N/2) Z^{mu nu} d theta d theta - (M_P^2/2) theta GB,
    Z = G.G, reproduces T1 iff N = 4 M_P^2/(3 gamma^2). Canonical phi = sqrt(N)(Lambda + kappa rho) theta gives
    S = -(1/(sqrt3 M_P)) [d phi - (1 - c_s)(d phi)_perp], c_s = Lambda/(Lambda + kappa rho)   (gamma drops out).
T5  Stueckelberg completion of the non-gradient piece: A -> A + d pi/m turns -(xi/2) S.A Ftilde into a pi-phi-gamma
    vertex with coefficient c = xi (1 - c_s)/(2 sqrt3 M_P m) (dimension -2).
T6  Tensor modes (2001.06373 Eqs. 44-48, (III)): the H'' coefficient of the tensor equation is A/(A^2 + B^2) and,
    with Ldot from Eq. 29 (exact), A = 0 at k_P = (Lambda + kappa rho)/(6 P) exactly (gamma-independent):
    k* = H/eps_S with eps_S = S_t/H = 2P/H in a flat universe.
"""
import sympy as sp

from common import check, summary

t = sp.symbols("t", real=True)
gam, kap, Lam = sp.symbols("gamma kappa Lambda", positive=True)
a = sp.Function("a", positive=True)(t)
H = sp.diff(a, t) / a
rho, p, pt = sp.Function("rho")(t), sp.Function("p")(t), sp.Function("ptilde")(t)

# T1
fried = {rho: (3 * H**2 - Lam) / kap}  # flat, leading order (torsion corrections O(gamma^2))
lhs30 = (Lam + kap * rho) * (Lam - kap * (rho + 3 * p) / 2) - Lam**2
rhs30 = 9 * ((Lam + kap * rho) / 3 * (sp.diff(a * pt, t) / a) + sp.Rational(2, 3) * (Lam - kap * (rho + 3 * p) / 2) * H * pt)
transport = sp.Rational(3, 2) * sp.diff(a**3 * (Lam + kap * rho) * 2 * pt, t) / a**3
cons = {sp.Derivative(rho, t): -3 * H * (rho + p)}
# pressure from the acceleration equation: Hdot = -kappa (rho + p)/2
p_sol = sp.solve(sp.Eq(sp.diff(H, t), -kap * (rho + p) / 2), p)[0]
GB = 24 * H**2 * (sp.diff(H, t) + H**2)
ok1 = sp.simplify((transport.subs(cons) - rhs30).subs(p, p_sol).subs(fried).doit()) == 0
ok2 = sp.simplify((lhs30 - sp.Rational(3, 8) * GB + Lam**2).subs(p, p_sol).subs(fried).doit()) == 0
check("T1 Eq. 30 at O(gamma^0) (P = gamma p~) is (3/2) a^-3 d_t[a^3 (L + k rho) 2 p~] = (3/8) GB - Lambda^2, "
      "i.e. nabla.(G.S) = (gamma/4) (GB - (8/3) Lambda^2)", ok1 and ok2)

# T2
P = sp.Function("P")(t)
eq = sp.diff(a**3 * 3 * H**2 * 2 * P, t) / a**3 - gam / 4 * (GB - sp.Rational(8, 3) * Lam**2)  # L + k rho = 3H^2
I3 = sp.Function("I3")(t)  # I3 = int a^3 dt
P_part = gam * H / 3 - gam * Lam**2 / 9 * I3 / (a**3 * H**2)
check("T2a P = gamma H/3 - (gamma Lambda^2/9) (int a^3 dt)/(a^3 H^2) is an exact particular solution for any flat a(t); "
      "the first term alone for Lambda = 0", sp.simplify(eq.subs(P, P_part).doit().subs(sp.Derivative(I3, t), a**3)) == 0)
aL = sp.exp(sp.sqrt(Lam / 3) * t)
check("T2c de Sitter (H^2 = Lambda/3, int a^3 dt = a^3/(3H)): particular P = 0, so on the tracker eps_S decays in the "
      "Lambda era (P -> P_i/a^3, 2001.06373)",
      sp.simplify((gam * sp.sqrt(Lam / 3) / 3 - gam * Lam**2 / 9 * (aL**3 / (3 * sp.sqrt(Lam / 3))) /
                   (aL**3 * Lam / 3))) == 0)
C = sp.Symbol("C")
check("T2b homogeneous solution P_h = C/(a^3 H^2): eps_h = 2 P_h/H ~ 1/(a H)^3 grows through radiation and matter eras",
      sp.simplify((sp.diff(a**3 * 3 * H**2 * 2 * P, t)).subs(P, C / (a**3 * H**2)).doit()) == 0)

# T3
H0, Or, Om = sp.symbols("H0 Omega_r Omega_m", positive=True)
ar = sp.Symbol("a", positive=True)
Hr = H0 * sp.sqrt(Or) / ar**2
check("T3a radiation: gamma H/3 = (gamma/(3 a^2)) H0 sqrt(Omega_r) (2001.06373 Eq. 27)",
      sp.simplify(gam * Hr / 3 - gam / (3 * ar**2) * H0 * sp.sqrt(Or)) == 0)
check("T3b matter: (gamma H/3)^2 = gamma^2 kappa rho/27 with kappa rho = 3 H^2 (2001.06373 Sec. VIII)",
      sp.simplify((gam * sp.Symbol("Hm") / 3)**2 - gam**2 * 3 * sp.Symbol("Hm")**2 / 27) == 0)
check("T3c radiation: P_h/P_tracker ~ a^3 (2001.06373 Eq. 28, delta_P = C a^3)",
      sp.simplify((1 / (ar**3 * Hr**2)) / Hr / ar**3) == sp.simplify(1 / (H0**3 * Or**sp.Rational(3, 2))))

# T4
MP, N = sp.symbols("M_P N", positive=True)
th = sp.Function("theta")(t)
G = Lam + kap * rho  # = -G_t^t
L = a**3 * (N / 2 * G**2 * sp.diff(th, t)**2 - MP**2 / 2 * th * GB - 3 * MP**2 / th)  # last: -(1/2k)(1/(4 theta)) eps eeee
EL = sp.euler_equations(L, [th], [t])[0].lhs
# transport in theta form: S_t = -(2/(3 gamma)) G thetadot = 2P  =>  d_t[a^3 G^2 thetadot] = -(3 gamma^2/8) a^3 GB
target = sp.diff(a**3 * G**2 * sp.diff(th, t), t) + sp.Rational(3, 8) * gam**2 * a**3 * (GB - 6 / th**2)  # (8/3) Lambda^2 = 6/theta^2
Nsol = sp.solve(sp.Eq(sp.simplify(EL - (-N) * target), 0), N)
check("T4a theta action normalisation N = 4 M_P^2/(3 gamma^2) reproduces the transport equation",
      any(sp.simplify(s - 4 * MP**2 / (3 * gam**2)) == 0 for s in Nsol), f"(N = {Nsol})")
Gs, cs = sp.symbols("G_s c_s", positive=True)
Nv = 4 * MP**2 / (3 * gam**2)
phidot, phigrad = sp.symbols("phidot phigrad")
thdot, thgrad = phidot / (sp.sqrt(Nv) * Gs), phigrad / (sp.sqrt(Nv) * Gs)
S_par = -sp.Rational(2, 3) / gam * Gs * thdot           # (kappa tau - Lambda)_t^t = -(Lambda + kappa rho)
S_perp = -sp.Rational(2, 3) / gam * (cs * Gs) * thgrad  # (kappa tau - Lambda)_i^j = -Lambda delta (dust)
check("T4b S = -(1/(sqrt3 M_P)) (phidot, c_s grad phi): gamma cancels; c_s = Lambda/(Lambda + kappa rho)",
      sp.simplify(S_par + phidot / (sp.sqrt(3) * MP)) == 0 and sp.simplify(S_perp + cs * phigrad / (sp.sqrt(3) * MP)) == 0)
check("T4c canonical kinetic term (N/2) G^2 thetadot^2 = phidot^2/2 and gradient term (N/2) Lambda^2 (grad theta)^2 "
      "= c_s^2 (grad phi)^2/2", sp.simplify(Nv / 2 * Gs**2 * thdot**2 - phidot**2 / 2) == 0
      and sp.simplify(Nv / 2 * (cs * Gs)**2 * thgrad**2 - cs**2 * phigrad**2 / 2) == 0)

# T5: -(xi/2) p_mu K^mu with the non-gradient piece p_i = (1 - c_s) d_i phi/(sqrt3 M_P) (sign irrelevant).
# Under A -> A + d sigma: delta(-(1/2) p.K) = (1/4) sigma (dp)_{nu mu} Ftilde^{mu nu} + total derivative (task 022 G1).
# (dp)_{t i} = (1-c_s) phi_ti/(sqrt3 M_P) * xi; sigma = pi/m.
xi, m = sp.symbols("xi m", positive=True)
coef = sp.Rational(1, 4) * 2 * xi * (1 - cs) / (sp.sqrt(3) * MP * m)  # two orderings (t,i),(i,t) in (dp).Ftilde
check("T5 pi-phi-gamma vertex coefficient c = xi (1 - c_s)/(2 sqrt3 M_P m)",
      sp.simplify(coef - xi * (1 - cs) / (2 * sp.sqrt(3) * MP * m)) == 0)

# T6: tensor equation structure from 2001.06373 Eqs. 44-48 and (III), N = 1
Hh = sp.Function("Hh")(t)
Af, Bf, kP, g, Hm = sp.Function("A")(t), sp.Function("B")(t), sp.Function("kP")(t), sp.Function("g")(t), sp.Function("Hm")(t)
den = Af**2 + Bf**2
x1, x2 = -kP * Hh, sp.diff(Hh, t) + (2 * Hm - g) * Hh
Bc = (Af * x1 + Bf * x2) / den
Ec = (-Bf * x1 + Af * x2) / den
cL, cr, cp = sp.symbols("c_L c_r c_p")
eq48 = sp.diff(Ec, t) + (Hm + g) * Ec - kP * Bc - cL * Hh
coef_dd = sp.simplify(sp.expand(eq48).coeff(sp.diff(Hh, t, 2)))
check("T6a coefficient of H'' in the tensor equation is A/(A^2 + B^2) (vanishes and flips sign where A = 0)",
      sp.simplify(coef_dd - Af / den) == 0)
gs, Ps, kPs, Ld = sp.symbols("g P k_P Ldot")
Ld29 = sp.solve(sp.Eq(Ps, 3 * Ld / Lam**2 * (gs * Ps + (Lam + kap * rho) / (6 * gam))), Ld)[0]
A44 = 1 - 3 * Ld / Lam**2 * (gs + kPs / gam)
kstar = sp.solve(sp.Eq(A44.subs(Ld, Ld29), 0), kPs)
check("T6b with Eq. 29 (exact), A = 0 at k_P* = (Lambda + kappa rho)/(6 P), independent of gamma",
      len(kstar) == 1 and sp.simplify(kstar[0] - (Lam + kap * rho) / (6 * Ps)) == 0, f"(k_P* = {kstar})")
Hs, eS = sp.symbols("H_s eps_S", positive=True)
check("T6c flat: Lambda + kappa rho = 3H^2 (+O(P^2)), eps_S = 2P/H  =>  k* = H/eps_S",
      sp.simplify((3 * Hs**2 / (6 * Ps)).subs(Ps, eS * Hs / 2) - Hs / eS) == 0)
summary()
