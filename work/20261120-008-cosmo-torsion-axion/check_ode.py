"""Checks of the field evolution and recombination (PASS/FAIL lines)."""
import numpy as np
import sympy as sp
from scipy.special import jv, gamma
from scipy.integrate import cumulative_simpson

from common import E, H0_EV, OMR, inp, _LNA, _A, TAU_TAB, tau_of_a
from field_evolution import evolve, phi_of_a, rho_today_over_phii2
from recombination import visibility


def report(name, ok, detail):
    print(f"{'PASS' if ok else 'FAIL'}: {name} -- {detail}")


# 1. conformal-time equation == cosmic-time Klein-Gordon (sympy)
tau, m = sp.symbols("tau m", positive=True)
a = sp.Function("a")(tau)
phi = sp.Function("phi")(tau)
# cosmic time: d/dt = (1/a) d/dtau; H = a'/a^2
H = sp.diff(a, tau) / a**2
phidot = sp.diff(phi, tau) / a
phiddot = sp.diff(phidot, tau) / a
kg = phiddot + 3 * H * phidot + m**2 * phi
conf = sp.diff(phi, tau, 2) + 2 * sp.diff(a, tau) / a * sp.diff(phi, tau) + a**2 * m**2 * phi
report("conformal eq. = a^2 x (phi_tt + 3H phi_t + m^2 phi)", sp.simplify(conf - a**2 * kg) == 0,
       "sympy simplify")

# 2. slow-roll (small-m) perturbative limit: 1 - phi(tau) = mu^2 int dtau a^-2 int dtau a^4
mu = 1e-35 / H0_EV
J = cumulative_simpson(_A**3 / E(_A), x=_LNA, initial=0.0)          # int a^4 dtau
I = cumulative_simpson(J / (_A**3 * E(_A)), x=_LNA, initial=0.0)     # int J a^-2 dtau
out = evolve(mu)
phi0 = out["sol"].y[1, -1]
pred0 = mu**2 * I[-1]
astar = 1 / (1 + inp("z_star"))
phis = phi_of_a(out, np.array([astar]))[0]
preds = mu**2 * np.interp(np.log(astar), _LNA, I)
r0 = abs((1 - phi0) / pred0 - 1)
report("slow-roll limit today (m = 1e-35 eV): beta/(g phi_i) = mu^2 [I(tau0) - I(tau*)]/2",
       r0 < 1e-3, f"ODE 1-phi0 = {1-phi0:.6e}, analytic {pred0:.6e}, rel diff {r0:.1e}")
bode = -(phi0 - phis) / 2
bana = mu**2 * (I[-1] - np.interp(np.log(astar), _LNA, I)) / 2
report("slow-roll beta at LSS-to-today", abs(bode / bana - 1) < 2e-3,
       f"ODE {bode:.6e}, analytic {bana:.6e}")

# 3. pure radiation background: exact solution phi = Gamma(5/4) (2/(mu t))^{1/4} J_{1/4}(mu t)
Erad = lambda aa: 1.0 / aa**2
mu_r = 1.0
o = evolve(mu_r, a_end=5.0, switch=1e12, Efun=Erad, a_init=1e-4, rtol=1e-11)
errs = []
for aa in (1.0, 2.0, 3.0, 4.5):
    t = aa**2 / 2
    exact = gamma(1.25) * (2 / (mu_r * t)) ** 0.25 * jv(0.25, mu_r * t)
    num = o["sol"].sol(aa)[1]          # tau = a in this background
    errs.append(abs(num - exact))
report("radiation era vs Bessel J_{1/4} solution", max(errs) < 1e-6, f"max |diff| = {max(errs):.1e}")

# 4. WKB continuation vs full ODE through recombination (m = 1e-25 eV)
mu_w = 1e-25 / H0_EV
ow = evolve(mu_w)                       # switches at m/H = 300 (before recombination)
of = evolve(mu_w, a_end=1 / 801.0, switch=1e12, rtol=1e-11)
aa = np.array([1 / 1301.0, 1 / 1091.0, 1 / 901.0])
# full ODE value at given a
taus = of["sol"].t
tq = np.interp(np.log(aa), np.log(of["sol"].y[0]), taus)
for _ in range(4):
    st = of["sol"].sol(tq)
    tq = tq - (st[0] - aa) / (st[0] ** 2 * E(st[0]))
full = of["sol"].sol(tq)[1]
from common import t_of_a_precise
wkb = phi_of_a(ow, aa, t_of=t_of_a_precise)
env = ow["A"] * (ow["a_s"] / aa) ** 1.5
d = np.max(np.abs(full - wkb) / env)
report("WKB (switch at m/H=300) vs full ODE, m=1e-25 eV, z=1300..900", d < 1e-2,
       f"max |diff|/envelope = {d:.1e}")

# 5. adiabatic invariant: rho a^3 constant after onset (m = 1e-29 eV, full ODE to a=1)
o5 = evolve(1e-29 / H0_EV, switch=1e12)
r_end = rho_today_over_phii2(o5)
o5b = evolve(1e-29 / H0_EV, a_end=0.5, switch=1e12)
r_half = rho_today_over_phii2(o5b) * 0.5**3
report("rho a^3 conserved after oscillation onset (m=1e-29 eV, a=0.5 vs 1)",
       abs(r_end / r_half - 1) < 2e-2, f"ratio {r_end / r_half:.4f}")

# 6. recombination: optical depth 1 at Planck z*
z, g, kap, x = visibility()
zk = np.interp(1.0, kap[::-1], z[::-1])
zs = inp("z_star")
report("Peebles+fudge: kappa(z)=1 vs Planck 2018 z*", abs(zk / zs - 1) < 3e-3, f"{zk:.1f} vs {zs:.2f}")
zp = z[np.argmax(g)]
report("visibility peak within 2% of z*", abs(zp / zs - 1) < 2e-2, f"peak z = {zp:.1f}")
fwhm = np.ptp(z[g > g.max() / 2])
report("visibility FWHM in z plausible (150-250)", 150 < fwhm < 250, f"FWHM = {fwhm:.0f}")
