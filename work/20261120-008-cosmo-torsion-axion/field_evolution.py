"""Homogeneous pseudoscalar in LCDM: phi'' + 2 (a'/a) phi' + a^2 m^2 phi = 0 (conformal time).

Spectator field (its energy density is not fed back into H). Integration in conformal
time (units 1/H0) until m/H = SWITCH, then the adiabatic (WKB) solution
phi = A (a_s/a)^{3/2} cos(m (t - t_s) + delta) in cosmic time.

beta = -(g/2) (phi_0 - phibar_LSS) with phibar_LSS = int g_vis(z) phi(z) dz (CONVENTIONS.md).
"""
import numpy as np
from scipy.integrate import solve_ivp

from common import (E, OMR, H0_EV, tau_of_a, t_of_a, t_of_a_precise, TAU_TODAY, RHO_CRIT0)
from recombination import visibility

SWITCH = 300.0
A_INIT = 1e-9


def evolve(mu, a_end=1.0, switch=SWITCH, Efun=E, a_init=A_INIT, rtol=1e-10):
    """Integrate with phi_i = 1. Returns dict with solution pieces.

    mu = m/H0. Efun(a) = H/H0. State y = (a, phi, dphi/dtau).
    """
    omr = OMR if Efun is E else Efun(1.0) ** 2   # only used for initial tau in radiation era
    tau_i = a_init / np.sqrt(omr) if Efun is E else a_init
    # slow-roll start: in radiation era phi = 1 - (mu t)^2/10 ... start well before m ~ H
    y0 = [a_init, 1.0, 0.0]

    def rhs(tau, y):
        a, p, dp = y
        aH = a * Efun(a)               # a'/a = a H  (H in H0 units)
        return [a * aH, dp, -2 * aH * dp - a * a * mu * mu * p]

    def ev_switch(tau, y):
        return mu - switch * Efun(y[0])
    ev_switch.terminal = True
    ev_switch.direction = 1

    def ev_end(tau, y):
        return y[0] - a_end
    ev_end.terminal = True

    tau_max = 10 * TAU_TODAY if Efun is E else 10 * a_end
    sol = solve_ivp(rhs, (tau_i, tau_max), y0, method="DOP853", rtol=rtol, atol=1e-14,
                    dense_output=True, events=[ev_switch, ev_end])
    out = {"sol": sol, "mu": mu, "switched": len(sol.t_events[0]) > 0}
    if out["switched"]:
        ts, ys = sol.t_events[0][0], sol.y_events[0][0]
        a_s, p, dp = ys
        H = Efun(a_s)
        pdot = dp / a_s                       # cosmic-time derivative
        Ac = p
        As = -(pdot + 1.5 * H * p) / mu
        out.update(a_s=a_s, A=np.hypot(Ac, As), delta=np.arctan2(As, Ac),
                   t_s=float(t_of_a_precise(a_s)[0]) if Efun is E else None)
    return out


def phi_of_a(out, a, t_of=None):
    """phi at scale factors a (array). Uses ODE before switch, WKB after."""
    a = np.atleast_1d(a)
    sol = out["sol"]
    res = np.empty_like(a)
    if out["switched"]:
        pre = a < out["a_s"]
    else:
        pre = np.ones_like(a, dtype=bool)
    if pre.any():
        # invert a(tau) on the ODE solution
        taus = sol.t
        avals = sol.y[0]
        tq = np.interp(np.log(a[pre]), np.log(avals), taus)
        # refine with one Newton step using a' = a^2 E
        for _ in range(3):
            st = sol.sol(tq)
            tq = tq - (st[0] - a[pre]) / (st[0] ** 2 * E(st[0]))
        res[pre] = sol.sol(tq)[1]
    if (~pre).any():
        tt = (t_of if t_of is not None else t_of_a)(a[~pre])
        res[~pre] = out["A"] * (out["a_s"] / a[~pre]) ** 1.5 * np.cos(
            out["mu"] * (tt - out["t_s"]) + out["delta"])
    return res


def envelope(out, a):
    """Oscillation amplitude at a (WKB) or |phi| if not oscillating."""
    if out["switched"] and a >= out["a_s"]:
        return out["A"] * (out["a_s"] / a) ** 1.5
    return abs(phi_of_a(out, np.array([a]))[0])


_VIS = None


def vis():
    global _VIS
    if _VIS is None:
        z, g, kap, x = visibility()
        _VIS = (z[::-1].copy(), g[::-1].copy())   # ascending z
    return _VIS


def lss_average(out, npts=None):
    """Visibility-weighted phi over recombination, with a phase-accurate t(a)."""
    z, g = vis()
    mu = out["mu"]
    if out["switched"] and out["a_s"] < 1 / (1 + z[0]):
        # number of oscillations across recombination ~ mu * Delta t
        dt = t_of_a(1 / (1 + z[0])) - t_of_a(1 / (1 + z[-1]))
        n = int(max(20001, 40 * mu * dt / (2 * np.pi)))
        zz = np.linspace(z[0], z[-1], n)
        gg = np.interp(zz, z, g)
        a = 1 / (1 + zz)
        # phase-accurate t: precise anchors + interpolation of the dense table difference
        anc = np.linspace(a.min(), a.max(), 200)
        t_anc = t_of_a_precise(anc)
        corr = t_anc - t_of_a(anc)
        tt = t_of_a(a) + np.interp(a, anc, corr)
        ph = out["A"] * (out["a_s"] / a) ** 1.5 * np.cos(mu * (tt - out["t_s"]) + out["delta"])
        return np.trapezoid(gg * ph, zz) / np.trapezoid(gg, zz)
    ph = phi_of_a(out, 1 / (1 + z))
    return np.trapezoid(g * ph, z) / np.trapezoid(g, z)


def rho_today_over_phii2(out):
    """rho_phi(a=1) / phi_i^2 in eV^2 (phi_i in eV)."""
    m = out["mu"] * H0_EV
    if out["switched"]:
        return 0.5 * m * m * (out["A"] * out["a_s"] ** 1.5) ** 2
    st = out["sol"].y[:, -1]
    a, p, dp = st
    pdot = dp / a * H0_EV
    return 0.5 * pdot**2 + 0.5 * m * m * p * p


def scan(masses_ev):
    rows = []
    _, _ = vis()
    for m in masses_ev:
        mu = m / H0_EV
        out = evolve(mu)
        lss = lss_average(out)
        z, g = vis()
        zpk = z[np.argmax(g)]
        env_lss = envelope(out, 1 / (1 + zpk))
        if out["switched"]:
            phi0 = np.nan
            env0 = out["A"] * out["a_s"] ** 1.5
        else:
            phi0 = out["sol"].y[1, -1]
            env0 = abs(phi0)
        rows.append(dict(m=m, phibar_lss=lss, env_lss=env_lss, phi0=phi0, env0=env0,
                         rho0=rho_today_over_phii2(out), switched=out["switched"],
                         a_s=out.get("a_s", np.nan)))
    return rows
