"""KK-02: initial-angle tuning of the Wilson-line axion (beta = k alpha Delta theta / 2 pi).

Potential: with a charged bulk fermion at M l ~ 220-280 (task 009) the Hosotani sum is dominated by n = 1,
higher harmonics are suppressed by exp(-M l) ~ 1e-96, so V = m^2 f^2 (1 - cos theta) about the minimum.
A field released from rest at theta_i (|theta_i| <= pi, uniform prior on the circle) rolls to the nearest
minimum, so |Delta theta| <= |theta_i| <= pi. Best case (full relaxation before today, frozen at z*):
|Delta theta| = |theta_i|. Required: |theta_i| >= theta_req = 2 pi beta/(k alpha);
P(uniform) = 1 - theta_req/pi; distance from the hilltop d = pi - theta_req.
Second check and overshoot: integrate theta'' + 3 H theta' + m^2 sin theta = 0 on a theta_i grid for m = 1-200 H0.
If the field is still in its first swing today (m ~ 5-30 H0, start near the hilltop), it can be caught on the far
side, so |Delta theta| = |theta_i| + |theta_0| > pi is possible. The joint probability uses theta_i uniform and
ln m uniform on [H0, H*] (Hosotani m ~ exp(-M l/2) with M l uniform); above 200 H0 full relaxation is used.
"""
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

from common import ACT, ALPHA, A_STAR, DEG, E, ESK, H0, H_STAR, check, done, info


def theta_req(beta_deg, k):
    return 2 * math.pi * beta_deg * DEG / (k * ALPHA)


def p_reach(beta_deg, k):
    """Fraction of uniform initial angles with |Delta theta| >= theta_req (best case, full relaxation)."""
    return max(0.0, 1 - theta_req(beta_deg, k) / math.pi)


def p_interval(lo, hi, k):
    return max(0.0, p_reach(max(lo, 0), k) - p_reach(hi, k))


def delta_theta(theta_i, mu):
    """Delta theta between z* and today for m = mu H0, released from rest at a = 1e-6 (frozen until H ~ m)."""
    def rhs(N, y):
        a = math.exp(N)
        e = E(a)
        return [y[1] / e, -3 * y[1] - mu**2 * math.sin(y[0]) / e]

    sol = solve_ivp(rhs, [math.log(1e-6), 0.0], [theta_i, 0.0], method="DOP853", rtol=1e-9, atol=1e-12,
                    dense_output=True)
    return sol.sol(0.0)[0] - sol.sol(math.log(A_STAR))[0]


MU_NUM_MAX = 200.0


def delta_theta_vec(th0, mu):
    """Vectorised delta_theta: all initial angles in one ODE system (same physics as delta_theta)."""
    n = len(th0)

    def rhs(N, y):
        a = math.exp(N)
        e = E(a)
        return np.concatenate([y[n:] / e, -3 * y[n:] - mu**2 * np.sin(y[:n]) / e])

    sol = solve_ivp(rhs, [math.log(1e-6), 0.0], np.concatenate([th0, np.zeros(n)]), method="DOP853",
                    rtol=1e-9, atol=1e-12, dense_output=True)
    return sol.sol(0.0)[:n] - sol.sol(math.log(A_STAR))[:n]


def main(nu=60, nlin=400, nmu=40):
    print("KK-02 Wilson-line axion: initial-angle tuning (candidate mechanism, ALP class)")
    bmax1 = ALPHA / 2 / DEG
    info(f"beta_max(k=1) = alpha/2 = {bmax1:.4f} deg")
    rows = []
    for lab, (b, s) in (("ACT DR6 2509.13654", ACT), ("ACT+Planck 2608.06480", ESK)):
        for k in (1, 2, 3):
            tc = theta_req(b, k) / math.pi
            tlo, thi = theta_req(b - s, k) / math.pi, theta_req(b + s, k) / math.pi
            pc = p_reach(b, k)
            p1 = p_interval(b - s, b + s, k)
            d_lo = 1 - tlo
            rows.append((lab, k, tc, tlo, thi, pc, p1, d_lo))
            info(f"{lab} beta = {b} +- {s} deg, k = {k}: theta_req = {tc:.3f} pi (1 sigma {tlo:.3f}-{thi:.3f} pi); "
                 f"P(beta >= central) = {pc:.3f}; P(beta in 1 sigma) = {p1:.3f}; "
                 f"1 sigma lower edge needs start within {d_lo:.3f} pi of the hilltop")
        kmin = math.ceil(2 * b * DEG / ALPHA)
        info(f"{lab}: smallest integer k reaching the central value = {kmin} (k in 6Z: k = 6, theta_req = "
             f"{theta_req(b, 6) / math.pi:.3f} pi)")
        nsig = (b - bmax1) / s
        info(f"{lab}: k = 1 ceiling {bmax1:.3f} deg sits {nsig:.2f} sigma below the central value")

    # checks on the closed forms
    check("H1 full relaxation: k=1 cannot reach the ACT central value (theta_req = 1.028 pi > pi)",
          abs(theta_req(ACT[0], 1) / math.pi - 1.0284) < 2e-3 and p_reach(ACT[0], 1) == 0.0)
    check("H2 full relaxation, k=1, ACT 1 sigma: P(beta in [0.141, 0.289] deg) = 1 - 0.673 = 0.327 (start within 0.33 pi of hilltop)",
          abs(p_interval(ACT[0] - ACT[1], ACT[0] + ACT[1], 1) - 0.327) < 2e-3)
    check("H3 full relaxation: k=1 cannot reach 2608.06480 even at its 1 sigma lower edge (0.220 deg > 0.209 deg)",
          p_reach(ESK[0] - ESK[1], 1) == 0.0)
    check("H4 2608.06480 central value needs k >= 2 (theta_req(k=2) = 0.663 pi)",
          math.ceil(2 * ESK[0] * DEG / ALPHA) == 2 and abs(theta_req(ESK[0], 2) / math.pi - 0.6625) < 2e-3)

    # independent numerical check and the overshoot channel: many trajectories integrated together
    u = np.concatenate([np.logspace(-9, -2.3, nu), np.linspace(0.005, 0.995, nlin)])  # u = 1 - theta_i/pi
    th = np.sort(math.pi * (1 - u))
    wts = np.gradient(th) / math.pi           # measure of each grid cell (uniform prior on (0, pi))
    targets = {"ACT>=central": ACT[0], "ACT>=lo1s": ACT[0] - ACT[1], "ESK>=central": ESK[0], "ESK>=lo1s": ESK[0] - ESK[1]}
    mus = np.logspace(0, math.log10(MU_NUM_MAX), nmu)
    tt = np.array([0.5, 0.9, 0.99]) * math.pi
    dv = delta_theta_vec(tt, 10.0)
    ds = np.array([delta_theta(t, 10.0) for t in tt])
    check("H0 vectorised and scalar integrators agree (m = 10 H0)", np.allclose(dv, ds, atol=1e-5))
    res = {}
    for mu in mus:
        dth = np.abs(delta_theta_vec(th, mu))
        res[mu] = dth
    # full relaxation check at the top of the numerical range
    dtop = res[mus[-1]]
    bulk = th < 0.9 * math.pi
    ratio = dtop[bulk] / th[bulk]
    info(f"m = {mus[-1]:.0f} H0: |Delta theta|/theta_i for theta_i < 0.9 pi in [{ratio.min():.4f}, {ratio.max():.4f}]")
    check(f"H5 m = {mus[-1]:.0f} H0 is fully relaxed away from the hilltop (|Delta theta| = theta_i to 3%)",
          abs(ratio - 1).max() < 0.03)
    pn = np.sum(wts * (dtop >= theta_req(ACT[0] - ACT[1], 1)))
    check(f"H6 m = {mus[-1]:.0f} H0: numerical P(k=1, beta >= ACT lower 1 sigma) = {pn:.3f} vs full-relaxation "
          f"{p_reach(ACT[0] - ACT[1], 1):.3f}", abs(pn - p_reach(ACT[0] - ACT[1], 1)) < 0.02)
    # overshoot: field caught on the far side of its first swing today -> |Delta theta| > theta_i
    mx = {mu: res[mu].max() / math.pi for mu in mus}
    best = max(mx, key=mx.get)
    info("max |Delta theta|/pi over theta_i, per m/H0: " + ", ".join(f"{mu:.1f}:{mx[mu]:.2f}" for mu in mus))
    check("H7 overshoot exists: for some m in the H0-H* window |Delta theta| > pi (k=1 can exceed 0.209 deg)",
          mx[best] > 1.0, f"max {mx[best]:.2f} pi at m = {best:.1f} H0")
    # joint probability: theta_i uniform, ln m uniform on [H0, H*] (x = M l uniform in the Hosotani mass)
    lnw = math.log(H_STAR / H0)
    out = {}
    for k in (1, 2):
        for lab, b in targets.items():
            tr = theta_req(b, k)
            pm = np.array([np.sum(wts * (res[mu] >= tr)) for mu in mus])
            num = np.trapezoid(pm, np.log(mus))
            ana = p_reach(b, k) * (lnw - math.log(MU_NUM_MAX))   # above MU_NUM_MAX: full relaxation
            pj = (num + ana) / lnw
            tail = max(pm[-1] - p_reach(b, k), 0) * (lnw - math.log(MU_NUM_MAX)) / lnw  # overshoot tail above 200 H0, upper bound
            okm = mus[pm > 0]
            rng = f"m = {okm.min():.1f}-{okm.max():.1f} H0 numerically" if len(okm) else "no m in grid"
            out[(k, lab)] = pj
            info(f"k = {k}, {lab} ({b:.3f} deg, theta_req = {tr / math.pi:.3f} pi): joint P = {pj:.3e}; "
                 f"(+ <= {tail:.1e} overshoot tail above {MU_NUM_MAX:.0f} H0); P(m) max = {pm.max():.3e}; reachable at {rng}")
            out[(k, lab, "tail")] = tail
    check("H8 k=1 joint probability of beta >= ACT central is small but nonzero (overshoot only)",
          0 < out[(1, "ACT>=central")] < 0.05)
    check("H9 k=1 joint probability of beta >= 2608.06480 central is below that for ACT central",
          out[(1, "ESK>=central")] < out[(1, "ACT>=central")])
    return rows, out


if __name__ == "__main__":
    import time
    t = time.time()
    _, o1 = main()
    print(f"[INFO] base grid took {time.time() - t:.0f} s; convergence run with 1.75x theta grid and 1.5x m grid")
    if "--noconv" not in sys.argv:
        _, o2 = main(nu=90, nlin=700, nmu=60)
        r = o2[(1, "ACT>=central")] / o1[(1, "ACT>=central")]
        check("H10 convergence: k=1 joint P(beta >= ACT central) changes < 10% on the finer grid", abs(r - 1) < 0.10,
              f"{o1[(1, 'ACT>=central')]:.2e} -> {o2[(1, 'ACT>=central')]:.2e}")
    sys.exit(done())
