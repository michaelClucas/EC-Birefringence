"""Checks R1-R3 (rotation from L = -1/4 Z (F^2 + c_odd F Ftilde)) and H1-H3 (Hosotani mass).

R: helicity modes on conformally flat FRW obey (Z A')' + Z k^2 A - h c_odd k Z' A = 0
(h = +-1; derived in derivation.md Sec. 4). The accumulated rotation, half the phase
difference of the two helicities, must equal (c_odd/2) Delta ln Z for any Z(tau).
H: V(theta) for a 5D Dirac fermion of charge q and mass M on a circle of circumference l,
V = (3/(pi^2 l^4)) sum_n (1 + n x + n^2 x^2/3) e^{-n x} cos(n q theta)/n^5, x = M l;
check against the direct KK-mode sum (Coleman-Weinberg with dimensional subtraction).
"""
import numpy as np
from scipy.integrate import solve_ivp, quad

from common import check, info


def rotation(c_odd, Zfun, dZfun, k, t0, t1):
    def rhs(t, y, h):
        A, P = y[0] + 1j * y[1], y[2] + 1j * y[3]   # P = Z A'
        Z, dZ = Zfun(t), dZfun(t)
        dA = P / Z
        dP = -Z * k**2 * A + h * c_odd * k * dZ * A
        return [dA.real, dA.imag, dP.real, dP.imag]
    out = {}
    for h in (+1, -1):
        Z0 = Zfun(t0)
        w0 = np.sqrt(k**2 - h * c_odd * k * dZfun(t0) / Z0)
        y0 = [1.0, 0.0, 0.0, -Z0 * w0]   # A = e^{-i w t}
        sol = solve_ivp(rhs, (t0, t1), y0, args=(h,), rtol=1e-10, atol=1e-12, dense_output=True)
        ts = np.linspace(t0, t1, 6001)
        A = sol.sol(ts)[0] + 1j * sol.sol(ts)[1]
        out[h] = np.unwrap(np.angle(A))[-1] - np.unwrap(np.angle(A))[0]
    return 0.5 * (out[-1] - out[+1])   # (phase_- - phase_+)/2 with phase = -int w


def main():
    k = 120.0
    cases = [(0.3, lambda t: 1 + 0.5 * np.tanh(t - 2), lambda t: 0.5 / np.cosh(t - 2) ** 2),
             (-1.2, lambda t: np.exp(-0.4 * t), lambda t: -0.4 * np.exp(-0.4 * t)),
             (2.0, lambda t: 2 + np.sin(t), lambda t: np.cos(t))]
    ok = True
    for c, Zf, dZf in cases:
        rot = rotation(c, Zf, dZf, k, 0.0, 4.0)
        pred = 0.5 * c * np.log(Zf(4.0) / Zf(0.0))
        info(f"rotation c_odd={c}: numeric {rot:.6e}, (c/2) Delta ln Z = {pred:.6e}")
        ok &= abs(rot + pred) < 3e-3 * abs(pred)
    check("R1 rotation = (c_odd/2) Delta ln Z in magnitude, one consistent sign, three Z(tau) profiles", ok)
    # 2608: Z = phi^3, c_odd = 2c -> beta = 3 c Delta ln phi; alpha ~ 1/Z -> |Delta alpha/alpha| = |beta|/|c|
    c, r = 0.7, 1.003
    beta = 0.5 * (2 * c) * np.log(r**3)
    dal = abs(r**-3 - 1)
    check("R2 2608 normalization: beta = 3c Delta ln phi and |Delta alpha/alpha| = |beta|/|c| to O(Delta)",
          abs(beta - 3 * c * np.log(r)) < 1e-15 and abs(dal - abs(beta) / c) < 1e-2 * dal)

    # H1-H3: Hosotani potential, closed form vs direct KK sum
    l, q = 1.0, 1

    def V_closed(th, x, nmax=200):
        n = np.arange(1, nmax + 1)
        return 3 / (np.pi**2 * l**4) * np.sum((1 + n * x + n**2 * x**2 / 3) * np.exp(-n * x) * np.cos(n * q * th) / n**5)

    def V_kk(th, x):
        # 4D Dirac fermion per KK level = -2 x complex scalar. Complex scalar with KK masses
        # m_k = (2 pi k + q th)/l: V_c = -int ds/(16 pi^2 s^3) e^{-s M^2} S(s), S = sum_k e^{-s m_k^2}.
        # theta-dependence: S(s,th) - S(s,0); direct KK sum for s > s*, Poisson-resummed form below.
        M = x / l
        ks = np.arange(-60, 61)
        ns = np.arange(1, 80)

        def dS(s):
            if s > 0.1 * l**2:
                m2 = ((2 * np.pi * ks + q * th) / l) ** 2
                m20 = ((2 * np.pi * ks) / l) ** 2
                return np.sum(np.exp(-s * m2)) - np.sum(np.exp(-s * m20))
            return (l / np.sqrt(4 * np.pi * s)) * 2 * np.sum(np.exp(-(ns * l) ** 2 / (4 * s)) * (np.cos(ns * q * th) - 1))

        def integrand(s):
            return 2 * np.exp(-s * M**2) * dS(s) / (16 * np.pi**2 * s**3)
        v1, _ = quad(integrand, 1e-4, 0.1 * l**2, limit=400)
        v2, _ = quad(integrand, 0.1 * l**2, np.inf, limit=400)
        return v1 + v2
    for x in (0.0, 2.0, 6.0):
        dv_closed = V_closed(1.1, x) - V_closed(0.0, x)
        dv_kk = V_kk(1.1, x)
        info(f"Hosotani x={x}: closed {dv_closed:.6e}, proper-time KK sum {dv_kk:.6e}")
        okh = abs(dv_closed - dv_kk) < 1e-4 * abs(dv_closed)
        check(f"H1 Hosotani closed form = proper-time KK sum at x = M l = {x}", okh)
    # mass: m^2 f^2 = V''(min); large x: n=1 term m^2 = q^2 g4^2 (x^2+3x+3) e^{-x}/(pi^2 l^2) (f = 1/(l g4))
    g4 = 0.3
    for x in (20.0, 60.0):
        h = 1e-3
        th0 = np.pi  # fermion potential minimum at theta = pi for q = 1
        Vpp = (V_closed(th0 + h, x) - 2 * V_closed(th0, x) + V_closed(th0 - h, x)) / h**2
        m2 = Vpp * (l * g4) ** 2
        m2a = q**2 * g4**2 * (x**2 + 3 * x + 3) * np.exp(-x) / (np.pi**2 * l**2)
        check(f"H2 m^2 = q^2 g4^2 (x^2+3x+3) e^-x/(pi^2 l^2) (n=1 term) at x = {x}", abs(m2 / m2a - 1) < 1e-4,
              f"ratio {m2 / m2a:.4f}")


if __name__ == "__main__":
    main()
