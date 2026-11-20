"""Item 5: anisotropic birefringence from inflationary fluctuations of a light spectator phi.

delta phi has a scale-invariant spectrum P_dphi = (H_inf/2 pi)^2 per ln k. Superhorizon modes
obey the same linear equation as the background, so delta phi_k(tau*) = (phibar_LSS/phi_i) x
delta phi_k,i; inside the horizon they decay (transfer T(k), massless approximation, computed
here). delta beta(n) = +(g/2) delta phi(tau*, D n) (observer term is a monopole).
  L(L+1)C_L/2pi = (g/2)^2 (H_inf/2pi)^2 (phibar_LSS/phi_i)^2 S_L,
  S_L = 2 L(L+1) int dk/k T(k)^2 j_L(k D)^2  (S_L -> 1 for T = 1).
With g fixed by beta: g = 2 beta / (F phi_i), F = |phi0 - phibar_LSS|/phi_i:
  A_CB = (beta * R / F)^2 (H_inf/(2 pi phi_i))^2 S_L,  R = phibar_LSS/phi_i.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import spherical_jn

from common import (E, tau_of_a, TAU_TODAY, inp, ACB_SPT, BETA_DEG, DEG, MP, GEV)
from variants import load_scan, required_g

BETA = BETA_DEG * DEG
Z_STAR = inp("z_star")
TAU_STAR = float(tau_of_a(1 / (1 + Z_STAR)))
D = TAU_TODAY - TAU_STAR


def transfer(k):
    """Massless mode delta'' + 2 aH delta' + k^2 delta = 0 from a tiny tau, delta_i = 1."""
    from common import a_of_tau
    ti = 1e-6 * TAU_STAR

    def rhs(t, y):
        a = a_of_tau(t)
        return [y[1], -2 * a * E(a) * y[1] - k * k * y[0]]
    s = solve_ivp(rhs, (ti, TAU_STAR), [1.0, 0.0], method="DOP853", rtol=1e-8, atol=1e-12)
    return s.y[0, -1]


def S_L(Ls, T_of_k=None, nk=40000):
    x = np.logspace(-2, np.log10(4000), nk)          # x = k D
    T = np.ones_like(x) if T_of_k is None else T_of_k(x / D)
    out = []
    for L in Ls:
        jl = spherical_jn(L, x)
        out.append(2 * L * (L + 1) * np.trapezoid(T**2 * jl**2 / x, x))
    return np.array(out)


def tabulated_transfer(nk=300):
    kk = np.logspace(np.log10(0.05 / D), np.log10(5000 / D), nk)
    TT = np.array([transfer(k) for k in kk])
    return lambda k: np.interp(np.log(k), np.log(kk), TT, left=1.0, right=0.0)


def main():
    H_MAX = np.pi * MP * np.sqrt(inp("r_BK18") * np.exp(inp("logA")) * 1e-10 / 2)
    print(f"tau* = {TAU_STAR:.4f}/H0, D = {D:.4f}/H0, H_inf,max (BK18) = {H_MAX/GEV:.2e} GeV")
    Tk = tabulated_transfer()
    Ls = [2, 10, 30, 50, 100, 200, 400]
    SL = S_L(Ls, Tk)
    print("S_L (horizon-entry suppression, massless):", ", ".join(f"L={L}: {s:.3f}" for L, s in zip(Ls, SL)))
    sc = load_scan()
    m, F, cap, phiE, greq, _ = required_g(sc)
    R = sc["phibar_lss"]
    print(f"SPT-3G A_CB < {ACB_SPT:.1e} (95%)")
    print(f"{'m [eV]':>10} {'R':>7} {'F':>7} {'(H/phi_i)_max (SL=1)':>21} {'phi_i [GeV]':>12} {'A_CB at H_max':>14}")
    for i in range(0, len(m), 4):
        if m[i] > 1.5e-27:
            continue
        rmax = 2 * np.pi * np.sqrt(ACB_SPT) * F[i] / (BETA * abs(R[i]))
        acb = (BETA * R[i] / F[i]) ** 2 * (H_MAX / (2 * np.pi * phiE[i])) ** 2
        print(f"{m[i]:10.3e} {R[i]:7.3f} {F[i]:7.3f} {rmax:21.3e} {phiE[i]/GEV:12.3e} {acb:14.3e}")
    print("A_CB(H_inf, phi_i) = (beta R/F)^2 (H_inf/(2 pi phi_i))^2 S_L; plateau case R = F = 1:"
          f" A_CB = {(BETA/(2*np.pi))**2:.3e} (H_inf/phi_i)^2")


if __name__ == "__main__":
    main()
