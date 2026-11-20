"""Hydrogen recombination (Peebles 3-level atom with the RECFAST fudge) and the
visibility function. Helium and reionization are neglected (see README).

visibility() returns z grid, g(z) normalised so that int g dz = 1 over [ZMIN, ZMAX].
SI units internally.
"""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

from common import inp, OMB, h, T0_K, H_ev, HBAR, HBARC

KB = inp("k_B")
ME = inp("m_e")
MH = inp("m_H")
EI = inp("E_ion_H")
SIGT = inp("sigma_T")
LAM2S = inp("Lambda_2s1s")
LYA = inp("lambda_Lya")
FUD = inp("recfast_F")
PA, PB, PC, PD = (inp(k) for k in ("pequignot_a", "pequignot_b", "pequignot_c", "pequignot_d"))
YP = inp("Y_P")
C_SI = 2.99792458e8
GN = inp("G_Newton")
EV_KG = inp("m_proton_kg")
H0_SI = h * 100.0 * 1e3 / 3.0856775814913673e22
RHO_CRIT_SI = 3 * H0_SI**2 / (8 * np.pi * GN)                    # kg/m^3
NH0 = (1 - YP) * OMB * RHO_CRIT_SI / (MH * EV_KG)                 # m^-3 today

ZMIN, ZMAX = 200.0, 2000.0


def alpha_B(T):
    t = T / 1e4
    return FUD * 1e-19 * PA * t**PB / (1 + PC * t**PD)        # m^3/s


def saha_x(z):
    T = T0_K * (1 + z)
    kT = KB * T
    nh = NH0 * (1 + z) ** 3
    # (m_e kT / 2 pi hbar^2)^{3/2} in m^-3: (m_e kT)^{3/2}/(2 pi)^{3/2}/(hbar c)^3
    S = (ME * kT / (2 * np.pi)) ** 1.5 / HBARC**3 * np.exp(-EI / kT) / nh
    return (-S + np.sqrt(S * S + 4 * S)) / 2


def H_si(z):
    return H_ev(1.0 / (1 + z)) / HBAR


def _rhs(z, y):
    x = min(max(y[0], 0.0), 1.0)
    T = T0_K * (1 + z)
    kT = KB * T
    nh = NH0 * (1 + z) ** 3
    aB = alpha_B(T)
    bB = aB * (ME * kT / (2 * np.pi)) ** 1.5 / HBARC**3 * np.exp(-EI / (4 * kT))
    H = H_si(z)
    K = LYA**3 / (8 * np.pi * H)
    C = (1 + K * LAM2S * nh * (1 - x)) / (1 + K * (LAM2S + bB) * nh * (1 - x))
    dxdt = -C * (aB * nh * x * x - bB * (1 - x) * np.exp(-0.75 * EI / kT))
    return [-dxdt / ((1 + z) * H)]


def ionization(n=20001):
    zs = 1700.0
    sol = solve_ivp(_rhs, (zs, ZMIN), [saha_x(zs)], method="LSODA", rtol=1e-8, atol=1e-12,
                    dense_output=True)
    z = np.linspace(ZMAX, ZMIN, n)
    x = np.where(z >= zs, np.minimum(saha_x(z), 1.0), sol.sol(np.minimum(z, zs))[0])
    return z, x


def visibility(n=20001):
    """z (descending from ZMAX to ZMIN), g(z) with int g dz = 1, optical depth kappa(z)."""
    z, x = ionization(n)
    ne = x * NH0 * (1 + z) ** 3
    dk_dz = ne * SIGT * C_SI / ((1 + z) * H_si(z))
    # kappa(z) = int_0^z ; z < ZMIN contributes ~ x_e(200) small; we add it approximately
    zz = z[::-1]
    kap_asc = cumulative_trapezoid(dk_dz[::-1], zz, initial=0.0)
    kap_low = dk_dz[-1] * ZMIN / 3.0          # rough tail below ZMIN (x_e frozen, ~(1+z)^2)
    kap = (kap_asc + kap_low)[::-1]
    g = dk_dz * np.exp(-kap)
    norm = np.trapezoid(g[::-1], zz)
    return z, g / norm, kap, x


if __name__ == "__main__":
    z, g, kap, x = visibility()
    print(f"peak of visibility at z = {z[np.argmax(g)]:.2f}")
    i = np.argmin(abs(kap - 1))
    print(f"kappa = 1 at z = {z[i]:.2f}; x_e(z=1000) = {np.interp(1000, z[::-1], x[::-1]):.4f}")
