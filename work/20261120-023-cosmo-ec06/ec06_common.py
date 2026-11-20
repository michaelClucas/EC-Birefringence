"""Shared pieces for task 20261120-023 (EC-06: homogeneous O(H) torsion background + CFJ photon coupling).

Reuses the flat LCDM background (photons + massless neutrinos) of task 008 (work/20261120-008-cosmo-torsion-axion/common.py).
All measured numbers come from data/.

Bookkeeping (CONVENTIONS.md): CFJ term -1/2 p_mu A_nu Ftilde^{mu nu}, beta = +1/2 int p_tau dtau (conformal time,
last scattering -> today). The orthonormal (hat) component with e^0hat = a dtau is p_0hat = p_tau / a, and dt = a dtau,
so beta = +1/2 int p_0hat dt = +1/2 int_0^{z_s} p_0hat(z) dz / ((1+z) H(z)).
Background shapes: p_0hat(z) = xi*eps * H(z) * (1+z)^n   (n = 0: tracks H)
                   p_0hat     = xi*eps * H0              (constant, k_AF-like)
Here xi*eps is one dimensionless number (p = xi S or xi T with S_0hat or T_0hat = eps H (1+z)^n).
"""
import csv
import os
import sys

import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data")
W008 = os.path.join(ROOT, "work", "20261120-008-cosmo-torsion-axion")
for p in (W008, ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

import common as c008  # noqa: E402  (task 008 background)
from common import E, H0_EV, OMB, OMM, OMR, OML, h  # noqa: E402,F401


def _load(name):
    with open(os.path.join(DATA, name), newline="", encoding="utf-8-sig") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


EC6 = _load("ec06_cosmo_inputs.csv")
COSMO = _load("cosmo_inputs.csv")
MEAS = _load("birefringence_measurements.csv")
PL18 = _load("planck2018_bestfit.csv")


def ec6(k):
    return float(EC6[k]["value"])


ZSTAR = float(COSMO["z_star"]["value"])
ZRE = ec6("z_re")
DZRE = ec6("dz_re")
ZHE2 = ec6("z_HeII")
DEG = np.pi / 180.0
BETA_ACT = float(MEAS["DPK25-beta-baseline"]["value"])

SHAPES = [("H", 0.0), ("H(1+z)^-1", -1.0), ("H(1+z)^-0.5", -0.5), ("H(1+z)^+0.5", 0.5), ("H(1+z)^+1", 1.0),
          ("const", None)]


def p_over_xieps(z, n):
    """p_0hat / (xi eps) in units of H0."""
    z = np.asarray(z, float)
    if n is None:
        return np.ones_like(z)
    return E(1.0 / (1.0 + z)) * (1.0 + z) ** n


def beta_unit(zs, n):
    """beta(z_s) in rad for xi*eps = 1: 1/2 int_0^zs p_0hat dz/((1+z)H) (numerical, cosmic-time route)."""
    f = lambda z: 0.5 * p_over_xieps(z, n) / ((1 + z) * E(1 / (1 + z)))
    return quad(f, 0.0, zs, epsabs=0, epsrel=1e-11, limit=400)[0]


def beta_unit_analytic(zs, n):
    """Analytic: n = 0 -> 1/2 ln(1+zs); n != 0 -> ((1+zs)^n - 1)/(2n); const -> 1/2 H0 t_L (t_L numerical, 1/H0)."""
    if n is None:
        return None
    if n == 0:
        return 0.5 * np.log1p(zs)
    return ((1 + zs) ** n - 1) / (2 * n)


def beta_unit_conformal(zs, n):
    """Independent route: beta = 1/2 int p_tau dtau with p_tau = a p_0hat, on the task-008 conformal-time table."""
    a_s = 1.0 / (1.0 + zs)
    lna = c008._LNA
    m = lna >= np.log(a_s)
    a = np.exp(lna[m])
    tau = c008.TAU_TAB[m]
    ptau = a * p_over_xieps(1 / a - 1, n)
    return 0.5 * np.trapezoid(ptau, tau)


def beta_curve(z, n):
    """beta(z_s) for an array of z_s (cumulative, rad, xi eps = 1)."""
    z = np.asarray(z, float)
    out = np.empty_like(z)
    prev, acc = 0.0, 0.0
    order = np.argsort(z)
    for i in order:
        acc += beta_unit(z[i], n) - beta_unit(prev, n) if z[i] > prev else 0.0
        prev = z[i]
        out[i] = acc
    return out


# ---------------- reionization (tanh, Planck/CAMB convention) -----------------------------------------------
INP008 = c008.INP
YP = float(INP008["Y_P"]["value"])
SIGT = float(INP008["sigma_T"]["value"])          # m^2
MH_EV = float(INP008["m_H"]["value"])
HBAR, MPC = c008.HBAR, c008.MPC
C_SI = 2.99792458e8


def xe_reion(z):
    """Free electrons per hydrogen (tanh in y = (1+z)^1.5), He I with H, He II at z = 3.5 (CAMB convention)."""
    fHe = YP / (4.0 * (1 - YP))  # n_He/n_H with m_He ~ 4 m_H (0.7% approximation)
    y, yre = (1 + z) ** 1.5, (1 + ZRE) ** 1.5
    dy = 1.5 * np.sqrt(1 + ZRE) * DZRE
    xe = (1 + fHe) / 2 * (1 + np.tanh((yre - y) / dy))
    xe += fHe / 2 * (1 + np.tanh((ZHE2 - z) / DZRE))
    return xe


def n_H0_m3():
    rho_b = OMB * c008.RHO_CRIT0  # eV^4
    nH_ev3 = (1 - YP) * rho_b / MH_EV
    return nH_ev3 / (c008.HBARC ** 3)  # eV^3 -> m^-3


def dtau_dz(z):
    """Thomson optical depth per unit z: n_e sigma_T c / ((1+z) H)."""
    H_si = H0_EV / HBAR * E(1 / (1 + z))
    return xe_reion(z) * n_H0_m3() * (1 + z) ** 3 * SIGT * C_SI / ((1 + z) * H_si)


def reion_visibility(zmax=40.0, n=8001):
    z = np.linspace(0.0, zmax, n)
    dt = dtau_dz(z)
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (dt[1:] + dt[:-1]) * np.diff(z))])
    g = dt * np.exp(-tau)
    return z, g, tau
