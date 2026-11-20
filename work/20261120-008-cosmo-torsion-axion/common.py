"""Shared inputs and LCDM background for task 20261120-008 (torsion axion, EC-03).

Units: natural units, energies in eV. In the field ODE, conformal time tau and
cosmic time t are in units of 1/H0 and the mass is mu = m/H0.
All measured numbers are loaded from data/ (no hardcoded measurements).
"""
import csv
import os

import numpy as np
from scipy.integrate import solve_ivp, cumulative_simpson, quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data")


def _load(name, key="id"):
    with open(os.path.join(DATA, name), newline="", encoding="utf-8-sig") as fh:
        return {row[key]: row for row in csv.DictReader(fh)}


INP = _load("cosmo_torsion_axion_inputs.csv")
BINP = _load("bounds_inputs.csv")
BOUNDS = _load("bounds.csv")
MEAS = _load("birefringence_measurements.csv")


def inp(k):
    return float(INP[k]["value"])


def binp(k):
    return float(BINP[k]["value"])


# ---- constants ---------------------------------------------------------------
HBAR = binp("hbar")            # eV s
HBARC = binp("hbarc")          # eV m
MPC = binp("Mpc")              # m
C_KMS = binp("c_kms")
MP = binp("MP_reduced")        # eV
GEV = binp("GeV")
CM3_TO_EV3 = (HBARC * 100.0) ** 3   # 1 cm^-3 -> eV^3 multiply by (hbar c in eV cm)^3
DEG = np.pi / 180.0

h = inp("H0") / 100.0
H0_EV = inp("H0") / (MPC / 1e3) * HBAR      # km/s/Mpc -> s^-1 -> eV
OMB = inp("ombh2") / h**2
OMC = inp("omch2") / h**2
OMM = OMB + OMC
T0_K = inp("T_CMB")
KB = inp("k_B")
T0_EV = KB * T0_K
NEFF = inp("N_eff")
# photon density: rho_g = (pi^2/15) T^4; neutrinos 7/8 (4/11)^{4/3} N_eff
RHO_CRIT0 = 3.0 * H0_EV**2 * MP**2          # eV^4
OMG = (np.pi**2 / 15.0) * T0_EV**4 / RHO_CRIT0
OMR = OMG * (1.0 + 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0) * NEFF)
OML = 1.0 - OMM - OMR
RHO_DM0 = OMC * RHO_CRIT0                   # eV^4
ALPHA = inp("alpha_em")

BETA_DEG = float(MEAS["DPK25-beta-baseline"]["value"])
BETA_ERR_DEG = float(MEAS["DPK25-beta-baseline"]["err_plus"])
ACB_SPT = float(MEAS["BAL25-ACB-SPT3G"]["value"])


def gev_inv_to_ev_inv(x):
    return x / GEV


def bound_upper(id_):
    """Numeric part of value_eV for a one-sided upper bound like '< 8e-22'."""
    s = BOUNDS[id_]["value_eV"].replace("<~", "").replace("<", "").strip().split()[0]
    return float(s)


def E(a):
    """H(a)/H0 for flat LCDM with radiation."""
    return np.sqrt(OMR / a**4 + OMM / a**3 + OML)


def H_ev(a):
    return H0_EV * E(a)


# ---- conformal and cosmic time tables ----------------------------------------
_LNA = np.linspace(np.log(1e-12), 0.0, 400001)
_A = np.exp(_LNA)
# dtau = da/(a^2 H) = dlna/(a H); dt = dlna/H  (units 1/H0)
_dtau = 1.0 / (_A * E(_A))
_dt = 1.0 / E(_A)
_tau0 = _A[0] / np.sqrt(OMR)          # radiation era: tau = a/sqrt(OMR)
_t0 = _A[0] ** 2 / (2 * np.sqrt(OMR))
TAU_TAB = _tau0 + cumulative_simpson(_dtau, x=_LNA, initial=0.0)
T_TAB = _t0 + cumulative_simpson(_dt, x=_LNA, initial=0.0)


def tau_of_a(a):
    return np.interp(np.log(a), _LNA, TAU_TAB)


def t_of_a(a):
    return np.interp(np.log(a), _LNA, T_TAB)


def a_of_tau(tau):
    return np.exp(np.interp(tau, TAU_TAB, _LNA))


TAU_TODAY = TAU_TAB[-1]
T_TODAY = T_TAB[-1]


def t_of_a_precise(a):
    """Cosmic time (1/H0) by direct quadrature, for phase-sensitive uses."""
    a = np.atleast_1d(a)
    return np.array([_t0 + quad(lambda x: 1.0 / (x * E(x)), _A[0], ai, epsabs=0, epsrel=1e-13,
                                limit=500)[0] for ai in a])
