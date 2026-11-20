"""Shared inputs and cosmology for task 20261120-007 (conventions of CONVENTIONS.md).

All measured numbers come from data/ (cosmo_inputs.csv, birefringence_measurements.csv)
or from work/20261120-005-bounds/derived_bounds.csv (lab |S_0| bound).
"""
import csv
import math
import os

from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
ZETA3 = 1.2020569031595942


def _rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def load():
    c = {k: float(r["value"]) for k, r in _rows(os.path.join(ROOT, "data", "cosmo_inputs.csv")).items()}
    errs = {k: r["err"] for k, r in _rows(os.path.join(ROOT, "data", "cosmo_inputs.csv")).items()}
    c["xi_nu_err"] = float(errs["xi_nu"])
    c["N_eff_obs_err"] = float(errs["N_eff_obs"])
    b = _rows(os.path.join(ROOT, "data", "birefringence_measurements.csv"))["DPK25-beta-baseline"]
    c["beta_deg"] = float(b["value"])
    c["beta_err_deg"] = float(b["err_plus"])
    d = _rows(os.path.join(ROOT, "work", "20261120-005-bounds", "derived_bounds.csv"))
    c["S0_lab_max"] = float(d["S0-boost-neutron"]["value_eV"])
    return c


class Cosmo:
    """Flat LCDM with photons + 3 massless-equivalent neutrino species in the radiation term."""

    def __init__(self, c):
        self.c = c
        self.hbar = c["hbar_eV_s"]
        self.H0 = c["H0"] / c["Mpc_km"] * self.hbar  # eV
        self.h = c["H0"] / 100.0
        self.MP = c["MP_reduced_eV"]
        self.kappa = 1.0 / self.MP**2
        self.T0 = c["k_B_eV_K"] * c["T_cmb_K"]  # eV
        self.Tnu0 = (4.0 / 11.0) ** (1.0 / 3.0) * self.T0
        self.rho_c = 3 * self.H0**2 * self.MP**2  # eV^4
        rho_g = math.pi**2 / 15 * self.T0**4
        self.Om_r = rho_g * (1 + 7 / 8 * (4 / 11) ** (4 / 3) * c["N_eff_SM"]) / self.rho_c
        self.Om_m = c["Omega_m"]
        self.Om_L = 1 - self.Om_m - self.Om_r
        self.a_star = 1 / (1 + c["z_star"])
        self.rho_dm0 = c["Omega_c_h2"] / self.h**2 * self.rho_c
        self.rho_b0 = c["Omega_b_h2"] / self.h**2 * self.rho_c

    def H(self, a):
        return self.H0 * math.sqrt(self.Om_r / a**4 + self.Om_m / a**3 + self.Om_L)

    def t_int(self, f, a1=None, a2=1.0):
        """int f(a) dt from a1 (default last scattering) to a2; dt = da/(a H)."""
        a1 = self.a_star if a1 is None else a1
        # integrate in ln a for accuracy
        return quad(lambda la: f(math.exp(la)) / self.H(math.exp(la)), math.log(a1), math.log(a2),
                    limit=400, epsabs=0, epsrel=1e-10)[0]


def beta_rad(c):
    return math.radians(c["beta_deg"])


def report(name, ok):
    print(("PASS" if ok else "FAIL") + ": " + name)
    return bool(ok)
