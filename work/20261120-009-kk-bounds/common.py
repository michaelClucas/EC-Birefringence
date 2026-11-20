"""Shared inputs for task 20261120-009 (KK bounds). Natural units, eV.

All measured numbers are loaded from data/. Nothing measured is hardcoded here.
"""
import csv
import os

import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data")


def _load(name, key="id"):
    with open(os.path.join(DATA, name), newline="", encoding="utf-8-sig") as fh:
        return {row[key]: row for row in csv.DictReader(fh)}


COSMO = _load("cosmo_torsion_axion_inputs.csv")
BINP = _load("bounds_inputs.csv")
BOUNDS = _load("bounds.csv")
MEAS = _load("birefringence_measurements.csv")
KKIN = _load("kk_bounds_inputs.csv")


def cin(k):
    return float(COSMO[k]["value"])


def binp(k):
    return float(BINP[k]["value"])


def kk(k):
    return float(KKIN[k]["value"]), float(KKIN[k]["err"])


def bound_num(id_):
    """Leading numeric part of value_eV (strips <, >, <~ and trailing text)."""
    s = BOUNDS[id_]["value_eV"]
    for t in ("<~", ">~", "<", ">"):
        s = s.replace(t, "")
    return float(s.strip().split()[0])


HBAR = binp("hbar")
HBARC = binp("hbarc")
MPC = binp("Mpc")
MP = binp("MP_reduced")
GEV = binp("GeV")
YR = binp("julian_year")
DEG = np.pi / 180.0

ALPHA = cin("alpha_em")
G4SQ = 4 * np.pi * ALPHA

h = cin("H0") / 100.0
H0 = cin("H0") / (MPC / 1e3) * HBAR
OMM = (cin("ombh2") + cin("omch2")) / h**2
T0 = cin("k_B") * cin("T_CMB")
NEFF_SM = cin("N_eff")
RHOC0 = 3 * H0**2 * MP**2
OMG = (np.pi**2 / 15) * T0**4 / RHOC0
RNU1 = 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)   # one neutrino species / photons (post e+e-)
OMR = OMG * (1 + RNU1 * NEFF_SM)
OML = 1 - OMM - OMR
Z_STAR = cin("z_star")
A_STAR = 1 / (1 + Z_STAR)

BETA_DEG = float(MEAS["DPK25-beta-baseline"]["value"])
BETA_ERR_DEG = float(MEAS["DPK25-beta-baseline"]["err_plus"])
BETA = BETA_DEG * DEG
BETA_ERR = BETA_ERR_DEG * DEG


def E(a):
    return np.sqrt(OMR / a**4 + OMM / a**3 + OML)


def H(a):
    return H0 * E(a)


H_STAR = H(A_STAR)

LOG = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    line = f"[{tag}] {name}" + (f"  ({detail})" if detail else "")
    print(line)
    LOG.append((name, ok))
    return ok


def info(msg):
    print(f"[INFO] {msg}")
