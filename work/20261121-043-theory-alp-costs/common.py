"""Shared loader for task 20261121-043 (Result C costs). Natural units, eV. All measured inputs from data/."""
import csv
import math
import os

import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data")


def load(name):
    with open(os.path.join(DATA, name), newline="", encoding="utf-8-sig") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


COSMO = load("cosmo_torsion_axion_inputs.csv")
BINP = load("bounds_inputs.csv")
BOUNDS = load("bounds.csv")
MEAS = load("birefringence_measurements.csv")
KKIN = load("kk_bounds_inputs.csv")
ALPC = load("alp_costs_inputs.csv")
EC06T = load("ec06_theory_inputs.csv")
MZIN = load("mz_inputs.csv")


def cin(k):
    return float(COSMO[k]["value"])


def binp(k):
    return float(BINP[k]["value"])


def bound_num(id_):
    s = BOUNDS[id_]["value_eV"]
    for t in ("<~", ">~", "<", ">"):
        s = s.replace(t, "")
    return float(s.strip().split()[0])


DEG = math.pi / 180
HBAR = binp("hbar")
HBARC = binp("hbarc")
MPC = binp("Mpc")
MP = binp("MP_reduced")
GEV = binp("GeV")
ALPHA = cin("alpha_em")

h = cin("H0") / 100
H0 = cin("H0") / (MPC / 1e3) * HBAR
OMM = (cin("ombh2") + cin("omch2")) / h**2
T0 = cin("k_B") * cin("T_CMB")
NEFF_SM = cin("N_eff")
RHOC0 = 3 * H0**2 * MP**2
OMG = (math.pi**2 / 15) * T0**4 / RHOC0
RNU1 = 7 / 8 * (4 / 11) ** (4 / 3)
OMR = OMG * (1 + RNU1 * NEFF_SM)
OML = 1 - OMM - OMR
A_STAR = 1 / (1 + cin("z_star"))


def E(a):
    return math.sqrt(OMR / a**4 + OMM / a**3 + OML)


H_STAR = H0 * E(A_STAR)

# measurements: (central, 1 sigma) in degrees
ACT = (float(MEAS["DPK25-beta-baseline"]["value"]), float(MEAS["DPK25-beta-baseline"]["err_plus"]))
ESK = (float(ALPC["ESK26-beta-joint"]["value"]), float(ALPC["ESK26-beta-joint"]["err"]))

NP = {"PASS": 0, "FAIL": 0}


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    NP[tag] += 1
    print(f"[{tag}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def info(msg):
    print(f"[INFO] {msg}")


def done():
    print(f"== {NP['PASS']} PASS, {NP['FAIL']} FAIL ==")
    return 0 if NP["FAIL"] == 0 else 1
