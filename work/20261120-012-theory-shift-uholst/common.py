"""Shared helpers for task 20261120-012 (shift symmetry, u-Holst with a hidden scalar).

Measured numbers come from data/. Nothing measured is hardcoded here.
"""
import csv
import os

import numpy as np

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
NEW = _load("shift_uholst_inputs.csv")


def cin(k):
    return float(COSMO[k]["value"])


def binp(k):
    return float(BINP[k]["value"])


def new(k):
    return float(NEW[k]["value"])


def bound_pair(id_):
    """(central, sigma) from a 'a +- b' value_eV entry of bounds.csv."""
    s = BOUNDS[id_]["value_eV"].replace("+-", " ").split()
    return float(s[0]), float(s[1])


def bound_num(id_):
    s = BOUNDS[id_]["value_eV"]
    for t in ("<~", ">~", "<", ">"):
        s = s.replace(t, "")
    return float(s.strip().split()[0])


DEG = np.pi / 180.0
MP = binp("MP_reduced")          # eV
HBAR = binp("hbar")              # eV s
HBARC = binp("hbarc")            # eV m
YR = binp("julian_year")         # s
MPC = binp("Mpc")
ALPHA = cin("alpha_em")
M_E = cin("m_e")
BETA_ACT = float(MEAS["DPK25-beta-baseline"]["value"]) * DEG   # rad
H0 = cin("H0") / (MPC / 1e3) * HBAR                             # eV
OM = (cin("ombh2") + cin("omch2")) / (cin("H0") / 100) ** 2

LOG = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  ({detail})" if detail else ""))
    LOG.append((name, ok))
    return ok


def info(msg):
    print(f"[INFO] {msg}")
