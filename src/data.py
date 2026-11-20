"""Load the cited numbers in data/. No measured number is hardcoded anywhere in src/."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
ARTIFACTS = ROOT / "artifacts"


def _read(name):
    with open(DATA / name, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def _num(s):
    return float(s) if s not in (None, "") else None


def measurements():
    """dict id -> row of data/birefringence_measurements.csv, numeric fields parsed."""
    out = {}
    for r in _read("birefringence_measurements.csv"):
        for k in ("value", "err_minus", "err_plus", "CL"):
            r[k] = _num(r[k])
        r["upper_limit_flag"] = int(r["upper_limit_flag"] or 0)
        out[r["id"]] = r
    return out


def planck2018_bestfit():
    return {r["id"]: float(r["value"]) for r in _read("planck2018_bestfit.csv")}


def act_dr6_bands():
    rows = _read("act_dr6_bands.csv")
    for r in rows:
        for k in ("nu_nominal_GHz", "nu_low_GHz", "nu_high_GHz"):
            r[k] = float(r[k])
    return {r["array_band"]: r for r in rows}
