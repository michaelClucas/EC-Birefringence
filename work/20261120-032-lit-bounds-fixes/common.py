"""Shared helpers for task 032: CSV loading by id and PASS/FAIL checks."""
import csv
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
_res = {"PASS": 0, "FAIL": 0}


def load(name, folder="data"):
    with open(os.path.join(ROOT, folder, name), newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def val(d, k, col="value"):
    return float(d[k][col])


def check(name, ok, msg=""):
    tag = "PASS" if ok else "FAIL"
    _res[tag] += 1
    print(f"{tag}: {name} {msg}")


def summary():
    print(f"== {_res['PASS']} PASS, {_res['FAIL']} FAIL")
    return _res["FAIL"]
