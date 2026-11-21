"""Shared helpers for task 20261121-036."""
import csv
import os

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
_count = {"PASS": 0, "FAIL": 0}


def check(name, ok, msg=""):
    tag = "PASS" if ok else "FAIL"
    _count[tag] += 1
    print(f"{tag}: {name} {msg}".rstrip())


def summary():
    print(f"== {_count['PASS']} PASS, {_count['FAIL']} FAIL ==")


def load(name):
    with open(os.path.join(ROOT, "data", name), newline="", encoding="utf-8-sig") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def val(table, key, col="value"):
    return float(table[key][col])
