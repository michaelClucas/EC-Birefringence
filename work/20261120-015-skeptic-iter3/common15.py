"""Loaders for task 015 (skeptic). Measured numbers come from data/; nothing measured is hardcoded."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
NPASS = [0, 0]


def table(name, key="id"):
    with open(DATA / name, newline="", encoding="utf-8-sig") as f:
        return {r[key]: r for r in csv.DictReader(f)}


def val(name, k):
    return float(table(name)[k]["value"])


def pair(id_):
    s = table("bounds.csv")[id_]["value_eV"].replace("+-", " ").split()
    return float(s[0]), float(s[1])


def check(label, ok, extra=""):
    NPASS[0 if ok else 1] += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({extra})" if extra else ""))


def info(s):
    print(f"[INFO] {s}")
