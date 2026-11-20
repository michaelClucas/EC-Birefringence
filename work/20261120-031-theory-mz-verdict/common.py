"""Shared helpers for task 20261120-031 (MZ / EC-06a verdict). Conventions: CONVENTIONS.md.

Frame indices A, B = 0..3, eta = diag(-1, 1, 1, 1), eps_{0123} = +1 so eps^{0123} = -1.
p-forms are stored as fully antisymmetric component arrays w[mu][nu]... with w = (1/p!) w_{mu nu ...} dx^mu ^ dx^nu ...
"""
import csv
import itertools
import os

import sympy as sp

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
R4 = range(4)
ETA = sp.diag(-1, 1, 1, 1)
_count = {"PASS": 0, "FAIL": 0}


def check(name, ok, msg=""):
    tag = "PASS" if ok else "FAIL"
    _count[tag] += 1
    print(f"{tag}: {name} {msg}".rstrip())


def summary():
    print(f"== {_count['PASS']} PASS, {_count['FAIL']} FAIL ==")


def load(name, key="id"):
    with open(os.path.join(ROOT, "data", name), newline="", encoding="utf-8") as f:
        return {r[key]: r for r in csv.DictReader(f)}


def eps_lo(*idx):
    return sp.LeviCivita(*idx)


def eps_up(*idx):
    return -sp.LeviCivita(*idx)


# ---------- forms in coordinate components ----------
def zero(p):
    if p == 1:
        return [sp.Integer(0)] * 4
    if p == 2:
        return [[sp.Integer(0)] * 4 for _ in R4]
    return [[[sp.Integer(0)] * 4 for _ in R4] for _ in R4]


def wedge11(a, b):
    return [[a[m] * b[n] - a[n] * b[m] for n in R4] for m in R4]


def wedge12(a, B):
    """(a ^ B)_{mnr} = a_m B_nr + a_n B_rm + a_r B_mn."""
    return [[[a[m] * B[n][r] + a[n] * B[r][m] + a[r] * B[m][n] for r in R4] for n in R4] for m in R4]


def d1(a, X):
    return [[sp.diff(a[n], X[m]) - sp.diff(a[m], X[n]) for n in R4] for m in R4]


def add(*forms, coeffs=None):
    coeffs = coeffs or [1] * len(forms)
    def rec(objs):
        if isinstance(objs[0], list):
            return [rec([o[i] for o in objs]) for i in range(len(objs[0]))]
        return sum(c * o for c, o in zip(coeffs, objs))
    return rec(list(forms))


def simp(form):
    if isinstance(form, list):
        return [simp(f) for f in form]
    return sp.simplify(form)


def is_zero(form):
    if isinstance(form, list):
        return all(is_zero(f) for f in form)
    return sp.simplify(form) == 0


def perms3():
    return list(itertools.permutations(range(3)))
