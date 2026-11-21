"""Task 044: run skeptic 039's three scripts and check that 039's local inputs agree with data/skeptic039_inputs.csv.

Usage (from the repo root, pinned, 1 thread):
  .venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261121-044-reproduce-039/run_checks.py
Each 039 script runs in its own process with cwd = the 039 folder (they locate data/ and inputs_039.csv from their
own __file__). A script fails on a nonzero exit, any FAIL: line, or a missing "== n PASS, 0 FAIL ==" summary.
Exit code 1 if anything fails. Writes no files.
"""
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
T039 = os.path.join(ROOT, "work", "20261121-039-skeptic-036-037")
SCRIPTS = ["rc_trace_indep.py", "endpoint_grain.py", "stueckelberg_cutoff.py"]
SUMMARY = re.compile(r"== (\d+) PASS, 0 FAIL ==")


def rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def main():
    bad = []
    npass = nfail = 0
    # C0: 039's local inputs equal the data/ copy (039's scripts cannot be edited to read data/)
    loc, dat = rows(os.path.join(T039, "inputs_039.csv")), rows(os.path.join(ROOT, "data", "skeptic039_inputs.csv"))
    same = set(loc) == set(dat) and all(float(loc[k]["value"]) == float(dat[k]["value"]) for k in loc)
    print(("PASS" if same else "FAIL") + f": C0 inputs_039.csv ({len(loc)} rows) equals data/skeptic039_inputs.csv "
          f"({len(dat)} rows) in ids and values")
    npass += same
    nfail += not same
    if not same:
        bad.append("inputs")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", PYTHONIOENCODING="utf-8")
    for s in SCRIPTS:
        t0 = time.time()
        r = subprocess.run([sys.executable, "-W", "ignore::SyntaxWarning", s], cwd=T039, env=env, capture_output=True,
                           text=True, encoding="utf-8", errors="replace")
        out = r.stdout + r.stderr
        p, f = len(re.findall(r"^PASS:", out, re.M)), len(re.findall(r"^FAIL:", out, re.M))
        m = SUMMARY.search(out)
        ok = r.returncode == 0 and f == 0 and m is not None and int(m.group(1)) == p
        npass += p
        nfail += f
        print(f"=== 039 {s}: {p} PASS, {f} FAIL, exit {r.returncode}, {time.time() - t0:.1f} s")
        print(out.rstrip())
        if not ok:
            bad.append(s)
    print(f"\nTOTAL: {npass} PASS, {nfail} FAIL" + (f"; failing: {bad}" if bad else ""))
    print("ALL PASS" if not bad else "SOME FAIL")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
