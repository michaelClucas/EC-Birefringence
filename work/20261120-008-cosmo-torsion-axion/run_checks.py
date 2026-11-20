"""Run all scripts of task 20261120-008 and summarize PASS/FAIL.

Usage (from the repo root, pinned to the program's CPUs):
  .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-008-cosmo-torsion-axion/run_checks.py > work/20261120-008-cosmo-torsion-axion/run.log
Single-threaded; each script runs in its own process (each < 3 min). Order matters:
beta_scan.py writes scan.csv, which variants.py and anisotropy.py read.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["beta_scan.py", "variants.py", "anisotropy.py",
           "check_ode.py", "check_numbers.py", "check_anisotropy.py"]


def main():
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    npass = nfail = 0
    bad = []
    for s in SCRIPTS:
        t0 = time.time()
        r = subprocess.run([sys.executable, s], cwd=HERE, env=env, capture_output=True, text=True)
        out = r.stdout + r.stderr
        p, f = out.count("PASS:"), out.count("FAIL:")
        npass += p
        nfail += f
        print(f"=== {s}: {p} PASS, {f} FAIL, exit {r.returncode}, {time.time() - t0:.1f} s")
        print(out.rstrip(), flush=True)
        if f or r.returncode:
            bad.append(s)
    print(f"\nTOTAL: {npass} PASS, {nfail} FAIL" + (f"; failing: {bad}" if bad else ""))
    print("ALL PASS" if not bad else "SOME FAIL")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
