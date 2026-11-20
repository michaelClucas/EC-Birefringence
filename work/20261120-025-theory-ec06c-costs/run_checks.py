"""Run all checks of task 20261120-025 and summarize PASS/FAIL.

From the repo root:
  .venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-025-theory-ec06c-costs/run_checks.py
Single-threaded; each script in its own process (about 4 min total). Writes sny_background.csv.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SCRIPTS = ["check_ft_lorentz.py", "check_stueckelberg.py", "ec06c_numbers.py", "sny_background.py", "check_sny_limit.py"]


def main():
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    npass = nfail = 0
    bad = []
    for s in SCRIPTS:
        t0 = time.time()
        r = subprocess.run([sys.executable, "-W", "ignore::SyntaxWarning", s], cwd=HERE, env=env,
                           capture_output=True, text=True)
        out = (r.stdout + r.stderr).replace(ROOT, "<repo>")
        p, f = out.count("PASS:"), out.count("FAIL:")
        npass += p
        nfail += f
        print(f"=== {s}: {p} PASS, {f} FAIL, exit {r.returncode}, {time.time() - t0:.1f} s")
        print(out.rstrip())
        if f or r.returncode:
            bad.append(s)
    print(f"\nTOTAL: {npass} PASS, {nfail} FAIL" + (f"; failing: {bad}" if bad else ""))
    print("ALL PASS" if not bad else "SOME FAIL")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
