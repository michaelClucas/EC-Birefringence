"""Run every script of task 20261120-023 and summarize PASS/FAIL.

Usage (from the repo root, pinned, 1 thread; total about 1 min):
  .venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261120-023-cosmo-ec06/run_checks.py > work/20261120-023-cosmo-ec06/run.log
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SCRIPTS = ["check_symbolic.py", "ec06_beta.py", "ec06_fit.py", "ec06_aniso.py", "ec06_kaf.py", "ec06_published_maps.py"]


def main():
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", ""))
    npass = nfail = 0
    bad = []
    for s in SCRIPTS:
        t0 = time.time()
        r = subprocess.run([sys.executable, "-W", "ignore::SyntaxWarning", s], cwd=HERE, env=env,
                           capture_output=True, text=True)
        out = (r.stdout + r.stderr).replace(ROOT, "<repo>")
        p = sum(1 for ln in out.splitlines() if ln.startswith("PASS:"))
        f = sum(1 for ln in out.splitlines() if ln.startswith("FAIL:"))
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
