"""Run all checks of task 20261121-036 and summarize PASS/FAIL.

From the repo root:
  .venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261121-036-theory-pgt0plus-heavy/run_checks.py
Single-threaded; each script in its own process (well under 3 min total). Output goes to run.log.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SCRIPTS = ["check_rc_trace.py", "check_scalar_mode.py", "endpoint_bounds.py"]


def main():
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    npass = nfail = 0
    lines = []
    for s in SCRIPTS:
        t0 = time.time()
        r = subprocess.run([sys.executable, "-W", "ignore::SyntaxWarning", s], cwd=HERE, env=env,
                           capture_output=True, text=True)
        out = (r.stdout + r.stderr).replace(ROOT, "<repo>")
        p, f = out.count("PASS:"), out.count("FAIL:")
        npass += p
        nfail += f + (r.returncode != 0)
        lines.append(f"=== {s}: {p} PASS, {f} FAIL, exit {r.returncode}, {time.time() - t0:.1f} s")
        lines.append(out.rstrip())
    lines.append(f"TOTAL: {npass} PASS, {nfail} FAIL")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(HERE, "run.log"), "w", encoding="utf-8") as fh:
        fh.write(text + "\n")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
