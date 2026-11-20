"""Run all checks of task 20261120-009 and write run.log.

From the repo root:
  .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-009-kk-bounds/run_checks.py
Each script runs in its own process (each under ~3 min on one core).
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["check_parity.py", "check_ec_uholst.py", "check_orbifold.py",
           "check_rotation_hosotani.py", "kk_bounds.py"]


def main():
    lines, npass, nfail = [], 0, 0
    for s in SCRIPTS:
        t = time.time()
        r = subprocess.run([sys.executable, s], cwd=HERE, capture_output=True, text=True)
        out = r.stdout + (r.stderr if r.returncode else "")
        lines.append(f"===== {s} ({time.time() - t:.0f} s, exit {r.returncode})")
        lines.append(out.rstrip())
        npass += out.count("[PASS]")
        nfail += out.count("[FAIL]") + (1 if r.returncode else 0)
    lines.append(f"===== TOTAL: {npass} PASS, {nfail} FAIL")
    text = "\n".join(lines) + "\n"
    print(text)
    with open(os.path.join(HERE, "run.log"), "w", encoding="utf-8") as fh:
        fh.write(text)
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
