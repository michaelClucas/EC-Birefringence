"""Run every check script of task 043 and write run.log (no absolute paths). Run from the repo root:
.venv/Scripts/python -W ignore::SyntaxWarning tools/pin.py .venv/Scripts/python -W ignore::SyntaxWarning work/20261121-043-theory-alp-costs/run_checks.py
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SCRIPTS = ["ec04_tuning.py", "kk02_hilltop.py", "kk03_kick.py", "ec06a_remnant.py"]


def main():
    lines, fails = [], 0
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for s in SCRIPTS:
        p = subprocess.run([sys.executable, "-W", "ignore::SyntaxWarning", os.path.join(HERE, s)], cwd=HERE,
                           capture_output=True, text=True, env=env)
        out = (p.stdout + p.stderr).replace(ROOT, "<repo>").replace(ROOT.replace("\\", "/"), "<repo>")
        lines.append(f"===== {s} (exit {p.returncode}) =====\n{out}")
        fails += p.returncode != 0 or "[FAIL]" in out
    npass = sum(l.count("[PASS]") for l in lines)
    nfail = sum(l.count("[FAIL]") for l in lines)
    lines.append(f"===== TOTAL: {npass} PASS, {nfail} FAIL; scripts failing: {fails} =====\n")
    with open(os.path.join(HERE, "run.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(lines[-1])
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
