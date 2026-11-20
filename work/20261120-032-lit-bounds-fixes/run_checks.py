"""Run all task-032 scripts (about 1 s each). Redirect to run.log."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
fails = 0
for s in ("verify_inputs.py", "hz_bounds.py", "kk04_shortfall.py", "ec02_yinkomatsu.py"):
    print(f"### {s}", flush=True)
    p = subprocess.run([sys.executable, "-W", "ignore::SyntaxWarning", os.path.join(HERE, s)], cwd=HERE,
                       capture_output=True, text=True)
    print(p.stdout + p.stderr, flush=True)
    fails += int(p.returncode != 0 or "FAIL:" in p.stdout)
print(f"### scripts with failures: {fails}")
sys.exit(1 if fails else 0)
