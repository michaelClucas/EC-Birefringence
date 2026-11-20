"""Run all checks for task 20261120-004-kk-reduction. Prints PASS/FAIL per check.

Run from the repo root, pinned:
  .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-004-kk-reduction/run_checks.py
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import check_kk_ricci, check_torsion, check_cs, check_forms  # noqa: E402

n_fail = 0
for mod in (check_kk_ricci, check_torsion, check_cs, check_forms):
    t0 = time.time()
    for name, ok, det in mod.run():
        n_fail += not ok
        print(("PASS" if ok else "FAIL"), "|", name, "|", det)
    print(f"  [{mod.__name__}: {time.time() - t0:.1f} s]")
print("ALL PASS" if n_fail == 0 else f"{n_fail} FAIL")
sys.exit(1 if n_fail else 0)
