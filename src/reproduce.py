"""Reproduce entry point: `python -m src.reproduce` (from the repo root, via tools/pin.py, 1 thread).

Steps:
  1. likelihood validation (task 010): published beta significances and SPT-3G results.
  2. claim numbers (task 014): regenerates every number in the "predicted beta" and "needed parameter vs bound"
     columns of CLAIMS.md by calling the work/ scripts (src/claims_numbers.py) and compares them with the expected
     values in src/claims_expected.csv, whose printed text is also checked to still occur in CLAIMS.md.
  3. candidate fits (task 014): EC-03', KK-02, generic ALP vs null on isotropic beta and SPT A_CB.
Options: --skip-validation (steps 2-3 only), --fast-scan (EC-03' scan on the 17-mass sub-grid instead of all 45).
Exit code is nonzero if any step fails or any claim number does not match.
"""
import importlib.util
import sys
import time

from . import claims_numbers as CN
from .data import ROOT

VALIDATE = ROOT / "work" / "20261120-010-data-likelihood" / "validate.py"
FIT = ROOT / "work" / "20261120-014-data-reproduce-fit" / "fit_candidates.py"


def _load(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def claim_numbers(full_scan):
    vals, times = CN.regenerate(full_scan=full_scan)
    for k, t in times.items():
        print(f"  {k}: {t:.1f} s")
    flags = {k: bool(v) for k, v in vals.items() if k.startswith("_") and str(v) in ("True", "False")}
    for k, v in flags.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k[1:]}")
    npts, diff = vals["_EC-03'.scan_points"], vals["_EC-03'.scan_vs_committed_maxreldiff"]
    print(f"  [INFO] EC-03' scan: {npts} masses, max rel. diff of beta/(g phi_i) vs committed "
          f"work/20261120-008 scan.csv = {diff:.1e}")
    bad = CN.print_table(CN.compare(vals))
    failed_flags = [k for k, v in flags.items() if not v]
    return vals, bool(bad or failed_flags)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    args = sys.argv[1:]
    failed = []
    t0 = time.time()
    if "--skip-validation" not in args:
        print("### likelihood validation (20261120-010)")
        if _load(VALIDATE).main():
            failed.append("likelihood validation")
    print("### claim numbers vs CLAIMS.md (20261120-014)")
    vals, bad = claim_numbers(full_scan="--fast-scan" not in args)
    if bad:
        failed.append("claim numbers")
    print("### candidate fits (20261120-014)")
    if _load(FIT).main(vals):
        failed.append("candidate fits")
    print(f"### reproduce: {'OK' if not failed else f'FAILED: {failed}'} ({time.time() - t0:.0f} s)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
