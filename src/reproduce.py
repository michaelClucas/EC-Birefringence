"""Reproduce entry point: `python -m src.reproduce` (from the repo root, via tools/pin.py).

Currently runs the likelihood validation (task 20261120-010). Each claim added to CLAIMS.md gets a step here.
Exit code is nonzero if any step fails.
"""
import importlib.util
import sys

from .data import ROOT

STEPS = [
    ("likelihood validation (20261120-010)", ROOT / "work" / "20261120-010-data-likelihood" / "validate.py"),
]


def _run(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.main()


def main():
    failed = []
    for label, path in STEPS:
        print(f"### {label}")
        if _run(path):
            failed.append(label)
    print("### reproduce:", "OK" if not failed else f"FAILED: {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
