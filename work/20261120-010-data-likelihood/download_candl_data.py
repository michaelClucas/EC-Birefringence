"""Fetch the public candl_data likelihood library (Balkenhol et al.) into artifacts/ (gitignored).

Pinned commit; the SPT-3G D1 BB lite files are about 0.5 MB. Run from the repo root:
    .venv/Scripts/python work/20261120-010-data-likelihood/download_candl_data.py
"""
import subprocess
import sys
from pathlib import Path

URL = "https://github.com/Lbalkenhol/candl_data.git"
COMMIT = "c4eb8fbc5b0eb259ba487494b3cff918eac5c47a"
ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "artifacts" / "candl_data"


def main():
    if not DEST.exists():
        DEST.parent.mkdir(parents=True, exist_ok=True)
        subprocess.check_call(["git", "clone", URL, str(DEST)])
    subprocess.check_call(["git", "-C", str(DEST), "fetch", "--depth", "1", "origin", COMMIT])
    subprocess.check_call(["git", "-C", str(DEST), "checkout", "--detach", COMMIT])
    y = DEST / "candl_data" / "SPT3G_D1_BB_lite_v0" / "SPT3G_D1_BB_lite.yaml"
    print("OK" if y.exists() else "MISSING", y.relative_to(ROOT))
    return 0 if y.exists() else 1


if __name__ == "__main__":
    sys.exit(main())
