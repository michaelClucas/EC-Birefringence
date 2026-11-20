"""Spot-check src/claims_expected.csv rows that come from 014 fit_candidates by running the source script.
Output csv is redirected to this folder (no write into work/20261120-014). Also recomputes the 014 Savage-Dickey
B10 with a common prior width for every model, to separate the prior (Occam) factor from the fit.
"""
import csv
import importlib.util
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("fit14", ROOT / "work" / "20261120-014-data-reproduce-fit" / "fit_candidates.py")
fit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fit)
fit.HERE = HERE / "artifacts_fit"
fit.HERE.mkdir(exist_ok=True)
fit.main()

rows = list(csv.DictReader(open(fit.HERE / "fit_results.csv", encoding="utf-8")))
get = lambda ds, m: next(r for r in rows if r["dataset"] == ds and r["model"].startswith(m))
exp = {r["key"]: r for r in csv.DictReader(open(ROOT / "src" / "claims_expected.csv", encoding="utf-8-sig"))}
for key, ds, m in (("penalty_cap27_AWP", "act_wmap_planck", "EC-03' 2.7%"), ("penalty_k1_AWP", "act_wmap_planck", "KK-02 k=1")):
    v = float(get(ds, m)["penalty_vs_ALP"])
    e = float(exp[key]["expected"])
    print(f"[{'PASS' if round(v, 1) == e else 'FAIL'}] claims_expected {key}: regenerated {v:.3f}, expected {e}")
b = float(get("act_wmap_planck", "EC-03' 2.7%")["beta_max_deg"])
print(f"[{'PASS' if round(b, 3) == float(exp['beta_max_cap27_deg']['expected']) else 'FAIL'}] beta_max_cap27_deg {b:.4f}")

# Savage-Dickey with the bounded models' physical cut but a common outer width: B10 = Z/(L(0) * width)
from src import likelihood as L  # noqa: E402
from scipy.stats import norm  # noqa: E402
for ds in ("act_dr6", "act_wmap_planck"):
    lk = L.IsotropicBeta(ds)
    mu, s = lk.mu, lk.sigma
    L0 = norm.pdf(0, mu, s)
    for B in (1.0, 0.418, 0.2091, 0.0998):
        frac = norm.cdf(B, mu, s) - norm.cdf(-B, mu, s)
        print(f"{ds}: uniform prior +-{B:.4f} deg: B10 = {frac / (2 * B * L0):.3g}  (likelihood mass kept {frac:.2f})")
