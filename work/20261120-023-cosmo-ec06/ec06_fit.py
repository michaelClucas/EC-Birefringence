"""Goal 3: fit xi*eps per background shape to the isotropic-beta likelihoods (src.likelihood), compare with a free-beta ALP.

beta(xi eps) = K_shape * xi eps is linear, so the fit is a re-labelling of the free-beta fit: same chi2_min, same Delta chi2.
Parameter counts: k_eff = 1 (xi eps, the only combination beta sees). k_full = 2 for EC-06 with a fixed shape (xi and eps
are separate parameters of the theory; beta fixes only their product; the shape n is fixed by the theory, 022) and 3 if n
is left free. Free-beta ALP: k_full = 2 (task 021). BIC with one compressed beta has ln n = 0: Delta BIC = Delta chi2.
Writes fit_results.csv.
"""
import csv
import math
import os

import numpy as np

from ec06_common import HERE, SHAPES, ZSTAR, DEG, beta_unit
from src.likelihood import IsotropicBeta
from src.compare import summary_1param_vs_null, aic

checks = []


def check(name, ok, info=""):
    checks.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}: {name} -- {info}")


out = []
for ds in ["act_dr6", "act_wmap_planck", "act_cosmoglobe"]:
    lk = IsotropicBeta(ds)
    alp = summary_1param_vs_null(lambda b: lk.chi2(b), 0.1, 1, -1.0, 1.0)
    out.append(dict(dataset=ds, model="free-beta ALP", shape="-", best=alp["best"], lo68=lk.mu - lk.sigma,
                    hi68=lk.mu + lk.sigma, units="deg", chi2_min=alp["chi2_min"], dchi2=alp["delta_chi2"],
                    sigma=alp["sigma"], dAIC_keff1=alp["delta_AIC"], dAIC_kfull=aic(alp["chi2_null"], 0) - aic(alp["chi2_min"], 2),
                    k_full=2, dBIC_n1=alp["delta_chi2"], B10=alp["B10"]))
    for name, n in SHAPES:
        K = beta_unit(ZSTAR, n) / DEG   # deg per unit xi eps
        s = summary_1param_vs_null(lambda x: lk.chi2(K * x), 0.1 / K, 1, -1.0 / K, 1.0 / K)
        lo, hi = (lk.mu - lk.sigma) / K, (lk.mu + lk.sigma) / K
        # 68% interval from Delta chi2 = 1 on a grid (independent of the analytic mu +- sigma)
        xs = np.linspace(s["best"] - 3 * lk.sigma / K, s["best"] + 3 * lk.sigma / K, 60001)
        c = lk.chi2(K * xs)
        inside = xs[c <= s["chi2_min"] + 1.0]
        check(f"F1 [{ds}, {name}] 68% interval from Delta chi2 = 1 grid = (mu +- sigma)/K",
              abs(inside[0] / lo - 1) < 1e-3 and abs(inside[-1] / hi - 1) < 1e-3, f"[{inside[0]:.4e}, {inside[-1]:.4e}]")
        check(f"F2 [{ds}, {name}] Delta chi2 equals free-beta ALP Delta chi2 (rel. 1e-5, Nelder-Mead tolerance)", abs(s["delta_chi2"] / alp["delta_chi2"] - 1) < 1e-5,
              f"{s['delta_chi2']:.4f} vs {alp['delta_chi2']:.4f}")
        check(f"F3 [{ds}, {name}] Savage-Dickey B10 equals ALP's (uniform beta in [-1,1] deg mapped to xi eps)",
              abs(s["B10"] / alp["B10"] - 1) < 1e-4, f"{s['B10']:.3f}")
        kf = 2
        out.append(dict(dataset=ds, model="EC-06 CFJ", shape=name, best=s["best"], lo68=lo, hi68=hi, units="xi*eps",
                        chi2_min=s["chi2_min"], dchi2=s["delta_chi2"], sigma=s["sigma"], dAIC_keff1=s["delta_AIC"],
                        dAIC_kfull=aic(s["chi2_null"], 0) - aic(s["chi2_min"], kf), k_full=kf, dBIC_n1=s["delta_chi2"],
                        B10=s["B10"]))

with open(os.path.join(HERE, "fit_results.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    for r in out:
        w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})

print("\ndataset | model | shape | best [68%] | dchi2 | sigma | dAIC(k=1) | dAIC(k_full=2) | B10")
for r in out:
    print(f"{r['dataset']} | {r['model']} | {r['shape']} | {r['best']:.4g} [{r['lo68']:.4g}, {r['hi68']:.4g}] {r['units']} | "
          f"{r['dchi2']:.3f} | {r['sigma']:.2f} | {r['dAIC_keff1']:.2f} | {r['dAIC_kfull']:.2f} | {r['B10']:.2f}")
r = [x for x in out if x["dataset"] == "act_dr6" and x["shape"] == "H"][0]
check("F4 ACT DR6, shape H: Delta chi2 = (0.215/0.074)^2", abs(r["dchi2"] - (0.215 / 0.074) ** 2) < 1e-4, f"{r['dchi2']:.4f}")
print(f"\n{sum(checks)} PASS, {len(checks) - sum(checks)} FAIL")
