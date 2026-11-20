"""Validate the likelihood layer (src/) against published numbers in data/. Prints PASS/FAIL/INFO lines.

Run from the repo root:
    .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-010-data-likelihood/validate.py
"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src import compare as C  # noqa: E402
from src import likelihood as L  # noqa: E402
from src.data import measurements  # noqa: E402

RESULTS = []


def report(status, label, msg):
    RESULTS.append((status, label))
    print(f"[{status}] {label}: {msg}")


def check(label, ok, msg):
    report("PASS" if ok else "FAIL", label, msg)


def published_sigma(row):
    m = re.search(r"([\d.]+) sigma", row["notes"])
    return float(m.group(1))


def iso_checks(M):
    for ds in ("act_dr6", "act_wmap_planck", "act_cosmoglobe"):
        lk = L.IsotropicBeta(ds)
        s = C.summary_1param_vs_null(lk.chi2, 0.2, n_data=1, prior_lo=-1.0, prior_hi=1.0)
        pub = published_sigma(M[lk.row_id])
        check(
            f"significance {ds}",
            round(s["sigma"], 1) == pub,
            f"beta = {lk.mu} +- {lk.sigma} deg, Delta chi2 = {s['delta_chi2']:.3f} -> {s['sigma']:.3f} sigma "
            f"(published {pub}); Delta AIC = {s['delta_AIC']:.2f}, Delta BIC(n=1) = {s['delta_BIC']:.2f}, "
            f"B10 (uniform beta in [-1,1] deg) = {s['B10']:.1f} (ln B10 = {s['lnB10']:.2f})",
        )
        B10a, _ = C.savage_dickey_gaussian(lk.mu, lk.sigma, -1.0, 1.0)
        check(f"Savage-Dickey numeric vs analytic {ds}", abs(s["B10"] / B10a - 1) < 1e-4, f"{s['B10']:.4f} vs {B10a:.4f}")

    # recompute the published combinations from components (independent data, as assumed by 2509.13654 Sec. VI)
    for comb, other in (("act_wmap_planck", "wmap_planck_pr4"), ("act_cosmoglobe", "cosmoglobe_wmap_planck")):
        mu, sig = L.inverse_variance_combine([L.IsotropicBeta("act_dr6"), L.IsotropicBeta(other)])
        pub = L.IsotropicBeta(comb)
        check(
            f"combination act_dr6 + {other}",
            abs(mu - pub.mu) < 1e-3 and abs(sig - pub.sigma) < 1e-3,
            f"{mu:.4f} +- {sig:.4f} vs published {pub.mu} +- {pub.sigma}",
        )
    try:
        L.joint("act_dr6", "act_wmap_planck")
        check("joint() refuses overlapping data", False, "no error raised")
    except ValueError as e:
        check("joint() refuses overlapping data", True, str(e))
    j = L.joint("act_dr6", "planck_pr4")
    s = C.summary_1param_vs_null(j.chi2, 0.2, n_data=j.n_data, prior_lo=-1.0, prior_hi=1.0)
    report("INFO", "joint act_dr6 x planck_pr4 (disjoint)", f"best beta = {s['best']:.3f} deg, {s['sigma']:.2f} sigma")
    lk = L.IsotropicBeta("planck_pr4")
    s = C.summary_1param_vs_null(lk.chi2, 0.2, n_data=1, prior_lo=-1.0, prior_hi=1.0)
    report("INFO", "planck_pr4 alone (ESK22-beta)", f"{lk.mu} +- {lk.sigma} deg, {s['sigma']:.2f} sigma, B10 = {s['B10']:.1f}")


def freq_checks(M):
    avail = sorted(k for k, r in M.items() if r["quantity"] in ("psi_alpha_plus_beta_deg", "psi_difference_deg",
                                                                 "beta0_150GHz_deg", "beta_powerlaw_index_n"))
    report("INFO", "per-frequency rows available", ", ".join(avail))

    bands = L.ACTPsiBands(alpha_prior=False)
    mu, sig = bands.achromatic_fit()
    pub = M["DPK25-psi-mean"]
    check("psi bands, stat only, weighted mean", abs(mu - pub["value"]) < 1.5e-3 and abs(sig - pub["err_plus"]) < 1.5e-3,
          f"{mu:.4f} +- {sig:.4f} vs published {pub['value']} +- {pub['err_plus']}")

    bands = L.ACTPsiBands()
    mu, sig = bands.achromatic_fit()
    ref = M["DPK25-beta-PA"]
    report("INFO", "psi bands + Pi^PA alpha prior, achromatic beta",
           f"{mu:.3f} +- {sig:.3f} deg (approximate reconstruction, EB only, no PA4 f220); "
           f"full analysis Table 2 Pi^PA (EB+TB): {ref['value']} +- {ref['err_plus']}")
    c2, best = C.fit_min(lambda b0, n: bands.chi2(lambda nu: L.beta_powerlaw(nu, b0, n)), [0.2, 0.0])
    c0, b0 = C.fit_min(lambda b0: bands.chi2(lambda nu: L.beta_powerlaw(nu, b0, 0.0)), [0.2])
    report("INFO", "psi bands power law beta0 (nu/150)^n",
           f"beta0 = {best[0]:.3f} deg, n = {best[1]:.2f}; chi2 {c2:.2f} vs achromatic {c0:.2f} "
           f"(Delta chi2 = {c0 - c2:.2f}, 1 extra param, {C.significance(c0 - c2, 1):.2f} sigma); "
           f"Delta AIC = {C.aic(c0, 1) - C.aic(c2, 2):.2f}, Delta BIC(n=4) = {C.bic(c0, 1, 4) - C.bic(c2, 2, 4):.2f}")
    d = L.PA5Difference(include_alpha_prior=False)
    dp = L.PA5Difference(include_alpha_prior=True)
    report("INFO", "PA5 f150 - f090 difference",
           f"{d.mu} +- {d.sigma:.3f} (stat), +- {dp.sigma:.3f} with Pi^PA alpha-difference prior; achromatic chi2 = "
           f"{d.chi2(lambda nu: 0.2 + 0 * nu):.2f} / {dp.chi2(lambda nu: 0.2 + 0 * nu):.2f}")
    pl = L.PlanckPowerLaw("all")
    report("INFO", "Planck PR4 power law (ESK22, all bands)", f"beta0 = {pl.b}, n = {pl.n}; chi2(n=0, beta0=0.29) = {pl.chi2(0.29, 0.0):.2f}")


def spt_checks(M):
    try:
        spt = L.SPTLite()
    except FileNotFoundError as e:
        report("FAIL", "SPT lite likelihood", str(e))
        return
    import candl
    import yaml

    like = candl.Like(L.SPT3G_YAML.as_posix(), feedback=False)
    tests = L.CANDL_DIR / "tests"
    cfg = yaml.safe_load(open(tests / "SPT3G_D1_BB_lite_test.yaml"))
    spec = np.loadtxt(tests / cfg["test_spectrum"])
    ll = float(like.log_like({"Dl": {"BB": spec[: like.ell_max - 1, 4]}}))
    check("candl SPT-3G D1 BB lite own test", abs(ll - cfg["test_logl"]) < 1e-6 * abs(cfg["test_logl"]),
          f"{ll:.6f} vs {cfg['test_logl']:.6f}")

    pts = [(0.0, 1.0, 0.0), (0.05, 0.7, 1e-4), (0.1, 1.2, -1e-4), (0.2, 0.5, 3e-4), (0.0, 0.9, 5e-5)]
    mx, off = spt.check_against_candl(pts)
    check("numpy offset-lognormal vs candl", mx < 1e-8, f"max |diff| = {mx:.2e} after constant {off:.4f}")

    s_cb = spt.cb_scale
    report("INFO", "CB template calibration",
           f"thin-LSS template x {s_cb:.3f} reproduces the published A_CB>=0, no-A_LT-prior limit by construction "
           f"(raw thin-LSS limit = {s_cb * M['BAL25-ACB-SPT3G']['value']:.3e}); Namikawa 2404.13771 abstract: exact spectrum "
           "about 2x below thin-LSS at l >~ 100, so a scale near 0.5 is expected")
    check("CB template scale consistent with thin-LSS overestimate (0.35-0.65)", 0.35 < s_cb < 0.65, f"{s_cb:.3f}")

    grid = np.linspace(0, 6e-4, 241)
    post = spt.acb_posterior(grid)
    ul = L.grid_upper_limit(post["acb"], post["post_acb"])
    pub = M["BAL25-ACB-SPT3G"]["value"]
    report("INFO", "SPT-3G A_CB 95% UL (no A_LT prior; calibration target)", f"{ul:.3e} vs published < {pub:.2e}")
    rul = L.grid_upper_limit(post["r"], post["post_r"])
    pub = M["BAL25-r-SPT3G"]["value"]
    check("SPT-3G r 95% UL (A_CB >= 0, no A_LT prior)", abs(rul / pub - 1) < 0.1, f"{rul:.3f} vs published < {pub}")

    m = re.search(r"N\(([\d.]+), ([\d.]+)\^2\)", M["BAL25-ACB-SPT3G-APS"]["dataset"])
    prior = (float(m.group(1)), float(m.group(2)))
    post_p = spt.acb_posterior(grid, alt_prior=prior)
    ulp = L.grid_upper_limit(post_p["acb"], post_p["post_acb"])
    pub = M["BAL25-ACB-SPT3G-APS"]["value"]
    check("SPT-3G A_CB 95% UL (A_LT APS prior)", abs(ulp / pub - 1) < 0.2,
          f"{ulp:.3e} vs published < {pub:.2e} (ratio {ulp / pub:.3f}; tol 20%; raw thin-LSS {ulp * s_cb:.3e})")

    post_n = spt.acb_posterior(np.linspace(-6e-4, 6e-4, 481), acb_nonneg=False)
    mu, sd = L.grid_mean_sd(post_n["acb"], post_n["post_acb"])
    pub = M["BAL25-ACB-SPT3G-neg"]
    check("SPT-3G A_CB negative allowed", abs(mu - pub["value"]) < 0.25 * pub["err_plus"] and abs(sd / pub["err_plus"] - 1) < 0.15,
          f"{mu:.3e} +- {sd:.3e} vs published {pub['value']:.2e} +- {pub['err_plus']:.2e} (tol 0.25 sigma, 15%)")
    ma, sa = L.grid_mean_sd(post_n["alt"], post_n["post_alt"])
    neg = post_n["acb"] < 0
    frac_neg = np.trapezoid(post_n["post_acb"][neg], post_n["acb"][neg]) / np.trapezoid(post_n["post_acb"], post_n["acb"])
    pub_frac = float(re.search(r"(\d+)% mass below zero", pub["notes"]).group(1)) / 100
    report("INFO", "negative-allowed extras", f"A_LT = {ma:.2f} +- {sa:.2f}; P(A_CB<0) = {frac_neg:.2f} (published {pub_frac})")
    ma, sa = L.grid_mean_sd(post["alt"], post["post_alt"])
    pub = M["BAL25-ALT-SPT3G"]
    check("SPT-3G A_LT (A_CB >= 0)", abs(ma - pub["value"]) < 0.25 * pub["err_plus"],
          f"{ma:.3f} +- {sa:.3f} vs published {pub['value']} +- {pub['err_plus']}")

    ap = L.SPTApprox()
    ula = ap.upper_limit()
    pub = M["BAL25-ACB-SPT3G"]["value"]
    check("SPT Gaussian-approx A_CB 95% UL (truncated at 0)", abs(ula / pub - 1) < 0.15, f"{ula:.3e} vs published < {pub:.2e}")


def main():
    M = measurements()
    print("== (a) isotropic beta ==")
    iso_checks(M)
    print("== (b) frequency dependence ==")
    freq_checks(M)
    print("== (c) SPT anisotropy ==")
    spt_checks(M)
    n_fail = sum(s == "FAIL" for s, _ in RESULTS)
    n_pass = sum(s == "PASS" for s, _ in RESULTS)
    print(f"== {n_pass} PASS, {n_fail} FAIL ==")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
