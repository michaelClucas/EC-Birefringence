"""Model comparison helpers. chi2 = -2 ln L throughout (additive constants dropped consistently)."""
import math

import numpy as np
from scipy import optimize, stats


def delta_chi2(chi2_model, chi2_null):
    """chi2_null - chi2_model (positive = model preferred)."""
    return chi2_null - chi2_model


def pvalue_from_delta_chi2(dchi2, dof):
    """Wilks: Delta chi2 ~ chi2(dof) under the null."""
    return float(stats.chi2.sf(max(dchi2, 0.0), dof))


def sigma_from_pvalue(p):
    """Two-sided Gaussian-equivalent significance."""
    return float(stats.norm.isf(p / 2))


def significance(dchi2, dof=1):
    return sigma_from_pvalue(pvalue_from_delta_chi2(dchi2, dof))


def aic(chi2_min, k):
    """AIC = chi2_min + 2k, k = number of free parameters fitted (state it explicitly at the call site)."""
    return chi2_min + 2 * k


def bic(chi2_min, k, n):
    """BIC = chi2_min + k ln n, n = number of data points in the likelihood (e.g. 1 for a compressed beta)."""
    return chi2_min + k * math.log(n)


def fit_min(chi2_fn, x0, bounds=None):
    """Minimise chi2 over parameters; returns (chi2_min, best_params)."""
    x0 = np.atleast_1d(np.asarray(x0, float))
    res = optimize.minimize(lambda x: float(chi2_fn(*x)), x0, bounds=bounds, method="L-BFGS-B" if bounds else "Nelder-Mead")
    return float(res.fun), res.x


def savage_dickey(loglike_fn, prior_lo, prior_hi, theta0=0.0, n=20001):
    """Bayes factor B_10 (1-parameter model with uniform prior on [prior_lo, prior_hi] vs nested null theta = theta0).
    Savage-Dickey: B_01 = p(theta0 | d) / pi(theta0). Returns B_10 = 1/B_01 and ln B_10."""
    if not prior_lo <= theta0 <= prior_hi:
        raise ValueError("theta0 must lie inside the prior")
    th = np.linspace(prior_lo, prior_hi, n)
    ll = np.array([float(loglike_fn(t)) for t in th]) if not _vectorises(loglike_fn) else np.asarray(loglike_fn(th), float)
    ll0 = float(loglike_fn(theta0))
    m = ll.max()
    Z = np.trapezoid(np.exp(ll - m), th)  # integral of likelihood over theta (unnormalised posterior)
    post0 = math.exp(ll0 - m) / Z
    prior0 = 1.0 / (prior_hi - prior_lo)
    B01 = post0 / prior0
    return 1.0 / B01, -math.log(B01)


def savage_dickey_gaussian(mu, sigma, prior_lo, prior_hi, theta0=0.0):
    """Analytic check of savage_dickey for a Gaussian likelihood N(mu, sigma)."""
    Z = sigma * math.sqrt(2 * math.pi) * (stats.norm.cdf((prior_hi - mu) / sigma) - stats.norm.cdf((prior_lo - mu) / sigma))
    post0 = math.exp(-0.5 * ((theta0 - mu) / sigma) ** 2) / Z
    B01 = post0 * (prior_hi - prior_lo)
    return 1.0 / B01, -math.log(B01)


def _vectorises(fn):
    try:
        out = np.asarray(fn(np.array([0.0, 1.0])))
        return out.shape == (2,)
    except Exception:
        return False


def summary_1param_vs_null(chi2_fn, x0, n_data, prior_lo, prior_hi, theta0=0.0):
    """Delta chi2, significance, Delta AIC, Delta BIC (null minus model; positive favours model), Savage-Dickey B_10.
    Model: 1 free parameter. Null: 0 free parameters (theta fixed at theta0)."""
    c_null = float(chi2_fn(theta0))
    c_min, best = fit_min(chi2_fn, [x0])
    d = delta_chi2(c_min, c_null)
    B10, lnB10 = savage_dickey(lambda t: -0.5 * chi2_fn(t), prior_lo, prior_hi, theta0)
    return {
        "best": float(best[0]),
        "chi2_min": c_min,
        "chi2_null": c_null,
        "delta_chi2": d,
        "sigma": significance(d, 1),
        "delta_AIC": aic(c_null, 0) - aic(c_min, 1),
        "delta_BIC": bic(c_null, 0, n_data) - bic(c_min, 1, n_data),
        "B10": B10,
        "lnB10": lnB10,
        "prior": (prior_lo, prior_hi),
    }
