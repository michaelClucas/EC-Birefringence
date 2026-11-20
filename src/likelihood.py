"""Likelihood layer for fitting candidate mechanisms to birefringence data.

Sign convention (CONVENTIONS.md): beta > 0 is a clockwise rotation on the sky (Minami-Komatsu / IAU), the
convention of arXiv:2509.13654. Every beta row in data/ is already in this convention. A model computing
beta = +1/2 int p_tau dtau (CFJ) or beta = -(g/2)(phi_0 - phi_LSS) feeds that number in directly.
All inputs and outputs here are in degrees for beta and rad^2 for A_CB.

(a) IsotropicBeta: Gaussian (or split-normal) likelihood on a published isotropic beta.
    Combining data sets: only published combined values are used. `joint()` refuses to multiply
    likelihoods whose underlying data overlap.
(b) Frequency dependence beta(nu) = beta0 (nu/nu0)^n:
    - PA5Difference: published PA5 f150 - f090 psi difference (2509.13654 Sec. IV).
    - ACTPsiBands: approximate reconstruction from the four published EB-only per-array-band
      psi_i = alpha_i + beta(nu_i) (stat only) with the optics-model alpha prior (Pi^PA) marginalised
      analytically. Approximate: no cross-band stat correlation (not published), no TB, no PA4 f220.
    - PlanckPowerLaw: Eskilt 2022 (2201.13347) beta0, n posteriors as independent split normals
      (their correlation is not published).
(c) SPT A_CB:
    - SPTLite: SPT-3G D1 BB lite (candl_data, offset log-normal) with r, A_LT, A_CB; templates from
      src/templates.py. Needs artifacts/ download (work/20261120-010-data-likelihood/download_candl_data.py).
    - SPTApprox: Gaussian in A_CB from the published A_CB = (-1.13 +- 0.94)e-4 (negative allowed).
"""
import math

import numpy as np

from .data import ARTIFACTS, act_dr6_bands, measurements

# --------------------------------------------------------------------------------------------- (a)

# name -> (row id in data/birefringence_measurements.csv, set of underlying data)
ISO_DATASETS = {
    "act_dr6": ("DPK25-beta-baseline", {"ACT_DR6"}),
    "act_wmap_planck": ("DPK25-comb-EK22", {"ACT_DR6", "WMAP", "PLANCK"}),
    "act_cosmoglobe": ("DPK25-comb-CG", {"ACT_DR6", "WMAP", "PLANCK"}),
    "planck_pr4": ("ESK22-beta", {"PLANCK"}),
    "wmap_planck_pr4": ("EK22-beta", {"WMAP", "PLANCK"}),
    "cosmoglobe_wmap_planck": ("CG23-beta-cited", {"WMAP", "PLANCK"}),
}


class IsotropicBeta:
    """-2 ln L = ((beta - mu)/sigma)^2, sigma = err_plus above the mean and err_minus below if split=True,
    else the mean of the two (default; matches how 2509.13654 combined EK22)."""

    n_data = 1

    def __init__(self, dataset="act_dr6", split=False):
        rid, self.components = ISO_DATASETS[dataset]
        row = measurements()[rid]
        assert row["quantity"] == "beta_iso_deg" and row["CL"] == 68, rid
        self.name, self.row_id = dataset, rid
        self.mu = row["value"]
        self.sig_lo, self.sig_hi = row["err_minus"], row["err_plus"]
        if not split:
            self.sig_lo = self.sig_hi = 0.5 * (self.sig_lo + self.sig_hi)

    def chi2(self, beta_deg):
        b = np.asarray(beta_deg, float)
        sig = np.where(b >= self.mu, self.sig_hi, self.sig_lo)
        return ((b - self.mu) / sig) ** 2

    def loglike(self, beta_deg):
        return -0.5 * self.chi2(beta_deg)

    @property
    def sigma(self):
        return 0.5 * (self.sig_lo + self.sig_hi)


class Joint:
    """Product of likelihoods with disjoint underlying data. Each member needs .components and .chi2(*args)."""

    def __init__(self, *likes):
        seen = set()
        for lk in likes:
            if seen & lk.components:
                raise ValueError(
                    f"{lk.name} shares data {sorted(seen & lk.components)} with another member; "
                    "use a published combined value instead of multiplying"
                )
            seen |= lk.components
        self.likes, self.components = likes, seen
        self.name = "+".join(lk.name for lk in likes)
        self.n_data = sum(lk.n_data for lk in likes)

    def chi2(self, *args):
        return sum(lk.chi2(*args) for lk in self.likes)

    def loglike(self, *args):
        return -0.5 * self.chi2(*args)


def joint(*names):
    return Joint(*[IsotropicBeta(n) for n in names])


def inverse_variance_combine(likes):
    """Published-style inverse-variance combination of independent Gaussian beta's (for validation only)."""
    w = np.array([1 / lk.sigma ** 2 for lk in likes])
    mu = np.array([lk.mu for lk in likes])
    return float(w @ mu / w.sum()), float(w.sum() ** -0.5)


# --------------------------------------------------------------------------------------------- (b)

def beta_powerlaw(nu, beta0, n, nu0=150.0):
    return beta0 * (np.asarray(nu, float) / nu0) ** n


class PA5Difference:
    """psi(PA5 f150) - psi(PA5 f090) = 0.21 +- 0.09 deg (stat only, EB only). Model: beta(nu_150) - beta(nu_090).
    include_alpha_prior adds the variance of alpha_150 - alpha_090 under the Pi^PA prior (2 s^2 (1 - rho))."""

    n_data = 1
    components = {"ACT_DR6"}
    name = "act_pa5_diff"

    def __init__(self, include_alpha_prior=True, freq="nominal"):
        m = measurements()
        row = m["DPK25-psi-PA5diff"]
        bands = act_dr6_bands()
        self.mu = row["value"]
        var = row["err_plus"] ** 2
        if include_alpha_prior:
            s = m["DPK25-alphaprior-PA5"]["value"]
            rho = m["DPK25-alphaprior-corr"]["value"]
            var += 2 * s * s * (1 - rho)
        self.sigma = math.sqrt(var)
        key = "nu_nominal_GHz" if freq == "nominal" else freq
        self.nu_hi = bands["PA5f150"][key]
        self.nu_lo = bands["PA5f090"][key]

    def chi2(self, beta_of_nu):
        d = beta_of_nu(self.nu_hi) - beta_of_nu(self.nu_lo)
        return ((d - self.mu) / self.sigma) ** 2

    def loglike(self, beta_of_nu):
        return -0.5 * self.chi2(beta_of_nu)


class ACTPsiBands:
    """psi_i = alpha_i + beta(nu_i) for PA5/PA6 f090/f150 (EB only, stat errors), alpha_i ~ N(0, C_alpha) with the
    Pi^PA prior (sigma per PA, correlation rho within a PA, independent across PAs). Marginalising alpha:
    psi ~ N(beta(nu), C_stat + C_alpha). Approximate (see module docstring)."""

    BANDS = ["PA5f090", "PA5f150", "PA6f090", "PA6f150"]
    components = {"ACT_DR6"}
    name = "act_psi_bands"
    n_data = 4

    def __init__(self, freq="nominal", alpha_corr=None, alpha_prior=True):
        m = measurements()
        bands = act_dr6_bands()
        key = "nu_nominal_GHz" if freq == "nominal" else freq
        self.psi = np.array([m[f"DPK25-psi-{b}"]["value"] for b in self.BANDS])
        stat = np.array([m[f"DPK25-psi-{b}"]["err_plus"] for b in self.BANDS])
        self.nu = np.array([bands[b][key] for b in self.BANDS])
        pa = [bands[b]["PA"] for b in self.BANDS]
        sa = np.array([m[f"DPK25-alphaprior-{p}"]["value"] for p in pa])
        rho = m["DPK25-alphaprior-corr"]["value"] if alpha_corr is None else alpha_corr
        corr = np.array([[1.0 if i == j else (rho if pa[i] == pa[j] else 0.0) for j in range(4)] for i in range(4)])
        self.cov = np.diag(stat ** 2) + (np.outer(sa, sa) * corr if alpha_prior else 0.0)
        self.icov = np.linalg.inv(self.cov)

    def chi2(self, beta_of_nu):
        d = self.psi - beta_of_nu(self.nu)
        return float(d @ self.icov @ d)

    def loglike(self, beta_of_nu):
        return -0.5 * self.chi2(beta_of_nu)

    def achromatic_fit(self):
        one = np.ones(4)
        F = one @ self.icov @ one
        return float(one @ self.icov @ self.psi / F), float(F ** -0.5)


class PlanckPowerLaw:
    """Eskilt 2022 beta_nu = beta0 (nu/150 GHz)^n, split-normal in beta0 and n separately (no published correlation)."""

    components = {"PLANCK"}
    n_data = 2
    VARIANTS = {"all": "all", "dust": "dust", "noFG": "noFG"}

    def __init__(self, variant="all"):
        m = measurements()
        b, n = m[f"ESK22-beta0-{variant}"], m[f"ESK22-n-{variant}"]
        self.name = f"planck_pr4_powerlaw_{variant}"
        self.b = (b["value"], b["err_minus"], b["err_plus"])
        self.n = (n["value"], n["err_minus"], n["err_plus"])
        self.nu0 = 150.0

    @staticmethod
    def _sn(x, p):
        mu, lo, hi = p
        return ((x - mu) / (hi if x >= mu else lo)) ** 2

    def chi2(self, beta0, n):
        return self._sn(beta0, self.b) + self._sn(n, self.n)

    def loglike(self, beta0, n):
        return -0.5 * self.chi2(beta0, n)


# --------------------------------------------------------------------------------------------- (c)

CANDL_DIR = ARTIFACTS / "candl_data" / "candl_data"
SPT3G_YAML = CANDL_DIR / "SPT3G_D1_BB_lite_v0" / "SPT3G_D1_BB_lite.yaml"


class SPTLite:
    """SPT-3G D1 BB lite, offset log-normal (Prince, Calabrese & Dunkley 2403.00085 Eq. 7, as in candl):
        -ln L = 1/2 (y - mu)^T C^-1 (y - mu) + sum y,   y = ln(W^T D - offset),
        D_l = r D^r + A_LT D^lens + A_CB D^CB.
    Vectorised re-implementation of candl's offset_lognorm_logl; `check_against_candl` compares to candl.

    D^CB: our template is the thin-LSS result (src/templates.py). Namikawa 2024 (2404.13771, abstract) finds the
    exact spectrum used by 2510.07928 is about 2x lower at l >~ 100. cb_scale multiplies our template:
      cb_scale="calibrated" (default): one number s fixed so that the A_CB >= 0, no-A_LT-prior 95% limit equals the
        published one (row BAL25-ACB-SPT3G); the A_LT-prior limit and the negative-allowed posterior are then
        independent checks (validate.py). s is cached in artifacts/templates/cb_scale.json.
      cb_scale=1.0: raw thin-LSS template."""

    components = {"SPT3G_D1"}
    name = "spt3g_d1_bb_lite"
    n_data = 5

    def __init__(self, templates=None, cb_scale="calibrated"):
        if not SPT3G_YAML.exists():
            raise FileNotFoundError(
                "candl_data not found; run work/20261120-010-data-likelihood/download_candl_data.py"
            )
        d = SPT3G_YAML.parent
        self.mean = np.loadtxt(d / "SPT3G_D1_BB_lite_offset_lognorm_mean.txt")
        self.offset = np.loadtxt(d / "SPT3G_D1_BB_lite_offset_lognorm_offsets.txt")
        self.cov = np.loadtxt(d / "SPT3G_D1_BB_lite_offset_lognorm_cov.txt")
        self.icov = np.linalg.inv(self.cov)
        win = np.loadtxt(d / "windows" / "BB_lxl_window_functions.txt")
        self.win_ell = win[:, 0].astype(int)  # 2..3001
        self.win = win[:, 1:]
        if templates is None:
            from .templates import spt_bb_templates

            templates = spt_bb_templates()
        self.t = templates
        sl = self.win_ell
        self.b_r = templates["dl_r"][sl] @ self.win
        self.b_lens = templates["dl_lens"][sl] @ self.win
        self.b_cb_raw = templates["dl_cb"][sl] @ self.win
        self.cb_scale = self._calibrated_scale() if cb_scale == "calibrated" else float(cb_scale)
        self.b_cb = self.cb_scale * self.b_cb_raw
        self.dl_cb = self.cb_scale * templates["dl_cb"]

    def _calibrated_scale(self):
        import json

        from .templates import CACHE

        f = CACHE.parent / "cb_scale.json"
        pub = measurements()["BAL25-ACB-SPT3G"]
        stamp = {"template": CACHE.name, "b_cb_raw": [float(v) for v in self.b_cb_raw], "ul_pub": pub["value"]}
        if f.exists():
            c = json.loads(f.read_text())
            if all(c.get(k) == v for k, v in stamp.items()):
                return c["scale"]
        self.b_cb = self.b_cb_raw
        post = self.acb_posterior(np.linspace(0, 6e-4, 241))
        ul_raw = grid_upper_limit(post["acb"], post["post_acb"], pub["CL"] / 100)
        stamp.update(scale=ul_raw / pub["value"], ul_raw=ul_raw)
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps(stamp, indent=1))
        return stamp["scale"]

    def binned(self, r, A_LT, A_CB):
        r, A_LT, A_CB = (np.asarray(v, float)[..., None] for v in (r, A_LT, A_CB))
        return r * self.b_r + A_LT * self.b_lens + A_CB * self.b_cb

    def loglike(self, r, A_LT, A_CB):
        b = self.binned(r, A_LT, A_CB) - self.offset
        ok = np.all(b > 0, axis=-1)
        y = np.log(np.where(b > 0, b, 1.0))
        dy = y - self.mean
        ll = -0.5 * np.einsum("...i,ij,...j->...", dy, self.icov, dy) - y.sum(-1)
        return np.where(ok, ll, -np.inf)

    def check_against_candl(self, points):
        """Return max |ln L_ours - ln L_candl| over points [(r, A_LT, A_CB), ...] (constant offset removed)."""
        import candl

        like = candl.Like(SPT3G_YAML.as_posix(), feedback=False)
        ours, theirs = [], []
        for r, a, c in points:
            dl = r * self.t["dl_r"] + a * self.t["dl_lens"] + c * self.dl_cb
            theirs.append(float(like.log_like({"Dl": {"BB": dl[2 : like.ell_max + 1]}})))
            ours.append(float(self.loglike(r, a, c)))
        diff = np.array(ours) - np.array(theirs)
        return float(np.max(np.abs(diff - diff.mean()))), float(diff.mean())

    def acb_posterior(self, acb_grid, r_max=1.0, alt_max=3.0, n_r=121, n_alt=241, alt_prior=None, acb_nonneg=True):
        """Marginal posterior of A_CB on acb_grid with uniform priors r in [0, r_max], A_LT in [0, alt_max]
        (as in 2510.07928 Sec. 3.1), optional Gaussian A_LT prior (mean, sigma), optional A_CB >= 0."""
        rg = np.linspace(0, r_max, n_r)
        ag = np.linspace(0, alt_max, n_alt)
        acb = np.asarray(acb_grid, float)
        if acb_nonneg:
            acb = acb[acb >= 0]
        R, A = np.meshgrid(rg, ag, indexing="ij")
        ll = np.stack([self.loglike(R, A, c) for c in acb], axis=-1)  # (n_r, n_alt, n_acb), loop keeps memory small
        if alt_prior is not None:
            ll = ll - 0.5 * ((A[..., None] - alt_prior[0]) / alt_prior[1]) ** 2
        ll -= np.max(ll)
        p = np.exp(ll)
        post_acb = np.trapezoid(np.trapezoid(p, rg, axis=0), ag, axis=0)
        post_alt = np.trapezoid(np.trapezoid(p, rg, axis=0), acb, axis=1)
        post_r = np.trapezoid(np.trapezoid(p, ag, axis=1), acb, axis=1)
        return {"acb": acb, "post_acb": post_acb, "alt": ag, "post_alt": post_alt, "r": rg, "post_r": post_r}


def grid_upper_limit(x, p, cl=0.95):
    c = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(x))])
    c /= c[-1]
    return float(np.interp(cl, c, x))


def grid_mean_sd(x, p):
    n = np.trapezoid(p, x)
    m = np.trapezoid(x * p, x) / n
    return float(m), float(np.sqrt(np.trapezoid((x - m) ** 2 * p, x) / n))


class SPTApprox:
    """Fallback: Gaussian in A_CB from the published negative-allowed posterior A_CB = (-1.13 +- 0.94)e-4
    (2510.07928 Sec. 3.2, no A_LT prior). Use with a prior A_CB >= 0 for the physical case."""

    components = {"SPT3G_D1"}
    name = "spt3g_acb_gauss_approx"
    n_data = 1

    def __init__(self):
        row = measurements()["BAL25-ACB-SPT3G-neg"]
        self.mu, self.sigma = row["value"], 0.5 * (row["err_minus"] + row["err_plus"])

    def chi2(self, A_CB):
        return ((np.asarray(A_CB, float) - self.mu) / self.sigma) ** 2

    def loglike(self, A_CB):
        return -0.5 * self.chi2(A_CB)

    def upper_limit(self, cl=0.95):
        """Upper limit with A_CB >= 0 prior (truncated Gaussian), analytic."""
        from scipy.stats import norm

        z0 = (0 - self.mu) / self.sigma
        p0 = norm.cdf(z0)
        return float(self.mu + self.sigma * norm.ppf(p0 + cl * (1 - p0)))
