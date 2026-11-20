"""Fit the live axion-class candidate mechanisms to the isotropic beta likelihood and SPT-3G A_CB.

Run from the repo root (1 thread, pinned):
    .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-014-data-reproduce-fit/fit_candidates.py

Models (beta in deg, sign convention of CONVENTIONS.md; every model allows either sign):
  null      beta = 0, A_CB = 0
  ALP       beta free, uniform prior [-1, 1] deg for Savage-Dickey (same prior as task 010)
  EC-03'(c) treated as a generic ALP with g = 8 alpha/(pi f) set by hand (skeptic 011: no torsion-generated g).
            Caps c: 2.7% (Planck 2018 + BOSS) and 6.9% (Planck 2018), 95% at m = 1e-32 eV, from
            data/ultralight_fraction_bounds.csv (2301.08361 Table 2); 5% (the 008 assumption); 12% (needed for
            0.215 deg; above the Planck 2018 limit at 1e-32 eV, labelled EXCL). No per-mass bound exists at 5-8e-33 eV.
            |beta| <= B(c) = B(5%) sqrt(c/5%), c = energy cap rho_phi/rho_DM. B(5%) from work/20261120-008 (via
            src.claims_numbers.ec03). beta is linear in phi_i at fixed m, so a uniform phi_i prior is uniform in beta.
            A_CB = (H_inf g_T F / 4 pi)^2 S_L at fixed beta (from 008 anisotropy.py with phi_i = 2 beta/(g_T F)).
  KK-02(k)  beta = k alpha Delta_theta / 2 pi, |Delta_theta| <= pi, so |beta| <= k alpha/2. Uniform Delta_theta prior.
            Anisotropy envelope A_CB <= (k alpha/2)^2 (theta fluctuating by O(pi), task 009).
Parameter counts: k_eff = parameters the beta data constrain (1 for every non-null model). k_full also counts the
parameters a model needs for beta to be what it is but which the beta data do not constrain:
  ALP: amplitude g phi_i and mass/timing m (2; task 021, after skeptic 015: a physical ALP needs m on the
       plateau and an initial amplitude just like EC-03' and KK-02, so all non-null models count 2)
  EC-03': phi_i and m (m must sit on the plateau 6e-33 .. 3e-31 eV); the cap is fixed per row
  KK-02: Delta_theta and the timing/mass parameter 2 pi L M (k is a discrete model label, not counted)
BIC uses n = 1 for a compressed beta (so ln n = 0 and BIC = chi2: no penalty, stated, not hidden) and n = 6 for
beta + the 5 SPT-3G bandpowers.
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src import compare as C  # noqa: E402
from src import likelihood as L  # noqa: E402
from src import claims_numbers as CN  # noqa: E402

HERE = Path(__file__).resolve().parent
DATASETS = ["act_dr6", "act_wmap_planck"]
ALP_PRIOR = 1.0  # deg


def cosmo_inputs():
    rows = {}
    for name in ("cosmo_torsion_axion_inputs.csv", "bounds_inputs.csv"):
        with open(ROOT / "data" / name, newline="", encoding="utf-8-sig") as f:
            rows.update({r["id"]: float(r["value"]) for r in csv.DictReader(f)})
    return rows


def ultralight_caps(omch2):
    """Omega_a/Omega_DM 95% caps at m = 1e-32 eV (task 013, 2301.08361 Table 2): Omega_a h^2 bound / Omega_c h^2."""
    out = {}
    with open(ROOT / "data" / "ultralight_fraction_bounds.csv", newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["id"] in ("ULF-R23-P18-m1e-32", "ULF-R23-P18BOSS-m1e-32"):
                assert r["units"] == "Omega_a h^2" and r["CL"] == "95"
                out[r["id"]] = float(r["value"].replace("<", "").strip()) / omch2
    return out


def spt_profile():
    """Return (A grid, Delta chi2_SPT(A) = chi2(A) - chi2(0)), chi2 profiled over r in [0,1], A_LT in [0,3]."""
    try:
        spt = L.SPTLite()
    except FileNotFoundError:
        g = L.SPTApprox()
        A = np.linspace(0, 1.5e-4, 301)
        return A, g.chi2(A) - g.chi2(0.0), g.name, g.n_data
    R, AL = np.meshgrid(np.linspace(0, 1, 201), np.linspace(0, 3, 301), indexing="ij")
    A = np.linspace(0, 1.5e-4, 301)
    c2 = np.array([-2 * np.max(spt.loglike(R, AL, a)) for a in A])
    return A, c2 - c2[0], spt.name, spt.n_data


def bounded_fit(lk, B):
    """Best fit of a Gaussian beta likelihood with |beta| <= B: clip the mean."""
    b = float(np.clip(lk.mu, -B, B))
    return b, float(lk.chi2(b))


def main(vals=None):
    """vals: output of src.claims_numbers.regenerate() (reused by src.reproduce); recomputed if None."""
    sys.stdout.reconfigure(encoding="utf-8")
    ci = cosmo_inputs()
    alpha = ci["alpha_em"]
    ec = vals if vals is not None else CN.ec03()
    kkv = vals if vals is not None else CN.kk()
    B5 = ec["EC-03'.beta_max_deg"]
    cap5 = ec["EC-03'.cap_at_best"]
    kk1 = kkv["KK-02.beta_max_k1_deg"]
    assert abs(kk1 / math.degrees(alpha / 2) - 1) < 1e-9, "KK-02 beta_max differs from k alpha/2"

    # EC-03' anisotropy: A_CB = (H_inf g_T F/4 pi)^2 S_L, independent of beta and cap at fixed m; S_L <= 1
    H_max = math.pi * ci["MP_reduced"] * math.sqrt(ci["r_BK18"] * math.exp(ci["logA"]) * 1e-10 / 2)  # eV
    acb_ec03 = (H_max * ec["_EC-03'.g_T_eVinv"] * ec["_EC-03'.F"] / (4 * math.pi)) ** 2

    models = [("ALP (free beta)", ALP_PRIOR, 1, 2, None)]  # k_full = 2 (amplitude, mass), task 021
    ulf = ultralight_caps(ci["omch2"])
    print("EC-03' caps from data/ultralight_fraction_bounds.csv (95%, m = 1e-32 eV; no per-mass bound at 5-8e-33 eV):")
    for k, v in ulf.items():
        print(f"  {k}: Omega_a/Omega_DM < {v:.4f}")
    p18 = ulf["ULF-R23-P18-m1e-32"]
    caps = [(ulf["ULF-R23-P18BOSS-m1e-32"], "P18+BOSS"), (cap5, "008 assumed"), (p18, "P18"), (0.12, "fit 0.215")]
    for cap, src in sorted(caps):
        B = B5 * math.sqrt(cap / cap5)
        phi_ok = ec["_EC-03'.phiE_GeV"] * math.sqrt(cap / cap5) < ec["_EC-03'.phi_excursion_cap_GeV"]
        assert phi_ok, "excursion cap pi f/N_q binds: beta no longer scales as sqrt(cap)"
        tag = " EXCL" if cap > p18 else ""  # above the Planck 2018 95% limit at 1e-32 eV
        models.append((f"EC-03' {cap:.1%} {src}{tag}", B, 1, 2, acb_ec03))
    for k in (1, 2):
        B = k * math.degrees(alpha / 2)
        models.append((f"KK-02 k={k}", B, 1, 2, (k * alpha / 2) ** 2))

    A, dspt, spt_name, spt_n = spt_profile()
    dspt_at = lambda a: float(np.interp(a, A, dspt)) if a <= A[-1] else float("inf")
    mono = bool(np.all(np.diff(dspt) >= -1e-9))
    print(f"SPT likelihood: {spt_name}; Delta chi2_SPT(A_CB) profiled over r, A_LT; monotone in A_CB >= 0: {mono}")
    print(f"  Delta chi2_SPT at A_CB = 1e-5, 5e-5, 1.2e-4: {dspt_at(1e-5):.3f}, {dspt_at(5e-5):.3f}, {dspt_at(1.2e-4):.3f}")
    print(f"  EC-03' A_CB,max = (H_inf,max g_T F/4pi)^2 = {acb_ec03:.2e} (H_inf,max = {H_max / 1e9:.2e} GeV, BK18)")

    out_rows = []
    for ds in DATASETS:
        lk = L.IsotropicBeta(ds)
        c_null = float(lk.chi2(0.0))
        print(f"\n=== {ds}: beta = {lk.mu} +- {lk.sigma} deg; chi2_null = {c_null:.3f}")
        hdr = (f"{'model':28s} {'beta_max':>8s} {'best':>6s} {'chi2':>6s} {'pen':>5s} {'dchi2':>6s} {'sigma':>5s} "
               f"{'dAIC_eff':>8s} {'dAIC_full':>9s} {'dBIC_n1':>7s} {'dBIC_n6':>7s} {'B10':>6s} "
               f"{'SPT dchi2@Amax':>14s}")
        print(hdr)
        c_alp = None
        for name, B, k_eff, k_full, acb_max in models:
            best, c_min = bounded_fit(lk, B)
            if c_alp is None:
                c_alp = c_min
            d = c_null - c_min
            B10, lnB10 = C.savage_dickey(lambda t: -0.5 * lk.chi2(t), -B, B)
            B10a, _ = C.savage_dickey_gaussian(lk.mu, lk.sigma, -B, B)
            assert abs(B10 / B10a - 1) < 1e-4
            # joint with SPT: A_CB in [0, A_max]; SPT prefers A_CB -> 0 (monotone), so the joint best fit has A_CB = 0
            A_best = 0.0
            c_joint = c_min + dspt_at(A_best)
            dj = c_null - c_joint
            spt_worst = dspt_at(acb_max) if acb_max is not None else float("nan")
            r = dict(dataset=ds, model=name, beta_max_deg=B, best_beta_deg=best, chi2_min=c_min, penalty_vs_ALP=c_min - c_alp,
                     dchi2_vs_null=d, sigma=C.significance(d, 1),
                     dAIC_keff=C.aic(c_null, 0) - C.aic(c_min, k_eff), dAIC_kfull=C.aic(c_null, 0) - C.aic(c_min, k_full),
                     dBIC_n1=C.bic(c_null, 0, 1) - C.bic(c_min, k_eff, 1),
                     dBIC_n6_joint=C.bic(c_null, 0, 1 + spt_n) - C.bic(c_joint, k_eff, 1 + spt_n),
                     dBIC_n6_joint_kfull=C.bic(c_null, 0, 1 + spt_n) - C.bic(c_joint, k_full, 1 + spt_n),
                     B10=B10, lnB10=lnB10, k_eff=k_eff, k_full=k_full, acb_max=acb_max if acb_max is not None else float("nan"),
                     spt_dchi2_at_acb_max=spt_worst, joint_dchi2_vs_null=dj)
            out_rows.append(r)
            bm = f"{B:8.3f}" if name != "ALP (free beta)" else "  (free)"
            print(f"{name:28s} {bm} {best:6.3f} {c_min:6.3f} {r['penalty_vs_ALP']:5.2f} {d:6.2f} {r['sigma']:5.2f} "
                  f"{r['dAIC_keff']:8.2f} {r['dAIC_kfull']:9.2f} {r['dBIC_n1']:7.2f} {r['dBIC_n6_joint']:7.2f} {B10:6.2f} "
                  f"{spt_worst:14.3g}")
        need = 100 * cap5 * (lk.mu / B5) ** 2
        print(f"EC-03' cap needed for the best fit {lk.mu} deg: {need:.1f}% rho_DM; KK-02 k=1 needs Delta_theta = "
              f"{lk.mu / kk1:.2f} pi")

    with open(HERE / "fit_results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        for r in out_rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print("\nwrote fit_results.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
