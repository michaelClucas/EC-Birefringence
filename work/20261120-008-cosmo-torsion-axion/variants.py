"""Items 1, 3, 4: couplings, QCD masses and beta for the EC-03 variants; required g(m)
for a generic ALP and for a torsion-fixed coupling without QCD mass.

Needs scan.csv from beta_scan.py. Writes g_required.png and variants_summary.csv.
All angles in radians internally.
"""
import csv
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import (HERE, inp, BOUNDS, bound_upper, GEV, CM3_TO_EV3, RHO_DM0, ALPHA, BETA_DEG,
                    BETA_ERR_DEG, DEG, MP, HBAR)

BETA = BETA_DEG * DEG
BETA_LO = (BETA_DEG - BETA_ERR_DEG) * DEG
RHO_LOC = inp("rho_loc_DM") * GEV * CM3_TO_EV3          # GeV/cm^3 -> eV^4
MAFA = inp("maxfa_QCD")
CQCD = inp("cgamma_QCD")
EN = inp("EoverN_universal")
F_CF = inp("f_CF")
NQ = inp("Nq_CF")
F_LM = inp("f_LM")
G_CHANDRA = bound_upper("A-Chandra")                      # eV^-1
Z_STAR = inp("z_star")
RAY = inp("rayleigh95_SOZ")
WASH = float(re.findall(r"[\d.eE+-]+", BOUNDS["U-washout-J0"]["value_eV"])[0])


def load_scan():
    with open(os.path.join(HERE, "scan.csv")) as fh:
        rows = list(csv.DictReader(fh))
    out = {k: np.array([float(r[k]) if r[k] not in ("True", "False") else r[k] == "True"
                        for r in rows]) for k in rows[0]}
    return out


def j0_inverse(y):
    """x with 1 - J0(x) = y (small y), by bisection."""
    from scipy.special import j0
    lo, hi = 0.0, 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if 1 - j0(mid) < y else (lo, mid)
    return 0.5 * (lo + hi)


def frac_cap(m):
    """Max rho_phi/rho_DM today from data/bounds.csv (U-frac-*), assumptions in README."""
    hl = bound_upper("U-frac-Hlozek")            # 0.05 for 1e-32 < m < 10^-25.5
    rog = BOUNDS["U-frac-Rogers23"]["value_eV"]  # "< 0.1 ; < 0.01"
    r_hi, r_lo = [float(x) for x in re.findall(r"<\s*([\d.]+)", rog)]
    cap = np.full_like(m, np.nan)
    for i, mm in enumerate(m):
        c = []
        if mm <= 1e-32:
            c.append(hl)                         # assumption: DE-like regime, same 5% cap
        if 1e-32 < mm < 10**-25.5:
            c.append(hl)
        if 1e-30 <= mm <= 1e-28:
            c.append(r_lo)
        if mm <= 1e-26:
            c.append(r_hi)
        if mm > 1e-26:
            c.append(r_hi)                       # assumption: 10% beyond Rogers' range
        cap[i] = min(c)
    return cap


def qcd_variants():
    print("== Item 1/3: torsion axion with QCD anomaly mass ==")
    cases = [("KR/BI, 1502.03694, f_a = f/N_q", F_CF / NQ),
             ("BI Nieh-Yan, 0911.2698, f_a = f", F_LM),
             ("KR/BI, 1502.03694, f_a = f (as in their Eq. 27)", F_CF)]
    rows = []
    for name, fa in cases:
        m = MAFA / fa
        g = ALPHA / (2 * np.pi * fa) * abs(EN + CQCD)
        g_noqcd = ALPHA / (2 * np.pi * fa) * EN
        phi_loc = np.sqrt(2 * RHO_LOC) / m                  # all local DM in this field
        b_loc = 0.5 * g * phi_loc                           # instantaneous amplitude
        period_s = 2 * np.pi * HBAR / m
        b_frozen = 0.5 * g * np.pi * fa                     # if it did not oscillate (no mass)
        b_frozen_noqcd = 0.5 * g_noqcd * np.pi * fa
        short = np.log10(BETA / b_loc)
        rows.append(dict(variant=name, f_a_GeV=fa / GEV, m_eV=m, g_GeVinv=g * GEV,
                         phi_loc_GeV=phi_loc / GEV, beta_loc_amp_deg=b_loc / DEG,
                         period_s=period_s, orders_short=short))
        print(f"{name}: f_a = {fa/GEV:.3e} GeV, m = {m:.3e} eV, g = {g*GEV:.3e} GeV^-1 "
              f"(E/N - 1.92 = {EN+CQCD:.3f}); local phi amplitude (all DM) = {phi_loc/GEV:.3e} GeV; "
              f"beta amplitude = {b_loc:.2e} rad = {b_loc/DEG:.2e} deg, period {period_s:.2e} s; "
              f"shortfall 10^{short:.1f}")
        print(f"    comparison: field frozen at pi f_a would give beta = {b_frozen/DEG:.3f} deg "
              f"(with QCD mixing), {b_frozen_noqcd/DEG:.3f} deg (E/N only)")
    # f_a independence of the oscillating amplitude: beta = alpha c sqrt(2 rho)/(4 pi m f_a)
    b_univ = ALPHA * abs(EN + CQCD) * np.sqrt(2 * RHO_LOC) / (4 * np.pi * MAFA)
    print(f"f_a-independent amplitude alpha|E/N-1.92| sqrt(2 rho_loc)/(4 pi m_a f_a) = {b_univ:.3e} rad")
    return rows, b_univ


def required_g(sc):
    m = sc["m"]
    F = 2 * sc["beta_maxabs"]                 # |phi0 - phibar_LSS| / phi_i (best phase)
    cap = frac_cap(m)
    phiE = np.sqrt(cap * RHO_DM0 / sc["rho0"])   # eV
    greq = 2 * BETA / (F * phiE)
    greq_lo = 2 * BETA_LO / (F * phiE)
    return m, F, cap, phiE, greq, greq_lo


def soz_window():
    """Reproduce 2512.21888 Eqs. (4), (6), (12) at m = 1e-24 eV, f_phi = 1 (local-field regime)."""
    m = 1e-24
    phi0 = np.sqrt(2 * RHO_LOC) / m
    g_req = 2 * BETA / (RAY * phi0)
    rho_star = RHO_DM0 * (1 + Z_STAR) ** 3
    x = j0_inverse(WASH)
    g_wash = x / (np.sqrt(2 * rho_star) / m)
    return phi0, g_req, g_wash, x


def main():
    qrows, b_univ = qcd_variants()
    sc = load_scan()
    m, F, cap, phiE, greq, greq_lo = required_g(sc)
    ok = m <= 1.5e-27                         # beyond this the LSS term is washed out
    print("\n== Item 4: generic ALP, homogeneous field, LSS + observer term ==")
    print(f"{'m [eV]':>10} {'F':>8} {'cap':>6} {'phi_i,max [GeV]':>15} {'g_min [GeV^-1]':>15} {'(g phi_i)_needed':>16}")
    for i in np.where(ok)[0][::3]:
        print(f"{m[i]:10.3e} {F[i]:8.3e} {cap[i]:6.3f} {phiE[i]/GEV:15.3e} {greq[i]*GEV:15.3e} {2*BETA/F[i]:16.3e}")
    i_best = np.argmin(np.where(ok, greq, np.inf))
    print(f"lowest required g: {greq[i_best]*GEV:.2e} GeV^-1 at m = {m[i_best]:.2e} eV "
          f"(phi_i = {phiE[i_best]/GEV:.2e} GeV = {phiE[i_best]/MP:.2f} M_P); Chandra bound "
          f"{G_CHANDRA*GEV:.1e} GeV^-1")
    allowed = ok & (greq < G_CHANDRA)
    print(f"m range with g_min < Chandra: {m[allowed].min():.2e} .. {m[allowed].max():.2e} eV")
    print(f"(g f)_min if phi_i <= pi f: 2 beta/(pi F) = {2*BETA/(np.pi*F[i_best]):.2e} (dimensionless)")

    phi0, g_req, g_wash, x = soz_window()
    print(f"\n2512.21888 regime (m=1e-24 eV, all DM): phi_0 = {phi0/GEV:.3e} GeV, "
          f"g_req = {g_req*GEV:.3e} GeV^-1, washout g_max = {g_wash*GEV:.3e} GeV^-1 "
          f"(x_max = {x:.4f}); ratio {g_req/g_wash:.2f}")

    # torsion-fixed coupling with photon anomaly but no QCD mass (assumption)
    print("\n== Item 4: torsion-fixed coupling g = alpha Qbar^2/(pi f) without a QCD mass ==")
    g_T = ALPHA / (2 * np.pi * (F_CF / NQ)) * EN          # = alpha * 8 / (pi f)
    tor = {}
    for lab, phicap in (("phi_i <= pi f/N_q", np.pi * F_CF / NQ), ("phi_i <= pi f", np.pi * F_CF)):
        bT = 0.5 * g_T * F * np.minimum(phicap, phiE)
        sel = ok & (bT >= BETA_LO)
        j = np.argmax(np.where(ok, bT, -1))
        rng = f"{m[sel].min():.2e} .. {m[sel].max():.2e} eV" if sel.any() else "none"
        print(f"g_T = {g_T*GEV:.3e} GeV^-1, {lab}: max beta = {bT[j]/DEG:.3f} deg at m = {m[j]:.2e} eV;"
              f" beta >= {BETA_LO/DEG:.3f} deg for m in {rng}")
        tor[lab] = bT

    # plot
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=110)
    mm = m[ok]
    ax.fill_between(mm, greq[ok] * GEV, G_CHANDRA * GEV, where=greq[ok] < G_CHANDRA,
                    color="C2", alpha=0.25, label=r"generic ALP: $\beta=0.215^\circ$ reachable")
    ax.loglog(mm, greq[ok] * GEV, "C2-", lw=2, label=r"$g_{\min}(m)$ (energy cap, best phase)")
    ax.loglog(mm, greq_lo[ok] * GEV, "C2:", lw=1, label=r"$g_{\min}$ for $\beta-1\sigma$")
    ax.axhline(G_CHANDRA * GEV, color="C3", ls="--", label="Chandra NGC 1275 (1907.05475)")
    ax.axhline(g_T * GEV, color="C0", ls="-.", label=r"torsion anomaly $g=8\alpha/(\pi f)$, no QCD mass")
    for r in qrows:
        ax.plot([r["m_eV"]], [r["g_GeVinv"]], "kx")
    ax.annotate("torsion axion with QCD mass\n" r"($\beta\approx 0$, washed out)",
                xy=(qrows[1]["m_eV"], qrows[1]["g_GeVinv"]), xytext=(1e-24, 1e-19), fontsize=7,
                arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.axvspan(1.5e-27, 1e-11, color="0.85", alpha=0.6, lw=0)
    ax.text(1e-21, 3e-11, "LSS term washed out;\nsingle-field local-DM\nregime: 2512.21888", fontsize=7)
    ax.set_xlim(1e-35, 1e-11)
    ax.set_ylim(1e-23, 1e-9)
    ax.set_xlabel("m [eV]")
    ax.set_ylabel(r"$g_{\phi\gamma}$ [GeV$^{-1}$]")
    ax.legend(fontsize=6.5, loc="upper left")
    ax.grid(alpha=0.3, which="both", lw=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "g_required.png"))

    with open(os.path.join(HERE, "variants_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(qrows[0].keys()))
        w.writeheader()
        for r in qrows:
            w.writerow({k: (f"{v:.4e}" if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote g_required.png, variants_summary.csv")


if __name__ == "__main__":
    main()
