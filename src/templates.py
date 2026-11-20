"""BB templates for the SPT-3G D1 BB lite likelihood: D_l^r (r=1), D_l^lens, D_l^CB (A_CB = 1 rad^2).

All D_l = l(l+1)C_l/2pi in muK^2, l = 0..LMAX_OUT.

D^r and D^lens: CAMB at the Planck 2018 best fit (data/planck2018_bestfit.csv), tensors with r=1, n_t=0,
pivot 0.05/Mpc.

D^CB: full-sky thin-last-scattering-surface result for a Gaussian rotation field alpha(n) with
L(L+1)C_L^aa/2pi = A_CB (L >= 1) acting on the unlensed E modes:
    C_l^BB = (1/pi) sum_{L,l'} (2L+1)(2l'+1) C_L^aa C_l'^EE (l L l'; 2 0 -2)^2 [1 + (-1)^(l+L+l')]/2
evaluated exactly with Wigner-d correlation functions on a Gauss-Legendre grid:
    C_l^BB = 1/(4 pi) int dx xi_a(x) [d^l_22(x) xi_+(x) + d^l_2-2(x) xi_-(x)],
    xi_a = sum (2L+1) C_L^aa P_L,  xi_+- = sum (2l'+1) C_l'^EE d^l'_2(+-2).
Balkenhol et al. (2510.07928) used the exact (non-thin-LSS) template of Lonappan et al. / Namikawa 2024,
which needs a modified CLASS build; the thin-LSS approximation is the documented difference.
"""
import numpy as np

from .data import ARTIFACTS, planck2018_bestfit

LMAX_OUT = 3001
LMAX_IN = 4000
CACHE = ARTIFACTS / "templates" / "spt_bb_templates_v1.npz"


def wigner_d_iter(m, mp, x, lmax):
    """Yield (l, d^l_{m mp}(x)) for l = max(|m|,|mp|)..lmax by the standard three-term recursion.
    Implemented for (m, mp) in {(0,0), (2,2), (2,-2)}."""
    if (m, mp) == (0, 0):
        prev, cur, l0 = np.zeros_like(x), np.ones_like(x), 0
    elif (m, mp) == (2, 2):
        prev, cur, l0 = np.zeros_like(x), ((1 + x) / 2) ** 2, 2
    elif (m, mp) == (2, -2):
        prev, cur, l0 = np.zeros_like(x), ((1 - x) / 2) ** 2, 2
    else:
        raise ValueError((m, mp))
    yield l0, cur
    for l in range(l0, lmax):
        a = (2 * l + 1) * (l * (l + 1) * x - m * mp)
        b = (l + 1) * np.sqrt(max((l * l - m * m) * (l * l - mp * mp), 0))
        c = l * np.sqrt(((l + 1) ** 2 - m * m) * ((l + 1) ** 2 - mp * mp)) if l > 0 else 1.0
        if l == 0:  # Legendre start: P_1 = x
            nxt = x * cur
        else:
            nxt = (a * cur - b * prev) / c
        prev, cur = cur, nxt
        yield l + 1, cur


def cb_bb_from_ee(cl_ee, cl_aa, lmax_out, npts=None):
    """C_l^BB (l=0..lmax_out) from rotation field spectrum cl_aa and E spectrum cl_ee (both C_l arrays from l=0)."""
    lmax_e = len(cl_ee) - 1
    lmax_a = len(cl_aa) - 1
    if npts is None:
        npts = (lmax_out + lmax_e + lmax_a) // 2 + 10
    x, w = np.polynomial.legendre.leggauss(npts)
    xi_a = np.zeros_like(x)
    for l, d in wigner_d_iter(0, 0, x, lmax_a):
        if cl_aa[l]:
            xi_a += (2 * l + 1) * cl_aa[l] * d
    xi_p = np.zeros_like(x)
    xi_m = np.zeros_like(x)
    for (l, dp), (_, dm) in zip(wigner_d_iter(2, 2, x, lmax_e), wigner_d_iter(2, -2, x, lmax_e)):
        xi_p += (2 * l + 1) * cl_ee[l] * dp
        xi_m += (2 * l + 1) * cl_ee[l] * dm
    out = np.zeros(lmax_out + 1)
    gp = w * xi_a * xi_p
    gm = w * xi_a * xi_m
    for (l, dp), (_, dm) in zip(wigner_d_iter(2, 2, x, lmax_out), wigner_d_iter(2, -2, x, lmax_out)):
        out[l] = (dp @ gp + dm @ gm) / (4 * np.pi)
    return out


def scale_invariant_cl_aa(lmax, A=1.0):
    L = np.arange(lmax + 1, dtype=float)
    cl = np.zeros(lmax + 1)
    cl[1:] = 2 * np.pi * A / (L[1:] * (L[1:] + 1))
    return cl


def dl_to_cl(dl):
    l = np.arange(len(dl), dtype=float)
    cl = np.zeros_like(dl)
    cl[2:] = 2 * np.pi * dl[2:] / (l[2:] * (l[2:] + 1))
    return cl


def cl_to_dl(cl):
    l = np.arange(len(cl), dtype=float)
    return l * (l + 1) * cl / (2 * np.pi)


def camb_spectra():
    import camb

    p = planck2018_bestfit()
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"], mnu=p["mnu"], tau=p["tau"])
    pars.InitPower.set_params(As=np.exp(p["logA"]) * 1e-10, ns=p["ns"], r=1.0, nt=0.0, pivot_tensor=0.05)
    pars.WantTensors = True
    pars.set_for_lmax(LMAX_IN, lens_potential_accuracy=1)
    res = camb.get_results(pars)
    ps = res.get_cmb_power_spectra(pars, CMB_unit="muK", lmax=LMAX_IN, raw_cl=False)
    return {
        "dl_ee_unlensed": ps["unlensed_scalar"][: LMAX_IN + 1, 1],
        "dl_bb_lens": ps["lensed_scalar"][: LMAX_IN + 1, 2],
        "dl_bb_r1": ps["tensor"][: LMAX_IN + 1, 2],
    }


def spt_bb_templates(force=False):
    """Return dict with ell, dl_r, dl_lens, dl_cb (per unit r, per A_LT=1, per A_CB=1 rad^2). Cached in artifacts/."""
    if CACHE.exists() and not force:
        z = np.load(CACHE)
        return {k: z[k] for k in z.files}
    s = camb_spectra()
    cl_ee = dl_to_cl(s["dl_ee_unlensed"])
    cl_bb_cb = cb_bb_from_ee(cl_ee, scale_invariant_cl_aa(LMAX_IN), LMAX_OUT)
    out = {
        "ell": np.arange(LMAX_OUT + 1),
        "dl_r": s["dl_bb_r1"][: LMAX_OUT + 1],
        "dl_lens": s["dl_bb_lens"][: LMAX_OUT + 1],
        "dl_cb": cl_to_dl(cl_bb_cb),
        "dl_ee_unlensed": s["dl_ee_unlensed"][: LMAX_OUT + 1],
    }
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, **out)
    return out
