"""Independent checks of the anisotropic-birefringence BB template (src/templates.py).

1. Wigner-d recursion: orthonormality int d^l_22 d^l'_22 dx = 2/(2l+1) delta_ll'.
2. Isotropic limit: only C_0^aa nonzero -> C_l^BB = C_0^aa C_l^EE / pi = 4 <alpha^2> C_l^EE.
3. Flat-sky convolution (independent formula, polar integration):
   C^BB(l) = 4 int d^2L/(2pi)^2 C^aa(L) C^EE(|l-L|) cos^2(2(phi_{l-L} - phi_l)), L >= 1,
   compared with the full-sky Wigner-d result at l = 200, 500, 1000, 2000 using the cached CAMB E spectrum.
   Both use C_L^aa only for L >= L0 = 30: for L ~ 1 the flat-sky integral does not approximate the discrete sum
   (sum_L (2L+1)/(2L(L+1)) vs int dL/(L+1) differ by gamma - 1/2 + ln 2 = 0.77, about 10% of the total), so the
   comparison tests the high-L kernel; the low-L part is covered by the exact isotropic-limit check.
Run: .venv/Scripts/python tools/pin.py .venv/Scripts/python work/20261120-010-data-likelihood/check_templates.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.templates import (cb_bb_from_ee, dl_to_cl, scale_invariant_cl_aa, spt_bb_templates,  # noqa: E402
                           wigner_d_iter)

fails = 0


def check(label, ok, msg):
    global fails
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {label}: {msg}")


x, w = np.polynomial.legendre.leggauss(300)
d = {l: v.copy() for l, v in wigner_d_iter(2, 2, x, 60)}
err = max(abs(w @ (d[a] * d[b]) - (2 / (2 * a + 1) if a == b else 0)) for a in range(2, 61, 7) for b in range(2, 61, 5))
check("Wigner-d orthogonality", err < 1e-12, f"max error {err:.1e}")

cl_ee = np.zeros(201)
cl_ee[2:] = 1 / np.arange(2, 201) ** 2.5
cl_aa = np.zeros(3)
cl_aa[0] = 0.7
bb = cb_bb_from_ee(cl_ee, cl_aa, 150)
r = bb[2:151] / (cl_ee[2:151] * 0.7 / np.pi)
check("isotropic limit", np.max(np.abs(r - 1)) < 1e-8, f"max |ratio-1| = {np.max(np.abs(r - 1)):.1e}")

t = spt_bb_templates()
cl_ee_camb = dl_to_cl(t["dl_ee_unlensed"])
lmax_e = len(cl_ee_camb) - 1
L0 = 30
cl_aa_cut = scale_invariant_cl_aa(lmax_e)
cl_aa_cut[:L0] = 0.0
cl_cb = cb_bb_from_ee(cl_ee_camb, cl_aa_cut, 2000)
Lg = np.unique(np.concatenate([np.linspace(L0 - 0.5, 100, 400), np.geomspace(100, 6000, 3000)]))
th = np.linspace(0, 2 * np.pi, 1024, endpoint=False)
for l in (200, 500, 1000, 2000):
    Lx = Lg[:, None] * np.cos(th)[None, :]
    Ly = Lg[:, None] * np.sin(th)[None, :]
    px, py = l - Lx, -Ly  # l along x
    lp = np.hypot(px, py)
    cee = np.interp(lp, np.arange(lmax_e + 1), cl_ee_camb, right=0.0)
    cos2 = np.cos(2 * np.arctan2(py, px)) ** 2
    caa = 2 * np.pi / (Lg * (Lg + 1))
    inner = (cee * cos2).mean(axis=1) * 2 * np.pi
    flat = 4 * np.trapezoid(Lg * caa * inner, Lg) / (2 * np.pi) ** 2
    rel = cl_cb[l] / flat - 1
    check(f"full-sky vs flat-sky l={l}", abs(rel) < 0.02, f"full {cl_cb[l]:.4e} flat {flat:.4e} (rel {rel:+.3%})")

print(f"D_l^CB (A_CB=1) at l=80,250,450: {t['dl_cb'][80]:.3e} {t['dl_cb'][250]:.3e} {t['dl_cb'][450]:.3e} muK^2")
sys.exit(1 if fails else 0)
