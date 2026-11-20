"""Skeptic 028: items 3 and 4, independent numbers.

W1  Lower edge of the Stueckelberg window from photons leaving compact objects (same mechanism as 025's
    solar-surface limit: moving-source term |dp| ~ xi v d^2Phi, conversion amplitude A ~ int|dp|dl/(4m)),
    for the Sun, a typical WD, Sirius B and a NICER neutron star; criteria A < 0.1 (025's choice) and A < 1.
W2  Plasma: conversion is vacuum-like only where the photon-pi detuning dk = w_p^2/(2w) satisfies dk r < 1.
    Sun: Baumbach-Allen corona, find r_vac for optical photons, integrate the mixing from there out.
    WD: hydrostatic exosphere scale height kT/(m_p g) vs R.
W3  Windows under the four choices (A < 0.1 or A < 1) x (upper edge Proca galactic 1e-26 eV or PDG 1e-18 eV).
C1  025's dim-7 gravitational-exchange cutoff (m M_P^2/xi)^(1/3), with and without a (4 pi)^(1/3) NDA factor.
C2  Lorentz-mode (eta) exchange: eta stiffness (from eta_sector.py) L = 4 M_P^2 f_TT H^2 (d_x eta)^2 =
    (2/3) alpha M_P^2 (d_x eta)^2 with alpha = f_TT Tbar (Tbar = 6H^2); pi couples to eta through
    -(xi/2m) pi (eta_xx - eta_tt) Ft^{tau x}. Integrating eta out gives a pi pi gamma gamma contact term with
    amplitude M ~ (3 xi^2/(8 alpha M_P^2 m^2)) E^4; partial-wave unitarity |M| < 8 pi gives
    Lambda_eta = (64 pi alpha M_P^2 m^2/(3 xi^2))^(1/4). Power counting, O(1) factors not tracked.
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
npass = nfail = 0


def check(name, ok, msg=""):
    global npass, nfail
    npass += bool(ok)
    nfail += (not ok)
    print(f"{'PASS' if ok else 'FAIL'}: {name} {msg}".rstrip())


def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        return {r["id"]: r for r in csv.DictReader(f)}


def val(d, k):
    return float(d[k]["value"])


bi = load(os.path.join(ROOT, "data", "bounds_inputs.csv"))
ci = load(os.path.join(ROOT, "data", "cosmo_inputs.csv"))
ec = load(os.path.join(ROOT, "data", "ec06c_costs_inputs.csv"))
meas = load(os.path.join(ROOT, "data", "birefringence_measurements.csv"))
loc = load(os.path.join(HERE, "inputs_028.csv"))

hbarc = val(bi, "hbarc")
c = val(bi, "c_kms") * 1e3
MP = val(bi, "MP_reduced")
v = val(bi, "v_sun_cmb") * 1e3 / c
zstar = val(ci, "z_star")
beta = math.radians(val(meas, "DPK25-beta-baseline"))
xi = beta / (1.5 * math.log(1 + zstar))
GMs, Rs = val(ec, "GM_sun"), val(ec, "R_sun")
m_proca, m_pdg = val(ec, "m_gamma_Proca_gal"), val(ec, "m_gamma_PDG")
print(f"xi = {xi:.3e}, v/c = {v:.3e}")

objs = {
    "Sun": (GMs, Rs),
    "WD 0.6 Msun": (GMs * val(loc, "M_WD_typ"), Rs * val(loc, "R_WD_typ")),
    "Sirius B": (GMs * val(loc, "M_SiriusB"), Rs * val(loc, "R_SiriusB")),
    "NS J0030": (GMs * val(loc, "M_J0030"), val(loc, "R_J0030") * 1e3),
}


def mmin(GM, r0, Amax):
    """025's estimate int|dp|dl ~ xi v Phi(r0)/r0 (the exact integral of xi v GM/r^3 from r0 is half this)."""
    Phi = GM / (r0 * c**2)
    return xi * v * Phi / r0 * hbarc / (4 * Amax)


res = {}
for name, (GM, R) in objs.items():
    res[name] = {A: mmin(GM, R, A) for A in (0.1, 1.0)}
    print(f"W1 {name:12s}: Phi = {GM / (R * c**2):.2e}, R = {R:.3e} m -> m_min = {res[name][0.1]:.1e} eV (A<0.1), "
          f"{res[name][1.0]:.1e} eV (A<1)")
check("W1a reproduces 025's solar value 6.6e-28 eV (A < 0.1)", abs(res["Sun"][0.1] / 6.6e-28 - 1) < 0.05)
check("W1b white dwarfs give a stronger lower edge than the Sun by > 100x", res["WD 0.6 Msun"][0.1] > 100 * res["Sun"][0.1])
check("W1c WD lower edge exceeds the Proca galactic upper limit 1e-26 eV even with A < 1",
      res["WD 0.6 Msun"][1.0] > m_proca and res["Sirius B"][1.0] > m_proca)

# W2 plasma
w_opt = 2.0  # eV, optical photon (choice)
wp2_per_ne = (3.7e-11)**2  # eV^2 per cm^-3: w_p = 3.7e-11 eV sqrt(n_e / cm^-3)


def ne_BA(xr):
    return val(loc, "ne_corona_model") * (2.99 * xr**-16 + 1.55 * xr**-6 + 0.036 * xr**-1.5)


xr = 1.0
while ne_BA(xr) * wp2_per_ne / (2 * w_opt) * (xr * Rs / hbarc) >= 1:
    xr *= 1.01
print(f"W2 Sun, optical: vacuum-like mixing starts at r = {xr:.1f} R_sun (n_e = {ne_BA(xr):.1e} cm^-3)")
m_sun_plasma = {A: mmin(GMs, xr * Rs, A) for A in (0.1, 1.0)}
print(f"W2 Sun with plasma (optical): m_min = {m_sun_plasma[0.1]:.1e} eV (A<0.1), {m_sun_plasma[1.0]:.1e} eV (A<1)")
check("W2a plasma lowers the solar edge but does not remove it (factor < 1e3)",
      1e-3 * res["Sun"][0.1] < m_sun_plasma[0.1] < res["Sun"][0.1])
GM, R = objs["WD 0.6 Msun"]
ne_thr = 2 * w_opt / (wp2_per_ne * R / hbarc)
kB, mp = 1.380649e-23, 1.67262192e-27  # CODATA exact / 2018
Hsc = kB * val(loc, "Teff_SiriusB") / (mp * GM / R**2)
nfold = math.log(1e17 / ne_thr)  # from a photospheric n_e ~ 1e17 cm^-3 (order of magnitude, choice)
print(f"W2 WD: vacuum-like needs n_e < {ne_thr:.1e} cm^-3; scale height {Hsc:.0f} m; reached after {nfold:.0f} "
      f"scale heights = {nfold * Hsc / 1e3:.1f} km above the photosphere (R = {R / 1e3:.0f} km)")
check("W2b WD: the plasma-detuned layer is < 1% of R, so vacuum mixing applies", nfold * Hsc < 0.01 * R)

print("W3 windows (lower edge = max over Sun-with-plasma, WD, Sirius B; NS not used):")
for A in (0.1, 1.0):
    lo = max(m_sun_plasma[A], res["WD 0.6 Msun"][A], res["Sirius B"][A])
    for up_name, up in (("Proca gal", m_proca), ("PDG", m_pdg)):
        state = f"open {lo:.0e} - {up:.0e} eV" if lo < up else f"CLOSED (lower {lo:.0e} > upper {up:.0e})"
        print(f"   A<{A}: {up_name:9s}: {state}")
lo01 = max(res["WD 0.6 Msun"][0.1], res["Sirius B"][0.1])
check("W3 with the Proca galactic limit the window is closed; with PDG only it is open above ~1e-24 eV",
      lo01 > m_proca and lo01 < m_pdg)
print(f"W3 neutron star (vacuum mixing, magnetosphere ignored): m_min = {res['NS J0030'][0.1]:.1e} eV (A<0.1)")

for m in (7e-28, m_proca, lo01, m_pdg):
    L = (m * MP**2 / xi)**(1 / 3)
    print(f"C1 m = {m:.0e} eV: (m M_P^2/xi)^(1/3) = {L / 1e9:.3g} GeV, x (4 pi)^(1/3) = "
          f"{L * (4 * math.pi)**(1 / 3) / 1e9:.3g} GeV")
check("C1 reproduces 025's 20-60 GeV at 7e-28 and 1e-26 eV",
      19e9 < (7e-28 * MP**2 / xi)**(1 / 3) < 25e9 and 50e9 < (1e-26 * MP**2 / xi)**(1 / 3) < 60e9)


def lam_eta(m, alpha):
    return (64 * math.pi * alpha * MP**2 * m**2 / (3 * xi**2))**0.25


for m in (7e-28, m_proca, lo01, m_pdg):
    print(f"C2 m = {m:.0e} eV: Lambda_eta = {lam_eta(m, 1):.2e} eV (alpha = 1), {lam_eta(m, 0.01):.2e} eV (alpha = 0.01)")
check("C2a Lambda_eta < 10 MeV for every m up to the PDG limit and alpha <= 1", lam_eta(m_pdg, 1.0) < 1e7)
check("C2b Lambda_eta < 1 keV everywhere in 025's window (m <= 1e-26 eV, alpha <= 1)", lam_eta(m_proca, 1.0) < 1e3)
m_need = (1e13)**2 * xi / (MP * math.sqrt(64 * math.pi / 3))
print(f"C2 m needed for Lambda_eta >= 10 TeV at alpha = 1: {m_need:.1e} eV (vs PDG {m_pdg:.0e} eV)")
check("C2c a collider-safe Lambda_eta needs m far above the PDG photon-mass limit", m_need > 1e3 * m_pdg)
print(f"== {npass} PASS, {nfail} FAIL ==")
