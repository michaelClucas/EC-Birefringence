"""Task 031 numbers for EC-06a (Magueijo-Zlosnik). All measured inputs from data/.

N1  tracker history eps_S(z) = 2P/H, P = gamma H/3 - (gamma Lambda^2/9) int(a^3 dt)/(a^3 H^2) (check_dynamics T2a),
    flat LCDM + radiation background (torsion corrections O(gamma^2)); xi gamma needed for beta.
N2  GW170817 speed on the tracker (026: c_T - 1 = 0.057 gamma^2 at z ~ 0.01, upper side 7e-16): xi >= 1.4e4 reproduced.
N3  tensor singular scale k* = H/eps_S (check_dynamics T6): the paper's own requirement k* >> k_LIGO
    (2001.06373 Eqs. 69-70) at f = 24 Hz gives eps_S(0) < H0/k_GW.
N4  off-tracker: general solution eps = eps_tr + eps_h0 (a H)_0^3/(a H)^3; ratio eps_h(z*)/eps_h(0);
    cancellation branch eps(0) = 0; GW path z <= z_GW; N_eff bound on gamma.
N5  strong coupling from theta (Lambda-field) exchange with a Stueckelberg photon:
    Lambda_SC = (16 pi)^(1/4) (2 sqrt3 M_P m/(xi (1 - c_s)))^(1/2); m needed for 10 TeV; 029's dim-7 arithmetic.
"""
import math

import numpy as np
from scipy.integrate import quad

from common import check, load, summary

cos = {k: float(v["value"]) for k, v in load("cosmo_inputs.csv").items()}
cos_err = {k: float(v["err"]) if v["err"] else 0.0 for k, v in load("cosmo_inputs.csv").items()}
th = {k: float(v["value"]) for k, v in load("ec06_theory_inputs.csv").items() if k.startswith("cT")}
mzi = {k: float(v["value"]) for k, v in load("mz_inputs.csv").items()}
ecc = {k: float(v["value"]) for k, v in load("ec06c_costs_inputs.csv").items() if k.startswith("m_gamma")}
meas = load("birefringence_measurements.csv")
beta = math.radians(float(meas["DPK25-beta-baseline"]["value"]))
MP = cos["MP_reduced_eV"]
hbar = cos["hbar_eV_s"]
h = cos["H0"] / 100
H0 = cos["H0"] / (cos["Mpc_km"]) / 1.0 * hbar  # km/s/Mpc -> s^-1 -> eV
Om = cos["Omega_m"]
Or = 2.469e-5 / h**2 * (cos["T_cmb_K"] / 2.7255) ** 4 * (1 + 0.2271 * cos["N_eff_SM"])
OL = 1 - Om - Or
zs = cos["z_star"]
a_s = 1 / (1 + zs)
print(f"inputs: beta = {math.degrees(beta)} deg, H0 = {H0:.4e} eV, Omega_m = {Om}, Omega_r = {Or:.3e}, z* = {zs}")


def E(a):
    return math.sqrt(Om / a**3 + Or / a**4 + OL)


# N1: tracker. Units H0 = 1. Lambda = 3 OL. int a^3 dt = int a^2 da/E (from a = 0).
Lam = 3 * OL


def I3(a):
    return quad(lambda b: b**2 / E(b), 0, a, limit=200)[0]


def eps_tr(a):  # per unit gamma
    return 2 / 3 - 2 / 9 * Lam**2 * I3(a) / (a**3 * E(a)**3)


eps0 = eps_tr(1.0)
eps_star = eps_tr(a_s)
J = quad(lambda a: eps_tr(a) / a, a_s, 1, limit=200)[0]  # int eps dz/(1+z) per gamma
xig = 2 * beta / J
const_eps_xieps = 2 * beta / math.log(1 + zs)
print(f"N1 tracker: eps_S/gamma = {eps_star:.4f} at z*, {eps0:.4f} today; int eps dz/(1+z) = {J:.3f} gamma "
      f"(constant 2/3 would give {2 / 3 * math.log(1 + zs):.3f}); xi*gamma needed = {xig:.3e}; "
      f"constant-eps xi*eps = {const_eps_xieps:.3e}")
check("N1 tracker: eps_S = 2 gamma/3 in radiation and matter eras, reduced by the Lambda term today",
      abs(eps_star - 2 / 3) < 1e-3 and eps0 < 2 / 3)

# N2: GW speed (tracker relation of 026): c_T - 1 = 0.057 gamma^2 (2001.06373 Eq. 65 at z = 0.01)
g_gw = math.sqrt(th["cT_GW170817_hi"] / 0.057)
xi_gw = xig / g_gw
print(f"N2 GW170817 speed: gamma <= {g_gw:.2e}, eps_S(z*) <= {2 / 3 * g_gw:.2e}, xi >= {xi_gw:.2e} (026: 1.4e4)")
check("N2 reproduces 026's xi >= 1.4e4 within 15% (Lambda-era eps reduction included here)",
      1.2e4 < xi_gw < 1.8e4)

# N3: k* condition
kGW = 2 * math.pi * mzi["GW170817_f_low"] * hbar
r = H0 / kGW
g_k = r / eps0
xi_k = xig / g_k
print(f"N3 k_GW(24 Hz) = {kGW:.3e} eV, H0/k_GW = {r:.2e}. Tracker at gamma = {g_gw:.1e}: k_GW/k* = "
      f"{eps0 * g_gw / r:.1e}. k* > k_GW needs eps_S(0) < {r:.1e}: gamma < {g_k:.1e}, xi >= {xi_k:.1e}")
check("N3 at the GW170817 speed bound the LIGO band sits >= 1e12 above k* (Eq. 72 used outside its regime)",
      eps0 * g_gw / r > 1e12)
check("N3 k* > k_GW on the tracker needs xi >= 1e16", xi_k > 1e16)

# N4: off-tracker
ratio_h = (1.0 * E(1.0))**3 / (a_s * E(a_s))**3  # eps_h(z*)/eps_h(0) = (aH)_0^3/(aH)_*^3
print(f"N4 homogeneous mode: eps_h(z*)/eps_h(0) = {ratio_h:.2e} (it grows toward today; no large-at-z* decaying branch)")
check("N4a homogeneous eps_h(z*)/eps_h(0) < 1e-3: off-tracker modes cannot make eps(z*) >> eps(0)", ratio_h < 1e-3)
zGW = cos["H0"] * mzi["GW170817_dL"] / 299792.458


def eps_cancel(a):  # per gamma, eps(0) = 0
    return eps_tr(a) - eps0 * (1.0 * E(1.0))**3 / (a * E(a))**3


Jc = quad(lambda a: eps_cancel(a) / a, a_s, 1, limit=200)[0]
zz = np.linspace(0, zGW, 21)
emax = max(abs(eps_cancel(1 / (1 + q))) for q in zz)
print(f"N4 cancellation branch (C tuned so eps(0) = 0): z_GW = {zGW:.4f}, max |eps|/gamma on z <= z_GW = {emax:.2e}, "
      f"beta integral {Jc:.3f} gamma (tracker {J:.3f})")
# N_eff bound on gamma (radiation-era deficit gamma^2/9, as in task 022)
f_nu = 0.2271
Neff_lo = cos["N_eff_obs"] - 2 * cos_err["N_eff_obs"]
dmax = f_nu * (cos["N_eff_SM"] - Neff_lo) / (1 + f_nu * cos["N_eff_SM"])
g_neff = math.sqrt(9 * dmax)
xi_best = 2 * beta / (Jc * g_neff)
g_k_c = r / emax
xi_k_c = 2 * beta / (Jc * g_k_c)
print(f"N4 N_eff >= {Neff_lo:.2f}: deficit <= {dmax:.3f}, gamma <= {g_neff:.2f}; best case if no GW bound applied: "
      f"xi >= {xi_best:.1e}. With the k* condition on z <= z_GW: gamma < {g_k_c:.1e}, xi >= {xi_k_c:.1e}")
check("N4b cancellation branch: with only N_eff, xi ~ 1e-3; with the k* condition, xi >= 1e14", 1e-4 < xi_best < 1e-2
      and xi_k_c > 1e14)

# N5: strong coupling
one_m_cs = Om / (Om + OL)  # cosmic mean today; ~1 in galaxies, in matter, and at z* (smallest value used)
m_list = [ecc["m_gamma_Proca_gal"], 1e-24, ecc["m_gamma_PDG"]]


def lam_sc(xi, m, omc):
    return (16 * math.pi)**0.25 * math.sqrt(2 * math.sqrt(3) * MP * m / (xi * omc))


for xi in (xi_best, xi_gw, xi_k):
    s = ", ".join(f"m = {m:.0e}: {lam_sc(xi, m, one_m_cs):.2e} eV" for m in m_list)
    print(f"N5 xi = {xi:.1e}: Lambda_SC (1 - c_s = {one_m_cs:.3f}) {s}")
E_safe = mzi["E_collider_safe"]
m_need = (E_safe / (16 * math.pi)**0.25)**2 * xi_gw * one_m_cs / (2 * math.sqrt(3) * MP)
m_need_best = (E_safe / (16 * math.pi)**0.25)**2 * xi_best * one_m_cs / (2 * math.sqrt(3) * MP)
print(f"N5 m needed for Lambda_SC = {E_safe:.0e} eV: {m_need:.1e} eV at xi = {xi_gw:.1e}; {m_need_best:.1e} eV at "
      f"xi = {xi_best:.1e} (PDG m < {ecc['m_gamma_PDG']:.0e} eV): {math.log10(m_need / ecc['m_gamma_PDG']):.1f} and "
      f"{math.log10(m_need_best / ecc['m_gamma_PDG']):.1f} orders above PDG")
check("N5a Lambda_SC < 10 keV for xi >= 1.4e4 at every m up to the PDG limit", lam_sc(xi_gw, ecc["m_gamma_PDG"], one_m_cs) < 1e4)
check("N5b even the best off-tracker case (xi ~ 1e-3, no GW bound) has Lambda_SC < 10 MeV at the PDG limit",
      lam_sc(xi_best, ecc["m_gamma_PDG"], one_m_cs) < 1e7)
check("N5c collider-safe Lambda_SC needs m >= 1e10 x PDG in every branch", min(m_need, m_need_best) > 1e10 * ecc["m_gamma_PDG"])
for m in (ecc["m_gamma_Proca_gal"], ecc["m_gamma_PDG"]):
    print(f"N5 029's dim-7 arithmetic (m M_P^2/xi)^(1/3) at xi = {xi_gw:.1e}, m = {m:.0e}: "
          f"{(m * MP**2 / xi_gw)**(1 / 3) / 1e9:.3g} GeV")
check("N5d 029's ~0.1 GeV reproduced as the dim-7 scale at m = 1e-26 eV, but the theta-exchange scale is lower",
      0.05e9 < (1e-26 * MP**2 / xi_gw)**(1 / 3) < 0.3e9 and lam_sc(xi_gw, 1e-26, one_m_cs) < 1e3)
g_alp = xi_gw / (math.sqrt(3) * MP) * 1e9
print(f"N5 gauge-invariant (gradient) piece = ALP coupling of canonical phi: g = xi/(sqrt3 M_P) = {g_alp:.1e} GeV^-1 "
      f"at xi = {xi_gw:.1e}")
summary()
