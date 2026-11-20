"""EC-06c repair costs in numbers (task 025). All measured inputs from data/.

N1  xi for the flat f(T) trace (eps_T = 3): beta = (3 xi/2) ln(1+z*)  (022 value 3.6e-4).
N2  Size of the U(1) breaking dp = xi dT on f(T) solutions (check_ft_lorentz.py L2f):
      static potential, cosmic frame:  |dT|_phys = |C(z)| H |grad Phi|,  C = -3 + 3 Hdot/H^2  (-7.5 in matter era)
      source moving at v w.r.t. the cosmic frame:  |dT|_phys ~ v |grad grad Phi| ~ v Phi / r^2
    Local sources (Sun, Earth, Galaxy) and cosmological perturbations.
N3  Repair (a) Stueckelberg photon: converted amplitude A ~ int|dp| dl/(4m) (check_stueckelberg.py S3).
    Lower limit on m from A < 0.1 (choice), from CFJ stability m > p/2, vs upper limits (PDG; Proca galactic).
    Naive EFT cutoff of the pi-matter-photon contact operator: Lambda ~ (m M_P^2/xi)^(1/3).
N4  Repair (b) gravitational CS (0907.5159): their consistency condition q RR~ = 16 pi G p FF~ is a field
    constraint, and it addresses diffeomorphism (stress-tensor) consistency for an external p, not U(1).
    Size check only: q ~ kappa^2 p vs LAGEOS 1/m_cs (m_cs = 2/q in their normalization) and GW kappa.
N5  Repair (c): (c1) ln(T_sc) completion: observer-endpoint dependence (ground vs L2) from L4 (T_sc = 6H^2 +
    2|grad Phi|^2); (c2) auxiliary scalar with -(M^2/2)(d theta - T)^2: M from CAST and from inflationary
    A_CB (A_CB = (xi H_inf/(4 pi M))^2).
"""
import math

from common import check, load, summary

bi = load("bounds_inputs.csv")
ci = load("cosmo_inputs.csv")
inp = load("ec06c_costs_inputs.csv")
bnd = load("bounds.csv")
meas = load("birefringence_measurements.csv")

hbarc = float(bi["hbarc"]["value"])  # eV m
c_kms = float(bi["c_kms"]["value"])
c_ms = c_kms * 1e3
Mpc = float(bi["Mpc"]["value"])
MP = float(bi["MP_reduced"]["value"])
H0_kms = float(ci["H0"]["value"])
Om = float(ci["Omega_m"]["value"])
zstar = float(ci["z_star"]["value"])
H0 = H0_kms / (Mpc / 1e3) * hbarc / c_ms  # s^-1 -> eV:  H[1/s] * hbar ; hbar = hbarc/c
H0_inv_m = H0_kms / (Mpc / 1e3) / c_ms  # H0/c in 1/m
beta = math.radians(float(meas["DPK25-beta-baseline"]["value"]))

# ---------------------------------------------------------------- N1
xi = beta / (1.5 * math.log(1 + zstar))
print(f"N1 xi (flat f(T), eps_T = 3) = {xi:.3e}")
check("N1 xi = 3.6e-4 (022)", abs(xi / 3.6e-4 - 1) < 0.02)


def Hz(z):  # flat LCDM, matter + Lambda (radiation omitted: H(z*) low by ~1.5 -> conservative for m_min)
    return H0 * math.sqrt(Om * (1 + z)**3 + 1 - Om)


def C_of(z):
    E2 = Om * (1 + z)**3 + 1 - Om
    hdot_over_h2 = -1.5 * Om * (1 + z)**3 / E2
    return -3 + 3 * hdot_over_h2


check("N2a C(matter era) = -7.5 (matches L2g: 15/tau = 7.5 aH)", abs(C_of(1e4) + 7.5) < 1e-3)
C0 = C_of(0)
print(f"N2 C(0) = {C0:.3f}, C(z*) = {C_of(zstar):.3f}")

# ---------------------------------------------------------------- N2 local sources
GMs, Rs = float(inp["GM_sun"]["value"]), float(inp["R_sun"]["value"])
GMe, Re = float(inp["GM_earth"]["value"]), float(inp["R_earth"]["value"])
au = float(inp["au"]["value"])
vc = float(inp["v_circ_MW"]["value"]) * 1e3
R0 = float(inp["R0_MW"]["value"]) * 1e3 * Mpc / 1e6
v_sun = float(bi["v_sun_cmb"]["value"]) * 1e3 / c_ms  # in units of c
src = {  # name: (Phi, r [m])
    "Sun surface": (GMs / (Rs * c_ms**2), Rs),
    "Sun at 1 au": (GMs / (au * c_ms**2), au),
    "Earth surface": (GMe / (Re * c_ms**2), Re),
    "Galaxy at Sun": (vc**2 / c_ms**2, R0),
}
A_max = 0.1  # choice: converted amplitude allowed
m_min_local = 0.0
rows = []
for name, (Phi, r) in src.items():
    gradPhi = Phi / r * hbarc  # eV (1/m -> eV)
    dp_static = xi * abs(C0) * H0 * gradPhi  # eV^2
    dp_moving = xi * v_sun * Phi / r**2 * hbarc**2  # eV^2
    # path integral over the source scale r: int |dp| dl ~ |dp| r
    I_static = dp_static * r / hbarc
    I_moving = dp_moving * r / hbarc
    mmin = max(I_static, I_moving) / (4 * A_max)
    m_min_local = max(m_min_local, mmin)
    rows.append((name, Phi, gradPhi, dp_static, dp_moving, mmin))
    print(f"N2 {name:14s}: Phi = {Phi:.2e}, |grad Phi| = {gradPhi:.2e} eV, |dp| static = {dp_static:.1e} eV^2, "
          f"moving (v = {v_sun * c_kms:.0f} km/s) = {dp_moving:.1e} eV^2 -> m_min(A<0.1) = {mmin:.1e} eV")
check("N2b moving-source term dominates the static (C H) term for every local source", all(r_[4] > r_[3] for r_ in rows))

# cosmological perturbations: per coherence patch a/k, amplitude xi |C| H Phi /(4m); random walk over N = k/(aH)
lnA = float(load("planck2018_bestfit.csv")["logA"]["value"])
As = math.exp(lnA) * 1e-10
Phi_prim = 0.6 * math.sqrt(As)
N_patch = 1e4  # generous choice (k/aH up to 1e4)
m_min_cos = 0.0
for z in (0.0, zstar):
    amp_m = xi * abs(C_of(z)) * Hz(z) * Phi_prim * math.sqrt(N_patch) / 4
    m_min_cos = max(m_min_cos, amp_m / A_max)
    print(f"N2 cosmological z = {z:6.0f}: H = {Hz(z):.2e} eV, Phi = {Phi_prim:.1e}: A*m = {amp_m:.1e} eV")
print(f"N3 m_min from conversion: local {m_min_local:.1e} eV, cosmological {m_min_cos:.1e} eV (vacuum, degenerate mixing)")

# CFJ stability: m > p/2 with p = xi eps H, eps = 3
m_stab_now, m_stab_star = xi * 3 * H0 / 2, xi * 3 * Hz(zstar) / 2
print(f"N3 CFJ tachyon removal m > p/2: today {m_stab_now:.1e} eV, at z* {m_stab_star:.1e} eV")
m_up_pdg = float(inp["m_gamma_PDG"]["value"])
m_up_proca = float(inp["m_gamma_Proca_gal"]["value"])
m_low = max(m_min_local, m_min_cos, m_stab_star)
print(f"N3 window: {m_low:.1e} eV < m_gamma < {m_up_pdg:.0e} eV (PDG) or < {m_up_proca:.0e} eV (Proca galactic field)")
check("N3a a photon-mass window exists below the PDG limit", m_low < m_up_pdg)
check("N3b the window survives the Proca galactic-field limit 1e-26 eV", m_low < m_up_proca, f"(lower end {m_low:.1e})")
for m in (m_low, m_up_proca, m_up_pdg):
    Lam = (m * MP**2 / xi)**(1 / 3)
    print(f"N3 naive cutoff (m M_P^2/xi)^(1/3) at m = {m:.0e} eV: {Lam:.1e} eV")
# superhorizon timelike-CFJ growth is bounded by exp(int p/2 dt) = exp((3 xi/2) N_efolds)
print(f"N3 max growth of the m = 0 tachyonic band over 60 e-folds: exp({1.5 * xi * 60:.3f})")

# ---------------------------------------------------------------- N4 gravitational CS sizes
lageos = float(inp["LAGEOS_mcs_inv"]["value"]) * 1e3  # m
kappa2 = 1 / MP**2  # eV^-2 (8 pi G)
p0 = xi * 3 * H0
q_nat = 2 * kappa2 * p0 * hbarc  # 16 pi G p -> length in m (16 pi G = 2 kappa^2)
q_max = 2 * lageos  # m_cs = 2/q, 1/m_cs <~ 1000 km
print(f"N4 q ~ 16 pi G p = {q_nat:.1e} m vs LAGEOS q < {q_max:.1e} m: margin {q_max / q_nat:.1e}")
check("N4a an O(kappa^2 p) gravitational CS coefficient is far below LAGEOS/GP-B", q_nat < 1e-50 * q_max)
k_gw = 2 * math.pi * 100 / c_ms  # 1/m at 100 Hz
kap_est = q_nat * H0_inv_m * k_gw  # 1/m
kap_max = float(inp["GW_kappa_birefringence"]["value"]) / (1e3 * Mpc)
print(f"N4 GW amplitude birefringence estimate q H0 k = {kap_est:.1e} /m vs {kap_max:.1e} /m (2101.11153)")
check("N4b GW birefringence from that coefficient is negligible", kap_est < 1e-50 * kap_max)

# ---------------------------------------------------------------- N5 (c1) ln T endpoint dependence
T0 = 6 * H0**2
Tstar = 6 * Hz(zstar)**2
g_ground = GMe / Re**2 / c_ms**2 * hbarc  # eV
g_L2 = GMs / au**2 / c_ms**2 * hbarc
lam = {"cosmological observer": T0, "ground (ACT)": T0 + 2 * g_ground**2, "L2 (Planck)": T0 + 2 * g_L2**2}
base = math.log(T0 / Tstar)
for k_, v in lam.items():
    print(f"N5 ln(T_obs/T_*) {k_:22s} = {math.log(v / Tstar):+.2f}")
ratio = math.log(lam["ground (ACT)"] / Tstar) / math.log(lam["L2 (Planck)"] / Tstar)
print(f"N5 beta(ground)/beta(L2) = {ratio:.2f}; sign relative to the cosmological-observer value: "
      f"{'flipped' if ratio * base < 0 or math.log(lam['ground (ACT)'] / Tstar) * base < 0 else 'same'}")
check("N5a ln T completion: local field dominates the observer endpoint (|ln(T_ground/T0)| > |ln(T0/T*)|)",
      abs(math.log(lam["ground (ACT)"] / T0)) > abs(base))
check("N5b ln T completion: ground and L2 observers differ by more than a factor 2", ratio > 2 or ratio < 0.5)

# (c2) auxiliary scalar
g_cast = float(bnd["A-CAST"]["value_eV"].replace("<", "").strip())  # eV^-1
M_cast = xi / g_cast
acb = float(meas["BAL25-ACB-SPT3G"]["value"])
M_over_Hinf = xi / (4 * math.pi * math.sqrt(acb))
print(f"N5 (c2) M > xi/g_CAST = {M_cast:.1e} eV = {M_cast / 1e9:.1e} GeV; M/H_inf > {M_over_Hinf:.1e} (SPT-3G A_CB)")
print(f"N5 (c2) M^2 T.T correction to the teleparallel action relative to M_P^2: {(M_cast / MP)**2:.1e} at the CAST minimum")
check("N5c (c2) both lower limits on M are far below M_P (correction to gravity negligible)",
      M_cast < 1e-6 * MP)
summary()
