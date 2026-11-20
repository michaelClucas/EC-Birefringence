"""Skeptic 026: independent numbers for EC-06 (beta formula, reionization share, f(T) xi, ln Q map,
MZ GW-speed bound, 1003.4823 k_AF conversion, A_CB order of magnitude). Written without importing 022/023 code.
Run from the repo root."""
import csv, math
import numpy as np
from scipy.integrate import quad

def load(path, key="id"):
    with open(path, newline="", encoding="utf-8") as f:
        return {r[key]: r for r in csv.DictReader(f)}

C = load("data/cosmo_inputs.csv"); PB = load("data/planck2018_bestfit.csv")
TH = load("data/ec06_theory_inputs.csv"); CO = load("data/ec06_cosmo_inputs.csv")
BM = load("data/birefringence_measurements.csv"); LZ = load("data/kaf_lowz_bounds.csv")
v = lambda d, k: float(d[k]["value"])

H0kms = v(C, "H0"); h = H0kms / 100; Om = v(C, "Omega_m"); zs = v(C, "z_star")
Tcmb = v(C, "T_cmb_K"); kB = v(C, "k_B_eV_K"); hbar_s = v(C, "hbar_eV_s"); hbarc = v(C, "hbar_c_eV_m")
Mpc_m = v(C, "Mpc_km") * 1e3; MP = v(C, "MP_reduced_eV"); Neff = v(C, "N_eff_SM")
beta_deg = float(BM["DPK25-beta-baseline"]["value"]) if "value" in BM["DPK25-beta-baseline"] else None
if beta_deg is None:
    beta_deg = float(list(BM["DPK25-beta-baseline"].values())[2])
beta = math.radians(beta_deg)
H0_eV = H0kms * 1e3 / Mpc_m * hbar_s                    # H0 in eV
rho_crit = 3 * H0_eV**2 * MP**2                          # eV^4
rho_g = math.pi**2 / 15 * (kB * Tcmb)**4
Og = rho_g / rho_crit; Orad = Og * (1 + 7/8 * (4/11)**(4/3) * Neff); OL = 1 - Om - Orad
E = lambda z: math.sqrt(Om*(1+z)**3 + Orad*(1+z)**4 + OL)
out = []
def rep(name, val, note=""):
    s = f"{name}: {val}" + (f"   [{note}]" if note else ""); print(s); out.append(s)

rep("ACT beta (deg, from data)", beta_deg)
# 1. beta for p_0hat = xi eps H: beta = 1/2 int p_0hat dt = 1/2 xi eps int H dt = (xi eps/2) int dz/(1+z)
num = quad(lambda z: 1/(1+z), 0, zs, epsrel=1e-12)[0]
# independent route: integrate in cosmic time t with H(t) on a grid of a, dt = da/(aH)
a = np.logspace(math.log10(1/(1+zs)), 0, 200001)
Ha = np.array([E(1/x-1) for x in a]); t_int = np.trapezoid(Ha / (a*Ha), a)  # int H dt
rep("int_0^z* dz/(1+z)", num, "analytic ln(1+z*) = %.6f" % math.log(1+zs))
rep("int H dt on a-grid", t_int)
xieps = 2*beta/math.log(1+zs); rep("xi*eps for ACT beta (p tracks H)", f"{xieps:.4e}")
# 2. reionization share
zre = v(CO, "z_re"); dz = v(CO, "dz_re")
rep("point share ln(1+z_re)/ln(1+z*)", round(math.log(1+zre)/math.log(1+zs), 4))
rep("share at z=10", round(math.log(11)/math.log(1+zs), 4))
# visibility-weighted, hydrogen-only tanh in y=(1+z)^1.5 (He omitted: sensitivity only)
yre = (1+zre)**1.5; dy = 1.5*(1+zre)**0.5*dz
xe = lambda z: 0.5*(1+math.tanh((yre-(1+z)**1.5)/dy))
dtau = lambda z: xe(z)*(1+z)**2/E(z)        # shape only; e^-tau ~ 1 within 5%
w = quad(dtau, 0, 40, limit=400, points=[zre])[0]
wavg = quad(lambda z: dtau(z)*math.log(1+z), 0, 40, limit=400, points=[zre])[0]/w
rep("visibility-avg share (H only, no He)", round(wavg/math.log(1+zs), 4))
# add HeII at z=3.5 (step, fHe ~ 0.08 model setting) as sensitivity
fHe = 0.08
xe2 = lambda z: (1+fHe)*xe(z) + fHe*0.5*(1+math.tanh((3.5-z)/dz))
d2 = lambda z: xe2(z)*(1+z)**2/E(z)
w2 = quad(d2, 0, 40, limit=400, points=[3.5, zre])[0]
rep("visibility-avg share (H+He, fHe=0.08 sensitivity)", round(quad(lambda z: d2(z)*math.log(1+z), 0, 40, limit=400, points=[3.5, zre])[0]/w2/math.log(1+zs), 4))
# 3. f(T): T_t = 3H exactly -> beta = (3 xi/2) ln(1+z*)
xiT = 2*beta/(3*math.log(1+zs)); rep("f(T) xi (eps_T = 3)", f"{xiT:.4e}")
# 4. ln Q (Q = 6H^2): p = c d ln Q -> beta per unit = (1/2) * 2 ln E(z*) = ln E(z*); trace per unit xi = 1.5 ln(1+z*)
rep("trace beta/xi = 1.5 ln(1+z*)", round(1.5*math.log(1+zs), 4))
rep("lnQ  beta per unit = ln E(z*)", round(math.log(E(zs)), 4), "same coefficient normalisation as 022 N7")
rep("ratio trace/lnQ - 1", round(1.5*math.log(1+zs)/math.log(E(zs)) - 1, 4))
rep("lnQ post-reion share ln E(z_re)/ln E(z*)", round(math.log(E(zre))/math.log(E(zs)), 4))
# 5. MZ GW speed. Matter tracker: omega^2 = [1+(1+3w)g^2/9]k^2 -> c_T-1 = g^2/18 (w=0) > 0: GW superluminal.
lo = float(TH["cT_GW170817_lo"]["value"]); hi = float(TH["cT_GW170817_hi"]["value"])
g_hi = math.sqrt(18*hi); g_lo = math.sqrt(18*abs(lo))
rep("gamma max using +7e-16 side (correct side, c_T>1)", f"{g_hi:.3e}", "eps_S = 2g/3 = %.2e, xi >= %.2e" % (2*g_hi/3, xieps/(2*g_hi/3)))
rep("gamma max using -3e-15 side (as in 022)", f"{g_lo:.3e}", "eps_S = %.2e, xi >= %.2e" % (2*g_lo/3, xieps/(2*g_lo/3)))
# Lambda-era first term of 2001.06373 Eq. 65 with p=0: 1 + r(1-r)/(1+r)^2 g^2, r = kappa rho/Lambda (assumed Lambda = 3H0^2 Omega_L)
for z in (0.0, 0.0098):
    r = Om*(1+z)**3/OL; cf = r*(1-r)/(1+r)**2/2
    rep(f"Eq.65 gamma^2 coefficient of c_T-1 at z={z}", round(cf, 4), "vs 1/18 = 0.0556 (normalisation of Lambda assumed)")
# 6. 1003.4823 conversion: constant k_AF: angle = k_AF * t_lookback
zeff = float(LZ["L-AFG10-UV"]["z_eff"])
tL = quad(lambda z: 1/((1+z)*E(z)), 0, zeff)[0] / H0_eV
ang = math.radians(abs(float(LZ["L-AFG10-UV"]["central"])) + 2*float(LZ["L-AFG10-UV"]["err_1sigma_or_limit"]))
rep("t_L(z_eff) eV^-1", f"{tL:.4e}", "file: 5.508e32")
rep("k_AF limit eV", f"{ang/tL:.3e}", "file: 1.65e-34; p0 = 2k")
rep("xi eps limit (p tracks H): 2 alpha/ln(1+z)", round(2*ang/math.log(1+zeff), 4), "file 0.136")
# 7. A_CB, kinematic assumption, BBKS P(k) normalised to sigma8 (independent of 023's EH code)
ns = v(PB, "ns"); s8 = v(CO, "sigma8"); Ob = v(C, "Omega_b_h2")/h**2
H0c = h/2997.92458            # Mpc^-1
Gam = Om*h*math.exp(-Ob*(1+math.sqrt(2*h)/Om))
def T(k):
    q = k/(h*Gam); return math.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**-0.25
def W8(x): return 3*(math.sin(x)-x*math.cos(x))/x**3
Pu = lambda k: k**ns*T(k)**2
sig2 = quad(lambda lk: math.exp(lk)**3*Pu(math.exp(lk))*W8(math.exp(lk)*8/h)**2/(2*math.pi**2), math.log(1e-5), math.log(50), limit=400)[0]
A = s8**2/sig2; P0 = lambda k: A*Pu(k)
Ela = lambda aa: math.sqrt(Om/aa**3 + 1 - Om)
def Dg(aa): return 2.5*Om*Ela(aa)*quad(lambda y: 1/(y*Ela(y))**3, 0, aa)[0]
D0 = Dg(1.0)
zg = np.linspace(0.01, 20, 400)
chi = np.array([quad(lambda z: 1/E(z), 0, zz)[0]/H0c for zz in zg])
Dz = np.array([Dg(1/(1+zz))/D0 for zz in zg])
fz = np.array([Om*(1+zz)**3/E(zz)**2 for zz in zg])**0.55
Hz = np.array([E(zz) for zz in zg])*H0c
xe_b = 2*beta/math.log(1+zs)
Wz = -(xe_b/6)*fz*Hz/(1+zg)*Dz     # delta beta = int dchi W delta_0
dchi = np.gradient(chi)
best = (0, 0)
for L in (2, 10, 100, 300, 770, 1500):
    CL = np.sum(dchi*Wz**2/chi**2*np.array([P0((L+0.5)/c) for c in chi]))
    DL = L*(L+1)*CL/(2*math.pi); rep(f"A_CB kinematic D_L at L={L}", f"{DL:.2e}")

