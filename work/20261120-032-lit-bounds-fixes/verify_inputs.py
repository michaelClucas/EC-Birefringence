"""Compare 028's recalled compact-star inputs (work/028/inputs_028.csv) with the opened values in
data/compact_star_inputs.csv, and run independent checks on the values that were not printed in a source."""
import math

from common import check, load, summary, val

new = load("compact_star_inputs.csv")
old = load("inputs_028.csv", folder="work/20261120-028-skeptic-ec06c-costs")
ec = load("ec06c_costs_inputs.csv")

print("id               028 recalled   opened        rel. diff")
for k in ("M_SiriusB", "R_SiriusB", "Teff_SiriusB", "M_WD_typ", "R_WD_typ", "M_J0030", "R_J0030", "ne_corona_model"):
    a, b = val(old, k), val(new, k)
    print(f"{k:16s} {a:<14.6g} {b:<13.6g} {(a - b) / b:+.3%}")
check("Sirius B mass, J0030 mass and radius, Baumbach-Allen normalisation as recalled",
      all(val(old, k) == val(new, k) for k in ("M_SiriusB", "M_J0030", "R_J0030", "ne_corona_model")))
dR = (val(old, "R_SiriusB") - val(new, "R_SiriusB")) / val(new, "R_SiriusB")
check("Sirius B radius 0.008098 (028) within 1 sigma of 1806.00061 Table 4 (but its arXiv ID and error bar are wrong)",
      abs(val(old, "R_SiriusB") - val(new, "R_SiriusB")) < val(new, "R_SiriusB", "err"), f"({dR:+.2%})")
# 028 W1: m_min ~ GM/R^2, so the Sirius B edge scales by (R_old/R_new)^2 at fixed M
f = (val(old, "R_SiriusB") / val(new, "R_SiriusB"))**2
print(f"   effect on 028 W1 Sirius B lower edge: x{f:.3f} with the opened radius (m_min ~ M/R^2)")
check("opened radius changes 028's Sirius B edge by < 5% (verdict unaffected)", abs(f - 1) < 0.05)

# Nauenberg (1972, ApJ 175 417; no arXiv) zero-temperature mass-radius relation, mu_e = 2 (theory input):
# R/R_sun = (0.0225/mu_e) sqrt(1 - (M/M3)^(4/3)) / (M/M3)^(1/3),  M3 = 5.816/mu_e^2 M_sun
mu_e = 2.0
M3 = 5.816 / mu_e**2


def R_nauenberg(M):
    x = M / M3
    return 0.0225 / mu_e * math.sqrt(1 - x**(4 / 3)) / x**(1 / 3)


Rsb = R_nauenberg(val(new, "M_SiriusB"))
print(f"   Nauenberg R(Sirius B mass) = {Rsb:.5f} R_sun vs opened {val(new, 'R_SiriusB'):.5f} (hot WD: measured larger)")
check("Nauenberg relation on Sirius B: within 5% of the opened radius", abs(Rsb / val(new, "R_SiriusB") - 1) < 0.05)
R06 = R_nauenberg(val(new, "M_WD_typ"))
print(f"   Nauenberg R(0.6 M_sun) = {R06:.5f} R_sun = {R06 * val(ec, 'R_sun') / 1e3:.0f} km; 028 used {val(old, 'R_WD_typ')}")
check("028's typical WD radius 0.0125 R_sun agrees with the zero-T relation to < 3%", abs(R06 / val(old, "R_WD_typ") - 1) < 0.03)

# Baumbach-Allen validity: Lemaire 1112.3850 Eq. (5) says 1.05 < x < 3; 028 W2a uses it at x ~ 19
n0 = val(new, "ne_corona_model")


def ne(x):
    return n0 * (2.99 * x**-16 + 1.55 * x**-6 + 0.036 * x**-1.5)


print(f"   n_e(BA) at x = 3: {ne(3):.2e} cm^-3; at x = 19 (028's solar mixing radius): {ne(19):.2e} cm^-3 (extrapolated)")
check("flag: 028's solar plasma radius (~19 R_sun) is outside the Baumbach-Allen validity range x < 3", 19 > 3)

# Retino et al. kg -> eV
kg_eV = 5.60958860e35  # eV/c^2 per kg (CODATA 2018; exact from c and e)
print(f"   Retino 2016: 1.4e-49 kg = {1.4e-49 * kg_eV:.2e} eV, 3.4e-51 kg = {3.4e-51 * kg_eV:.2e} eV")
check("data row m_gamma_Retino = 1.4e-49 kg in eV", abs(val(new, "m_gamma_Retino") / (1.4e-49 * kg_eV) - 1) < 0.01)
summary()
