"""KK-04: radion Brans-Dicke gamma - 1 versus Cassini, in orders of magnitude."""
import math

from common import check, load, summary, val

kk = load("kk_bounds_inputs.csv")
g, s = val(kk, "gamma_minus_1_Cassini"), val(kk, "gamma_minus_1_Cassini", "err")
b1 = s                      # |gamma - 1| < 2.3e-5, as the program quotes it (the 1 sigma error)
b2 = abs(g) + 2 * s         # conservative 2 sigma: |gamma - 1| < 6.7e-5
# theory: omega_BD = 0 radion, y-independent (bulk) matter: gamma = (1+w)/(2+w) = 1/2 (task 009);
# brane-localised matter alpha0^2 = 1/6: |gamma - 1| = 2 alpha0^2/(1+alpha0^2) = 2/7 (029 item 6, hand estimate)
pred = {"bulk matter (omega = 0)": 0.5, "brane matter (alpha0^2 = 1/6)": 2 / 7}
for name, p in pred.items():
    print(f"{name:30s}: |gamma-1| = {p:.4f}; orders short vs {b1:.1e}: {math.log10(p / b1):.2f}; "
          f"vs 2 sigma {b2:.1e}: {math.log10(p / b2):.2f}")
o = math.log10(0.5 / b1)
check("KK-04 bulk: 4.3 orders (029's ~4.3)", abs(o - 4.3) < 0.05, f"({o:.3f})")
check("KK-04 brane: still > 3.5 orders even at 2 sigma", math.log10((2 / 7) / b2) > 3.5)
summary()
