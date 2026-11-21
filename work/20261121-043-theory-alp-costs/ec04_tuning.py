"""EC-04: tuning of the PGT 0- torsion mode mass (task 007, Sec. C).

L = (A/2) S.S + (w/2)(div S)^2 + ...: the 0- mode has m^2 = A/w with A = 3/(2 kappa) = 1.5 M_P^2 (minimal EC).
w is the dimensionless coefficient of (div S)^2. m ~ H0 (needed to roll after recombination and give beta today)
requires w = 1.5 (M_P/H0)^2; an O(1) w gives m ~ M_P.
"""
import math
import sys

from common import H0, MP, check, done, info


def main():
    print("EC-04 PGT 0- mode: mass tuning (candidate mechanism, ALP class with hand-added CFJ coupling)")
    A = 1.5 * MP**2
    w = {mu: A / (mu * H0) ** 2 for mu in (0.3, 1.0, 5.0)}
    for mu, v in w.items():
        info(f"m = {mu} H0: w = A/m^2 = {v:.2e}")
    info(f"(M_P/H0)^2 = {(MP / H0) ** 2:.2e}; O(1) w gives m = sqrt(1.5) M_P = {math.sqrt(1.5) * MP / 1e9:.2e} GeV")
    check("T1 w(m = H0) = 1.5 (M_P/H0)^2 = 4.3e120 (task 007 quotes ~4e120)", abs(w[1.0] / 4.3e120 - 1) < 0.05)
    check("T2 thawing range m = 0.3-5 H0 spans w = 1.7e119 - 4.8e121", 1e119 < w[5.0] < 2e119 and 4e121 < w[0.3] < 5e121)
    return w


if __name__ == "__main__":
    main()
    sys.exit(done())
