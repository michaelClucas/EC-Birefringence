"""EC-05: fermion-loop CFJ term from axial torsion, with the minimal-EC torsion source. Derivation Sec. B.

Coefficient conversions to our CFJ normalisation L = -(1/2) p_mu A_nu Ftilde^{mu nu} and our torsion
coupling L = (3/4) S_mu J5^mu (work/20261120-003):
  Dobado-Maroto astro-ph/9706044: coupling -(1/8) psibar gamma5 Sslash_DM psi with S_DM = eps T = 6 S (Eq. 4),
      L_eff = -(alpha/24 pi) sum Q^2 eps^{mu nu a b} F_{mu nu} A_a S_DM_b (Eq. 11), dimensional regularisation, m_f >> momenta
  Jackiw-Kostelecky hep-ph/9901358: coupling -b_mu psibar gamma^mu gamma5 psi (Eq. 1, e absorbed in A),
      L_k = (1/2) k_mu eps^{mu nu a b} F_{a b} A_nu, k = (3/16 pi^2) e^2 b (Eq. 16, "nonperturbative");
      Coleman-Glashow routing and Pauli-Villars give k = 0 (Eq. 17 and text); perturbatively undetermined.
Only magnitudes are compared (sign conventions for eps differ between the papers).
"""
import math
import sys

import sympy as sp

from common import Cosmo, beta_rad, load, report
from ec02_nogo import dn_closed, dneff
from scipy.optimize import brentq


def main():
    ok = True
    al, Q2, S = sp.symbols("alpha Q2 S", positive=True)
    e2 = 4 * sp.pi * al
    # DM: eps^{mu nu a b} F_{mu nu} = 2 Ftilde^{ab}  =>  L = -(alpha/12pi) Q2 Ftilde^{ab} A_a S_DM_b
    #     compare -(1/2) p_b A_a Ftilde^{ba} = +(1/2) p_b A_a Ftilde^{ab}  =>  |p| = (alpha/6pi) Q2 |S_DM|
    p_DM = al / (6 * sp.pi) * Q2 * 6 * S
    # JK: (1/2) k_mu eps^{mu nu a b} F_{ab} A_nu = k_mu Ftilde^{mu nu} A_nu  =>  |p| = 2|k|, k = 3 e^2 b/(16 pi^2), |b| = (3/4)|S|
    p_JK = 2 * 3 * e2 / (16 * sp.pi**2) * Q2 * sp.Rational(3, 4) * S
    # coupling map check: (1/8) S_DM = (3/4) S  <=>  S_DM = 6 S
    ok &= report("DM coupling (1/8) S_DM equals ours (3/4) S for S_DM = 6 S", sp.Rational(1, 8) * 6 == sp.Rational(3, 4))
    rDM = sp.simplify(p_DM / S)
    rJK = sp.simplify(p_JK / S)
    print(f"|p|/|S|:  Dobado-Maroto = {rDM},  Jackiw-Kostelecky = {rJK},  Coleman-Glashow / Pauli-Villars = 0")
    ok &= report("JK / DM = 9/8", sp.simplify(rJK / rDM - sp.Rational(9, 8)) == 0)
    ok &= report("JK value equals 003's |p| = 2 C e^2 (3/4)|S| with C = 3/16pi^2",
                 sp.simplify(rJK - 2 * sp.Rational(3, 16) / sp.pi**2 * e2 * sp.Rational(3, 4) * Q2) == 0)

    # gradient background: loop function for a derivative coupling to a massive fermion,
    # B1(tau) = 1 - tau f(tau)^2, f = arcsin(1/sqrt(tau)), tau = 4 m_f^2/q^2 (Bauer-Neubert-Thamm 1708.00443)
    t = sp.symbols("tau", positive=True)
    u = sp.symbols("u", positive=True)  # u = 1/tau
    B1 = 1 - sp.asin(sp.sqrt(u)) ** 2 / u
    ser = sp.series(B1, u, 0, 3).removeO()
    ok &= report("B1 = -1/(3 tau) - 8/(45 tau^2) + ...  (decoupling ~ q^2/(12 m_f^2))",
                 sp.simplify(ser - (-u / 3 - sp.Rational(8, 45) * u**2)) == 0)

    c = load()
    cos = Cosmo(c)
    beta = beta_rad(c)
    a_ = c["alpha_em"]
    rnum = {"JK (max)": float(rJK.subs({al: a_, Q2: c["SumQ2_SM"]})),
            "DM dim-reg": float(rDM.subs({al: a_, Q2: c["SumQ2_SM"]}))}
    xi95 = c["xi_nu"] + 2 * c["xi_nu_err"]
    dN95 = c["N_eff_obs"] + 2 * c["N_eff_obs_err"] - c["N_eff_SM"]
    xi_gen = brentq(lambda x: dneff(x) - dN95 / 2, 0, 5)
    Tnu = lambda a: cos.Tnu0 / a
    src = {"nu-A massless, xi95": lambda a: 3 * dn_closed(Tnu(a), xi95),
           "nu-C generous N_eff": lambda a: dn_closed(Tnu(a), xi95) + 2 * dn_closed(Tnu(a), xi_gen)}
    print(f"\nsum Q^2 = {c['SumQ2_SM']:.0f} (all SM charged fermions, no decoupling: generous)")
    for name, f in src.items():
        I = cos.t_int(f)
        for rn, r in rnum.items():
            bmax = r * cos.kappa / 4 * I
            print(f" {name:22s} {rn:11s}: |p|/|S| = {r:.4f}, max beta = {math.degrees(bmax):.1e} deg, "
                  f"shortfall {math.log10(beta / bmax):.1f} orders, Lambda_eff = M_P sqrt(2/r) = {cos.MP * math.sqrt(2 / r) / 1e9:.1e} GeV")
    # needed S for a constant background (003 number) through the loop
    t0 = cos.t_int(lambda a: 1.0)
    for rn, r in rnum.items():
        print(f" constant S needed through the loop ({rn}): |S^0hat| = {2 * beta / t0 / r:.1e} eV")
    # Reig-type bound on any universal J5.K operator (dual-axion completion + helioscope)
    xi_max = c["helioscope_gag_gae_GeV"] * 1e-9 / (cos.kappa * c["m_e_eV"])
    print(f" Reig/helioscope: universal J5.A Ftilde operator needs Lambda > {cos.MP * math.sqrt(2 / xi_max) / 1e9:.1e} GeV;"
          f" EC-05 has Lambda = {cos.MP * math.sqrt(2 / rnum['JK (max)']) / 1e9:.1e} GeV (allowed, but far too weak)")
    # physical suppression for a gradient (FRW) background with m_f >> H
    Hs = cos.H(cos.a_star)
    print(f" gradient S on FRW, CG/PV prescription or decoupling: extra factor ~ H_LSS^2/(12 m_e^2) = {Hs**2 / (12 * c['m_e_eV']**2):.1e}")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
