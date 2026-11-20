"""Derived columns of data/kaf_lowz_bounds.csv (task 024).

For each row, from the published value (columns central, err_1sigma_or_limit, kind) compute:
  kAF_T_limit_eV  : bound on |(k_AF)_T| for a constant isotropic CFJ coefficient, in eV
  p0_limit_eV     : same as our CFJ p_0 = 2 |(k_AF)_T| (CONVENTIONS.md, literature/bounds.md Sec. 2)
  ec06_xi_eps_limit: bound on xi*eps for an EC-06 background p_0hat = xi*eps*H, for which
                     beta(z) = (xi*eps/2) ln(1+z)   (beta = 1/2 Int p_0hat dt, Int H dt = ln(1+z))
Rules:
  measurements (meas_*): limit = |central| + 2 sigma ("~95%", Gaussian 1 sigma assumed)
  angle -> constant k_AF: |(k_AF)_T| = alpha / t_L(z_eff)   (KM 0905.0031 Eq. 158 with constant coefficient)
  h0-scaled values use h from data/cosmo_inputs.csv (Planck 2018).
Also prints the CMB-side requirement xi*eps from data/birefringence_measurements.csv (ACT DR6 beta) for comparison.

Usage (repo root): python work/20261120-024-lit-ec06-novelty/kaf_lowz.py [--write]
Without --write it recomputes and checks the columns already in the file (exit 1 on mismatch > 1%).
"""
import csv, math, sys

ROOT = '.'
def load(path, key='id'):
    with open(f'{ROOT}/{path}', newline='', encoding='utf-8') as f:
        return {r[key]: r for r in csv.DictReader(f)}

cosmo = load('data/cosmo_inputs.csv')
units = load('data/bounds_inputs.csv')
meas = load('data/birefringence_measurements.csv')

H0_kms = float(cosmo['H0']['value'])
Om = float(cosmo['Omega_m']['value'])
z_star = float(cosmo['z_star']['value'])
h = H0_kms / 100.0
hbar = float(units['hbar']['value'])              # eV s
Mpc_km = float(units['Mpc']['value']) / 1e3
GeV = float(units['GeV']['value'])
H0_eV = H0_kms / Mpc_km * hbar                     # eV

def E(z):
    return math.sqrt(Om * (1 + z) ** 3 + (1 - Om))

def lookback_eV_inv(z, n=20000):
    # t_L = Int_0^z dz / ((1+z) H(z)), flat LCDM, radiation neglected (error < 1e-3 for z < 4)
    if z <= 0:
        return 0.0
    s, dz = 0.0, z / n
    for i in range(n):
        zz = (i + 0.5) * dz
        s += dz / ((1 + zz) * E(zz))
    return s / H0_eV

def derive(r):
    kind = r['kind']
    c = float(r['central'] or 0)
    e = float(r['err_1sigma_or_limit'] or 0)
    zeff = float(r['z_eff']) if r['z_eff'] else None
    k = p = x = None
    note = ''
    if kind == 'meas_angle_deg':
        a = math.radians(abs(c) + 2 * e)
        if zeff:
            tL = lookback_eV_inv(zeff)
            k = a / tL
            p = 2 * k
            x = 2 * a / math.log(1 + zeff)
            note = (f'alpha_lim = |{c}|+2x{e} deg = {math.degrees(a):.2f} deg; t_L(z={zeff}) = {tL:.3e} eV^-1 '
                    f'(H0={H0_kms}, Om={Om}); k = alpha/t_L, p0 = 2k; xi*eps = 2 alpha/ln(1+z)')
        else:
            note = f'alpha_lim = {math.degrees(a):.2f} deg; no z_eff, no conversion'
    elif kind in ('upper_p0_GeVh', 'meas_p0_GeVh', 'meas_p0_eVh'):
        scale = (GeV if 'GeV' in kind else 1.0) * h
        v = e if kind.startswith('upper') else abs(c) + 2 * e
        p = v * scale
        k = p / 2
        note = (f'p0_lim = {v:.3g} x {"1e9 x " if "GeV" in kind else ""}h ({h:.4f}) eV'
                f'{"" if kind.startswith("upper") else " (|central|+2 sigma)"}; k_AF = p0/2')
    elif kind == 'upper_kAF_GeV':
        k = e * GeV
        p = 2 * k
        note = 'k_AF x 1e9 eV/GeV; p0 = 2 k_AF'
    elif kind == 'upper_spatial_GeV':
        note = f'spatial components only: < {e * GeV:.1e} eV; no bound on (k_AF)_T'
    return k, p, x, note

def fmt(v):
    return '' if v is None else f'{v:.3g}'

path = f'{ROOT}/data/kaf_lowz_bounds.csv'
with open(path, newline='', encoding='utf-8') as f:
    rd = csv.DictReader(f)
    fields, rows = rd.fieldnames, list(rd)

write = '--write' in sys.argv
bad = 0
print(f'H0 = {H0_eV:.4e} eV, h = {h}, t_L(z*) = {lookback_eV_inv(z_star):.4e} eV^-1')
for r in rows:
    k, p, x, note = derive(r)
    new = {'kAF_T_limit_eV': fmt(k), 'p0_limit_eV': fmt(p), 'ec06_xi_eps_limit': fmt(x), 'conversion_note': note}
    if write:
        r.update(new)
    else:
        for col in ('kAF_T_limit_eV', 'p0_limit_eV', 'ec06_xi_eps_limit'):
            old, nv = r[col], new[col]
            ok = (old == nv == '') or (old and nv and abs(float(old) / float(nv) - 1) < 0.01)
            if not ok:
                bad += 1
                print('MISMATCH', r['id'], col, old, nv)
    print(f"{r['id']:18s} k_AF<{new['kAF_T_limit_eV']:>9s} eV  p0<{new['p0_limit_eV']:>9s} eV  xi*eps<{new['ec06_xi_eps_limit']:>7s}")

if write:
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

# CMB-side requirement for EC-06, for comparison
b = math.radians(float(meas['DPK25-beta-baseline']['value']))
xe = 2 * b / math.log(1 + z_star)
print(f'EC-06 needs xi*eps = {xe:.3g} for ACT beta = {math.degrees(b)} deg (z* = {z_star})')
for zz in (0.47, 0.811, 2.80):
    print(f'  that xi*eps gives beta(z={zz}) = {math.degrees(xe / 2 * math.log(1 + zz)):.4f} deg')
# constant k_AF: CMB-implied value vs low-z
kc = b / lookback_eV_inv(z_star)
print(f'constant k_AF fitting ACT beta: |(k_AF)_T| = {kc:.3g} eV (p0 = {2 * kc:.3g} eV)')
print('PASS' if bad == 0 else f'FAIL ({bad} mismatches)')
sys.exit(1 if bad else 0)
