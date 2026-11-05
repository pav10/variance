"""varlp_cert.py -- produce and rigorously verify auxiliary functions
    f_I(t) = dist(t, I)^2 - sum_i c_i log|Q_i(t)|  >=  m_I   (all real t)
for a cover of [0, 1/2] by intervals I.  The c_i are rounded down to rationals, and the
inequality is checked on all of R:
  * on a compact window by interval arithmetic (mpmath.iv), subdividing cells until the
    lower bound  dist_lo^2 - sum c_i log(max|Q_i(cell)|)  is >= m_I;
  * outside it by the explicit tail bound  (|t|-1/2)^2 - sum c_i (log A_i + deg_i log|t|),
    which is increasing there and is checked at the window ends.
Usage: python3 varlp_cert.py POLYS.json WIDTH NGRID MARGIN OUT.json
"""
import sys, json
from fractions import Fraction
import numpy as np
from mpmath import iv, mp, mpf, log as mlog
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from varlp import solve

iv.prec = 120
mp.prec = 120
LO, HI = -7.0, 7.5

def ivF(f):
    f = Fraction(f)
    return iv.mpf(f.numerator) / iv.mpf(f.denominator)

def iv_poly(coeffs, X):
    r = iv.mpf(0)
    for a in coeffs:
        r = r * X + a
    return r

def cell_lower(I, cs, polys, a, b):
    s0, s1 = I
    # exact lower bound of dist(t, I)^2 on [a, b]
    if b < s0:   dlo = (Fraction(s0) - b) ** 2
    elif a > s1: dlo = (a - Fraction(s1)) ** 2
    else:        dlo = Fraction(0)
    X = iv.mpf([ivF(a).a, ivF(b).b])
    tot = ivF(dlo)
    for c, p in zip(cs, polys):
        v = iv_poly(p, X)
        M = max(abs(v.a), abs(v.b))          # upper bound of |Q(t)| on the cell
        if M == 0: continue
        tot -= ivF(c) * iv.log(iv.mpf([M, M]))
    return tot.a                             # rigorous lower bound (as mpf interval endpoint)

def verify(I, cs, polys, m, cell=Fraction(1, 64), maxdepth=30):
    ncell = 0
    stack = []
    a = Fraction(LO)
    while a < HI:
        stack.append((a, min(a + cell, Fraction(HI)), 0)); a += cell
    worst = None
    while stack:
        a, b, dep = stack.pop()
        lb = cell_lower(I, cs, polys, a, b)
        ncell += 1
        if lb >= ivF(m).b:
            continue
        if dep >= maxdepth:
            return False, (float(a), float(b), float(lb)), ncell
        mid = (a + b) / 2
        stack.append((a, mid, dep + 1)); stack.append((mid, b, dep + 1))
    # tails: |t| >= 7: (|t| - 1/2)^2 - sum c (log A + deg log|t|) increasing if 2(|t|-1/2)|t| > sum c deg
    s = sum(c * (len(p) - 1) for c, p in zip(cs, polys))
    for T in (Fraction(-LO), Fraction(HI)):
        if 2 * (T - Fraction(1, 2)) * T <= s: return False, ('tail-mono', float(T)), ncell
        A = [sum(abs(x) for x in p) for p in polys]
        g = ivF((T - Fraction(1, 2)) ** 2) - sum((ivF(c) * (iv.log(iv.mpf(A_i)) + (len(p) - 1) * iv.log(ivF(T)))
                                             for c, p, A_i in zip(cs, polys, A)), iv.mpf(0))
        if g.a < ivF(m).b: return False, ('tail', float(T), float(g.a)), ncell
    return True, None, ncell

if __name__ == '__main__':
    polys = json.load(open(sys.argv[1]))
    width = Fraction(sys.argv[2]); ngrid = int(sys.argv[3]); margin = Fraction(sys.argv[4])
    out = sys.argv[5]
    T = np.linspace(LO, HI, ngrid)
    certs = []; overall = None
    s0 = Fraction(0)
    while s0 < Fraction(1, 2):
        s1 = min(s0 + width, Fraction(1, 2))
        mval, c = solve((float(s0), float(s1)), polys, T)
        used = [i for i in range(len(c)) if c[i] > 1e-7]
        cs = [Fraction(int(c[i] * 10**7), 10**7) for i in used]          # rounded down
        ps = [polys[i] for i in used]
        mcert = Fraction(int((mval - float(margin)) * 10**5), 10**5)
        ok, info, ncell = verify((s0, s1), cs, ps, mcert)
        print(f"I=[{float(s0):.4f},{float(s1):.4f}] LP={mval:.5f} certified m={float(mcert):.5f} ok={ok} "
              f"polys={len(ps)} cells={ncell} {info if not ok else ''}", flush=True)
        certs.append({'I': [str(s0), str(s1)], 'm': str(mcert), 'ok': ok,
                      'c': [str(x) for x in cs], 'Q': ps})
        if ok: overall = mcert if overall is None else min(overall, mcert)
        else: overall = Fraction(0)
        s0 = s1
    json.dump({'bound': str(overall), 'certs': certs}, open(out, 'w'))
    print('CERTIFIED BOUND:', float(overall), overall)
