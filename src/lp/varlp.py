"""varlp.py -- Smyth-type auxiliary functions for the variance of totally real
algebraic integers.

For an interval I = [s0, s1] (the possible values of tau(x) after x -> +-x + k we may take
tau(x) in [0, 1/2]) find c_i >= 0 and m such that

    f(t) = dist(t, I)^2 - sum_i c_i log|Q_i(t)|  >=  m      for all real t.

Then every totally real algebraic integer x with tau(x) in I that is not a root of any Q_i
has Var(x) >= m: average f over the conjugates and use |N(Q_i(x))| >= 1.
"""
import numpy as np, sys, json
from scipy.optimize import linprog

def poly_logabs(coeffs, t):
    # coeffs: highest degree first (numpy convention)
    v = np.polyval(np.array(coeffs, dtype=float), t)
    with np.errstate(divide='ignore'):
        return np.log(np.abs(v))

def solve(I, polys, tgrid):
    s0, s1 = I
    d2 = np.where(tgrid < s0, (s0 - tgrid) ** 2, np.where(tgrid > s1, (tgrid - s1) ** 2, 0.0))
    L = np.array([poly_logabs(p, tgrid) for p in polys])          # k x T
    # drop grid points too close to a root (log -> -inf makes the constraint vacuous)
    ok = np.all(np.isfinite(L), axis=0) & np.all(L > -30, axis=0)
    L = L[:, ok]; d2 = d2[ok]
    k = len(polys)
    # variables: c_1..c_k, m ; maximise m  s.t.  m + sum c_i L_i(t) <= d2(t)
    A = np.hstack([L.T, np.ones((L.shape[1], 1))])
    res = linprog(c=np.r_[np.zeros(k), -1.0], A_ub=A, b_ub=d2,
                  bounds=[(0, None)] * k + [(None, None)], method='highs')
    if not res.success:
        return None, None
    return res.x[-1], res.x[:k]

if __name__ == '__main__':
    polys = json.load(open(sys.argv[1]))
    T = np.linspace(-7, 7.5, int(sys.argv[3]) if len(sys.argv) > 3 else 60001)
    width = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    worst = 9
    for s0 in np.arange(0, 0.5, width):
        m, c = solve((s0, min(s0 + width, 0.5)), polys, T)
        worst = min(worst, m)
        used = [(i, round(ci, 4)) for i, ci in enumerate(c) if ci > 1e-6]
        print(f"I=[{s0:.3f},{min(s0+width,0.5):.3f}]  m={m:.5f}  used={used}")
    print("min over intervals:", worst)
