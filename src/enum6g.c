/* [enum6g.c: enum6.c of pav10/2O with one change: the quadratic filter F3 is skipped when the
 *  environment variable ENUM_NOF3 is set (used for G in regime nu_2 < 2 nu, where F3 is not valid)]
 *
 * enum6.c -- enumerate the possible shortest vectors x_1 of the variance
 * lattice of a totally real sextic field K admitting a universal Z-form,
 * in the case K = Q(x_1).
 *
 * Output: every monic f = x^6 + a1 x^5 + ... + a6 in Z[x] such that
 *   (E1) f has only real roots,
 *   (E2) a1 in {0,-1,-2,-3}              (normalisation x -> +-x + k),
 *   (E3) S := sum (theta_i - tau)^2 <= SMAX = 15 + 3 sqrt(21)
 *        (i.e. Var = S/6 <= (5+sqrt21)/2, Lemma "nu bound"),
 * and whose roots pass the filters (with V = S/6, tau = -a1/6):
 *   (F1) downspread  tau - (ceil(theta_min) - 1) >= V,
 *   (F2) upspread    (floor(theta_max) + 1) - tau >= V,
 *   (F3) no monic P in Z[t] of degree 2 with P(theta_i) > 0 for all i and
 *        P(tau) < 0.
 * Each filter expresses tau(alpha) >= nu(K) = V for a specific totally
 * positive non-rational alpha in Z[x_1] (x_1 - n, -x_1 - n', P(x_1)).
 * All comparisons are made with safety margins in the direction that keeps
 * polynomials; exact re-verification of everything is done in PARI/GP.
 *
 * Enumeration: Rolle/interlacing ("Robinson") method.  With
 *   g_m(x) = sum_{i=0}^m a_i C(6-i, m-i) x^{m-i}   (so g_6 = f, g_m' = (7-m) g_{m-1}),
 * f real-rooted implies every g_m real-rooted with roots interlacing those of
 * g_{m-1}; given a_0..a_{m-1} this confines a_m (the constant term of g_m)
 * to an explicit interval.  The enumeration is exhaustive.
 *
 * usage: enum6 a1 a2 [smax [all|strong|why]]   (the subtree with the given a1,a2)
 * output lines: a1 a2 a3 a4 a5 a6
 * stderr: statistics.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef long double ld;
#ifndef DEG
#define DEG 6
#endif
#define N DEG
#define EPS 1e-9L

static long a[N + 1];
static ld Cb[N + 1][N + 1];
static ld Rt[N + 1][N + 1];      /* Rt[m][0..m-1]: roots of g_m, ascending */
static ld LB, UB, SMAX, TAU, VAR;
/* ---- split modes (sharpened rigidity), see notes: "Case A pruning" ------------------
 * m := min(3/2 V, V + 3/4).  With n = ceil(theta_1) - 1 and tau_- = tau - n:
 *  if tau_- < m then x_1 - n = y^2 with y in O_K (square dichotomy), tau(y)^2 <= tau_- - V,
 *     sum y^2 = d tau_-, sum y^4 = S + d tau_-^2;           -> mode YL (enumerate y)
 *  symmetrically n' - x_1 = y^2 if n' - tau < m (n' = floor(theta_d) + 1)  -> mode YU
 *  otherwise theta_1 <= floor(tau - m) + 1 and theta_d >= ceil(tau + m) - 1 -> mode X
 * The three modes together cover every admissible x_1.
 * NU is the budget threshold nu(K) = Var(x_1) (= VAR except in the y-enumeration). */
static int YMODE = 0, XSPLIT = 0;
static long YN = 0;                 /* rung: x_1 = YN + y^2 (YL) or YN - y^2 (YU) */
static long long YP4 = 0;           /* pinned power sum p_4(y) */
static ld TAUX = 0, VARX = 0, NU = 0, XA = 0, XB = 0;
static long long n_rung = 0, n_ybud = 0;
static long long n_leaf = 0, n_out = 0, n_f1 = 0, n_f2 = 0, n_f3 = 0, n_ns = 0, n_rig = 0;
static int NOFILTER = 0, STRONG = 0, WHY = 0, SHARP = 0, NOF3 = 0;   /* NOF3: env ENUM_NOF3 (G, regime 2) */
static long A3K = 0, A3N = 0;          /* optional split of a subtree by a3 mod A3N (env ENUM_A3MOD=k/n) */
/* moment pruning (only when filters are on): F1/F2 force theta_1 - tau <= -AA and
 * theta_d - tau >= BB, hence every even centred power sum M_k >= AA^k + BB^k.
 * M_k depends only on a_1..a_k, so whole subtrees are cut at level k = 4, 6. */
static ld AA = 0, BB = 0, MIDR = -1;   /* MIDR: radius of the window for theta_2..theta_{d-1} */
static long long n_prune = 0;

static ld g(int m, ld x)
{
    ld v = 0;
    for (int i = 0; i <= m; i++) v = v * x + (ld)a[i] * Cb[N - i][m - i];
    return v;
}

static ld findroot(int m, ld lo, ld hi)
{
    ld flo = g(m, lo), fhi = g(m, hi);
    if (flo == 0) return lo;
    if (fhi == 0) return hi;
    if ((flo > 0) == (fhi > 0)) return fabsl(flo) < fabsl(fhi) ? lo : hi; /* tangency */
    ld x = 0.5L * (lo + hi);
    for (int it = 0; it < 300; it++) {
        ld fx = g(m, x);
        if (fx == 0) return x;
        if ((fx > 0) == (flo > 0)) { lo = x; flo = fx; } else { hi = x; fhi = fx; }
        ld d = (ld)(N - m + 1) * g(m - 1, x);
        ld xn = (d != 0) ? x - fx / d : 0.5L * (lo + hi);
        if (!(xn > lo && xn < hi)) xn = 0.5L * (lo + hi);
        if (fabsl(xn - x) <= 1e-18L * (1 + fabsl(x))) return xn;
        x = xn;
        if (hi - lo <= 1e-18L * (1 + fabsl(x))) return x;
    }
    return x;
}

/* F3: is there a monic integer quadratic with both roots in the gap
 * (th[k], th[k+1]) that contains tau, and with tau strictly between them? */
static int quad_violation(const ld *th)
{
    int k = -1;
    for (int i = 0; i + 1 < N; i++)
        if (th[i] < TAU && TAU < th[i + 1]) { k = i; break; }
    if (k < 0) return 0;               /* tau (numerically) at a root: keep */
    ld l = th[k], r = th[k + 1];
    long bmin = (long)ceill(-2 * r), bmax = (long)floorl(-2 * l);
    for (long b = bmin; b <= bmax; b++) {
        ld lo1 = -l * l - b * l, lo2 = -r * r - b * r;
        ld lo = lo1 > lo2 ? lo1 : lo2;
        ld up = -TAU * TAU - b * TAU;
        /* need integer c with lo < c < up; declare violation only with margin */
        long c = (long)floorl(lo + 1e-7L) + 1;
        if ((ld)c > lo + 1e-7L && (ld)c < up - 1e-7L) return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------------
 * STRONG mode (optional): order-level tests inside Z[x_1], done in C.
 *  NS : some y = c x^2 + b x  or  y = x^3 + a x^2 + b x  (small integers) has
 *       Var(y) < Var(x_1): then x_1 is not a shortest vector of Lambda_K.
 *  RIG: rigidity for the downspread elements alpha_- = x_1 - n, alpha_+ = -x_1 + n'
 *       (they generate K when x_1 does): if tau(alpha) < 3/2 Var(x_1) there must be
 *       kappa = 0,3 mod 4 and z in O_K with z^2 = 4 alpha - kappa and
 *       4 Var(x_1) <= Var(z) <= 4 tau(alpha);  here z = sum of the conjugates with
 *       signs, so some sign pattern of +-sqrt(4 sigma_j(alpha) - kappa) must have
 *       integral elementary symmetric functions.  Range of kappa: see notes, Lemma 5.1
 *       (0 < 4 min alpha - kappa <= max(d^2 tau/(d-1), Var(alpha)/Var(x_1))).
 * All tests exclude only on a clear numerical violation.
 * ---------------------------------------------------------------------- */
static ld varv(const ld *v)
{
    ld m = 0, s2 = 0;
    for (int j = 0; j < N; j++) m += v[j];
    m /= N;
    for (int j = 0; j < N; j++) s2 += (v[j] - m) * (v[j] - m);
    return s2 / N;
}

static int ns_violation(const ld *th)
{
    ld p1[N], p2[N], p3[N], y[N];
    ld tol = 1e-9L * (1 + VAR);
    for (int j = 0; j < N; j++) { p1[j] = th[j]; p2[j] = th[j] * th[j]; p3[j] = p2[j] * th[j]; }
    /* covariances */
    ld m1 = 0, m2 = 0, m3 = 0;
    for (int j = 0; j < N; j++) { m1 += p1[j]; m2 += p2[j]; m3 += p3[j]; }
    m1 /= N; m2 /= N; m3 /= N;
    ld V11 = 0, V22 = 0, V33 = 0, C21 = 0, C31 = 0, C32 = 0;
    for (int j = 0; j < N; j++) {
        ld d1 = p1[j] - m1, d2 = p2[j] - m2, d3 = p3[j] - m3;
        V11 += d1 * d1; V22 += d2 * d2; V33 += d3 * d3;
        C21 += d2 * d1; C31 += d3 * d1; C32 += d3 * d2;
    }
    V11 /= N; V22 /= N; V33 /= N; C21 /= N; C31 /= N; C32 /= N;
    /* y = c x^2 + b x, c = 1..3 */
    for (int c = 1; c <= 3; c++) {
        ld bs = -c * C21 / V11;
        for (long b = (long)floorl(bs) - 1; b <= (long)ceill(bs) + 1; b++) {
            ld v = c * c * V22 + 2 * c * b * C21 + (ld)b * b * V11;
            if (v < VAR - tol) {
                /* confirm directly from the values; y is non-rational since deg x_1 >= 3 */
                for (int j = 0; j < N; j++) y[j] = c * p2[j] + b * p1[j];
                ld vy = varv(y);
                if (vy < VAR - tol && vy > tol) return 1;
            }
        }
    }
    /* y = x^3 + a x^2 + b x (non-rational only when deg x_1 >= 4):
     * real minimiser of the quadratic in (a,b), then a box */
    ld det = V22 * V11 - C21 * C21;
    if (N >= 4 && det > 0) {
        ld as = (-C32 * V11 + C31 * C21) / det, bs = (-C31 * V22 + C32 * C21) / det;
        for (long aa = (long)floorl(as) - 2; aa <= (long)ceill(as) + 2; aa++)
            for (long bb = (long)floorl(bs) - 2; bb <= (long)ceill(bs) + 2; bb++) {
                ld v = V33 + 2 * aa * C32 + 2 * bb * C31 + (ld)aa * aa * V22 + 2 * (ld)aa * bb * C21 + (ld)bb * bb * V11;
                if (v < VAR - tol) {
                    for (int j = 0; j < N; j++) y[j] = p3[j] + aa * p2[j] + bb * p1[j];
                    ld vy = varv(y);
                    if (vy < VAR - tol && vy > tol) return 1;
                }
            }
    }
    return 0;
}

/* av[0..N-1]: conjugates of a totally positive alpha generating K, tau(alpha) < 3/2 VAR.
 * Returns 1 if the rigidity-by-squares condition can hold (keep), 0 if it certainly fails. */
static int sq_ok(const ld *av)
{
    ld mn = av[0], ta = 0;
    for (int j = 0; j < N; j++) { if (av[j] < mn) mn = av[j]; ta += av[j]; }
    ta /= N;
    ld va = varv(av);
    ld u1 = (ld)N * N * ta / (N - 1), u2 = va / VAR;
    ld umax = (u1 > u2 ? u1 : u2) + 1;
    long kmin = (long)floorl(4 * mn - umax) - 1, kmax = (long)ceill(4 * mn);
    if (SHARP) {
        /* sharpened rigidity: kappa = 4b - c^2 = Q(e)Q(r) - B(e,r)^2 >= 0 (Cauchy-Schwarz),
         * and kappa + tau(z)^2 <= 4 (tau(alpha) - nu)  (from Var(z) >= 4 nu) */
        kmin = 0;
        long k2 = (long)floorl(4 * (ta - VAR) + 1e-9L);
        if (k2 < kmax) kmax = k2;
    }
    ld s[N], c[N + 1];
    for (long kap = kmin; kap <= kmax; kap++) {
        if ((ld)kap >= 4 * mn - 1e-12L || kap > kmax) break;
        long r = ((kap % 4) + 4) % 4;
        if (r == 1 || r == 2) continue;
        ld rhs = 4 * ta - kap - 4 * VAR;         /* (e1/N)^2 <= rhs  <=>  Var(z) >= 4 VAR */
        if (rhs < -1e-9L) continue;
        for (int j = 0; j < N; j++) s[j] = sqrtl(4 * av[j] - kap);
        for (long pat = 0; pat < (1L << (N - 1)); pat++) {
            ld e1 = s[0];
            for (int j = 1; j < N; j++) e1 += ((pat >> (j - 1)) & 1) ? -s[j] : s[j];
            if (fabsl(e1 - roundl(e1)) > 1e-6L) continue;
            ld q = (e1 / N) * (e1 / N);
            if (q > rhs + 1e-9L) continue;          /* Var(z) < 4 VAR */
            if (q < -(ld)kap - 1e-9L) continue;     /* Var(z) > 4 tau(alpha) */
            /* all elementary symmetric functions integral? */
            for (int k = 0; k <= N; k++) c[k] = 0;
            c[0] = 1;
            for (int j = 0; j < N; j++) {
                ld zj = (j == 0 || !((pat >> (j - 1)) & 1)) ? s[j] : -s[j];
                for (int k = j + 1; k >= 1; k--) c[k] -= zj * c[k - 1];
            }
            int ok = 1;
            for (int k = 1; k <= N; k++)
                if (fabsl(c[k] - roundl(c[k])) > 1e-5L * (1 + fabsl(c[k]) * 1e-12L)) { ok = 0; break; }
            if (ok) return 1;
        }
    }
    return 0;
}

/* returns 0 (pass), 1 (alpha_- = x - nm fails), 2 (alpha_+ = np - x fails) */
static int rig_violation(const ld *th, long *nout)
{
    ld av[N];
    long nm = (long)ceill(th[0]) - 1, np = (long)floorl(th[N - 1]) + 1;
    if (TAU - nm < 1.5L * VAR - 1e-9L) {
        for (int j = 0; j < N; j++) av[j] = th[j] - nm;
        if (!sq_ok(av)) { *nout = nm; return 1; }
    }
    if (np - TAU < 1.5L * VAR - 1e-9L) {
        for (int j = 0; j < N; j++) av[j] = np - th[j];
        if (!sq_ok(av)) { *nout = np; return 2; }
    }
    return 0;
}

static int USEX = 0;              /* y-mode: printpoly prints the x-polynomial XCO[] */
static long long XCO[N + 1];

static void printpoly(const char *tag, long extra)
{
    if (tag) printf("%s %ld ", tag, extra);
    if (USEX) { for (int i = 1; i <= N; i++) printf("%lld%c", XCO[i], i == N ? '\n' : ' '); return; }
    for (int i = 1; i <= N; i++) printf("%ld%c", a[i], i == N ? '\n' : ' ');
}

static int xfilters(const ld *th)      /* x-context filters; returns 1 if th passes */
{
    ld dm = TAU - (ceill(th[0]) - 1);
    if (dm < VAR - 1e-7L) { n_f1++; return 0; }
    ld dp = (floorl(th[N - 1]) + 1) - TAU;
    if (dp < VAR - 1e-7L) { n_f2++; return 0; }
    if (N >= 3 && quad_violation(th)) { n_f3++; return 0; }
    if (STRONG) {
        long nn = 0;
        int rv;
        if (N >= 3 && ns_violation(th)) { n_ns++; if (WHY) printpoly("NS", 0); return 0; }
        if ((rv = rig_violation(th, &nn))) { n_rig++; if (WHY) printpoly(rv == 1 ? "RIGM" : "RIGP", nn); return 0; }
    }
    return 1;
}

/* leaf of the y-enumeration: th = conjugates of y (ascending) */
static void yleaf(ld *th)
{
    /* y itself is a non-rational element of O_K: its downspreads are >= nu(K) = NU */
    if (TAU - (ceill(th[0]) - 1) < NU - 1e-7L || (floorl(th[N - 1]) + 1) - TAU < NU - 1e-7L) {
        n_ybud++;
        if (WHY) printpoly(YMODE == 1 ? "YBL" : "YBU", YN);     /* y-polynomial, rung n */
        return;
    }
    ld xs[N], mn = 1e30L;
    for (int j = 0; j < N; j++) {
        ld q = th[j] * th[j];
        if (q < mn) mn = q;
        xs[j] = (YMODE == 1) ? YN + q : YN - q;
    }
    if (mn > 1 + 1e-9L) { n_rung++; return; }          /* wrong rung (another rung covers it) */
    for (int i = 1; i < N; i++)                          /* sort ascending */
        for (int j = i; j > 0 && xs[j] < xs[j - 1]; j--) { ld t = xs[j]; xs[j] = xs[j - 1]; xs[j - 1] = t; }
    /* exact coefficients of the x-polynomial from the integer y-polynomial g = a[]:
     * g(t) = E(t^2) + t O(t^2), prod (u - y_j^2) = (-1)^N (E(u)^2 - u O(u)^2) =: G(u);
     * lower: X(t) = G(t - n);  upper: X(t) = (-1)^N G(n - t). */
    __int128 E[N + 1], O[N + 1], G[N + 1], X[N + 1];
    for (int k = 0; k <= N; k++) E[k] = O[k] = G[k] = X[k] = 0;      /* index = power of u */
    for (int i = 0; i <= N; i++) {
        int pw = N - i;                                             /* power of t */
        if (pw % 2 == 0) E[pw / 2] += a[i]; else O[(pw - 1) / 2] += a[i];
    }
    for (int p = 0; p <= N; p++) for (int q = 0; q <= N; q++) {
        if (p + q <= N) G[p + q] += E[p] * E[q];
        if (p + q + 1 <= N) G[p + q + 1] -= O[p] * O[q];
    }
    if (N % 2) for (int k = 0; k <= N; k++) G[k] = -G[k];
    /* substitute u = s*t + c  with (s,c) = (1,-n) lower, (-1,n) upper; Horner in t */
    long sg = (YMODE == 1) ? 1 : -1, cc = (YMODE == 1) ? -YN : YN;
    for (int k = N; k >= 0; k--) {
        /* X := X * (sg t + cc) + G[k] */
        __int128 T[N + 1];
        for (int j = 0; j <= N; j++) T[j] = 0;
        for (int j = 0; j < N; j++) { T[j + 1] += X[j] * sg; T[j] += X[j] * cc; }
        for (int j = 0; j <= N; j++) X[j] = T[j];
        X[0] += G[k];
    }
    if (YMODE == 2 && N % 2) for (int k = 0; k <= N; k++) X[k] = -X[k];
    if (X[N] != 1) { fprintf(stderr, "WARNING: x-polynomial not monic\n"); return; }
    /* cross-check against the numerical roots */
    {
        ld c[N + 1];
        for (int k = 0; k <= N; k++) c[k] = 0;
        c[0] = 1;
        for (int j = 0; j < N; j++)
            for (int k = j + 1; k >= 1; k--) c[k] -= xs[j] * c[k - 1];
        for (int k = 1; k <= N; k++)
            if (fabsl(c[k] - (ld)X[N - k]) > 1e-4L * (1 + fabsl(c[k])))
                fprintf(stderr, "NOTE: numerical/exact mismatch at coefficient %d (%Lg vs %Lg)\n", k, c[k], (ld)X[N - k]);
    }
    for (int k = 1; k <= N; k++) XCO[k] = (long long)X[N - k];
    ld st = TAU, sv = VAR;
    TAU = TAUX; VAR = VARX; USEX = 1;
    int ok = xfilters(xs);
    TAU = st; VAR = sv;
    if (ok) { n_out++; printpoly(WHY ? "OUT" : NULL, 0); }
    USEX = 0;
}

static void leaf(void)
{
    ld br[N + 2], th[N];
    n_leaf++;
    br[0] = LB;
    for (int i = 0; i < N - 1; i++) br[i + 1] = Rt[N - 1][i];
    br[N] = UB;
    for (int i = 0; i < N; i++) th[i] = findroot(N, br[i], br[i + 1]);
    if (YMODE) { yleaf(th); return; }
    if (NOFILTER) { n_out++; for (int i = 1; i <= N; i++) printf("%ld%c", a[i], i == N ? '\n' : ' '); return; }
    /* F1: downspread.  theta_min is never an integer for irreducible f. */
    ld dm = TAU - (ceill(th[0]) - 1);
    if (dm < NU - 1e-7L) { n_f1++; return; }
    ld dp = (floorl(th[N - 1]) + 1) - TAU;
    if (dp < NU - 1e-7L) { n_f2++; return; }
    if (XSPLIT && (th[0] > XA + 1e-9L || th[N - 1] < XB - 1e-9L)) { n_rung++; return; }
    if (!NOF3 && N >= 3 && quad_violation(th)) { n_f3++; return; }   /* P(x_1) non-rational needs deg x_1 >= 3 */
    if (STRONG) {
        long nn = 0; int rv;
        if (N >= 3 && ns_violation(th)) { n_ns++; if (WHY) printpoly("NS", 0); return; }  /* c x^2 + b x non-rational needs deg >= 3 */
        if ((rv = rig_violation(th, &nn))) { n_rig++; if (WHY) printpoly(rv == 1 ? "RIGM" : "RIGP", nn); return; }
    }
    n_out++;
    printpoly(WHY ? "OUT" : NULL, 0);
}

/* centred power sum M_k = sum_j (theta_j - tau)^k from a_1..a_k (Newton's identities) */
static ld centred_moment(int k)
{
    ld e[N + 1], p[N + 1];
    e[0] = 1;
    for (int i = 1; i <= k; i++) e[i] = (i % 2 ? -1 : 1) * (ld)a[i];
    p[0] = N;
    for (int j = 1; j <= k; j++) {
        ld v = (j % 2 ? 1 : -1) * j * e[j];
        for (int i = 1; i < j; i++) v += (i % 2 ? 1 : -1) * e[i] * p[j - i];
        p[j] = v;
    }
    ld M = 0, bin = 1, t = -TAU;
    for (int i = 0; i <= k; i++) {
        /* C(k,i) p_i t^(k-i) */
        ld tp = 1;
        for (int r = 0; r < k - i; r++) tp *= t;
        M += bin * p[i] * tp;
        bin = bin * (k - i) / (i + 1);
    }
    return M;
}

static void rec(int m)
{
    ld lo = -1e30L, hi = 1e30L;
    a[m] = 0;
    for (int i = 0; i < m - 1; i++) {
        ld hv = g(m, Rt[m - 1][i]);
        int k = m - (i + 1);
        if (k % 2 == 0) { if (-hv > lo) lo = -hv; }
        else            { if (-hv < hi) hi = -hv; }
    }
    long amin = (long)ceill(lo - EPS * (1 + fabsl(lo)));
    long amax = (long)floorl(hi + EPS * (1 + fabsl(hi)));
    if (YMODE && m == 4) {
        /* Newton: p4 = e1 p3 - e2 p2 + e3 p1 - 4 e4 with e_i = (-1)^i a_i; p4 is pinned */
        long long e1 = -a[1], e2 = a[2], e3 = -a[3];
        long long p1 = e1, p2 = e1 * p1 - 2 * e2, p3 = e1 * p2 - e2 * p1 + 3 * e3;
        long long num = e1 * p3 - e2 * p2 + e3 * p1 - YP4;
        if (num % 4 != 0) return;
        long a4 = (long)(num / 4);
        if (a4 < amin || a4 > amax) return;
        amin = amax = a4;
    }
    for (long am = amin; am <= amax; am++) {
        if (m == 3 && A3N > 0 && ((am % A3N) + A3N) % A3N != A3K) continue;   /* ENUM_A3MOD=k/n */
        a[m] = am;
        if (m == N) { leaf(); continue; }
        if (!NOFILTER && ((N > 4 && m == 4) || (N > 6 && m == 6)) && AA > 0 && BB > 0) {
            ld Mk = centred_moment(m), bound = powl(AA, m) + powl(BB, m);
            if (Mk < bound * (1 - 1e-12L) - 1e-9L) { n_prune++; continue; }
        }
        ld br[N + 2];
        br[0] = LB;
        for (int i = 0; i < m - 1; i++) br[i + 1] = Rt[m - 1][i];
        br[m] = UB;
        for (int i = 0; i < m; i++) Rt[m][i] = findroot(m, br[i], br[i + 1]);
        /* middle-window pruning: the i-th root of f^(d-m) lies in [theta_i, theta_{i+d-m}];
         * for 2 <= i <= m-1 this is inside [theta_2, theta_{d-1}], and F1/F2 force all middle
         * roots into tau +- MIDR with MIDR^2 = S - AA^2 - BB^2. */
        if (N > 3 && MIDR >= 0 && m >= 3) {
            int bad = 0;
            for (int i = 1; i <= m - 2; i++)
                if (Rt[m][i] < TAU - MIDR - 1e-9L || Rt[m][i] > TAU + MIDR + 1e-9L) { bad = 1; break; }
            if (bad) { n_prune++; continue; }
        }
        rec(m + 1);
    }
    a[m] = 0;
}

/* run the y-enumeration for all forced rungs; mode 1 = lower (x = n + y^2), 2 = upper (x = n - y^2) */
static int yrun(int mode, ld SX, ld mm)
{
    long a1x = a[1], a2x = a[2];
    YMODE = mode; TAUX = TAU; VARX = VAR; NU = VAR;
    long long SXN = (long long)(N - 1) * a1x * a1x - 2LL * N * a2x;      /* N * S_x, exact */
    /* rungs: T = tau - n (lower) or n - tau (upper) in [V, mm) */
    long nlo = (long)floorl(-4 * (SX + 10)) , nhi = -nlo;
    for (long n = nlo; n <= nhi; n++) {
        ld T = (mode == 1) ? TAUX - n : n - TAUX;
        if (T < VARX - 1e-12L || T >= mm + 1e-12L) continue;
        long long p2 = (mode == 1) ? (long long)(-a1x) - (long long)N * n : (long long)N * n + a1x;   /* = N T */
        long long p4N = SXN + p2 * p2;                  /* N * p4 */
        if (p4N % N != 0) continue;
        YP4 = p4N / N; YN = n;
        long smax = (long)floorl(N * sqrtl(fmaxl(T - VARX, 0)) + 1e-9L);
        for (long sy = 0; sy <= smax; sy++) {
            if ((sy * sy - p2) % 2 != 0) continue;
            a[1] = -sy; a[2] = (sy * sy - p2) / 2;
            ld Sy = (ld)p2 - (ld)sy * sy / N;              /* centred T2 of y */
            if (Sy <= 0) continue;
            TAU = (ld)sy / N; VAR = Sy / N;
            AA = TAU - floorl(TAU - NU + 1 + 1e-12L);
            BB = ceill(TAU + NU - 1 - 1e-12L) - TAU;
            MIDR = -1;
            if (NU <= 1) AA = BB = 0;
            if (AA > 0 && BB > 0) {
                ld r2 = Sy - AA * AA - BB * BB;
                if (r2 < -1e-9L) continue;
                MIDR = sqrtl(r2 > 0 ? r2 : 0) * (1 + 1e-12L) + 1e-12L;
            }
            ld rad = sqrtl((N - 1) * Sy / N);
            LB = TAU - rad - 0.01L; UB = TAU + rad + 0.01L;
            Rt[1][0] = TAU;
            ld br[3] = {LB, TAU, UB};
            for (int i = 0; i < 2; i++) Rt[2][i] = findroot(2, br[i], br[i + 1]);
            rec(3);
        }
    }
    fprintf(stderr, "a1=%ld a2=%ld S=%.4Lf mode=%s leaves=%lld F1=%lld F2=%lld F3=%lld NS=%lld RIG=%lld ybud=%lld rung=%lld out=%lld pruned=%lld\n",
            a1x, a2x, SX, mode == 1 ? "YL" : "YU", n_leaf, n_f1, n_f2, n_f3, n_ns, n_rig, n_ybud, n_rung, n_out, n_prune);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s a1 a2\n", argv[0]); return 1; }
    for (int n = 0; n <= N; n++)
        for (int k = 0; k <= n; k++)
            Cb[n][k] = (k == 0 || k == n) ? 1 : Cb[n - 1][k - 1] + Cb[n - 1][k];
    /* d*nu_max with d*nu >= 2(nu-1)^2:  nu_max = ((d+4)+sqrt((d+4)^2-16))/4 */
    SMAX = (ld)N * ((N + 4) + sqrtl((ld)(N + 4) * (N + 4) - 16)) / 4 + 1e-9L;
    { const char *e = getenv("ENUM_A3MOD"); if (e) sscanf(e, "%ld/%ld", &A3K, &A3N); }
    if (getenv("ENUM_NOF3")) NOF3 = 1;
    a[0] = 1;
    a[1] = atol(argv[1]);
    a[2] = atol(argv[2]);
    /* optional: argv[3] = SMAX override (0 = default),
     *           argv[4] = "all" (no filters) or "strong" (add the NS/RIG tests) */
    if (argc > 3 && strtold(argv[3], NULL) > 0) SMAX = strtold(argv[3], NULL);
    if (argc > 4) {
        if (!strcmp(argv[4], "strong")) STRONG = 1;
        else if (!strcmp(argv[4], "why")) STRONG = WHY = 1;          /* print the reason for each exclusion */
        else if (!strcmp(argv[4], "sharp")) STRONG = SHARP = 1;      /* strong + kappa >= 0 */
        else if (!strcmp(argv[4], "whysharp")) STRONG = SHARP = WHY = 1;
        else NOFILTER = 1;
    }
    ld S = (ld)(N - 1) * a[1] * a[1] / N - 2 * (ld)a[2];
    if (S <= 0 || S > SMAX) return 0;
    TAU = -(ld)a[1] / N;
    VAR = S / N;
    NU = VAR;
    if (argc > 5) {
        ld mm = 1.5L * VAR < VAR + 0.75L ? 1.5L * VAR : VAR + 0.75L;
        if (!strcmp(argv[5], "X")) {
            XSPLIT = 1;
            XA = floorl(TAU - mm + 1e-12L) + 1;
            XB = ceill(TAU + mm - 1e-12L) - 1;
        } else if (!strcmp(argv[5], "YL") || !strcmp(argv[5], "YU")) {
            return yrun(argv[5][1] == 'L' ? 1 : 2, S, mm);
        }
    }
    ld rad = sqrtl((N - 1) * S / N);
    /* integers L = floor(tau - V + 1), U = ceil(tau + V - 1) (rounded conservatively) */
    AA = TAU - floorl(TAU - VAR + 1 + 1e-12L);
    BB = ceill(TAU + VAR - 1 - 1e-12L) - TAU;
    if (VAR <= 1) AA = BB = 0;
    if (XSPLIT) {
        if (TAU - XA > AA) AA = TAU - XA;
        if (XB - TAU > BB) BB = XB - TAU;
    }
    if (!NOFILTER && AA > 0 && BB > 0) {
        ld r2 = S - AA * AA - BB * BB;
        if (r2 < -1e-9L) {       /* no leaf can pass F1 and F2 */
            fprintf(stderr, "a1=%ld a2=%ld S=%.4Lf leaves=0 F1=0 F2=0 F3=0 NS=0 RIG=0 out=0 pruned=1 (empty by F1/F2)\n", a[1], a[2], S);
            return 0;
        }
        MIDR = sqrtl(r2 > 0 ? r2 : 0) * (1 + 1e-12L) + 1e-12L;
    }
    LB = TAU - rad - 0.01L;
    UB = TAU + rad + 0.01L;
    Rt[1][0] = TAU;
    if (N == 2) {
        leaf();
    } else {
        /* roots of g_2 */
        ld br[3] = {LB, TAU, UB};
        for (int i = 0; i < 2; i++) Rt[2][i] = findroot(2, br[i], br[i + 1]);
        rec(3);
    }
    fprintf(stderr, "a1=%ld a2=%ld S=%.4Lf%s leaves=%lld F1=%lld F2=%lld F3=%lld NS=%lld RIG=%lld rung=%lld out=%lld pruned=%lld\n",
            a[1], a[2], S, XSPLIT ? " mode=X" : "", n_leaf, n_f1, n_f2, n_f3, n_ns, n_rig, n_rung, n_out, n_prune);
    return 0;
}
