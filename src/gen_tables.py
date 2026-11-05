"""gen_tables.py -- write the LaTeX tables of notes/G_beyond5.tex from results/."""
import collections, glob, json, os, re
from fractions import Fraction
R = 'results'; G = 'notes/gen'
os.makedirs(G, exist_ok=True)

def sci(n):
    m, e = f"{n:.2e}".split('e')
    return f"${m}\\cdot10^{{{int(e)}}}$"

def pol(s):
    return s.replace('*', '')

def census_rows(d):
    fs = sorted(glob.glob(f'{R}/census/deg{d}_N*.txt'))
    fs = [f for f in fs if 'survivors' not in f and 'part' not in f]
    if not fs: return None, None
    rows = [l.rstrip('\n').split('|') for l in open(fs[0])]
    surv = {}
    for f in glob.glob(f'{R}/census/deg{d}_survivors_*.txt'):
        for l in open(f):
            p = l.rstrip('\n').split('|'); surv[p[0]] = p
    return rows, surv

def tex_census():
    lines = [r'\begin{tabular}{rrrrrrrrr}', r'\toprule',
             r'$d$ & fields & $\max d_K$ & U & T & N & left & $h^+=h$ & max $N$\\', r'\midrule']
    thm = []
    for d in range(3, 9):
        rows, surv = census_rows(d)
        if rows is None: continue
        st = collections.Counter(r[6] for r in rows)
        # survivors finished with larger bound
        left = [r for r in rows if r[6] == 'S']
        nleft = 0; members = []
        for r in left:
            p = surv.get(r[0])
            if p is None or p[2].startswith('SURVIVOR'):
                nleft += 1; members.append(r[0])
            else:
                st[{'unit': 'U', 'trace': 'T', 'norm': 'N'}[p[2]]] += 1
        err = st.get('E', 0)
        hp = sum(1 for r in rows if r[3] == '1')
        nmax = max([int(r[8]) for r in rows if r[6] == 'N'] + [0])
        for r in left:
            p = surv.get(r[0])
            if p is not None and p[2] == 'norm':
                nmax = max(nmax, int(p[4].split('=')[1]))
        dmax = max(int(r[0]) for r in rows)
        lines.append(f"{d} & {len(rows)} & {sci(dmax)} & {st['U']} & {st['T']} & {st['N']} & {nleft}"
                     + (f" ({','.join(members)})" if members else '') + (f" + {err} err" if err else '')
                     + f" & {100*hp/len(rows):.1f}\\% & {nmax if nmax else '--'}\\\\")
        thm.append((d, len(rows), dmax, nleft, members, err))
    lines += [r'\bottomrule', r'\end{tabular}']
    open(f'{G}/census_table.tex', 'w').write('\n'.join(lines) + '\n')
    return thm

def tex_unram():
    out = [r'\begin{tabular}{rll}', r'\toprule', r'$d_K$ & witness $\alpha$ (minimal polynomial) & $\Nm_{K/\Q}\alpha$\\', r'\midrule']
    for f in ['g2o_d6.txt']:
        for l in open(f'{R}/{f}'):
            p = l.rstrip('\n').split('|')
            if p[2] != '2unram=1': continue
            out.append(f"{p[0]} & ${pol(p[4])}$ & {p[5][2:]}\\\\")
    out.append(r'\midrule')
    for l in open(f'{R}/g2o_d8B.txt'):
        p = l.rstrip('\n').split('|')
        if p[0] in ('442050625', '656505625'):
            out.append(f"{p[0]} & ${pol(p[4])}$ & {p[5][2:]}\\\\")
    out += [r'\bottomrule', r'\end{tabular}']
    open(f'{G}/unram_table.tex', 'w').write('\n'.join(out) + '\n')

def tex_spectrum(vmax=Fraction(19, 10)):
    vals = collections.defaultdict(lambda: [set(), 0])
    for f in glob.glob(f'{R}/spectrum/spec_d*.txt'):
        for l in open(f):
            p = l.split('|'); v = Fraction(p[1])
            if v < vmax: vals[v][0].add(int(p[0])); vals[v][1] += 1
    out = [r'\begin{tabular}{llll}', r'\toprule', r'$\Var$ & decimal & degrees & \#\\', r'\midrule']
    for v in sorted(vals):
        out.append(f"${v.numerator}/{v.denominator}$ & {float(v):.4f} & {','.join(map(str, sorted(vals[v][0])))} & {vals[v][1]}\\\\")
    out += [r'\bottomrule', r'\end{tabular}']
    open(f'{G}/spectrum_table.tex', 'w').write('\n'.join(out) + '\n')

if __name__ == '__main__':
    print(tex_census()); tex_unram(); tex_spectrum()

import numpy as np, math

def tex_census_theorem(thm):
    by = {d: (n, dm, nl, mem, err) for d, n, dm, nl, mem, err in thm}
    hi = [d for d in (6, 7, 8) if d in by]
    parts = []
    for d in hi:
        n, dm, nl, mem, err = by[d]
        parts.append(f"{n} fields of degree {d} ($d_K\\le{sci(dm)[1:-1]}$)")
    t = (r"Every one of the " + ", ".join(parts) +
         r" in the LMFDB lists has an explicit element $\alpha\in\OK^{\tp}$ with $2\alpha\notin\SoS$,"
         r" certified as in \S\ref{sec:algo}; none of these fields lies in $\Gcl$. In degrees "
         + ", ".join(str(d) for d in (3, 4, 5) if d in by) +
         r" the fields without a witness are exactly the members of \eqref{eq:G5} in the lists.")
    bad = [d for d in hi if by[d][2] or by[d][4]]
    if bad:
        t += r" \emph{(Exceptions: " + ", ".join(f"degree {d}: {by[d][2]} unresolved, {by[d][4]} errors" for d in bad) + ".)}"
    open(f'{G}/census_theorem.tex', 'w').write(t + '\n')

GAMMA = {3: 2/math.sqrt(3), 4: 2**(1/3), 5: math.sqrt(2), 6: 8**(1/5), 7: (64/3)**(1/6), 8: 64**(1/7)}

def tex_stats():
    out = []
    out.append(r"\begin{table}[t]\caption{Normalised shape minimum $\hat\nu$ in the census: quantiles,"
               r" the Hermite bound $\gamma_{d-1}$, medians by type of first witness, and the share of"
               r" fields with $h^+=h$ (equivalently $U_K^+=U_K^2$).}\label{tab:shape}")
    out.append(r"\begin{tabular}{rrrrrrrrrrr}\toprule")
    out.append(r"$d$ & $q_{0.01}$ & $q_{0.1}$ & $q_{0.5}$ & $q_{0.9}$ & max & $\gamma_{d-1}$ & U & T & N & $h^+=h$\\\midrule")
    cub = None
    for d in range(3, 9):
        rows, _ = census_rows(d)
        if rows is None: continue
        shp = np.array([float(r[13]) for r in rows if r[13]])
        if len(shp) == 0: continue
        med = {k: np.median([float(r[13]) for r in rows if r[6] == k and r[13]]) if any(r[6] == k for r in rows) else float('nan') for k in 'UTN'}
        hp = sum(1 for r in rows if r[3] == '1') / len(rows)
        q = np.quantile(shp, [0.01, 0.1, 0.5, 0.9])
        fm = lambda x: '--' if x != x else f"{x:.2f}"
        out.append(f"{d} & " + " & ".join(f"{x:.2f}" for x in q) + f" & {shp.max():.2f} & {GAMMA[d]:.2f} & "
                   + " & ".join(fm(med[k]) for k in 'UTN') + f" & {100*hp:.1f}\\%\\\\")
        if d == 3: cub = (rows, shp)
    out.append(r"\bottomrule\end{tabular}\end{table}")
    out.append("")
    out.append(r"""Table~\ref{tab:shape} summarises the shapes. Three features stand out.
(a) In the prime degrees $5$ and $7$ the shape minima concentrate in a band around $0.7$--$1.0$, well
below the Hermite bound, while in degrees $4$, $6$, $8$ there is a heavy lower tail: these are the
imprimitive fields, whose minimal vector lies in a subfield of small variance and whose shape is
therefore skewed. (b) In degree $6$, witnesses of type T (small trace, i.e.\ $t(K)$ small
compared with $\nu$ and $\nu_2$: Corollary~\ref{cor:dich}(4)) occur for well-rounded shapes and
witnesses of type N for skewed ones: the medians of $\hat\nu$ are $0.77$ and $0.37$. This is the
budget--cost mechanism seen statistically: a round shape forces $\nu$ to be large compared with the
cheapest totally positive element, while a skewed (imprimitive) shape inherits a small $\nu$ from a
subfield and only a norm-type witness is left. In degree $4$ the distinction is not visible.
(c) Between one half and two thirds of the fields have $U_K^{+}\neq U_K^{2}$; in odd degree all of
these are excluded by the unit criterion, since then $(2)$ is not the square of an ideal and $2u$
can never be a square. The proportion of fields with $h^{+}=h$ that we observe ($31$--$48\%$ in
degrees $4$--$7$) is the kind of statistic predicted by the heuristics of Dummit and Voight for
signature ranks of units \cite{DV}; we have not compared the numbers with their predictions.""")
    if cub is not None:
        rows, shp = cub
        D = np.array([int(r[0]) for r in rows if r[13]], dtype=float)
        out.append("")
        lines = [r"\begin{table}[t]\caption{Cubic fields: empirical distribution of $\hat\nu=1/\operatorname{Im}z$"
                 r" (the shape $z$ in the modular fundamental domain) against the hyperbolic-measure"
                 r" prediction $P(\hat\nu<s)=3s/\pi$ $(s\le1)$ of Bhargava--Harron equidistribution.}\label{tab:bh}",
                 r"\begin{tabular}{rrrrr}\toprule $s$ & $3s/\pi$ & all & $d_K\le10^5$ & $d_K>10^6$\\\midrule"]
        for sv in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
            a = np.mean(shp < sv); b = np.mean(shp[D <= 1e5] < sv); c = np.mean(shp[D > 1e6] < sv)
            lines.append(f"{sv:.1f} & {3*sv/math.pi:.3f} & {a:.3f} & {b:.3f} & {c:.3f}\\\\")
        lines.append(r"\bottomrule\end{tabular}\end{table}")
        out += lines
        out.append(r"""For cubic fields the shape is a point $z$ of the modular curve and
$\hat\nu=1/\operatorname{Im}z$, so Bhargava--Harron equidistribution predicts $P(\hat\nu<s)=3s/\pi$ for
$s\le1$ (Table~\ref{tab:bh}). The deep cusp is, however, arithmetically forbidden: a generator $x$ of a
cubic field has $d_K\le\disc(x)\le(3\Var x)^{3}$, so $\nu(K)\ge d_K^{1/3}/3$ and
$\hat\nu\ge\sqrt3\,d_K^{-1/6}$. This lower cut-off is what makes the set of cubic fields with bounded
$\nu$ finite, although equidistribution alone would allow $\asymp\sqrt X$ of them below discriminant
$X$; it is the shape-theoretic form of the effective Northcott theorem of \cite{PB}.""")
    open(f'{G}/stats_body.tex', 'w').write('\n'.join(out) + '\n')

def tex_cyclo():
    rows = []
    for f in sorted(glob.glob(f'{R}/cyclo*.log')):
        for l in open(f):
            if not l.startswith('m='): continue
            p = dict(x.split('=', 1) for x in l.strip().split('|') if '=' in x)
            typ = [x for x in l.strip().split('|') if x in ('unit', 'trace', 'norm')]
            rows.append((int(p['m']), int(p['d']), typ[0] if typ else 'none', p.get('N', ''), p.get('Tr', '')))
    out = [r"\begin{center}\begin{tabular}{rrrrrl}\toprule $m$ & $m_0$ & $[F_m:\Q]$ & type & $\Nm\alpha$ & $\Tr\alpha$\\\midrule"]
    for m, d, t, N, T in sorted(rows):
        out.append(f"{m} & {m//4} & {d} & {t} & {N} & {T}\\\\")
    out.append(r"\bottomrule\end{tabular}\end{center}")
    open(f'{G}/cyclo_table.tex', 'w').write('\n'.join(out) + '\n')

if __name__ == '__main__':
    thm = tex_census(); tex_census_theorem(thm); tex_stats(); tex_cyclo()

def tex_regime1():
    out = []
    # cubic Case B
    if os.path.exists(f'{R}/cubicB.txt'):
        rows = [l.rstrip('\n').split('|') for l in open(f'{R}/cubicB.txt')]
        excl = sum(1 for r in rows if 'EXCL' in r[3])
        H = [r for r in rows if r[3].startswith('H(F)')]
        forced = [r for r in rows if r[3].startswith('forced')]
        nK = sum(int(r[3].split()[1]) for r in forced)
        surv = [r for r in forced if 'SURVIVOR' in r[4]]
        dmax = max(int(r[0]) for r in rows)
        out.append(r"\begin{proposition}\label{prop:cubicB}"
                   r"Let $K\in\Gcl$ be sextic, with a minimal vector of $\Lambda_K$ in a cubic subfield $F$."
                   r" Then $\nu(F)=\nu(K)\le9.58$ (both rungs cost at least $\frac12\nu$, and the spread $s$ of the"
                   r" minimal vector satisfies $s^{2}\le6\nu$), hence $d_F\le(3\nu)^{3}\le23{,}700$, and $t(F)\ge\frac12\nu(F)$. Among the "
                   + f"{len(rows)}" + r" cubic fields with $\nu(F)\le9.58$ (all have $d_F\le" + f"{dmax}" + r"$), "
                   + f"{excl}" + r" violate $t(F)\ge\frac12\nu(F)$; for "
                   + f"{len(forced)}" + r" of them some $\alpha\in\OF^{\tp}$ with $\tau(\alpha)<\nu(F)$ has no half-square in $F$, which forces"
                   r" $K=F(\sqrt{2\alpha-c})$ for an admissible $c$; the " + f"{nK}" + r" fields so obtained "
                   + (r"all have certified witnesses" if not surv else f"leave {len(surv)} survivors") +
                   r". The remaining cubic fields, for which $K$ is not determined, are exactly "
                   + ", ".join(f"$d_F={r[0]}$" for r in H) + r" (that is, $\Q(\zeta_7)^{+}$, $\Q(\zeta_9)^{+}$ and $\Q(\rho)$).\end{proposition}")
    if os.path.exists(f'{R}/regime1.txt'):
        rows = [l.rstrip('\n').split('|') for l in open(f'{R}/regime1.txt')]
        ram = [r for r in rows if r[2] == '2unram=0']
        surv = [r for r in rows if 'SURVIVOR' in r[3] or 'ERR' in r[3]]
        out.append(r"\begin{proposition}\label{prop:reg1}"
                   r"No sextic $K\in\Gcl$ is in regime~$1$ ($\nu_2\ge2\nu$) with a minimal vector generating $K$."
                   r"\end{proposition}" if not surv else r"\begin{proposition}\label{prop:reg1}[partial]\end{proposition}")
        out.append(r"\begin{proof}By the discussion above $K$ is generated by a root of one of the $9{,}785$ polynomials that"
                   r" \cite{TwoO} excludes only by rigidity. They define " + f"{len(rows)}" + r" distinct fields ("
                   + f"{len(ram)}" + r" with $2$ ramified), and each has a certified witness (file \texttt{results/regime1.txt})."
                   + (r"" if not surv else f" Exceptions: {len(surv)}.") + r"\end{proof}")
    open(f'{G}/regime1.tex', 'w').write('\n\n'.join(out) + '\n')

if __name__ == '__main__':
    tex_regime1()

def tex_inherited():
    rows = [l.rstrip('\n').split('|') for l in open(f'{R}/census/deg6_N300.txt')]
    lab = [('5/4', r'$\sqrt5\in K$ ($\nu=\frac54$)'), ('14/9', r'$\Q(\zeta_7)^+\subseteq K$ ($\nu=\frac{14}9$)'),
           ('2', r'$\nu=2$'), ('20/9', r'$\nu=\frac{20}9$'), ('3', r'$\nu=3$')]
    out = [r"\begin{table}[t]\caption{Sextic fields of the census whose $\nu$ is that of a member of $\Gcl$"
           r" (by Theorem~\ref{thm:smallnu} the first two rows are exactly the fields containing $\Q(\sqrt5)$, resp.~"
           r"$\Q(\zeta_7)^+$ but not $\sqrt5$): number of fields, types of first witness, and the norms of"
           r" the type-N witnesses.}\label{tab:inherited}",
           r"\begin{tabular}{lrrrrrr}\toprule family & fields & U & T & N & median $N$ & max $N$\\\midrule"]
    for k, name in lab:
        v = [r for r in rows if r[5] == k]
        typ = collections.Counter(r[6] for r in v)
        W = [int(r[8]) for r in v if r[6] == 'N']
        out.append(f"{name} & {len(v)} & {typ['U']} & {typ['T']} & {typ['N']} & "
                   + (f"{int(np.median(W))} & {max(W)}" if W else "-- & --") + r"\\")
    out.append(r"\bottomrule\end{tabular}\end{table}")
    out.append(r"""Table~\ref{tab:inherited} isolates the \emph{inherited} families, the sextic fields that contain
a member of $\Gcl$ (or have the same $\nu$ as one). They are the extreme cases of the census: they
pass the unit and trace tests far more often than average, and the fields of largest discriminant in
the $\Q(\zeta_7)^{+}$ family need the largest norm witnesses of all sextic fields ($\Nm\alpha=43$ at
$d_K\approx9\cdot10^{7}$). This is consistent with the reduction of \S\ref{sec:ramified}, where the
cases left open are exactly extensions of $\Q(\sqrt2),\Q(\sqrt3),\Q(\sqrt5),\Q(\zeta_7)^{+},
\Q(\zeta_9)^{+},\Q(\rho)$.""")
    return '\n'.join(out)

if __name__ == '__main__':
    open(f'{G}/stats_body.tex', 'a').write('\n' + tex_inherited() + '\n')

def tex_regime2X():
    fs = sorted(glob.glob(f'{R}/regime2X/status*.txt'))
    if not fs: return ''
    rows = [l.rstrip('\n').split('|', 1) for f in fs for l in open(f)]
    cnt = collections.Counter(r[1].split()[0] for r in rows)
    fields = {}
    surv = [r for r in rows if r[1].startswith('SURVIVOR')]
    t = (r"\begin{proposition}\label{prop:reg2X}No sextic $K\in\Gcl$ has a minimal vector $x_1$ generating $K$"
         r" with both rungs regular.\end{proposition}" if not surv else
         r"\begin{proposition}\label{prop:reg2X}[partial: survivors remain]\end{proposition}")
    t += (r"\begin{proof}If both rungs are regular, the downspread inequalities hold and $\nu\le4.79$, so the"
          r" enumeration of \cite{TwoO} without its quadratic filter (which needs $\tau\ge\nu$ for \emph{all}"
          r" totally positive $\alpha$, not available in regime~$2$) contains the minimal polynomial of $x_1$."
          r" This enumeration (\texttt{src/enum6g.c}, \texttt{src/run\_regime2X.sh}) has $178{,}438$ leaves,"
          r" exactly as in \cite{TwoO}, and outputs $42{,}232$ polynomials. Of the "
          + f"{len(rows)}" + r" irreducible ones, " +
          ", ".join(f"{v} are excluded by {k}" for k, v in sorted(cnt.items()) if k != 'SURVIVOR') +
          r", where UNR means $2$ unramified (Theorem~\ref{thm:unram}), R1 means $\nu_2\ge2\nu$"
          r" (Proposition~\ref{prop:reg1}), NS that $x_1$ is not minimal in $\OK$, BC0 an element with"
          r" $\tau<\min(\nu,\frac12\nu_2)$ (a witness by Corollary~\ref{cor:dich}(4)), HS an element with"
          r" $\tau<\nu$ and no half-square (Corollary~\ref{cor:dich}(2)), and W a certified witness."
          + ("" if not surv else f" {len(surv)} survive.") + r"\end{proof}")
    return t

if __name__ == '__main__':
    open(f'{G}/regime1.tex', 'a').write('\n' + tex_regime2X() + '\n')
