\\ gtest.gp -- exact tools for the class
\\      G = { K totally real : 2 O_K^+ is contained in Sigma O_K^2 }.
\\
\\ Conventions.  A field is given by a monic integral polynomial; internally the
\\ number field lives in the variable w (so that polynomials over K can use x).
\\ Elements of O_K are column vectors on the integral basis nf.zk (zk[1] = 1).
\\
\\ Everything that is *decided* here is decided exactly:
\\   * total positivity: an algebraic integer gamma != 0 has |sigma_j(gamma)| >= H^(1-d),
\\     H = house(gamma); with realprecision 100 and |error| < 1e-80 a numerical sign is
\\     therefore exact as long as H^(d-1) < 1e60.  Ties (|value| < 1e-60) fall back to
\\     Sturm sequences on the characteristic polynomial.
\\   * sums of squares: sosdecide() enumerates ALL x in O_K with x^2 <= beta
\\     (Fincke-Pohst on the weighted form sum_j sigma_j(x)^2/sigma_j(beta) <= d, then an
\\     exact check), and searches all multisets with a memo on the remainder.  The search
\\     is a decision procedure: every summand of a representation is in that finite set.
\\ Floating point is only used to locate candidates and to order the search.

default(realprecision, 100);
\\ run gp with  -D parisizemax=2G  (setting it inside a read file aborts the read)

\\ ------------------------------------------------------------------ field record
\\ K = [nf, bnf or 0, d, M, fu (coords), sgnfu (F2 matrix), T (trace form), Gvar]
gfield(P, wantbnf = 1) =
{
  my(Q = subst(P, variable(P), 'w), nf, bnf = 0, d = poldegree(P), M, fu = [], S = 0, T, G);
  if (wantbnf, bnf = bnfinit(Q, 1); nf = bnf.nf, nf = nfinit(Q));
  if (nf.sign[1] != d, error("gfield: not totally real"));
  M = nf[5][1];
  if (wantbnf,
    fu = apply(u -> nfalgtobasis(nf, u), bnf.fu);
    fu = apply(u -> if (type(u) == "t_COL", u, u * unitvec(d, 1)), fu);
    fu = lllunits(nf, fu, d);
    \\ sign matrix over F2: columns = sign vectors of -1, fu_1, ..., fu_{d-1}
    S = matrix(d, d, i, j, if (j == 1, 1, signexact(nf, M, fu[j-1], i) < 0)));
  if (nf.roots != vecsort(nf.roots), error("gfield: nf.roots not sorted"));
  T = matrix(d, d, i, j, nfelttrace(nf, nfeltmul(nf, unitvec(d, i), unitvec(d, j))));
  \\ d^2 Var on the basis zk[2..d] of Lambda = O/Z
  G = matrix(d-1, d-1, i, j, d*T[i+1,j+1] - T[1,i+1]*T[1,j+1]);
  [nf, bnf, d, M, fu, S, T, G, vecmax(abs(M))];
}
\\ LLL-reduce a basis of fundamental units w.r.t. the log embedding, so that the units
\\ used later are balanced (bnfinit may return units with huge coefficients)
lllunits(nf, fu, d) =
{
  my(r = #fu, th, Lm, U, nfu);
  if (r == 0, return(fu));
  localprec(600);
  th = real(polroots(nf.pol));
  Lm = matrix(d, r, j, i, log(abs(subst(lift(nfbasistoalg(nf, fu[i])), 'w, th[j]))));
  U = qflll(Lm);
  nfu = vector(r, i, my(u = vector(d)~); u[1] = 1;
    for (j = 1, r, if (U[j,i], u = nfeltmul(nf, u, nfeltpow(nf, fu[j], U[j,i]))));
    if (type(u) == "t_COL", u, u * concat([1], vector(d-1))~));
  nfu;
}
gnf(K) = K[1];
gbnf(K) = K[2];
gdeg(K) = K[3];
unitvec(d, i) = my(v = vector(d)~); v[i] = 1; v;

emb(K, x) = K[4] * x;
tr(K, x) = nfelttrace(K[1], x);
nrm(K, x) = nfeltnorm(K[1], x);
\\ PARI returns rational results of nf arithmetic as scalars: force columns
col(K, a) = if (type(a) == "t_COL", a, a * unitvec(K[3], 1));
mul(K, x, y) = col(K, nfeltmul(K[1], x, y));
sq(K, x) = col(K, nfeltmul(K[1], x, x));
pw(K, x, n) = col(K, nfeltpow(K[1], x, n));
alg(K, x) = nfbasistoalg(K[1], x);
bas(K, a) = col(K, nfalgtobasis(K[1], a));
cpoly(K, x) = charpoly(Mod(lift(nfbasistoalg(K[1], x)), K[1].pol), 'x);
minpol(K, x) = minpoly(Mod(lift(nfbasistoalg(K[1], x)), K[1].pol), 'x);
isrational(K, x) = (x[2..#x] == vector(#x-1)~);

\\ Rigorous numerics.  nf.M is computed to ~100 digits, so the computed value of
\\ sigma_j(x) = (M x)_j has absolute error < embtol(K,x) = (|x|_1 + 1) max|M| 10^-90.
\\ A sign is read off numerically only when |value| exceeds that bound; otherwise
\\ the decision is exact (Sturm sequence of the characteristic polynomial).
embtol(K, x) = (normlp(x, 1) + 1) * K[9] * 1e-90;
\\ exact sign test: 1 if x is totally positive, 0 otherwise (x a nonzero integer of K)
istp(K, x) =
{
  my(v = emb(K, x), m = vecmin(v), e = embtol(K, x));
  if (m > e, return(1));
  if (m < -e, return(0));
  my(cp = cpoly(K, x));
  if (polcoef(cp, 0) == 0, return(0));
  polsturm(cp, [-oo, 0]) == 0;
}
\\ certified high-precision embeddings: returns [v, err] with |sigma_j(x) - v_j| <= err,
\\ raising the precision until every |v_j| > err (x != 0).  The roots are computed by
\\ polroots at precision pr (error << 10^(-pr+10)); the evaluation error is bounded by
\\ |a|_1 * d * (1+R)^d * 10^(-pr+10), a = polynomial of x, R = max |root|.
embhp(K, x) =
{
  my(nf = K[1], P = nf.pol, d = K[3], a = lift(nfbasistoalg(nf, x)), A, pr = 120);
  A = normlp(Vec(a), 1) + 1;
  while (1,
    localprec(pr);
    my(r = vecsort(real(polroots(P))), R = vecmax(abs(r)), err, v);
    \\ match the ordering of nf.roots (sorted as well)
    err = A * d * (1 + R)^d * 10^(-pr + 10);
    v = vector(d, j, subst(a, variable(P), r[j]));
    if (vecmin(abs(v)) > err, return([v, err]));
    pr *= 2);
}
\\ exact sign of sigma_i(x) for x != 0, embeddings ordered as nf.roots (ascending)
signexact(nf, M, x, i) =
{
  my(v = (M[i,] * x), e = (normlp(x, 1) + 1) * vecmax(abs(M)) * 1e-90);
  if (abs(v) > e, return(sign(v)));
  my(K0 = [nf, 0, poldegree(nf.pol)]);
  sign(embhp(K0, x)[1][i]);
}
\\ totally >= 0 (allows x = 0)
istnn(K, x) = if (x == vector(#x)~, 1, istp(K, x));

\\ ------------------------------------------------------------------ squares in K
\\ is a (coords) a square in K?  returns a square root (coords) or 0
sqrtK(K, a) =
{
  my(nf = K[1], r);
  if (a == vector(#a)~, return(a));
  r = nfroots(nf, 'x^2 - nfbasistoalg(nf, a));
  if (#r == 0, 0, col(K, nfalgtobasis(nf, r[1])));
}

\\ ------------------------------------------------------------------ sums of squares
\\ all x != 0 (one of +-x) in O_K with beta - x^2 totally >= 0
sqcands(K, beta) =
{
  my(d = K[3], M = K[4], v = emb(K, beta), G, q, Xc = List(), Xs = List(), pr, ok = 0);
  if (!istp(K, beta), error("sqcands: beta not totally positive"));
  \\ the weighted Gram sum_j m_j m_j^t / sigma_j(beta) is ill-conditioned when beta is
  \\ unbalanced: work at a precision adapted to the spread of the conjugates
  pr = 100 + 3 * ceil(log(normlp(beta, 1) + 2) / log(10)) + 3 * ceil(d * log(K[9] + 2) / log(10));
  while (!ok,
    localprec(pr);
    my(nfh = nfnewprec(K[1]), Mh = nfh[5][1], vh);
    vh = Mh * beta;
    if (vecmin(vh) <= 0, pr *= 2; next);
    G = Mh~ * matdiagonal(vector(d, j, 1/vh[j])) * Mh;
    G = (G + G~)/2;
    q = iferr(my(U = qflllgram(G)); U * qfminim(U~ * G * U, d + 1e-20, , 2)[3], E, 0);
    if (q === 0, pr *= 2, ok = 1));
  for (k = 1, #q,
    my(c = q[,k], s = sq(K, c));
    if (istnn(K, beta - s), listput(Xc, c); listput(Xs, s)));
  [Vec(Xc), Vec(Xs)];
}

\\ decide beta in Sigma O^2.  Returns [1, [x_1,...,x_N]] or [0, #candidates].
sosdecide(K, beta) =
{
  my(C, Xc, Xs, Xe, n, trs, perm, sqmap = Map(), memo = Map(), res);
  C = sqcands(K, beta); Xc = C[1]; Xs = C[2]; n = #Xc;
  trs = vector(n, i, tr(K, Xs[i]));
  perm = vecsort(trs, , 5);   \\ indices, decreasing trace of x^2
  Xc = vector(n, i, Xc[perm[i]]); Xs = vector(n, i, Xs[perm[i]]);
  Xe = vector(n, i, emb(K, Xs[i]));
  my(Xt = vector(n, i, embtol(K, Xs[i])));
  for (i = 1, n, mapput(sqmap, Xs[i], i));
  res = sosrec(K, Xs, Xe, Xt, sqmap, memo, beta);
  if (type(res) != "t_VEC", return([0, n]));
  res = vector(#res, k, Xc[res[k]]);
  \\ certify the representation exactly
  if (sum(k = 1, #res, sq(K, res[k])) != beta, error("sosdecide: bad representation"));
  [1, res];
}

\\ returns a vector of candidate indices whose squares sum to rem, or 0
sosrec(K, Xs, Xe, Xt, sqmap, memo, rem) =
{
  my(i, r2, v, e, res, n = #Xs, d = K[3], z = vector(d)~);
  if (rem == z, return([]));
  if (mapisdefined(sqmap, rem, &i), return([i]));
  if (mapisdefined(memo, rem), return(0));
  v = emb(K, rem); e = embtol(K, rem);
  \\ a candidate is skipped only if x^2 <= rem fails provably (numeric margin)
  \\ two squares first (cheap), then recursion
  for (j = 1, n,
    if (vecmin(v - Xe[j]) < -(e + Xt[j]), next);
    r2 = rem - Xs[j];
    if (mapisdefined(sqmap, r2, &i), return([j, i])));
  for (j = 1, n,
    if (vecmin(v - Xe[j]) < -(e + Xt[j]), next);
    r2 = rem - Xs[j];
    if (r2 == z, return([j]));
    if (!istp(K, r2), next);
    res = sosrec(K, Xs, Xe, Xt, sqmap, memo, r2);
    if (type(res) == "t_VEC", return(concat([j], res))));
  mapput(memo, rem, 0);
  0;
}

\\ ------------------------------------------------------------------ units
\\ F2 kernel of the sign matrix -> totally positive units modulo squares
\\ returns list of units (coords) forming an F2-basis of U^+/U^2
tpunits(K) =
{
  my(S = K[6], d = K[3], ker, L = List());
  ker = lift(matker(Mod(S, 2)));
  for (k = 1, #ker,
    my(e = ker[,k], u = unitvec(d, 1));
    if (e[1] % 2, u = -u);
    for (i = 2, d, if (e[i] % 2, u = mul(K, u, K[5][i-1])));
    listput(L, u));
  Vec(L);
}

\\ all products of subsets of a list (F2-span), coords
spanprod(K, B) =
{
  my(L = List());
  listput(L, unitvec(K[3], 1));
  for (i = 1, #B, my(m = #L); for (k = 1, m, listput(L, mul(K, L[k], B[i]))));
  Vec(L);
}

\\ unit criterion: returns 0 if passed, else a totally positive unit u with
\\ neither u nor 2u a square (then 2u is not a sum of squares: Lemma unitdich).
unitcrit(K) =
{
  my(B = tpunits(K), L = spanprod(K, B));
  for (k = 2, #L,
    my(u = L[k]);
    if (sqrtK(K, u) === 0 && sqrtK(K, 2*u) === 0, return(u)));
  0;
}
\\ h^+/h = |U^+/U^2|
hplusoverh(K) = 2^#tpunits(K);

\\ ------------------------------------------------------------------ generators
\\ log|sigma_j(x)| for a list of elements, at high precision (columns = elements).
\\ Only used to choose representatives, never to decide anything.
logemb(K, L) =
{
  my(d = K[3], E = vector(#L, i, embhp(K, L[i])[1]));
  matrix(d, #L, j, i, log(abs(E[i][j])));
}
\\ reduce a totally positive alpha modulo squares of units, to make its trace small
unitreduce(K, a) =
{
  my(d = K[3], fu = K[5], r = d - 1, Lg, la, c, best, bt);
  if (r == 0, return(a));
  Lg = 2 * logemb(K, fu);
  la = logemb(K, [a])[,1];
  la = la - vecsum(la)/d * vector(d, j, 1)~;
  c = matsolve(Lg~*Lg, -Lg~*la);
  c = round(c);
  best = a; bt = tr(K, a);
  my(apply_c = cc -> my(b = a); for (i = 1, r, if (cc[i], b = mul(K, b, pw(K, fu[i], 2*cc[i])))); b);
  my(b = apply_c(c), tb = tr(K, b));
  if (tb < bt, best = b; bt = tb);
  \\ local search
  my(improved = 1, cc = c);
  while (improved, improved = 0;
    for (i = 1, r, for (s = -1, 1, if (s == 0, next);
      my(c2 = cc); c2[i] += s;
      my(b2 = apply_c(c2), t2 = tr(K, b2));
      if (t2 < bt, best = b2; bt = t2; cc = c2; improved = 1))));
  best;
}

\\ totally positive generators of the ideal I modulo U^2 (empty if not narrowly principal)
tpgens(K, I) =
{
  my(bnf = K[2], d = K[3], r, g, s, e, a, B, L);
  r = bnfisprincipal(bnf, I);
  if (r[1] != vector(#r[1])~, return([]));
  g = col(K, r[2]);
  s = vector(d, i, signexact(K[1], K[4], g, i) < 0)~;
  e = matsolvemod(K[6], 2, s);
  if (e === 0, return([]));
  a = g;
  if (e[1] % 2, a = -a);
  for (i = 2, d, if (e[i] % 2, a = mul(K, a, K[5][i-1])));
  if (!istp(K, a), error("tpgens: sign bug"));
  B = tpunits(K);
  L = spanprod(K, B);
  apply(u -> unitreduce(K, mul(K, a, u)), L);
}

\\ ------------------------------------------------------------------ witness search
\\ search for alpha in O^+ with 2 alpha not a sum of squares, by increasing norm.
\\ returns [type, alpha, norm, info] or 0.  type: "unit" or "norm".
gwitness(K, Nmax, verbose = 0) =
{
  my(u = unitcrit(K), nf = K[1], IL);
  if (type(u) == "t_COL", return(["unit", u, 1, 0]));
  IL = ideallist(nf, Nmax);
  for (N = 2, Nmax,
    my(Is = IL[N]);
    for (k = 1, #Is,
      my(gens = tpgens(K, Is[k]));
      for (j = 1, #gens,
        my(a = gens[j], r = sosdecide(K, 2*a));
        if (verbose, print("N=", N, " alpha=", alg(K, a), " -> ", r[1]));
        if (r[1] == 0, return(["norm", a, N, r[2]])))));
  0;
}

\\ ------------------------------------------------------------------ variance invariants
\\ nu(K) as [nu, vectors in zk coords (first coord 0)]
nuK(K) =
{
  my(d = K[3], q = qfminim(K[8]), vs);
  vs = vector(#q[3], k, concat([0], q[3][,k]~)~);
  [q[2]/d^2, vs];
}
\\ dyadic square-root radical r = prod p^ceil(e/2)
dyadicrad(K) =
{
  my(nf = K[1], P = idealprimedec(nf, 2), r = 1);
  for (i = 1, #P, r = idealmul(nf, r, idealpow(nf, P[i], ceil(P[i].e/2))));
  r;
}
\\ nu_2(K) = min Var on (Z + r) \ Z
nu2K(K) =
{
  my(d = K[3], r = dyadicrad(K), B, H, Gs);
  B = idealhnf(K[1], r);
  B = mathnf(matrix(d-1, d, i, j, B[i+1, j]));
  Gs = B~ * K[8] * B;
  qfminim(Gs)[2] / d^2;
}
Var(K, x) = my(d = K[3], y = x[2..d]); (y~ * K[8] * y) / d^2;
tau(K, x) = tr(K, x) / K[3];

\\ t(K) = min tau over O^+ \ Z, by enumerating Lambda up to d^2 Var <= d^2 (d-1) T^2,
\\ with the optimal integer shift for each class.  T = an upper bound (e.g. from tsearch0)
tK(K, T) =
{
  my(d = K[3], R = d^2*(d-1)*T^2, q, best = T, arg = 0);
  q = qfminim(K[8], R + 1e-9, , 2)[3];
  for (k = 1, #q,
    for (s = -1, 1, if (s == 0, next);
      my(w = s * concat([0], q[,k]~)~, v = emb(K, w), n = floor(-vecmin(v)) + 1, a = w, ta);
      a[1] += n;
      if (!istp(K, a), a[1] += 1);
      if (!istp(K, a), error("tK: shift bug"));
      \\ smallest shift: decrease while still totally positive
      while (1, my(b = a); b[1] -= 1; if (istp(K, b), a = b, break));
      ta = tau(K, a);
      if (ta < best || (ta == best && arg === 0), best = ta; arg = a)));
  [best, arg];
}

\\ ------------------------------------------------------------------ small trace
\\ all totally positive non-rational alpha with tau(alpha) <= T (one per element, not up to
\\ units).  Uses Var(alpha) <= (d-1) tau(alpha)^2 and the optimal integer shift per class.
tpos_small(K, T) =
{
  my(d = K[3], R = d^2*(d-1)*T^2, q, L = List());
  q = qfminim(K[8], R + 1e-9, , 2)[3];
  for (k = 1, #q,
    for (s = -1, 1, if (s == 0, next);
      my(w = s * concat([0], q[,k]~)~, v = emb(K, w), a = w);
      a[1] += floor(-vecmin(v)) + 1;
      while (!istp(K, a), a[1] += 1);
      while (1, my(b = a); b[1] -= 1; if (istp(K, b), a = b, break));
      while (tau(K, a) <= T, listput(L, a); a[1] += 1)));
  L = Vec(L);
  vecsort(L, a -> tau(K, a));
}

\\ certification of a witness, independent of how it was found: alpha is totally
\\ positive (Sturm sequence of its characteristic polynomial, no floating point) and
\\ the decision procedure, rerun from scratch, finds no representation of 2 alpha.
certify(K, a) =
{
  my(cp = cpoly(K, a));
  if (polcoef(cp, 0) == 0 || polsturm(cp, [-oo, 0]) != 0, error("certify: not totally positive"));
  if (sosdecide(K, 2*a)[1] != 0, error("certify: 2 alpha IS a sum of squares"));
  1;
}

\\ witness search: units, then small trace (tau <= Tmax), then small norm (<= Nmax)
gwitness2(K, Nmax, Tmax) =
{
  my(u = unitcrit(K), L, w = 0);
  if (type(u) == "t_COL", w = ["unit", u, 1, 0]);
  if (w === 0 && Tmax > 0,
    L = tpos_small(K, Tmax);
    for (i = 1, #L,
      my(r = sosdecide(K, 2*L[i]));
      if (r[1] == 0, w = ["trace", L[i], abs(nrm(K, L[i])), r[2]]; break)));
  if (w === 0, w = gwitness(K, Nmax));
  if (type(w) == "t_VEC", certify(K, w[2]));
  w;
}
