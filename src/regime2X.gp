\\ regime2X.gp -- sextic Case A, regime nu_2 < 2 nu, both rungs of x_1 regular (nu <= 4.79).
\\ Candidates: output of src/run_regime2X.sh (F1, F2 on, F3 off).  Field-level filters:
\\   UNR  2 unramified (Theorem unram)          R1   nu_2 >= 2 nu (regime 1: Proposition reg1)
\\   NS   nu(K) < Var(x_1): x_1 not minimal
\\   BC0  some alpha in O_K^+ \ Z has tau(alpha) < min(nu, nu_2/2)   (Corollary dich(1),(4))
\\   HS   some alpha with tau(alpha) < nu has no half-square 2 alpha - c in O_K^2 (Corollary dich(2))
\\   W    a certified witness (if all filters pass)
read("src/gtest.gp");
varpol(f) = my(d = poldegree(f), a1 = polcoef(f, d-1), a2 = polcoef(f, d-2)); ((d-1)*a1^2/d - 2*a2)/d;
fieldtest2(g) =
{
  my(D = nfdisc(g), K, nu, nu2, L);
  if (D % 2, return("UNR"));
  K = gfield(g, 0); nu = nuK(K)[1]; nu2 = nu2K(K);
  if (nu2 >= 2*nu, return(Str("R1 nu=", nu)));
  L = tpos_small(K, nu - 1/10^9);
  for (i = 1, #L,
    my(a = L[i], ta = tau(K, a), m = vecmin(emb(K, a)), ok = 0);
    if (ta >= nu, next);
    if (ta < min(nu, nu2/2), return(Str("BC0 nu=", nu)));     \\ alpha is then a witness (Cor. dich(4))
    \\ half-square: 2 alpha - c = y^2 with 0 <= c < 2 min sigma(alpha) and c <= 2 tau(alpha) - nu_2
    for (c = 0, min(floor(2*m - 1e-30), floor(2*ta - nu2)), if (type(sqrtK(K, 2*a - c*unitvec(6, 1))) == "t_COL", ok = 1; break));
    if (!ok, return(Str("HS nu=", nu))));
  my(KB = gfield(g), w = gwitness2(KB, 300, 3));
  if (w === 0, Str("SURVIVOR nu=", nu), Str("W ", w[1], " N=", w[3], " nu=", nu));
}
r2x(infile, outfile, part, nparts) =
{
  my(L = readstr(infile), seen = Map(), k = 0);
  for (i = 1, #L,
    my(c = eval(Str("[", strjoin(strsplit(L[i], " "), ","), "]")), f, g, V, res);
    if (#c != 6, next);
    f = Pol(concat([1], c), 'x);
    if (!polisirreducible(f), next);
    k++; if ((k - 1) % nparts != part, next);
    V = varpol(f); g = polredabs(f);
    my(key = [g, V]);
    if (mapisdefined(seen, key, &res), write(outfile, Str(f, "|", res)); next);
    \\ NS is decided per (field, Var(x_1)); the remaining tests only depend on the field
    my(nuK0 = nuK(gfield(g, 0))[1]);
    res = if (nuK0 < V, "NS", fieldtest2(g));
    mapput(seen, key, res);
    write(outfile, Str(f, "|", res)));
}
