\\ cubicB.gp -- sextic K in G whose minimal vector lies in a cubic subfield F ([K:F] = 2).
\\ Necessary: nu(F) = nu(K) <= 9.58 (rungs cost >= nu/2, spread^2 <= 6 nu), so d_F <= (3 nu)^3;
\\ t(F) >= nu(F)/2; and every alpha in O_F^+ \ Z with tau(alpha) < nu(F) has 2 alpha - c a square
\\ in K (0 <= c < 2 min sigma(alpha)).  If some such alpha has no half-square in F, then
\\ K = F(sqrt(2 alpha - c)) for an admissible c: finitely many K, tested for G-witnesses.
read("src/gtest.gp");
NUB = 9.58;
cubB(infile, outfile, Dmax) =
{
  my(L = readstr(infile));
  for (i = 2, #L,
    my(p = strsplit(L[i], "|"), co = eval(Str("[", p[2], "]")), P = Pol(Vecrev(concat(co, [1])), 'x), D = nfdisc(P));
    if (D > Dmax, break);
    my(F = gfield(P, 0), nu = nuK(F)[1]);
    if (nu > NUB, next);
    if (D > (3*nu)^3, next);              \\ a generator x has d_F <= disc(x) <= (3 Var x)^3
    my(T = tpos_small(F, nu - 1/10^6), bad = 0, forced = 0, Ks = List());
    if (#T && tau(F, T[1]) < nu/2, write(outfile, Str(D, "|", P, "|nu=", nu, "|EXCL t<nu/2")); next);
    \\ find an alpha with no half-square in F; record the admissible quadratic extensions
    for (k = 1, #T,
      my(a = T[k], m = vecmin(emb(F, a)), cs = List(), has = 0);
      for (c = 0, floor(2*m - 1e-30),
        my(b = 2*a - c*unitvec(3, 1));
        if (type(sqrtK(F, b)) == "t_COL", has = 1; break);
        listput(cs, b));
      if (!has, forced = 1;
        for (j = 1, #cs, listput(Ks, polredbest(rnfequation(F[1], 'x^2 - alg(F, cs[j])))));
        break));
    if (!forced, write(outfile, Str(D, "|", P, "|nu=", nu, "|H(F) holds: K not determined")); next);
    \\ test the forced sextic fields
    my(res = List());
    for (j = 1, #Ks,
      my(Q = Ks[j]);
      if (polsturm(Q) != 6, listput(res, "notTR"); next);
      my(K = gfield(Q), w = gwitness2(K, 300, 3));
      listput(res, if (w === 0, Str("SURVIVOR ", Q), Str(nfdisc(Q), ":", w[1], ":N=", w[3]))));
    write(outfile, Str(D, "|", P, "|nu=", nu, "|forced ", #Ks, "|", Vec(res))));
}
