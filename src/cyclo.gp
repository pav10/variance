\\ cyclo.gp -- G-witnesses for real cyclotomic fields F_m = Q(zeta_m)^+.
read("src/gtest.gp");
\\ minimal polynomial of 2cos(2 pi/m)
realcyc(m) = my(R = charpoly(Mod('z + 'z^(m-1), polcyclo(m, 'z)), 'x)); factor(R)[1,1];
cyctest(m, Nmax, Tmax) =
{
  my(P = realcyc(m), Q = polredbest(P), t0 = getabstime(), K, w, res);
  K = gfield(Q);
  w = gwitness2(K, Nmax, Tmax);
  res = Str("m=", m, "|d=", poldegree(P), "|hplus/h=", 2^#tpunits(K), "|",
    if (w === 0, Str("SURVIVOR(N<=", Nmax, ",tau<=", Tmax, ")"),
      Str(w[1], "|N=", w[3], "|Tr=", tr(K, w[2]), "|minpoly=", minpol(K, w[2]))),
    "|", getabstime() - t0, "ms");
  print(res);
  res;
}
