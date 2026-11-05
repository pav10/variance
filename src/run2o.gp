\\ run2o.gp -- G-witnesses for every field that the 2O pipeline (lifting problem) left
\\ for its maximal-order stage in degrees 6 and 8.  Input files: data/2o/*.
read("src/gtest.gp");
g2o(file, out, Nmax, Tmax, part = 0, nparts = 1) =
{
  my(L = readvec(file), seen = Map(), k = 0);
  for (i = 1, #L,
    my(P = L[i][1], D, K, w, t0 = getabstime(), lin);
    D = nfdisc(P);
    if (mapisdefined(seen, D), next); mapput(seen, D, 1);
    k++; if ((k-1) % nparts != part, next);
    K = gfield(polredbest(P));
    w = gwitness2(K, Nmax, Tmax);
    lin = Str(D, "|", K[1].pol, "|2unram=", (D % 2 != 0), "|",
      if (w === 0, "SURVIVOR", Str(w[1], "|", minpol(K, w[2]), "|N=", w[3])), "|", getabstime() - t0, "ms");
    write(out, lin));
}
