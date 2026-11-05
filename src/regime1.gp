\\ regime1.gp -- sextic Case A, regime nu_2 >= 2 nu: the candidates that the lifting-problem
\\ pipeline excludes only by rigidity (not available for G when 2 ramifies) are tested at the
\\ level of the field: distinct fields, G-witness search.
read("src/gtest.gp");
reg1(infile, outfile, part, nparts, Nmax, Tmax) =
{
  my(L = readstr(infile), seen = Map(), k = 0);
  for (i = 1, #L,
    my(f = eval(L[i]), g = polredabs(f), D);
    if (mapisdefined(seen, g), next); mapput(seen, g, 1);
    k++; if ((k - 1) % nparts != part, next);
    D = nfdisc(g);
    my(K = gfield(g), t0 = getabstime(), w = iferr(alarm(600, gwitness2(K, Nmax, Tmax)), E, Str("ERR ", errname(E))));
    write(outfile, Str(D, "|", g, "|2unram=", D % 2 != 0, "|",
      if (type(w) == "t_STR", w, if (w === 0, "SURVIVOR", Str(w[1], "|", minpol(K, w[2]), "|N=", w[3]))),
      "|", getabstime() - t0, "ms")));
}
