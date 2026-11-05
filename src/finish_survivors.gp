\\ finish_survivors.gp -- rerun the census survivors of a degree with a larger norm bound
read("src/gtest.gp");
finish(infile, outfile, Nmax, Tmax) =
{
  my(L = readstr(infile));
  for (i = 1, #L,
    my(p = strsplit(L[i], "|"));
    if (p[7] != "S", next);
    my(co = eval(p[2]), P = Pol(Vecrev(concat(co, [1])), 'x), K = gfield(P), w = gwitness2(K, Nmax, Tmax));
    write(outfile, Str(p[1], "|", P, "|", if (w === 0, Str("SURVIVOR(N<=", Nmax, ")"),
      Str(w[1], "|", minpol(K, w[2]), "|N=", w[3], "|Tr=", tr(K, w[2]))))));
}
