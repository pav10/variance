\\ spectrum.gp -- irreducible totally real algebraic integers of small variance.
\\ Reads results/spectrum/raw_dD.txt (lines "a1 ... aD") and keeps irreducible f.
varpol(f) = my(d = poldegree(f), a1 = polcoef(f, d-1), a2 = polcoef(f, d-2)); ((d-1)*a1^2/d - 2*a2)/d;
spec(D, V0, out) =
{
  my(L = readstr(Str("results/spectrum/raw_d", D, ".txt")), R = List());
  for (i = 1, #L,
    my(c = eval(Str("[", strjoin(strsplit(L[i], " "), ","), "]")), f);
    if (#c != D, next);
    f = Pol(concat([1], c), 'x);
    if (!polisirreducible(f), next);
    my(v = varpol(f));
    if (v <= V0, listput(R, [v, f, nfdisc(f), polgalois(f)[1]])));
  R = vecsort(Vec(R), 1);
  for (i = 1, #R, write(out, Str(D, "|", R[i][1], "|", R[i][1]*1., "|", R[i][2], "|disc=", R[i][3])));
  #R;
}
