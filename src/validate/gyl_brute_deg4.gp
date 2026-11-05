rd(file) = my(L = readstr(file), R = List()); for (i = 1, #L, my(c = eval(Str("[", strjoin(strsplit(L[i], " "), ","), "]"))); if (#c == 4, listput(R, Pol(concat([1], c), 'x)))); Vec(R);
expected(f) =
{
  my(d = 4, th = vecsort(real(polroots(f))), a1 = polcoef(f, d-1), a2 = polcoef(f, d-2), ta = -a1/d, nu = ((d-1)*a1^2/d - 2*a2)/d, n, T, Q = subst(f, 'x, 'w));
  n = ceil(th[1]) - 1; T = ta - n;
  if ((floor(th[d]) + 1) - ta < nu/2 - 1e-9, return(0));
  my(r1 = (sqrt(2*d) + sqrt(2*d + 12))/3, NUONE = r1^2);
  if (nu > NUONE && (floor(th[d]) + 1) - ta >= nu - 1e-9, return(0));
  for (c = 0, 1,
    if (T < (nu + c)/2 || T >= nu, next);
    my(r = nfroots(Q, 'x^2 - (2*('w - n) - c)));
    for (k = 1, #r,
      my(y = lift(r[k]), cp = charpoly(Mod(y, Q), 'x), yt = vecsort(real(polroots(cp))), ty = -polcoef(cp, d-1)/d, vy = vecsum(apply(z -> (z - ty)^2, yt))/d);
      if (ty < 0, next);
      if (vy < nu - 1e-9, next);
      if (ty - (ceil(yt[1]) - 1) < nu/2 - 1e-9 || (floor(yt[d]) + 1) - ty < nu/2 - 1e-9, next);
      return(1)));
  0;
}
A = rd("/tmp/claude-0/gval/all.txt"); G = rd("/tmp/claude-0/gval/gyl.txt");
E = Set(select(f -> polisirreducible(f) && expected(f), A));
GG = Set(select(f -> polisirreducible(f), G));
print("expected ", #E, "  GYL ", #GG, "  only-expected ", #setminus(E, GG), "  only-GYL ", #setminus(GG, E));
print(setminus(E, GG)); print(setminus(GG, E));
