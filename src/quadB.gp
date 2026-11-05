\\ quadB.gp -- quadratic fields F that can contain the minimal vector of a sextic K in G
\\ (Case B, [K:F] = 3 odd, so squares of K lying in F are squares in F).  Necessary condition:
\\ every alpha in O_F^+ \ Z with tau(alpha) < nu(F) has 2 alpha - c in O_F^2 for some integer
\\ 0 <= c < 2 min sigma(alpha)   (Corollary dich(2) with nu(K) = nu(F)).
read("src/gtest.gp");
quadok(D) =
{
  my(K = gfield(x^2 - D, 0), nu = nuK(K)[1], L = tpos_small(K, nu - 1/1000));
  for (i = 1, #L,
    my(a = L[i], ta = tau(K, a), ok = 0, m = vecmin(emb(K, a)));
    if (ta >= nu, next);
    for (c = 0, floor(2*m - 1e-30), if (type(sqrtK(K, 2*a - c*unitvec(2,1))) == "t_COL", ok = 1; break));
    if (!ok, return([0, alg(K, a), ta, nu])));
  [1, nu];
}
{forsquarefree(D = 2, 2000, my(r = quadok(D[1])); if (r[1], print("D=", D[1], " passes, nu=", r[2])));}
