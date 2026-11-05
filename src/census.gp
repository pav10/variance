\\ census.gp -- witness census over a list of totally real fields.
\\ Input file: lines  "index|c0,c1,...,c_{n-1}|mono"  (rvisser7 / LMFDB format, coefficients
\\ ascending without the leading 1).  Output: one line per field
\\   disc|coeffs|d|hplus/h|dyadic type|nu|status|witness minpoly|norm|time_ms
\\ status: U = unit witness, N = small-norm witness, S = survivor (no witness up to Nmax).
\\ Run:  gp -q -D parisizemax=2G census.gp  after defining  INFILE, OUTFILE, SHARD, NSH, NMAX.

read("src/gtest.gp");
TMAX = 0;   \\ trace stage bound (0 = off); set before calling census()
TLIM = 120; \\ per-field time limit in seconds (status E = error or time-out)

dyadtype(K) =
{
  my(P = idealprimedec(K[1], 2));
  strjoin(vector(#P, i, Str(P[i].e, "^", P[i].f)), ",");
}

censusone(P, Nmax) =
{
  my(t0 = getabstime(), K, w, st, mp = "", nm = "", wt = "", nu, nu2, hp, D, d, h, R, shp);
  K = gfield(P);
  D = nfdisc(P); d = poldegree(P);
  hp = 2^#tpunits(K);
  nu = nuK(K)[1]; nu2 = nu2K(K);
  h = K[2].no; R = K[2].reg;
  \\ normalised shape minimum: nu divided by det(Var|Lambda)^(1/(d-1)) = (D/d^d)^(1/(d-1))
  shp = nu / (D / d^d)^(1/(d-1));
  w = gwitness2(K, Nmax, TMAX);
  if (w === 0, st = "S",
    st = if (w[1] == "unit", "U", if (w[1] == "trace", "T", "N"));
    mp = minpol(K, w[2]); nm = w[3]; wt = tr(K, w[2]));
  [D, d, hp, dyadtype(K), nu, st, mp, nm, getabstime() - t0, nu2, h, R, shp, wt];
}

census(infile, outfile, shard, nsh, Nmax) =
{
  my(L = readstr(infile), k = 0);
  for (i = 2, #L,
    my(s = L[i], parts, co, P, r);
    if (#s == 0, next);
    parts = strsplit(s, "|");
    co = eval(Str("[", parts[2], "]"));
    k++;
    if ((k - 1) % nsh != shard, next);
    P = Pol(Vecrev(concat(co, [1])), 'x);
    r = iferr(alarm(TLIM, censusone(P, Nmax)), E, [nfdisc(P), poldegree(P), "", "", "", "E", Str(errname(E)), "", 0]);
    if (#r < 14, r = concat(r, vector(14 - #r, k, "")));
    write(outfile, Str(r[1], "|", Vecrev(Vec(P))[1..#co], "|", r[2], "|", r[3], "|", r[4], "|",
      r[5], "|", r[6], "|", r[7], "|", r[8], "|", r[9], "|", r[10], "|", r[11], "|",
      if (r[12] === "", "", strprintf("%.6f", r[12])), "|", if (r[13] === "", "", strprintf("%.6f", r[13])), "|", r[14])));
}
