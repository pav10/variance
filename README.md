# variance — sums of squares, the class G, and the trace variance

Work on the class

    G = { K totally real : 2 O_K^+ ⊆ Σ O_K^2 },

whose members of degree ≤ 5 are known (Kala–Yatsyna: Q, Q(√2), Q(√3), Q(√5), Q(√2,√5),
Q(ζ20)+, Q(ζ7)+, Q(ρ) with ρ³ = 4ρ + 2), using the trace variance Var(x) = τ(x²) − τ(x)².
The write-up is [`notes/G_beyond5.pdf`](notes/G_beyond5.pdf) (source `notes/G_beyond5.tex`).

## Main results (see the notes for precise statements)

* **Exact witnesses.** `src/gtest.gp` decides `2α ∈ ΣO_K²` exactly and searches for exceptional
  elements (unit criterion, small trace, small norm); every witness is certified without floating
  point or class-group input.
* **Census.** All LMFDB totally real fields of degree 3–8 in the lists used by the `sails` project:
  in degrees ≤ 5 exactly the Kala–Yatsyna members survive; in degrees 6, 7, 8 every field has a
  certified witness.
* **2 unramified, degrees 6–8.** For fields with 2 unramified, membership in G implies all the
  conditions used in the degree 6–8 lifting-problem classification of `pav10/2O`; with explicit
  G-witnesses for every field that classification leaves, no field of degree 6, 7, 8 with 2
  unramified lies in G.
* **Bottom of the variance spectrum.** All totally real algebraic integers with Var ≤ 2 in degrees
  ≤ 10, and a Smyth-type auxiliary-function certificate (interval arithmetic) showing Var(x) ≥ 1.69195
  for every other totally real algebraic integer. Hence the values below 1.69195 are exactly
  0, 5/4, 14/9, 27/16, and ν(K) < 1.69195 iff K contains Q(√5) (ν = 5/4) or Q(ζ7)+ (ν = 14/9).
* Shapes/statistics, witnesses for residual real cyclotomic fields (F₃₆, F₄₄, ...), and a reduction
  of the 2-ramified sextic problem.

## Layout

```
src/gtest.gp        toolkit (decision procedure, witnesses, invariants)   [gp -D parisizemax=2G]
src/census.gp       census driver; src/runcensus.sh DEG NMAX [NSH]  (TMAX env = trace bound)
src/finish_survivors.gp  rerun survivors with a larger norm bound
src/run2o.gp        G-witnesses for the fields of pav10/2O (data/2o/)
src/regime1.gp      ramified sextic Case A, regime nu_2 >= 2 nu
src/quadB.gp        quadratic subfields in Case B
src/cyclo.gp        real cyclotomic fields
src/enum6.c         Robinson enumerator (copied from pav10/2O), src/varspec.sh, src/spectrum.gp
src/lp/             auxiliary-function LP (varlp.py) and its interval-arithmetic certificate
src/gen_tables.py   LaTeX tables of the notes from results/
data/fields/        LMFDB totally real field lists (rvisser7/indecomposables mirror)
data/2o/            field lists produced by the pav10/2O pipeline
results/            census, witnesses, spectrum, LP certificates
notes/              write-up
```

Requirements: PARI/GP ≥ 2.15 (with galdata), gcc, python3 with numpy, scipy, mpmath; pdflatex for the notes.
