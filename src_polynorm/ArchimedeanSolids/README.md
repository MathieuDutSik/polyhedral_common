The Platonic and Archimedean solids
===================================

Exact vertex coordinates of the 5 Platonic and 13 Archimedean solids, in
the input format of `POLYNORM_BetkeHenk` (one vertex `(1, v)` per row):

  * rational: tetrahedron, cube, octahedron, truncated tetrahedron,
    cuboctahedron, truncated octahedron (`arith = gmp`);
  * over `Q(sqrt 2)` with `x = sqrt 2`: truncated cube, rhombicuboctahedron,
    truncated cuboctahedron (`arith = Qsqrt2`);
  * over `Q(sqrt 5)` with `x = sqrt 5`: dodecahedron, icosahedron,
    icosidodecahedron, truncated dodecahedron, truncated icosahedron,
    rhombicosidodecahedron, truncated icosidodecahedron (`arith = Qsqrt5`);
  * the snub cube over `Q(t)`, `t^3 = t^2 + t + 1` (the tribonacci constant),
    and the snub dodecahedron over `Q(xi)`, `xi^3 - 2 xi = phi` (degree 6),
    with the field descriptions `snub_cube.field` and
    `snub_dodecahedron.field` (`arith = RealAlgebraic=<file>.field`).

The files are produced by `python3 generate_solids.py <outdir>` (sympy and
mpmath), which checks that every edge has the same length and every vertex
the expected degree.

The densest lattice packings, as computed by Betke and Henk
(arXiv:math/9909172), are reproduced by for example

    POLYNORM_BetkeHenk Qsqrt5 icosahedron.ext stdout 8

and the table is in `../NOTES.md`.
