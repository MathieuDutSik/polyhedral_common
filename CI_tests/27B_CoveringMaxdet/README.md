# Covering density optimization over iso-Delaunay domains

This test exercises the determinant maximization ("MAXDET") formulation of the
sphere covering problem implemented in `src_delaunay/CoveringMaxdet.h` and
reached through the `FileCoveringOptimum` query of
`LATT_AnalysisIsoDelaunay`.

Over one iso-Delaunay domain the covering optimization is the convex problem

```
minimize   -log det Q          subject to   Q in the L-type cone,  mu(Q) <= 1
```

where `mu(Q) <= 1` is a linear matrix inequality (Proposition 6.1 of the
reference below): one `(d+1) x (d+1)` circumradius block per orbit of Delone
polytopes, every entry linear in `Q`. The header solves it with a barrier
path-following method that depends on nothing outside this repository.

## What is checked

**Lattice domains.** The full pipeline — enumerate the iso-Delaunay domains of
the classic `GL_d(Z)` T-space with `LATT_SerialLattice_IsoDelaunayDomain`
(dumping each of them via `PrefixIsoDelaunayDomains`), then optimize the
covering density over each — has to reproduce the classical optima:

| dimension | domains | best covering density | attained by |
| --------- | ------- | --------------------- | ----------- |
| 3         | 1       | 1.4635030689668180 = 5 pi sqrt(5) / 24 | `A_3^*` |
| 4         | 3       | 1.7655285081493524    | `A_4^*` |
| 5         | 222     | 2.1242859089916246    | `A_5^*` |

At every optimum the squared covering radius has to be 1, which is the
normalization the LMI imposes and therefore a check that the constraint is
active where it should be.

**Periodic domains.** The same optimization is run over the iso-Delaunay
domains of the periodic point set `Z^3 + {0, (1/3,1/3,1/3)}`, enumerated by
`LATT_SerialPeriodic_IsoDelaunayDomain`. The optima there are not published
constants, so what is checked is what is known independently:

* the squared covering radius is again 1 at every optimum;
* the point density factor is `2 / 3^3` — a periodic set of `m` cosets with
  common denominator `N` needs the lattice density formula multiplied by
  `m / N^d`, and getting that factor wrong is the easiest way to make the
  whole periodic computation silently meaningless;
* no domain beats `A_3^*`, the best lattice covering of dimension 3 and
  conjecturally the best covering of dimension 3 altogether.

## Reference

Mathieu Dutour Sikirić, Achill Schürmann, Frank Vallentin, *A generalization of
Voronoi's reduction theory and its application*, Duke Mathematical Journal 142
(2008) 127--164, arXiv:math/0601084 — Section 6.2 for the determinant
maximization formulation, and the 5-dimensional iso-Delaunay classification
(222 domains).

Note: the 5-dimensional enumeration alone takes a few minutes, so this is a
long-running test, appropriate for the bi-monthly CI schedule.
