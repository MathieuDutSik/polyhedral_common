# src_non_lattice_packing — status and issues

Double-precision search-and-certify tools for m-periodic sphere packings,
after Andreanov–Kallus (DCG 2019) and Schuermann's periodic-extremeness
theory. The exposition of the theory is in
`Paper_non_lattice/NonLatticePacking.tex`; this file records what the
implementation does, what is validated, and what is missing.

## What works (validated)

* `evaluate`: packing density `phi = V_n (lambda/2)^n m / sqrt(det Q)` with
  `lambda^2` from an exact bounded enumeration of the pair vectors
  `v + c_t - c_s`. Reproduces to machine precision:
  - fcc (m=1): `phi = 0.7404804896931 = pi/sqrt(18)`, 12 minimal vectors;
  - bcc as `Z^3 + {0,(1/2,1/2,1/2)}` (m=2): `phi = 0.6801747615878`;
  - `D_4`: `pi^2/16`, 24 minimal vectors; `D_5`: `pi^2 sqrt(2)/30`, 40.
* `certify`: the first-order extremeness test in the `(Q,C)` chart — the LP
  asking for a direction with `(n/2)(min_j g_j . d)/lambda^2 > (1/2) gdet . d`
  over the active pair-norm gradients. fcc, bcc, `D_4`, `D_5` all certify
  `extreme=1` with full-rank active families and `floating=0`, matching their
  known perfection/eutaxy. The rate is recomputed from the returned direction
  (the LP is fully degenerate, every active row has zero right-hand side, and
  the simplex objective cannot be trusted — the same failure mode diagnosed
  and fixed in the covering certificate).
* `descend`: staged-softmin L-BFGS maximization of `log phi` on a frozen
  short-vector list, re-enumerated per round, scale pinned between stages.
  Correctly HOLDS local optima (bcc stays at 0.6801747: it is extreme, so no
  descent path to fcc exists — that is the theory, not a failure).

## Bugs found and fixed during construction (recorded because instructive)

1. **Coset collision = free density.** Excluding all near-zero vectors from
   the minimum let two cosets collapse onto each other, doubling `m` at no
   cost in `lambda` (`phi = 1.047 > 0.74` in d=3, an impossibility). Only the
   trivial same-coset zero vector may be excluded; a near-zero cross-coset
   vector is a genuine collision whose tiny norm must dominate the minimum.
2. **Unreduced cosets defeat a zero-centered enumeration box.** Nothing
   constrains the cosets to a fundamental domain, and once a descent drifts
   them far from the origin every cross-coset vector lies outside the box:
   `lambda` silently becomes the pure-lattice minimum and the collapse cheat
   reopens. The enumeration now reduces each coset difference into
   `[-1/2,1/2)^n` (exact — the integer shift is folded into the stored
   integer part), and the descent wraps `C` into `[0,1)^n` each round.

## Known issues / future work

* **Global search is weak from random seeds.** The softmin objective is
  driven only by the few near-minimal vectors, and plain multistart lands at
  poor local optima (`phi ~ 0.3` in `n=3, m=2` after 5 starts). The remedy,
  as for the covering, is basin hopping with cheap kicks on top of this
  descent; not yet wired for packing.
* **The fluid/floating diagnostic is untested on a true fluid family.** The
  first genuine fluid packing (Andreanov–Kallus) is the 9-dimensional fluid
  diamond `D_9^+`; no low-dimensional test case exists because fluidity
  needs `mu > 2 rho` at an extreme lattice, first satisfied at `D_9`.
  Validating `floating > 0` on `D_9^+(b)` is a natural test to add.
* **No exact vertex enumeration.** The generalized Voronoi algorithm
  (enumerate the vertices of the Ryshkov-like polyhedron `R(lambda)` up to
  the symmetry group, solve the rank condition on faces, test algebraic
  extremeness exactly) is what produces Andreanov–Kallus's complete d<=5,
  m=2 classification. It needs exact dual-description machinery over the
  `(n+m-1)`-dimensional form space — the natural next implementation, using
  the existing polyhedral_common exact infrastructure. For m>2, additionally,
  Theorem 2.6 (local optimum => algebraically extreme or fluid) is only
  conjectural, and the search must walk the `m(m+1)/2`-dimensional faces of
  `R(lambda)` rather than its edges.
* Double precision only; a candidate record found here must be re-verified
  exactly (rational or algebraic reconstruction), as for the coverings.
