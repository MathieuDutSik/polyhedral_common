All 0/1 solutions of A x = b
============================

Test data for the enumeration of every `x` in `{0,1}^n` satisfying

    A x = b

with `A` a nonnegative integer matrix and `b` an integer vector. The
enumeration program is not written yet; this directory holds the two
reference instances, converted to the polyhedral_common format, and
the complete list of their solutions.

Files
-----

For each instance `k` in `{1, 2}`:

* `Problem<k>.matrix` -- the matrix `A`, in the matrix format read by
  `ReadMatrixFile`: a first line `nbRow nbCol` followed by the rows.
* `Problem<k>.rhs` -- the right hand side `b`, in the vector format read
  by `ReadVectorFile`: a first line `nbRow` followed by the entries.
* `Problem<k>.solutions` -- the complete set of solutions, as a matrix with
  one solution per row, so `n_sol x n`. The rows are sorted
  lexicographically, so that the file can be compared as it is with a
  sorted output of the program.

|            | `m` (rows) | `n` (columns) | solutions | entries of `A` | reference time |
|------------|------------|---------------|-----------|----------------|----------------|
| Problem1   | 13         | 406           | 56        | 0 to 360       | 6 min 38 s     |
| Problem2   | 113        | 276           | 6         | 0 to 48        | 3.94 s         |

The right hand side is `b = (30, ..., 30, 1140)` for the first
instance and `b = (10, ..., 10, 2080)` for the second. The solutions of
Problem1 have 10 to 13 nonzero coordinates, those of Problem2 have
exactly 50.

Provenance
----------

The instances come from `problem1.in` / `problem2.in` of an earlier
solver, whose own format was a header line `m n 1`, a blank line and
then `m` lines holding the `n` coefficients of a row followed by its
right hand side; the solutions were the `.out` files, one line of `n`
characters `0` or `1` per solution. The conversion only splits `A`
from `b` and rewrites the solutions as a matrix; no value was changed.

Every solution recorded here was checked to be 0/1 and to satisfy
`A x = b` exactly, reading the files back through `ReadMatrixFile` and
`ReadVectorFile`. The reference times above are those of the earlier
solver (`problem1.time`, `problem2.time`), kept as an order of
magnitude for what the instances cost; they are not a target.

The solution lists are believed to be complete, that being the point of
the instances, but completeness was not re-established independently
here -- it is what the new code is meant to confirm.

Status of the enumeration program
---------------------------------

`src_milp/LATT_ZeroOneSolutions`, built by `make -C src_milp`,
enumerates the solutions with
the branch and bound of `src_milp/zero_one_solution.h`:

    cd CI_tests/08B_AllZeroOneSolutions
    ../../src_milp/LATT_ZeroOneSolutions gmp Problem1.matrix Problem1.rhs ZeroOne out

* **Problem1 is solved**: the 56 solutions, complete, in 68 s and
  5328493 nodes, the output being the same set as `Problem1.solutions`.
  The reference solver took 6 min 38 s on it.
* **Problem2 is solved by the lattice method**, in 4 min 57 s and
  2284675431 nodes, the 6 solutions being the recorded ones. The
  branch and bound of `zero_one_solution.h` does not close its tree,
  its rows being sums of about ten small coefficients equal to 10, so
  the bound propagation has too much slack. Problem1 too is solved by
  the lattice, in 59 s and 181278175 nodes.

What it took, and what did not work
-----------------------------------

The lattice is `ker_Z([A | -b])` with `Q(z,s) = ||2z - s 1||^2 + s^2`,
built by `MILP_ZeroOneLattice`; the 0/1 solutions are exactly its
vectors of norm n+1 with s = +-1. Homogenizing is what removes the
coset, so there is no Hermite normal form of a particular solution
with 700 digits and no rational centre with a huge denominator.

What did not work was enumerating that lattice by norm alone. The
solutions have norm 277 where the lattice minimum is about 40, so it
is not a short vector search but a search through a dense region seven
times the minimum: the widest level of a plain Fincke-Pohst holds
about 10^23 nodes, and about 10^17 even for a perfectly flat profile
of the same determinant, which is the best any reduction could give.
BKZ was therefore not worth writing, and `... profile` reports those
numbers for any instance.

What works is not enumerating a ball at all. A solution has *every*
ambient coordinate equal to +-1, far stronger than having norm n+1,
and `zero_one_enum.h` carries that into the tree, after the enumeration
of A. Wassermann's solvediophant:

* Determined coordinates. With `first_nonzero[l]` the least j with
  B[j][l] nonzero, every vector of the span of b_0, ..., b_{i-1} has
  its l-th coordinate zero for i = first_nonzero[l], so the projection
  P_i v already has the final value of the coordinate l: it is settled
  at level i and must be +-1 there. Half of the 277 coordinates are
  settled by level 50 of 165, and only 10 wait until the last level.
* Hoelder. P_i being an orthogonal projection,
  cs_i = ||P_i v||^2 = <P_i v, v> <= ||v||_infinity ||P_i v||_1, so
  cs_i <= ||w_i||_1 holds along every solution. On Problem2 this fires
  1663318138 times against 58828807 for the coordinates, so it carries
  most of the load.
* Dual bounds, |t_i| <= min(||d_i||_1, sqrt(Fd) ||d_i||_2). Cheap, but
  measured on these instances they alone leave 10^224 and 10^354
  candidates, so they are a complement and not a lever.

The ball of radius sqrt(n+1) exceeds the cube [-1,1]^{n+1} by about
10^86 in dimension 277, which is the whole distance between the 10^17
above and what the search actually costs.

The Gram-Schmidt and the tests are in double precision, as in
solvediophant, with a tolerance eps; every vector reaching the bottom
is verified exactly in integers, so a reported solution is certain
while completeness holds up to that tolerance. Beware that the
solutions have squared norm exactly equal to the bound, so the
accumulated norm lands on it and rounding pushes it barely over: the
early return on a negative remaining radius needs the same tolerance
as the pruning test, or solutions are silently lost. That was the one
real bug, and it cost 8 of the 56 solutions of a random test case.

The CI test
-----------

    cd CI_tests/08B_AllZeroOneSolutions
    ../../gap.sh < TestZeroOneSolutions.g

driven by `.github/workflows/ci_08B_zero_one_solutions.yml`, on day 8
of the even months, `ci_08A_cone_int` holding the odd ones. It runs

* seven small cases under the three arithmetics of the branch and
  bound, and the same seven through the lattice enumeration, which is
  an independent implementation;
* the node budget, which has to stop the branch and bound on Problem1;
* Problem1 by the branch and bound, 5328493 nodes;
* Problem1 by the lattice, 181278175 nodes;
* Problem2 by the lattice, 2284675431 nodes. It is out of reach of the
  branch and bound.

About seven minutes, five of them Problem2.

The small cases are checked against a brute force enumeration of the
2^n vectors done in GAP, so their expected answer comes neither from
the branch and bound nor from the lattice. For every case, each
returned vector is checked to be 0/1 and to satisfy `A x = b`, and the
set of returned solutions is compared with the expected one. Problem1
and Problem2 are compared against their recorded solutions, and
Problem1 being solved by both methods is a cross-check between two
independent implementations.

The node counts above are not asserted, only the solution sets are:
the lattice enumeration works in double precision, so its node counts
depend on the machine while its answers do not.

