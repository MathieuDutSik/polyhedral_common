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
* **Problem2 is not solved by that method.** Its rows are of the form
  "a sum of about ten small coefficients equals 10", so the bound
  propagation has too much slack to prune and the search tree does not
  close. It is the harder instance for a combinatorial search even
  though the reference solver, which is lattice based, does it in
  3.94 s.

What was measured on Problem2, for whoever picks this up:

* Bound propagation alone reaches 300000 nodes without a solution.
* The exact linear programming relaxation is cheap at the root
  (0.7 s for the 778 x 277 program) but its bounds are loose: it only
  gives 44 <= sum_j x_j <= 64, where the solutions have 50 ones.
* Adding the 112 rows of an LLL reduced basis of the row lattice of
  `[A|b]` as extra constraints (coefficients then at most 7) does not
  close the tree either.
* The lattice formulation works and is the promising direction, but
  the enumeration is the obstacle. With `x0` an integral solution and
  `K` a basis of the kernel, `x` is 0/1 exactly when
  `|| 2x - 1 ||^2 <= n`, so the solutions are the lattice points of a
  ball and the count of points equals the count of solutions. The
  kernel has dimension 164, LLL brings its Gram from 10^37 down to
  norms 10 to 74 in 2.5 s, and Babai round-off brings the centre from
  a squared norm of 10^711 down to 2420. But the exact Fincke-Pohst of
  `CVPSolver` does not finish in dimension 164 even at a fifth of the
  radius. Closing this needs the stronger pruning of a dedicated
  solvediophant style enumeration.

The CI test
-----------

    cd CI_tests/08B_AllZeroOneSolutions
    ../../gap.sh < TestZeroOneSolutions.g

driven by `.github/workflows/ci_08B_zero_one_solutions.yml`, on day 8
of the even months, `ci_08A_cone_int` holding the odd ones. It runs
seven small cases under the three arithmetics `gmp`, `gmp_boost` and
`multi_boost`, then the node budget, then Problem1. The whole test is
a little over a minute, almost all of it Problem1.

The small cases are checked against a brute force enumeration of the
2^n vectors done in GAP, so their expected answer does not come from
the program under test. For every case, small or not, each returned
vector is checked to be 0/1 and to satisfy `A x = b`, and the set of
returned solutions is compared with the expected one. Problem1 is
compared against `Problem1.solutions`.

Problem2 is not part of the test, for the reason below.
