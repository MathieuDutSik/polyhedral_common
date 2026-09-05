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

Still to do
-----------

The GAP driver, the access point and the workflow `ci_08B_...` come
with the enumeration program. Day 8 is free for an even month schedule,
`ci_08A_cone_int` already running in odd months.
