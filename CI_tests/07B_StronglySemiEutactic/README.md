Strongly semi-eutactic forms
============================

This section exercises `src_perfect/IsStronglySemiEutactic.cpp`.

A positive definite quadratic form `A` is *strongly semi-eutactic* when

    A^{-1} = mu * sum_{v in S} v v^T

for some `mu > 0` and some subset `S` of `Min_h(A)`, the set of
shortest vectors of `A` with one representative kept per antipodal
pair `{v,-v}`. The strongly eutactic forms are the special case
`S = Min_h(A)`. The program answers the question by a branch and bound
search in exact arithmetic; a positive answer comes with the
certificate `(mu, S)`, a negative answer is the exhaustion of the
search tree, and the search can also stop undecided when the node
budget `max_node` is reached.

Running the test
----------------

    cd CI_tests/07B_StronglySemiEutactic
    ../../gap.sh < TestStronglySemiEutactic.g

The binary is built by `make -C src_perfect -f Makefile_eutactic`.
The whole test takes about half a minute.

`ListCases.g` holds the cases as a list of records and
`TestStronglySemiEutactic.g` runs them through the access point
`test_strongly_semi_eutactic` of `../access_points.g`. For each case
the reported minimum, number of shortest vectors, resolved/unresolved
status and yes/no answer are compared with the expected ones, and for
a positive answer also `mu` and `|S|`. The printed subset is checked
to have `|S|` distinct entries within range, and the certificate is
checked against the identity

    n = <A^{-1}, A> = mu * sum_{v in S} v^T A v = mu * |S| * min

obtained by pairing the eutaxy relation with `A` itself.

The cases
---------

* **The 34 strongly semi-eutactic perfect forms of dimension 9**
  (`Perfect9_01` ... `Perfect9_34`). This is the complete list among
  the 2237251040 perfect forms of dimension 9, obtained by an
  exhaustive scan of the netCDF database of the paper "The complete
  classification of perfect forms in dimension 9". Among them are
  `A_9`, `A_9^2`, `A_9^5`, `D_9`, `Q_99`, `Q_129` and `Lambda_9`.
  Ten of them are strongly eutactic (`S = Min_h`), the other 24 have
  a proper subset `S`, which is what makes them the interesting part
  of the test.
  The expected `mu` and `|S|` are those of the scan, and each of the
  34 certificates was re-verified exactly (the listed vectors are of
  minimal norm and `A^{-1} = mu * sum_{v in S} v v^T` holds) before
  being recorded here. The program reproduces all 34 of them.

* **Classical strongly eutactic lattices**: `A_3`, `D_4`, `Z^4`,
  `A_4^*`, `E_6^*`, `E_7^*`. For those the minimal vectors were
  enumerated independently and `A^{-1} = mu * sum_{v in Min_h} v v^T`
  was checked, so the expected `mu` and `|S| = |Min_h|` do not come
  from the program.

* **`E_8 + A_1`**: strongly semi-eutactic with `|S| = 9` out of the
  121 antipodal pairs, that is with a subset far smaller than
  `Min_h`. Here the expected `mu` and `|S|` are recorded from the
  program, the run-time checks above being what guards them.

* **Forms that are not strongly semi-eutactic**. The orthogonal sums
  `A_2 + <2>`, `A_3 + <2>` and `E_6 + <2>` are refuted by an argument
  that does not need the program: `A^{-1}` is block diagonal and its
  last entry `1/2` is nonzero, so the last basis vector belongs to
  `S` and `mu = 1/2`; the first block then requires
  `sum_{v in S} v v^T = 2 * B^{-1}` for `B` the block, which is not
  integral since `det(A_2) = 3`, `det(A_3) = 4` and `det(E_6) = 3`.
  Likewise for `E_7 + A_2`: the `A_2` block admits `mu = 1/3` only,
  and `3 * E_7^{-1}` is not integral since `det(E_7) = 2`. That last
  case is also the most substantial refutation of the test, 843 nodes
  of search. `E_6 + A_3` and `A_8 + A_1` are refuted too.

* **`E7_plus_A2_budget100`**: the same `E_7 + A_2`, run with
  `max_node = 100`, which is less than the 843 nodes the search
  needs. This exercises the unresolved outcome of the program.

Arithmetics
-----------

Every case is run with the `gmp` arithmetic. The cases carrying
`sweep:=true` in `ListCases.g` are run with `gmp_boost` and
`multi_boost` as well, and must give the same answer.
