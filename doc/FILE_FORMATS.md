Input and Output File Formats
=============================

Almost every command-line program in this project exchanges data as plain
text files holding **matrices**. This document specifies those formats exactly,
so that inputs can be written by hand or generated from another system, and
outputs can be parsed reliably.

The reference implementation of the readers and writers is
`basic_common_cpp/src_matrix/MAT_MatrixFund.h` (`ReadMatrix`, `WriteMatrix`,
`ReadListMatrix`, `WriteMatrixGAP`, ...). If this document and that code ever
disagree, the code is authoritative.


The native matrix format
------------------------

This is the format read by `ReadMatrix` / `ReadMatrixFile` and written by
`WriteMatrix`. It is what every program means by "a matrix file".

Layout:

```
<nbRow> <nbCol>
<entry(0,0)> <entry(0,1)> ... <entry(0,nbCol-1)>
<entry(1,0)> <entry(1,1)> ... <entry(1,nbCol-1)>
...
```

* The **first two numbers** are the row count and the column count, both
  non-negative integers. Negative values are rejected.
* The **remaining `nbRow * nbCol` numbers** are the matrix entries in
  **row-major order** (all of row 0, then all of row 1, ...).

Parsing is done purely with the C++ stream operator `>>`, so the format is
**entirely whitespace-insensitive**: spaces, tabs and newlines are
interchangeable as separators. The line structure shown above is a convention
for readability only — the entire matrix may equally be written on a single
line, or one entry per line. The reader stops once it has consumed
`2 + nbRow * nbCol` tokens; anything after that in the file is ignored.

### Entry syntax depends on the arithmetic

Each entry is parsed by the `operator>>` of the number type selected on the
command line (the `arith` argument of most programs). In practice:

* **Integers** (`mpz_class`, `int64_t`, ...): a plain signed integer, e.g.
  `-3`, `0`, `42`.
* **Rationals** (`mpq_class`, `cpp_rational`, `safe_rational`, ...): either a
  plain integer, or a fraction `p/q` with no surrounding spaces, e.g. `3/4`,
  `-7/2`, `5`. The fraction is expected in lowest terms on output; on input any
  `p/q` is accepted and normalised.
* **Real-algebraic / quadratic fields** (`Qsqrt2`, `Qsqrt5`,
  `RealAlgebraic=...`): see the arithmetic-specific notes at the top of the
  program's usage message; these types have their own token syntax.

### Example

A 3x4 matrix with rational entries:

```
3 4
 1 0 0 0
 0 1/2 0 0
 -1 0 3 7/3
```

The leading space before each entry on output is cosmetic; on input it is
irrelevant.

### Geometric convention

Most programs interpret a matrix geometrically as a family of vectors, one per
**row**, in homogeneous coordinates (the first column is the constant term).
The precise meaning — vertices of a polytope, inequalities of a cone, a Gram
matrix, a generator family — is fixed by each program and stated in the
per-directory `README.md` and in the program's own usage message (run it with
no arguments to print that message).


The list-of-matrices format
---------------------------

Read by `ReadListMatrix` / `ReadListMatrixFile`, written by `WriteListMatrix`.
Used wherever a program consumes or produces several matrices (for example a
list of group generators, or a family of forms).

Layout:

```
<n_mat>
<matrix 0 in the native format>
<matrix 1 in the native format>
...
```

* The **first number** is the count of matrices.
* It is followed by that many matrices, **each in the native matrix format
  above** (i.e. each begins with its own `<nbRow> <nbCol>` header). The
  matrices need not share dimensions.

The same whitespace-insensitivity applies throughout.


Output formats for other systems
--------------------------------

Several programs can emit their result in a form directly consumable by another
system, usually selected by a `choice` argument (e.g. `GAP`, `CPP`, `Number`).

### GAP

`WriteMatrixGAP` emits a matrix as a nested GAP list:

```
[ [ 1, 0, 0, 0 ],
[ 0, 1/2, 0, 0 ],
[ -1, 0, 3, 7/3 ] ]
```

Rows are separated by `,` and, unless the "one line" variant is used, a
newline. When a program writes such output to a file it wraps it as a complete
GAP statement:

```
return [ [ ... ], ... ];
```

so the file can be read directly with GAP's `ReadAsFunction` / `Read`. A list
of matrices (`WriteListMatrixGAP`) is emitted as a GAP list of such matrices,
again optionally wrapped in `return ...;`.

### Python

`WriteListMatrixPYTHON` emits the analogous nested Python-list structure, for
consumption from Python tooling.


Vectors
-------

Vectors follow the same "count, then entries" philosophy as matrices and are
parsed with the same stream mechanism. When a program expects a vector rather
than a matrix, its usage message states so explicitly; when in doubt a vector
of length `n` can be supplied as a `1 x n` matrix in the native format.


Practical notes
---------------

* **Discovering the exact I/O of one program.** Run the program with no
  arguments (or the wrong number of arguments). Every command-line program
  prints a usage message naming its input files, its `arith` choices, and its
  output modes.
* **Files vs. standard streams.** Most programs take an input file path and an
  optional output file path; when the output path is omitted the result is
  written to standard error / standard output as documented in the usage
  message.
* **Choosing the arithmetic.** Prefer the cheapest arithmetic that is safe for
  your data. `safe_rational` uses machine integers and aborts cleanly on
  overflow; `mpq_class` / `cpp_rational` are unbounded but slower. This choice
  is made per run via the `arith` argument, never hard-coded in the file.
