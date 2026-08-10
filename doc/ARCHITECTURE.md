Architecture Overview
=====================

This document is a map of the code: what the main building blocks are, how they
stack, and where to start reading. It is meant to orient someone who did not
write the code. For the meaning of the individual programs in a directory, see
that directory's `README.md`; for file formats, see `doc/FILE_FORMATS.md`; for
the debugging / timing conventions, see `DEVELOPER.md`.


What this project is
--------------------

A collection of C++ tools for polytopes, lattices, and quadratic forms, aimed at
solving *record-size* problems. Two design decisions shape everything else:

* **Header-only, templated on the arithmetic.** Nearly all the logic lives in
  `.h` files and is templated on a number type `T` (chosen per run: machine
  rationals, GMP, boost multiprecision, quadratic / real-algebraic fields, ...).
  The `.cpp` files are thin command-line drivers that pick a `T`, read input,
  call into the headers, and write output.
* **Symmetry is exploited pervasively.** Groups (permutation and matrix) are not
  an afterthought; canonical forms, stabilisers, and orbit enumeration are
  woven through the core so that highly symmetric problems become tractable.

Because the code is header-only, the directory "dependencies" below are really
*include-path availability*, not a link-time DAG. The graph is intentionally not
acyclic — in particular `src_poly` and `src_group` include each other and form a
single tightly-coupled core.


The foundation: two git submodules
-----------------------------------

Everything rests on two upstream projects vendored as git submodules. **Do not
edit files inside them in place** — changes are overwritten on the next
submodule update. To change them, edit the upstream working copy, push there,
then run `./update.sh` to advance the pointer (see `CLAUDE.md`).

* **`basic_common_cpp/`** — the general-purpose C++ foundation, exposed to the
  rest of the tree through the symlinked directories `src_basic` (utilities,
  files, threading, interrupts), `src_number` (the arithmetic types),
  `src_matrix` (dense/sparse matrices, `MyMatrix`, the matrix I/O of
  `doc/FILE_FORMATS.md`), `src_comb` (bitsets, hashing, combinatorics), and
  `src_graph` (graphs).
* **`permutalib/`** — the permutation-group library (orbits, stabilisers,
  backtrack search) used wherever symmetry is applied.


The layers
----------

From the bottom up. Each layer may use anything below it.

```
  ┌──────────────────────────────────────────────────────────────────────┐
  │  Domain solvers (the "record problem" applications)                    │
  │  src_perfect  src_delaunay  src_lorentzian  src_indefinite             │
  │  src_ctype    src_copos      (+ WIP: src_rankin, src_robust_covering,  │
  │                               src_k_coverings, src_single_delaunay,    │
  │                               src_poincare_polyhedron)                 │
  ├──────────────────────────────────────────────────────────────────────┤
  │  Enumeration backbone:  src_enum_schemes   src_dualdesc                │
  │  (orbit / adjacency / cell-complex enumeration, serial + MPI)          │
  ├──────────────────────────────────────────────────────────────────────┤
  │  Lattice core:  src_latt   src_short        (+ src_polydecomp)         │
  ├──────────────────────────────────────────────────────────────────────┤
  │  Polyhedral & symmetry core:  src_poly  <——>  src_group                │
  ├──────────────────────────────────────────────────────────────────────┤
  │  Foundation-only helpers:  src_isotropy  src_sparse_solver             │
  │                            src_polarization  src_polygen               │
  ├──────────────────────────────────────────────────────────────────────┤
  │  Foundation (submodules): basic_common_cpp (src_basic/number/matrix/   │
  │                           comb/graph)  +  permutalib                    │
  └──────────────────────────────────────────────────────────────────────┘
```

### Foundation-only helpers

These depend on the submodules alone and are reused everywhere:

* **`src_isotropy`** — testing whether a quadratic form is isotropic and finding
  an isotropic vector.
* **`src_sparse_solver`** — sparse solver for systems of linear equations.
* **`src_polarization`**, **`src_polygen`** — polarization and polytope
  generation helpers.

### Polyhedral & symmetry core

The heart of the system, and mutually dependent by design:

* **`src_poly`** — polyhedral computations: dual description (facets <-> vertices
  via cddlib / lrs / beneath-and-beyond), linear programming, redundancy
  removal, integral points. See `src_poly/README.md`.
* **`src_group`** — automorphism groups, canonical forms, and (in)equivalence of
  polytopes and weighted structures (`WeightMatrix`, `PolytopeEquiStab`), built
  on `permutalib` and `nauty`/`bliss`.

### Lattice core

* **`src_latt`** — fundamental lattice computations: canonicalization of
  positive-definite forms and shortest-vector enumeration. See
  `src_latt/README.md`.
* **`src_short`** — short-vector related computations.
* **`src_polydecomp`** — polytope-decomposition helpers used by the Delaunay
  machinery.

### Enumeration backbone

The generic engine that turns "enumerate all orbits of X adjacent to Y" into a
reusable, parallelisable computation:

* **`src_enum_schemes`** — the abstract **adjacency scheme**, **cell complex**,
  and **equivariant decomposition** frameworks, each with a serial and an
  **MPI** variant (`POLY_AdjacencyScheme.h` / `POLY_MPI_AdjacencyScheme.h`,
  `POLY_CellComplex.h` / `POLY_MPI_CellComplex.h`).
* **`src_dualdesc`** — dual description at scale: the symmetry-quotiented
  **recursive dual description** (`POLY_RecursiveDualDesc.h`), serial and MPI
  drivers, and a disk-backed databank (`Databank*.h`) for computations larger
  than memory.

### Domain solvers

Each targets a specific mathematical problem and is built on all of the above.
Note they also build on *each other*: `src_lorentzian`, `src_ctype`, and
`src_copos` all use `src_perfect`; `src_delaunay` and `src_robust_covering` use
`src_polydecomp` / `src_polygen`.

* **`src_perfect`** — perfect forms.
* **`src_delaunay`** — Delaunay polytopes and the space of Delaunay
  tessellations (iso-Delaunay / L-types).
* **`src_lorentzian`** — Vinberg / Edgewalk algorithms for hyperbolic forms and
  perfect hyperbolic forms.
* **`src_indefinite`** — reduction of indefinite forms (LLL-based, isotropic).
* **`src_ctype`** — C-types (used to enumerate all C-types in dimension 6).
* **`src_copos`** — copositivity and strict copositivity.
* **Work in progress:** `src_rankin` (Rankin constants), `src_single_delaunay`,
  `src_poincare_polyhedron` (Poincaré polyhedron theorem for tilings),
  `src_robust_covering` (robust covering density), `src_k_coverings`.

### Interoperability

* **`src_export_oscar`**, **`src_export_wasm`** — bridges to OSCAR (Julia) and
  WebAssembly. See also the top-level `julia/`, `PythonScript/`, and GAP tooling.


Cross-cutting concerns
----------------------

* **Arithmetic as a template parameter.** The number type is chosen at the
  command line, never hard-coded. The same header serves machine-integer,
  GMP, and quadratic/real-algebraic runs. The exact types (`mpq_class`,
  `mpz_class`) are unbounded and always correct; the kernels transparently
  attempt the hot inner computations over machine integers first (the
  `TryInt64` deferred-overflow path) and fall back to the exact type on
  overflow, so machine-integer speed no longer requires an unsafe type
  choice.
* **Symmetry.** Canonicalization and stabiliser computations (`src_group`,
  `permutalib`) let the domain solvers enumerate one representative per orbit
  instead of the full set — the difference between feasible and infeasible on
  the target problems.
* **Parallelism.** Distributed runs use **MPI** and live in `src_dualdesc` and
  the MPI variants of `src_enum_schemes` (and are exercised by `src_ctype`).
  Everything else is serial; the two models share the same headers.
* **Diagnostics.** Output is gated behind named `DEBUG_<CONTEXT>`,
  `SANITY_CHECK_<CONTEXT>`, and `TIMINGS_<CONTEXT>` macros — see `DEVELOPER.md`
  and `CLAUDE.md` for the conventions.


The build
---------

Each `src_*` directory is independently buildable: its `Makefile` declares its
include paths (the `CLOCAL` line — the authoritative dependency list for that
module) and its targets. There is also a top-level `CMakeLists.txt` covering the
tree. External dependencies are GMP, Boost, Eigen, cddlib, and nauty (plus
optional FLINT / bliss / MPI); the submodules supply the rest.


Where to start reading
----------------------

1. `README.md` (top level) for the one-paragraph tour of each directory.
2. `doc/FILE_FORMATS.md` to be able to construct an input and read an output.
3. The `README.md` of the directory you care about, then run its program with no
   arguments to see the usage message.
4. For the mathematics behind a module, follow the literature references in that
   `README.md`.
5. To understand an *enumeration* (perfect forms, Delaunay, C-types...), read
   `src_enum_schemes/POLY_AdjacencyScheme.h` first — it is the shared skeleton
   that each domain solver specialises.
