# Joint (Q, c) descent for periodic covering densities

Numerical local optimization of the covering density of a 2-point periodic
set `Z^5 + {0, c}` over the **joint** space of the form `Q` and the coset
position `c`. Unlike the exact iso-Delaunay machinery of `src_delaunay/`,
`c` is a continuous variable here — no denominator, and no requirement that
`-I` preserve the point set — so this searches the full 2-point family, at
the price of being numerical rather than exact.

* `evaluator.py` — covering density via a genuine periodic Delaunay
  triangulation (scipy/qhull) in the metric `Q`. The enumeration box is
  LLL-reduced and grows automatically when a cell touches its shell, since a
  too-small box makes qhull fill the star with spurious cells whose
  circumradii over-estimate the covering radius (`rng=2` sufficed for a
  well-reduced form, `rng=3` was needed for walk optima; the growth is
  iterative and releases the previous triangulation — a recursive version
  cost 14 GB per process). Validated: it reproduces the certified
  2.160060765053 to 12 digits.
* `descend.py` — outer loop of (re-tessellate, SLSQP on
  `(cholesky(Q), c)` with the cell list fixed), the constraints being the
  squared circumradii of all translation classes of cells. Each round's
  result is verified against a fresh tessellation, so a stale cell list
  cannot produce a fake improvement.
* `multistart.py` — random restarts of the descent.

Requires `numpy` and `scipy` (`python3 -m venv venv && venv/bin/pip install
numpy scipy`). Everything is specific to dimension 5 through the constant
`N = 5` in `evaluator.py` and the hypersimplex constant `KAPPA5`.

## What it established (September 2026)

Seeding the descent at the per-family flip-walk optima of
`doc/COVERING_OPTIMIZATION.md` and freeing `c`:

| seed (family, walk value) | joint descent reaches | c moved |
| --- | --- | --- |
| N=6, 2.199698 | **2.160060765** | 0.25 |
| N=5, 2.517200 | 2.241145 | 0.42 |
| N=5, 2.507244 | 2.314943 | 0.36 |
| N=5, 2.510016 | 2.454742 | 0.09 |

The N=6 optimum drains into the same attractor as the N=4 family: the
certified `Q(sqrt2)` point at density 2.160060765053 is a local optimum of
the full continuum problem (576 of 864 cell classes simultaneously active,
no joint descent direction), reached from independent families. The N=5
walk values were 0.2–0.27 above their continuum limits — for that family
the rational denominators genuinely were the binding constraint. About a
hundred random multistarts landed only in shallower basins (2.5–2.7).
No configuration below 2.160060765053 was found; the record to beat is
`Theta(A_5^*) = 2.124285909`.
