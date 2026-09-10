# Joint (Q, c) descent for periodic covering densities

Numerical local optimization of the covering density of a 2-point periodic
set `Z^5 + {0, c}` over the **joint** space of the form `Q` and the coset
position `c`. Unlike the exact iso-Delaunay machinery of `src_delaunay/`,
`c` is a continuous variable here — no denominator, and no requirement that
`-I` preserve the point set — so this searches the full 2-point family, at
the price of being numerical rather than exact.

* `evaluator.py` — covering density via a genuine periodic Delaunay
  triangulation (scipy/qhull) in the metric `Q`, for any number of cosets
  `m`. The point cloud is the Q-ball `|x|_Q <= R` enumerated through an
  LLL-reduced basis, and every kept cell is required to have its whole
  circumsphere inside the ball (`|center|_Q + radius <= R`), which makes it
  provably a Delone cell of the full periodic set; if a star cell fails the
  test, `R` grows. Earlier versions enumerated a hypercube in the reduced
  basis, which both wastes a factor ~50 in volume in dimension 5 and never
  stabilizes for skewed bases where genuine cells reach the shell along
  near-cancelling combinations. Cells are anchored on any lattice-0 vertex,
  not only coset-0 ones — the set is not translation-symmetric by a coset,
  so cells whose vertices all lie in nonzero cosets are genuine classes.
  Validated: reproduces the certified m=2 value 2.160060765053 to 12 digits
  and the exact C++ value 2.3765773479972 of an m=3 configuration digit for
  digit.
* `descend.py` — outer loop of (re-tessellate, SLSQP on
  `(cholesky(Q), c)` with the cell list fixed), the constraints being the
  squared circumradii of all translation classes of cells. Each round's
  result is verified against a fresh tessellation, so a stale cell list
  cannot produce a fake improvement.
* `multistart.py` — random restarts of the descent.

* `m3_campaign.py` — the m=3 exploration: seeded at the m=3 flip-walk
  optimum, at the m=2 optimum plus a random third coset, and at random.

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

For `m = 3` the continuum campaign found the attractor 2.3398 below the
flip-walk value 2.3766, and — pointedly — seeding at the m=2 optimum plus a
random third coset always stalls near `1.5 x 2.16006`: the third point's
density factor is never recovered locally. Note this is a statement about
the descents, not a principle -- an optimal covering can have cells below
the covering radius (the dimension 6 record `L^c_6` does), and the
2.160060765053 optimum itself has one slack orbit.
