# Perfect domain of E7

The perfect domain (Voronoi cone) of the perfect form E7 is the cone spanned by
the rank one forms `v v^T`, with `v` running over the 63 pairs of minimal
vectors of the E7 root lattice. As E7 is a perfect form those rank one forms
span the whole space of symmetric 7x7 matrices, so the cone is full dimensional
of dimension 28.

All rank one forms of minimal vectors satisfy `<E7, v v^T> = min(E7)`, so the
rays lie in an affine hyperplane and the cone is the homogenization of a
polytope. `PerfectE7.ext` is written in the usual "first column equal to 1"
format: 63 vertices, 28 columns, obtained by trading the coordinate `X(0,0)`
against the constant functional `<E7, X> / min(E7)`. It is regenerated with

```
python3 ../GeneratePerfectCone.py E7 PerfectE7.ext
../../../src_group/GRP_LinPolytope_Automorphism rational PerfectE7.ext Oscar PerfectE7.grp
```

`PerfectE7.grp` is `Aut(E7) / {+-1}`, of order 1451520.

The dual description has **157 orbits of facets**, 79900912 facets in total,
with incidences ranging from 27 (simplicial facets, 91 orbits) to 46.

The run takes about 12 seconds. The two heuristics that matter:

* `Split.heu` sends everything of incidence at most 45 to the direct dual
  description. The single facet orbit of incidence 46 costs about 170 seconds
  when handled directly by cdd, and a couple of seconds when the recursive
  adjacency decomposition is applied to it instead.
* `DualDesc_heu.ts` routes `delta <= 1` (simplices and near simplices) to
  `small_polytopes`. 91 of the 157 facet orbits are simplices, so this is the
  bulk of the direct dual descriptions performed.
