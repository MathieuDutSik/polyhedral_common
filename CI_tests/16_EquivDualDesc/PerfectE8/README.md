# Perfect domain of E8

Same construction as `../PerfectE7`, one dimension up: the cone spanned by the
rank one forms `v v^T` of the 120 pairs of minimal vectors of the E8 root
lattice, full dimensional in the 36 dimensional space of symmetric 8x8
matrices. `PerfectE8.ext` has 120 vertices and 36 columns, `PerfectE8.grp` is
`Aut(E8) / {+-1}`, of order 348364800. Regenerate both with

```
python3 ../GeneratePerfectCone.py E8 PerfectE8.ext
../../../src_group/GRP_LinPolytope_Automorphism rational PerfectE8.ext Oscar PerfectE8.grp
```

This is the computation the recursive adjacency decomposition method was
originally built for: the dual description has 25075566937584 facets in 83092
orbits. With the heuristics below it runs from scratch in **2h17m** on one
core. It is **not part of the CI**: the entry stays commented out in
`TestCases.g`, next to CUT_K8 and CUT_K55, and is meant to be run by hand.
The shape of the run is: everything found in the first half hour, 83091 of the
83092 orbits done within two hours, and the last twenty minutes spent on the
single orbit 229, of incidence 75 and stabilizer 23040.

`NumericalType = "integer"` is set. It makes no difference to the result (same
83092 orbits, identical as a set) and none to the speed either: 2h23m against
the 2h17m of `rational`, one run each, which is inside what this machine
resolves. The reason is that the kernels already run on integers in the
rational setting -- both the reverse search and the normaliz port scale each
row into the underlying ring once and then work over TryInt64 -- so the
setting only removes the per-subproblem conversion in the surrounding layer.

* `input.nml` is the entry used by `TestCases.g`, saving disabled.
* `PerfectE8_saving.nml` has `Saving = T` for both the polyhedral database and
  the bank and a `max_runtime`, so the run checkpoints and is restarted by the
  loop in `run_benchmark.sh` until the `orbits` file appears.

## What the facets look like

`POLY_sampling_facets rational lp_cdd:iter_300 PerfectE8.ext` finds facets of
incidence 42 up to 75 (the rank of a facet is 35), with stabilizers of order 2
for the incidence 42 ones up to 23040 for the incidence 75 ones. Since the
total facet count divided by the number of orbits is 3.0e8, just below
`|GRP| = 3.5e8`, almost every one of the 83092 orbits has a stabilizer of order
1 or 2, so the bulk of the work is on facets of incidence in the low forties.

Those are found quickly: all 83092 orbits are *found* within about ninety
minutes. What costs is *finishing* them, that is running the adjacency
decomposition on each so that the Balinski criterion can certify that no facet
is missing. The last few dozen orbits, of incidence 47 and above, dominate the
run.

## Where the thresholds come from

Cost of a *direct* dual description of one facet (`POLY_dual_description
mpq_class <prog> CPP`), measured on facets actually met in the run:

| incidence | delta | ridges | cdd | lrs |
|-----------|-------|--------------|--------|--------|
| 42 | 7 | 2096 | 0.15s | 0.07s |
| 45 | 10 | 14961 | 0.24s | 0.55s |
| 47 | 12 | 39529 | 0.7s | |
| 48 | 13 | 117924 | 4.1s | 7.7s |
| 49 | 14 | 185187 | 7.9s | |
| 50 | 15 | 325110 | 25.9s | 27.5s |
| 51 | 16 | 1821563 | >8182s | 2362s |
| 52 | 17 | 2456177 | >8561s | 2400s |
| 54 | 19 | 696621 | 1234s | |

Hence:

* `Split.heu` decides on `delta` *and* on the size of the stabilizer. `delta`
  rather than `incidence` makes the rule scale free, since the same rule then
  applies at every level of the recursion, where the rank drops. But `delta`
  alone is not enough, and getting the rule wrong is expensive in both
  directions:

  - Too eager to split. At `delta <= 11` the incidence 47 to 50 facets went to
    the recursion and cost about half an hour each, against the second or so of
    the table above: an 8x loss on the tail of the run.
  - Too eager to compute directly. The stabilizer is what the recursion trades
    on, since it divides the number of ridges that have to be flipped. The nine
    facet orbits that are left at the end of the run are

    | orbit | incidence | delta | stabilizer |
    |-------|-----------|-------|------------|
    | 170   | 54 | 19 | 24 |
    | 749   | 54 | 19 | 120 |
    | 1010  | 54 | 19 | 216 |
    | 236   | 57 | 22 | 48 |
    | 152   | 58 | 23 | 240 |
    | 5292  | 60 | 25 | 7200 |
    | 216   | 66 | 31 | 2304 |
    | 552   | 70 | 35 | 5040 |
    | 229   | 75 | 40 | 23040 |

    and orbit 5292 has a stabilizer of 7200, so its ridge work is divided by
    7200 by the recursion. Sent to a direct lrs instead it ran for more than two
    hours without finishing, while a standalone lrs on the incidence 57 one,
    of stabilizer 48, takes 1555 seconds. A rule on `delta` alone puts those two
    on the same side, which is wrong.
* Above `delta = 15` the table stops being monotone in the incidence: the
  incidence 51 and 52 facets are harder for cdd than the incidence 54 one, and
  cdd does not finish on them at all while lrs does, three times faster.
  `DualDesc_heu.ts` therefore uses two programs and no sampling: lrs below
  `delta = 16`, normaliz above it, and normaliz whenever the incidence exceeds
  44 whatever the delta. That last guard matters because the hinge was measured
  at rank 35, and deeper in the recursion the rank has dropped so `delta` alone
  understates the size of the output. Leaving cdd in the sampler as a third
  option cost a factor of five on the whole run: the sampler state is not
  persisted across the checkpoint restarts, so cdd was re-explored on giant
  subpolytopes every time, taking 17 of the 45 direct dual descriptions above
  200000 facets at 2x to 10x the normaliz time.
* `Bank.heu` / `CheckBank.heu` bank and query at `delta >= 12`, that is the
  subpolytopes that are expensive enough to be worth a canonical form lookup.
* `InitFacet.heu` uses `lp_cdd_min` so the enumeration starts on a facet of
  minimum incidence (42) rather than on one of the incidence 75 ones that plain
  `lp_cdd` overwhelmingly returns.
* `OrbitSplit.heu` forces `canonic`. The default sends `groupsize_sma < 100` to
  the exhaustive method, and the stabilizers here are of order 1 or 2, so the
  default applies exactly where it must not: on the splitting of a bank entry
  with `|BigGRP|=737280`, `|SmaGRP|=2` and 60140 orbits, the exhaustive method
  does not finish while `canonic` takes 85 seconds. That case is kept as
  `CI_tests/DoubleCosets/DBL/DoubleCoset_n48_big737280_sma2_vf60140_idx0` if it
  is regenerated; it is far outside the rest of that collection, where every
  file has at most 414 orbits and an index of at most 16.

The canonicalization method is left at the default (`ChoiceCanonicalizationFile`
unset, that is `groupsize < 500 store`, `canonic` otherwise). Measured over
equal wall clock at the top level, `canonic` treated 752 orbits, `guess` 371 and
`canonic_initial_triv` 2, so the default is the right choice here despite most
facet stabilizers being tiny.

The Balinski based `AdvancedTerminationCriterion` is on: it lets the enumeration
stop once the facets found are provably all of them, well before every one of
the 83092 orbits has been through the adjacency decomposition.

Note that the checked in `Makefile_serial_dualdesc` compiles with
`-DDEBUG_DUAL_DESC`, which prints two lines per direct dual description. On a
run that performs tens of millions of them that logging is not free; build with
`make -f Makefile_serial_dualdesc COPTIONS=-DGMP_POOL` for the benchmark.
