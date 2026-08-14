# TODO

## Integer arithmetic of the dual-description kernels (normaliz / lrs / bb / dd)

Background: the kernels currently run their machine-integer fast path over the
deferred-overflow TryCarryInt64 (exact per-operation detection, exact-ring
fallback). Micro-benchmarks (Apple Silicon) put the per-operation checking
cost at 1.6x-3x over raw int64_t on the hot loop shapes (scalar product, FM
combination, Bareiss row); end-to-end the arithmetic share is smaller, so the
observed whole-run overhead against a raw-int64 implementation is ~1.1x-1.4x.
The reference Normaliz long-long mode is NOT rigorous (unchecked hand-unrolled
scalar products; range checks applied to results after products that may
already have wrapped, which is UB); we do not want to copy that. The items
below are the rigorous ways to close the gap.

* **Hadamard-certified raw int64 for dimension <= 15.** For a submatrix with
  entries bounded by G, every Bareiss intermediate is a product of two
  consecutive minors, so Had_k(G) * Had_{k+1}(G) < 2^63 with
  Had_k(G) = (G*sqrt(k))^k certifies the whole elimination a priori --
  no runtime checks at all. With G = 1 this holds up to dimension ~15
  (fails at 16: ~1.4e20). Same regime covers scalar products and FM
  combinations (facet coefficients bounded by Had_{d-1}(G)). Applicable to
  all four backends on low-dimensional (sub)problems; useless at d ~ 98,
  where the Hadamard bound is ~1e96 while the measured coefficients on the
  P_9_6 hard family are <= 96.

* **Precondition-gated raw int64_t for high dimension.** Runtime-gated, but
  rigorous (contrast Normaliz): verify the precondition BEFORE running a
  block raw, never check results after possibly-wrapped operations.
  Invariants: generator entries <= G (checked at input conversion), stored
  facet coefficients <= B (checked at the number_hyperplane choke point,
  after gcd reduction), with B = sqrt(2^62 / (dim * G)); then every
  intermediate of a scalar product (<= dim*G*B) and of an FM combination
  (<= 2*dim*G*B^2) is provably wrap-free. For Bareiss eliminations track the
  max workspace entry E per step (one compare per write) and require
  2*E^2 < 2^63 before each row update; violation falls back per call
  (rank: to the comparison test), per object, or per run (to TryCarryInt64,
  then the exact ring). Measured headroom on the P_9_6 hard family:
  facet coefficients <= 96 against a budget of ~2.2e8 (the
  NMZ_COEFF_STATS instrumentation in POLY_DualDesc_normaliz.h measures
  this; the max-Bareiss-entry counterpart is still to be measured).
  Expected gain: the arithmetic share only, est. 5-15% on the hard family.

* **NEON widening multiplies for the classification loop.** AArch64 NEON has
  no 64x64 vector multiply (and AVX2 none either; vpmullq is AVX-512DQ), so
  int64 loops cannot vectorize on current targets -- this is why
  TrySimdInt64 cannot beat TryCarryInt64 there (measured: it is 10-30%
  slower; kept selectable via POLY_NORMALIZ_TRY_SIMD / POLY_LRS_TRY_SIMD for
  future ISAs). But 32x32->64 widening multiplies (smull/smlal) DO exist:
  when the tracked bounds certify coefficients < 2^31 (P_9_6 family: <= 96),
  storing facet normals as int32 makes the scalar products of the
  classification loop vectorizable 2-4 wide with 64-bit accumulation.
  To be pursued after the precondition-gated int64 mode exists, since it
  reuses the same bound bookkeeping.
