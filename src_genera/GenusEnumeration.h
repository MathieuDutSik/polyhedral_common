// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSENUMERATION_H_
#define SRC_GENERA_GENUSENUMERATION_H_

// clang-format off
#include "Boost_bitset.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "LatticeStabEquiCan.h"
#include "LatticeRootDecomposition.h"
#include "LatticeAutomorphismVinberg.h"
#include "LatticePleskenSouvignier.h"
#include "InvariantVectorFamily.h"
#include "ClassicLLL.h"
#include "Positivity.h"
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_GENUS_ENUMERATION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_GENUS_ENUMERATION
#endif

#ifdef TIMINGS
#define TIMINGS_GENUS_ENUMERATION
#endif

#ifdef TRACK_INFO
#define TRACK_INFO_GENUS_ENUMERATION
#endif

/*
  Enumeration of the isometry classes of a genus of even positive definite
  lattices, by Kneser p-neighbours.

  The scheme is NOT the adjacency-decomposition scheme of src_perfect and
  src_delaunay, and deliberately does not use their DataFunc machinery. There
  the orbits are joined by a geometric adjacency relation and the enumeration
  terminates when the orbit list is closed under it. Here every lattice of the
  genus is adjacent to a fixed lattice through some p-neighbour, so there is no
  useful notion of a facet, and the termination criterion is arithmetic: the
  Smith-Minkowski-Siegel mass

      mass(G) = sum over the classes L of the genus of 1 / |Aut(L)| .

  The enumeration stops when the accumulated sum reaches the mass supplied on
  input. That also certifies completeness independently of the enumeration, so
  a bug in the neighbour construction cannot silently return a short list: it
  shows up as a mass that never closes.

  Two facts govern the design.

  (1) The p-neighbour graph is connected on a SPINOR genus, not on a genus. A
  genus with several spinor genera therefore needs one seed lattice per spinor
  genus; with a single seed the enumeration converges to a proper subset and
  the mass falls short. The program reports that rather than pretending to be
  complete, and accepts several seeds on input.

  (2) The deduplication is by canonical form (src_latt/LatticeStabEquiCan.h),
  so testing whether a new neighbour is already known is a hash lookup rather
  than an isometry test against every class found so far. That is the point of
  doing this here: the pairwise scheme is quadratic in the class number.
 */

/*
  How the neighbour primes are scheduled.

  OutPrime: one prime for the whole run -- the one given, or the smallest
  prime NOT dividing the determinant. At such a prime the neighbour graph
  is connected on the spinor genus (Kneser), so the enumeration terminates
  with the full list; the cost is that an even determinant forces p >= 5
  and the projective space grows as p^{n-1}.

  Iterate: primes 2, 3, 5, 7, ... in order, including those dividing the
  determinant. A p-neighbour is in the genus for EVERY p (Schulze-Pillot,
  properties (i)-(iv): same determinant, integral, even, same genus), so
  every class found is genuine at any prime; what is lost at p | det is
  only the guarantee that one prime connects the graph. The scheduler
  always works at the smallest prime that still has unexpanded classes, so
  classes discovered at an expensive prime are expanded at the cheap
  primes first, and a new prime enters only when every known class has
  been expanded at every prime in play and the mass has not closed.
  Deterministic, and termination is inherited from OutPrime: the schedule
  eventually reaches a prime not dividing the determinant.
 */
enum class NeighborStrategy { OutPrime, Iterate };

// The description of the genus. `prime` is the prime used for the neighbour
// step; 0 means "choose the smallest prime not dividing 2*det". With the
// Iterate strategy the prime is ignored.
template <typename T> struct GenusSpec {
  int rank;
  T det;
  int prime;
  NeighborStrategy strategy;
};

/*
  The accumulated result of an enumeration.

  Three different kinds of number meet here and none of them is the others.
  The Gram matrices are over T. The group orders are counts, so they are of
  the integer type the group itself counts with, Tgroup::Tint, and not of the
  type the matrices happen to have. The mass is a sum of 1 / |Aut(L)| and its
  target is read as a fraction, so it lives in the field over the orders.
 */
template <typename T, typename Tgroup> struct GenusEnumerationResult {
  using TintGroup = typename Tgroup::Tint;
  using Tmass = typename overlying_field<TintGroup>::field_type;
  std::vector<MyMatrix<T>> ListGram;
  std::vector<TintGroup> ListAutOrder;
  Tmass accumulated_mass;
  Tmass target_mass;
  bool complete;
  // The primes the neighbour steps used, in the order they came into play,
  // and how many new classes the neighbours at each of them produced. The
  // seeds are not counted. With the out_prime strategy there is one entry.
  std::vector<int> ListPrimeUsed;
  std::vector<size_t> ListNewClassPerPrime;
};

// Inverse of a modulo the prime p, by extended Euclid. p is small.
inline int InverseModP(int a, int p) {
  int a_red = ((a % p) + p) % p;
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (a_red == 0) {
    std::cerr << "GENUS: InverseModP called on a multiple of p\n";
    throw TerminalException{1};
  }
#endif
  int r0 = p, r1 = a_red, s0 = 0, s1 = 1;
  while (r1 != 0) {
    int q = r0 / r1;
    int r2 = r0 - q * r1;
    int s2 = s0 - q * s1;
    r0 = r1;
    r1 = r2;
    s0 = s1;
    s1 = s2;
  }
  return ((s0 % p) + p) % p;
}

template <typename Tint> int ResidueModP(Tint const &x, int p) {
  Tint p_T(p);
  Tint r = x % p_T;
  int r_int = UniversalScalarConversion<int, Tint>(r);
  return ((r_int % p) + p) % p;
}

// The p-adic valuation of a non-zero integer.
template <typename Tint> int ValuationP(Tint const &x, int p) {
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (x == 0) {
    std::cerr << "GENUS: ValuationP called on zero\n";
    throw TerminalException{1};
  }
#endif
  Tint p_T(p);
  Tint work = x;
  if (work < 0) {
    work = -work;
  }
  int val = 0;
  while (work % p_T == 0) {
    work /= p_T;
    val++;
  }
  return val;
}

/*
  The points of the projective space P^{n-1}(F_p), indexed WITHOUT being
  enumerated. There are (p^n - 1) / (p - 1) of them: 16383 for n = 14 and
  p = 2, but 1.5 x 10^9 at p = 5, where materialising the list (as an
  earlier version did) is itself a hundred gigabytes. The representatives
  are the vectors whose first non-zero coordinate is 1; the index runs
  through the blocks ipos = 0, ..., n-1 (position of that coordinate), each
  block of size p^{n-ipos-1}, the free coordinates read little-endian
  (v[ipos+1] least significant). ProjectiveLineOfIndex and
  IndexOfProjectiveLine are mutually inverse bijections for this order.
 */
inline size_t NumberProjectiveLines(int n, int p) {
  size_t const max_val = std::numeric_limits<size_t>::max();
  size_t p_s = static_cast<size_t>(p);
  size_t total = 0;
  size_t power = 1;
  for (int k = 0; k < n; k++) {
    if (total > max_val - power) {
      std::cerr << "GENUS: (p^n - 1)/(p - 1) overflows size_t for p=" << p
                << " n=" << n << "\n";
      throw TerminalException{1};
    }
    total += power;
    if (k + 1 < n) {
      if (power > max_val / p_s) {
        std::cerr << "GENUS: p^" << (k + 1) << " overflows size_t for p=" << p
                  << "\n";
        throw TerminalException{1};
      }
      power *= p_s;
    }
  }
  return total;
}

inline std::vector<int> ProjectiveLineOfIndex(int n, int p, size_t index) {
  size_t p_s = static_cast<size_t>(p);
  // block size p^{n-1} of the block ipos = 0
  size_t block = 1;
  for (int i = 0; i < n - 1; i++) {
    block *= p_s;
  }
  size_t work = index;
  int ipos = 0;
  while (work >= block) {
    work -= block;
    block /= p_s;
    ipos++;
#ifdef SANITY_CHECK_GENUS_ENUMERATION
    if (ipos == n) {
      std::cerr << "GENUS: ProjectiveLineOfIndex, index out of range\n";
      throw TerminalException{1};
    }
#endif
  }
  std::vector<int> v(n, 0);
  v[ipos] = 1;
  for (int i = ipos + 1; i < n; i++) {
    v[i] = static_cast<int>(work % p_s);
    work /= p_s;
  }
  return v;
}

inline size_t IndexOfProjectiveLine(std::vector<int> const &v, int p) {
  size_t p_s = static_cast<size_t>(p);
  int n = static_cast<int>(v.size());
  int ipos = 0;
  while (ipos < n && ((v[ipos] % p) + p) % p == 0) {
    ipos++;
  }
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (ipos == n) {
    std::cerr << "GENUS: IndexOfProjectiveLine called on the zero vector\n";
    throw TerminalException{1};
  }
#endif
  int scal = InverseModP(v[ipos], p);
  // start of the block: sum_{k < ipos} p^{n-k-1}
  size_t block = 1;
  for (int i = 0; i < n - 1; i++) {
    block *= p_s;
  }
  size_t index = 0;
  for (int k = 0; k < ipos; k++) {
    index += block;
    block /= p_s;
  }
  // offset inside the block, little-endian in the normalised coordinates
  size_t mult = 1;
  for (int i = ipos + 1; i < n; i++) {
    int val = ((((v[i] % p) + p) % p) * scal) % p;
    index += mult * static_cast<size_t>(val);
    mult *= p_s;
  }
  return index;
}

// Kept only for the sanity check of the indexing, as an INDEPENDENT
// enumeration of the same order (small cases only); the enumeration itself
// never materialises the projective space.
inline std::vector<std::vector<int>> EnumerateProjectiveLines(int n, int p) {
  std::vector<std::vector<int>> ListLines;
  for (int ipos = 0; ipos < n; ipos++) {
    int n_free = n - ipos - 1;
    int64_t n_case = 1;
    for (int i = 0; i < n_free; i++) {
      n_case *= p;
    }
    for (int64_t idx = 0; idx < n_case; idx++) {
      std::vector<int> v(n, 0);
      v[ipos] = 1;
      int64_t work = idx;
      for (int i = ipos + 1; i < n; i++) {
        v[i] = static_cast<int>(work % p);
        work /= p;
      }
      ListLines.push_back(v);
    }
  }
  return ListLines;
}

#ifdef SANITY_CHECK_GENUS_ENUMERATION
// The indexing against the independent enumeration: exhaustive when the
// space is small, and for a large space the roundtrip index -> line ->
// index on the block boundaries and on a stride through the range.
inline void CheckProjectiveLineIndexing(int n, int p) {
  size_t n_line = NumberProjectiveLines(n, p);
  if (n_line <= 100000) {
    std::vector<std::vector<int>> ListLines = EnumerateProjectiveLines(n, p);
    if (ListLines.size() != n_line) {
      std::cerr << "GENUS: NumberProjectiveLines disagrees with the "
                   "enumeration\n";
      throw TerminalException{1};
    }
    for (size_t i = 0; i < n_line; i++) {
      if (ProjectiveLineOfIndex(n, p, i) != ListLines[i] ||
          IndexOfProjectiveLine(ListLines[i], p) != i) {
        std::cerr << "GENUS: projective line indexing broken at index " << i
                  << "\n";
        throw TerminalException{1};
      }
    }
    return;
  }
  auto roundtrip = [&](size_t i) -> void {
    std::vector<int> v = ProjectiveLineOfIndex(n, p, i);
    if (IndexOfProjectiveLine(v, p) != i) {
      std::cerr << "GENUS: projective line roundtrip broken at index " << i
                << "\n";
      throw TerminalException{1};
    }
  };
  size_t block_start = 0;
  size_t block = n_line - (n_line - 1) / static_cast<size_t>(p);
  // block = p^{n-1}; the boundaries of every block and their neighbours
  for (int ipos = 0; ipos < n; ipos++) {
    roundtrip(block_start);
    if (block_start + block - 1 < n_line) {
      roundtrip(block_start + block - 1);
    }
    block_start += block;
    block /= static_cast<size_t>(p);
  }
  size_t stride = 1 + n_line / 10000;
  for (size_t i = 0; i < n_line; i += stride) {
    roundtrip(i);
  }
}
#endif

/*
  The p-neighbour of the lattice with Gram matrix G attached to an isotropic
  point v of P^{n-1}(F_p):

      N = { y in L : (y, v) = 0 mod p }  +  Z (v/p) .

  For N to be an even lattice one needs v_p((v,v)) >= e + 2 with e = v_p(2),
  that is (v,v) = 0 mod p^2 for odd p and mod 8 for p = 2. A point that is
  isotropic mod p only gives valuation e + 1 in general, and is corrected here
  by v <- v - a p u for a suitable u; that correction is the whole of the p = 2
  difficulty, and is invisible for odd p.

  Returns nothing when v lies in the radical mod p, or when the valuation is
  too small to be corrected (which can only happen for p = 2).
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GenusNeighbor(MyMatrix<Tint> const &G, std::vector<int> const &v_line, int p,
              [[maybe_unused]] std::ostream &os) {
  int n = G.rows();
  MyVector<Tint> v(n);
  for (int i = 0; i < n; i++) {
    v(i) = Tint(v_line[i]);
  }
  // w = G v, the linear form y -> (y, v)
  auto compute_w = [&](MyVector<Tint> const &x) -> MyVector<Tint> {
    MyVector<Tint> w(n);
    for (int i = 0; i < n; i++) {
      Tint sum(0);
      for (int j = 0; j < n; j++) {
        sum += G(i, j) * x(j);
      }
      w(i) = sum;
    }
    return w;
  };
  MyVector<Tint> w = compute_w(v);
  int jpiv = -1;
  for (int i = 0; i < n; i++) {
    if (ResidueModP(w(i), p) != 0) {
      jpiv = i;
      break;
    }
  }
  if (jpiv == -1) {
    // v lies in the radical of the form mod p: no neighbour from it
    return {};
  }
  auto norm_of = [&](MyVector<Tint> const &x, MyVector<Tint> const &wx) -> Tint {
    Tint sum(0);
    for (int i = 0; i < n; i++) {
      sum += x(i) * wx(i);
    }
    return sum;
  };
  int e = (p == 2) ? 1 : 0;
  Tint nrm = norm_of(v, w);
  int val = ValuationP(nrm, p);
  if (val <= e) {
    return {};
  }
  if (val == e + 1) {
    // correct v so that the valuation rises to at least e + 2. With
    // a = nrm / (2 p (v, e_jpiv)) taken mod p, the two leading terms of
    // (v - a p e_jpiv, v - a p e_jpiv) = nrm - 2 a p (v, e_jpiv) + a^2 p^2 (..)
    // cancel. Both nrm and 2 p (v,e_jpiv) have valuation e + 1, so the ratio
    // is a p-adic unit and a is well defined.
    Tint den = Tint(2) * Tint(p) * w(jpiv);
    Tint p_T(p);
    Tint num_red = nrm;
    Tint den_red = den;
    for (int i = 0; i < e + 1; i++) {
      num_red /= p_T;
      den_red /= p_T;
    }
    int a = (ResidueModP(num_red, p) *
             InverseModP(ResidueModP(den_red, p), p)) % p;
    v(jpiv) -= Tint(a) * Tint(p);
    w = compute_w(v);
    nrm = norm_of(v, w);
#ifdef SANITY_CHECK_GENUS_ENUMERATION
    if (ValuationP(nrm, p) < e + 2) {
      std::cerr << "GENUS: the correction of the isotropic vector failed, "
                << "valuation is " << ValuationP(nrm, p) << " and should be at "
                << "least " << (e + 2) << "\n";
      throw TerminalException{1};
    }
#endif
  }
  /*
    A basis of L_v = { y : (y, v) = 0 mod p }: with c_i = w_i / w_jpiv mod p,
    the vectors e_i - c_i e_jpiv for i != jpiv, together with p e_jpiv. That is
    a triangular basis of determinant p, so the index is p as it must be.
   */
  int inv_pivot = InverseModP(ResidueModP(w(jpiv), p), p);
  MyMatrix<Tint> Bsub(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Bsub(i, j) = Tint(0);
    }
  }
  for (int i = 0; i < n; i++) {
    if (i == jpiv) {
      Bsub(i, jpiv) = Tint(p);
    } else {
      int c = (ResidueModP(w(i), p) * inv_pivot) % p;
      Bsub(i, i) = Tint(1);
      Bsub(i, jpiv) = -Tint(c);
    }
  }
  /*
    N = L_v + Z (v/p), so p N = p L_v + Z v. Take a Z-basis of the lattice
    generated by the rows of p * Bsub together with v, then divide by p.
   */
  MyMatrix<Tint> Agen(n + 1, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Agen(i, j) = Tint(p) * Bsub(i, j);
    }
  }
  for (int j = 0; j < n; j++) {
    Agen(n, j) = v(j);
  }
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_zb;
#endif
  MyMatrix<Tint> Bp = GetZbasis(Agen);
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: neighbour GetZbasis|=" << time_zb << "\n";
#endif
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (Bp.rows() != n) {
    std::cerr << "GENUS: the neighbour basis has " << Bp.rows()
              << " rows instead of " << n << "\n";
    throw TerminalException{1};
  }
#endif
  // Gram of the neighbour: (Bp/p) G (Bp/p)^T = Bp G Bp^T / p^2
  MyMatrix<Tint> Prod = Bp * G * Bp.transpose();
  Tint p2 = Tint(p) * Tint(p);
  MyMatrix<T> GramN(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
#ifdef SANITY_CHECK_GENUS_ENUMERATION
      if (Prod(i, j) % p2 != 0) {
        std::cerr << "GENUS: the neighbour Gram matrix is not integral at ("
                  << i << "," << j << ")\n";
        throw TerminalException{1};
      }
#endif
      GramN(i, j) = UniversalScalarConversion<T, Tint>(Prod(i, j) / p2);
    }
  }
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  for (int i = 0; i < n; i++) {
    T two(2);
    if (ResInt(GramN(i, i), two) != T(0)) {
      std::cerr << "GENUS: the neighbour lattice is not even at " << i << "\n";
      throw TerminalException{1};
    }
  }
  if (DeterminantMat(GramN) !=
      UniversalScalarConversion<T, Tint>(DeterminantMat(G))) {
    std::cerr << "GENUS: the neighbour has a different determinant\n";
    throw TerminalException{1};
  }
#endif
  /*
    Reduce before returning. The basis coming out of GetZbasis and divided by p
    is arbitrary and can be very skewed, and everything downstream -- the
    canonical form and the automorphism group -- starts by enumerating short
    vectors, which is where a skewed basis is punishing. Without this the
    enumeration stalls on the first neighbour.
   */
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_lll;
#endif
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(GramN, os);
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: neighbour LLL|=" << time_lll << "\n";
#endif
  return rec.GramMatRed;
}

// The scheme by which classes are told apart and their automorphism groups
// computed.
//
//   CanonicalForm      -- a canonical Gram matrix per neighbour, recognising
//                         a known class by a dictionary lookup; automorphism
//                         groups by the Vinberg decomposition.
//   PleskenSouvignier  -- scheme A of the isomorphism strategy of Hecke.jl:
//                         neighbours bucketed by a cheap isometry invariant
//                         and, within a bucket, told apart by the
//                         Plesken-Souvignier isometry test on the Z-spanning
//                         family; automorphism groups by the same engine.
//   MatrixGroup        -- scheme B: same bucketing, but the isometry and the
//                         automorphism group come from the MatrixGroup
//                         algorithm on the full-rank (not Z-spanning)
//                         iterated-shortest family.
enum class GenusScheme { CanonicalForm, PleskenSouvignier, MatrixGroup };

/*
  The automorphism data of one lattice: the order, needed for the mass, and
  matrix generators, needed to cut the projective points into orbits. Both are
  derived from a single invariant vector family, since extracting it is the
  expensive part and computing them separately would do it twice.
 */
// The order is a count, so it is of the type the group counts with.
template <typename T, typename Tint, typename Tgroup> struct LatticeAutInfo {
  typename Tgroup::Tint order;
  std::vector<MyMatrix<Tint>> ListGenMat;
};

template <typename T, typename Tint, typename Tgroup>
LatticeAutInfo<T, Tint, Tgroup>
GetLatticeAutInfo(MyMatrix<T> const &GramMat, GenusScheme scheme,
                  std::ostream &os) {
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time;
#endif
  LatticeAutInfo<T, Tint, Tgroup> info;
  if (scheme == GenusScheme::CanonicalForm) {
    /*
      The automorphism group through the Vinberg decomposition Aut(L) =
      W(R) rtimes Aut(L, rho): the Weyl group W(R) of the root system is not
      enumerated, its order is read off the Dynkin type and its generators
      are the reflections in the simple roots; only the small residual
      Aut(L, rho) is searched. The order and the generators are exactly
      what the enumeration needs -- the order for the mass, the generators
      for the orbits of projective points.
     */
    VinbergAutom<Tint, Tgroup> vin =
        ComputeAutomorphismVinberg<T, Tint, Tgroup>(GramMat, os);
    info.order = vin.order;
    info.ListGenMat = vin.ListGen;
#ifdef DEBUG_GENUS_ENUMERATION
    os << "GENUS: |Aut|=" << vin.order << " = |W(R)|=" << vin.weyl_order
       << " x residual=" << vin.residual_order << ", " << info.ListGenMat.size()
       << " generators\n";
#endif
  } else if (scheme == GenusScheme::PleskenSouvignier) {
    // Scheme A: the Plesken-Souvignier engine on the Z-spanning family. The
    // order is the product of the stabilizer-chain orbit lengths.
    std::vector<MyMatrix<T>> ListMat{GramMat};
    PleskenSouvignierAutomResult<Tint, Tgroup> res =
        PleskenSouvignierLatticeAutomorphism<T, Tint, Tgroup>(ListMat, os);
    info.order = res.order;
    info.ListGenMat = res.ListGen;
  } else {
    // Scheme B: the MatrixGroup algorithm on the full-rank family. The order
    // is read off the permutation action of the generators on an invariant
    // family.
    info.ListGenMat =
        IteratedShortestAutomorphismGroup<T, Tint, Tgroup>(GramMat, os);
    MyMatrix<Tint> fam = root_iterated_shortest<T, Tint>(GramMat, os);
    info.order = OrderFromGens<Tint, Tgroup>(fam, info.ListGenMat);
  }
#ifdef DEBUG_GENUS_ENUMERATION
  if (scheme != GenusScheme::CanonicalForm) {
    os << "GENUS: |Aut|=" << info.order << ", " << info.ListGenMat.size()
       << " generators\n";
  }
#endif
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: GetLatticeAutInfo|=" << time << "\n";
#endif
  return info;
}

/*
  The canonical Gram matrix, used as the key of the dictionary of classes.

  ComputeCanonicalForm chooses the invariant vector family, which is what
  governs the cost here: the weight matrix is quadratic in the size of the
  family and the canonical labelling worse than quadratic, so a family that
  is one shell too large is fatal. On TestData/SlowCanonic/slow_canonic_2 the
  shells give 16942 vectors, which is a graph with 3.8 billion edges that
  nauty cannot allocate, against 238 for V_cv.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<T> GenusCanonicalGram(MyMatrix<T> const &GramMat, std::ostream &os) {
  // The Weyl-chamber canonical form (Brandhorst-Hofmann, Prop. 7.6): the
  // characteristic vectors in one Weyl chamber and the simple roots, a
  // family spanning the lattice, so neither the Weyl group nor the glue
  // step enters the computation. About 2x faster than the root-decomposed
  // form on the 813 classes of the determinant-351 genus. Every class of
  // the enumeration uses this same function, so the deduplication stays
  // consistent; the keys are not comparable with those of another method.
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_can;
#endif
  MyMatrix<Tint> B =
      ComputeCanonicalFormWeylChamber<T, Tint, Tgroup>(GramMat, os);
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: GenusCanonicalGram|=" << time_can << "\n";
#endif
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  return B_T * GramMat * B_T.transpose();
}

/*
  The projective points to visit, one per orbit of the automorphism group,
  as INDICES in the order of ProjectiveLineOfIndex. Without the orbits the
  rank-14 enumeration would canonicalise 16383 neighbours per lattice at
  p = 2; the orbits cut that by the order of the group acting.

  Nothing of the size of the projective space is materialised: the images
  under the generators are computed point by point through the index
  bijection, and the visited set is one bit per point (190 MB at p = 5,
  n = 14, where the point count is 1.5 x 10^9 and the earlier
  materialisation of the lines, the index map and one permutation table
  per generator was the memory explosion of 2026-10-25).

  Generators acting trivially on the projective space -- the scalar
  matrices mod p, in particular -Id, which is a generator of every
  automorphism group here -- are discarded first. When nothing remains the
  group acts trivially, every point is its own representative, and the
  selection says so through all = true rather than materialising
  0, ..., n_line - 1.
 */
struct ProjectiveLineSelection {
  // true: every index in [0, n_line) is a representative; ListRepr empty.
  bool all;
  size_t n_line;
  std::vector<size_t> ListRepr;
  size_t size() const { return all ? n_line : ListRepr.size(); }
  size_t at(size_t i) const { return all ? i : ListRepr[i]; }
};

template <typename Tint>
bool IsScalarModP(MyMatrix<Tint> const &eGen, int p) {
  int n = eGen.rows();
  int diag = ResidueModP(eGen(0, 0), p);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int res = ResidueModP(eGen(i, j), p);
      if (i == j && res != diag) {
        return false;
      }
      if (i != j && res != 0) {
        return false;
      }
    }
  }
  return true;
}

template <typename Tint>
ProjectiveLineSelection
ProjectiveLineOrbitRepresentatives(int n, int p,
                                   std::vector<MyMatrix<Tint>> const &ListGen) {
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  CheckProjectiveLineIndexing(n, p);
#endif
  size_t n_line = NumberProjectiveLines(n, p);
  // the generators as reduced integer matrices mod p, scalars dropped
  std::vector<std::vector<int>> ListGenRed;
  for (auto &eGen : ListGen) {
    if (!IsScalarModP(eGen, p)) {
      std::vector<int> gred(n * n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          gred[i * n + j] = ResidueModP(eGen(i, j), p);
        }
      }
      ListGenRed.push_back(gred);
    }
  }
  if (ListGenRed.empty()) {
    return {true, n_line, {}};
  }
  auto image = [&](std::vector<int> const &v,
                   std::vector<int> const &gred) -> size_t {
    std::vector<int> w(n, 0);
    for (int j = 0; j < n; j++) {
      int64_t sum = 0;
      for (int i = 0; i < n; i++) {
        sum += static_cast<int64_t>(v[i]) * gred[i * n + j];
      }
      w[j] = static_cast<int>(sum % p);
    }
    return IndexOfProjectiveLine(w, p);
  };
  Face status(n_line);
  std::vector<size_t> ListRepr;
  std::vector<size_t> ListActive;
  for (size_t i_line = 0; i_line < n_line; i_line++) {
    if (status[i_line]) {
      continue;
    }
    ListRepr.push_back(i_line);
    status[i_line] = 1;
    ListActive.push_back(i_line);
    while (!ListActive.empty()) {
      size_t pos = ListActive.back();
      ListActive.pop_back();
      std::vector<int> v = ProjectiveLineOfIndex(n, p, pos);
      for (auto &gred : ListGenRed) {
        size_t img = image(v, gred);
        if (!status[img]) {
          status[img] = 1;
          ListActive.push_back(img);
        }
      }
    }
  }
  return {false, n_line, std::move(ListRepr)};
}

/*
  The smallest prime not dividing the determinant. Taking 2 whenever the
  determinant is odd matters a great deal: the number of projective points is
  (p^n - 1)/(p - 1), so in rank 16 the choice p = 2 gives 65535 of them while
  p = 3 gives 21523360, and in rank 14 it is 16383 against 2391484. The p = 2
  case needs the valuation correction in GenusNeighbor, which is exactly why
  that correction is there.

  Connectivity of the p-neighbour graph on a spinor genus is the reason for
  excluding the primes dividing the determinant. It is not relied upon blindly:
  if the graph fails to be connected the accumulated mass simply stops short of
  the target and the program reports the enumeration as incomplete.
 */
template <typename T> int ChooseNeighborPrime(T const &det) {
  std::vector<int> ListPrime{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (auto &p : ListPrime) {
    T p_T(p);
    // Divisibility, not an exact division: over a ring "det / p" truncates
    // and would say that every prime divides the determinant. ResInt is
    // defined over both.
    if (ResInt(det, p_T) != T(0)) {
      return p;
    }
  }
  std::cerr << "GENUS: no suitable neighbouring prime found for det=" << det
            << "\n";
  throw TerminalException{1};
}

/*
  A cheap isometry invariant: the multiset of norms of the iterated-shortest
  vector family. It is an invariant of the lattice and, unlike a theta
  series, cheap to build even on the root-rich classes where the norm-ball
  explodes. Two lattices with different invariants are certainly not
  isometric, so only the equal-invariant ones need the isometry test.
 */
template <typename T, typename Tint>
std::map<T, int> GenusIsoInvariant(MyMatrix<T> const &GramMat,
                                   std::ostream &os) {
  MyMatrix<Tint> SHV = IteratedShortestVectorFamilyHalf<T, Tint>(GramMat, os);
  std::map<T, int> inv;
  for (int i = 0; i < SHV.rows(); i++) {
    MyVector<Tint> v = GetMatrixRow(SHV, i);
    inv[EvaluationQuadForm<T, Tint>(GramMat, v)]++;
  }
  return inv;
}

/*
  The enumeration itself. Starts from the seed lattices, forms p-neighbours,
  keeps the new ones (recognised either by canonical form or by the
  isomorphism strategy), and stops when the accumulated mass reaches the
  target.
 */
template <typename T, typename Tint, typename Tgroup>
GenusEnumerationResult<T, Tgroup>
GenusEnumeration(std::vector<MyMatrix<T>> const &ListSeed,
                 typename GenusEnumerationResult<T, Tgroup>::Tmass const
                     &TotalMass,
                 int prime, GenusScheme scheme, NeighborStrategy strategy,
                 std::ostream &os) {
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  using TintGroup = typename Tgroup::Tint;
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_total;
#endif
  int n = ListSeed[0].rows();
  GenusEnumerationResult<T, Tgroup> result;
  result.target_mass = TotalMass;
  result.accumulated_mass = Tmass(0);
  result.complete = false;
  // The classes are held in canonical form, so recognising an already known
  // one is a dictionary lookup rather than a comparison against every class
  // found so far. That is the whole point of canonicalising: the pairwise
  // scheme is quadratic in the class number.
  // CanonicalForm strategy: classes held by canonical Gram, so a known one
  // is a dictionary lookup. Isomorphism strategy: classes bucketed by the
  // cheap invariant, a bucket holding the indices of the classes with that
  // invariant.
  std::map<MyMatrix<T>, size_t> MapCanonic;
  std::map<std::map<T, int>, std::vector<size_t>> MapInvariant;
  std::vector<MyMatrix<T>> ListGramWork;
  // The automorphism data is kept alongside: it is the expensive part, it is
  // needed twice (the order for the mass, the generators for the orbits of
  // projective points), and recomputing it in the neighbour loop doubles the
  // cost of the whole enumeration. The generators also prune the isometry
  // test of the isomorphism strategy.
  std::vector<LatticeAutInfo<T, Tint, Tgroup>> ListAutInfo;
  // Diagnostic counters for the discriminating power of the cheap isometry
  // invariant: how many isometry tests are run, and how many return a match.
  // tests_per_positive close to 1 means the invariant is almost perfect (a
  // duplicate is found with a single test); larger means the buckets carry
  // non-isometric collisions that cost extra tests.
  [[maybe_unused]] size_t n_iso_test = 0;
  [[maybe_unused]] size_t n_iso_positive = 0;

  // Records a new class and its automorphism data; returns its index.
  auto f_add = [&](MyMatrix<T> const &StoreGram,
                   LatticeAutInfo<T, Tint, Tgroup> const &info) -> size_t {
    size_t idx = result.ListGram.size();
    result.ListGram.push_back(StoreGram);
    result.ListAutOrder.push_back(info.order);
    result.accumulated_mass +=
        Tmass(1) / UniversalScalarConversion<Tmass, TintGroup>(info.order);
    ListGramWork.push_back(StoreGram);
    ListAutInfo.push_back(info);
#ifdef DEBUG_GENUS_ENUMERATION
    os << "GENUS: class " << result.ListGram.size() << ", |Aut|=" << info.order
       << ", mass " << result.accumulated_mass << " / " << TotalMass << "\n";
#endif
    return idx;
  };

  auto f_insert = [&](MyMatrix<T> const &GramMat) -> bool {
    if (scheme == GenusScheme::CanonicalForm) {
      MyMatrix<T> GramCan = GenusCanonicalGram<T, Tint, Tgroup>(GramMat, os);
      if (MapCanonic.count(GramCan) == 1) {
        return false;
      }
      LatticeAutInfo<T, Tint, Tgroup> info =
          GetLatticeAutInfo<T, Tint, Tgroup>(GramCan, scheme, os);
      MapCanonic[GramCan] = f_add(GramCan, info);
      return true;
    }
    // Isomorphism schemes: bucket by the cheap invariant, then run the
    // isometry test of the chosen engine against the same-invariant classes.
    std::map<T, int> inv = GenusIsoInvariant<T, Tint>(GramMat, os);
    std::vector<size_t> &bucket = MapInvariant[inv];
    for (size_t idx : bucket) {
      bool isamorphic;
      if (scheme == GenusScheme::PleskenSouvignier) {
        // Scheme A: the Plesken-Souvignier isometry test on the Z-spanning
        // family, passing the known class's automorphism generators to prune
        // the search.
        std::vector<MyMatrix<T>> LM1{GramMat};
        std::vector<MyMatrix<T>> LM2{result.ListGram[idx]};
        isamorphic = PleskenSouvignierLatticeIsometry<T, Tint>(
                         LM1, LM2, ListAutInfo[idx].ListGenMat, os)
                         .has_value();
      } else {
        // Scheme B: the MatrixGroup isometry test on the full-rank family.
        isamorphic = IteratedShortestEquivalence<T, Tint, Tgroup>(
                         GramMat, result.ListGram[idx], os)
                         .has_value();
      }
      n_iso_test++;
      if (isamorphic) {
        n_iso_positive++;
#ifdef DEBUG_GENUS_ISO_STATS
        if (n_iso_positive % 100 == 0) {
          os << "GENUS_ISO: tests=" << n_iso_test
             << " positive=" << n_iso_positive << " ratio="
             << (double(n_iso_test) / double(n_iso_positive)) << "\n";
        }
#endif
        return false;
      }
    }
    LatticeAutInfo<T, Tint, Tgroup> info =
        GetLatticeAutInfo<T, Tint, Tgroup>(GramMat, scheme, os);
    bucket.push_back(f_add(GramMat, info));
    return true;
  };

  for (auto &eSeed : ListSeed) {
    if (result.accumulated_mass == TotalMass) {
      break;
    }
    (void)f_insert(eSeed);
  }
  // Expands one known class at one prime: the orbits of projective points,
  // the neighbour of each representative, insertion of the new classes.
  // Stops early when the mass closes.
  auto expand_at = [&](size_t pos, int p, size_t i_prime) -> void {
    MyMatrix<T> GramMat = ListGramWork[pos];
    std::vector<MyMatrix<Tint>> const &ListGenMat = ListAutInfo[pos].ListGenMat;
    MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(GramMat);
#ifdef TIMINGS_GENUS_ENUMERATION
    MicrosecondTime time_orb;
#endif
    ProjectiveLineSelection ListRepr =
        ProjectiveLineOrbitRepresentatives<Tint>(n, p, ListGenMat);
#ifdef TIMINGS_GENUS_ENUMERATION
    os << "|GENUS: ProjectiveLineOrbits|=" << time_orb << "\n";
#endif
#ifdef DEBUG_GENUS_ENUMERATION
    os << "GENUS: lattice " << (pos + 1) << " has " << ListRepr.size()
       << " orbits of projective points at p=" << p << "\n";
#endif
    // the representatives are indices; the line vectors are instantiated
    // one at a time, so nothing of the size of the projective space exists
    for (size_t i_orbit = 0; i_orbit < ListRepr.size(); i_orbit++) {
#ifdef DEBUG_GENUS_ENUMERATION
      os << "GENUS: lattice " << (pos + 1) << " orbit " << (i_orbit + 1)
         << " / " << ListRepr.size() << "\n";
#endif
      std::vector<int> eLine =
          ProjectiveLineOfIndex(n, p, ListRepr.at(i_orbit));
      std::optional<MyMatrix<T>> opt =
          GenusNeighbor<T, Tint>(G_int, eLine, p, os);
      if (opt) {
#ifdef DEBUG_GENUS_ENUMERATION
        os << "GENUS:   neighbour built, testing against known classes\n";
#endif
#ifdef TRACK_INFO_GENUS_ENUMERATION
        // With GENUS_NEIGHBOUR_DUMP set, the Gram matrix of every neighbour
        // is written to that file before it is canonicalized, each write
        // replacing the last. When a canonicalization dies (time, memory),
        // the file holds the Gram matrix that killed it, which is otherwise
        // unrecoverable: the neighbour is a function of the whole run state.
        static const char *dump_file = std::getenv("GENUS_NEIGHBOUR_DUMP");
        if (dump_file) {
          WriteMatrixFile(std::string(dump_file), *opt);
        }
#endif
        if (f_insert(*opt)) {
          result.ListNewClassPerPrime[i_prime]++;
        }
        if (result.accumulated_mass == TotalMass) {
          return;
        }
      }
    }
  };

  if (strategy == NeighborStrategy::OutPrime) {
    result.ListPrimeUsed.push_back(prime);
    result.ListNewClassPerPrime.push_back(0);
    size_t pos = 0;
    while (pos < ListGramWork.size() &&
           result.accumulated_mass != TotalMass) {
      expand_at(pos, prime, 0);
      pos++;
    }
  } else {
    /*
      The iterate strategy: see NeighborStrategy. ListPrimeUsed holds the
      primes in play in increasing order; next_expand[i] is the number of
      classes, in discovery order, already expanded at ListPrimeUsed[i].
      The invariant exits: the mass closes (complete), or the accumulated
      mass passes the target (impossible for genuine neighbours, so it
      would witness a programming error and stops the run as incomplete
      rather than spinning), or NumberProjectiveLines overflows at an
      absurd prime (an exception, equally a witness that something is
      wrong with the genus data).
     */
    std::vector<int> &ListPrimeUsed = result.ListPrimeUsed;
    std::vector<size_t> next_expand;
    auto next_prime = [](int q) -> int {
      int cand = q + 1;
      while (true) {
        bool is_prime = cand >= 2;
        for (int d = 2; d * d <= cand; d++) {
          if (cand % d == 0) {
            is_prime = false;
            break;
          }
        }
        if (is_prime) {
          return cand;
        }
        cand++;
      }
    };
    while (result.accumulated_mass != TotalMass) {
      if (result.accumulated_mass > TotalMass) {
        break;
      }
      size_t i_work = ListPrimeUsed.size();
      for (size_t i = 0; i < ListPrimeUsed.size(); i++) {
        if (next_expand[i] < ListGramWork.size()) {
          i_work = i;
          break;
        }
      }
      if (i_work == ListPrimeUsed.size()) {
        int p_new =
            ListPrimeUsed.empty() ? 2 : next_prime(ListPrimeUsed.back());
        ListPrimeUsed.push_back(p_new);
        result.ListNewClassPerPrime.push_back(0);
        next_expand.push_back(0);
#ifdef DEBUG_GENUS_ENUMERATION
        os << "GENUS: iterate strategy brings in prime " << p_new << "\n";
#endif
        continue;
      }
      expand_at(next_expand[i_work], ListPrimeUsed[i_work], i_work);
      next_expand[i_work]++;
    }
  }
  result.complete = (result.accumulated_mass == TotalMass);
#ifdef DEBUG_GENUS_ISO_STATS
  os << "GENUS_ISO: FINAL isometry_tests=" << n_iso_test
     << " positive=" << n_iso_positive << " tests_per_positive="
     << (n_iso_positive > 0 ? double(n_iso_test) / double(n_iso_positive) : 0.0)
     << "\n";
#endif
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: GenusEnumeration|=" << time_total << "\n";
#endif
  return result;
}

//
// Input and output
//

/*
  The genus description, as key/value lines:

      rank 14
      det 243
      prime 3        (optional, 0 or absent means "choose automatically")

  It carries no information the enumeration needs beyond the prime -- the
  p-neighbour construction preserves the genus by itself -- but it lets the
  seed lattices be checked against the genus they are supposed to represent,
  which catches a mismatched input immediately rather than as a mass that
  never closes.
 */
template <typename T>
GenusSpec<T> ReadGenusSpecFile(std::string const &file_name) {
  if (!FILE_IsExistingFile(file_name)) {
    std::cerr << "GENUS: the file " << file_name << " does not exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  GenusSpec<T> spec;
  spec.rank = -1;
  spec.det = T(0);
  spec.prime = 0;
  spec.strategy = NeighborStrategy::OutPrime;
  std::string key;
  while (is >> key) {
    if (key == "rank") {
      is >> spec.rank;
    } else if (key == "det") {
      is >> spec.det;
    } else if (key == "prime") {
      is >> spec.prime;
    } else if (key == "strategy") {
      std::string strat;
      is >> strat;
      if (strat == "out_prime") {
        spec.strategy = NeighborStrategy::OutPrime;
      } else if (strat == "iterate") {
        spec.strategy = NeighborStrategy::Iterate;
      } else {
        std::cerr << "GENUS: unrecognised strategy \"" << strat << "\" in "
                  << file_name << ". Allowed: out_prime, iterate\n";
        throw TerminalException{1};
      }
    } else {
      std::cerr << "GENUS: unrecognised key \"" << key << "\" in "
                << file_name << ". Allowed: rank, det, prime, strategy\n";
      throw TerminalException{1};
    }
  }
  if (spec.rank <= 0) {
    std::cerr << "GENUS: the rank is missing or not positive in " << file_name
              << "\n";
    throw TerminalException{1};
  }
  if (spec.det <= T(0)) {
    std::cerr << "GENUS: the determinant is missing or not positive in "
              << file_name << "\n";
    throw TerminalException{1};
  }
  return spec;
}

// The mass, as "num den" or as a single integer. A fraction, so it is read
// over the field and not over the ring the Gram matrices live in.
template <typename Tmass> Tmass ReadMassFile(std::string const &file_name) {
  if (!FILE_IsExistingFile(file_name)) {
    std::cerr << "GENUS: the file " << file_name << " does not exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  Tmass num, den;
  if (!(is >> num)) {
    std::cerr << "GENUS: failed to read the mass from " << file_name << "\n";
    throw TerminalException{1};
  }
  if (!(is >> den)) {
    den = Tmass(1);
  }
  if (den == Tmass(0)) {
    std::cerr << "GENUS: the denominator of the mass is zero in " << file_name
              << "\n";
    throw TerminalException{1};
  }
  return num / den;
}

template <typename T, typename Tgroup>
void WriteGenusEnumerationResult(
    std::ostream &os, GenusEnumerationResult<T, Tgroup> const &result,
    std::string const &OutFormat, int prime) {
  if (OutFormat == "Summary") {
    os << "prime = " << prime << "\n";
    os << "primes used (new classes found):";
    for (size_t i = 0; i < result.ListPrimeUsed.size(); i++) {
      os << " " << result.ListPrimeUsed[i] << " ("
         << result.ListNewClassPerPrime[i] << ")";
    }
    os << "\n";
    os << "class number = " << result.ListGram.size() << "\n";
    os << "accumulated mass = " << result.accumulated_mass << "\n";
    os << "target mass = " << result.target_mass << "\n";
    os << "complete = " << (result.complete ? "true" : "false") << "\n";
    if (!result.complete) {
      os << "The accumulated mass is short of the target. Either the genus "
         << "has several spinor genera and one seed per spinor genus is "
         << "needed, or the prime is not suitable.\n";
    }
    return;
  }
  if (OutFormat == "CPP") {
    os << result.ListGram.size() << "\n";
    for (auto &eGram : result.ListGram) {
      WriteMatrix(os, eGram);
    }
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(ListGram:=[";
    for (size_t i = 0; i < result.ListGram.size(); i++) {
      if (i > 0) {
        os << ",\n";
      }
      WriteMatrixGAP(os, result.ListGram[i]);
    }
    os << "],\n ListAutOrder:=[";
    for (size_t i = 0; i < result.ListAutOrder.size(); i++) {
      if (i > 0) {
        os << ", ";
      }
      os << result.ListAutOrder[i];
    }
    os << "],\n accumulated_mass:=" << result.accumulated_mass
       << ",\n target_mass:=" << result.target_mass
       << ",\n complete:=" << (result.complete ? "true" : "false")
       << ",\n prime:=" << prime << ",\n ListPrimeUsed:=[";
    for (size_t i = 0; i < result.ListPrimeUsed.size(); i++) {
      os << (i > 0 ? ", " : "") << result.ListPrimeUsed[i];
    }
    os << "],\n ListNewClassPerPrime:=[";
    for (size_t i = 0; i < result.ListNewClassPerPrime.size(); i++) {
      os << (i > 0 ? ", " : "") << result.ListNewClassPerPrime[i];
    }
    os << "]);\n";
    return;
  }
  std::cerr << "GENUS: failed to find a matching entry for OutFormat="
            << OutFormat << "\n";
  throw TerminalException{1};
}

#endif  //  SRC_GENERA_GENUSENUMERATION_H_
