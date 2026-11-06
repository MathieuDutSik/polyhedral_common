// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSENUMERATION_H_
#define SRC_GENERA_GENUSENUMERATION_H_

// clang-format off
#include "Boost_bitset.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "LatticeStabEquiCan.h"
#include "LatticeMethods.h"
#include "LatticeRootDecomposition.h"
#include "LatticeAutomorphismVinberg.h"
#include "LatticePleskenSouvignier.h"
#include "InvariantVectorFamily.h"
#include "ClassicLLL.h"
#include "Positivity.h"
#include "Namelist.h"
#include "Basic_random.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <random>
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

  RandomWalk: most neighbours of the full expansion are classes already
  known, so each new class is first treated by random 2-neighbours. A class
  with |Aut| = g gets the number b, the smallest integer with g < 2^b, and
  the bound max(1, ceil(hits_factor * b)); random projective lines mod 2 are
  drawn (uniformly, so an orbit is hit in proportion to its size and no
  orbit computation is needed) and their neighbours inserted, until that
  many of them are known classes. A new class never stops the walk. A line
  need not give a neighbour, and a class may have no 2-neighbour at all, so
  the draws that give none are bounded by failed_draws_factor * b.
  The walk is at p = 2 by default. With walk_prime_max > 2, every draw takes
  instead a prime uniformly among those of [walk_prime_min, walk_prime_max]
  not dividing 2 det, and an isotropic line at that prime: a neighbour costs
  the same at any prime, at such a prime the neighbour graph is connected on
  the spinor genus (Kneser), which p = 2 is not, and the neighbours at a
  large prime are spread much more widely over the genus.
  The pending classes are treated by decreasing |Aut|. If the mass closes in
  this first phase the enumeration is over. Otherwise the full expansion
  follows, at 2, then 3, 5, ... as in Iterate, and every new class it finds
  is again treated by the random walk first. In both phases the class with
  the largest |Aut| goes first. Deterministic for a given seed, and
  terminating since the second phase is Iterate.
 */
enum class NeighborStrategy { OutPrime, Iterate, RandomWalk };

/*
  The order of the full expansions at a prime. Default (order aut in the
  namelist): the one of the schedule, discovery order for OutPrime and
  Iterate, largest |Aut| first for RandomWalk. Deficit, the default of the
  namelist: the class with the largest Hecke deficit, the number
  of its points giving a neighbour that the classes already expanded do not
  account for (see GenusEnumerationState), ties going to the larger |Aut|.
  When the mass has not closed, the missing classes are neighbours of the
  classes with a positive deficit only, so those go first.
 */
enum class GenusExpansionOrder { Default, Deficit };

// The parameters of the random walk phase of the RandomWalk strategy.
struct GenusRandomWalkOptions {
  uint64_t seed = 1;
  double hits_factor = 0.25;
  size_t failed_draws_factor = 16;
  int walk_prime_min = 2;
  int walk_prime_max = 2;
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
  // The cost in neighbours: how many the full expansion built, and how many
  // the random walk built and how many new classes those gave. Whether the
  // mass closed during the first, random walk only, phase.
  size_t n_neighbor_full = 0;
  // The expansions started, and those the Hecke bookkeeping skipped.
  size_t n_expansion = 0;
  size_t n_expansion_skipped = 0;
  size_t n_neighbor_random = 0;
  size_t n_new_random = 0;
  bool closed_in_random_walk = false;
  // Whether the enumeration was stopped by max_runtime_second.
  bool stopped_by_runtime = false;
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
/*
  The neighbour before reduction: Bp, whose rows divided by p are a basis of
  N in the coordinates of L, and the Gram matrix of that basis. Split from
  GenusNeighbor so that the enumeration of lattices with an isometry
  (GenusEnumerationGroup.h) can carry the isometry over to N; GenusNeighbor
  does the same work as before.
 */
template <typename T, typename Tint> struct GenusNeighborBasis {
  MyMatrix<Tint> Bp;
  MyMatrix<T> GramN;
};


// w = G x, the linear form y -> (y, x).
template <typename Tint>
MyVector<Tint> GenusLinearForm(MyMatrix<Tint> const &G, MyVector<Tint> const &x) {
  int n = G.rows();
  MyVector<Tint> w(n);
  for (int i = 0; i < n; i++) {
    Tint sum(0);
    for (int j = 0; j < n; j++) {
      sum += G(i, j) * x(j);
    }
    w(i) = sum;
  }
  return w;
}

template <typename Tint>
Tint GenusNormOf(MyVector<Tint> const &x, MyVector<Tint> const &wx) {
  Tint sum(0);
  for (int i = 0; i < x.size(); i++) {
    sum += x(i) * wx(i);
  }
  return sum;
}

// What the neighbour construction needs of the point v: the lift, its linear
// form, a pivot of that form mod p, and the norm with its valuation.
template <typename Tint> struct GenusNeighborLine {
  MyVector<Tint> v;
  MyVector<Tint> w;
  int jpiv;
  Tint nrm;
  int val;
};

/*
  The point v of P^{n-1}(F_p) gives a neighbour exactly when it is not in the
  radical of the form mod p and v_p((v,v)) > e, e = v_p(2): the valuation
  e + 1 is corrected by GenusNeighbor, a smaller one is not. This test is
  cheap, and it is what the Hecke bookkeeping counts the lines with.
 */
template <typename Tint>
std::optional<GenusNeighborLine<Tint>>
GenusNeighborLineData(MyMatrix<Tint> const &G, std::vector<int> const &v_line,
                      int p) {
  int n = G.rows();
  MyVector<Tint> v(n);
  for (int i = 0; i < n; i++) {
    v(i) = Tint(v_line[i]);
  }
  MyVector<Tint> w = GenusLinearForm(G, v);
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
  int e = (p == 2) ? 1 : 0;
  Tint nrm = GenusNormOf(v, w);
  int val = ValuationP(nrm, p);
  if (val <= e) {
    return {};
  }
  return GenusNeighborLine<Tint>{v, w, jpiv, nrm, val};
}

template <typename T, typename Tint>
std::optional<GenusNeighborBasis<T, Tint>>
GenusNeighborKernel(MyMatrix<Tint> const &G, std::vector<int> const &v_line,
                    int p, [[maybe_unused]] std::ostream &os) {
  int n = G.rows();
  std::optional<GenusNeighborLine<Tint>> opt_line =
      GenusNeighborLineData(G, v_line, p);
  if (!opt_line) {
    return {};
  }
  MyVector<Tint> v = opt_line->v;
  MyVector<Tint> w = opt_line->w;
  int jpiv = opt_line->jpiv;
  Tint nrm = opt_line->nrm;
  int val = opt_line->val;
  int e = (p == 2) ? 1 : 0;
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
    w = GenusLinearForm(G, v);
    nrm = GenusNormOf(v, w);
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
  return GenusNeighborBasis<T, Tint>{std::move(Bp), std::move(GramN)};
}

template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GenusNeighbor(MyMatrix<Tint> const &G, std::vector<int> const &v_line, int p,
              [[maybe_unused]] std::ostream &os) {
  std::optional<GenusNeighborBasis<T, Tint>> opt =
      GenusNeighborKernel<T, Tint>(G, v_line, p, os);
  if (!opt) {
    return {};
  }
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
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(opt->GramN, os);
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: neighbour LLL|=" << time_lll << "\n";
#endif
  return rec.GramMatRed;
}

/*
  How the classes are told apart and their automorphism groups computed, by
  the methods of src_latt/LatticeMethods.h.

  by_canonical true: a canonical Gram matrix per neighbour (canonical), and a
  known class is a dictionary lookup. by_canonical false: the neighbours are
  bucketed by a cheap isometry invariant and, within a bucket, told apart by
  the isometry test (isomorphism), as in the isomorphism strategy of
  Hecke.jl; the Plesken-Souvignier test is pruned by the automorphism group
  of the known class. The automorphism groups (automorphism) give the orders
  for the mass and the generators for the orbits of projective points.

  The defaults are the Weyl-chamber methods (Brandhorst-Hofmann, Prop. 7.6)
  for both: the characteristic vectors in one Weyl chamber and the simple
  roots, a family spanning the lattice. On the classes of the genera of
  determinant 30 (rank 19) and 210 (rank 15), the canonical form takes 46 ms
  against 94 for weyl, 249 for vinberg and 4.8 s for plesken_souvignier, and
  the automorphism group 38 and 16 ms against 73 and 35 for vinberg, the
  next fastest.
 */
struct GenusMethods {
  bool by_canonical = true;
  LatticeStabCanIsomMethod canonical = LatticeStabCanIsomMethod::WeylChamber;
  LatticeStabCanIsomMethod isomorphism =
      LatticeStabCanIsomMethod::PleskenSouvignier;
  LatticeStabCanIsomMethod automorphism = LatticeStabCanIsomMethod::WeylChamber;
};

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
GetLatticeAutInfo(MyMatrix<T> const &GramMat, GenusMethods const &methods,
                  std::ostream &os) {
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time;
#endif
  std::vector<MyMatrix<T>> ListMat{GramMat};
  IntAutomorphismGroup<Tint, Tgroup> aut =
      LatticeAutomorphismGroup<T, Tint, Tgroup>(ListMat, methods.automorphism,
                                                os);
  LatticeAutInfo<T, Tint, Tgroup> info;
  info.order = aut.order;
  info.ListGenMat = std::move(aut.ListGen);
#ifdef DEBUG_GENUS_ENUMERATION
  os << "GENUS: |Aut|=" << info.order << ", " << info.ListGenMat.size()
     << " generators\n";
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
MyMatrix<T> GenusCanonicalGram(MyMatrix<T> const &GramMat,
                               LatticeStabCanIsomMethod method,
                               std::ostream &os) {
  // Every class of the enumeration uses this same method, so the
  // deduplication stays consistent; the keys are not comparable with those
  // of another method.
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_can;
#endif
  MyMatrix<Tint> B = LatticeCanonicalBasis<T, Tint, Tgroup>(GramMat, method, os);
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
  // the number of points of the orbit of ListRepr[i]
  std::vector<size_t> ListOrbitSize;
  size_t size() const { return all ? n_line : ListRepr.size(); }
  size_t at(size_t i) const { return all ? i : ListRepr[i]; }
  size_t orbit_size(size_t i) const { return all ? 1 : ListOrbitSize[i]; }
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
    return {true, n_line, {}, {}};
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
  std::vector<size_t> ListOrbitSize;
  std::vector<size_t> ListActive;
  for (size_t i_line = 0; i_line < n_line; i_line++) {
    if (status[i_line]) {
      continue;
    }
    ListRepr.push_back(i_line);
    status[i_line] = 1;
    ListActive.push_back(i_line);
    size_t orbit_size = 0;
    while (!ListActive.empty()) {
      size_t pos = ListActive.back();
      ListActive.pop_back();
      orbit_size++;
      std::vector<int> v = ProjectiveLineOfIndex(n, p, pos);
      for (auto &gred : ListGenRed) {
        size_t img = image(v, gred);
        if (!status[img]) {
          status[img] = 1;
          ListActive.push_back(img);
        }
      }
    }
    ListOrbitSize.push_back(orbit_size);
  }
  return {false, n_line, std::move(ListRepr), std::move(ListOrbitSize)};
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
  The state of an enumeration, shared by the two neighbour schedules: the
  classes found so far with their automorphism data, the structure that
  recognises a known class, and the accumulated mass. The two operations on
  it, inserting a lattice and expanding a known class at one prime, are the
  same for every schedule; what differs is only the order of the expansions,
  see GenusEnumerationOutPrime and GenusEnumerationIterate.
 */
template <typename T, typename Tint, typename Tgroup>
struct GenusEnumerationState {
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  using TintGroup = typename Tgroup::Tint;
  int n;
  GenusMethods methods;
  Tmass const &TotalMass;
  std::ostream &os;
  GenusEnumerationResult<T, Tgroup> result;
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
  size_t n_iso_test = 0;
  size_t n_iso_positive = 0;
  /*
    The Hecke bookkeeping. At a prime p, write T(M, L) for the number of
    points of P^{n-1}(F_p) whose neighbour of M is isometric to L. Counting
    the neighbouring pairs up to isometry gives

        T(M, L) / |Aut M| = T(L, M) / |Aut L| ,

    so the expansions already done tell part of the neighbours of every other
    class: hecke_in[i][L] is the sum of T(M, L) / |Aut M| over the classes
    M != L expanded at the i-th prime. When |Aut L| hecke_in[i][L] is the
    number of points of L that give a neighbour, every neighbour of L is a
    class already expanded, and the expansion of L can be skipped: it would
    find nothing new, and its neighbours towards the classes expanded later
    are zero by the same identity, so the counts stay exact. hecke_row[i][M]
    is the row T(M, .) of an expanded class, kept for the checks only.
   */
  bool hecke;
  GenusExpansionOrder order;
  std::vector<std::vector<Tmass>> hecke_in;
  // hecke_npoint[i]: the number of points giving a neighbour at the i-th
  // prime, -1 while not computed. Whether a point gives a neighbour depends
  // on L / p^2 L only, so on L (x) Z_p, which is the same for all the
  // classes of the genus: one number per prime.
  std::vector<int64_t> hecke_npoint;
#ifdef SANITY_CHECK_GENUS_ENUMERATION_EXTENSIVE
  std::vector<std::map<size_t, std::map<size_t, size_t>>> hecke_row;
#endif

  // The deficit order needs the bookkeeping, so it switches it on.
  // The time limit: 0 for none. When it is reached the enumeration stops as
  // if it were over, and the classes found so far are the result, with
  // complete false.
  int max_runtime_second = 0;
  std::chrono::steady_clock::time_point time_start =
      std::chrono::steady_clock::now();

  GenusEnumerationState(int _n, GenusMethods const &_methods,
                        Tmass const &_TotalMass,
                        bool _hecke, GenusExpansionOrder _order,
                        std::ostream &_os)
      : n(_n), methods(_methods), TotalMass(_TotalMass), os(_os),
        hecke(_hecke || _order == GenusExpansionOrder::Deficit),
        order(_order) {
    result.target_mass = TotalMass;
    result.accumulated_mass = Tmass(0);
    result.complete = false;
  }

  bool is_mass_closed() const { return result.accumulated_mass == TotalMass; }

  // The mass is closed, or the time limit is reached.
  bool is_finished() {
    if (is_mass_closed()) {
      return true;
    }
    if (max_runtime_second > 0 && !result.stopped_by_runtime) {
      auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                         std::chrono::steady_clock::now() - time_start)
                         .count();
      if (elapsed >= max_runtime_second) {
        result.stopped_by_runtime = true;
      }
    }
    return result.stopped_by_runtime;
  }

  // Records a new class and its automorphism data; returns its index.
  size_t add(MyMatrix<T> const &StoreGram,
             LatticeAutInfo<T, Tint, Tgroup> const &info) {
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
  }

  // Inserts the lattice if its class is not known yet; returns whether it
  // was new.
  bool insert(MyMatrix<T> const &GramMat) {
    return insert_index(GramMat).first;
  }

  // The same, with the index of the class of the lattice, new or known.
  std::pair<bool, size_t> insert_index(MyMatrix<T> const &GramMat) {
    if (methods.by_canonical) {
      MyMatrix<T> GramCan =
          GenusCanonicalGram<T, Tint, Tgroup>(GramMat, methods.canonical, os);
      // One search of the map: the entry is created here and its index
      // filled in once the class is added.
      auto [iter, is_new] = MapCanonic.try_emplace(GramCan, 0);
      if (!is_new) {
        return {false, iter->second};
      }
      LatticeAutInfo<T, Tint, Tgroup> info =
          GetLatticeAutInfo<T, Tint, Tgroup>(iter->first, methods, os);
      iter->second = add(iter->first, info);
      return {true, iter->second};
    }
    // By isometry: bucket by the cheap invariant, then run the isometry
    // test against the same-invariant classes, with the automorphism
    // generators of the known class.
    std::map<T, int> inv = GenusIsoInvariant<T, Tint>(GramMat, os);
    std::vector<size_t> &bucket = MapInvariant[inv];
    for (size_t idx : bucket) {
      std::vector<MyMatrix<T>> LM1{GramMat};
      std::vector<MyMatrix<T>> LM2{result.ListGram[idx]};
      bool isamorphic = LatticeIsometry<T, Tint, Tgroup>(
                            LM1, LM2, ListAutInfo[idx].ListGenMat,
                            methods.isomorphism, os)
                            .has_value();
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
        return {false, idx};
      }
    }
    LatticeAutInfo<T, Tint, Tgroup> info =
        GetLatticeAutInfo<T, Tint, Tgroup>(GramMat, methods, os);
    size_t idx_new = add(GramMat, info);
    bucket.push_back(idx_new);
    return {true, idx_new};
  }

  // The number of points of the class pos giving a neighbour at p, from the
  // orbits ListRepr of its automorphism group.
  size_t count_points(int p, ProjectiveLineSelection const &ListRepr,
                      MyMatrix<Tint> const &G_int) const {
    size_t n_point = 0;
    for (size_t i_orbit = 0; i_orbit < ListRepr.size(); i_orbit++) {
      std::vector<int> eLine = ProjectiveLineOfIndex(n, p, ListRepr.at(i_orbit));
      if (GenusNeighborLineData(G_int, eLine, p)) {
        n_point += ListRepr.orbit_size(i_orbit);
      }
    }
    return n_point;
  }

  // The number of points giving a neighbour at the i-th prime p, the same
  // for all the classes; computed on the class pos the first time.
  size_t n_point_genus(size_t pos, int p, size_t i_prime) {
    if (hecke_npoint.size() <= i_prime) {
      hecke_npoint.resize(i_prime + 1, -1);
    }
    if (hecke_npoint[i_prime] < 0) {
      MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(ListGramWork[pos]);
      ProjectiveLineSelection ListRepr = ProjectiveLineOrbitRepresentatives<Tint>(
          n, p, ListAutInfo[pos].ListGenMat);
      hecke_npoint[i_prime] =
          static_cast<int64_t>(count_points(p, ListRepr, G_int));
    }
    return static_cast<size_t>(hecke_npoint[i_prime]);
  }

  // The Hecke deficit of the class pos at p: its points giving a neighbour
  // that the classes expanded at p do not account for.
  Tmass deficit(size_t pos, int p, size_t i_prime) {
    size_t n_point = n_point_genus(pos, p, i_prime);
    Tmass in(0);
    if (hecke_in.size() > i_prime && hecke_in[i_prime].size() > pos) {
      in = hecke_in[i_prime][pos];
    }
    Tmass order_mass =
        UniversalScalarConversion<Tmass, TintGroup>(result.ListAutOrder[pos]);
    return Tmass(n_point) - in * order_mass;
  }

  // The class to expand next among the unexpanded ones at the i-th prime p,
  // which are listed in increasing order: default_pick (a position in that
  // list, the choice of the schedule) with the default order, the one of
  // largest deficit with the deficit order.
  size_t pick(std::vector<size_t> const &unexpanded, int p, size_t i_prime,
              size_t default_pick) {
    if (order == GenusExpansionOrder::Default) {
      return unexpanded[default_pick];
    }
    size_t best = unexpanded[0];
    Tmass best_def = deficit(best, p, i_prime);
    for (size_t i = 1; i < unexpanded.size(); i++) {
      size_t pos = unexpanded[i];
      Tmass def = deficit(pos, p, i_prime);
      if (def > best_def ||
          (def == best_def &&
           result.ListAutOrder[pos] > result.ListAutOrder[best])) {
        best = pos;
        best_def = def;
      }
    }
    return best;
  }

  // Expands one known class at one prime: the orbits of projective points,
  // the neighbour of each representative, insertion of the new classes.
  // Stops early when the mass closes. i_prime is the index of p in
  // result.ListPrimeUsed, for the count of the new classes. f_new is called
  // on the index of every new class, as soon as it is inserted.
  void expand_at(size_t pos, int p, size_t i_prime) {
    expand_at(pos, p, i_prime, [](size_t) -> void {});
  }

  template <typename Fnew>
  void expand_at(size_t pos, int p, size_t i_prime, Fnew f_new) {
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
    result.n_expansion++;
    Tmass order_mass =
        UniversalScalarConversion<Tmass, TintGroup>(result.ListAutOrder[pos]);
    [[maybe_unused]] bool would_skip = false;
    if (hecke) {
      if (hecke_in.size() <= i_prime) {
        hecke_in.resize(i_prime + 1);
      }
      hecke_in[i_prime].resize(ListGramWork.size(), Tmass(0));
      // the points of the lattice that give a neighbour
      size_t n_point = n_point_genus(pos, p, i_prime);
#ifdef SANITY_CHECK_GENUS_ENUMERATION
      size_t n_point_here = count_points(p, ListRepr, G_int);
      if (n_point_here != n_point) {
        std::cerr << "GENUS: lattice " << (pos + 1) << " has " << n_point_here
                  << " points giving a neighbour at p=" << p
                  << " against " << n_point << " for the genus\n";
        throw TerminalException{1};
      }
#endif
      Tmass accounted = hecke_in[i_prime][pos] * order_mass;
#ifdef SANITY_CHECK_GENUS_ENUMERATION
      if (accounted > Tmass(n_point)) {
        std::cerr << "GENUS: the expanded classes account for " << accounted
                  << " neighbours of lattice " << (pos + 1) << " at p=" << p
                  << " which has only " << n_point << "\n";
        throw TerminalException{1};
      }
#endif
      if (accounted == Tmass(n_point)) {
#ifdef DEBUG_GENUS_ENUMERATION
        os << "GENUS: lattice " << (pos + 1) << " skipped at p=" << p
           << ", its " << n_point << " neighbours are all expanded classes\n";
#endif
#ifdef SANITY_CHECK_GENUS_ENUMERATION_EXTENSIVE
        // the expansion is done anyway, and checked below
        would_skip = true;
#else
        result.n_expansion_skipped++;
        return;
#endif
      }
    }
    // the row T(pos, .) of the expansion, for the bookkeeping
    std::map<size_t, size_t> row;
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
        result.n_neighbor_full++;
        std::pair<bool, size_t> pair = insert_index(*opt);
        if (hecke) {
          row[pair.second] += ListRepr.orbit_size(i_orbit);
        }
        if (pair.first) {
          result.ListNewClassPerPrime[i_prime]++;
          f_new(pair.second);
        }
        if (is_finished()) {
          return;
        }
      }
    }
    if (hecke) {
      std::vector<Tmass> &vin = hecke_in[i_prime];
      vin.resize(ListGramWork.size(), Tmass(0));
      for (auto &kv : row) {
        if (kv.first != pos) {
          vin[kv.first] += Tmass(kv.second) / order_mass;
        }
      }
#ifdef SANITY_CHECK_GENUS_ENUMERATION_EXTENSIVE
      std::cerr << "GENUS: SANITY_CHECK_GENUS_ENUMERATION_EXTENSIVE, the Hecke "
                   "identity for lattice "
                << (pos + 1) << " at p=" << p << "\n";
      if (hecke_row.size() <= i_prime) {
        hecke_row.resize(i_prime + 1);
      }
      auto &rows = hecke_row[i_prime];
      auto get = [](std::map<size_t, size_t> const &m, size_t k) -> size_t {
        auto iter = m.find(k);
        return iter == m.end() ? 0 : iter->second;
      };
      for (auto &kv : rows) {
        size_t other = kv.first;
        Tmass order_other = UniversalScalarConversion<Tmass, TintGroup>(
            result.ListAutOrder[other]);
        Tmass lhs = Tmass(get(row, other)) / order_mass;
        Tmass rhs = Tmass(get(kv.second, pos)) / order_other;
        if (lhs != rhs) {
          std::cerr << "GENUS: Hecke identity fails between lattices "
                    << (pos + 1) << " and " << (other + 1) << " at p=" << p
                    << ": T=" << get(row, other) << " |Aut|=" << order_mass
                    << " against T=" << get(kv.second, pos)
                    << " |Aut|=" << order_other << "\n";
          throw TerminalException{1};
        }
      }
      if (would_skip) {
        for (auto &kv : row) {
          if (kv.first != pos && rows.count(kv.first) == 0) {
            std::cerr << "GENUS: lattice " << (pos + 1) << " would have been "
                      << "skipped at p=" << p << " but has the unexpanded "
                      << "neighbour " << (kv.first + 1) << "\n";
            throw TerminalException{1};
          }
        }
        result.n_expansion_skipped++;
      }
      rows[pos] = row;
      std::cerr << "GENUS: SANITY_CHECK_GENUS_ENUMERATION_EXTENSIVE, done\n";
#endif
    }
  }
};

/*
  The out_prime schedule: see NeighborStrategy. Every known class, in
  discovery order (or by deficit, see GenusExpansionOrder), is expanded at
  the one prime until the mass closes.
 */
template <typename T, typename Tint, typename Tgroup>
void GenusEnumerationOutPrime(GenusEnumerationState<T, Tint, Tgroup> &state,
                              int prime) {
  auto &result = state.result;
  result.ListPrimeUsed.push_back(prime);
  result.ListNewClassPerPrime.push_back(0);
  std::vector<bool> is_expanded;
  while (!state.is_finished()) {
    size_t n_class = state.ListGramWork.size();
    is_expanded.resize(n_class, false);
    std::vector<size_t> unexpanded;
    for (size_t pos = 0; pos < n_class; pos++) {
      if (!is_expanded[pos]) {
        unexpanded.push_back(pos);
      }
    }
    if (unexpanded.empty()) {
      break;
    }
    size_t pos = state.pick(unexpanded, prime, 0, 0);
    is_expanded[pos] = true;
    state.expand_at(pos, prime, 0);
  }
}

// The smallest prime larger than q, for the schedules that bring in the
// primes one after the other.
inline int GenusNextPrime(int q) {
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
}

/*
  The iterate schedule: see NeighborStrategy. ListPrimeUsed holds the
  primes in play in increasing order; is_expanded[i][pos] tells whether the
  class pos has been expanded at ListPrimeUsed[i]. At a prime the classes go
  in discovery order (or by deficit, see GenusExpansionOrder).
  The invariant exits: the mass closes (complete), or the accumulated
  mass passes the target (impossible for genuine neighbours, so it
  would witness a programming error and stops the run as incomplete
  rather than spinning), or NumberProjectiveLines overflows at an
  absurd prime (an exception, equally a witness that something is
  wrong with the genus data).
 */
template <typename T, typename Tint, typename Tgroup>
void GenusEnumerationIterate(GenusEnumerationState<T, Tint, Tgroup> &state) {
  auto &result = state.result;
  std::vector<int> &ListPrimeUsed = result.ListPrimeUsed;
  std::vector<std::vector<bool>> is_expanded;
  while (!state.is_finished()) {
    if (result.accumulated_mass > state.TotalMass) {
      break;
    }
    size_t n_class = state.ListGramWork.size();
    size_t i_work = ListPrimeUsed.size();
    std::vector<size_t> unexpanded;
    for (size_t i = 0; i < ListPrimeUsed.size(); i++) {
      is_expanded[i].resize(n_class, false);
      for (size_t pos = 0; pos < n_class; pos++) {
        if (!is_expanded[i][pos]) {
          unexpanded.push_back(pos);
        }
      }
      if (!unexpanded.empty()) {
        i_work = i;
        break;
      }
    }
    if (i_work == ListPrimeUsed.size()) {
      int p_new = ListPrimeUsed.empty() ? 2 : GenusNextPrime(ListPrimeUsed.back());
      ListPrimeUsed.push_back(p_new);
      result.ListNewClassPerPrime.push_back(0);
      is_expanded.push_back(std::vector<bool>(n_class, false));
#ifdef DEBUG_GENUS_ENUMERATION
      state.os << "GENUS: iterate strategy brings in prime " << p_new << "\n";
#endif
      continue;
    }
    int p = ListPrimeUsed[i_work];
    size_t pos = state.pick(unexpanded, p, i_work, 0);
    is_expanded[i_work][pos] = true;
    state.expand_at(pos, p, i_work);
  }
}

// The smallest b with g < 2^b, for a class with |Aut| = g: the scale of the
// number of known classes the random walk must hit from that class before
// it counts as treated.
template <typename TintGroup>
size_t GenusRandomWalkLogOrder(TintGroup const &order) {
  size_t b = 0;
  TintGroup pow(1);
  while (pow <= order) {
    pow *= 2;
    b++;
  }
  return b;
}

/*
  The random_walk schedule: see NeighborStrategy. The classes waiting for
  their random walk are in pending, the full expansion is scheduled as in
  GenusEnumerationIterate except that, at a given prime, the unexpanded class
  with the largest |Aut| goes first. A class found by either phase is put in
  pending and the random walks are run before the full expansion goes on.
 */
template <typename T, typename Tint, typename Tgroup>
void GenusEnumerationRandomWalk(GenusEnumerationState<T, Tint, Tgroup> &state,
                                GenusRandomWalkOptions const &rw) {
  auto &result = state.result;
  int n = state.n;
  std::mt19937_64 engine(rw.seed);
  // The primes of the walk. At p = 2 the lines are drawn uniformly and most
  // of them give a neighbour; at an odd prime only about one line in p is
  // isotropic, so the isotropic ones are drawn directly.
  std::vector<int> walk_primes;
  if (rw.walk_prime_max <= 2) {
    walk_primes.push_back(2);
  } else {
    T det = DeterminantMat(state.ListGramWork[0]);
    int p_start = std::max(rw.walk_prime_min, 3);
    for (int p = GenusNextPrime(p_start - 1); p <= rw.walk_prime_max;
         p = GenusNextPrime(p)) {
      if (ResInt(det, T(p)) != T(0)) {
        walk_primes.push_back(p);
      }
    }
  }
  if (walk_primes.empty()) {
    std::cerr << "GENUS: no prime of [" << rw.walk_prime_min << ", "
              << rw.walk_prime_max << "] for the random walk does not divide "
              << "2 det\n";
    throw TerminalException{1};
  }
  size_t n_line_2 = NumberProjectiveLines(n, 2);
  // A random isotropic point mod p, by rejection: about p draws of O(n^2)
  // each, nothing next to the canonical form of the neighbour, and exactly
  // uniform over the isotropic points. A form without isotropic point mod p
  // (rank <= 2) gives nothing after 64 p tries.
  auto draw_line = [&](MyMatrix<Tint> const &G_int,
                       int p) -> std::optional<std::vector<int>> {
    if (p == 2) {
      return ProjectiveLineOfIndex(n, 2, random_index(engine, n_line_2));
    }
    std::vector<int64_t> G_mod(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        G_mod[i * n + j] = ResidueModP(G_int(i, j), p);
      }
    }
    int64_t p64 = p;
    std::vector<int> v(n);
    for (int i_try = 0; i_try < 64 * p; i_try++) {
      bool is_zero = true;
      for (int i = 0; i < n; i++) {
        v[i] = random_int(engine, 0, p - 1);
        if (v[i] != 0) {
          is_zero = false;
        }
      }
      if (is_zero) {
        continue;
      }
      int64_t q = 0;
      for (int i = 0; i < n; i++) {
        int64_t row = 0;
        for (int j = 0; j < n; j++) {
          row = (row + G_mod[i * n + j] * v[j]) % p64;
        }
        q = (q + row * v[i]) % p64;
      }
      if (q == 0) {
        return v;
      }
    }
    return {};
  };
  // The position, among the indices in vec, of the class of largest |Aut|,
  // the first discovered among the ties.
  auto get_max_aut = [&](std::vector<size_t> const &vec) -> size_t {
    size_t i_best = 0;
    for (size_t i = 1; i < vec.size(); i++) {
      if (result.ListAutOrder[vec[i]] > result.ListAutOrder[vec[i_best]]) {
        i_best = i;
      }
    }
    return i_best;
  };
  std::vector<size_t> pending;
  auto random_walk = [&](size_t pos) -> void {
    MyMatrix<Tint> G_int =
        UniversalMatrixConversion<Tint, T>(state.ListGramWork[pos]);
    size_t log_order = GenusRandomWalkLogOrder(result.ListAutOrder[pos]);
    size_t bound = std::max(
        size_t(1), static_cast<size_t>(std::ceil(rw.hits_factor *
                                                 static_cast<double>(log_order))));
    size_t max_failed = rw.failed_draws_factor * log_order;
    size_t n_hit = 0;
    size_t n_failed = 0;
    while (n_hit < bound && n_failed < max_failed) {
      // with one prime no draw, so that the p = 2 walk keeps its sequence
      int p = (walk_primes.size() == 1)
                  ? walk_primes[0]
                  : walk_primes[random_index(engine, walk_primes.size())];
      std::optional<std::vector<int>> opt_line = draw_line(G_int, p);
      if (!opt_line) {
        n_failed++;
        continue;
      }
      std::optional<MyMatrix<T>> opt =
          GenusNeighbor<T, Tint>(G_int, *opt_line, p, state.os);
      if (!opt) {
        n_failed++;
        continue;
      }
      result.n_neighbor_random++;
      if (state.insert(*opt)) {
        result.n_new_random++;
        pending.push_back(state.ListGramWork.size() - 1);
        if (state.is_finished()) {
          return;
        }
      } else {
        n_hit++;
      }
    }
#ifdef DEBUG_GENUS_ENUMERATION
    state.os << "GENUS: random walk on lattice " << (pos + 1) << " |Aut|="
             << result.ListAutOrder[pos] << " bound=" << bound
             << " failed=" << n_failed
             << " hits=" << n_hit << "\n";
#endif
  };
  auto run_pending = [&]() -> void {
    while (!pending.empty() && !state.is_finished()) {
      size_t i_best = get_max_aut(pending);
      size_t pos = pending[i_best];
      pending.erase(pending.begin() + i_best);
      random_walk(pos);
    }
  };
  for (size_t pos = 0; pos < state.ListGramWork.size(); pos++) {
    pending.push_back(pos);
  }
  run_pending();
  if (state.is_finished()) {
    result.closed_in_random_walk = state.is_mass_closed();
    return;
  }
#ifdef DEBUG_GENUS_ENUMERATION
  state.os << "GENUS: random walk phase over with "
           << state.ListGramWork.size() << " classes, full expansion\n";
#endif
  // The full expansion: is_expanded[i][pos] tells whether the class pos has
  // been expanded at ListPrimeUsed[i].
  std::vector<int> &ListPrimeUsed = result.ListPrimeUsed;
  std::vector<std::vector<bool>> is_expanded;
  while (!state.is_finished()) {
    if (result.accumulated_mass > state.TotalMass) {
      break;
    }
    size_t n_class = state.ListGramWork.size();
    size_t i_work = ListPrimeUsed.size();
    std::vector<size_t> unexpanded;
    for (size_t i = 0; i < ListPrimeUsed.size(); i++) {
      is_expanded[i].resize(n_class, false);
      for (size_t pos = 0; pos < n_class; pos++) {
        if (!is_expanded[i][pos]) {
          unexpanded.push_back(pos);
        }
      }
      if (!unexpanded.empty()) {
        i_work = i;
        break;
      }
    }
    if (i_work == ListPrimeUsed.size()) {
      int p_new = ListPrimeUsed.empty() ? 2 : GenusNextPrime(ListPrimeUsed.back());
      ListPrimeUsed.push_back(p_new);
      result.ListNewClassPerPrime.push_back(0);
      is_expanded.push_back(std::vector<bool>(n_class, false));
#ifdef DEBUG_GENUS_ENUMERATION
      state.os << "GENUS: random walk strategy brings in prime " << p_new
               << "\n";
#endif
      continue;
    }
    size_t pos = state.pick(unexpanded, ListPrimeUsed[i_work], i_work,
                            get_max_aut(unexpanded));
    is_expanded[i_work][pos] = true;
    // A new class found by the expansion gets its random walk at once,
    // before the expansion goes on.
    state.expand_at(pos, ListPrimeUsed[i_work], i_work,
                    [&](size_t idx) -> void {
                      pending.push_back(idx);
                      run_pending();
                    });
  }
}

/*
  The enumeration itself. Starts from the seed lattices, forms p-neighbours
  on the schedule of the strategy, keeps the new ones (recognised either by
  canonical form or by the isomorphism strategy), and stops when the
  accumulated mass reaches the target.
 */
template <typename T, typename Tint, typename Tgroup>
GenusEnumerationResult<T, Tgroup>
GenusEnumeration(std::vector<MyMatrix<T>> const &ListSeed,
                 typename GenusEnumerationResult<T, Tgroup>::Tmass const
                     &TotalMass,
                 int prime, GenusMethods const &methods,
                 NeighborStrategy strategy,
                 GenusRandomWalkOptions const &rw, bool hecke,
                 GenusExpansionOrder order, int max_runtime_second,
                 std::ostream &os) {
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_total;
#endif
  int n = ListSeed[0].rows();
  GenusEnumerationState<T, Tint, Tgroup> state(n, methods, TotalMass, hecke,
                                               order, os);
  state.max_runtime_second = max_runtime_second;
  for (auto &eSeed : ListSeed) {
    if (state.is_finished()) {
      break;
    }
    (void)state.insert(eSeed);
  }
  if (strategy == NeighborStrategy::OutPrime) {
    GenusEnumerationOutPrime(state, prime);
  } else if (strategy == NeighborStrategy::Iterate) {
    GenusEnumerationIterate(state);
  } else {
    GenusEnumerationRandomWalk(state, rw);
  }
  GenusEnumerationResult<T, Tgroup> &result = state.result;
  result.complete = state.is_mass_closed();
#ifdef DEBUG_GENUS_ISO_STATS
  os << "GENUS_ISO: FINAL isometry_tests=" << state.n_iso_test
     << " positive=" << state.n_iso_positive << " tests_per_positive="
     << (state.n_iso_positive > 0
             ? double(state.n_iso_test) / double(state.n_iso_positive)
             : 0.0)
     << "\n";
#endif
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: GenusEnumeration|=" << time_total << "\n";
#endif
  return std::move(result);
}

//
// Input and output
//

// The method of an enumeration, read from the namelist. `prime` is the prime
// of the out_prime strategy; 0 means "choose the smallest prime not dividing
// det". With the Iterate and RandomWalk strategies the prime is ignored.
struct GenusSpec {
  int prime;
  NeighborStrategy strategy;
  GenusMethods methods;
  GenusRandomWalkOptions rw;
  bool hecke;
  GenusExpansionOrder order;
};

/*
  The namelist of GENUS_Enumerate. The genus is given by its seed lattices,
  one per spinor genus, and its mass; the rest chooses the method.
 */
inline FullNamelist NAMELIST_GetStandard_GENUS_ENUMERATE() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  {
    std::map<std::string, std::string> ListStringValues_doc;
    ListStringValues_doc["arithmetic"] = "Default: gmp\n\
The arithmetic of the Gram matrices";
    ListStringValues_doc["FileLattice"] = "Default: unset\n\
One or more Gram matrices in the ListMatrix format, the seeds. One seed per\n\
SPINOR genus is needed: the p-neighbour graph is connected on a spinor\n\
genus, not on a genus, so a single seed may reach only part of the genus";
    ListStringValues_doc["FileMass"] = "Default: unset\n\
The mass of the genus, as \"num den\" or a single integer. The enumeration\n\
stops when sum 1/|Aut| reaches it, and that equality certifies completeness";
    ListBlock["DATA"].setListStringValues_doc(ListStringValues_doc);
  }
  // METHOD
  {
    std::map<std::string, std::string> ListStringValues_doc;
    std::map<std::string, std::string> ListIntValues_doc;
    std::map<std::string, std::string> ListBoolValues_doc;
    ListStringValues_doc["strategy"] = "Default: out_prime\n\
The schedule of the neighbour steps, see NeighborStrategy:\n\
out_prime: one prime, the given one or the smallest not dividing det\n\
iterate: the primes 2, 3, 5, ... in order, including those dividing det,\n\
  cheapest first, until the mass closes\n\
random_walk: every new class is first treated by random neighbours (the\n\
  RANDOM_WALK block), then the expansion of iterate";
    ListIntValues_doc["prime"] = "Default: 0\n\
The prime of out_prime, 0 for the smallest prime not dividing det. Ignored\n\
by iterate and random_walk";
    ListStringValues_doc["recognition"] = "Default: canonical\n\
How a known class is recognised, see GenusMethods: canonical (by the\n\
canonical form of canonical_method) or isomorphism (by the isometry test\n\
of isomorphism_method within the buckets of a cheap invariant)";
    ListStringValues_doc["canonical_method"] =
        "Default: weyl_chamber\nThe method of the canonical form, with "
        "recognition canonical. The methods are:" +
        DocLatticeStabCanIsomMethods();
    ListStringValues_doc["isomorphism_method"] =
        "Default: plesken_souvignier\nThe method of the isometry test, with "
        "recognition isomorphism. The methods are:" +
        DocLatticeStabCanIsomMethods();
    ListStringValues_doc["automorphism_method"] =
        "Default: weyl_chamber\nThe method of the automorphism group of every "
        "class. The methods are:" +
        DocLatticeStabCanIsomMethods();
    ListBoolValues_doc["hecke"] = "Default: F\n\
Whether to skip the expansions whose neighbours are all expanded classes,\n\
by the Hecke bookkeeping of GenusEnumerationState";
    ListStringValues_doc["order"] = "Default: deficit\n\
The order of the full expansions at a prime, see GenusExpansionOrder:\n\
aut: discovery order for out_prime and iterate, largest |Aut| first for\n\
  random_walk\n\
deficit: the class with the largest Hecke deficit first (implies hecke)";
    ListBlock["METHOD"].setListStringValues_doc(ListStringValues_doc);
    ListBlock["METHOD"].setListIntValues_doc(ListIntValues_doc);
    ListBlock["METHOD"].setListBoolValues_doc(ListBoolValues_doc);
  }
  // RANDOM_WALK
  {
    std::map<std::string, std::string> ListIntValues_doc;
    std::map<std::string, std::string> ListDoubleValues_doc;
    ListIntValues_doc["seed"] = "Default: 1\n\
The seed of the random walk";
    ListDoubleValues_doc["hits_factor"] = "Default: 0.25\n\
A class with |Aut| = g, and b the smallest integer with g < 2^b, is treated\n\
when max(1, ceil(hits_factor * b)) of its random neighbours are known\n\
classes";
    ListIntValues_doc["failed_draws_factor"] = "Default: 16\n\
The draws giving no neighbour are bounded by failed_draws_factor * b";
    ListIntValues_doc["walk_prime_min"] = "Default: 2\n\
The smallest prime of the random walk";
    ListIntValues_doc["walk_prime_max"] = "Default: 2\n\
The largest prime of the random walk: 2 for the walk at p = 2, otherwise\n\
the primes of [walk_prime_min, walk_prime_max] not dividing 2 det";
    ListBlock["RANDOM_WALK"].setListIntValues_doc(ListIntValues_doc);
    ListBlock["RANDOM_WALK"].setListDoubleValues_doc(ListDoubleValues_doc);
  }
  // SYSTEM
  {
    std::map<std::string, std::string> ListStringValues_doc;
    ListStringValues_doc["OutFormat"] = "Default: GAP\n\
GAP: a record with the Gram matrices, the group orders, the mass and the\n\
  counters of the enumeration\n\
CPP: the list of Gram matrices\n\
Summary: the class number, the mass check and the counters";
    ListStringValues_doc["OutFile"] = "Default: stderr\n\
The file of the output, or stderr, or stdout";
    std::map<std::string, std::string> ListIntValues_doc;
    ListIntValues_doc["max_runtime_second"] = "Default: 0\n\
The maximum runtime of the enumeration in seconds (0 = unlimited). When it\n\
is reached the classes found so far are written, with complete false";
    ListBlock["SYSTEM"].setListStringValues_doc(ListStringValues_doc);
    ListBlock["SYSTEM"].setListIntValues_doc(ListIntValues_doc);
  }
  return FullNamelist(ListBlock);
}

/*
  The method part of the namelist. The genus itself is not described: the
  p-neighbour construction preserves it, and the seeds are checked against
  each other (same rank, same determinant) by GENUS_Enumerate.
 */
inline GenusSpec ReadGenusSpec(FullNamelist const &eFull) {
  SingleBlock const &BlockMETHOD = eFull.get_block("METHOD");
  SingleBlock const &BlockRW = eFull.get_block("RANDOM_WALK");
  GenusSpec spec;
  spec.prime = BlockMETHOD.get_int("prime");
  std::string strat = BlockMETHOD.get_string("strategy");
  if (strat == "out_prime") {
    spec.strategy = NeighborStrategy::OutPrime;
  } else if (strat == "iterate") {
    spec.strategy = NeighborStrategy::Iterate;
  } else if (strat == "random_walk") {
    spec.strategy = NeighborStrategy::RandomWalk;
  } else {
    std::cerr << "GENUS: unrecognised strategy \"" << strat
              << "\". Allowed: out_prime, iterate, random_walk\n";
    throw TerminalException{1};
  }
  std::string recognition = BlockMETHOD.get_string("recognition");
  if (recognition == "canonical") {
    spec.methods.by_canonical = true;
  } else if (recognition == "isomorphism") {
    spec.methods.by_canonical = false;
  } else {
    std::cerr << "GENUS: unrecognised recognition \"" << recognition
              << "\". Allowed: canonical, isomorphism\n";
    throw TerminalException{1};
  }
  spec.methods.canonical = ParseLatticeStabCanIsomMethod(
      BlockMETHOD.get_string("canonical_method"));
  spec.methods.isomorphism = ParseLatticeStabCanIsomMethod(
      BlockMETHOD.get_string("isomorphism_method"));
  spec.methods.automorphism = ParseLatticeStabCanIsomMethod(
      BlockMETHOD.get_string("automorphism_method"));
  spec.hecke = BlockMETHOD.get_bool("hecke");
  std::string order = BlockMETHOD.get_string("order");
  if (order == "aut") {
    spec.order = GenusExpansionOrder::Default;
  } else if (order == "deficit") {
    spec.order = GenusExpansionOrder::Deficit;
  } else {
    std::cerr << "GENUS: unrecognised order \"" << order
              << "\". Allowed: aut, deficit\n";
    throw TerminalException{1};
  }
  int seed = BlockRW.get_int("seed");
  int failed_draws_factor = BlockRW.get_int("failed_draws_factor");
  if (seed < 0 || failed_draws_factor < 0) {
    std::cerr << "GENUS: seed and failed_draws_factor must be nonnegative\n";
    throw TerminalException{1};
  }
  spec.rw.seed = static_cast<uint64_t>(seed);
  spec.rw.hits_factor = BlockRW.get_double("hits_factor");
  spec.rw.failed_draws_factor = static_cast<size_t>(failed_draws_factor);
  spec.rw.walk_prime_min = BlockRW.get_int("walk_prime_min");
  spec.rw.walk_prime_max = BlockRW.get_int("walk_prime_max");
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
    os << "neighbours built: full expansion " << result.n_neighbor_full
       << ", random walk " << result.n_neighbor_random << " ("
       << result.n_new_random << " new classes)\n";
    os << "expansions: " << result.n_expansion << " ("
       << result.n_expansion_skipped << " skipped by the Hecke bookkeeping)\n";
    os << "closed in the random walk phase = "
       << (result.closed_in_random_walk ? "true" : "false") << "\n";
    os << "stopped by max_runtime_second = "
       << (result.stopped_by_runtime ? "true" : "false") << "\n";
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
       << ",\n stopped_by_runtime:="
       << (result.stopped_by_runtime ? "true" : "false")
       << ",\n prime:=" << prime << ",\n ListPrimeUsed:=[";
    for (size_t i = 0; i < result.ListPrimeUsed.size(); i++) {
      os << (i > 0 ? ", " : "") << result.ListPrimeUsed[i];
    }
    os << "],\n ListNewClassPerPrime:=[";
    for (size_t i = 0; i < result.ListNewClassPerPrime.size(); i++) {
      os << (i > 0 ? ", " : "") << result.ListNewClassPerPrime[i];
    }
    os << "],\n n_neighbor_full:=" << result.n_neighbor_full
       << ",\n n_neighbor_random:=" << result.n_neighbor_random
       << ",\n n_new_random:=" << result.n_new_random
       << ",\n n_expansion:=" << result.n_expansion
       << ",\n n_expansion_skipped:=" << result.n_expansion_skipped
       << ",\n closed_in_random_walk:="
       << (result.closed_in_random_walk ? "true" : "false") << ");\n";
    return;
  }
  std::cerr << "GENUS: failed to find a matching entry for OutFormat="
            << OutFormat << "\n";
  throw TerminalException{1};
}

#endif  //  SRC_GENERA_GENUSENUMERATION_H_
