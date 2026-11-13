// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSENUMERATIONODD_H_
#define SRC_GENERA_GENUSENUMERATIONODD_H_

// clang-format off
#include "GenusEnumeration.h"
#include "LocalDensity.h"
#include "GenusMass.h"
#include "LocalDensityOverlattice.h"
#include <cstdint>
#include <unordered_map>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_GENUS_ENUMERATION_ODD
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_GENUS_ENUMERATION_ODD
#endif

/*
  The enumeration of a genus of odd integral positive definite lattices
  through their even sublattices, which keeps the prime 2 and the root
  systems of the even enumeration.

  For an odd lattice L, x -> b(x, x) mod 2 is additive, so the even
  sublattice L_0 = { x : b(x, x) even } is the kernel of a linear form mod 2
  and has index 2. It is determined by L, so an isometry L -> L' maps L_0
  onto L'_0, and L_0 (x) Z_q is the even sublattice of L (x) Z_q: the genus
  of L_0 depends only on that of L. Hence
  * the even sublattices of the lattices of the genus G of L are the
    classes of the genus G_0 of L_0;
  * the lattices of G with even sublattice isometric to a class M of G_0 are
    the odd overlattices N of index 2 of M that lie in G, and two of them
    are isometric iff they are in one Aut(M)-orbit, an isometry between
    them preserving M;
  * Aut(N) is the stabilizer of N in Aut(M), of order |Aut M| / |orbit|.
  An index 2 overlattice is N = M + Z x/2 with x in M, x/2 in M^#, that is
  M x = 0 mod 2 (x in the radical of the form mod 2, as a row vector), and
  N is odd iff b(x/2, x/2) is odd, that is x M x^T = 4 mod 8; N depends only
  on x mod 2M. So the classes of G are obtained, class M by class M of G_0,
  as the orbits of Aut(M) on these x mod 2 that give a lattice of G (same
  2-adic symbol as L, the other primes being unchanged), with no canonical
  form. The enumeration of G_0 is certified by its own mass, and the masses
  1/|Aut N| must add up to the mass of G, a second certificate.
 */

// A basis (rows) of the even sublattice of the odd integral form GramMat.
template <typename T, typename Tint>
MyMatrix<Tint> EvenSublatticeBasis(MyMatrix<T> const &GramMat) {
  int n = GramMat.rows();
  int j_odd = -1;
  for (int i = 0; i < n; i++) {
    if (!IsInteger(GramMat(i, i) / T(2))) {
      j_odd = i;
      break;
    }
  }
  if (j_odd < 0) {
    std::cerr << "GENUS_ODD: the form is even, it has no even sublattice of "
                 "index 2\n";
    throw TerminalException{1};
  }
  // b(x, x) = sum_i a_ii x_i mod 2: the e_i with a_ii even, e_i + e_j for
  // the other i, and 2 e_j, with a_jj odd.
  MyMatrix<Tint> B = ZeroMatrix<Tint>(n, n);
  for (int i = 0; i < n; i++) {
    if (i == j_odd) {
      B(i, i) = 2;
    } else if (IsInteger(GramMat(i, i) / T(2))) {
      B(i, i) = 1;
    } else {
      B(i, i) = 1;
      B(i, j_odd) = 1;
    }
  }
  return B;
}

template <typename T, typename Tint>
MyMatrix<T> EvenSublatticeGram(MyMatrix<T> const &GramMat) {
  MyMatrix<T> B = UniversalMatrixConversion<T, Tint>(
      EvenSublatticeBasis<T, Tint>(GramMat));
  return B * GramMat * B.transpose();
}

// The 2-adic symbol of a lattice of the genus G of an odd lattice, taken on
// the even lattice with twice its form, so that the even code applies.
template <typename T, typename Tint>
std::vector<int> OddTwoAdicKey(MyMatrix<T> const &GramMat) {
  MyMatrix<mpq_class> G2 =
      UniversalMatrixConversion<mpq_class, T>(GramMat) * mpq_class(2);
  return LocalSymbolKey<mpq_class, mpz_class>(G2, 2);
}

/*
  The classes of G above one class M of G_0, appended to result: the orbits
  of the automorphism group of M (generators ListGen, order aut_order) on
  the x mod 2 described above.
 */
template <typename T, typename Tint, typename Tgroup>
void GenusOddOverlattices(MyMatrix<T> const &GramM,
                          std::vector<MyMatrix<Tint>> const &ListGen,
                          typename Tgroup::Tint const &aut_order,
                          std::vector<int> const &key_odd,
                          GenusEnumerationResult<T, Tgroup> &result,
                          [[maybe_unused]] std::ostream &os) {
  using TintGroup = typename Tgroup::Tint;
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  int n = GramM.rows();
  if (n > 64) {
    std::cerr << "GENUS_ODD: rank " << n << " above 64\n";
    throw TerminalException{1};
  }
  MyMatrix<Tint> G = UniversalMatrixConversion<Tint, T>(GramM);
  // The rows of G mod 2 as masks, and a basis of the radical of G mod 2
  // (the x with x G = 0 mod 2) by elimination over F_2.
  std::vector<uint64_t> col(n, 0);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (ResidueModP(G(i, j), 2) != 0) {
        col[j] |= uint64_t(1) << i;
      }
    }
  }
  // Solve the system col[j] . x = 0 for all j: row reduce the n equations.
  std::vector<uint64_t> eqs = col;
  std::vector<int> pivot_of_row;
  std::vector<uint64_t> red;
  for (uint64_t e : eqs) {
    for (size_t k = 0; k < red.size(); k++) {
      if ((e >> pivot_of_row[k]) & 1) {
        e ^= red[k];
      }
    }
    if (e != 0) {
      int piv = __builtin_ctzll(e);
      for (size_t k = 0; k < red.size(); k++) {
        if ((red[k] >> piv) & 1) {
          red[k] ^= e;
        }
      }
      red.push_back(e);
      pivot_of_row.push_back(piv);
    }
  }
  std::vector<bool> is_pivot(n, false);
  for (int piv : pivot_of_row) {
    is_pivot[piv] = true;
  }
  std::vector<uint64_t> radical;
  for (int f = 0; f < n; f++) {
    if (is_pivot[f]) {
      continue;
    }
    uint64_t x = uint64_t(1) << f;
    for (size_t k = 0; k < red.size(); k++) {
      if ((red[k] >> f) & 1) {
        x |= uint64_t(1) << pivot_of_row[k];
      }
    }
    radical.push_back(x);
  }
  int r = static_cast<int>(radical.size());
  if (r > 24) {
    std::cerr << "GENUS_ODD: a radical mod 2 of dimension " << r
              << ", too large to list\n";
    throw TerminalException{1};
  }
  // The candidates: the nonzero x of the radical with x G x^T = 4 mod 8.
  auto norm_mod8 = [&](uint64_t x) -> int64_t {
    int64_t s = 0;
    for (int i = 0; i < n; i++) {
      if (!((x >> i) & 1)) {
        continue;
      }
      for (int j = 0; j < n; j++) {
        if ((x >> j) & 1) {
          s += ResidueModP(G(i, j), 8);
        }
      }
    }
    return ((s % 8) + 8) % 8;
  };
  std::vector<uint64_t> cand;
  for (uint64_t c = 1; c < (uint64_t(1) << r); c++) {
    uint64_t x = 0;
    for (int k = 0; k < r; k++) {
      if ((c >> k) & 1) {
        x ^= radical[k];
      }
    }
    if (norm_mod8(x) == 4) {
      cand.push_back(x);
    }
  }
#ifdef DEBUG_GENUS_ENUMERATION_ODD
  os << "GENUS_ODD: radical mod 2 of dimension " << r << ", " << cand.size()
     << " candidate vectors\n";
#endif
  if (cand.empty()) {
    return;
  }
  std::unordered_map<uint64_t, size_t> index_of;
  for (size_t i = 0; i < cand.size(); i++) {
    index_of[cand[i]] = i;
  }
  // The generators mod 2, acting on row vectors: x -> x g, the XOR of the
  // rows of g over the bits of x.
  std::vector<std::vector<uint64_t>> ListRows;
  for (auto &g : ListGen) {
    std::vector<uint64_t> rows(n, 0);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (ResidueModP(g(i, j), 2) != 0) {
          rows[i] |= uint64_t(1) << j;
        }
      }
    }
    ListRows.push_back(std::move(rows));
  }
  auto image = [&](uint64_t x, std::vector<uint64_t> const &rows) -> uint64_t {
    uint64_t w = 0;
    for (int i = 0; i < n; i++) {
      if ((x >> i) & 1) {
        w ^= rows[i];
      }
    }
    return w;
  };
  std::vector<bool> seen(cand.size(), false);
  for (size_t i0 = 0; i0 < cand.size(); i0++) {
    if (seen[i0]) {
      continue;
    }
    seen[i0] = true;
    std::vector<size_t> active{i0};
    size_t orbit_size = 0;
    while (!active.empty()) {
      size_t i = active.back();
      active.pop_back();
      orbit_size++;
      for (auto &rows : ListRows) {
        uint64_t y = image(cand[i], rows);
        auto iter = index_of.find(y);
#ifdef SANITY_CHECK_GENUS_ENUMERATION_ODD
        if (iter == index_of.end()) {
          std::cerr << "GENUS_ODD: an automorphism maps a candidate outside "
                       "the candidates\n";
          throw TerminalException{1};
        }
#endif
        if (!seen[iter->second]) {
          seen[iter->second] = true;
          active.push_back(iter->second);
        }
      }
    }
    // The overlattice of the representative: e_k replaced by x/2, k a bit
    // of x, a basis since e_k = 2 (x/2) - sum_{i != k} x_i e_i.
    uint64_t x = cand[i0];
    int k = __builtin_ctzll(x);
    MyMatrix<T> B = IdentityMat<T>(n);
    for (int j = 0; j < n; j++) {
      B(k, j) = ((x >> j) & 1) ? T(1) / T(2) : T(0);
    }
    MyMatrix<T> GramN = B * GramM * B.transpose();
    if (OddTwoAdicKey<T, Tint>(GramN) != key_odd) {
      continue;
    }
    TintGroup order_n = aut_order / TintGroup(orbit_size);
#ifdef SANITY_CHECK_GENUS_ENUMERATION_ODD
    if (order_n * TintGroup(orbit_size) != aut_order) {
      std::cerr << "GENUS_ODD: an orbit of size " << orbit_size
                << " does not divide |Aut| = " << aut_order << "\n";
      throw TerminalException{1};
    }
#endif
    result.ListGram.push_back(GramN);
    result.ListAutOrder.push_back(order_n);
    result.accumulated_mass +=
        Tmass(1) / UniversalScalarConversion<Tmass, TintGroup>(order_n);
  }
}

/*
  The classes of the genus of the odd seeds ListSeed, through the genus of
  their even sublattices. TotalMassOdd is the mass of the odd genus; the
  mass of the even genus is computed. The other arguments are those of
  GenusEnumeration, which enumerates the even genus; the counters of the
  result are those of that enumeration.
 */
template <typename T, typename Tint, typename Tgroup>
GenusEnumerationResult<T, Tgroup> GenusEnumerationOdd(
    std::vector<MyMatrix<T>> const &ListSeed,
    typename GenusEnumerationResult<T, Tgroup>::Tmass const &TotalMassOdd,
    int prime, GenusMethods const &methods, NeighborStrategy strategy,
    GenusRandomWalkOptions const &rw, bool hecke,
    GenusOrderOptions const &order, int max_runtime_second, std::ostream &os) {
  using Tmass = typename GenusEnumerationResult<T, Tgroup>::Tmass;
  std::vector<MyMatrix<T>> ListSeedEven;
  for (auto &eSeed : ListSeed) {
    ListSeedEven.push_back(EvenSublatticeGram<T, Tint>(eSeed));
  }
  Tmass TotalMassEven = GenusMass<mpq_class, mpz_class>(
      UniversalMatrixConversion<mpq_class, T>(ListSeedEven[0]));
#ifdef DEBUG_GENUS_ENUMERATION_ODD
  os << "GENUS_ODD: the even sublattice genus has the mass " << TotalMassEven
     << "\n";
#endif
  GenusEnumerationResult<T, Tgroup> resEven = GenusEnumeration<T, Tint, Tgroup>(
      ListSeedEven, TotalMassEven, prime, methods, strategy, rw, hecke, order,
      max_runtime_second, os);
  GenusEnumerationResult<T, Tgroup> result = resEven;
  result.ListGram.clear();
  result.ListAutOrder.clear();
  result.accumulated_mass = Tmass(0);
  result.target_mass = TotalMassOdd;
  std::vector<int> key_odd = OddTwoAdicKey<T, Tint>(ListSeed[0]);
  for (size_t i = 0; i < resEven.ListGram.size(); i++) {
    MyMatrix<T> const &GramM = resEven.ListGram[i];
    std::vector<MyMatrix<T>> ListMatM{GramM};
    IntAutomorphismGroup<Tint, Tgroup> aut =
        LatticeAutomorphismGroup<T, Tint, Tgroup>(ListMatM,
                                                   methods.automorphism, os);
#ifdef SANITY_CHECK_GENUS_ENUMERATION_ODD
    if (aut.order != resEven.ListAutOrder[i]) {
      std::cerr << "GENUS_ODD: the automorphism group of class " << i
                << " has order " << aut.order << " against "
                << resEven.ListAutOrder[i] << " in the enumeration\n";
      throw TerminalException{1};
    }
#endif
    GenusOddOverlattices<T, Tint, Tgroup>(GramM, aut.ListGen, aut.order,
                                          key_odd, result, os);
  }
  result.complete =
      resEven.complete && result.accumulated_mass == TotalMassOdd;
  return result;
}

// clang-format off
#endif  // SRC_GENERA_GENUSENUMERATIONODD_H_
// clang-format on
