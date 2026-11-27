// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSINDEFINITE_H_
#define SRC_GENERA_GENUSINDEFINITE_H_

// clang-format off
#include "GenusSpinor.h"
#include "GenusEnumeration.h"
#include "FiniteFieldProjective.h"
#include "Indefinite_LLL.h"
#include "SignatureSymmetric.h"
#include <map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
# define DEBUG_GENUS_INDEFINITE
#endif

#ifdef SANITY_CHECK
# define SANITY_CHECK_GENUS_INDEFINITE
#endif

/*
  The classes of an indefinite genus of lattices over Z of rank >= 3.

  By Eichler's theorem (strong approximation for the spin group) the
  classes of an indefinite genus of rank >= 3 are its spinor genera, and
  the isometry classes are its improper spinor genera: the cosets of the
  improper spinor kernel K in A = prod_{p | 2d} Z_p^x / Z_p^{x2}
  (GenusSpinor.h, with Delta(-1) in K since the real place is isotropic).
  An l-neighbor, l a prime not dividing 2d, moves a lattice from the coset
  c to c + Delta(l). So the representatives are found by a search on the
  group Q = A / K: from the given lattice, the l-neighbors of every lattice
  found, for l = 3, 5, 7, ..., until every element of Q has a lattice.
  Unlike taking one neighbor of the given lattice per generator of Q, this
  reaches the sums of the generators too.

  The neighbors are those of GenusEnumeration.h, of the even lattice: L for
  an even L, L(2) (the form doubled) for an odd one, as l is odd. Each new
  lattice is reduced by the indefinite LLL of Indefinite_LLL.h, so that the
  entries do not grow along the search.
 */

template <typename T> std::pair<int, int> GenusIndefiniteSignature(MyMatrix<T> const &G) {
  DiagSymMat<T> diag = DiagonalizeSymmetricMatrix(G, std::cerr);
  return {diag.nbPlus, diag.nbMinus};
}

template <typename T, typename Tint> struct GenusIndefiniteResult {
  // the number of classes, |A / K|
  int n_class;
  // the Gram matrices of the representatives, the given lattice first
  std::vector<MyMatrix<Tint>> l_gram;
  // the primes l used
  std::vector<int> l_prime;
};

/*
  An l-neighbor of the lattice of Gram matrix G (even), or nothing when no
  line of P^{n-1}(F_l) among the first max_lines gives one.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>> GenusIndefiniteNeighbor(MyMatrix<Tint> const &G,
                                                   int l, size_t max_lines,
                                                   std::ostream &os) {
  int n = G.rows();
  size_t n_line = NumberProjectiveLines(n, l);
  for (size_t idx = 0; idx < std::min(n_line, max_lines); idx++) {
    std::vector<int> v = ProjectiveLineOfIndex(n, l, idx);
    std::optional<GenusNeighborBasis<T, Tint>> opt =
        GenusNeighborKernel<T, Tint>(G, v, l, os);
    if (opt) {
      return opt->GramN;
    }
  }
  return {};
}

template <typename T, typename Tint>
GenusIndefiniteResult<T, Tint>
GenusIndefiniteRepresentatives(MyMatrix<Tint> const &G, std::ostream &os) {
  MyMatrix<T> G_T = UniversalMatrixConversion<T, Tint>(G);
  GenusSymbol s = GenusSymbolOfGram<T, Tint>(G_T);
  SpinorGenusData data = GenusSpinorData<T, Tint>(s, G_T, true);
  SpinorSubspace const &K = data.kernel();
  GenusIndefiniteResult<T, Tint> res;
  res.n_class = data.n_spinor_genera();
  Tint det = T_abs(UniversalScalarConversion<Tint, T>(DeterminantMat(G_T)));
  // the even lattice: G, or 2 G for an odd one (the neighbors at odd l of
  // L(2) are the L'(2), L' the neighbors of L)
  bool even = s.even;
  Tint scale = even ? Tint(1) : Tint(2);
  // the lattices found, by their coset of K (reduced by the echelon basis)
  std::map<uint64_t, size_t> coset_index;
  coset_index[K.reduce(0)] = 0;
  res.l_gram.push_back(G);
  std::vector<uint64_t> l_coset{K.reduce(0)};
  auto reduce = [&](MyMatrix<T> const &M) -> MyMatrix<Tint> {
    ResultReduction<T, Tint> red = IndefiniteReduction<T, Tint>(M, os);
    return UniversalMatrixConversion<Tint, T>(red.Mred);
  };
  for (int l = 3; res.l_gram.size() < static_cast<size_t>(res.n_class); l += 2) {
    bool is_prime = true;
    for (int q = 3; q * q <= l; q += 2) {
      if (l % q == 0) {
        is_prime = false;
      }
    }
    if (!is_prime || ResInt(det, Tint(l)) == 0) {
      continue;
    }
    uint64_t delta = K.reduce(SpinorDelta(data.A, Tint(l)));
    if (delta == 0) {
      continue;
    }
    res.l_prime.push_back(l);
    // every lattice found so far gives the coset shifted by Delta(l)
    size_t n_now = res.l_gram.size();
    for (size_t i = 0; i < n_now; i++) {
      uint64_t c = K.reduce(l_coset[i] ^ delta);
      if (coset_index.count(c)) {
        continue;
      }
      MyMatrix<Tint> Ge = scale * res.l_gram[i];
      std::optional<MyMatrix<T>> opt =
          GenusIndefiniteNeighbor<T, Tint>(Ge, l, 1000000, os);
#ifdef SANITY_CHECK_GENUS_INDEFINITE
      if (!opt) {
        std::cerr << "GINDEF: no " << l << "-neighbor found, but an "
                  << "indefinite lattice of rank >= 3 has isotropic lines "
                     "modulo l\n";
        throw TerminalException{1};
      }
#endif
      MyMatrix<T> M = *opt / UniversalScalarConversion<T, Tint>(scale);
      MyMatrix<Tint> M_red = reduce(M);
#ifdef SANITY_CHECK_GENUS_INDEFINITE
      if (GenusSymbolOfGram<T, Tint>(UniversalMatrixConversion<T, Tint>(M_red)) != s) {
        std::cerr << "GINDEF: a neighbor is not in the genus\n";
        throw TerminalException{1};
      }
#endif
      coset_index[c] = res.l_gram.size();
      res.l_gram.push_back(M_red);
      l_coset.push_back(c);
#ifdef DEBUG_GENUS_INDEFINITE
      os << "GINDEF: the coset " << c << " by a " << l << "-neighbor, "
         << res.l_gram.size() << " of " << res.n_class << " classes\n";
#endif
    }
#ifdef SANITY_CHECK_GENUS_INDEFINITE
    if (l > 100000) {
      std::cerr << "GINDEF: the primes up to 100000 do not generate A / K\n";
      throw TerminalException{1};
    }
#endif
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_GENUSINDEFINITE_H_
// clang-format on
