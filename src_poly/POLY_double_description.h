// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DOUBLE_DESCRIPTION_H_
#define SRC_POLY_POLY_DOUBLE_DESCRIPTION_H_

// clang-format off
#include "Boost_bitset.h"
#include "MAT_Matrix.h"
#include "MAT_NonUniqueRescale.h"
#include "POLY_SimplexClarkson.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_DOUBLE_DESCRIPTION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_DOUBLE_DESCRIPTION
#endif

#ifdef TIMINGS
#define TIMINGS_DOUBLE_DESCRIPTION
#endif

/*
  A from-scratch implementation of the double description method.

  Credit where it is due: the algorithm follows
  --- K. Fukuda, A. Prodon, "Double Description Method Revisited",
      Combinatorics and Computer Science, LNCS 1120 (1996), 91--111,
  --- the cddlib implementation of Komei Fukuda (dd_DoubleDescription in
      cddcore.c), from which the earlier in-tree C++ translation in
      POLY_cddlib.h also derives,
  --- with design hints from the conversion procedure of the Parma
      Polyhedra Library (integer coefficients with gcd normalization and
      cached scalar products).

  What differs from the cddlib design:
  --- No division: the whole computation runs in a ring (e.g. mpz_class).
      The only combination used, w = (a.r+) r- - (a.r-) r+, is division
      free, and the rays are reduced by their content (an exact division)
      to keep the entries small. For a field type T the public entry
      points scale the rows once to the underlying ring, in the same way
      as POLY_lrslib.h, the fraction-free simplex and the beneath-beyond
      method.
  --- Standard containers: the rays live in a std::vector instead of the
      linked list of cddlib, and the zero sets are Face bitsets.
  --- The zero sets are maintained incrementally over the processed rows
      only, so no masking with the set of added halfspaces is needed in
      the adjacency test, and the final zero sets are exactly the facet
      incidences.

  The conventions: the input EXT is a m x d matrix whose rows generate a
  cone; the output is the dual description, i.e. the facets of that cone,
  which are the extreme rays of { x : EXT x >= 0 }. DualDescription_incd
  returns for each facet the set of rows lying on it, matching
  cdd::DualDescription_incd; DualDescription returns the facet vectors.

  The insertion order of the rows is the lexicographic-minimum order,
  the robust default of Fukuda-Prodon and cddlib.
 */

namespace double_desc {

// Reduce the vector by its content, an exact division that keeps the
// entries small along the computation. The gcd is seeded at zero and
// computed with an early exit, and is positive so that no sign flip
// occurs.
template <typename T> void NormalizeRayContent(MyVector<T> &V) {
  int n = V.size();
  T eGCD(0);
  for (int i = 0; i < n; i++) {
    eGCD = GcdPair(eGCD, V(i));
    if (eGCD == 1) {
      return;
    }
  }
  if constexpr (is_totally_ordered<T>::value) {
    if (eGCD < 0) {
      eGCD = -eGCD;
    }
  }
  if (eGCD == 0 || eGCD == 1) {
    return;
  }
  for (int i = 0; i < n; i++) {
    V(i) = V(i) / eGCD;
  }
}

/*
  Greedy selection of a maximal set of linearly independent rows among the
  listed candidates, in the given order, by division-free elimination.
  Each candidate row is reduced against the echelon rows accumulated so
  far using cross-multiplication only; the reduced rows are content
  normalized to keep the entries small. Usable over a ring.
 */
template <typename T>
std::vector<int> SelectIndependentRows(MyMatrix<T> const &M,
                                       std::vector<int> const &candidates,
                                       size_t const &max_rank) {
  int nbCol = M.cols();
  std::vector<MyVector<T>> echelon;
  std::vector<int> pivcol;
  std::vector<int> selected;
  for (auto &iRow : candidates) {
    MyVector<T> V = GetMatrixRow(M, iRow);
    for (size_t k = 0; k < echelon.size(); k++) {
      int c = pivcol[k];
      if (V(c) != 0) {
        T coef1 = echelon[k](c);
        T coef2 = V(c);
        for (int j = 0; j < nbCol; j++) {
          V(j) = V(j) * coef1 - coef2 * echelon[k](j);
        }
        NormalizeRayContent(V);
      }
    }
    int c_piv = -1;
    for (int j = 0; j < nbCol; j++) {
      if (V(j) != 0) {
        c_piv = j;
        break;
      }
    }
    if (c_piv >= 0) {
      echelon.push_back(V);
      pivcol.push_back(c_piv);
      selected.push_back(iRow);
      if (selected.size() == max_rank) {
        return selected;
      }
    }
  }
  return selected;
}

// Greedy selection of a maximal independent set of columns, done by the
// division-free row selection on the transposed matrix.
template <typename T>
std::vector<int> SelectIndependentColumns(MyMatrix<T> const &M) {
  MyMatrix<T> Mtr = TransposedMat(M);
  int nbRow = Mtr.rows();
  std::vector<int> candidates(nbRow);
  for (int i = 0; i < nbRow; i++) {
    candidates[i] = i;
  }
  return SelectIndependentRows(Mtr, candidates, nbRow);
}

// A ray of the intermediate cone: its coordinates, the set of processed
// rows it is tight on, and the cached scalar product with the row being
// inserted.
template <typename T> struct RayEntry {
  MyVector<T> coord;
  Face zero;
  T val;
};

/*
  The kernel of the double description method over the (ring) type T.
  Precondition: EXT has full column rank d. Zero rows are harmless: they
  are tight on every ray. Returns the list of extreme rays of
  { x : EXT x >= 0 } together with their tight sets, which are the facets
  of the cone generated by the rows of EXT with their incidences.
 */
template <typename T>
std::vector<RayEntry<T>> DoubleDescription_Kernel(MyMatrix<T> const &EXT,
                                                  [[maybe_unused]] std::ostream &os) {
  int m = EXT.rows();
  int d = EXT.cols();
#ifdef TIMINGS_DOUBLE_DESCRIPTION
  MicrosecondTime time;
#endif
  // The lexicographic-minimum insertion order.
  std::vector<int> order(m);
  for (int i = 0; i < m; i++) {
    order[i] = i;
  }
  std::sort(order.begin(), order.end(), [&](int i1, int i2) -> bool {
    for (int j = 0; j < d; j++) {
      if (EXT(i1, j) != EXT(i2, j)) {
        return EXT(i1, j) < EXT(i2, j);
      }
    }
    return i1 < i2;
  });
  // The initial basis: the first d independent rows in the order.
  std::vector<int> basis = SelectIndependentRows(EXT, order, d);
#ifdef SANITY_CHECK_DOUBLE_DESCRIPTION
  if (static_cast<int>(basis.size()) != d) {
    std::cerr << "DOUBLEDESC: The matrix does not have full column rank: "
              << "rank=" << basis.size() << " d=" << d << "\n";
    throw TerminalException{1};
  }
#endif
  // The initial rays: the columns of the scaled inverse of the basis
  // matrix B, with B r_j = |det B| e_j, computed fraction-free.
  MyMatrix<T> B(d, d);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      B(i, j) = EXT(basis[i], j);
    }
  }
  std::vector<RayEntry<T>> rays;
  for (int j = 0; j < d; j++) {
    MyVector<T> e_j = ZeroVector<T>(d);
    e_j(j) = T(1);
    std::optional<std::pair<MyVector<T>, T>> opt =
        FractionFreeSolveSquare(B, e_j);
#ifdef SANITY_CHECK_DOUBLE_DESCRIPTION
    if (!opt) {
      std::cerr << "DOUBLEDESC: The basis matrix is singular, which "
                   "contradicts the independent row selection\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> r = opt->first;
    if (opt->second < 0) {
      for (int k = 0; k < d; k++) {
        r(k) = -r(k);
      }
    }
    NormalizeRayContent(r);
    Face zero(m);
    for (auto &iRow : basis) {
      if (iRow != basis[j]) {
        zero[iRow] = 1;
      }
    }
    rays.push_back({std::move(r), std::move(zero), T(0)});
  }
#ifdef DEBUG_DOUBLE_DESCRIPTION
  os << "DOUBLEDESC: initial basis built, d=" << d << " m=" << m << "\n";
#endif
  // The main loop: insert the remaining rows in the order.
  Face processed(m);
  for (auto &iRow : basis) {
    processed[iRow] = 1;
  }
  for (auto &iRow : order) {
    if (processed[iRow] == 1) {
      continue;
    }
    processed[iRow] = 1;
    int n_neg = 0;
    for (auto &ray : rays) {
      T eSum(0);
      for (int j = 0; j < d; j++) {
        eSum += EXT(iRow, j) * ray.coord(j);
      }
      ray.val = eSum;
      if (eSum < 0) {
        n_neg++;
      } else {
        if (eSum == 0) {
          ray.zero[iRow] = 1;
        }
      }
    }
    if (n_neg == 0) {
      continue;
    }
    // The new rays from the adjacent pairs (positive, negative). The
    // common tight set bitset is hoisted out of the loops (the allocation
    // per pair dominates the profile otherwise) and the tight set sizes
    // are precomputed: a set of size c can only be contained in a tight
    // set of size >= c. The subset scan is restricted to the rays tight
    // on the least populated row of the common tight set: any ray whose
    // tight set contains the common set is in particular tight on every
    // one of its rows, so the restriction loses nothing.
    size_t n_rays = rays.size();
    std::vector<RayEntry<T>> new_rays;
    // A packed word representation of the tight sets for the pair loop:
    // the bitwise operations run on raw words without per-call overhead.
    size_t n_words = (m + 63) / 64;
    std::vector<uint64_t> zw(n_rays * n_words, 0);
    std::vector<size_t> zero_count(n_rays);
    std::vector<std::vector<int>> row_tight(m);
    std::vector<int> Rpos, Rneg;
    for (size_t i = 0; i < n_rays; i++) {
      zero_count[i] = rays[i].zero.count();
      boost::dynamic_bitset<>::size_type pos = rays[i].zero.find_first();
      while (pos != boost::dynamic_bitset<>::npos) {
        row_tight[pos].push_back(i);
        zw[i * n_words + pos / 64] |= uint64_t(1) << (pos % 64);
        pos = rays[i].zero.find_next(pos);
      }
      if (rays[i].val > 0) {
        Rpos.push_back(i);
      } else {
        if (rays[i].val < 0) {
          Rneg.push_back(i);
        }
      }
    }
    std::vector<uint64_t> common(n_words);
    for (auto &i1 : Rpos) {
      uint64_t const *w1 = zw.data() + i1 * n_words;
      for (auto &i2 : Rneg) {
        uint64_t const *w2 = zw.data() + i2 * n_words;
        // The combinatorial adjacency test of Fukuda-Prodon: the common
        // tight set must have at least d-2 elements and must not be
        // contained in the tight set of any other ray. The scan runs over
        // the rays tight on the least populated row of the common set:
        // any containing ray is tight on every one of its rows, so the
        // restriction loses nothing.
        size_t count_common = 0;
        for (size_t k = 0; k < n_words; k++) {
          common[k] = w1[k] & w2[k];
          count_common += __builtin_popcountll(common[k]);
        }
        if (count_common + 2 < static_cast<size_t>(d)) {
          continue;
        }
        size_t best_size = std::numeric_limits<size_t>::max();
        int best_row = -1;
        for (size_t k = 0; k < n_words; k++) {
          uint64_t w = common[k];
          while (w != 0) {
            int pos = 64 * k + __builtin_ctzll(w);
            w &= w - 1;
            if (row_tight[pos].size() < best_size) {
              best_size = row_tight[pos].size();
              best_row = pos;
            }
          }
        }
        bool adjacent = true;
        for (auto &i3 : row_tight[best_row]) {
          if (i3 == i1 || i3 == i2 || zero_count[i3] < count_common) {
            continue;
          }
          uint64_t const *w3 = zw.data() + i3 * n_words;
          bool subset = true;
          for (size_t k = 0; k < n_words; k++) {
            if ((common[k] & ~w3[k]) != 0) {
              subset = false;
              break;
            }
          }
          if (subset) {
            adjacent = false;
            break;
          }
        }
        if (!adjacent) {
          continue;
        }
        // The division-free combination w = (a.r1) r2 - (a.r2) r1, which
        // is tight on the new row and feasible for all processed rows.
        MyVector<T> w(d);
        for (int j = 0; j < d; j++) {
          w(j) = rays[i1].val * rays[i2].coord(j) -
                 rays[i2].val * rays[i1].coord(j);
        }
        NormalizeRayContent(w);
        // The exact tight set, evaluated over the processed rows.
        Face zero(m);
        boost::dynamic_bitset<>::size_type kRow = processed.find_first();
        while (kRow != boost::dynamic_bitset<>::npos) {
          T eSum(0);
          for (int j = 0; j < d; j++) {
            eSum += EXT(kRow, j) * w(j);
          }
          if (eSum == 0) {
            zero[kRow] = 1;
          }
#ifdef SANITY_CHECK_DOUBLE_DESCRIPTION
          if (eSum < 0) {
            std::cerr << "DOUBLEDESC: The new ray is infeasible at the "
                         "processed row "
                      << kRow << " which contradicts its construction\n";
            throw TerminalException{1};
          }
#endif
          kRow = processed.find_next(kRow);
        }
        new_rays.push_back({std::move(w), std::move(zero), T(0)});
      }
    }
    // Remove the infeasible rays and append the new ones.
    size_t pos = 0;
    for (size_t i1 = 0; i1 < n_rays; i1++) {
      if (!(rays[i1].val < 0)) {
        if (pos != i1) {
          rays[pos] = std::move(rays[i1]);
        }
        pos++;
      }
    }
    rays.resize(pos);
    for (auto &ray : new_rays) {
      rays.push_back(std::move(ray));
    }
#ifdef DEBUG_DOUBLE_DESCRIPTION
    os << "DOUBLEDESC: row " << iRow << " inserted, n_rays=" << rays.size()
       << "\n";
#endif
  }
#ifdef SANITY_CHECK_DOUBLE_DESCRIPTION
  // Each output ray must be feasible everywhere with the exact tight set
  // and the tight set must have rank d-1 (a facet).
  for (auto &ray : rays) {
    std::vector<int> tight;
    for (int i = 0; i < m; i++) {
      T eSum(0);
      for (int j = 0; j < d; j++) {
        eSum += EXT(i, j) * ray.coord(j);
      }
      if (eSum < 0) {
        std::cerr << "DOUBLEDESC: An output ray is infeasible at row " << i
                  << "\n";
        throw TerminalException{1};
      }
      if ((eSum == 0) != (ray.zero[i] == 1)) {
        std::cerr << "DOUBLEDESC: The tight set of an output ray is "
                     "incorrect at row "
                  << i << "\n";
        throw TerminalException{1};
      }
      if (eSum == 0) {
        tight.push_back(i);
      }
    }
    std::vector<int> sel = SelectIndependentRows(EXT, tight, d);
    if (static_cast<int>(sel.size()) != d - 1) {
      std::cerr << "DOUBLEDESC: An output ray has a tight set of rank "
                << sel.size() << " instead of " << (d - 1) << "\n";
      throw TerminalException{1};
    }
  }
#endif
#ifdef TIMINGS_DOUBLE_DESCRIPTION
  os << "|DOUBLEDESC: DoubleDescription_Kernel|=" << time << "\n";
#endif
  return rays;
}

// Assemble the public output from the kernel rays: the facet vectors are
// embedded back into the unselected columns with zeros (yielding valid
// facet inequalities of the original cone) and converted to T; the tight
// sets are the facet incidences and need no conversion.
template <typename T, typename Twork>
std::pair<MyMatrix<T>, vectface>
BuildOutput(std::vector<RayEntry<Twork>> const &rays,
            std::vector<int> const &colsel, int const &m, int const &d) {
  int n_facet = rays.size();
  int d_red = colsel.size();
  MyMatrix<T> FAC = ZeroMatrix<T>(n_facet, d);
  vectface vf(m);
  for (int i_facet = 0; i_facet < n_facet; i_facet++) {
    for (int j = 0; j < d_red; j++) {
      FAC(i_facet, colsel[j]) =
          UniversalScalarConversion<T, Twork>(rays[i_facet].coord(j));
    }
    vf.push_back(rays[i_facet].zero);
  }
  return {std::move(FAC), std::move(vf)};
}

/*
  The general entry point. For a field type T the rows are scaled once to
  the underlying ring and the kernel runs there; for a ring type the
  kernel runs directly. A column reduction is applied first when the
  matrix does not have full column rank: the incidences are not affected.
 */
template <typename T>
std::pair<MyMatrix<T>, vectface> DualDescription_pair(MyMatrix<T> const &EXT,
                                                      std::ostream &os) {
  int m = EXT.rows();
  int d = EXT.cols();
  if constexpr (is_ring_field<T>::value && !std::is_floating_point_v<T>) {
    using Tring = typename underlying_ring<T>::ring_type;
    MyMatrix<Tring> EXTring(m, d);
    for (int iRow = 0; iRow < m; iRow++) {
      FractionVectorRing<T> fr =
          RemoveFractionVectorPlusCoeffRing(GetMatrixRow(EXT, iRow));
      AssignMatrixRow(EXTring, iRow, fr.TheVect);
    }
    std::vector<int> colsel = SelectIndependentColumns(EXTring);
    if (static_cast<int>(colsel.size()) < d) {
      MyMatrix<Tring> EXTred = SelectColumn(EXTring, colsel);
      std::vector<RayEntry<Tring>> rays = DoubleDescription_Kernel(EXTred, os);
      return BuildOutput<T, Tring>(rays, colsel, m, d);
    }
    std::vector<RayEntry<Tring>> rays = DoubleDescription_Kernel(EXTring, os);
    return BuildOutput<T, Tring>(rays, colsel, m, d);
  } else {
    std::vector<int> colsel = SelectIndependentColumns(EXT);
    if (static_cast<int>(colsel.size()) < d) {
      MyMatrix<T> EXTred = SelectColumn(EXT, colsel);
      std::vector<RayEntry<T>> rays = DoubleDescription_Kernel(EXTred, os);
      return BuildOutput<T, T>(rays, colsel, m, d);
    }
    std::vector<RayEntry<T>> rays = DoubleDescription_Kernel(EXT, os);
    return BuildOutput<T, T>(rays, colsel, m, d);
  }
}

// The facet incidences of the cone generated by the rows of EXT: one Face
// per facet listing the rows lying on it. Matches the conventions of
// cdd::DualDescription_incd and lrs::DualDescription_incd.
template <typename T>
vectface DualDescription_incd(MyMatrix<T> const &EXT, std::ostream &os) {
  return DualDescription_pair(EXT, os).second;
}

// The facet vectors of the cone generated by the rows of EXT, i.e. the
// extreme rays of { x : EXT x >= 0 }. Matches cdd::DualDescription.
template <typename T>
MyMatrix<T> DualDescription(MyMatrix<T> const &EXT, std::ostream &os) {
  return DualDescription_pair(EXT, os).first;
}

// clang-format off
}  // namespace double_desc
#endif  // SRC_POLY_POLY_DOUBLE_DESCRIPTION_H_
// clang-format on
