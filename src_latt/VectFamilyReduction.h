// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_VECTFAMILYREDUCTION_H_
#define SRC_LATT_VECTFAMILYREDUCTION_H_
// clang-format off
#include "ClassicLLL.h"
#include "LatticeReduction.h"
#include "MAT_Matrix.h"
#include "norms.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_VECT_FAMILY_REDUCTION
#endif

#ifdef DISABLE_DEBUG_VECT_FAMILY_REDUCTION
#undef DEBUG_VECT_FAMILY_REDUCTION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_VECT_FAMILY_REDUCTION
#endif

#ifdef TIMINGS
#define TIMINGS_VECT_FAMILY_REDUCTION
#endif

/*
  Reduction of a family of vectors, that is a change of coordinates in the
  ambient space making the coefficients of the family small.

  This is not the same problem as reducing a lattice basis, and the difference
  decides which reduction to use. What a consumer of the output cares about is
  the size of the integers it will have to compute with; for the dual
  description in particular the relevant quantity is the Hadamard estimate
  sqr_estimate_facet_coefficients of norms.h, which bounds the coefficients of
  the facets that will be produced. That quantity depends on ALL the vectors,
  symmetrically. A reduction aiming at one short vector, which is what LLL
  does, optimises the wrong thing here; Seysen's measure, which penalises any
  basis vector being long and is invariant under exchanging the lattice with
  its dual, is a much closer match to it.

  So rather than pick a reducer on general grounds, the "best" method runs the
  candidates and keeps whichever actually minimises the estimate, the input
  itself being among the candidates so that the result is never worse than
  what was handed in. The cost is a small multiple of one reduction and is
  negligible against a dual description.

  This header sits in src_latt and not beside the other reducers in
  src_isotropy because it uses all of them, and BKZ and slide reduction need
  the shortest-vector enumerator of Shvec_exact.h. The directory dependency
  runs src_latt -> src_isotropy; a consumer of everything belongs at the top of
  that stack.
 */

/*
  The measure a reduction is judged by here. Lower is better. The facet
  estimate is the primary criterion, being what the dual description will pay;
  the L1 norm breaks its ties, which are common on small families.
 */
template <typename T> struct VectFamilyQuality {
  T facet_estimate;
  T l1_norm;
  bool operator<(VectFamilyQuality<T> const &other) const {
    if (facet_estimate != other.facet_estimate) {
      return facet_estimate < other.facet_estimate;
    }
    return l1_norm < other.l1_norm;
  }
};

template <typename T>
VectFamilyQuality<T> ComputeVectFamilyQuality(MyMatrix<T> const &M) {
  return {sqr_estimate_facet_coefficients(M), L1_norm_mat(M)};
}

/*
  The methods that reduce through a single Gram-matrix reduction. "best" is
  handled separately, being a search over these.
 */
inline std::vector<std::string> VectFamilyReductionSingleMethods() {
  return LatticeReductionSingleMethods();
}

inline std::vector<std::string> VectFamilyReductionMethods() {
  std::vector<std::string> methods = VectFamilyReductionSingleMethods();
  methods.push_back("best");
  return methods;
}

template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>>
ReduceVectorFamilySingle(MyMatrix<T> const &M, std::string const &method,
                         std::ostream &os) {
  using Tint = typename underlying_ring<T>::ring_type;
  auto f_reduce = [&](MyMatrix<T> const &G,
                      std::ostream &os_i) -> LLLreduction<T, Tint> {
    return LatticeReducedGeneral<T, Tint>(G, method, os_i);
  };
  return ReduceVectorFamilyKernel(M, f_reduce, os);
}

/*
  Every reduction here needs the column Gram matrix to be positive definite,
  which asks that the family span the ambient space. Checked once and reported
  in those terms: a rank deficient family is a modelling error upstream, not a
  numerical accident, and the message from a singular elimination deep inside
  a reducer would not say so.
 */
template <typename T>
void CheckVectFamilyFullColumnRank(MyMatrix<T> const &M) {
  int nbCol = M.cols();
  int rnk = RankMat(M);
  if (rnk != nbCol) {
    std::cerr << "VECT_FAMILY_REDUCTION: the family has rank " << rnk
              << " in dimension " << nbCol
              << ", so the Gram matrix of the columns is singular and no "
                 "reduction of it is defined. Restrict the family to the "
                 "subspace it spans first.\n";
    throw TerminalException{1};
  }
}

/*
  The best-of search. Returns the winning method as well, so that a caller
  running many families can see which reduction is earning its place and
  eventually stop searching.
 */
template <typename T> struct VectFamilyReductionResult {
  MyMatrix<T> Mred;
  MyMatrix<T> Pmat;
  std::string method;
  VectFamilyQuality<T> quality;
};

template <typename T>
VectFamilyReductionResult<T>
ReduceVectorFamilyBest(MyMatrix<T> const &M, std::ostream &os) {
#ifdef TIMINGS_VECT_FAMILY_REDUCTION
  MicrosecondTime time;
#endif
  // The identity is a candidate, so the search can never return something
  // worse than its input.
  MyMatrix<T> Pmat_id = IdentityMat<T>(M.cols());
  VectFamilyReductionResult<T> best{M, Pmat_id, "none",
                                    ComputeVectFamilyQuality(M)};
  for (auto &method : VectFamilyReductionSingleMethods()) {
    std::pair<MyMatrix<T>, MyMatrix<T>> pair =
        ReduceVectorFamilySingle(M, method, os);
    VectFamilyQuality<T> quality = ComputeVectFamilyQuality(pair.first);
#ifdef DEBUG_VECT_FAMILY_REDUCTION
    os << "VECT_FAMILY_REDUCTION: method " << method
       << " facet_estimate=" << quality.facet_estimate
       << " L1=" << quality.l1_norm << "\n";
#endif
    if (quality < best.quality) {
      best = {std::move(pair.first), std::move(pair.second), method,
              std::move(quality)};
    }
  }
#ifdef DEBUG_VECT_FAMILY_REDUCTION
  os << "VECT_FAMILY_REDUCTION: best method is " << best.method << "\n";
#endif
#ifdef TIMINGS_VECT_FAMILY_REDUCTION
  os << "VECT_FAMILY_REDUCTION: ReduceVectorFamilyBest took " << time << "\n";
#endif
  return best;
}

template <typename T>
VectFamilyReductionResult<T>
ReduceVectorFamilyGeneral(MyMatrix<T> const &M, std::string const &method,
                          std::ostream &os) {
  CheckVectFamilyFullColumnRank(M);
  if (method == "best") {
    return ReduceVectorFamilyBest(M, os);
  }
  // No check against a fixed list here: the parametrised methods make the set
  // of valid names infinite, and LatticeReducedGeneral rejects what it cannot
  // parse with a message naming the shapes it accepts.
  std::pair<MyMatrix<T>, MyMatrix<T>> pair =
      ReduceVectorFamilySingle(M, method, os);
  VectFamilyQuality<T> quality = ComputeVectFamilyQuality(pair.first);
  return {std::move(pair.first), std::move(pair.second), method,
          std::move(quality)};
}

// clang-format off
#endif  // SRC_LATT_VECTFAMILYREDUCTION_H_
// clang-format on
