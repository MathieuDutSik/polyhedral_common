// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_POLYTOPEEQUISTABMETHODS_H_
#define SRC_GROUP_POLYTOPEEQUISTABMETHODS_H_

// clang-format off
#include "PolytopeEquiStab.h"
#include "PolytopePleskenSouvignier.h"
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  The automorphism group, the canonical form and the isomorphism of
  polytopes (linear equivalence of the rows) by a chosen method:
  * standard: the subset method above THRESHOLD_USE_SUBSET_SCHEME_STAB
    (automorphism) or THRESHOLD_USE_SUBSET_SCHEME_CANONIC (canonical form)
    rows, the graph of the whole weight matrix below.
  * graph: the graph of the whole weight matrix, quadratic in the number of
    rows.
  * subset: the subset method (with the individualization for the
    automorphism group, see f_for_stab_heuristic).
  * plesken_souvignier: the Plesken-Souvignier backtrack, see
    PolytopePleskenSouvignier.h. Only for rational coordinates.
  The isomorphism is computed from the canonical orderings of both sides,
  which have to come from the same method.
 */
enum class LinPolytopeMethod { standard, graph, subset, plesken_souvignier };

inline std::optional<LinPolytopeMethod>
LinPolytopeMethodFromString(std::string const &method) {
  if (method == "standard") {
    return LinPolytopeMethod::standard;
  }
  if (method == "graph") {
    return LinPolytopeMethod::graph;
  }
  if (method == "subset") {
    return LinPolytopeMethod::subset;
  }
  if (method == "plesken_souvignier") {
    return LinPolytopeMethod::plesken_souvignier;
  }
  return {};
}

inline void LinPolytopeMethodPrintDescriptions(std::ostream &os) {
  os << "standard           : subset method for many rows, graph of the "
        "weight matrix otherwise (default)\n";
  os << "graph              : graph of the whole weight matrix\n";
  os << "subset             : subset method\n";
  os << "plesken_souvignier : Plesken-Souvignier backtrack (rational "
        "coordinates only)\n";
}

// The Plesken-Souvignier method needs integral coordinates, which a real
// algebraic field does not give.
template <typename T> void LinPolytopeMethodCheckPleskenSouvignier() {
  if constexpr (!is_implementation_of_Q<
                    typename overlying_field<T>::field_type>::value) {
    std::cerr << "LinPolytopeMethod: plesken_souvignier is only available "
                 "for rational coordinates\n";
    throw TerminalException{1};
  }
}

template <typename T, typename Tgroup>
Tgroup LinPolytope_Automorphism_Method(MyMatrix<T> const &EXT,
                                       LinPolytopeMethod method,
                                       std::ostream &os) {
  if (method == LinPolytopeMethod::standard) {
    return LinPolytope_Automorphism<T, Tgroup>(EXT, os);
  }
  if (method == LinPolytopeMethod::plesken_souvignier) {
    LinPolytopeMethodCheckPleskenSouvignier<T>();
    if constexpr (is_implementation_of_Q<
                      typename overlying_field<T>::field_type>::value) {
      return LinPolytope_Automorphism_PleskenSouvignier<T, Tgroup>(EXT, os);
    }
  }
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  MyMatrix<T> Qmat = GetQmatrix(EXTred, os);
  if (method == LinPolytopeMethod::subset) {
    return GroupFromStabGeneratorsOrder<Tgroup>(
        LinPolytope_Automorphism_GramMat_GensOrder_Heuristic<T, Tgroup>(
            EXTred, Qmat, os),
        EXT.rows());
  }
  size_t max_poss_val =
      weightmatrix_get_nb(IsSymmetricMatrix(Qmat), EXTred.rows());
  auto f_dispatch = [&]<typename Tidx_value>() -> Tgroup {
    return GroupFromStabGeneratorsOrder<Tgroup>(
        LinPolytope_Automorphism_GramMat_GensOrder_Kernel<T, Tgroup,
                                                          Tidx_value>(
            EXTred, Qmat, os),
        EXT.rows());
  };
  return call_with_smallest_unsigned(
      max_poss_val, "LinPolytope_Automorphism_Method", f_dispatch);
}

template <typename T, typename Tgroup>
std::vector<typename Tgroup::Telt::Tidx>
LinPolytope_CanonicOrdering_Method(MyMatrix<T> const &EXT,
                                   LinPolytopeMethod method,
                                   std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  if (method == LinPolytopeMethod::plesken_souvignier) {
    LinPolytopeMethodCheckPleskenSouvignier<T>();
    if constexpr (is_implementation_of_Q<
                      typename overlying_field<T>::field_type>::value) {
      return LinPolytope_CanonicOrdering_PleskenSouvignier<T, Tgroup>(EXT, os)
          .first;
    }
  }
  size_t threshold = THRESHOLD_USE_SUBSET_SCHEME_CANONIC;
  if (method == LinPolytopeMethod::graph) {
    threshold = std::numeric_limits<size_t>::max();
  }
  if (method == LinPolytopeMethod::subset) {
    threshold = 0;
  }
  return LinPolytope_CanonicOrdering<T, Tidx>(EXT, threshold, os);
}

// The canonical form of the column reduction of EXT: CanonicalizeOrderedMatrix
// needs a matrix of full column rank, and the linear class of the rows does
// not depend on the reduction.
template <typename T, typename Tgroup>
MyMatrix<T> LinPolytope_CanonicForm_Method(MyMatrix<T> const &EXT,
                                           LinPolytopeMethod method,
                                           std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  std::vector<Tidx> CanonicOrd =
      LinPolytope_CanonicOrdering_Method<T, Tgroup>(EXTred, method, os);
  MyMatrix<T> EXTreord(EXTred.rows(), EXTred.cols());
  for (int iRow = 0; iRow < EXTred.rows(); iRow++) {
    EXTreord.row(iRow) = EXTred.row(CanonicOrd[iRow]);
  }
  return CanonicalizeOrderedMatrix(EXTreord);
}

template <typename T, typename Tgroup>
std::optional<std::vector<typename Tgroup::Telt::Tidx>>
LinPolytope_Isomorphism_Method(MyMatrix<T> const &EXT1,
                               MyMatrix<T> const &EXT2,
                               LinPolytopeMethod method, std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  using Tfield = typename overlying_field<T>::field_type;
  if (EXT1.rows() != EXT2.rows()) {
    return {};
  }
  // On the column reductions, which IsomorphismFromCanonicReord needs of
  // full column rank; the linear equivalences of the rows are the same.
  MyMatrix<T> EXTred1 = ColumnReduction(EXT1);
  MyMatrix<T> EXTred2 = ColumnReduction(EXT2);
  if (EXTred1.cols() != EXTred2.cols()) {
    return {};
  }
  std::vector<Tidx> CanonicReord1 =
      LinPolytope_CanonicOrdering_Method<T, Tgroup>(EXTred1, method, os);
  std::vector<Tidx> CanonicReord2 =
      LinPolytope_CanonicOrdering_Method<T, Tgroup>(EXTred2, method, os);
  std::optional<std::pair<std::vector<Tidx>, MyMatrix<Tfield>>> IsoInfo =
      IsomorphismFromCanonicReord<T, Tfield, Tidx>(
          EXTred1, EXTred2, CanonicReord1, CanonicReord2, os);
  if (!IsoInfo) {
    return {};
  }
  return IsoInfo->first;
}

// clang-format off
#endif  // SRC_GROUP_POLYTOPEEQUISTABMETHODS_H_
// clang-format on
