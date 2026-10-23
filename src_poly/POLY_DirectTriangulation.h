// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DIRECTTRIANGULATION_H_
#define SRC_POLY_POLY_DIRECTTRIANGULATION_H_

// clang-format off
#include "POLY_DualDesc_reverse_search.h"
#include "POLY_DualDesc_beneath_and_beyond.h"
#include "POLY_RedundancyElimination.h"
#include <algorithm>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_DIRECT_TRIANGULATION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_DIRECT_TRIANGULATION
#endif

/*
  The triangulation of a polytope (or of a pointed cone) given by its
  vertices (the rows of EXT, in homogeneous coordinates), with a choice of
  method, as POLY_DirectDualDesc.h does for the dual description:
  --- lrs: the reverse search triangulation (rev_search::GetTriangulationDet_f),
      a lexicographic triangulation enumerated with bounded memory;
  --- beneath_beyond: the placing triangulation of the beneath-and-beyond
      method (beneath_beyond::GetTriangulationDet_f).
  The two methods give in general different triangulations of the same
  polytope. Each simplex is returned as the sorted list of the indices of its
  vertices (rows of EXT), with the signed determinant of
  SelectRow(EXT, simplex) for this order.
 */
enum class TriangulationProgram {
  lrs,
  beneath_beyond,
};

template <>
struct std::formatter<TriangulationProgram>
    : std::formatter<std::string_view> {
  static constexpr std::string_view name(TriangulationProgram prog) {
    switch (prog) {
    case TriangulationProgram::lrs:
      return "lrs";
    case TriangulationProgram::beneath_beyond:
      return "beneath_beyond";
    }
    return "unknown";
  }
  template <typename FormatContext>
  auto format(TriangulationProgram prog, FormatContext &ctx) const {
    return std::formatter<std::string_view>::format(name(prog), ctx);
  }
};

// The same names as for the dual description programs.
inline std::optional<TriangulationProgram>
triangulation_program_from_string_opt(std::string const &prog) {
  if (prog == "lrs" || prog == "rs" || prog == "reverse_search")
    return TriangulationProgram::lrs;
  if (prog == "beneath_beyond" || prog == "bb")
    return TriangulationProgram::beneath_beyond;
  return {};
}

inline TriangulationProgram
triangulation_program_from_string(std::string const &prog) {
  std::optional<TriangulationProgram> opt =
      triangulation_program_from_string_opt(prog);
  if (!opt) {
    std::cerr << "TRIG: ERROR: unknown triangulation program prog=" << prog
              << ", the possibilities are lrs (rs, reverse_search) and "
                 "beneath_beyond (bb)\n";
    throw TerminalException{1};
  }
  return *opt;
}

// Sorts the simplex, returning the sign of the permutation applied: the
// determinant of the sorted rows is that sign times the one of the rows in
// the original order.
inline int sort_simplex_sign(std::vector<int> &simplex) {
  int sign = 1;
  size_t len = simplex.size();
  for (size_t i = 1; i < len; i++) {
    int val = simplex[i];
    size_t j = i;
    while (j > 0 && simplex[j - 1] > val) {
      simplex[j] = simplex[j - 1];
      j--;
      sign = -sign;
    }
    simplex[j] = val;
  }
  return sign;
}

// Beneath-and-beyond pivots and decides orientations by exact zero and sign
// tests, so it needs exact arithmetic: with double the rounding residues are
// taken as nonzero, and the row selection goes out of bounds.
[[noreturn]] inline void error_beneath_beyond_floating_point() {
  std::cerr << "TRIG: the beneath_beyond method requires an exact arithmetic, "
               "use lrs for floating point\n";
  throw TerminalException{1};
}

/*
  The input of the triangulation dispatchers has to be made of distinct
  extreme rays: a triangulation of a repeated or a non-extreme row does not
  make sense. Exact arithmetic only, as for the check of the determinants.
 */
template <typename T>
void check_direct_triangulation_input(
    [[maybe_unused]] MyMatrix<T> const &EXT,
    [[maybe_unused]] std::string const &context,
    [[maybe_unused]] std::ostream &os) {
  if constexpr (!std::is_floating_point_v<T>) {
#ifdef SANITY_CHECK_DIRECT_TRIANGULATION
    check_distinct_nonzero_rays(EXT, context);
#endif
#ifdef SANITY_CHECK_DIRECT_TRIANGULATION_EXTENSIVE
    std::cerr << "TRIG: SANITY_CHECK_DIRECT_TRIANGULATION_EXTENSIVE, checking "
                 "that the rows are extreme rays in "
              << context << "\n";
    check_all_extreme_rays(EXT, context, os);
    std::cerr << "TRIG: SANITY_CHECK_DIRECT_TRIANGULATION_EXTENSIVE, done\n";
#endif
  }
}

// The triangulation streamed: f_trig_det(simplex, det) is called for each
// simplex, simplex the sorted indices of its vertices and det the signed
// determinant of SelectRow(EXT, simplex).
template <typename T, typename Ftrig_det>
void DirectTriangulationDet_f(MyMatrix<T> const &EXT,
                              TriangulationProgram prog, std::ostream &os,
                              Ftrig_det f_trig_det) {
  check_direct_triangulation_input(EXT, "DirectTriangulationDet_f", os);
#ifdef DEBUG_DIRECT_TRIANGULATION
  os << "TRIG: DirectTriangulationDet_f, prog=" << std::format("{}", prog)
     << " |EXT|=" << EXT.rows() << " / " << EXT.cols() << "\n";
#endif
  std::vector<int> sorted;
  auto f_sorted = [&](std::vector<int> const &simplex, T const &det) -> void {
    sorted = simplex;
    int sign = sort_simplex_sign(sorted);
#ifdef SANITY_CHECK_DIRECT_TRIANGULATION
    // Exact arithmetic only: in floating point the two computations differ
    // by rounding.
    if constexpr (!std::is_floating_point_v<T>) {
      T det_check = DeterminantMat(SelectRow(EXT, sorted));
      T det_sorted = (sign > 0) ? det : -det;
      if (det_check != det_sorted) {
        std::cerr << "TRIG: the determinant " << det_sorted
                  << " of the simplex is not the one of its vertices "
                  << det_check << ", prog=" << std::format("{}", prog) << "\n";
        throw TerminalException{1};
      }
    }
#endif
    if (sign > 0) {
      f_trig_det(sorted, det);
    } else {
      T det_neg = -det;
      f_trig_det(sorted, det_neg);
    }
  };
  switch (prog) {
  case TriangulationProgram::lrs:
    rev_search::GetTriangulationDet_f<T>(EXT, f_sorted);
    return;
  case TriangulationProgram::beneath_beyond:
    if constexpr (std::is_floating_point_v<T>) {
      error_beneath_beyond_floating_point();
    } else {
      beneath_beyond::GetTriangulationDet_f<T>(EXT, os, f_sorted);
    }
    return;
  }
  std::cerr << "TRIG: ERROR in DirectTriangulationDet_f: no handler for "
            << std::format("{}", prog) << "\n";
  throw TerminalException{1};
}

template <typename T>
std::vector<std::pair<std::vector<int>, T>>
DirectTriangulationDet(MyMatrix<T> const &EXT, TriangulationProgram prog,
                       std::ostream &os) {
  std::vector<std::pair<std::vector<int>, T>> l_trig;
  auto f_trig_det = [&](std::vector<int> const &simplex, T const &det) -> void {
    l_trig.push_back({simplex, det});
  };
  DirectTriangulationDet_f<T>(EXT, prog, os, f_trig_det);
  return l_trig;
}

template <typename T>
std::vector<std::vector<int>> DirectTriangulation(MyMatrix<T> const &EXT,
                                                  TriangulationProgram prog,
                                                  std::ostream &os) {
  std::vector<std::vector<int>> l_trig;
  auto f_trig_det = [&](std::vector<int> const &simplex,
                        [[maybe_unused]] T const &det) -> void {
    l_trig.push_back(simplex);
  };
  DirectTriangulationDet_f<T>(EXT, prog, os, f_trig_det);
  return l_trig;
}

template <typename T>
std::vector<std::pair<std::vector<int>, T>>
DirectTriangulationDet(MyMatrix<T> const &EXT, std::string const &prog,
                       std::ostream &os) {
  return DirectTriangulationDet<T>(EXT, triangulation_program_from_string(prog),
                                   os);
}

template <typename T>
std::vector<std::vector<int>> DirectTriangulation(MyMatrix<T> const &EXT,
                                                  std::string const &prog,
                                                  std::ostream &os) {
  return DirectTriangulation<T>(EXT, triangulation_program_from_string(prog),
                                os);
}

// GAP-format string of a triangulation: the simplices with their vertex
// indices shifted to 1-based.
inline std::string
StringTriangulationGAP(std::vector<std::vector<int>> const &l_trig) {
  std::string ret = "[";
  for (size_t i = 0; i < l_trig.size(); i++) {
    if (i > 0)
      ret += ",";
    ret += "[";
    for (size_t j = 0; j < l_trig[i].size(); j++) {
      if (j > 0)
        ret += ",";
      ret += std::to_string(l_trig[i][j] + 1);
    }
    ret += "]";
  }
  ret += "]";
  return ret;
}

// The volume of a polytope given by its vertices (rows of EXT, first
// coordinate 1): the sum of the absolute values of the determinants of the
// simplices of a triangulation, divided by dim factorial. The value is in the
// field of fractions of T, since the volume of a lattice polytope is in
// general not an integer.
template <typename T>
typename overlying_field<T>::field_type
DirectVolumePolytope(MyMatrix<T> const &EXT, TriangulationProgram prog,
                     std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  T sum_det(0);
  auto f_trig_det = [&]([[maybe_unused]] std::vector<int> const &simplex,
                        T const &det) -> void { sum_det += T_abs(det); };
  DirectTriangulationDet_f<T>(EXT, prog, os, f_trig_det);
  int dim = EXT.cols() - 1;
  Tfield fact(1);
  for (int u = 2; u <= dim; u++)
    fact *= Tfield(u);
  return UniversalScalarConversion<Tfield, T>(sum_det) / fact;
}

// The triangulation (simplices as sorted vertex indices) together with the
// facets, as incidences over the rows of EXT. Reverse search obtains both
// from one enumeration. Beneath-and-beyond computes them in two runs: its
// placing triangulation does not keep track of the facets.
template <typename T>
std::pair<std::vector<std::vector<int>>, vectface>
DirectTriangulationFacet(MyMatrix<T> const &EXT, TriangulationProgram prog,
                         std::ostream &os) {
  check_direct_triangulation_input(EXT, "DirectTriangulationFacet", os);
#ifdef DEBUG_DIRECT_TRIANGULATION
  os << "TRIG: DirectTriangulationFacet, prog=" << std::format("{}", prog)
     << " |EXT|=" << EXT.rows() << " / " << EXT.cols() << "\n";
#endif
  switch (prog) {
  case TriangulationProgram::lrs: {
    std::pair<std::vector<std::vector<int>>, vectface> pair =
        rev_search::GetTriangulationFacet(EXT);
    for (auto &simplex : pair.first)
      std::sort(simplex.begin(), simplex.end());
    return pair;
  }
  case TriangulationProgram::beneath_beyond: {
    if constexpr (std::is_floating_point_v<T>) {
      error_beneath_beyond_floating_point();
    } else {
      std::vector<std::vector<int>> l_trig = DirectTriangulation(EXT, prog, os);
      vectface vf = beneath_beyond::DualDescription_incd(EXT, os);
      return {std::move(l_trig), std::move(vf)};
    }
  }
  }
  std::cerr << "TRIG: ERROR in DirectTriangulationFacet: no handler for "
            << std::format("{}", prog) << "\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_POLY_POLY_DIRECTTRIANGULATION_H_
// clang-format on
