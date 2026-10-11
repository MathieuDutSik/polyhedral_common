// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DIRECTTRIANGULATION_H_
#define SRC_POLY_POLY_DIRECTTRIANGULATION_H_

// clang-format off
#include "POLY_DualDesc_reverse_search.h"
#include "POLY_DualDesc_beneath_and_beyond.h"
#include <format>
#include <optional>
#include <string>
#include <string_view>
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
      method (POLY_DualDescription_BeneathBeyondTriangulationDet_f).
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

// The triangulation streamed: f_trig_det(simplex, det) is called for each
// simplex, simplex the sorted indices of its vertices and det the signed
// determinant of SelectRow(EXT, simplex).
template <typename T, typename Ftrig_det>
void DirectTriangulationDet_f(MyMatrix<T> const &EXT,
                              TriangulationProgram prog, std::ostream &os,
                              Ftrig_det f_trig_det) {
#ifdef DEBUG_DIRECT_TRIANGULATION
  os << "TRIG: DirectTriangulationDet_f, prog=" << std::format("{}", prog)
     << " |EXT|=" << EXT.rows() << " / " << EXT.cols() << "\n";
#endif
  std::vector<int> sorted;
  auto f_sorted = [&](std::vector<int> const &simplex, T const &det) -> void {
    sorted = simplex;
    int sign = sort_simplex_sign(sorted);
#ifdef SANITY_CHECK_DIRECT_TRIANGULATION
    T det_check = DeterminantMat(SelectRow(EXT, sorted));
    T det_sorted = (sign > 0) ? det : -det;
    if (det_check != det_sorted) {
      std::cerr << "TRIG: the determinant " << det_sorted
                << " of the simplex is not the one of its vertices "
                << det_check << ", prog=" << std::format("{}", prog) << "\n";
      throw TerminalException{1};
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
    POLY_DualDescription_BeneathBeyondTriangulationDet_f<T>(EXT, os, f_sorted);
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

// clang-format off
#endif  // SRC_POLY_POLY_DIRECTTRIANGULATION_H_
// clang-format on
