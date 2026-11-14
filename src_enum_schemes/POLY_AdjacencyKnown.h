// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ENUM_SCHEMES_POLY_ADJACENCYKNOWN_H_
#define SRC_ENUM_SCHEMES_POLY_ADJACENCYKNOWN_H_

// clang-format off
#include "Boost_bitset.h"
#include "MAT_MatrixFund.h"
#include <optional>
#include <unordered_map>
#include <vector>
// clang-format on

/*
  Helpers for the enumerations whose adjacencies cross the facets of a
  polytope (Delaunay polytopes, perfect domains) and which provide the
  adjacencies known from the other side (AdjacencyKnown of
  POLY_AdjacencyScheme.h).

  Such an adjacency is given by a facet of the object. The adjacency known
  from the other side is the facet shared with the object it comes from, so
  it is found from the vertices of that facet (FaceOfRows), and it replaces
  the computation of the adjacency across the facets of its orbit
  (MatchKnownFacets).
 */

// The rows of EXT among the rows of M, as a face of EXT. Every row of M has
// to be a row of EXT, otherwise there is no face.
template <typename T>
std::optional<Face> FaceOfRows(MyMatrix<T> const &EXT, MyMatrix<T> const &M) {
  ContainerMatrix<T> cont(EXT);
  Face f(EXT.rows());
  int n_row = M.rows();
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = GetMatrixRow(M, i_row);
    std::optional<size_t> opt = cont.GetIdx_v(V);
    if (!opt) {
      return {};
    }
    f[*opt] = 1;
  }
  return f;
}

/*
  l_rep are the representatives of the orbits of the facets under GRP, and
  l_known the facets of the known adjacencies. For each representative, the
  index in l_known of a facet in its orbit, if there is one: the known
  adjacency then replaces the computation for that orbit. A known facet is
  used at most once. The orbits are identified by canonical images.
 */
template <typename Tgroup>
std::vector<std::optional<size_t>>
MatchKnownFacets(Tgroup const &GRP, vectface const &l_rep,
                 std::vector<Face> const &l_known) {
  std::unordered_map<Face, size_t> map_known;
  size_t n_known = l_known.size();
  for (size_t u = 0; u < n_known; u++) {
    map_known.emplace(GRP.OptCanonicalImage(l_known[u]), u);
  }
  std::vector<std::optional<size_t>> l_match;
  for (auto &eRep : l_rep) {
    if (map_known.empty()) {
      l_match.push_back({});
      continue;
    }
    auto iter = map_known.find(GRP.OptCanonicalImage(eRep));
    if (iter == map_known.end()) {
      l_match.push_back({});
    } else {
      l_match.push_back(iter->second);
      map_known.erase(iter);
    }
  }
  return l_match;
}

// clang-format off
#endif  // SRC_ENUM_SCHEMES_POLY_ADJACENCYKNOWN_H_
// clang-format on
