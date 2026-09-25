// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_BASIC_H_
#define SRC_POLYNORM_POLYNORM_BASIC_H_

// clang-format off
#include "GRP_GroupFct.h"
#include "PolytopeEquiStabInt.h"
#include "POLY_DirectDualDesc.h"
#include "POLY_Fundamental.h"
#include "POLY_PolytopeInt.h"
#include <algorithm>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_BASIC
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYNORM_BASIC
#endif

#ifdef TIMINGS
#define TIMINGS_POLYNORM_BASIC
#endif

/*
  Polyhedral norms.

  A polytope P in R^n with 0 in its interior defines the gauge
     ||x||_P = min{t >= 0 : x in tP} = max_j ell_j(x)
  where the facets of P are written as P = {x : ell_j(x) <= 1, j in [m]}.
  The gauge is a norm when P is centrally symmetric and only an asymmetric
  distance function otherwise; all of the code below works in the asymmetric
  case.

  The two lattice problems of this directory are, for the lattice Z^n:
  --- Packing: the largest alpha such that the translates alpha P + z,
      z in Z^n, have pairwise disjoint interiors.
  --- Covering: the smallest mu such that the translates mu P + z,
      z in Z^n, cover R^n.
  Both are invariant under translations of P, so P is recentered at the
  isobarycenter of its vertices, which puts 0 in the interior and makes the
  symmetry group of the problem the group of A in GL_n(Z) with PA = P
  (vectors are row vectors, matrices act on the right).

  Conventions of this file:
  --- EXT is the vertex matrix in homogeneous coordinates: row (1, v).
  --- Lmat is the m x n matrix of the linear forms ell_j, so that
      P = {x : Lmat x <= 1}. The gauge is the maximum of the rows on x.
  --- The difference body D = P - P is the unit ball of the symmetrization
      ||x||_D = max(||x||_P, ||-x||_P) only when P is symmetric; in general
      it is a genuine convex body of its own and LmatDiff holds its facets.
 */

// The facet normals of P as linear forms: from the vertex incidences of the
// facets, each facet inequality is normalized so that its constant term is
// 1, which is possible since 0 is interior. Row j of the returned matrix is
// ell_j with P = {x : ell_j(x) <= 1}.
template <typename T>
MyMatrix<T> PolyNorm_LinearForms(MyMatrix<T> const &EXT, vectface const &vf) {
  int n = EXT.cols() - 1;
  int m = vf.size();
  MyMatrix<T> Lmat(m, n);
  for (int j = 0; j < m; j++) {
    Face f = vf[j];
    MyVector<T> w = FindFacetInequality(EXT, f);
    // w(0) + sum_k w(k) x_k >= 0 on P, and 0 is interior so w(0) > 0.
#ifdef SANITY_CHECK_POLYNORM_BASIC
    if (w(0) <= 0) {
      std::cerr << "POLYNORM: The origin is not in the interior of the "
                << "polytope, facet j=" << j << " has w(0)=" << w(0) << "\n";
      throw TerminalException{1};
    }
#endif
    for (int k = 0; k < n; k++) {
      Lmat(j, k) = -w(k + 1) / w(0);
    }
  }
  return Lmat;
}

// The gauge ||x||_P = max_j ell_j(x). For x = 0 this is 0 and for x != 0 it
// is positive since 0 is interior.
template <typename T>
T PolyNorm_Gauge(MyMatrix<T> const &Lmat, MyVector<T> const &x) {
  int m = Lmat.rows();
  int n = Lmat.cols();
  T val(0);
  for (int j = 0; j < m; j++) {
    T scal(0);
    for (int k = 0; k < n; k++) {
      scal += Lmat(j, k) * x(k);
    }
    if (j == 0 || scal > val) {
      val = scal;
    }
  }
  return val;
}

// The index of a facet attaining the gauge, together with the gauge.
template <typename T>
std::pair<T, int> PolyNorm_GaugeArgmax(MyMatrix<T> const &Lmat,
                                       MyVector<T> const &x) {
  int m = Lmat.rows();
  int n = Lmat.cols();
  T val(0);
  int idx = -1;
  for (int j = 0; j < m; j++) {
    T scal(0);
    for (int k = 0; k < n; k++) {
      scal += Lmat(j, k) * x(k);
    }
    if (j == 0 || scal > val) {
      val = scal;
      idx = j;
    }
  }
  return {val, idx};
}

// The lattice points of t P for the polytope P = {x : Lmat x <= 1}: the
// inequalities are t - ell_j(x) >= 0.
template <typename T, typename Tint>
std::vector<MyVector<Tint>> PolyNorm_LatticePoints(MyMatrix<T> const &Lmat,
                                                   T const &t,
                                                   std::ostream &os) {
  int m = Lmat.rows();
  int n = Lmat.cols();
  MyMatrix<T> FAC(m, n + 1);
  for (int j = 0; j < m; j++) {
    FAC(j, 0) = t;
    for (int k = 0; k < n; k++) {
      FAC(j, k + 1) = -Lmat(j, k);
    }
  }
  return GetListIntegralPoint<T, Tint>(FAC, os);
}

// One symmetry of the problem: the matrix A in GL_n(Z) with PA = P and the
// permutation of the facets that it induces. The permutation of the
// vertices is the one that the group computation works with; the matrix is
// what acts on lattice points.
template <typename Tint> struct PolyNormSymm {
  MyMatrix<Tint> A;
  std::vector<int> facet_perm;
};

// The recentered polytope with everything the packing and covering
// computations need.
template <typename T, typename Tint, typename Tgroup> struct PolyNormData {
  int n;
  // The isobarycenter of the input vertices: the recentered polytope is
  // P_input - center.
  MyVector<T> center;
  // The recentered vertices, homogeneous coordinates.
  MyMatrix<T> EXT;
  // The facets of P as linear forms, P = {x : Lmat x <= 1}.
  MyMatrix<T> Lmat;
  // The vertex incidences of the facets, indexed as the rows of Lmat.
  vectface ListIncd;
  // The facets of the difference body D = P - P as linear forms.
  MyMatrix<T> LmatDiff;
  // The group {A in GL_n(Z) : PA = P} as permutations of the vertices.
  Tgroup GRP;
  // Its generators as matrices.
  std::vector<PolyNormSymm<Tint>> ListGen;
  // All its elements as matrices when the group is small enough (see
  // PolyNorm_MaxGroupSizeElements), which is what the canonical forms and
  // the stabilizer computations of the covering search use. Unset when the
  // group is too large: the search then only uses the orbits of the facets
  // computed from the generators.
  std::optional<std::vector<PolyNormSymm<Tint>>> ListElt;
};

// Above this order the elements of the symmetry group are not listed.
#ifndef POLYNORM_MAX_GROUP_SIZE_ELEMENTS
#define POLYNORM_MAX_GROUP_SIZE_ELEMENTS 50000
#endif
const size_t PolyNorm_MaxGroupSizeElements = POLYNORM_MAX_GROUP_SIZE_ELEMENTS;

template <typename T, typename Tint, typename Telt>
PolyNormSymm<Tint>
PolyNorm_SymmFromPerm(MyMatrix<T> const &EXT, vectface const &ListIncd,
                      std::unordered_map<Face, int> const &map_incd,
                      Telt const &ePerm) {
  int n = EXT.cols() - 1;
  MyMatrix<T> M = FindTransformation<T, Telt>(EXT, EXT, ePerm);
#ifdef SANITY_CHECK_POLYNORM_BASIC
  // The polytope is recentered at its isobarycenter, so the affine map has
  // no translation part.
  if (M(0, 0) != 1) {
    std::cerr << "POLYNORM: The transformation does not preserve the "
              << "homogeneous coordinate\n";
    throw TerminalException{1};
  }
  for (int k = 0; k < n; k++) {
    if (M(0, k + 1) != 0) {
      std::cerr << "POLYNORM: The transformation has a translation part, "
                << "the polytope should be centered\n";
      throw TerminalException{1};
    }
  }
#endif
  MyMatrix<T> A_T(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A_T(i, j) = M(i + 1, j + 1);
    }
  }
#ifdef SANITY_CHECK_POLYNORM_BASIC
  if (!IsIntegralMatrix(A_T)) {
    std::cerr << "POLYNORM: The transformation is not integral\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> A = UniversalMatrixConversion<Tint, T>(A_T);
  // The vertex k is mapped to the vertex ePerm(k), so the facet with
  // incidence f is mapped to the facet with incidence ePerm(f).
  size_t n_vert = EXT.rows();
  int m = ListIncd.size();
  std::vector<int> facet_perm(m);
  for (int j = 0; j < m; j++) {
    Face f = ListIncd[j];
    Face g(n_vert);
    boost::dynamic_bitset<>::size_type k = f.find_first();
    while (k != boost::dynamic_bitset<>::npos) {
      g[ePerm.at(k)] = 1;
      k = f.find_next(k);
    }
    auto iter = map_incd.find(g);
#ifdef SANITY_CHECK_POLYNORM_BASIC
    if (iter == map_incd.end()) {
      std::cerr << "POLYNORM: The image of a facet is not a facet\n";
      throw TerminalException{1};
    }
#endif
    facet_perm[j] = iter->second;
  }
  return {std::move(A), std::move(facet_perm)};
}

// The polytope is given by its vertices in homogeneous coordinates. It
// must be full dimensional.
template <typename T, typename Tint, typename Tgroup>
PolyNormData<T, Tint, Tgroup> PolyNorm_BuildData(MyMatrix<T> const &EXT_input,
                                                 std::ostream &os) {
  using Telt = typename Tgroup::Telt;
#ifdef TIMINGS_POLYNORM_BASIC
  MicrosecondTime time;
#endif
  int n_vert = EXT_input.rows();
  int n = EXT_input.cols() - 1;
  if (n_vert == 0) {
    std::cerr << "POLYNORM: The polytope has no vertices\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < n_vert; i++) {
    if (EXT_input(i, 0) != 1) {
      std::cerr << "POLYNORM: The vertices should be in homogeneous "
                << "coordinates with first entry 1, row i=" << i
                << " has " << EXT_input(i, 0) << "\n";
      throw TerminalException{1};
    }
  }
  if (RankMat(EXT_input) != n + 1) {
    std::cerr << "POLYNORM: The polytope should be full dimensional\n";
    throw TerminalException{1};
  }
  PolyNormData<T, Tint, Tgroup> data;
  data.n = n;
  // Only the distinct vertices are kept: the symmetry group is computed
  // from the vertex set (the automorphism code does not accept repeated
  // rows) and the isobarycenter must be that of the vertices, so that it is
  // fixed by the symmetries. A row is a vertex when the facets containing
  // it have normals spanning a hyperplane.
  MyMatrix<T> EXT_vert = EXT_input;
  {
    std::unordered_set<MyVector<T>> set_row;
    std::vector<int> ListDistinct;
    for (int i = 0; i < n_vert; i++) {
      MyVector<T> v = GetMatrixRow(EXT_input, i);
      if (set_row.count(v) == 0) {
        set_row.insert(v);
        ListDistinct.push_back(i);
      }
    }
    if (int(ListDistinct.size()) < n_vert) {
#ifdef DEBUG_POLYNORM_BASIC
      os << "POLYNORM: " << n_vert - ListDistinct.size()
         << " input rows are repeated and are dropped\n";
#endif
      EXT_vert = SelectRow(EXT_input, ListDistinct);
      n_vert = ListDistinct.size();
    }
  }
  {
    vectface vf_input = DirectFacetComputationIncidence(EXT_vert, "cdd", os);
    int m_input = vf_input.size();
    MyMatrix<T> FAC_input(m_input, n + 1);
    for (int j = 0; j < m_input; j++) {
      MyVector<T> w = FindFacetInequality(EXT_vert, vf_input[j]);
      for (int k = 0; k <= n; k++) {
        FAC_input(j, k) = w(k);
      }
    }
    std::vector<int> ListVert;
    for (int i = 0; i < n_vert; i++) {
      std::vector<int> ListFacet;
      for (int j = 0; j < m_input; j++) {
        if (vf_input[j][i] == 1) {
          ListFacet.push_back(j);
        }
      }
      MyMatrix<T> FAC_sel = SelectRow(FAC_input, ListFacet);
      if (RankMat(FAC_sel) == n) {
        ListVert.push_back(i);
      }
    }
    if (int(ListVert.size()) < n_vert) {
#ifdef DEBUG_POLYNORM_BASIC
      os << "POLYNORM: " << n_vert - ListVert.size()
         << " input rows are not vertices and are dropped\n";
#endif
      EXT_vert = SelectRow(EXT_vert, ListVert);
      n_vert = ListVert.size();
    }
  }
  // Recentering at the isobarycenter of the vertices.
  MyVector<T> center = ZeroVector<T>(n);
  for (int i = 0; i < n_vert; i++) {
    for (int k = 0; k < n; k++) {
      center(k) += EXT_vert(i, k + 1);
    }
  }
  center /= T(n_vert);
  data.center = center;
  MyMatrix<T> EXT(n_vert, n + 1);
  for (int i = 0; i < n_vert; i++) {
    EXT(i, 0) = 1;
    for (int k = 0; k < n; k++) {
      EXT(i, k + 1) = EXT_vert(i, k + 1) - center(k);
    }
  }
  data.EXT = EXT;
  // The facets of P.
  data.ListIncd = DirectFacetComputationIncidence(EXT, "cdd", os);
  data.Lmat = PolyNorm_LinearForms(EXT, data.ListIncd);
#ifdef TIMINGS_POLYNORM_BASIC
  os << "|POLYNORM: facets of P|=" << time << "\n";
#endif
  // The difference body D = P - P, whose vertices are among the
  // differences of vertices.
  std::unordered_set<MyVector<T>> set_diff;
  std::vector<MyVector<T>> list_diff;
  for (int i = 0; i < n_vert; i++) {
    for (int j = 0; j < n_vert; j++) {
      MyVector<T> v(n);
      for (int k = 0; k < n; k++) {
        v(k) = EXT(i, k + 1) - EXT(j, k + 1);
      }
      if (set_diff.count(v) == 0) {
        set_diff.insert(v);
        list_diff.push_back(v);
      }
    }
  }
  int n_diff = list_diff.size();
  MyMatrix<T> EXTdiff(n_diff, n + 1);
  for (int i = 0; i < n_diff; i++) {
    EXTdiff(i, 0) = 1;
    for (int k = 0; k < n; k++) {
      EXTdiff(i, k + 1) = list_diff[i](k);
    }
  }
  vectface vf_diff = DirectFacetComputationIncidence(EXTdiff, "cdd", os);
  data.LmatDiff = PolyNorm_LinearForms(EXTdiff, vf_diff);
#ifdef TIMINGS_POLYNORM_BASIC
  os << "|POLYNORM: facets of P - P|=" << time << "\n";
#endif
#ifdef DEBUG_POLYNORM_BASIC
  os << "POLYNORM: n=" << n << " n_vert=" << n_vert
     << " m=" << data.Lmat.rows() << " n_diff=" << n_diff
     << " m_diff=" << data.LmatDiff.rows() << "\n";
#endif
  // The symmetry group. The integral automorphism code takes an integral
  // vertex matrix; a uniform scaling of the homogeneous rows changes
  // neither the affine maps nor their integrality.
  MyMatrix<T> EXT_scal = RemoveFractionMatrix(EXT);
  MyMatrix<Tint> EXT_int = UniversalMatrixConversion<Tint, T>(EXT_scal);
#ifdef DEBUG_POLYNORM_BASIC
  os << "POLYNORM: EXT_int=\n";
  WriteMatrix(os, EXT_int);
#endif
  data.GRP = LinPolytopeIntegral_Automorphism<Tint, Tgroup>(EXT_int, os);
#ifdef TIMINGS_POLYNORM_BASIC
  os << "|POLYNORM: LinPolytopeIntegral_Automorphism|=" << time << "\n";
#endif
  std::unordered_map<Face, int> map_incd;
  int m = data.ListIncd.size();
  for (int j = 0; j < m; j++) {
    map_incd[data.ListIncd[j]] = j;
  }
  for (auto &eGen : data.GRP.GeneratorsOfGroup()) {
    data.ListGen.push_back(PolyNorm_SymmFromPerm<T, Tint, Telt>(
        EXT, data.ListIncd, map_incd, eGen));
  }
  if (data.GRP.size() <= PolyNorm_MaxGroupSizeElements) {
    std::vector<PolyNormSymm<Tint>> ListElt;
    for (auto &eElt : data.GRP.get_all_element()) {
      ListElt.push_back(PolyNorm_SymmFromPerm<T, Tint, Telt>(
          EXT, data.ListIncd, map_incd, eElt));
    }
    data.ListElt = std::move(ListElt);
  }
#ifdef DEBUG_POLYNORM_BASIC
  os << "POLYNORM: |GRP|=" << data.GRP.size()
     << " n_gen=" << data.ListGen.size()
     << " elements listed=" << data.ListElt.has_value() << "\n";
#endif
#ifdef TIMINGS_POLYNORM_BASIC
  os << "|POLYNORM: group elements as matrices|=" << time << "\n";
#endif
  return data;
}

// Representatives of the orbits of the facets under the symmetry group,
// as indices into the rows of Lmat.
template <typename T, typename Tint, typename Tgroup>
std::vector<int>
PolyNorm_FacetOrbitRepresentatives(PolyNormData<T, Tint, Tgroup> const &data) {
  std::unordered_map<Face, int> map_incd;
  int m = data.ListIncd.size();
  for (int j = 0; j < m; j++) {
    map_incd[data.ListIncd[j]] = j;
  }
  vectface vf_orb = OrbitSplittingSet(data.ListIncd, data.GRP);
  std::vector<int> ListRep;
  for (auto &f : vf_orb) {
    ListRep.push_back(map_incd.at(f));
  }
  std::sort(ListRep.begin(), ListRep.end());
  return ListRep;
}

// The orbit of each facet under the group generated by the facet
// permutations of a list of symmetries: the orbits are numbered in the
// order of their smallest facet and the entry j is the number of the orbit
// of the facet j.
template <typename Tint>
std::vector<int>
PolyNorm_FacetOrbitIds(std::vector<PolyNormSymm<Tint>> const &ListSymm,
                       int m) {
  std::vector<int> orbit_id(m, -1);
  int n_orbit = 0;
  for (int j = 0; j < m; j++) {
    if (orbit_id[j] == -1) {
      std::vector<int> orbit{j};
      orbit_id[j] = n_orbit;
      size_t pos = 0;
      while (pos < orbit.size()) {
        int u = orbit[pos];
        for (auto &eSymm : ListSymm) {
          int v = eSymm.facet_perm[u];
          if (orbit_id[v] == -1) {
            orbit_id[v] = n_orbit;
            orbit.push_back(v);
          }
        }
        pos++;
      }
      n_orbit++;
    }
  }
  return orbit_id;
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_BASIC_H_
// clang-format on
