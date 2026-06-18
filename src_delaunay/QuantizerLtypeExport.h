// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Export of a (simplex-Delaunay) L-type to JSON for downstream symbolic
// polynomial construction in Julia. See polyhedral_common/julia/ for the
// consumer.
//
// The exported JSON contains:
//   - dim                                  lattice dimension n
//   - tspace.basis                         T-space basis (symmetric matrices)
//   - cone_facets_tspace                   linear forms positive in the cone's
//                                          relative interior, in T-space coords
//   - delaunay_simplices                   star(0, 𝓓): for each Delaunay
//                                          n-simplex containing the origin,
//                                          its n non-zero lattice vertices
//   - triangulation                        boundary (n-1)-simplices of DV(0),
//                                          each as a list of n star indices.
//                                          Julia adds the apex 0 to get the
//                                          n-simplex {0, c_{i_1}, …, c_{i_n}}
//                                          on which to integrate ‖x‖²_Q.
//   - test_point_Q                         the Gram matrix from the
//                                          IsoDelaunayDomain (interior of cone)
//
// Assumes every Delaunay polytope in IDD.DT is an n-simplex.
//
// Algorithm:
//   1. BFS in the Delaunay adjacency graph to enumerate star(0, 𝓓).
//   2. Solve for the numerical circumcenter c_i of each star member D_i,
//      using IDD.GramMat as Q_0 (a representative quadratic form interior
//      to the L-type cone).
//   3. Triangulate the polytope DV(0) = conv(c_1, …, c_N) via lrs.
//   4. Extract its boundary (n-1)-simplices: each (n-1)-facet of a
//      triangulation n-simplex appears in either two n-simplices (interior)
//      or one (boundary). Keep the singletons.
//   5. Cone facets via the existing ComputeDefiningIneqIsoDelaunayDomain.
//   6. Emit JSON.
//
// clang-format off
#ifndef SRC_DELAUNAY_QUANTIZERLTYPEEXPORT_H_
#define SRC_DELAUNAY_QUANTIZERLTYPEEXPORT_H_

#include "IsoDelaunayDomains.h"
#include "POLY_lrslib.h"
#include <fstream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>
// clang-format on

namespace quantizer_export {

// ---------------------------------------------------------------------------
// Combinatorial helpers (star(0, 𝓓) enumeration)
// ---------------------------------------------------------------------------
//
// Per the project's convention, Delaunay-vertex matrices use the rational
// type `T` even though entries are integer-valued (homogeneous coord = 1,
// lattice coords = integers). We therefore template these helpers on `T`.

template <typename T>
MyMatrix<T> TranslateToOrigin(MyMatrix<T> const &EXT, int v_idx) {
  int nv = EXT.rows();
  int nc = EXT.cols();
  MyMatrix<T> out(nv, nc);
  for (int j = 0; j < nv; j++) {
    out(j, 0) = T(1);
    for (int k = 1; k < nc; k++) {
      out(j, k) = EXT(j, k) - EXT(v_idx, k);
    }
  }
  return out;
}

template <typename T>
bool ContainsOrigin(MyMatrix<T> const &EXT) {
  int nv = EXT.rows();
  int nc = EXT.cols();
  for (int j = 0; j < nv; j++) {
    bool is_zero = true;
    for (int k = 1; k < nc; k++) {
      if (EXT(j, k) != 0) {
        is_zero = false;
        break;
      }
    }
    if (is_zero) {
      return true;
    }
  }
  return false;
}

template <typename T>
std::vector<std::vector<T>> CanonicalKey(MyMatrix<T> const &EXT) {
  int nv = EXT.rows();
  int nc = EXT.cols();
  std::vector<std::vector<T>> rows;
  rows.reserve(nv);
  for (int j = 0; j < nv; j++) {
    std::vector<T> v(nc - 1);
    for (int k = 1; k < nc; k++) {
      v[k - 1] = EXT(j, k);
    }
    rows.push_back(std::move(v));
  }
  std::sort(rows.begin(), rows.end());
  return rows;
}

template <typename T, typename Tgroup>
std::vector<MyMatrix<T>>
EnumerateStarOf0(DelaunayTesselation<T, Tgroup> const &DT, std::ostream &os) {
  std::set<std::vector<std::vector<T>>> seen;
  std::vector<MyMatrix<T>> star;
  struct State {
    int i_orb;
    MyMatrix<T> EXT;
  };
  std::queue<State> bfs;

  int n_orbits = DT.l_dels.size();
  for (int i_orb = 0; i_orb < n_orbits; i_orb++) {
    MyMatrix<T> const &EXT_orb = DT.l_dels[i_orb].EXT;
    int nv = EXT_orb.rows();
    for (int v = 0; v < nv; v++) {
      MyMatrix<T> EXT_t = TranslateToOrigin(EXT_orb, v);
      auto key = CanonicalKey(EXT_t);
      if (seen.insert(key).second) {
        bfs.push({i_orb, EXT_t});
        star.push_back(EXT_t);
      }
    }
  }

  while (!bfs.empty()) {
    State s = bfs.front();
    bfs.pop();
    Delaunay_Entry<T, Tgroup> const &entry = DT.l_dels[s.i_orb];
    int nc = s.EXT.cols();
    MyVector<T> t_shift(nc - 1);
    for (int k = 1; k < nc; k++) {
      t_shift(k - 1) = s.EXT(0, k) - entry.EXT(0, k);
    }
    for (Delaunay_AdjO<T> const &adj : entry.ListAdj) {
      MyMatrix<T> nbr_in_orbframe = DT.l_dels[adj.iOrb].EXT * adj.eBigMat;
      MyMatrix<T> nbr_EXT(nbr_in_orbframe.rows(), nc);
      for (int j = 0; j < nbr_in_orbframe.rows(); j++) {
        nbr_EXT(j, 0) = T(1);
        for (int k = 1; k < nc; k++) {
          nbr_EXT(j, k) = nbr_in_orbframe(j, k) + t_shift(k - 1);
        }
      }
      if (!ContainsOrigin(nbr_EXT)) {
        continue;
      }
      auto key = CanonicalKey(nbr_EXT);
      if (seen.insert(key).second) {
        bfs.push({adj.iOrb, nbr_EXT});
        star.push_back(nbr_EXT);
      }
    }
  }

  os << "QuantExport: enumerated |star(0,D)| = " << star.size() << "\n";
  return star;
}

// ---------------------------------------------------------------------------
// Numerical circumcenter under a sample Gram matrix Q_0
// ---------------------------------------------------------------------------

// Solve for the circumcenter c of the simplex D (with origin among its
// vertices) in the metric Q_0.
//
// Setup: v_0 = 0, v_1, …, v_n are the vertices of D. The circumcenter satisfies
//   |c - v_k|²_{Q_0} = |c|²_{Q_0}   ⇔   2 v_k^T Q_0 c = |v_k|²_{Q_0}.
// Let M have rows v_k^T (k=1..n) and b_k = ½ v_k^T Q_0 v_k. Then
//   c = Q_0^{-1} · M^{-1} · b.
template <typename T>
MyVector<T> CircumcenterUnder(MyMatrix<T> const &D, MyMatrix<T> const &Q0) {
  int nv = D.rows();
  int n = D.cols() - 1;
  int zero_row = -1;
  for (int j = 0; j < nv; j++) {
    bool is_zero = true;
    for (int k = 1; k <= n; k++) {
      if (D(j, k) != 0) {
        is_zero = false;
        break;
      }
    }
    if (is_zero) {
      zero_row = j;
      break;
    }
  }
  if (zero_row < 0) {
    std::cerr << "CircumcenterUnder: simplex does not contain origin\n";
    throw TerminalException{1};
  }

  MyMatrix<T> M(n, n);
  MyVector<T> b(n);
  int row = 0;
  T half(1);
  half /= T(2);
  for (int j = 0; j < nv; j++) {
    if (j == zero_row) {
      continue;
    }
    MyVector<T> v(n);
    for (int k = 0; k < n; k++) {
      v(k) = D(j, k + 1);
    }
    for (int k = 0; k < n; k++) {
      M(row, k) = v(k);
    }
    T q(0);
    for (int a = 0; a < n; a++) {
      for (int c = 0; c < n; c++) {
        q += v(a) * Q0(a, c) * v(c);
      }
    }
    b(row) = half * q;
    row++;
  }

  MyMatrix<T> Minv = Inverse(M);
  MyMatrix<T> Q0inv = Inverse(Q0);
  MyVector<T> Qc = Minv * b;   // = Q_0 c
  MyVector<T> c = Q0inv * Qc;
  return c;
}

// ---------------------------------------------------------------------------
// Boundary extraction from the lrs triangulation of DV(0)
// ---------------------------------------------------------------------------

// Given a triangulation of an n-polytope into n-simplices (each represented
// as a Face = boolean vector over N vertex indices), return the list of
// boundary (n-1)-simplices: those (n-1)-faces appearing in exactly one
// n-simplex.
inline std::vector<std::vector<int>>
ExtractBoundarySimplices(vectface const &trig, std::ostream &os) {
  std::map<std::vector<int>, int> facet_count;
  for (size_t t = 0; t < trig.size(); t++) {
    Face const &T_simp = trig[t];
    std::vector<int> verts;
    boost::dynamic_bitset<>::size_type pos = T_simp.find_first();
    while (pos != boost::dynamic_bitset<>::npos) {
      verts.push_back(static_cast<int>(pos));
      pos = T_simp.find_next(pos);
    }
    // Each (n-1)-facet of this n-simplex is obtained by dropping one vertex.
    int nverts = static_cast<int>(verts.size());
    for (int skip = 0; skip < nverts; skip++) {
      std::vector<int> facet;
      facet.reserve(nverts - 1);
      for (int i = 0; i < nverts; i++) {
        if (i != skip) {
          facet.push_back(verts[i]);
        }
      }
      facet_count[facet] += 1;
    }
  }
  std::vector<std::vector<int>> boundary;
  for (auto const &kv : facet_count) {
    if (kv.second == 1) {
      boundary.push_back(kv.first);
    } else if (kv.second != 2) {
      std::cerr << "QuantExport: ill-formed triangulation, facet shared by "
                << kv.second << " simplices (expected 1 or 2)\n";
      throw TerminalException{1};
    }
  }
  os << "QuantExport: boundary (n-1)-simplices = " << boundary.size() << "\n";
  return boundary;
}

// ---------------------------------------------------------------------------
// JSON emission
// ---------------------------------------------------------------------------

template <typename T>
void WriteRationalJSON(std::ostream &out, T const &val) {
  // For now, the matrices we emit have integer or simple-rational entries.
  // Boost / mpq_class's operator<< produces a fractional form like "1/2"
  // which is not valid JSON. We sidestep by producing a string with explicit
  // numerator / denominator handling only when needed.
  std::ostringstream tmp;
  tmp << val;
  std::string s = tmp.str();
  // If the string contains '/', convert to JSON array [num, den] to keep
  // strict JSON conformance. Otherwise emit as plain number.
  auto pos = s.find('/');
  if (pos == std::string::npos) {
    out << s;
  } else {
    std::string num = s.substr(0, pos);
    std::string den = s.substr(pos + 1);
    out << "[" << num << ", " << den << "]";
  }
}

template <typename T, typename Tint, typename Tgroup>
void WriteQuantizerLtypeJSON(IsoDelaunayDomain<T, Tint, Tgroup> const &IDD,
                             LinSpaceMatrix<T> const &LinSpa,
                             std::string const &output_path,
                             std::ostream &os) {
  int n = LinSpa.n;

  std::vector<FullAdjInfo<T>> facets =
      ComputeDefiningIneqIsoDelaunayDomain<T, Tgroup>(IDD.DT, LinSpa.ListLineMat,
                                                     os);

  std::vector<MyMatrix<T>> star = EnumerateStarOf0<T, Tgroup>(IDD.DT, os);
  int N = static_cast<int>(star.size());

  MyMatrix<T> EXT_DV(N, n + 1);
  for (int i = 0; i < N; i++) {
    MyVector<T> c = CircumcenterUnder<T>(star[i], IDD.GramMat);
    EXT_DV(i, 0) = T(1);
    for (int k = 0; k < n; k++) {
      EXT_DV(i, k + 1) = c(k);
    }
  }

  vectface trig = lrs::GetTriangulation(EXT_DV);
  os << "QuantExport: DV(0) triangulated into " << trig.size() << " n-simplices\n";

  std::vector<std::vector<int>> boundary = ExtractBoundarySimplices(trig, os);

  std::ofstream out(output_path);
  out << "{\n";
  out << "  \"schema_version\": 1,\n";
  out << "  \"dim\": " << n << ",\n";

  out << "  \"tspace\": {\n";
  out << "    \"basis\": [\n";
  for (size_t i = 0; i < LinSpa.ListMat.size(); i++) {
    MyMatrix<T> const &M = LinSpa.ListMat[i];
    out << "      [";
    for (int r = 0; r < M.rows(); r++) {
      if (r > 0) out << ", ";
      out << "[";
      for (int c = 0; c < M.cols(); c++) {
        if (c > 0) out << ", ";
        WriteRationalJSON(out, M(r, c));
      }
      out << "]";
    }
    out << "]";
    if (i + 1 < LinSpa.ListMat.size()) out << ",";
    out << "\n";
  }
  out << "    ]\n";
  out << "  },\n";

  out << "  \"cone_facets_tspace\": [\n";
  for (size_t i = 0; i < facets.size(); i++) {
    MyVector<T> const &v = facets[i].eIneq;
    out << "    [";
    for (int k = 0; k < v.size(); k++) {
      if (k > 0) out << ", ";
      WriteRationalJSON(out, v(k));
    }
    out << "]";
    if (i + 1 < facets.size()) out << ",";
    out << "\n";
  }
  out << "  ],\n";

  out << "  \"delaunay_simplices\": [\n";
  for (int i = 0; i < N; i++) {
    MyMatrix<T> const &D = star[i];
    out << "    {\"vertices\": [";
    bool first = true;
    for (int j = 0; j < D.rows(); j++) {
      bool is_zero = true;
      for (int k = 1; k < D.cols(); k++) {
        if (D(j, k) != 0) {
          is_zero = false;
          break;
        }
      }
      if (is_zero) continue;
      if (!first) out << ", ";
      first = false;
      out << "[";
      for (int k = 1; k < D.cols(); k++) {
        if (k > 1) out << ", ";
        // Lattice coords are integer-valued; convert through Tint for clean
        // JSON output (otherwise mpq_class prints as "n/1").
        Tint val = UniversalScalarConversion<Tint, T>(D(j, k));
        out << val;
      }
      out << "]";
    }
    out << "]}";
    if (i + 1 < N) out << ",";
    out << "\n";
  }
  out << "  ],\n";

  out << "  \"triangulation\": [\n";
  for (size_t i = 0; i < boundary.size(); i++) {
    out << "    {\"delaunay_indices\": [";
    for (size_t k = 0; k < boundary[i].size(); k++) {
      if (k > 0) out << ", ";
      out << boundary[i][k];
    }
    out << "]}";
    if (i + 1 < boundary.size()) out << ",";
    out << "\n";
  }
  out << "  ],\n";

  out << "  \"test_point_Q\": [";
  for (int r = 0; r < IDD.GramMat.rows(); r++) {
    if (r > 0) out << ", ";
    out << "[";
    for (int c = 0; c < IDD.GramMat.cols(); c++) {
      if (c > 0) out << ", ";
      WriteRationalJSON(out, IDD.GramMat(r, c));
    }
    out << "]";
  }
  out << "]\n";
  out << "}\n";
  out.close();

  os << "QuantExport: wrote " << output_path << "\n";
}

}  // namespace quantizer_export

#endif  // SRC_DELAUNAY_QUANTIZERLTYPEEXPORT_H_
