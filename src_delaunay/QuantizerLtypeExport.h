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

// Canonicalize a Delaunay-vertex matrix modulo Z^n translation: translate so
// the lex-smallest non-homogeneous row becomes the origin, then sort the rows.
template <typename T>
std::vector<std::vector<T>> CanonicalKeyModZn(MyMatrix<T> const &EXT) {
  int nv = EXT.rows();
  int nc = EXT.cols();
  // Find lex-min row by lattice coords (columns 1..).
  int min_row = 0;
  for (int j = 1; j < nv; j++) {
    bool smaller = false;
    for (int k = 1; k < nc; k++) {
      if (EXT(j, k) < EXT(min_row, k)) {
        smaller = true;
        break;
      }
      if (EXT(j, k) > EXT(min_row, k)) {
        break;
      }
    }
    if (smaller) {
      min_row = j;
    }
  }
  std::vector<std::vector<T>> rows;
  rows.reserve(nv);
  for (int j = 0; j < nv; j++) {
    std::vector<T> v(nc - 1);
    for (int k = 1; k < nc; k++) {
      v[k - 1] = EXT(j, k) - EXT(min_row, k);
    }
    rows.push_back(std::move(v));
  }
  std::sort(rows.begin(), rows.end());
  return rows;
}

// Lift an n×n matrix g (acting on lattice rows by right multiplication) to an
// (n+1)×(n+1) matrix that fixes the homogeneous column. Result is a block
// matrix [[1, 0], [0, g]] when rows are (1, v).
template <typename T>
MyMatrix<T> LiftToHomog(MyMatrix<T> const &g) {
  int n = g.rows();
  MyMatrix<T> out(n + 1, n + 1);
  for (int i = 0; i < n + 1; i++) {
    for (int j = 0; j < n + 1; j++) {
      out(i, j) = T(0);
    }
  }
  out(0, 0) = T(1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      out(i + 1, j + 1) = g(i, j);
    }
  }
  return out;
}

// Enumerate Delaunay simplices containing the origin. Implementation follows
// the orbit-of-shapes approach:
//   1. Compute the lattice automorphism group Aut(L) via LINSPA_ComputeStabilizer.
//   2. Lift each n×n generator to a (n+1)×(n+1) homogeneous matrix.
//   3. For each Delaunay orbit rep, BFS the Aut(L)-orbit of its vertex matrix,
//      canonicalizing modulo Z^n translation, to collect all distinct Delaunay
//      shapes that occur in star(0, 𝓓).
//   4. For each shape, translate by -v for each vertex v to put it at origin;
//      the resulting matrix is one member of star(0, 𝓓).
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>>
EnumerateStarOf0(DelaunayTesselation<T, Tgroup> const &DT,
                 LinSpaceMatrix<T> const &LinSpa,
                 MyMatrix<T> const &Q0, std::ostream &os) {
  std::vector<MyMatrix<T>> autL_gens =
      LINSPA_ComputeStabilizer<T, Tint, Tgroup>(LinSpa, Q0, os);
  os << "QuantExport: Aut(L) has " << autL_gens.size() << " generator(s)\n";

  std::vector<MyMatrix<T>> autL_homog;
  autL_homog.reserve(autL_gens.size());
  for (MyMatrix<T> const &g : autL_gens) {
    autL_homog.push_back(LiftToHomog(g));
  }

  // BFS the orbit of each Delaunay orbit rep under Aut(L), modulo Z^n.
  std::set<std::vector<std::vector<T>>> seen_shapes;
  std::vector<MyMatrix<T>> shape_reps;
  for (Delaunay_Entry<T, Tgroup> const &entry : DT.l_dels) {
    auto k0 = CanonicalKeyModZn(entry.EXT);
    if (seen_shapes.insert(k0).second) {
      shape_reps.push_back(entry.EXT);
    }
  }
  {
    std::queue<MyMatrix<T>> q;
    for (MyMatrix<T> const &S : shape_reps) {
      q.push(S);
    }
    while (!q.empty()) {
      MyMatrix<T> X = q.front();
      q.pop();
      for (MyMatrix<T> const &g_h : autL_homog) {
        MyMatrix<T> gX = X * g_h;
        auto key = CanonicalKeyModZn(gX);
        if (seen_shapes.insert(key).second) {
          q.push(gX);
          shape_reps.push_back(gX);
        }
      }
    }
  }
  os << "QuantExport: |distinct shapes mod Z^n| = " << shape_reps.size() << "\n";

  // For each shape, translate so each vertex becomes the origin.
  std::set<std::vector<std::vector<T>>> seen_star;
  std::vector<MyMatrix<T>> star;
  for (MyMatrix<T> const &D : shape_reps) {
    int nv = D.rows();
    for (int v = 0; v < nv; v++) {
      MyMatrix<T> D_t = TranslateToOrigin(D, v);
      auto key = CanonicalKey(D_t);
      if (seen_star.insert(key).second) {
        star.push_back(D_t);
      }
    }
  }
  os << "QuantExport: enumerated |star(0,D)| = " << star.size() << "\n";
  return star;
}

// ---------------------------------------------------------------------------
// Cone-facet inequality orbit expansion under Aut(L)
// ---------------------------------------------------------------------------
//
// The classic-T-space basis (TSPACE_canonical_get_list_matrices) is:
//   M_{aa} = E_{aa}                       (diagonal entries)
//   M_{ab} = E_{ab} + E_{ba}              (off-diagonal entries, a < b)
//
// A linear inequality V · a > 0 on T-space coords corresponds to the symmetric
// "dual" matrix X (with V[idx(a,a)] = X[a,a] and V[idx(a,b)] = 2·X[a,b] for
// a < b). Aut(L) acts on Q ∈ S^n by Q → g^T Q g; the dual action on X is
// X → g X g^T.

template <typename T>
MyVector<T> InequalityFromSymMatClassic(MyMatrix<T> const &X) {
  int n = X.rows();
  int d = n * (n + 1) / 2;
  MyVector<T> V(d);
  int idx = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      T mult = (i == j) ? T(1) : T(2);
      V(idx) = mult * X(i, j);
      idx++;
    }
  }
  return V;
}

template <typename T>
MyMatrix<T> SymMatFromInequalityClassic(MyVector<T> const &V, int n) {
  MyMatrix<T> X(n, n);
  int idx = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      T mult = (i == j) ? T(1) : T(2);
      X(i, j) = V(idx) / mult;
      X(j, i) = X(i, j);
      idx++;
    }
  }
  return X;
}

template <typename T>
std::vector<MyVector<T>> OrbitInequalityUnderAutL(
    MyVector<T> const &V_start,
    std::vector<MyMatrix<T>> const &autL_gens,
    int n) {
  auto to_key = [](MyVector<T> const &v) -> std::vector<T> {
    std::vector<T> out(v.size());
    for (int i = 0; i < v.size(); i++) {
      out[i] = v(i);
    }
    return out;
  };
  std::set<std::vector<T>> seen;
  std::vector<MyVector<T>> result;
  std::queue<MyVector<T>> q;
  MyVector<T> V0 = ScalarCanonicalizationVector(V_start);
  if (seen.insert(to_key(V0)).second) {
    q.push(V0);
    result.push_back(V0);
  }
  while (!q.empty()) {
    MyVector<T> V = q.front();
    q.pop();
    MyMatrix<T> X = SymMatFromInequalityClassic(V, n);
    for (MyMatrix<T> const &g : autL_gens) {
      MyMatrix<T> X_new = g * X * g.transpose();
      MyVector<T> V_new = InequalityFromSymMatClassic(X_new);
      V_new = ScalarCanonicalizationVector(V_new);
      if (seen.insert(to_key(V_new)).second) {
        q.push(V_new);
        result.push_back(V_new);
      }
    }
  }
  return result;
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

  std::vector<FullAdjInfo<T>> facets_raw =
      ComputeDefiningIneqIsoDelaunayDomain<T, Tgroup>(IDD.DT, LinSpa.ListLineMat,
                                                     os);
  os << "QuantExport: ComputeDefiningIneq returned " << facets_raw.size()
     << " orbit rep(s) of inequalities\n";

  // Orbit-expand inequalities under Aut(L) using matrix_in_t_space (which
  // accounts for the two transposes — Q-space action vs functional action,
  // row vs column action). Then dedup and filter redundant.
  std::vector<MyMatrix<T>> autL_gens_for_ineq =
      LINSPA_ComputeStabilizer<T, Tint, Tgroup>(LinSpa, IDD.GramMat, os);
  os << "QuantExport: Aut(L) has " << autL_gens_for_ineq.size()
     << " generator(s) for inequality expansion\n";
  std::vector<MyMatrix<T>> tspace_action;
  tspace_action.reserve(autL_gens_for_ineq.size());
  for (MyMatrix<T> const &g : autL_gens_for_ineq) {
    tspace_action.push_back(matrix_in_t_space(g, LinSpa));
  }
  auto to_key = [](MyVector<T> const &v) -> std::vector<T> {
    std::vector<T> out(v.size());
    for (int i = 0; i < v.size(); i++) {
      out[i] = v(i);
    }
    return out;
  };
  std::set<std::vector<T>> facet_seen;
  std::vector<MyVector<T>> facets_all;
  std::queue<MyVector<T>> ineq_bfs;
  for (FullAdjInfo<T> const &fai : facets_raw) {
    MyVector<T> v0 = RemoveFractionVector(fai.eIneq);
    if (facet_seen.insert(to_key(v0)).second) {
      ineq_bfs.push(v0);
      facets_all.push_back(v0);
    }
  }
  while (!ineq_bfs.empty()) {
    MyVector<T> v = ineq_bfs.front();
    ineq_bfs.pop();
    for (MyMatrix<T> const &MatSpace : tspace_action) {
      MyVector<T> v_img = MatSpace * v;
      MyVector<T> v_img_red = RemoveFractionVector(v_img);
      if (facet_seen.insert(to_key(v_img_red)).second) {
        ineq_bfs.push(v_img_red);
        facets_all.push_back(v_img_red);
      }
    }
  }
  os << "QuantExport: after Aut(L) orbit-expansion, |all inequalities| = "
     << facets_all.size() << "\n";

  // Filter to non-redundant (= actual cone facets) using lrs-based redundancy
  // detection, matching f_adj's pattern.
  MyMatrix<T> FAC(facets_all.size(), facets_all[0].size());
  for (size_t i = 0; i < facets_all.size(); i++) {
    for (int j = 0; j < facets_all[i].size(); j++) {
      FAC(i, j) = facets_all[i](j);
    }
  }
  MyMatrix<T> FAC_extend = AddFirstZeroColumn(FAC);
  std::vector<int> ListIrred = get_non_redundant_indices(FAC_extend, os);
  std::vector<MyVector<T>> facets;
  for (int idx : ListIrred) {
    facets.push_back(facets_all[idx]);
  }
  os << "QuantExport: irredundant cone facets = " << facets.size() << "\n";

  std::vector<MyMatrix<T>> star =
      EnumerateStarOf0<T, Tint, Tgroup>(IDD.DT, LinSpa, IDD.GramMat, os);
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
    MyVector<T> const &v = facets[i];
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

  out << "  \"autL_generators\": [\n";
  for (size_t i = 0; i < autL_gens_for_ineq.size(); i++) {
    MyMatrix<T> const &g = autL_gens_for_ineq[i];
    out << "    [";
    for (int r = 0; r < g.rows(); r++) {
      if (r > 0) out << ", ";
      out << "[";
      for (int c = 0; c < g.cols(); c++) {
        if (c > 0) out << ", ";
        WriteRationalJSON(out, g(r, c));
      }
      out << "]";
    }
    out << "]";
    if (i + 1 < autL_gens_for_ineq.size()) out << ",";
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
