// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_IGUSA_IGUSA_FACE_VERTICES_H_
#define SRC_IGUSA_IGUSA_FACE_VERTICES_H_

// clang-format off
#include "igusa_facet_enum.h"
#include <functional>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_IGUSA_FACE_VERTICES
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_IGUSA_FACE_VERTICES
#endif

#ifdef TIMINGS
#define TIMINGS_IGUSA_FACE_VERTICES
#endif

/*
  The vertices of the face of P = conv(I_n) of a valid inequality
  tr(F X) >= b, from the list of the orbits of vertices of P.

  By the structure theorem of the faces, F = B^T F0 B with F0 positive
  definite of rank r and the face is the lift of the bounded face
  Phi0 = conv{ Y in I_r : tr(F0 Y) = b } of P_r. The vertices of Phi0 are
  vertices of P_r, so they are the forms X = g A g^T with A in the list of
  the vertex representatives of P_r, g in GL_r(Z) and tr(F0 X) = b.

  The theory of Voronoi gives F0 = sum_v lambda_v v v^T, lambda_v > 0, over
  minimal vectors v of the minimizer of tr(F0 X) over the Ryshkov
  polyhedron, with sum_v lambda_v = mu. Then
      tr(F0 X) = sum_v lambda_v A[x_v],   x_v = g^T v,
  and A[x_v] is an integer >= 1, so with the slack s = b - mu:
      sum_v lambda_v (A[x_v] - 1) = s.
  The v span R^r. For r independent support vectors v_1, ..., v_r the
  images x_{v_i} determine g, and every other support vector u = sum c_i v_i
  has x_u = sum c_i x_{v_i}. The x_{v_i} are searched by backtracking among
  the vectors of A of norm at most 1 + s / lambda_{v_i}; each support vector
  determined by the assigned ones is charged lambda_u (A[x_u] - 1) to the
  slack. Since g and g k, k in Aut(A), give the same X, each x_{v_i} is
  taken up to the stabilizer in Aut(A) of the x_{v_j}, j < i.
 */

template <typename T> struct IgusaFaceVertices {
  // The rank r of F, a basis B of W (r x n) and F0 with F = B^T F0 B
  int rank;
  MyMatrix<T> B;
  MyMatrix<T> F0;
  T ryshkov_min;
  T slack;
  int n_support;
  // The vertices of Phi0 (r x r) and the index of their representative
  std::vector<MyMatrix<T>> l_vertex;
  std::vector<size_t> l_type;
  // The dimension of Phi0, the affine dimension of its vertices (-1 if none)
  int dim_face;
  // The number of nodes of the backtracking
  size_t n_node;
};

// The affine dimension of a list of symmetric matrices, -1 if it is empty.
template <typename T>
int igusa_affine_dimension(std::vector<MyMatrix<T>> const &l_mat) {
  if (l_mat.empty()) {
    return -1;
  }
  int n = l_mat[0].rows();
  int len = n * (n + 1) / 2;
  MyMatrix<T> M(l_mat.size(), len + 1);
  for (size_t i = 0; i < l_mat.size(); i++) {
    MyVector<T> V = SymmetricMatrixToVector(l_mat[i]);
    M(i, 0) = T(1);
    for (int j = 0; j < len; j++) {
      M(i, j + 1) = V(j);
    }
  }
  return RankMat(M) - 1;
}

// The generators h of Aut(A) as maps x -> M x on the vectors, A[M x] = A[x].
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> igusa_vector_automorphisms(MyMatrix<T> const &A,
                                                       std::ostream &os) {
  std::vector<MyMatrix<Tint>> l_h =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(A, os);
  std::vector<MyMatrix<Tint>> l_M;
  for (auto &h : l_h) {
    MyMatrix<T> h_T = UniversalMatrixConversion<T, Tint>(h);
    if (h_T * A * h_T.transpose() == A) {
      l_M.push_back(h.transpose());
      continue;
    }
    if (h_T.transpose() * A * h_T == A) {
      l_M.push_back(h);
      continue;
    }
    std::cerr << "IGUSA_FACE_VERTICES: not an automorphism of A\n";
    throw TerminalException{1};
  }
  return l_M;
}

/*
  The order of the basis vectors among the support: each new vector is the
  one that determines the most support vectors together with the previous
  ones, ties broken by the largest weight (fewest candidates).
 */
template <typename T, typename Tint>
std::vector<size_t> igusa_support_basis_order(
    std::vector<std::pair<MyVector<Tint>, T>> const &l_weight, int r) {
  int n_w = l_weight.size();
  std::vector<size_t> l_sel;
  auto rank_with = [&](std::vector<size_t> const &l_idx) -> int {
    MyMatrix<T> M(l_idx.size(), r);
    for (size_t i = 0; i < l_idx.size(); i++) {
      for (int j = 0; j < r; j++) {
        M(i, j) = UniversalScalarConversion<T, Tint>(l_weight[l_idx[i]].first(j));
      }
    }
    return RankMat(M);
  };
  auto n_in_span = [&](std::vector<size_t> const &l_idx) -> int {
    int rk = l_idx.size();
    int cnt = 0;
    for (int u = 0; u < n_w; u++) {
      std::vector<size_t> l_ext = l_idx;
      l_ext.push_back(u);
      if (rank_with(l_ext) == rk) {
        cnt++;
      }
    }
    return cnt;
  };
  while (static_cast<int>(l_sel.size()) < r) {
    int best = -1;
    int best_cnt = -1;
    for (int u = 0; u < n_w; u++) {
      std::vector<size_t> l_ext = l_sel;
      l_ext.push_back(u);
      if (rank_with(l_ext) != static_cast<int>(l_ext.size())) {
        continue;
      }
      int cnt = n_in_span(l_ext);
      if (cnt > best_cnt ||
          (cnt == best_cnt && l_weight[u].second > l_weight[best].second)) {
        best = u;
        best_cnt = cnt;
      }
    }
    if (best == -1) {
      std::cerr << "IGUSA_FACE_VERTICES: the support does not span R^r\n";
      throw TerminalException{1};
    }
    l_sel.push_back(best);
  }
  return l_sel;
}

/*
  The vertices of the bounded face { Y in conv(I_r) : tr(F0 Y) = b } that
  are in the orbits of the forms of ListRepr, for F0 positive definite and
  tr(F0 Y) >= b valid. space is the full space of the symmetric r x r
  matrices.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<std::pair<MyMatrix<T>, size_t>>
igusa_bounded_face_vertices(IgusaSpace<T, Tint> &space, MyMatrix<T> const &F0,
                            T const &b,
                            std::vector<MyMatrix<T>> const &ListRepr,
                            IgusaFaceVertices<T> &info, std::ostream &os) {
  int r = F0.rows();
  IgusaVoronoiData<T, Tint> vd = igusa_voronoi_data(space, F0, os);
  std::vector<std::pair<MyVector<Tint>, T>> const &l_weight = vd.l_weight;
  int n_w = l_weight.size();
  T slack = b - vd.ryshkov_min;
  info.ryshkov_min = vd.ryshkov_min;
  info.slack = slack;
  info.n_support = n_w;
  info.n_node = 0;
  std::vector<std::pair<MyMatrix<T>, size_t>> l_ret;
  if (slack < 0) {
    return l_ret;
  }
  std::vector<size_t> l_basis = igusa_support_basis_order<T, Tint>(l_weight, r);
  // V: the basis vectors as columns, and the coefficients of the support
  // vectors in that basis, with the level where they get determined.
  MyMatrix<T> V(r, r);
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < r; j++) {
      V(j, i) = UniversalScalarConversion<T, Tint>(l_weight[l_basis[i]].first(j));
    }
  }
  MyMatrix<T> Vinv = Inverse(V);
  std::vector<MyVector<T>> l_coeff(n_w);
  std::vector<std::vector<int>> l_level_supp(r);
  for (int u = 0; u < n_w; u++) {
    MyVector<T> v_T = UniversalVectorConversion<T, Tint>(l_weight[u].first);
    l_coeff[u] = Vinv * v_T;
    int level = 0;
    for (int i = 0; i < r; i++) {
      if (l_coeff[u](i) != 0) {
        level = i;
      }
    }
    l_level_supp[level].push_back(u);
  }
#ifdef DEBUG_IGUSA_FACE_VERTICES
  os << "IGUSA_FACE_VERTICES: r=" << r << " |support|=" << n_w
     << " slack=" << slack << " supports by level:";
  for (int i = 0; i < r; i++) {
    os << " " << l_level_supp[i].size();
  }
  os << "\n";
#endif
  T lambda_min = l_weight[l_basis[0]].second;
  for (auto &idx : l_basis) {
    if (l_weight[idx].second < lambda_min) {
      lambda_min = l_weight[idx].second;
    }
  }
  T max_norm = T(1) + slack / lambda_min;
  std::set<MyMatrix<T>> set_found;
  for (size_t i_repr = 0; i_repr < ListRepr.size(); i_repr++) {
    MyMatrix<T> const &A = ListRepr[i_repr];
    // The vectors of A of norm at most max_norm, both signs, no zero
    std::set<MyVector<Tint>> set_x;
    for (auto &x : computeLevel_GramMat<T, Tint>(A, max_norm, os)) {
      if (IsZeroVector(x)) {
        continue;
      }
      set_x.insert(x);
      set_x.insert(MyVector<Tint>(-x));
    }
    std::vector<MyVector<Tint>> l_x(set_x.begin(), set_x.end());
    std::vector<T> l_norm;
    for (auto &x : l_x) {
      l_norm.push_back(EvaluationQuadForm<T, Tint>(A, x));
    }
    // Aut(A) as a permutation group on l_x, which it preserves
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    std::vector<MyMatrix<Tint>> l_M =
        igusa_vector_automorphisms<T, Tint, Tgroup>(A, os);
    std::map<MyVector<Tint>, Tidx> map_x;
    for (size_t i = 0; i < l_x.size(); i++) {
      map_x[l_x[i]] = i;
    }
    std::vector<Telt> l_perm;
    for (auto &M : l_M) {
      std::vector<Tidx> eList(l_x.size());
      for (size_t i = 0; i < l_x.size(); i++) {
        eList[i] = map_x.at(MyVector<Tint>(M * l_x[i]));
      }
      l_perm.push_back(Telt(eList));
    }
    // The stabilizer in Aut(A) of the vectors assigned before level k
    std::vector<Tgroup> l_grp(r + 1);
    l_grp[0] = Tgroup(l_perm, l_x.size());
    int n_x = l_x.size();
    // The backtracking: X_assigned holds the images x_{v_i} of the basis
    std::vector<MyVector<T>> l_assigned(r);
    std::function<void(int, T const &)> f_rec = [&](int k,
                                                    T const &rem) -> void {
      info.n_node++;
      if (k == r) {
        if (rem != 0) {
          return;
        }
        // g^T = Xm V^{-1}, Xm with the columns x_{v_i}
        MyMatrix<T> Xm(r, r);
        for (int i = 0; i < r; i++) {
          for (int j = 0; j < r; j++) {
            Xm(j, i) = l_assigned[i](j);
          }
        }
        MyMatrix<T> gT = Xm * Vinv;
        if (!IsIntegralMatrix(gT) || T_abs(DeterminantMat(gT)) != 1) {
          return;
        }
        MyMatrix<T> g = gT.transpose();
        MyMatrix<T> X = g * A * gT;
#ifdef SANITY_CHECK_IGUSA_FACE_VERTICES
        if (frobenius_inner(F0, X) != b) {
          std::cerr << "IGUSA_FACE_VERTICES: the form is not incident\n";
          throw TerminalException{1};
        }
#endif
        if (set_found.count(X) == 0) {
          set_found.insert(X);
          l_ret.push_back({X, i_repr});
        }
        return;
      }
      T const &lambda_k = l_weight[l_basis[k]].second;
      T bound_k = T(1) + rem / lambda_k;
      auto try_vector = [&](size_t i_x) -> void {
        if (l_norm[i_x] > bound_k) {
          return;
        }
        l_assigned[k] = UniversalVectorConversion<T, Tint>(l_x[i_x]);
        T rem2 = rem;
        for (auto &u : l_level_supp[k]) {
          MyVector<T> x_u = ZeroVector<T>(r);
          for (int i = 0; i <= k; i++) {
            if (l_coeff[u](i) != 0) {
              x_u += l_coeff[u](i) * l_assigned[i];
            }
          }
          if (!IsIntegralVector(x_u)) {
            return;
          }
          T norm = EvaluationQuadForm<T, T>(A, x_u);
          if (norm < 1) {
            return;
          }
          rem2 -= l_weight[u].second * (norm - 1);
          if (rem2 < 0) {
            return;
          }
        }
        f_rec(k + 1, rem2);
      };
      /*
        g and g k, k in Aut(A), give the same form and k acts on all the
        x_v at once, so x_{v_k} is taken up to the stabilizer of the
        vectors assigned before.
       */
      Tgroup const &grp = l_grp[k];
      if (grp.size() == 1) {
        l_grp[k + 1] = grp;
        for (int i_x = 0; i_x < n_x; i_x++) {
          try_vector(i_x);
        }
        return;
      }
      std::vector<int> orbit = igusa_orbit_index(grp, n_x);
      int i_orb_next = 0;
      for (int i_x = 0; i_x < n_x; i_x++) {
        if (orbit[i_x] != i_orb_next) {
          continue;
        }
        i_orb_next++;
        if (l_norm[i_x] > bound_k) {
          continue;
        }
        l_grp[k + 1] = grp.Stabilizer_OnPoints(static_cast<Tidx>(i_x));
        try_vector(i_x);
      }
    };
#ifdef TIMINGS_IGUSA_FACE_VERTICES
    MicrosecondTime time;
#endif
    f_rec(0, slack);
#ifdef TIMINGS_IGUSA_FACE_VERTICES
    os << "|IGUSA_FACE_VERTICES: representative " << i_repr
       << " |vectors|=" << l_x.size() << " |Aut(A)|=" << l_grp[0].size()
       << " n_node=" << info.n_node << "|=" << time
       << "\n";
#endif
  }
  return l_ret;
}

/*
  The vertices of the face of tr(F X) >= b, F positive semidefinite, as
  vertices of the bounded face Phi0 of P_r. ListRepr are the vertex
  representatives of P_r, r the rank of F.
 */
template <typename T, typename Tint, typename Tgroup>
IgusaFaceVertices<T>
igusa_face_vertices(IgusaParameters<T> const &params, MyMatrix<T> const &F,
                    T const &b, std::vector<MyMatrix<T>> const &ListRepr,
                    std::ostream &os) {
  IgusaFaceVertices<T> info;
  IgusaImageFactorization<T> fact = igusa_image_factorization(F);
  info.B = fact.B;
  info.F0 = fact.F0;
  info.rank = fact.B.rows();
  for (auto &A : ListRepr) {
    if (A.rows() != info.rank) {
      std::cerr << "IGUSA_FACE_VERTICES: the vertex representatives should "
                   "be of size the rank of F, rank="
                << info.rank << " size=" << A.rows() << "\n";
      throw TerminalException{1};
    }
  }
  LinSpaceMatrix<T> LinSpa_r = ComputeCanonicalSpace<T>(info.rank);
  IgusaSpace<T, Tint> space_r =
      build_igusa_space<T, Tint>(LinSpa_r, params, os);
  igusa_insert_initial_cuts(space_r);
  std::vector<std::pair<MyMatrix<T>, size_t>> l_pair =
      igusa_bounded_face_vertices<T, Tint, Tgroup>(space_r, info.F0, b,
                                                   ListRepr, info, os);
  for (auto &pair : l_pair) {
    info.l_vertex.push_back(pair.first);
    info.l_type.push_back(pair.second);
  }
  info.dim_face = igusa_affine_dimension(info.l_vertex);
  return info;
}

inline FullNamelist NAMELIST_GetStandard_IGUSA_FACE_VERTICES() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  {
    std::map<std::string, std::string> ListStringValues_doc;
    std::map<std::string, std::string> ListBoolValues_doc;
    ListStringValues_doc["arithmetic"] =
        doc_rational_integer_arithmetics_by_library();
    ListStringValues_doc["IlpMethod"] = doc_IlpMethod();
    ListStringValues_doc["OutFile"] = doc_StderrStdoutFile("stderr");
    ListStringValues_doc["FileFacet"] = "The file containing the matrix F\n\
(n x n, positive semidefinite) of the face tr(F X) >= rhs";
    ListStringValues_doc["FacetRhs"] = "The right hand side rhs of the face,\n\
a rational number";
    ListStringValues_doc["FileListVertex"] = "The file containing the\n\
representatives of the orbits of vertices of P_r, r the rank of F, in\n\
ListMatrix format";
    ListBoolValues_doc["CompareIncidence"] = "Default: F\n\
For F positive definite: also compute the incidence as IGUSA_FacetIncidence\n\
does and compare";
    ListBlock["DATA"].setListStringValues_doc(ListStringValues_doc);
    ListBlock["DATA"].setListBoolValues_doc(ListBoolValues_doc);
  }
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_IGUSA_IGUSA_FACE_VERTICES_H_
// clang-format on
