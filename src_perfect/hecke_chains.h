// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_HECKE_CHAINS_H_
#define SRC_PERFECT_HECKE_CHAINS_H_

// clang-format off
#include "hecke_common.h"
// clang-format on

/*
  The chains of the perfect form complex for the Hecke operators: the
  elementary operations, the filter of the span of a chain, and the walk
  through the perfect cones to the cone containing a point.
 */

//
// Chains: elementary operations
//

template <typename T, typename Tint>
std::vector<PerfectFaceEntry<T, Tint>>
chain_right_multiply(std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                     T const &scal, MyMatrix<Tint> const &M) {
  std::vector<PerfectFaceEntry<T, Tint>> chain_ret;
  for (auto &fe : chain) {
    T value = fe.value * scal;
    MyMatrix<Tint> Mprod = fe.M * M;
    PerfectFaceEntry<T, Tint> fe_new{fe.iOrb, value, std::move(Mprod)};
    chain_ret.emplace_back(std::move(fe_new));
  }
  return chain_ret;
}

// The vectors of all the cells of a chain.
template <typename T, typename Tint, typename Tgroup>
MyMatrix<Tint>
chain_vector_family(int const &index,
                    std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                    FullComplexEnumeration<T, Tint, Tgroup> const &fce) {
  int n = fce.pctdi.LinSpa.n;
  std::vector<MyVector<Tint>> l_vect;
  for (auto &fe : chain) {
    MyMatrix<Tint> EXT = fce.levels[index].l_faces[fe.iOrb].EXT * fe.M;
    for (int i_row = 0; i_row < EXT.rows(); i_row++) {
      l_vect.push_back(GetMatrixRow(EXT, i_row));
    }
  }
  return MatrixFromVectorFamilyDim(n, l_vect);
}

/*
  The filter selecting the cells whose vectors belong to the span of the
  vectors of the chain. The span is saturated by the matrices of ListComm
  (e.g. the multiplication by omega for Hermitian T-spaces): the cells at
  infinity of the complex live in subspaces stable under them and the
  plain rational span of the vectors would contain no cell at all. If the
  span is everything then no filtering.
 */
template <typename T, typename Tint, typename Tgroup>
std::pair<std::function<bool(MyMatrix<Tint> const &)>, bool>
get_subspace_filter(int const &index,
                    std::vector<PerfectFaceEntry<T, Tint>> const &chain,
                    FullComplexEnumeration<T, Tint, Tgroup> const &fce,
                    [[maybe_unused]] std::ostream &os) {
  int n = fce.pctdi.LinSpa.n;
  MyMatrix<Tint> V_int = chain_vector_family(index, chain, fce);
  MyMatrix<T> V = UniversalMatrixConversion<T, Tint>(V_int);
  std::vector<MyMatrix<T>> const &ListComm = fce.pctdi.LinSpa.ListComm;
  int rnk = RankMat(V);
  while (true) {
    std::vector<MyVector<T>> l_vect;
    for (int i_row = 0; i_row < V.rows(); i_row++) {
      l_vect.push_back(GetMatrixRow(V, i_row));
    }
    for (auto &C : ListComm) {
      MyMatrix<T> Vimg = V * C;
      for (int i_row = 0; i_row < Vimg.rows(); i_row++) {
        l_vect.push_back(GetMatrixRow(Vimg, i_row));
      }
    }
    MyMatrix<T> V_new = MatrixFromVectorFamilyDim(n, l_vect);
    int rnk_new = RankMat(V_new);
    if (rnk_new == rnk) {
      break;
    }
    V = V_new;
    rnk = rnk_new;
  }
#ifdef DEBUG_HECKE_OPERATORS
  os << "HECKE: get_subspace_filter |V|=" << V.rows() << " rnk=" << rnk << " n=" << n << "\n";
#endif
  if (rnk == n) {
    std::function<bool(MyMatrix<Tint> const &)> f_all = [](MyMatrix<Tint> const &) -> bool { return true; };
    return {f_all, false};
  }
  // The rows e of Equa satisfy V e = 0.
  MyMatrix<T> Equa = NullspaceTrMat(V);
  MyMatrix<T> EquaT = Equa.transpose();
  std::function<bool(MyMatrix<Tint> const &)> f_sub = [EquaT](MyMatrix<Tint> const &EXT) -> bool {
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
    MyMatrix<T> prod = EXT_T * EquaT;
    return IsZeroMatrix(prod);
  };
  return {f_sub, true};
}

//
// Walking through the top dimensional cones to reach a cell
//

/*
  The generic search of a cell (fce_face_search) goes through the
  computation of a positive definite form with the cell as set of
  shortest vectors, which is not adapted to the cells at infinity. Here
  we instead walk from a starting cone along a segment towards the
  barycenter of the cell: at every step the exit facet of the segment is
  computed and the adjacent cone is entered. The walk ends with a cone
  whose closure contains the barycenter, which forces the cell to be a
  face of that cone. The starting point is chosen randomly in the
  interior of the starting cone so that the segment is generic.

  All computations are done in the frame of the reference cone: the cone
  EXT_i M is mapped back to EXT_i by w -> w M^{-1} on the vectors, and a
  vector w corresponds to the point of coordinates (w ListMat_k w^T)_k
  of the T-space.
 */
template <typename T, typename Tint, typename Tgroup> struct ConeWalker {
private:
  FullComplexEnumeration<T, Tint, Tgroup> const &fce;
  std::vector<MyMatrix<T>> l_scal;
  std::vector<std::vector<MyVector<T>>> ll_facet_fct;
  std::ostream &os;
  MyVector<T> get_coords(MyMatrix<T> const &W,
                         std::vector<T> const &weights) const {
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    int n_mat = ListMat.size();
    MyVector<T> coords = ZeroVector<T>(n_mat);
    int n_row = W.rows();
    for (int i_row = 0; i_row < n_row; i_row++) {
      MyVector<T> V = GetMatrixRow(W, i_row);
      for (int k = 0; k < n_mat; k++) {
        coords(k) += weights[i_row] * EvaluationQuadForm<T, T>(ListMat[k], V);
      }
    }
    return coords;
  }
public:
  ConeWalker(FullComplexEnumeration<T, Tint, Tgroup> const &_fce, std::ostream &_os)
      : fce(_fce), os(_os) {
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    for (auto &cone : fce.pctdi.l_perfect) {
      MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(ListMat, cone.EXT);
      l_scal.push_back(ScalMat);
      std::vector<MyVector<T>> l_fct;
      for (auto &adj : cone.full_adjacencies(os)) {
        MyMatrix<T> ScalFacet = SelectRow(ScalMat, adj.f_ext);
        MyMatrix<T> NSP = NullspaceTrMat(ScalFacet);
        if (NSP.rows() != 1) {
          std::cerr << "HECKE: The facet should define a functional of dimension 1, |NSP|=" << NSP.rows() << "\n";
          throw TerminalException{1};
        }
        MyVector<T> fct = GetMatrixRow(NSP, 0);
        // Orientation: positive on the vectors outside of the facet.
        int n_row = ScalMat.rows();
        T sign_val(0);
        for (int i_row = 0; i_row < n_row; i_row++) {
          if (adj.f_ext[i_row] == 0) {
            MyVector<T> V = GetMatrixRow(ScalMat, i_row);
            T val = ScalarProduct(V, fct);
            if (val != 0) {
              sign_val = val;
              break;
            }
          }
        }
        if (sign_val == 0) {
          std::cerr << "HECKE: Failed to orient the facet functional\n";
          throw TerminalException{1};
        }
        if (sign_val < 0) {
          fct = -fct;
        }
        l_fct.push_back(fct);
      }
      ll_facet_fct.push_back(l_fct);
    }
  }
  /*
    A top dimensional cone (i_perfect, M) whose closure contains the
    barycenter of the rays of the (rational) vectors of P.
   */
  /*
    Whether the closure of the cone (i_perfect, M) contains the barycenter
    of the rays of the vectors of P. If yes, returns the number of facets
    containing it.
   */
  std::optional<int> n_facet_containing(int const &i_perfect, MyMatrix<Tint> const &M,
                                        MyMatrix<T> const &P) const {
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
    MyMatrix<T> Minv = Inverse(M_T);
    std::vector<T> weights_p(P.rows(), T(1));
    MyVector<T> cp = get_coords(P * Minv, weights_p);
    int n_tight = 0;
    for (auto &fct : ll_facet_fct[i_perfect]) {
      T b = ScalarProduct(fct, cp);
      if (b < 0) {
        return {};
      }
      if (b == 0) {
        n_tight += 1;
      }
    }
    return n_tight;
  }
  TopPerfectCone<Tint> find_containing_cone(MyMatrix<Tint> const &EXT_target) const {
    MyMatrix<T> P = UniversalMatrixConversion<T, Tint>(EXT_target);
    return find_containing_cone_T(P);
  }
  TopPerfectCone<Tint> find_containing_cone_T(MyMatrix<T> const &P) const {
    int n = fce.pctdi.LinSpa.n;
    std::vector<T> weights_p(P.rows(), T(1));
    MyMatrix<T> Q0 = UniversalMatrixConversion<T, Tint>(fce.pctdi.l_perfect[0].EXT);
    size_t n_attempt = 0;
    uint64_t seed = 12345;
    auto next_random = [&]() -> uint64_t {
      seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
      return seed >> 33;
    };
    while (true) {
      n_attempt += 1;
      if (n_attempt > 20) {
        throw HeckeFailure{"find_containing_cone: the walk cycled in 20 attempts"};
      }
      std::vector<T> weights_q(Q0.rows());
      for (int i_row = 0; i_row < Q0.rows(); i_row++) {
        int w_val = 1 + static_cast<int>(next_random() % 1000);
        weights_q[i_row] = T(w_val);
      }
      int i_cone = 0;
      MyMatrix<Tint> M = IdentityMat<Tint>(n);
      std::unordered_set<MyMatrix<Tint>> set_visited;
      bool cycle = false;
      HeckeBudget budget;
      while (true) {
        if (set_visited.size() > budget.max_walk_steps) {
          throw HeckeFailure{"find_containing_cone: walk longer than " + std::to_string(budget.max_walk_steps) + " steps"};
        }
        MyMatrix<Tint> EXT_cone = fce.pctdi.l_perfect[i_cone].EXT * M;
        MyMatrix<Tint> key = tot_set(EXT_cone);
        if (!set_visited.insert(key).second) {
          cycle = true;
          break;
        }
        MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
        MyMatrix<T> Minv = Inverse(M_T);
        MyVector<T> cq = get_coords(Q0 * Minv, weights_q);
        MyVector<T> cp = get_coords(P * Minv, weights_p);
        std::optional<T> best_t;
        int best_f = -1;
        std::vector<MyVector<T>> const &l_fct = ll_facet_fct[i_cone];
        int n_facet = l_fct.size();
        for (int i_f = 0; i_f < n_facet; i_f++) {
          T a = ScalarProduct(l_fct[i_f], cq);
          T b = ScalarProduct(l_fct[i_f], cp);
          if (b < 0) {
            T t_f = a / (a - b);
            if (!best_t || t_f < *best_t) {
              best_t = t_f;
              best_f = i_f;
            }
          }
        }
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: walk i_cone=" << i_cone << " best_f=" << best_f;
        if (best_t) {
          os << " best_t=" << *best_t;
        }
        os << " cp=" << StringVectorGAP(cp) << " cq=" << StringVectorGAP(cq) << "\n";
#endif
        if (best_f == -1) {
          return {i_cone, M};
        }
        sing_adj<Tint> const &adj = fce.pctdi.l_perfect[i_cone].full_adjacencies(os)[best_f];
        M = adj.eMat * M;
        i_cone = adj.jCone;
      }
      if (cycle) {
#ifdef DEBUG_HECKE_OPERATORS
        os << "HECKE: find_containing_cone, cycle detected, new attempt\n";
#endif
      }
    }
  }
  /*
    The vertex of the complex on the ray of the functional
    Q -> sum_w Q[w] for the rational vectors w of W. The vectors of the
    vertex are the vectors of the containing cone whose functional is
    proportional to it (they are not necessarily proportional to the
    rows of W over Z, e.g. for Hermitian T-spaces).
   */
  FceFaceSearch<Tint> search_vertex(MyMatrix<T> const &W) const {
    TopPerfectCone<Tint> tpc = find_containing_cone_T(W);
    int i_cone = tpc.i_perfect;
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(tpc.M);
    MyMatrix<T> Minv = Inverse(M_T);
    std::vector<T> weights(W.rows(), T(1));
    MyVector<T> p = get_coords(W * Minv, weights);
    int n_mat = p.size();
    int l_nz = -1;
    for (int k = 0; k < n_mat; k++) {
      if (p(k) != 0) {
        l_nz = k;
        break;
      }
    }
    if (l_nz == -1) {
      std::cerr << "HECKE: The functional of the vertex should be nonzero\n";
      throw TerminalException{1};
    }
    MyMatrix<T> const &ScalMat = l_scal[i_cone];
    int n_row = ScalMat.rows();
    Face f_ext(n_row);
    for (int i_row = 0; i_row < n_row; i_row++) {
      bool is_prop = true;
      for (int k = 0; k < n_mat; k++) {
        if (ScalMat(i_row, k) * p(l_nz) != p(k) * ScalMat(i_row, l_nz)) {
          is_prop = false;
          break;
        }
      }
      if (is_prop) {
        f_ext[i_row] = 1;
      }
    }
    if (f_ext.count() == 0) {
      std::cerr << "HECKE: No vector of the containing cone is on the ray of the vertex\n";
      throw TerminalException{1};
    }
    triple<Tint> t{static_cast<size_t>(i_cone), f_ext, tpc.M};
    triple<Tint> t_can = canonicalize_triple(fce.pctdi.l_perfect, t, os);
    int vlevel = fce.levels.size() - 1;
    int n_orb = fce.levels[vlevel].l_faces.size();
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      std::vector<triple<Tint>> const &l_triple = fce.levels[vlevel].l_faces[iOrb].l_triple;
      std::optional<MyMatrix<Tint>> opt = test_triple_in_listtriple(fce.pctdi.l_perfect, l_triple, t_can, os);
      if (opt) {
        return {vlevel, iOrb, *opt};
      }
    }
    std::cerr << "HECKE: Failed to find the vertex among the orbits\n";
    throw TerminalException{1};
  }
  /*
    The smallest face of the complex whose relative interior contains the
    barycenter of the rays of the rational vectors of W: the walk gives a
    cone containing it, the tight facets of the cone cut out the face.
   */
  FceFaceSearch<Tint> search_face_containing(MyMatrix<T> const &W) const {
    TopPerfectCone<Tint> tpc = find_containing_cone_T(W);
    int i_cone = tpc.i_perfect;
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(tpc.M);
    MyMatrix<T> Minv = Inverse(M_T);
    std::vector<T> weights(W.rows(), T(1));
    MyVector<T> p = get_coords(W * Minv, weights);
    std::vector<sing_adj<Tint>> const &l_adj = fce.pctdi.l_perfect[i_cone].full_adjacencies(os);
    std::vector<MyVector<T>> const &l_fct = ll_facet_fct[i_cone];
    int n_row = fce.pctdi.l_perfect[i_cone].EXT.rows();
    Face f_ext(n_row);
    f_ext.set();
    for (size_t i_f = 0; i_f < l_fct.size(); i_f++) {
      T val = ScalarProduct(l_fct[i_f], p);
      if (val == 0) {
        f_ext &= l_adj[i_f].f_ext;
      }
    }
    triple<Tint> t{static_cast<size_t>(i_cone), f_ext, tpc.M};
    triple<Tint> t_can = canonicalize_triple(fce.pctdi.l_perfect, t, os);
    MyMatrix<Tint> EXT_face = SelectRow(fce.pctdi.l_perfect[i_cone].EXT, f_ext);
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    int n_mat = ListMat.size();
    MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(ListMat, EXT_face);
    int index = n_mat - RankMat(ScalMat);
    int n_orb = fce.levels[index].l_faces.size();
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      std::vector<triple<Tint>> const &l_triple = fce.levels[index].l_faces[iOrb].l_triple;
      std::optional<MyMatrix<Tint>> opt = test_triple_in_listtriple(fce.pctdi.l_perfect, l_triple, t_can, os);
      if (opt) {
        return {index, iOrb, *opt};
      }
    }
    std::cerr << "HECKE: Failed to find the face among the orbits of level " << index << "\n";
    throw TerminalException{1};
  }
  /*
    The position of a cell in the complex, by walking to a cone
    containing it.
   */
  FceFaceSearch<Tint> search_cell(MyMatrix<Tint> const &EXT_cell) const {
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: search_cell EXT_cell=" << StringMatrixGAP(EXT_cell) << "\n";
#endif
    TopPerfectCone<Tint> tpc = find_containing_cone(EXT_cell);
    MyMatrix<Tint> EXT_cone = fce.pctdi.l_perfect[tpc.i_perfect].EXT * tpc.M;
#ifdef DEBUG_HECKE_OPERATORS
    os << "HECKE: search_cell i_perfect=" << tpc.i_perfect << " M=" << StringMatrixGAP(tpc.M) << " EXT_cone=" << StringMatrixGAP(EXT_cone) << "\n";
#endif
    ContainerMatrix<Tint> cont(EXT_cone);
    Face f_ext(EXT_cone.rows());
    for (int i_row = 0; i_row < EXT_cell.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(EXT_cell, i_row);
      std::optional<size_t> opt = cont.GetIdx_v(V);
      if (!opt) {
        std::cerr << "HECKE: The vector of the cell is not a vector of the containing cone\n";
        throw TerminalException{1};
      }
      f_ext[*opt] = 1;
    }
    triple<Tint> t{static_cast<size_t>(tpc.i_perfect), f_ext, tpc.M};
    triple<Tint> t_can = canonicalize_triple(fce.pctdi.l_perfect, t, os);
    std::vector<MyMatrix<T>> const &ListMat = fce.pctdi.LinSpa.ListMat;
    int n_mat = ListMat.size();
    MyMatrix<T> ScalMat = get_scal_mat<T, Tint>(ListMat, EXT_cell);
    int index = n_mat - RankMat(ScalMat);
    int n_orb = fce.levels[index].l_faces.size();
    for (int iOrb = 0; iOrb < n_orb; iOrb++) {
      std::vector<triple<Tint>> const &l_triple = fce.levels[index].l_faces[iOrb].l_triple;
      std::optional<MyMatrix<Tint>> opt = test_triple_in_listtriple(fce.pctdi.l_perfect, l_triple, t_can, os);
      if (opt) {
        return {index, iOrb, *opt};
      }
    }
    std::cerr << "HECKE: Failed to find the cell among the orbits of level " << index << "\n";
    throw TerminalException{1};
  }
};

// clang-format off
#endif  // SRC_PERFECT_HECKE_CHAINS_H_
// clang-format on
