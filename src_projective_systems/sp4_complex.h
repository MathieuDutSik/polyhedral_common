// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PROJECTIVE_SYSTEMS_SP4_COMPLEX_H_
#define SRC_PROJECTIVE_SYSTEMS_SP4_COMPLEX_H_

// clang-format off
#include "hecke_operators.h"
#include "PolytopeEquiStab.h"
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SP4_COMPLEX
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SP4_COMPLEX
#endif

#ifdef TIMINGS
#define TIMINGS_SP4_COMPLEX
#endif

/*
  The cell complex of Sp_4(Z) of the projective systems (port of
  ProjectiveSystem.g of MyPolyhedral, the functions used by the driver
  HeckeOperatorSp4/TheBigComputation.g).

  The cells are the configurations of vectors of Z^4 (taken up to sign)
  whose cones of rank one forms v v^T meet the Siegel space: the
  intersection of the Voronoi decomposition of the positive definite
  forms of rank 4 with the Siegel space, of dimension 6. The cells of
  the top level are the perfect cones of A_4 (10 vectors) and D_4 (12
  vectors), the cells of the other levels are the subsets of their
  vectors that are equivalent under Sp_4(Z) to one of the listed
  representative configurations (the projective systems), down to a
  single vector. A cell whose vectors span Z^4 meets the interior of the
  cone of the positive definite forms ("well rounded"), the others are at
  infinity.

  The levels are indexed as in src_perfect: the top cells at the index 0
  and the faces of a cell of the index k at the index k + 1 (the GAP code
  has the reverse order, its level L is the index 7 - L here). The
  orientation of a cell is combinatorial: the sign character of its
  stabilizer, determined from the signs of the boundary (the sign of the
  permutation of the two vectors for the 2-vector cells), as in
  SP4Z_GetCellComplex.

  The symplectic form is the standard J of SYMPL_GetSymplecticForm(2):
  J[i][2 + i] = 1, J[2 + i][i] = -1.
 */

namespace sp4 {

using Tint = mpz_class;
using T = mpq_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, mpz_class>;

inline MyMatrix<Tint> symplectic_form() {
  MyMatrix<Tint> J = ZeroMatrix<Tint>(4, 4);
  for (int i = 0; i < 2; i++) {
    J(i, 2 + i) = 1;
    J(2 + i, i) = -1;
  }
  return J;
}

// Whether g is in Sp_4(Z): integral with g J g^T = J.
template <typename Tf> bool is_in_sp4(MyMatrix<Tf> const &g) {
  if (!IsIntegralMatrix(g)) {
    return false;
  }
  MyMatrix<Tf> J = UniversalMatrixConversion<Tf, Tint>(symplectic_form());
  return g * J * g.transpose() == J;
}

// The canonical form of a vector up to sign: the first nonzero entry is
// positive.
inline MyVector<Tint> sign_normalize(MyVector<Tint> const &v) {
  for (int i = 0; i < v.size(); i++) {
    if (v(i) != 0) {
      if (v(i) < 0) {
        return -v;
      }
      return v;
    }
  }
  return v;
}

// The key of a configuration taken up to the signs of its vectors and to
// their order: the sorted list of the sign normalized vectors.
using ConfigKey = std::vector<int64_t>;

inline ConfigKey config_key(MyMatrix<Tint> const &EXT) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  std::vector<std::vector<int64_t>> l_row;
  for (int i = 0; i < n_row; i++) {
    MyVector<Tint> v = sign_normalize(GetMatrixRow(EXT, i));
    std::vector<int64_t> w(n_col);
    for (int j = 0; j < n_col; j++) {
      w[j] = v(j).get_si();
    }
    l_row.emplace_back(std::move(w));
  }
  std::sort(l_row.begin(), l_row.end());
  ConfigKey key;
  for (auto &w : l_row) {
    key.insert(key.end(), w.begin(), w.end());
  }
  return key;
}

// The position of v or -v among the rows of EXT.
inline std::optional<int> position_antipodal(MyMatrix<Tint> const &EXT, MyVector<Tint> const &v) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  for (int i = 0; i < n_row; i++) {
    bool same = true, opposite = true;
    for (int j = 0; j < n_col; j++) {
      if (EXT(i, j) != v(j)) {
        same = false;
      }
      if (EXT(i, j) != -v(j)) {
        opposite = false;
      }
    }
    if (same || opposite) {
      return i;
    }
  }
  return {};
}

// The permutation of the rows of EXT (up to sign) induced by M, if M
// maps the configuration to itself.
inline std::optional<std::vector<Tidx>> permutation_antipodal(MyMatrix<Tint> const &EXT, MyMatrix<Tint> const &M) {
  MyMatrix<Tint> EXTimg = EXT * M;
  std::vector<Tidx> l_pos;
  for (int i = 0; i < EXT.rows(); i++) {
    std::optional<int> opt = position_antipodal(EXT, GetMatrixRow(EXTimg, i));
    if (!opt) {
      return {};
    }
    l_pos.push_back(*opt);
  }
  return l_pos;
}

// The configurations of EXT * M for M given by the subset of rows.
inline MyMatrix<Tint> select_rows(MyMatrix<Tint> const &EXT, std::vector<int> const &l_row) {
  MyMatrix<Tint> M(l_row.size(), EXT.cols());
  for (size_t i = 0; i < l_row.size(); i++) {
    for (int j = 0; j < EXT.cols(); j++) {
      M(i, j) = EXT(l_row[i], j);
    }
  }
  return M;
}

// (EXT; -EXT)
inline MyMatrix<Tint> totality(MyMatrix<Tint> const &EXT) {
  int n_row = EXT.rows();
  MyMatrix<Tint> EXTtot(2 * n_row, EXT.cols());
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < EXT.cols(); j++) {
      EXTtot(i, j) = EXT(i, j);
      EXTtot(n_row + i, j) = -EXT(i, j);
    }
  }
  return EXTtot;
}

/*
  The representative configurations of the levels (SP4Z_GetList_Faces and
  SP4Z_GetList_Faces_Coho), in the order of the indices here: the top
  cells A_4 and D_4 first, the single vector last. The configurations of
  GAP are multiplied by the permutation matrix exchanging the first two
  coordinates, which turns the form of SP4Z_GetList_Faces_Coho into J.
 */
inline std::vector<std::vector<MyMatrix<Tint>>> get_representative_configurations() {
  auto f_conf = [](std::vector<std::vector<int>> const &l_vect) -> MyMatrix<Tint> {
    MyMatrix<Tint> M(l_vect.size(), 4);
    for (size_t i = 0; i < l_vect.size(); i++) {
      // The exchange of the first two coordinates.
      M(i, 0) = l_vect[i][1];
      M(i, 1) = l_vect[i][0];
      M(i, 2) = l_vect[i][2];
      M(i, 3) = l_vect[i][3];
    }
    return M;
  };
  std::vector<std::vector<int>> id{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  auto f_ext = [&](std::vector<std::vector<int>> const &l_add) -> MyMatrix<Tint> {
    std::vector<std::vector<int>> l = id;
    l.insert(l.end(), l_add.begin(), l_add.end());
    return f_conf(l);
  };
  MyMatrix<Tint> face0_1 = f_ext({{1, 1, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, -1}, {0, 1, 0, 1}, {1, 1, 1, 0}, {1, 0, 1, -1}});
  MyMatrix<Tint> face0_2 = f_ext({{1, 1, 0, 0}, {1, 0, -1, 0}, {0, 1, 0, -1}, {0, 0, 1, -1},
                                  {1, 1, -1, 0}, {1, 1, 0, -1}, {1, 0, -1, 1}, {0, 1, 1, -1}});
  MyMatrix<Tint> face1_1 = f_ext({{1, 1, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, -1}, {0, 1, 0, 1}, {1, 0, 1, -1}});
  MyMatrix<Tint> face1_2 = f_ext({{1, 1, 0, 0}, {0, 1, 0, -1}, {0, 0, 1, -1}, {1, 1, -1, 0}, {1, 0, -1, 1}});
  MyMatrix<Tint> face2_1 = f_ext({{1, 0, 0, 1}, {0, 1, 1, 0}});
  MyMatrix<Tint> face2_2 = f_ext({{1, 1, 0, 0}, {1, 0, 0, -1}, {0, 1, 0, 1}, {1, 0, 1, -1}});
  MyMatrix<Tint> face2_3 = f_ext({{1, 1, 0, 0}, {1, 0, 1, 0}, {0, 1, 0, 1}, {1, 0, 1, -1}});
  MyMatrix<Tint> face3_1 = f_ext({{1, 0, 0, 1}});
  MyMatrix<Tint> face3_2 = f_ext({{1, 0, 1, 0}, {0, 1, 0, 1}});
  MyMatrix<Tint> face3_3 = f_ext({{1, 1, 0, 0}, {0, 1, 0, 1}, {1, 0, 1, -1}});
  MyMatrix<Tint> face4_1 = f_ext({});
  MyMatrix<Tint> face4_2 = f_conf({{1, 0, 0, 0}, {0, 1, 0, 0}, {1, 1, 0, 0}});
  MyMatrix<Tint> face5_1 = f_conf({{1, 0, 0, 0}, {0, 1, 0, 0}});
  MyMatrix<Tint> face6_1 = f_conf({{1, 0, 0, 0}});
  return {{face0_1, face0_2}, {face1_1, face1_2}, {face2_1, face2_2, face2_3}, {face3_1, face3_2, face3_3},
          {face4_1, face4_2}, {face5_1}, {face6_1}};
}

/*
  The elements of Sp_4(Z) preserving a configuration of rank 4, as the
  list of all of them (the group is finite) and the corresponding
  permutations of the rows of (EXT; -EXT) (SymplecticStabilizer).
 */
struct FiniteStabilizer {
  std::vector<MyMatrix<Tint>> l_elt;
  std::vector<MyMatrix<Tint>> l_gens;
};

inline FiniteStabilizer symplectic_stabilizer(MyMatrix<Tint> const &EXT, std::ostream &os) {
  MyMatrix<Tint> EXTtot = totality(EXT);
  MyMatrix<T> EXTtot_T = UniversalMatrixConversion<T, Tint>(EXTtot);
  Tgroup GRP = LinPolytope_Automorphism<T, Tgroup>(EXTtot_T, os);
  FiniteStabilizer fs;
  std::vector<Telt> l_perm_gens;
  for (auto &elt : GRP) {
    MyMatrix<T> M_T = FindTransformation(EXTtot_T, EXTtot_T, elt);
    if (is_in_sp4(M_T)) {
      fs.l_elt.push_back(UniversalMatrixConversion<Tint, T>(M_T));
    }
  }
  // A generating set: the elements not in the group generated by the
  // previous ones (the groups are small).
  std::set<MyMatrix<Tint>> set_generated;
  MyMatrix<Tint> id = IdentityMat<Tint>(4);
  set_generated.insert(id);
  for (auto &M : fs.l_elt) {
    if (set_generated.count(M) > 0) {
      continue;
    }
    fs.l_gens.push_back(M);
    std::vector<MyMatrix<Tint>> l_new(set_generated.begin(), set_generated.end());
    size_t pos = 0;
    while (pos < l_new.size()) {
      for (auto &g : fs.l_gens) {
        MyMatrix<Tint> h = l_new[pos] * g;
        if (set_generated.insert(h).second) {
          l_new.push_back(h);
        }
      }
      pos += 1;
    }
  }
  return fs;
}

// The subsets of k elements of [0, n).
inline std::vector<std::vector<int>> combinations(int const &n, int const &k) {
  std::vector<std::vector<int>> l_comb;
  std::vector<int> sel(k);
  std::function<void(int, int)> f_rec = [&](int start, int depth) -> void {
    if (depth == k) {
      l_comb.push_back(sel);
      return;
    }
    for (int i = start; i < n; i++) {
      sel[depth] = i;
      f_rec(i + 1, depth + 1);
    }
  };
  f_rec(0, 0);
  return l_comb;
}

// The positions of rows of EXT forming a basis of its span (greedy).
inline std::vector<int> basis_rows(MyMatrix<Tint> const &EXT) {
  std::vector<int> l_idx;
  int rnk = 0;
  for (int i = 0; i < EXT.rows(); i++) {
    std::vector<int> l_try = l_idx;
    l_try.push_back(i);
    if (RankMat(select_rows(EXT, l_try)) > rnk) {
      l_idx = l_try;
      rnk += 1;
    }
  }
  return l_idx;
}

/*
  The elements M of Sp_4(Z) with eFace1 * M contained in eFace2 (up to
  the signs), for eFace1 of rank 4, one per image (ListEmbedding). With
  B a basis of rows of eFace1 and S its image, M = B^{-1} S is in Sp_4(Q)
  iff S J S^T = B J B^T: the rows of S are chosen among +-eFace2 one at a
  time with the symplectic pairings of B (a backtracking, instead of all
  the ordered choices of GAP), then M must be integral and map eFace1 into
  +-eFace2.
 */
inline std::vector<MyMatrix<Tint>> list_embedding(MyMatrix<Tint> const &eFace1, MyMatrix<Tint> const &eFace2) {
  MyMatrix<Tint> J = symplectic_form();
  std::vector<int> l_basis = basis_rows(eFace1);
  MyMatrix<Tint> B = select_rows(eFace1, l_basis);
  MyMatrix<Tint> GB = B * J * B.transpose();
  MyMatrix<T> B_inv = Inverse(UniversalMatrixConversion<T, Tint>(B));
  MyMatrix<Tint> eFace2tot = totality(eFace2);
  int nb2tot = eFace2tot.rows();
  // The symplectic pairings of the rows of +-eFace2.
  MyMatrix<Tint> G2 = eFace2tot * J * eFace2tot.transpose();
  std::vector<MyMatrix<Tint>> l_sol;
  std::set<ConfigKey> set_img;
  std::vector<int> sel(4);
  std::function<void(int)> f_rec = [&](int k) -> void {
    if (k == 4) {
      MyMatrix<Tint> S = select_rows(eFace2tot, sel);
      MyMatrix<T> M_T = B_inv * UniversalMatrixConversion<T, Tint>(S);
      if (!IsIntegralMatrix(M_T)) {
        return;
      }
      MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
      MyMatrix<Tint> EXTimg = eFace1 * M;
      for (int i = 0; i < EXTimg.rows(); i++) {
        if (!position_antipodal(eFace2, GetMatrixRow(EXTimg, i))) {
          return;
        }
      }
      if (set_img.insert(config_key(EXTimg)).second) {
        l_sol.push_back(M);
      }
      return;
    }
    for (int c = 0; c < nb2tot; c++) {
      bool is_ok = true;
      for (int l = 0; l < k; l++) {
        if (sel[l] == c || G2(c, sel[l]) != GB(k, l)) {
          is_ok = false;
          break;
        }
      }
      if (is_ok) {
        sel[k] = c;
        f_rec(k + 1);
      }
    }
  };
  f_rec(0);
  return l_sol;
}

// The matrices M with eFace1 contained in eFace2 * M, one per eFace2 * M
// (ListEmbeddingRev).
inline std::vector<MyMatrix<Tint>> list_embedding_rev(MyMatrix<Tint> const &eFace1, MyMatrix<Tint> const &eFace2) {
  std::vector<MyMatrix<Tint>> l_ret;
  std::set<ConfigKey> set_key;
  for (auto &M : list_embedding(eFace1, eFace2)) {
    MyMatrix<Tint> M_inv = Inverse(M);
    if (set_key.insert(config_key(eFace2 * M_inv)).second) {
      l_ret.push_back(M_inv);
    }
  }
  return l_ret;
}

// Whether the rows of B are pairwise orthogonal for J.
inline bool is_isotropic(MyMatrix<Tint> const &B) {
  MyMatrix<Tint> J = symplectic_form();
  for (int i = 0; i < B.rows(); i++) {
    for (int j = i + 1; j < B.rows(); j++) {
      MyVector<Tint> vi = GetMatrixRow(B, i);
      MyVector<Tint> vj = GetMatrixRow(B, j);
      if (vi.dot(J * vj) != 0) {
        return false;
      }
    }
  }
  return true;
}

// Whether the rows of B (independent) span a primitive sublattice of Z^4:
// the gcd of the maximal minors is 1.
inline bool is_primitive_basis(MyMatrix<Tint> const &B) {
  int r = B.rows();
  Tint g = 0;
  for (auto &cols : combinations(4, r)) {
    MyMatrix<Tint> sub(r, r);
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < r; j++) {
        sub(i, j) = B(i, cols[j]);
      }
    }
    Tint det = DeterminantMat(sub);
    g = GcdPair(g, det);
  }
  return g == 1 || g == -1;
}

/*
  The subsets of the rows of eFace2 (a top cell) that are images of eFace1
  by an element of Sp_4(Z) (SP4Z_GetAllEmbeddings), as sorted lists of
  positions. For eFace1 of rank 4, by list_embedding. Otherwise eFace1 spans
  a Lagrangian subspace: a basis of its span (from its rows) is mapped to
  rows of +-eFace2 spanning a primitive isotropic sublattice, and the
  linear map must map eFace1 into eFace2 (up to sign). Sp_4(Z) is
  transitive on the bases of the primitive isotropic sublattices of a
  given rank, so such a map extends to Sp_4(Z).
 */
inline std::vector<std::vector<int>> get_all_embeddings(MyMatrix<Tint> const &eFace1, MyMatrix<Tint> const &eFace2) {
  std::set<std::vector<int>> set_img;
  auto f_image_set = [&](MyMatrix<Tint> const &EXTimg) -> std::optional<std::vector<int>> {
    std::vector<int> l_pos;
    for (int i = 0; i < EXTimg.rows(); i++) {
      std::optional<int> opt = position_antipodal(eFace2, GetMatrixRow(EXTimg, i));
      if (!opt) {
        return {};
      }
      l_pos.push_back(*opt);
    }
    std::sort(l_pos.begin(), l_pos.end());
    return l_pos;
  };
  int rnk = RankMat(eFace1);
  if (rnk == 4) {
    for (auto &M : list_embedding(eFace1, eFace2)) {
      std::optional<std::vector<int>> opt = f_image_set(eFace1 * M);
      set_img.insert(*opt);
    }
    return std::vector<std::vector<int>>(set_img.begin(), set_img.end());
  }
  MyMatrix<T> eFace1_T = UniversalMatrixConversion<T, Tint>(eFace1);
  MyMatrix<T> basis = RowReduction(eFace1_T);
  MyMatrix<Tint> basis_int = UniversalMatrixConversion<Tint, T>(basis);
  if (!is_isotropic(basis_int) || !is_primitive_basis(basis_int)) {
    std::cerr << "SP4: get_all_embeddings, eFace1 should be a basis of a primitive Lagrangian sublattice\n";
    throw TerminalException{1};
  }
  // The coordinates of the rows of eFace1 in the basis.
  MyMatrix<T> coord(eFace1.rows(), rnk);
  for (int i = 0; i < eFace1.rows(); i++) {
    std::optional<MyVector<T>> opt = SolutionMat(basis, GetMatrixRow(eFace1_T, i));
    MyVector<T> c = unfold_opt(opt, "SP4: the rows should be in the span of the basis");
    for (int j = 0; j < rnk; j++) {
      coord(i, j) = c(j);
    }
  }
  MyMatrix<Tint> eFace2tot = totality(eFace2);
  int nb2tot = eFace2tot.rows();
  for (auto &pset : combinations(nb2tot, rnk)) {
    std::vector<int> perm = pset;
    do {
      MyMatrix<Tint> img = select_rows(eFace2tot, perm);
      if (RankMat(img) != rnk || !is_isotropic(img) || !is_primitive_basis(img)) {
        continue;
      }
      MyMatrix<T> EXTimg_T = coord * UniversalMatrixConversion<T, Tint>(img);
      if (!IsIntegralMatrix(EXTimg_T)) {
        continue;
      }
      std::optional<std::vector<int>> opt = f_image_set(UniversalMatrixConversion<Tint, T>(EXTimg_T));
      if (opt) {
        set_img.insert(*opt);
      }
    } while (std::next_permutation(perm.begin(), perm.end()));
  }
  return std::vector<std::vector<int>>(set_img.begin(), set_img.end());
}

/*
  The neighbors of the representative top cells across their ridges
  (SP4Z_GetFlippingInfo): for the top cell iPerf, the ridge (of the index
  1) given by eSet, the neighbor is the top cell jPerf mapped by eMatB.
 */
struct FlipInfo {
  int iRidge;
  std::vector<int> eSet;
  int jPerf;
  MyMatrix<Tint> eMatB;
};

inline std::vector<std::vector<FlipInfo>>
get_flipping_info(std::vector<std::vector<MyMatrix<Tint>>> const &l_conf, [[maybe_unused]] std::ostream &os) {
  std::vector<MyMatrix<Tint>> const &l_perf = l_conf[0];
  std::vector<MyMatrix<Tint>> const &l_ridge = l_conf[1];
  std::vector<std::vector<FlipInfo>> ll_adj;
  for (size_t iPerf = 0; iPerf < l_perf.size(); iPerf++) {
    MyMatrix<Tint> const &ePerf = l_perf[iPerf];
    std::vector<FlipInfo> l_adj;
    for (size_t iRidge = 0; iRidge < l_ridge.size(); iRidge++) {
      for (auto &eMat : list_embedding(l_ridge[iRidge], ePerf)) {
        MyMatrix<Tint> eFace1 = l_ridge[iRidge] * eMat;
        std::vector<int> eSet;
        for (int i = 0; i < eFace1.rows(); i++) {
          eSet.push_back(*position_antipodal(ePerf, GetMatrixRow(eFace1, i)));
        }
        std::sort(eSet.begin(), eSet.end());
        ConfigKey key_perf = config_key(ePerf);
        std::vector<std::pair<int, MyMatrix<Tint>>> l_cont;
        for (size_t jPerf = 0; jPerf < l_perf.size(); jPerf++) {
          for (auto &eMatB : list_embedding_rev(eFace1, l_perf[jPerf])) {
            if (config_key(l_perf[jPerf] * eMatB) != key_perf) {
              l_cont.push_back({jPerf, eMatB});
            }
          }
        }
        if (l_cont.size() != 1) {
          std::cerr << "SP4: get_flipping_info, the ridge should have one neighbor, |l_cont|=" << l_cont.size()
                    << "\n";
          throw TerminalException{1};
        }
        l_adj.push_back({static_cast<int>(iRidge), eSet, l_cont[0].first, l_cont[0].second});
      }
    }
#ifdef DEBUG_SP4_COMPLEX
    os << "SP4: iPerf=" << iPerf << " |l_adj|=" << l_adj.size() << "\n";
#endif
    ll_adj.emplace_back(std::move(l_adj));
  }
  return ll_adj;
}

/*
  A top cell of the complex of Sp_4(Z): the representative top cell iPerf
  mapped by g, with its stabilizer in Sp_4(Z) (the elements, finite), the
  permutation group on the rows of (EXT; -EXT) and its neighbors.
 */
struct TopCell {
  int iPerf;
  MyMatrix<Tint> g;
  MyMatrix<Tint> EXT;
  std::vector<MyMatrix<Tint>> l_stab_gens;
  Tgroup perm_grp;
  struct Adj {
    int iCell;
    MyMatrix<Tint> eEquiv;
    std::vector<int> eSet;
  };
  std::vector<Adj> l_adj;
  // Rows of (EXT; -EXT) forming a basis, and the inverse of that basis.
  std::vector<int> l_basis;
  MyMatrix<T> basis_inv;
  MyMatrix<Tint> EXTtot;
  // The matrix of a permutation of the rows of (EXT; -EXT) in perm_grp:
  // the basis rows go to their images.
  MyMatrix<Tint> matrix_of(Telt const &elt) const {
    std::vector<int> l_img;
    for (auto &i : l_basis) {
      l_img.push_back(elt.at(i));
    }
    MyMatrix<T> M_T = basis_inv * UniversalMatrixConversion<T, Tint>(select_rows(EXTtot, l_img));
    return UniversalMatrixConversion<Tint, T>(M_T);
  }
};

// The permutation of the rows of (EXT; -EXT) induced by M.
inline Telt permutation_totality(MyMatrix<Tint> const &EXT, MyMatrix<Tint> const &M) {
  MyMatrix<Tint> EXTtot = totality(EXT);
  MyMatrix<Tint> EXTimg = EXTtot * M;
  int n = EXTtot.rows();
  std::vector<Tidx> l_pos(n);
  for (int i = 0; i < n; i++) {
    MyVector<Tint> v = GetMatrixRow(EXTimg, i);
    int pos = -1;
    for (int k = 0; k < n; k++) {
      if (GetMatrixRow(EXTtot, k) == v) {
        pos = k;
      }
    }
    if (pos == -1) {
      std::cerr << "SP4: permutation_totality, M does not preserve the configuration\n";
      throw TerminalException{1};
    }
    l_pos[i] = pos;
  }
  return Telt(l_pos);
}

/*
  The orbits of top cells under Sp_4(Z) (SP4Z_GetTopCells with the group
  Sp_4(Z) itself), by flipping across the ridges, and the generators of
  Sp_4(Z) found on the way (the stabilizers and the equivalences).
 */
struct TopCellResult {
  std::vector<TopCell> l_top;
  std::vector<MyMatrix<Tint>> l_gens_G;
};

inline TopCellResult get_top_cells(std::vector<std::vector<MyMatrix<Tint>>> const &l_conf,
                                   std::vector<std::vector<FlipInfo>> const &ll_adj, std::ostream &os) {
  std::vector<MyMatrix<Tint>> const &l_perf = l_conf[0];
  std::vector<FiniteStabilizer> l_stab;
  for (auto &ePerf : l_perf) {
    l_stab.push_back(symplectic_stabilizer(ePerf, os));
  }
  TopCellResult res;
  std::set<MyMatrix<Tint>> set_gens;
  MyMatrix<Tint> id = IdentityMat<Tint>(4);
  auto f_insert_gen = [&](MyMatrix<Tint> const &M) -> void {
    if (M != id && set_gens.insert(M).second) {
      res.l_gens_G.push_back(M);
    }
  };
  // An equivalence g1^{-1} s g2 (s in the stabilizer of the representative)
  // between two top cells of the same type: here all of them are in Sp_4(Z).
  auto f_equivalence = [&](TopCell const &c1, int const &iPerf, MyMatrix<Tint> const &g) -> MyMatrix<Tint> {
    if (c1.iPerf != iPerf) {
      return MyMatrix<Tint>(0, 0);
    }
    return Inverse(c1.g) * l_stab[iPerf].l_elt[0] * g;
  };
  struct InsertResult {
    int iCell;
    MyMatrix<Tint> eEquiv;
  };
  auto f_insert = [&](int const &iPerf, MyMatrix<Tint> const &g) -> InsertResult {
    for (size_t iCell = 0; iCell < res.l_top.size(); iCell++) {
      MyMatrix<Tint> test = f_equivalence(res.l_top[iCell], iPerf, g);
      if (test.rows() > 0) {
        f_insert_gen(test);
        return {static_cast<int>(iCell), test};
      }
    }
    TopCell cell;
    cell.iPerf = iPerf;
    cell.g = g;
    cell.EXT = l_perf[iPerf] * g;
    MyMatrix<Tint> g_inv = Inverse(g);
    std::vector<Telt> l_perm;
    for (auto &s : l_stab[iPerf].l_gens) {
      MyMatrix<Tint> s_conj = g_inv * s * g;
      cell.l_stab_gens.push_back(s_conj);
      f_insert_gen(s_conj);
      l_perm.push_back(permutation_totality(cell.EXT, s_conj));
    }
    Telt id_perm(2 * cell.EXT.rows());
    cell.perm_grp = Tgroup(l_perm, id_perm);
    cell.EXTtot = totality(cell.EXT);
    cell.l_basis = basis_rows(cell.EXTtot);
    cell.basis_inv = Inverse(UniversalMatrixConversion<T, Tint>(select_rows(cell.EXTtot, cell.l_basis)));
    res.l_top.emplace_back(std::move(cell));
    return {static_cast<int>(res.l_top.size()) - 1, id};
  };
  f_insert(0, id);
  size_t pos = 0;
  while (pos < res.l_top.size()) {
    int iPerf = res.l_top[pos].iPerf;
    MyMatrix<Tint> g = res.l_top[pos].g;
    std::vector<TopCell::Adj> l_adj;
    for (auto &adj : ll_adj[iPerf]) {
      InsertResult ins = f_insert(adj.jPerf, adj.eMatB * g);
      l_adj.push_back({ins.iCell, ins.eEquiv, adj.eSet});
    }
    res.l_top[pos].l_adj = l_adj;
    pos += 1;
  }
#ifdef DEBUG_SP4_COMPLEX
  os << "SP4: |l_top|=" << res.l_top.size() << " |l_gens_G|=" << res.l_gens_G.size() << "\n";
#endif
  return res;
}

// A cell given as the subset eSet of the rows of the top cell iTopCell,
// mapped by eMat: its vectors are l_top[iTopCell].EXT{eSet} * eMat.
struct Triple {
  int iTopCell;
  std::vector<int> eSet;
  MyMatrix<Tint> eMat;
};

/*
  The equivalence of two triples (SP4Z_TestEquivalenceTriple): an element
  M of Sp_4(Z) with EXT1 * M = EXT2 (up to the signs) obtained from the
  stabilizer of their top cell.
 */
inline std::optional<MyMatrix<Tint>> test_equivalence_triple(std::vector<TopCell> const &l_top, Triple const &t1,
                                                             Triple const &t2) {
  if (t1.iTopCell != t2.iTopCell || t1.eSet.size() != t2.eSet.size()) {
    return {};
  }
  TopCell const &cell = l_top[t1.iTopCell];
  int n_ext = cell.EXT.rows();
  Face f1(2 * n_ext), f2(2 * n_ext);
  for (auto &i : t1.eSet) {
    f1[i] = 1;
    f1[n_ext + i] = 1;
  }
  for (auto &i : t2.eSet) {
    f2[i] = 1;
    f2[n_ext + i] = 1;
  }
  std::optional<Telt> opt = cell.perm_grp.RepresentativeAction_OnSets(f1, f2);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> M = Inverse(t1.eMat) * cell.matrix_of(*opt) * t2.eMat;
  return M;
}

/*
  An orbit of cells: its triples (the ways the cell appears in the top
  cells, up to the equivalence), the generators of its stabilizer, its
  vectors, its boundary (the faces, with the elements mapping the
  representative of their orbit to them and the signs) and its
  orientation character: the stabilizer acts on the vectors (up to sign)
  by the permutations perm_gens, with the signs sign_gens; the group
  perm_or on the vectors and two extra points, generated by the
  permutations extended by the transposition of the extra points for the
  sign -1, gives the sign of any element of the stabilizer.
 */
struct CellOrbit {
  std::vector<Triple> l_triple;
  std::vector<MyMatrix<Tint>> l_stab_gens;
  MyMatrix<Tint> EXT;
  std::vector<BoundEntry<Tint>> l_bound;
  std::vector<std::vector<int>> l_set_m1;
  Tgroup perm_or;
  Tgroup perm_grp;
};

// The sign of an element g of the stabilizer of the cell.
inline int orbit_orientation(CellOrbit const &orbit, MyMatrix<Tint> const &g) {
  std::optional<std::vector<Tidx>> opt = permutation_antipodal(orbit.EXT, g);
  if (!opt) {
    std::cerr << "SP4: orbit_orientation, the element does not stabilize the cell\n";
    throw TerminalException{1};
  }
  std::vector<Tidx> l_pos = *opt;
  int n = orbit.EXT.rows();
  l_pos.push_back(n);
  l_pos.push_back(n + 1);
  Telt elt(l_pos);
  if (orbit.perm_or.isin(elt)) {
    return 1;
  }
  return -1;
}

// The sign changes of the faces (UntangleAllSigns): each codimension two
// face occurs in two faces, which gives a relation between their signs.
struct SignOccurrence {
  int sign;
  int idx;
};

inline std::vector<int> untangle_all_signs(std::vector<std::vector<SignOccurrence>> const &l_case, int const &n_var) {
  std::vector<int> v_sign(n_var, 0);
  std::vector<int> assigned(n_var, 0);
  v_sign[0] = 1;
  assigned[0] = 1;
  std::vector<int> used(l_case.size(), 0);
  while (true) {
    int n_unused = 0;
    for (size_t i = 0; i < l_case.size(); i++) {
      std::vector<SignOccurrence> const &c = l_case[i];
      if (c.size() != 2) {
        std::cerr << "SP4: untangle_all_signs, a face of codimension two occurs " << c.size() << " times\n";
        throw TerminalException{1};
      }
      int sum = assigned[c[0].idx] + assigned[c[1].idx];
      if (used[i] == 0 && sum > 0) {
        n_unused += 1;
        if (sum == 2) {
          if (c[0].sign * v_sign[c[0].idx] + c[1].sign * v_sign[c[1].idx] != 0) {
            std::cerr << "SP4: untangle_all_signs, inconsistent signs\n";
            throw TerminalException{1};
          }
        } else {
          int ref = (assigned[c[0].idx] == 1) ? c[0].idx : c[1].idx;
          int oth = (assigned[c[0].idx] == 1) ? c[1].idx : c[0].idx;
          v_sign[oth] = -c[0].sign * c[1].sign * v_sign[ref];
          assigned[oth] = 1;
        }
        used[i] = 1;
      }
    }
    if (n_unused == 0) {
      break;
    }
  }
  for (auto &s : v_sign) {
    if (s != 1 && s != -1) {
      std::cerr << "SP4: untangle_all_signs, a sign is not determined\n";
      throw TerminalException{1};
    }
  }
  return v_sign;
}

/*
  The cell complex of Sp_4(Z) (SP4Z_GetCellComplex): the orbits of cells
  of each level with their stabilizers, boundaries and orientations, and
  for each top cell its faces of each level (the triples with the orbit
  and the element mapping the representative to them).
 */
struct FaceIndex {
  std::vector<int> eSet;
  int iOrbit;
  MyMatrix<Tint> eEquiv;
};

struct Sp4CellComplex {
  std::vector<std::vector<MyMatrix<Tint>>> l_conf;
  std::vector<TopCell> l_top;
  std::vector<MyMatrix<Tint>> l_gens_G;
  // levels[index][iOrbit]
  std::vector<std::vector<CellOrbit>> levels;
  // face_index[index][iTopCell]
  std::vector<std::vector<std::vector<FaceIndex>>> face_index;
};

inline Sp4CellComplex compute_sp4_cell_complex(std::ostream &os) {
#ifdef TIMINGS_SP4_COMPLEX
  MicrosecondTime time;
#endif
  Sp4CellComplex cx;
  cx.l_conf = get_representative_configurations();
  int n_levels = cx.l_conf.size();
  std::vector<std::vector<FlipInfo>> ll_adj = get_flipping_info(cx.l_conf, os);
  TopCellResult tcr = get_top_cells(cx.l_conf, ll_adj, os);
  cx.l_top = std::move(tcr.l_top);
  cx.l_gens_G = std::move(tcr.l_gens_G);
#ifdef TIMINGS_SP4_COMPLEX
  os << "|SP4: top cells|=" << time << "\n";
#endif
  std::vector<TopCell> const &l_top = cx.l_top;
  int n_top = l_top.size();
  // The saturation of a triple: the triples of the same cell in the top
  // cells containing it, and the stabilizer.
  auto f_saturation = [&](Triple const &t0) -> CellOrbit {
    CellOrbit orbit;
    orbit.EXT = select_rows(l_top[t0.iTopCell].EXT, t0.eSet) * t0.eMat;
    std::set<MyMatrix<Tint>> set_gens;
    auto f_insert_gen = [&](MyMatrix<Tint> const &M) -> void {
#ifdef SANITY_CHECK_SP4_COMPLEX
      if (!permutation_antipodal(orbit.EXT, M)) {
        std::cerr << "SP4: a generator of the stabilizer does not stabilize the cell\n";
        throw TerminalException{1};
      }
#endif
      if (set_gens.insert(M).second) {
        orbit.l_stab_gens.push_back(M);
      }
    };
    TopCell const &cell0 = l_top[t0.iTopCell];
    int n_ext = cell0.EXT.rows();
    Face f(2 * n_ext);
    for (auto &i : t0.eSet) {
      f[i] = 1;
      f[n_ext + i] = 1;
    }
    MyMatrix<Tint> eMat_inv = Inverse(t0.eMat);
    for (auto &elt : cell0.perm_grp.Stabilizer_OnSets(f).GeneratorsOfGroup()) {
      f_insert_gen(eMat_inv * cell0.matrix_of(elt) * t0.eMat);
    }
    std::vector<int> status;
    auto f_insert_triple = [&](Triple const &t) -> void {
      for (auto &t_old : orbit.l_triple) {
        std::optional<MyMatrix<Tint>> opt = test_equivalence_triple(l_top, t, t_old);
        if (opt) {
          f_insert_gen(*opt);
          return;
        }
      }
      orbit.l_triple.push_back(t);
      status.push_back(0);
    };
    f_insert_triple(t0);
    size_t pos = 0;
    while (pos < orbit.l_triple.size()) {
      Triple t = orbit.l_triple[pos];
      for (auto &adj : l_top[t.iTopCell].l_adj) {
        if (!std::includes(adj.eSet.begin(), adj.eSet.end(), t.eSet.begin(), t.eSet.end())) {
          continue;
        }
        MyMatrix<Tint> fMat = adj.eEquiv * t.eMat;
        MyMatrix<Tint> EXTimg = l_top[adj.iCell].EXT * fMat;
        std::vector<int> fSet;
        for (int i = 0; i < orbit.EXT.rows(); i++) {
          std::optional<int> opt = position_antipodal(EXTimg, GetMatrixRow(orbit.EXT, i));
          if (!opt) {
            std::cerr << "SP4: the cell is not in the adjacent top cell\n";
            throw TerminalException{1};
          }
          fSet.push_back(*opt);
        }
        std::sort(fSet.begin(), fSet.end());
        f_insert_triple({adj.iCell, fSet, fMat});
      }
      pos += 1;
    }
    return orbit;
  };
  cx.levels = std::vector<std::vector<CellOrbit>>(n_levels);
  cx.face_index = std::vector<std::vector<std::vector<FaceIndex>>>(n_levels);
  for (int index = 0; index < n_levels; index++) {
    std::vector<CellOrbit> &l_orbit = cx.levels[index];
    auto f_insert_orbit = [&](Triple const &t) -> std::pair<int, MyMatrix<Tint>> {
      for (size_t iOrbit = 0; iOrbit < l_orbit.size(); iOrbit++) {
        for (auto &t_orb : l_orbit[iOrbit].l_triple) {
          std::optional<MyMatrix<Tint>> opt = test_equivalence_triple(l_top, t_orb, t);
          if (opt) {
            return {static_cast<int>(iOrbit), *opt};
          }
        }
      }
      l_orbit.push_back(f_saturation(t));
      return {static_cast<int>(l_orbit.size()) - 1, IdentityMat<Tint>(4)};
    };
    cx.face_index[index] = std::vector<std::vector<FaceIndex>>(n_top);
    for (int iTopCell = 0; iTopCell < n_top; iTopCell++) {
      MyMatrix<Tint> const &ePerf = cx.l_conf[0][l_top[iTopCell].iPerf];
      for (auto &eFace : cx.l_conf[index]) {
        for (auto &eSet : get_all_embeddings(eFace, ePerf)) {
          Triple t{iTopCell, eSet, IdentityMat<Tint>(4)};
          std::pair<int, MyMatrix<Tint>> pair = f_insert_orbit(t);
          cx.face_index[index][iTopCell].push_back({eSet, pair.first, pair.second});
        }
      }
    }
#ifdef DEBUG_SP4_COMPLEX
    os << "SP4: index=" << index << " |l_orbit|=" << l_orbit.size() << "\n";
#endif
  }
#ifdef TIMINGS_SP4_COMPLEX
  os << "|SP4: orbits of cells|=" << time << "\n";
#endif
  /*
    The boundaries, from the single vectors (the last index) upwards, with
    the orientation characters.
   */
  for (int index = n_levels - 1; index >= 0; index--) {
    for (auto &orbit : cx.levels[index]) {
      Triple const &t0 = orbit.l_triple[0];
      std::vector<std::vector<Tidx>> l_perm_gens;
      for (auto &g : orbit.l_stab_gens) {
        l_perm_gens.push_back(*permutation_antipodal(orbit.EXT, g));
      }
      int n_ext = orbit.EXT.rows();
      std::vector<int> l_sign_gens(orbit.l_stab_gens.size(), 1);
      if (index < n_levels - 1) {
        // The faces: the cells of the next index contained in the cell, in
        // its top cell.
        std::vector<int> l_iface;
        std::vector<MyMatrix<Tint>> l_elt;
        for (auto &fi : cx.face_index[index + 1][t0.iTopCell]) {
          if (std::includes(t0.eSet.begin(), t0.eSet.end(), fi.eSet.begin(), fi.eSet.end())) {
            l_iface.push_back(fi.iOrbit);
            l_elt.push_back(fi.eEquiv);
            std::vector<int> set_m1;
            for (auto &x : fi.eSet) {
              set_m1.push_back(std::lower_bound(t0.eSet.begin(), t0.eSet.end(), x) - t0.eSet.begin());
            }
            orbit.l_set_m1.push_back(set_m1);
          }
        }
        int n_face = l_iface.size();
        std::vector<int> l_sign;
        if (index == n_levels - 2) {
          if (n_face != 2) {
            std::cerr << "SP4: a cell of two vectors should have two faces\n";
            throw TerminalException{1};
          }
          l_sign = {1, -1};
        } else {
          // The faces of codimension two, each occurring in two faces.
          std::vector<std::vector<int>> l_set_m2;
          std::vector<MyMatrix<Tint>> l_elt_m2;
          std::vector<std::vector<SignOccurrence>> l_case;
          for (int i_face = 0; i_face < n_face; i_face++) {
            CellOrbit const &orbit_m1 = cx.levels[index + 1][l_iface[i_face]];
            for (auto &eb : orbit_m1.l_bound) {
              MyMatrix<Tint> add_elt = eb.M * l_elt[i_face];
              CellOrbit const &orbit_m2 = cx.levels[index + 2][eb.iOrb];
              MyMatrix<Tint> EXTimg = orbit_m2.EXT * add_elt;
              std::vector<int> set_img;
              for (int i = 0; i < EXTimg.rows(); i++) {
                set_img.push_back(*position_antipodal(orbit.EXT, GetMatrixRow(EXTimg, i)));
              }
              std::sort(set_img.begin(), set_img.end());
              auto iter = std::find(l_set_m2.begin(), l_set_m2.end(), set_img);
              if (iter == l_set_m2.end()) {
                l_set_m2.push_back(set_img);
                l_elt_m2.push_back(add_elt);
                l_case.push_back({{eb.sign, i_face}});
              } else {
                size_t pos = iter - l_set_m2.begin();
                MyMatrix<Tint> quot = l_elt_m2[pos] * Inverse(add_elt);
                int mul_sign = orbit_orientation(orbit_m2, quot);
                l_case[pos].push_back({eb.sign * mul_sign, i_face});
              }
            }
          }
          l_sign = untangle_all_signs(l_case, n_face);
        }
        for (int i_face = 0; i_face < n_face; i_face++) {
          orbit.l_bound.push_back({l_iface[i_face], l_sign[i_face], l_elt[i_face]});
        }
        // The signs of the generators of the stabilizer.
        for (size_t i_gen = 0; i_gen < orbit.l_stab_gens.size(); i_gen++) {
          MyMatrix<Tint> const &g = orbit.l_stab_gens[i_gen];
          if (index == n_levels - 2) {
            // The sign of the permutation of the two vectors.
            l_sign_gens[i_gen] = (l_perm_gens[i_gen][0] == 0) ? 1 : -1;
            continue;
          }
          std::vector<int> const &set1 = orbit.l_set_m1[0];
          MyMatrix<Tint> EXT1img = select_rows(orbit.EXT, set1) * g;
          std::vector<int> set1_img;
          for (int i = 0; i < EXT1img.rows(); i++) {
            set1_img.push_back(*position_antipodal(orbit.EXT, GetMatrixRow(EXT1img, i)));
          }
          std::sort(set1_img.begin(), set1_img.end());
          auto iter = std::find(orbit.l_set_m1.begin(), orbit.l_set_m1.end(), set1_img);
          if (iter == orbit.l_set_m1.end()) {
            std::cerr << "SP4: the image of a face is not a face\n";
            throw TerminalException{1};
          }
          size_t pos = iter - orbit.l_set_m1.begin();
          MyMatrix<Tint> quot = orbit.l_bound[0].M * g * Inverse(orbit.l_bound[pos].M);
          CellOrbit const &orbit_m1 = cx.levels[index + 1][orbit.l_bound[0].iOrb];
          int sign = orbit_orientation(orbit_m1, quot);
          l_sign_gens[i_gen] = orbit.l_bound[0].sign * orbit.l_bound[pos].sign * sign;
        }
      }
      // The orientation group: the permutations extended by the sign.
      std::vector<Telt> l_ext_perm, l_perm;
      for (size_t i_gen = 0; i_gen < orbit.l_stab_gens.size(); i_gen++) {
        std::vector<Tidx> v = l_perm_gens[i_gen];
        l_perm.push_back(Telt(v));
        if (l_sign_gens[i_gen] == 1) {
          v.push_back(n_ext);
          v.push_back(n_ext + 1);
        } else {
          v.push_back(n_ext + 1);
          v.push_back(n_ext);
        }
        l_ext_perm.push_back(Telt(v));
      }
      orbit.perm_or = Tgroup(l_ext_perm, Telt(n_ext + 2));
      orbit.perm_grp = Tgroup(l_perm, Telt(n_ext));
    }
  }
#ifdef TIMINGS_SP4_COMPLEX
  os << "|SP4: boundaries|=" << time << "\n";
#endif
  return cx;
}

/*
  The complex of Sp_4(Z) as step B of src_perfect sees it (see
  PerfectComplexView in hecke_homology.h).
 */
struct Sp4ComplexView {
  Sp4CellComplex const &cx;
  int dim() const { return 4; }
  int n_levels() const { return cx.levels.size(); }
  int n_orbit(int const &index) const { return cx.levels[index].size(); }
  std::vector<MyMatrix<Tint>> const &stabilizer_generators(int const &index, int const &iOrb) const {
    return cx.levels[index][iOrb].l_stab_gens;
  }
  int orientation(int const &index, int const &iOrb, MyMatrix<Tint> const &g) const {
    return orbit_orientation(cx.levels[index][iOrb], g);
  }
  bool is_well_rounded(int const &index, int const &iOrb) const {
    return RankMat(cx.levels[index][iOrb].EXT) == 4;
  }
  // The stabilizer acts on the vectors (up to sign) with kernel {+-1}.
  int64_t stabilizer_order(int const &index, int const &iOrb) const {
    return 2 * cx.levels[index][iOrb].perm_grp.size().get_si();
  }
  std::vector<BoundEntry<Tint>> const &boundary(int const &index, int const &iOrb) const {
    return cx.levels[index][iOrb].l_bound;
  }
  std::vector<MyMatrix<Tint>> ambient_generators([[maybe_unused]] std::ostream &os) const { return cx.l_gens_G; }
  bool in_ambient_group(MyMatrix<T> const &g) const { return is_in_sp4(g); }
  bool has_cells_at_infinity() const { return true; }
};

}  // namespace sp4

// clang-format off
#endif  // SRC_PROJECTIVE_SYSTEMS_SP4_COMPLEX_H_
// clang-format on
