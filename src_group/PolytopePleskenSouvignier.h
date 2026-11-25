// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_POLYTOPEPLESKENSOUVIGNIER_H_
#define SRC_GROUP_POLYTOPEPLESKENSOUVIGNIER_H_

// clang-format off
#include "PleskenSouvignier.h"
#include "ClassicLLL.h"
#include "MAT_MatrixInt.h"
#include "WeightMatrix.h"
#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYTOPE_PLESKEN_SOUVIGNIER
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYTOPE_PLESKEN_SOUVIGNIER
#endif

#ifdef TIMINGS
#define TIMINGS_POLYTOPE_PLESKEN_SOUVIGNIER
#endif

/*
  The linear automorphism group and the canonical form of a polytope, a
  cone or any finite family of vectors EXT, by the polytope mode of the
  Plesken-Souvignier backtrack (see PleskenSouvignier.h). The group is the
  one of the g in GL_n(Q) with EXT g = P EXT, P a permutation matrix, acting
  on the rows: the group of LinPolytope_Automorphism, computed without the
  weight matrix of the scalar products of the rows, whose graph is quadratic
  in their number (4.9 GB for the 4577 rays of a ridge of the local cone of
  E6).

  The reduction to the integral setting of the backtrack:
  * EXT is column reduced to its rank n and scaled to integers by a single
    multiplier, which commutes with the linear maps.
  * Every automorphism preserves L, the Z-span of the rows, so it is
    integral in a Z-basis of L, in which the rows are expressed.
  * The form is S^{-1}, S = V^T V, scaled to an integral form, as
    GetQmatrix does: from V g = P V comes g^T S g = S, that is
    g S^{-1} g^T = S^{-1}. The basis of L is then LLL reduced for it, as
    the lattice wrappers of src_latt/LatticePleskenSouvignier.h do.
  The group order is the product of the orbit lengths along the
  stabilizer chain of the searched basis, the pointwise stabilizer of a
  basis being trivial.
 */
template <typename Tint> struct PolytopePleskenSouvignierData {
  // The rows, in the LLL reduced basis of L.
  MyMatrix<Tint> V;
  // The invariant form in that basis, positive definite.
  MyMatrix<Tint> M;
};

template <typename T>
PolytopePleskenSouvignierData<typename underlying_ring<
    typename overlying_field<T>::field_type>::ring_type>
PolytopePleskenSouvignierPrepare(MyMatrix<T> const &EXT, std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  using Tint = typename underlying_ring<Tfield>::ring_type;
  MyMatrix<Tfield> EXTred =
      UniversalMatrixConversion<Tfield, T>(ColumnReduction(EXT));
  MyMatrix<Tfield> EXTint = UniversalMatrixConversion<Tfield, Tint>(
      RemoveFractionMatrixPlusCoeffRing(EXTred).TheMat);
  MyMatrix<Tfield> Zbasis = GetZbasis(EXTint);
  MyMatrix<Tfield> V1 = EXTint * Inverse(Zbasis);
  MyMatrix<Tfield> S = V1.transpose() * V1;
  MyMatrix<Tfield> Qmat = UniversalMatrixConversion<Tfield, Tint>(
      RemoveFractionMatrixPlusCoeffRing(Inverse(S)).TheMat);
  LLLreduction<Tfield, Tint> rec = LLLreducedBasis<Tfield, Tint>(Qmat, os);
  // A row v in the reduced basis is v P^{-1}, and the form P Q P^T.
  MyMatrix<Tfield> PmatInv =
      Inverse(UniversalMatrixConversion<Tfield, Tint>(rec.Pmat));
  return {UniversalMatrixConversion<Tint, Tfield>(V1 * PmatInv),
          UniversalMatrixConversion<Tint, Tfield>(rec.GramMatRed)};
}

template <typename Tint, typename Tgroup>
StabGeneratorsOrder<typename Tgroup::Telt::Tidx, typename Tgroup::Tint>
PolytopePleskenSouvignierAutomorphism(
    PolytopePleskenSouvignierData<Tint> const &data, std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  MyMatrix<Tint> const &V = data.V;
  int nbRow = V.rows();
  PleskenSouvignierAutomResult<Tint, Tgroup> result =
      PleskenSouvignierAutomorphism<Tint, Tgroup>(
          {data.M}, V, PSAutomorphismOptions<Tint>{.polytope = true}, os);
  std::unordered_map<MyVector<Tint>, Tidx> MapRow;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MapRow[GetMatrixRow(V, iRow)] = static_cast<Tidx>(iRow);
  }
  std::vector<std::vector<Tidx>> ListGen;
  for (auto &eGen : result.ListGen) {
    MyMatrix<Tint> Vimg = V * eGen;
    std::vector<Tidx> eList(nbRow);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      auto iter = MapRow.find(GetMatrixRow(Vimg, iRow));
#ifdef SANITY_CHECK_POLYTOPE_PLESKEN_SOUVIGNIER
      if (iter == MapRow.end()) {
        std::cerr << "PPS: a generator maps a row outside the family\n";
        throw TerminalException{1};
      }
#endif
      eList[iRow] = iter->second;
    }
    ListGen.push_back(std::move(eList));
  }
#ifdef DEBUG_POLYTOPE_PLESKEN_SOUVIGNIER
  os << "PPS: nbRow=" << nbRow << " n=" << V.cols()
     << " |ListGen|=" << ListGen.size() << " order=" << result.order
     << " depth=" << result.depth_used << " nodes=" << result.nodes_used
     << "\n";
#endif
  return {std::move(ListGen), result.order};
}

template <typename T, typename Tgroup>
StabGeneratorsOrder<typename Tgroup::Telt::Tidx, typename Tgroup::Tint>
LinPolytope_Automorphism_PleskenSouvignier_GensOrder(MyMatrix<T> const &EXT,
                                                     std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  using Tint = typename underlying_ring<Tfield>::ring_type;
#ifdef TIMINGS_POLYTOPE_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  PolytopePleskenSouvignierData<Tint> data =
      PolytopePleskenSouvignierPrepare(EXT, os);
#ifdef TIMINGS_POLYTOPE_PLESKEN_SOUVIGNIER
  os << "|PPS: preparation|=" << time << "\n";
#endif
  StabGeneratorsOrder<typename Tgroup::Telt::Tidx, typename Tgroup::Tint>
      gens_order = PolytopePleskenSouvignierAutomorphism<Tint, Tgroup>(data, os);
#ifdef TIMINGS_POLYTOPE_PLESKEN_SOUVIGNIER
  os << "|PPS: automorphism|=" << time << "\n";
#endif
  return gens_order;
}

template <typename T, typename Tgroup>
Tgroup LinPolytope_Automorphism_PleskenSouvignier(MyMatrix<T> const &EXT,
                                                  std::ostream &os) {
  return GroupFromStabGeneratorsOrder<Tgroup>(
      LinPolytope_Automorphism_PleskenSouvignier_GensOrder<T, Tgroup>(EXT,
                                                                      os),
      EXT.rows());
}

/*
  The canonical ordering of the rows from the automorphism group GRP (a
  permutation group of the rows): CanonicOrd[i] is the row put at position
  i, and for two families equivalent under GL_n(Q) the reordered families
  are equivalent row by row, so that CanonicalizeOrderedMatrix gives the
  same matrix for both (see LinPolytope_CanonicForm).

  It is a search over the ordered bases (x_1, ..., x_n) made of rows, with
  the sequence of orbit computations of the automorphism group:
  * At a node (x_1, ..., x_k), the rows independent of x_1, ..., x_k are
    split into classes by their key (the norm, then the products with
    x_1, ..., x_k, in that order). Everything here is intrinsic: the form
    is invariant up to a positive scalar, which preserves the order of the
    keys.
  * The class chosen is the one with the fewest orbits of the pointwise
    stabilizer H of x_1, ..., x_k in GRP, then the smallest, then the one
    of smallest key. These are intrinsic as well.
  * The search branches on one row of each H-orbit of the class: an
    element of H maps the subtree of a row to the one of another row of
    its orbit, and the leaves of the two give the same matrix below.
  * At a leaf, the matrix of the coordinates of the rows in the basis
    x_1, ..., x_n, which does not depend on the coordinates of the input,
    has its rows sorted; the canonical form is the smallest such matrix
    over the leaves, and the ordering is the one sorting it.
  The number of leaves is the product over the levels of the number of
  orbits of the chosen class, 1 for a family whose group is transitive
  enough or whose keys separate the rows.
 */
template <typename Tint, typename Tgroup>
std::vector<typename Tgroup::Telt::Tidx>
PolytopePleskenSouvignierCanonicOrdering(
    PolytopePleskenSouvignierData<Tint> const &data, Tgroup const &GRP,
    [[maybe_unused]] std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  using Tfield = typename overlying_field<Tint>::field_type;
  MyMatrix<Tint> const &V = data.V;
  int m = V.rows();
  int n = V.cols();
  MyMatrix<Tint> W = V * data.M;
  auto scal = [&](int i, int j) -> Tint {
    Tint sum(0);
    for (int k = 0; k < n; k++) {
      sum += V(i, k) * W(j, k);
    }
    return sum;
  };
  MyMatrix<Tfield> Vf = UniversalMatrixConversion<Tfield, Tint>(V);
  std::optional<MyMatrix<Tfield>> best;
  std::vector<Tidx> best_ord;
  [[maybe_unused]] size_t n_leaf = 0;
  // The rows sorted lexicographically, then compared with the best leaf.
  auto f_leaf = [&](std::vector<int> const &x) -> void {
    n_leaf++;
    MyMatrix<Tfield> B(n, n);
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        B(i, k) = Vf(x[i], k);
      }
    }
    MyMatrix<Tfield> C = Vf * Inverse(B);
    std::vector<Tidx> ord(m);
    for (int i = 0; i < m; i++) {
      ord[i] = static_cast<Tidx>(i);
    }
    auto cmp_rows = [&](MyMatrix<Tfield> const &A, int iA,
                        MyMatrix<Tfield> const &B2, int iB) -> int {
      for (int k = 0; k < n; k++) {
        if (A(iA, k) < B2(iB, k)) {
          return -1;
        }
        if (B2(iB, k) < A(iA, k)) {
          return 1;
        }
      }
      return 0;
    };
    std::sort(ord.begin(), ord.end(), [&](Tidx a, Tidx b) -> bool {
      return cmp_rows(C, a, C, b) < 0;
    });
    MyMatrix<Tfield> Csort(m, n);
    for (int i = 0; i < m; i++) {
      Csort.row(i) = C.row(ord[i]);
    }
    bool is_better = !best.has_value();
    if (!is_better) {
      for (int i = 0; i < m; i++) {
        int test = cmp_rows(Csort, i, *best, i);
        if (test != 0) {
          is_better = test < 0;
          break;
        }
      }
    }
    if (is_better) {
      best = std::move(Csort);
      best_ord = std::move(ord);
    }
  };
  // The orbit index of every row under the generators of H.
  auto get_orbit_index = [&](Tgroup const &H) -> std::vector<int> {
    std::vector<int> orb_idx(m, -1);
    std::vector<typename Tgroup::Telt> LGen = H.GeneratorsOfGroup();
    int n_orb = 0;
    for (int i = 0; i < m; i++) {
      if (orb_idx[i] != -1) {
        continue;
      }
      std::vector<int> orb{i};
      orb_idx[i] = n_orb;
      size_t pos = 0;
      while (pos < orb.size()) {
        for (auto &eGen : LGen) {
          int img = eGen.at(static_cast<Tidx>(orb[pos]));
          if (orb_idx[img] == -1) {
            orb_idx[img] = n_orb;
            orb.push_back(img);
          }
        }
        pos++;
      }
      n_orb++;
    }
    return orb_idx;
  };
  std::vector<std::vector<Tint>> key0(m);
  for (int i = 0; i < m; i++) {
    key0[i].push_back(scal(i, i));
  }
  std::function<void(std::vector<int> &, std::vector<std::vector<Tint>> const &,
                     Tgroup const &)>
      f_rec = [&](std::vector<int> &x,
                  std::vector<std::vector<Tint>> const &key,
                  Tgroup const &H) -> void {
    int k = x.size();
    if (k == n) {
      f_leaf(x);
      return;
    }
    // The rows independent of x: v z^T != 0 for a z of the kernel.
    std::vector<uint8_t> indep(m, 1);
    if (k > 0) {
      MyMatrix<Tint> Prefix(k, n);
      for (int i = 0; i < k; i++) {
        Prefix.row(i) = V.row(x[i]);
      }
      MyMatrix<Tint> NSP = NullspaceTrMat(Prefix);
#ifdef SANITY_CHECK_POLYTOPE_PLESKEN_SOUVIGNIER
      if (NSP.rows() != n - k || !IsZeroMatrix(MyMatrix<Tint>(
                                     Prefix * NSP.transpose()))) {
        std::cerr << "PPS: wrong kernel of the prefix\n";
        throw TerminalException{1};
      }
#endif
      MyMatrix<Tint> Prod = V * NSP.transpose();
      for (int i = 0; i < m; i++) {
        indep[i] = 0;
        for (int j = 0; j < n - k; j++) {
          if (Prod(i, j) != 0) {
            indep[i] = 1;
            break;
          }
        }
      }
    }
    std::map<std::vector<Tint>, std::vector<int>> classes;
    for (int i = 0; i < m; i++) {
      if (indep[i]) {
        classes[key[i]].push_back(i);
      }
    }
    std::vector<int> orb_idx = get_orbit_index(H);
    std::vector<int> const *chosen = nullptr;
    size_t chosen_norb = 0;
    // The classes are visited by increasing key, so a strict improvement
    // keeps the smallest key among the ties.
    for (auto &kv : classes) {
      std::vector<int> const &cls = kv.second;
      std::vector<int> list_orb;
      for (auto &i : cls) {
        list_orb.push_back(orb_idx[i]);
      }
      std::sort(list_orb.begin(), list_orb.end());
      size_t norb = std::unique(list_orb.begin(), list_orb.end()) -
                    list_orb.begin();
      if (chosen == nullptr || norb < chosen_norb ||
          (norb == chosen_norb && cls.size() < chosen->size())) {
        chosen = &cls;
        chosen_norb = norb;
      }
    }
    std::vector<uint8_t> orb_done(m, 0);
    for (auto &i : *chosen) {
      if (orb_done[orb_idx[i]]) {
        continue;
      }
      orb_done[orb_idx[i]] = 1;
      std::vector<std::vector<Tint>> key_next = key;
      for (int j = 0; j < m; j++) {
        key_next[j].push_back(scal(j, i));
      }
      x.push_back(i);
      f_rec(x, key_next, H.Stabilizer_OnPoints(static_cast<Tidx>(i)));
      x.pop_back();
    }
  };
  std::vector<int> x;
  f_rec(x, key0, GRP);
#ifdef DEBUG_POLYTOPE_PLESKEN_SOUVIGNIER
  os << "PPS: canonic ordering, m=" << m << " n=" << n
     << " |GRP|=" << GRP.size() << " leaves=" << n_leaf << "\n";
#endif
  return best_ord;
}

/*
  The canonical ordering of the rows of EXT (see
  PolytopePleskenSouvignierCanonicOrdering) with the generators of the
  automorphism group, as GetGroupCanonicalizationVector_Kernel gives them.
 */
template <typename T, typename Tgroup>
std::pair<std::vector<typename Tgroup::Telt::Tidx>,
          std::vector<std::vector<typename Tgroup::Telt::Tidx>>>
LinPolytope_CanonicOrdering_PleskenSouvignier(MyMatrix<T> const &EXT,
                                              std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  using Tint = typename underlying_ring<Tfield>::ring_type;
#ifdef TIMINGS_POLYTOPE_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  PolytopePleskenSouvignierData<Tint> data =
      PolytopePleskenSouvignierPrepare(EXT, os);
  StabGeneratorsOrder<typename Tgroup::Telt::Tidx, typename Tgroup::Tint>
      gens_order = PolytopePleskenSouvignierAutomorphism<Tint, Tgroup>(data, os);
#ifdef TIMINGS_POLYTOPE_PLESKEN_SOUVIGNIER
  os << "|PPS: preparation and automorphism|=" << time << "\n";
#endif
  Tgroup GRP = GroupFromStabGeneratorsOrder<Tgroup>(gens_order, EXT.rows());
  std::vector<typename Tgroup::Telt::Tidx> CanonicOrd =
      PolytopePleskenSouvignierCanonicOrdering<Tint, Tgroup>(data, GRP, os);
#ifdef TIMINGS_POLYTOPE_PLESKEN_SOUVIGNIER
  os << "|PPS: canonic ordering|=" << time << "\n";
#endif
  return {std::move(CanonicOrd), std::move(gens_order.ListGen)};
}

// clang-format off
#endif  // SRC_GROUP_POLYTOPEPLESKENSOUVIGNIER_H_
// clang-format on
