// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_APPROXAUTOEQUIV_H_
#define SRC_RANKIN_APPROXAUTOEQUIV_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include "PolytopeEquiStabInt.h"
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_APPROX_AUTO_EQUIV
#endif

#ifdef DISABLE_DEBUG_APPROX_AUTO_EQUIV
#undef DEBUG_APPROX_AUTO_EQUIV
#endif

/*
  Approximate automorphism / equivalence for positive definite Gram matrices
  given up to some tolerance "tol".

  The strategy mirrors the exact code in src_latt/LatticeStabEquiCan.h:
   * Enumerate the short vectors with a safety margin around the minimum so
     that vectors whose norm differs by less than tol from an included vector
     are also enumerated.
   * Build the matrix of scalar products and identify two scalar values that
     differ by less than tol.
   * Compute the automorphism group of the resulting weighted graph.
   * Lift each permutation generator to a linear map on the short vectors and
     verify that the map is an approximate automorphism of the Gram matrix.
 */



/*
  Compute from the diagonal the entries.
 */
template <typename T, typename Tint>
std::pair<MyMatrix<Tint>, T> GetApproxShortVectors_tol(MyMatrix<T> const &eG,
                                                       std::optional<T> const& max_norm_opt,
                                                       T const &tol,
                                                       std::ostream &os) {
  int dim = eG.rows();
  // Maximum of the LLL reduced form.
  T max_norm = [&]() -> T {
    if (max_norm_opt) {
      return *max_norm_opt + 2 * tol;
    } else {
      return GetMaxNorm<T,Tint>(eG, os);
    }
  }();
#ifdef DEBUG_APPROX_AUTO_EQUIV
  os << "AAE: arithmetic minimum max_norm=" << max_norm << " tol=" << tol << "\n";
#endif

  CVPSolver<T, Tint> solver(eG, os);
  // We enumerate every lattice vector v (in the central-symmetry half) of norm
  // at most max_norm, not only those that sit on the exact shell, so that the
  // "+ 2 * tol" buffer above (used when max_norm comes from another lattice)
  // actually captures the extra vectors instead of giving an empty set.
  std::vector<MyVector<Tint>> ListVect = solver.at_most_norm_vectors(max_norm);
#ifdef DEBUG_APPROX_AUTO_EQUIV
  os << "AAE: |short vectors|=" << ListVect.size() << " (centrally symmetric)\n";
#endif
  MyMatrix<Tint> M = MatrixFromVectorFamilyDim(dim, ListVect);
  return {M, max_norm};
}






// Build the symmetric weight matrix of scalar products v_i^T eG v_j on the
// rows of SHV, identifying values that differ by less than tol.
template <typename T, typename Tint>
WeightMatrix<true, T, uint32_t> GetWeightMatrix_tol(MyMatrix<T> const &eG,
                                                   MyMatrix<Tint> const &SHV,
                                                   T const &tol,
                                                   std::ostream &os) {
  size_t n_row = SHV.rows();
  std::vector<MyVector<T>> ListV(n_row);
  std::vector<MyVector<T>> ListGV(n_row);
  for (size_t i = 0; i < n_row; i++) {
    MyVector<Tint> v_int = GetMatrixRow(SHV, static_cast<int>(i));
    MyVector<T> v_T = UniversalVectorConversion<T, Tint>(v_int);
    ListV[i] = v_T;
    ListGV[i] = eG * v_T;
  }
  std::vector<T> ListVal;
  auto get_idx = [&](T const &val) -> uint32_t {
    for (size_t u = 0; u < ListVal.size(); u++) {
      T diff = ListVal[u] - val;
      if (T_abs(diff) < tol) {
        return static_cast<uint32_t>(u);
      }
    }
    uint32_t pos = static_cast<uint32_t>(ListVal.size());
    ListVal.push_back(val);
    return pos;
  };
  WeightMatrix<true, T, uint32_t> WMat(n_row, os);
  for (size_t i_row = 0; i_row < n_row; i_row++) {
    for (size_t j_row = 0; j_row <= i_row; j_row++) {
      T scal = ListGV[i_row].dot(ListV[j_row]);
      uint32_t pos = get_idx(scal);
      WMat.intDirectAssign(i_row, j_row, pos);
    }
  }
  WMat.SetWeight(ListVal);
  return WMat;
}

// Check if M is an approximate automorphism of eG using the convention
// M eG M^T = eG up to tolerance (matches LatticeStabEquiCan.h).
template <typename T, typename Tint>
bool IsApproximateAutomorphism(MyMatrix<T> const &eG, MyMatrix<Tint> const &M,
                               T const &tol) {
  MyMatrix<T> M_T = UniversalMatrixConversion<T,Tint>(M);
  MyMatrix<T> diff = M_T * eG * M_T.transpose() - eG;
  int n = diff.rows();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (T_abs(diff(i, j)) >= tol) {
        return false;
      }
    }
  }
  return true;
}

// Reorder the weights of WMat2 so that they align with the weights of WMatRef
// within tolerance tol. Returns false if the two weight lists cannot be put in
// approximate bijection (different sizes or some weight has no match within
// tol). The semantics mirror the exact RenormalizeWeightMatrix in
// WeightMatrix.h, but pair weights using T_abs(diff) < tol instead of equality.
template <typename T>
bool RenormalizeWeightMatrix_tol(
    WeightMatrix<true, T, uint32_t> const &WMatRef,
    WeightMatrix<true, T, uint32_t> &WMat2, T const &tol) {
  if (WMatRef.rows() != WMat2.rows()) {
    return false;
  }
  std::vector<T> const &ListWeightRef = WMatRef.GetWeight();
  std::vector<T> const &ListWeight2 = WMat2.GetWeight();
  size_t nbEnt = ListWeightRef.size();
  if (nbEnt != ListWeight2.size()) {
    return false;
  }
  std::vector<uint32_t> gListRev(nbEnt);
  std::vector<bool> matched(nbEnt, false);
  for (size_t i = 0; i < nbEnt; i++) {
    bool found = false;
    for (size_t j = 0; j < nbEnt; j++) {
      if (matched[j]) {
        continue;
      }
      T diff = ListWeight2[j] - ListWeightRef[i];
      if (T_abs(diff) < tol) {
        gListRev[j] = static_cast<uint32_t>(i);
        matched[j] = true;
        found = true;
        break;
      }
    }
    if (!found) {
      return false;
    }
  }
  WMat2.ReorderingOfWeights(gListRev);
  return true;
}

// Compute generators of the approximate automorphism group of eG.
// Each returned generator M satisfies M * eG * M^T approx eG within tol,
// matching the convention used in src_latt/LatticeStabEquiCan.h.
template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<MyMatrix<Tint>>>
ApproximateAutomorphismGroup(MyMatrix<T> const &eG, T const &tol,
                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tgr = GraphListAdj;
  std::optional<T> max_norm_opt;
  std::pair<MyMatrix<Tint>,T> pair = GetApproxShortVectors_tol<T, Tint>(eG, max_norm_opt, tol, os);
  MyMatrix<Tint> const& SHV = pair.first;
  int n_row = SHV.rows();
  if (n_row == 0) {
    std::cerr << "AAE: empty short vector set\n";
    throw TerminalException{1};
  }
  WeightMatrix<true, T, uint32_t> WMat =
      GetWeightMatrix_tol<T, Tint>(eG, SHV, tol, os);
  WMat.ReorderingSetWeight();
  Tgroup GRP =
      GetStabilizerWeightMatrix<T, Tgr, Tgroup, uint32_t>(WMat, os);
#ifdef DEBUG_APPROX_AUTO_EQUIV
  os << "AAE: |GRP|=" << GRP.size() << " on " << n_row << " short vectors\n";
#endif
  std::vector<MyMatrix<Tint>> ListGen;
  for (auto &ePerm : GRP.GeneratorsOfGroup()) {
    std::optional<MyMatrix<Tint>> opt =
        FindTransformation<Tint, Telt>(SHV, SHV, ePerm);
    if (!opt) {
#ifdef DEBUG_APPROX_AUTO_EQUIV
      os << "AAE: skipping permutation: no linear lift exists\n";
#endif
      return {};
    }
    MyMatrix<Tint> const &M = *opt;
    if (!IsApproximateAutomorphism(eG, M, tol)) {
#ifdef DEBUG_APPROX_AUTO_EQUIV
      os << "AAE: The matrix is not an approximate automorphism\n";
#endif
      return {};
    }
    ListGen.push_back(M);
  }
#ifdef DEBUG_APPROX_AUTO_EQUIV
  os << "AAE: returning " << ListGen.size() << " approximate generators\n";
#endif
  return ListGen;
}

// Check that the integer matrix M is an approximate equivalence from eG1 to
// eG2, i.e. that M * eG1 * M^T approx eG2 within tol on every entry (same
// convention as in src_latt/LatticeStabEquiCan.h).
template <typename T, typename Tint>
bool IsApproximateEquivalence(MyMatrix<T> const &eG1, MyMatrix<T> const &eG2,
                              MyMatrix<Tint> const &M, T const &tol) {
  MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
  MyMatrix<T> diff = M_T * eG1 * M_T.transpose() - eG2;
  int n = diff.rows();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (T_abs(diff(i, j)) >= tol) {
        return false;
      }
    }
  }
  return true;
}

// Approximate arithmetic equivalence between two positive definite Gram
// matrices eG1 and eG2 up to tolerance tol. On success returns an integer
// matrix M with M * eG1 * M^T approx eG2 within tol (same convention as
// ArithmeticEquivalence in src_latt/LatticeStabEquiCan.h). Returns std::nullopt
// if no such equivalence is detected. There is no "honest failure" report (in
// contrast to the automorphism path): inability to find an equivalence is
// reported as "not equivalent".
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
ApproximateEquivalence(MyMatrix<T> const &eG1, MyMatrix<T> const &eG2,
                       T const &tol, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  std::optional<T> max_norm_opt1;
  std::pair<MyMatrix<Tint>, T> pair1 =
      GetApproxShortVectors_tol<T, Tint>(eG1, max_norm_opt1, tol, os);
  MyMatrix<Tint> const &SHV1 = pair1.first;
  T const &max_norm = pair1.second;
  std::optional<T> max_norm_opt2(max_norm);
  std::pair<MyMatrix<Tint>, T> pair2 =
      GetApproxShortVectors_tol<T, Tint>(eG2, max_norm_opt2, tol, os);
  MyMatrix<Tint> const &SHV2 = pair2.first;
#ifdef DEBUG_APPROX_AUTO_EQUIV
  os << "AAE: |SHV1|=" << SHV1.rows() << " |SHV2|=" << SHV2.rows() << "\n";
#endif
  if (SHV1.rows() != SHV2.rows()) {
    return {};
  }
  WeightMatrix<true, T, uint32_t> WMat1 =
      GetWeightMatrix_tol<T, Tint>(eG1, SHV1, tol, os);
  WMat1.ReorderingSetWeight();
  WeightMatrix<true, T, uint32_t> WMat2 =
      GetWeightMatrix_tol<T, Tint>(eG2, SHV2, tol, os);
  WMat2.ReorderingSetWeight();
  bool test = RenormalizeWeightMatrix_tol<T>(WMat1, WMat2, tol);
  if (!test) {
#ifdef DEBUG_APPROX_AUTO_EQUIV
    os << "AAE: tolerance-aware renormalization failed\n";
#endif
    return {};
  }
  std::optional<Telt> opt_perm =
      TestEquivalenceWeightMatrix_norenorm_perm<T, Telt, uint32_t>(WMat2, WMat1,
                                                                   os);
  if (!opt_perm) {
#ifdef DEBUG_APPROX_AUTO_EQUIV
    os << "AAE: no graph isomorphism between the weighted matrices\n";
#endif
    return {};
  }
  std::optional<MyMatrix<Tint>> opt_mat =
      FindTransformation<Tint, Telt>(SHV2, SHV1, *opt_perm);
  if (!opt_mat) {
#ifdef DEBUG_APPROX_AUTO_EQUIV
    os << "AAE: no linear lift of the permutation\n";
#endif
    return {};
  }
  MyMatrix<Tint> const &M = *opt_mat;
  if (!IsApproximateEquivalence(eG1, eG2, M, tol)) {
#ifdef DEBUG_APPROX_AUTO_EQUIV
    os << "AAE: lift is not an approximate equivalence\n";
#endif
    return {};
  }
  return M;
}

// clang-format off
#endif  // SRC_RANKIN_APPROXAUTOEQUIV_H_
// clang-format on
