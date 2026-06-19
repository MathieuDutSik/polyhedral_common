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

// Enumerate short vectors of eG up to norm M + 100 * tol where M is the
// arithmetic minimum of eG. The 100 * tol safety margin guarantees that for
// every vector v of norm X kept in the family, the family also contains every
// vector whose norm lies in [X - tol, X + tol] (the property may only fail at
// the boundary, in which case the bound is enlarged until it holds).
template <typename T, typename Tint>
MyMatrix<Tint> GetApproxShortVectors_tol(MyMatrix<T> const &eG, T const &tol,
                                         std::ostream &os) {
  int dim = eG.rows();
  Tshortest<T, Tint> rec_shv = T_ShortestVector<T, Tint>(eG, os);
  T M = rec_shv.min;
#ifdef DEBUG_APPROX_AUTO_EQUIV
  os << "AAE: arithmetic minimum M=" << M << " tol=" << tol << "\n";
#endif
  T bound = M + T(100) * tol;
  CVPSolver<T, Tint> solver(eG, os);
  std::vector<MyVector<Tint>> ListVect;
  while (true) {
    ListVect = solver.at_most_norm_vectors(bound);
    // Check the [X - tol, X + tol] inclusion property for boundary vectors.
    T max_norm(0);
    bool first = true;
    for (auto const &V : ListVect) {
      MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
      T norm = V_T.dot(eG * V_T);
      if (first || norm > max_norm) {
        max_norm = norm;
        first = false;
      }
    }
    if (first) {
      break;
    }
    if (max_norm + tol <= bound) {
      break;
    }
    bound = max_norm + T(2) * tol;
#ifdef DEBUG_APPROX_AUTO_EQUIV
    os << "AAE: enlarging short vector bound to " << bound << "\n";
#endif
  }
  std::vector<MyVector<Tint>> Filtered;
  Filtered.reserve(2 * ListVect.size());
  for (auto &V : ListVect) {
    if (IsZeroVector(V)) {
      continue;
    }
    Filtered.push_back(V);
    MyVector<Tint> Vneg = -V;
    Filtered.push_back(Vneg);
  }
#ifdef DEBUG_APPROX_AUTO_EQUIV
  os << "AAE: |short vectors|=" << Filtered.size() << " (centrally symmetric)\n";
#endif
  return MatrixFromVectorFamilyDim(dim, Filtered);
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
template <typename T>
bool IsApproximateAutomorphism(MyMatrix<T> const &eG, MyMatrix<T> const &M,
                               T const &tol) {
  MyMatrix<T> diff = M * eG * M.transpose() - eG;
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

// Compute generators of the approximate automorphism group of eG.
// Each returned generator M satisfies M * eG * M^T approx eG within tol,
// matching the convention used in src_latt/LatticeStabEquiCan.h.
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<T>>
ApproximateAutomorphismGroup(MyMatrix<T> const &eG, T const &tol,
                             std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tgr = GraphListAdj;
  MyMatrix<Tint> SHV = GetApproxShortVectors_tol<T, Tint>(eG, tol, os);
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
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  std::vector<MyMatrix<T>> ListGen;
  for (auto &ePerm : GRP.GeneratorsOfGroup()) {
    std::optional<MyMatrix<T>> opt =
        FindTransformationGeneral<T, Telt>(SHV_T, SHV_T, ePerm);
    if (!opt) {
#ifdef DEBUG_APPROX_AUTO_EQUIV
      os << "AAE: skipping permutation: no linear lift exists\n";
#endif
      continue;
    }
    MyMatrix<T> const &M = *opt;
    if (!IsApproximateAutomorphism(eG, M, tol)) {
#ifdef DEBUG_APPROX_AUTO_EQUIV
      os << "AAE: skipping lift: not an approximate automorphism\n";
#endif
      continue;
    }
    ListGen.push_back(M);
  }
#ifdef DEBUG_APPROX_AUTO_EQUIV
  os << "AAE: returning " << ListGen.size() << " approximate generators\n";
#endif
  return ListGen;
}

// clang-format off
#endif  // SRC_RANKIN_APPROXAUTOEQUIV_H_
// clang-format on
