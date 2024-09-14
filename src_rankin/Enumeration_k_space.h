// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_ENUMERATION_K_SPACE_H_
#define SRC_RANKIN_ENUMERATION_K_SPACE_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include <utility>
#include <unordered_set>
#include <vector>
// clang-format on

// Returns the Hermite constant at the n-th power.
// That is the mmaximum of min(A)^n / det(A)
// --
// So we would have min(A)^n / det(A) <= H(n)
// which gets us min(A)^n <= H(n) det(A)
template <typename T> T GetUpperBoundHermitePower(int n) {
  if (n == 1) {
    return 1;
  }
  if (n == 2) {
    return T(4) / 3;
  }
  if (n == 3) {
    return 2;
  }
  if (n == 4) {
    return 4;
  }
  if (n == 5) {
    return 8;
  }
  if (n == 6) {
    return T(64) / 3;
  }
  if (n == 7) {
    return 64;
  }
  if (n == 8) {
    return 256;
  }
  int h = n * (n - 1) / 2;
  return (4 / 3) ^ h;
}

/*
  orthogonal projector means self adjoint.
  that is that we have
  <x, p(y)> = < p(x), y>
  We have the scalar product <x,y> = X G Y^t
  and so we get
  <p(x), y> = XP G Y^T
  <x, p(y)> = XG P^T Y
  so we have PG = G P^T
  ---
  See the code of __GetOrthogonalProjector in
  SublatticeEnumeration.g
  Though the single formula needs to be checked.
*/
template <typename T, typename Tint>
MyMatrix<T> GetOrthogonalProjector(MyMatrix<T> const &TheGramMat,
                                   MyMatrix<Tint> const &TheSubBasis) {
  int n = TheGramMat.rows();
  if (n != TheSubBasis.cols()) {
    std::cerr << "The matrix size does not match\n";
    throw TerminalException{1};
  }
  MyMatrix<T> TheSubBasis_T = UniversalMatrixConversion<T, Tint>(TheSubBasis);
  MyMatrix<T> eGram = TheSubBasis_T * TheGramMat * TheSubBasis_T.transpose();
  MyMatrix<T> eGramInv = Inverse(eGram);

  MyMatrix<T> TheProj =
      TheSubBasis_T.transpose() * eGramInv * TheSubBasis_T * TheGramMat;
#ifdef DEBUG_RANKIN
  MyMatrix<T> TheProjSqr = TheProj * TheProj;
  if (TheProj * TheProj != TheProj) {
    std::cerr << "The matrix is not a projector\n";
    throw TerminalException{1};
  }
  if (TheProj * TheGramMat != TheGramMat * TheProj) {
    std::cerr << "The obtained projector is not self adjoint\n";
    throw TerminalException{1};
  }
#endif
  return TheProj;
}

template <typename T, typename Tint>
std::pair<MyMatrix<T>, MyMatrix<Tint>>
GetOrthogonalProjector_dim1(MyMatrix<T> const &TheGramMat,
                            MyVector<Tint> const &eVect, std::ostream &os) {
  os << "RNK: GetOrthogonalProjector_dim1, step 1\n";
  int n = TheGramMat.rows();
  os << "RNK: GetOrthogonalProjector_dim1, step 2 n=" << n << "\n";
  MyMatrix<Tint> eVect_M = MatrixFromVector(eVect);
  os << "RNK: GetOrthogonalProjector_dim1, step 3\n";
  MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect);
  os << "RNK: GetOrthogonalProjector_dim1, step 4 eVect_T="
     << StringVector(eVect_T) << "\n";
  MyVector<T> eVect_T_TheGramMat = eVect_T.transpose() * TheGramMat;
  os << "RNK: GetOrthogonalProjector_dim1, step 5 eVect_T_GramMat="
     << StringVector(eVect_T_TheGramMat) << "\n";
  T rNorm = eVect_T_TheGramMat.dot(eVect_T);
  os << "RNK: GetOrthogonalProjector_dim1, step 6 rNorm=" << rNorm << "\n";
  MyMatrix<Tint> TheCompl = SubspaceCompletionInt(eVect_M, n);
  os << "RNK: GetOrthogonalProjector_dim1, step 7\n";
  MyMatrix<T> TheProj(n - 1, n);
  os << "RNK: GetOrthogonalProjector_dim1, step 8\n";
  for (int i = 0; i < n - 1; i++) {
    MyVector<Tint> fVect = GetMatrixRow(TheCompl, i);
    MyVector<T> fVect_T = UniversalVectorConversion<T, Tint>(fVect);
    T scal = (eVect_T_TheGramMat.dot(fVect_T)) / rNorm;
    MyVector<T> rVect = fVect_T - eVect_T * scal;
    AssignMatrixRow(TheProj, i, rVect);
  }
  os << "RNK: GetOrthogonalProjector_dim1, step 9\n";
  return {TheProj, TheCompl};
}

template <typename T, typename Tint>
T UpperBoundRankinMinimalDeterminant(MyMatrix<T> const &TheGramMat, int k,
                                     std::ostream &os) {
  os << "RNK: UpperBoundRankinMinimalDeterminant k=" << k << "\n";
  T_shvec_info<T, Tint> SHVmin = computeMinimum_GramMat<T, Tint>(TheGramMat);
  T rNorm = SHVmin.minimum;
  os << "RNK: UpperBoundRankinMinimalDeterminant rNorm=" << rNorm << "\n";
  if (k == 1) {
    return rNorm;
  }
  MyVector<Tint> eVect = SHVmin.short_vectors[0];
  os << "RNK: UpperBoundRankinMinimalDeterminant eVect=" << StringVector(eVect)
     << "\n";
  MyMatrix<T> TheProj =
      GetOrthogonalProjector_dim1(TheGramMat, eVect, os).first;
  os << "RNK: UpperBoundRankinMinimalDeterminant TheProj\n";
  MyMatrix<T> ReducedGramMat = TheProj * TheGramMat * TheProj.transpose();
  os << "RNK: UpperBoundRankinMinimalDeterminant ReducedGramMat\n";
  T upper =
      UpperBoundRankinMinimalDeterminant<T, Tint>(ReducedGramMat, k - 1, os);
  return rNorm * upper;
}

template <typename T> T get_multiple_value(MyMatrix<T> const &M) {
  FractionMatrix<T> fr = RemoveFractionMatrixPlusCoeff(M);
  int n = M.rows();
  auto get_mat_mult = [&]() -> T {
    for (int i = 0; i < n; i++) {
      T quot = fr.TheMat(i, i) / 2;
      if (!IsInteger(quot)) {
        return 1;
      }
    }
    return 2;
  };
  return get_mat_mult() / fr.TheMult;
}

// Find the upper bound
// --- For floating point types, use the power function
// --- For exact types, like GMP we can use the denominators of M
//     to get to the bound
// We want to find the maximum feasible value of m such that m^k <= C
template <typename T>
T MaxKBound(T const &C, int const &k, MyMatrix<T> const &M) {
  double expo = 1.0 / static_cast<double>(k);
  if constexpr (std::is_same_v<T, double>) {
    return std::pow(C, expo);
  }
  if constexpr (!std::is_same_v<T, double>) {
    // Assumed to be exact arithmetic.
    // pretty inefficient but actually fine.
    T delta = get_multiple_value(M);
    auto is_corr = [&](T const &val) -> T {
      T ret = val;
      for (int i = 1; i < k; i++) {
        ret *= val;
      }
      return ret <= C;
    };
    T val = 0;
    while (true) {
      T val_new = val + delta;
      if (is_corr(val_new)) {
        val = val_new;
      } else {
        return val;
      }
    }
  }
}

template <typename T, typename Tint> struct VectorProjection {
  MyMatrix<T> TheProj;
  MyMatrix<Tint> TheCompl;
  MyVector<Tint> eV;
  T rNorm;
  MyMatrix<T> ReducedGramMat;
};

template <typename T, typename Tint>
VectorProjection<T, Tint> GetVectorProjection(MyMatrix<T> const &TheGramMat,
                                              MyVector<Tint> const &eV,
                                              std::ostream &os) {
  os << "RNK: Beginning of GetVectorProjection\n";
  if (TheGramMat.rows() != eV.size()) {
    std::cerr << "|TheGramMat|=" << TheGramMat.rows() << " |eV|=" << eV.size()
              << "\n";
    std::cerr << "Error: They should be equal\n";
    throw TerminalException{1};
  }
  MyVector<T> eV_T = UniversalVectorConversion<T, Tint>(eV);
  MyVector<T> eV_T_gram = eV_T.transpose() * TheGramMat;
  T rNorm = eV_T_gram.dot(eV_T);
  std::pair<MyMatrix<T>, MyMatrix<Tint>> pair =
      GetOrthogonalProjector_dim1(TheGramMat, eV, os);
  MyMatrix<T> const &TheProj = pair.first;
  MyMatrix<Tint> const &TheCompl = pair.second;
  MyMatrix<T> ReducedGramMat = TheProj * TheGramMat * TheProj.transpose();
  return {std::move(TheProj), std::move(TheCompl), eV, std::move(rNorm),
          std::move(ReducedGramMat)};
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtendSublattice(VectorProjection<T, Tint> const &vp,
                                MyMatrix<Tint> const &eLatt, std::ostream &os) {
  MyMatrix<Tint> ePart = eLatt * vp.TheCompl;
  MyMatrix<Tint> fLatt = ConcatenateMatVec(ePart, vp.eV);
  os << "RNK: ExtendSublattice, fLatt=\n";
  WriteMatrix(os, fLatt);
  MyMatrix<Tint> gLatt = ComputeRowHermiteNormalForm_second(fLatt);
  os << "RNK: ExtendSublattice, gLatt=\n";
  WriteMatrix(os, gLatt);
  return gLatt;
}

// The function that returns the Rankin k-minimum.
// It should work both for floating point types and exact types.
// If the threshold need to be used, then that should be in the MaxDet.
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>> Rankin_k_level(MyMatrix<T> const &A, int const &k,
                                           T const &MaxDet, std::ostream &os) {
  if (k == 1) {
    T bound = MaxDet;
    T_shvec_info<T, Tint> SHVmin = computeLevel_GramMat<T, Tint>(A, bound);
    std::vector<MyMatrix<Tint>> RetList;
    for (auto &eV : SHVmin.short_vectors) {
      MyMatrix<Tint> M = MatrixFromVector(eV);
      RetList.push_back(M);
    }
    return RetList;
  }
  // We use the HermiteNormalForm
  std::unordered_set<MyMatrix<Tint>> set_subspaces;
  // We are now using the Hermite constant to get a bound on the minimum
  // That is we have min(A)^k <= H(n) * MaxDet
  T upper = GetUpperBoundHermitePower<T>(k) * MaxDet;
  T bound = MaxKBound(upper, k, A);
  T_shvec_info<T, Tint> SHVmin = computeLevel_GramMat<T, Tint>(A, bound);
  for (auto &eV : SHVmin.short_vectors) {
    VectorProjection<T, Tint> vp = GetVectorProjection(A, eV, os);
    T TheAskDet = MaxDet / vp.rNorm;
    std::vector<MyMatrix<Tint>> SpecEnum =
        Rankin_k_level<T, Tint>(vp.ReducedGramMat, k - 1, TheAskDet, os);
    for (auto &eLatt : SpecEnum) {
      MyMatrix<Tint> gLatt = ExtendSublattice(vp, eLatt, os);
#ifdef DEBUG_RANKIN
      MyMatrix<T> gLatt_T = UniversalMatrixConversion<T, Tint>(gLatt);
      MyMatrix<T> eProdMat = gLatt_T * TheGramMat * gLatt_T.transpose();
      T eDet = DeterminantMat(eProdMat);
      if (eDet > MaxDet) {
        std::cerr << "We have eDet=" << eDet << " MaxDet=" << MaxDet << "\n";
        std::cerr << "Incoherent result\n";
        throw TerminalException{1};
      }
#endif
      set_subspaces.insert(gLatt);
    }
  }
  std::vector<MyMatrix<Tint>> vec_subspaces;
  for (auto &mat : set_subspaces) {
    vec_subspaces.push_back(mat);
  }
  return vec_subspaces;
}

template <typename T, typename Tint> struct ResultKRankinMin {
  T min;
  std::vector<MyMatrix<Tint>> l_spaces;
};

template <typename T, typename Tint>
ResultKRankinMin<T, Tint> Rankin_k_minimum(MyMatrix<T> const &A, int const &k,
                                           T const &tol, std::ostream &os) {
  if (k == 1) {
    T bound = UpperBoundRankinMinimalDeterminant<T, Tint>(A, 1, os);
    T bound_search = bound * (1 + tol);
    os << "RNK: k=" << k << " bound=" << bound
       << " bound_search=" << bound_search << "\n";
    T_shvec_info<T, Tint> SHVmin =
        computeLevel_GramMat<T, Tint>(A, bound_search);
    os << "RNK: |SHVmin.short_vectors|=" << SHVmin.short_vectors.size() << "\n";
    std::vector<MyMatrix<Tint>> RetList;
    for (auto &eV : SHVmin.short_vectors) {
      MyMatrix<Tint> M = MatrixFromVector(eV);
      RetList.push_back(M);
    }
    return {bound, RetList};
  }
  // We use the HermiteNormalForm
  T DetMin = 0;
  std::unordered_set<MyMatrix<Tint>> set_subspaces;
  auto f_insert = [&](MyMatrix<Tint> const &gLatt) -> void {
    MyMatrix<T> gLatt_T = UniversalMatrixConversion<T, Tint>(gLatt);
    MyMatrix<T> eProdMat = gLatt_T * A * gLatt_T.transpose();
    T eDet = DeterminantMat(eProdMat);
    os << "RNK: f_insert eDet=" << eDet << "\n";
    if (set_subspaces.size() == 0) {
      DetMin = eDet;
      os << "RNK: Now DetMin=" << DetMin << " Case 1\n";
      set_subspaces.insert(gLatt);
    } else {
      if (eDet < DetMin * (1 - tol)) {
        set_subspaces.clear();
        DetMin = eDet;
        os << "RNK: Now DetMin=" << DetMin << " Case 2\n";
        set_subspaces.insert(gLatt);
      } else {
        if (eDet <= DetMin * (1 + tol)) {
          set_subspaces.insert(gLatt);
          os << "RNK: Now |set_subspaces|=" << set_subspaces.size() << "\n";
        }
      }
    }
  };
  // We are now using the Hermite constant to get a bound on the minimum
  // That is we have min(A)^k <= H(n) * MaxDet
  os << "RNK: Rankin_k_minimum, step 1\n";
  T MaxDet = UpperBoundRankinMinimalDeterminant<T, Tint>(A, k, os);
  os << "RNK: Rankin_k_minimum, step 2 MaxDet=" << MaxDet << "\n";
  os << "RNK: Rankin_k_minimum, k=" << k
     << " HermitePower=" << GetUpperBoundHermitePower<T>(k) << "\n";
  T upper = GetUpperBoundHermitePower<T>(k) * MaxDet;
  os << "RNK: Rankin_k_minimum, step 3 upper=" << upper << "\n";
  T bound = MaxKBound(upper, k, A);
  os << "RNK: Rankin_k_minimum, step 4 bound=" << bound << "\n";
  T_shvec_info<T, Tint> SHVmin = computeLevel_GramMat<T, Tint>(A, bound);
  os << "RNK: Rankin_k_minimum, step 5\n";
  for (auto &eV : SHVmin.short_vectors) {
    os << "RNK: eV, step 1\n";
    VectorProjection<T, Tint> vp = GetVectorProjection(A, eV, os);
    os << "RNK: eV, step 2\n";
    T TheAskDet = MaxDet / vp.rNorm;
    os << "RNK: eV, step 3\n";
    std::vector<MyMatrix<Tint>> SpecEnum =
        Rankin_k_level<T, Tint>(vp.ReducedGramMat, k - 1, TheAskDet, os);
    os << "RNK: eV, step 4\n";
    for (auto &eLatt : SpecEnum) {
      MyMatrix<Tint> gLatt = ExtendSublattice(vp, eLatt, os);
      f_insert(gLatt);
    }
    os << "RNK: eV, step 5\n";
  }
  os << "RNK: Rankin_k_minimum, step 6\n";
  std::vector<MyMatrix<Tint>> vec_subspaces;
  for (auto &mat : set_subspaces) {
    vec_subspaces.push_back(mat);
  }
  os << "RNK: Rankin_k_minimum, step 7\n";
  return {DetMin, vec_subspaces};
}

// clang-format off
#endif  // SRC_RANKIN_ENUMERATION_K_SPACE_H_
// clang-format on
