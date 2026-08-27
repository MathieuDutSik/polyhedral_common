// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ENUMERATION_K_SPACE_H_
#define SRC_LATT_ENUMERATION_K_SPACE_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include "LatticeStabEquiCan.h"
#include <utility>
#include <unordered_set>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ENUMERATION_K_SPACE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ENUMERATION_K_SPACE
#endif

// Returns the Hermite constant at the n-th power.
// That is the maximum of min(A)^n / det(A)
// --
// So we would have min(A)^n / det(A) <= H(n)
// which gets us min(A)^n <= H(n) det(A)
// --
// The constants are known exactly up to n=8 and for n=9.
// The case n=9 follows from the enumeration of the perfect
// lattices in dimension 9 which gives gamma_9 = 2 and so
// H(9) = 2^9 = 512.
// For higher n, we use the bound H(n) <= (4/3)^(n(n-1)/2)
// obtained from Hermite's inequality.
template <typename T> T GetUpperBoundHermitePower(int n) {
  if (n == 1) {
    return T(1);
  }
  if (n == 2) {
    return T(4) / T(3);
  }
  if (n == 3) {
    return T(2);
  }
  if (n == 4) {
    return T(4);
  }
  if (n == 5) {
    return T(8);
  }
  if (n == 6) {
    return T(64) / T(3);
  }
  if (n == 7) {
    return T(64);
  }
  if (n == 8) {
    return T(256);
  }
  if (n == 9) {
    return T(512);
  }
  int h = n * (n - 1) / 2;
  T base = T(4) / T(3);
  T result(1);
  for (int i = 0; i < h; i++)
    result *= base;
  return result;
}

/*
  Returns the matrix P of the orthogonal projection onto the span of
  TheSubBasis. With B = TheSubBasis and G = TheGramMat we set
  P = B^T (B G B^T)^{-1} B G
  which acts on column vectors; a row vector x projects as x P^T,
  which fixes the row span of B and maps to it. This is the transpose
  of the matrix built by __GetOrthogonalProjector in
  SublatticeEnumeration.g (whose rows are the images of the basis
  vectors, acting on row vectors as x M).
  ---
  An orthogonal projector is self adjoint, that is
  <p(x), y> = <x, p(y)>
  With the scalar product <x,y> = x G y^T and p(x) = x P^T we get
  <p(x), y> = x P^T G y^T
  <x, p(y)> = x G P y^T
  so the condition is P^T G = G P, that is G P is symmetric.
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
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
  if (TheProj * TheProj != TheProj) {
    std::cerr << "The matrix is not a projector\n";
    throw TerminalException{1};
  }
  if (TheProj.transpose() * TheGramMat != TheGramMat * TheProj) {
    std::cerr << "The obtained projector is not self adjoint\n";
    throw TerminalException{1};
  }
#endif
  return TheProj;
}

template <typename T, typename Tint>
std::pair<MyMatrix<T>, MyMatrix<Tint>>
GetOrthogonalProjector_dim1(MyMatrix<T> const &TheGramMat,
                            MyVector<Tint> const &eVect,
                            [[maybe_unused]] std::ostream &os) {
  int n = TheGramMat.rows();
  MyMatrix<Tint> eVect_M = MatrixFromVector(eVect);
  MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect);
  MyVector<T> eVect_T_TheGramMat = eVect_T.transpose() * TheGramMat;
  T rNorm = eVect_T_TheGramMat.dot(eVect_T);
#ifdef DEBUG_ENUMERATION_K_SPACE
  os << "ENUM_K_SPACE: GetOrthogonalProjector_dim1 n=" << n
     << " rNorm=" << rNorm << "\n";
#endif
  MyMatrix<Tint> TheCompl = SubspaceCompletionInt(eVect_M, n);
  MyMatrix<T> TheProj(n - 1, n);
  for (int i = 0; i < n - 1; i++) {
    MyVector<Tint> fVect = GetMatrixRow(TheCompl, i);
    MyVector<T> fVect_T = UniversalVectorConversion<T, Tint>(fVect);
    T scal = (eVect_T_TheGramMat.dot(fVect_T)) / rNorm;
    MyVector<T> rVect = fVect_T - eVect_T * scal;
    AssignMatrixRow(TheProj, i, rVect);
  }
  return {TheProj, TheCompl};
}

template <typename T, typename Tint>
T UpperBoundRankinMinimalDeterminant(MyMatrix<T> const &TheGramMat, int k,
                                     std::ostream &os) {
  Tshortest<T, Tint> RecMin = T_ShortestVector<T, Tint>(TheGramMat, os);
  T rNorm = RecMin.min;
#ifdef DEBUG_ENUMERATION_K_SPACE
  os << "ENUM_K_SPACE: UpperBoundRankinMinimalDeterminant k=" << k
     << " rNorm=" << rNorm << "\n";
#endif
  if (k == 1) {
    return rNorm;
  }
  MyVector<Tint> eVect = GetMatrixRow(RecMin.SHV, 0);
  MyMatrix<T> TheProj =
      GetOrthogonalProjector_dim1(TheGramMat, eVect, os).first;
  MyMatrix<T> ReducedGramMat = TheProj * TheGramMat * TheProj.transpose();
  T upper =
      UpperBoundRankinMinimalDeterminant<T, Tint>(ReducedGramMat, k - 1, os);
  return rNorm * upper;
}

template <typename T> T get_multiple_value(MyMatrix<T> const &M) {
  FractionMatrix<T> fr = RemoveFractionMatrixPlusCoeff(M);
  int n = M.rows();
  auto get_mat_mult = [&]() -> T {
    for (int i = 0; i < n; i++) {
      T quot = fr.TheMat(i, i) / T(2);
      if (!IsInteger(quot)) {
        return T(1);
      }
    }
    return T(2);
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
    auto is_corr = [&](T const &val) -> bool {
      T ret = val;
      for (int i = 1; i < k; i++) {
        ret *= val;
      }
      return ret <= C;
    };
    T val(0);
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

// A k-dimensional sublattice is represented by its row Hermite Normal
// Form which provides a canonical representative and so makes the
// set based computations possible.
template <typename Tint>
MyMatrix<Tint> CanonicalizeSublatticeBasis(MyMatrix<Tint> const &M) {
  return ComputeRowHermiteNormalForm_second(M);
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
                                MyMatrix<Tint> const &eLatt,
                                [[maybe_unused]] std::ostream &os) {
  MyMatrix<Tint> ePart = eLatt * vp.TheCompl;
  MyMatrix<Tint> fLatt = ConcatenateMatVec(ePart, vp.eV);
  MyMatrix<Tint> gLatt = ComputeRowHermiteNormalForm_second(fLatt);
#ifdef DEBUG_ENUMERATION_K_SPACE
  os << "ENUM_K_SPACE: ExtendSublattice, gLatt=\n";
  WriteMatrix(os, gLatt);
#endif
  return gLatt;
}

// The function that returns the k-dimensional subspaces of determinant at
// most MaxDet. A subspace is represented by the saturated sublattice
// L cap V which is encoded by its row Hermite Normal Form. The saturation
// forces us to consider only primitive vectors in the enumeration: any
// saturated sublattice has its minimal vector primitive in the full lattice.
// It should work both for floating point types and exact types.
// If the threshold need to be used, then that should be in the MaxDet.
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>> Rankin_k_level(MyMatrix<T> const &A, int const &k,
                                           T const &MaxDet, std::ostream &os) {
  if (k == 1) {
    T bound = MaxDet;
    std::vector<MyVector<Tint>> short_vectors =
        computeLevel_GramMat<T, Tint>(A, bound, os);
    std::vector<MyMatrix<Tint>> RetList;
    for (auto &eV : short_vectors) {
      if (IsVectorPrimitive(eV)) {
        MyMatrix<Tint> M = MatrixFromVector(eV);
        RetList.push_back(CanonicalizeSublatticeBasis(M));
      }
    }
    return RetList;
  }
  // We use the HermiteNormalForm
  std::unordered_set<MyMatrix<Tint>> set_subspaces;
  // We are now using the Hermite constant to get a bound on the minimum
  // That is we have min(A)^k <= H(n) * MaxDet
  T upper = GetUpperBoundHermitePower<T>(k) * MaxDet;
  T bound = MaxKBound(upper, k, A);
  std::vector<MyVector<Tint>> short_vectors =
      computeLevel_GramMat<T, Tint>(A, bound, os);
  for (auto &eV : short_vectors) {
    if (!IsVectorPrimitive(eV)) {
      continue;
    }
    VectorProjection<T, Tint> vp = GetVectorProjection(A, eV, os);
    T TheAskDet = MaxDet / vp.rNorm;
    std::vector<MyMatrix<Tint>> SpecEnum =
        Rankin_k_level<T, Tint>(vp.ReducedGramMat, k - 1, TheAskDet, os);
    for (auto &eLatt : SpecEnum) {
      MyMatrix<Tint> gLatt = ExtendSublattice(vp, eLatt, os);
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
      MyMatrix<T> gLatt_T = UniversalMatrixConversion<T, Tint>(gLatt);
      MyMatrix<T> eProdMat = gLatt_T * A * gLatt_T.transpose();
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
#ifdef DEBUG_ENUMERATION_K_SPACE
    os << "ENUM_K_SPACE: k=" << k << " bound=" << bound
       << " bound_search=" << bound_search << "\n";
#endif
    std::vector<MyVector<Tint>> short_vectors =
        computeLevel_GramMat<T, Tint>(A, bound_search, os);
#ifdef DEBUG_ENUMERATION_K_SPACE
    os << "ENUM_K_SPACE: |short_vectors|=" << short_vectors.size() << "\n";
#endif
    std::vector<MyMatrix<Tint>> RetList;
    for (auto &eV : short_vectors) {
      if (IsVectorPrimitive(eV)) {
        MyMatrix<Tint> M = MatrixFromVector(eV);
        RetList.push_back(CanonicalizeSublatticeBasis(M));
      }
    }
    return {bound, RetList};
  }
  // We use the HermiteNormalForm
  T DetMin(0);
  std::unordered_set<MyMatrix<Tint>> set_subspaces;
  auto f_insert = [&](MyMatrix<Tint> const &gLatt) -> void {
    MyMatrix<T> gLatt_T = UniversalMatrixConversion<T, Tint>(gLatt);
    MyMatrix<T> eProdMat = gLatt_T * A * gLatt_T.transpose();
    T eDet = DeterminantMat(eProdMat);
    if (set_subspaces.empty()) {
      DetMin = eDet;
      set_subspaces.insert(gLatt);
    } else {
      if (eDet < DetMin * (1 - tol)) {
        set_subspaces.clear();
        DetMin = eDet;
        set_subspaces.insert(gLatt);
      } else {
        if (eDet <= DetMin * (1 + tol)) {
          set_subspaces.insert(gLatt);
        }
      }
    }
#ifdef DEBUG_ENUMERATION_K_SPACE
    os << "ENUM_K_SPACE: f_insert eDet=" << eDet << " DetMin=" << DetMin
       << " |set_subspaces|=" << set_subspaces.size() << "\n";
#endif
  };
  // We are now using the Hermite constant to get a bound on the minimum
  // That is we have min(A)^k <= H(n) * MaxDet
  T MaxDet = UpperBoundRankinMinimalDeterminant<T, Tint>(A, k, os);
  T upper = GetUpperBoundHermitePower<T>(k) * MaxDet;
  T bound = MaxKBound(upper, k, A);
#ifdef DEBUG_ENUMERATION_K_SPACE
  os << "ENUM_K_SPACE: Rankin_k_minimum, k=" << k << " MaxDet=" << MaxDet
     << " upper=" << upper << " bound=" << bound << "\n";
#endif
  std::vector<MyVector<Tint>> short_vectors =
      computeLevel_GramMat<T, Tint>(A, bound, os);
  for (auto &eV : short_vectors) {
    if (!IsVectorPrimitive(eV)) {
      continue;
    }
    VectorProjection<T, Tint> vp = GetVectorProjection(A, eV, os);
    T TheAskDet = MaxDet / vp.rNorm;
    std::vector<MyMatrix<Tint>> SpecEnum =
        Rankin_k_level<T, Tint>(vp.ReducedGramMat, k - 1, TheAskDet, os);
    for (auto &eLatt : SpecEnum) {
      MyMatrix<Tint> gLatt = ExtendSublattice(vp, eLatt, os);
      f_insert(gLatt);
    }
  }
  std::vector<MyMatrix<Tint>> vec_subspaces;
  for (auto &mat : set_subspaces) {
    vec_subspaces.push_back(mat);
  }
  return {DetMin, vec_subspaces};
}

//
// The orbit enumeration code.
// This is a translation of parts of SublatticeEnumeration.g
//

// Computes the orbit of the sublattice eLatt under the group
// generated by ListGen (acting on the right).
// This is a translation of ComputeOrbitSublattice.
template <typename Tint>
std::vector<MyMatrix<Tint>>
OrbitSublattice(std::vector<MyMatrix<Tint>> const &ListGen,
                MyMatrix<Tint> const &eLatt) {
  std::unordered_set<MyMatrix<Tint>> set_latt;
  std::vector<MyMatrix<Tint>> l_latt;
  auto f_insert = [&](MyMatrix<Tint> const &fLatt) -> void {
    if (set_latt.count(fLatt) == 1)
      return;
    set_latt.insert(fLatt);
    l_latt.push_back(fLatt);
  };
  f_insert(CanonicalizeSublatticeBasis(eLatt));
  size_t pos = 0;
  while (pos < l_latt.size()) {
    MyMatrix<Tint> wLatt = l_latt[pos];
    pos++;
    for (auto &eGen : ListGen) {
      MyMatrix<Tint> eImgLatt = wLatt * eGen;
      MyMatrix<Tint> fLatt = CanonicalizeSublatticeBasis(eImgLatt);
      f_insert(fLatt);
    }
  }
  return l_latt;
}

template <typename Tint> struct SublatticeOrbit {
  MyMatrix<Tint> representative;
  size_t orbit_size;
};

// Splits a list of sublattices (assumed closed under the group
// generated by ListGen) into orbits.
template <typename Tint>
std::vector<SublatticeOrbit<Tint>>
OrbitSplittingSublattices(std::vector<MyMatrix<Tint>> const &l_latt,
                          std::vector<MyMatrix<Tint>> const &ListGen) {
  std::unordered_set<MyMatrix<Tint>> set_remain;
  for (auto &eLatt : l_latt) {
    set_remain.insert(CanonicalizeSublatticeBasis(eLatt));
  }
  std::vector<SublatticeOrbit<Tint>> l_orbit;
  while (!set_remain.empty()) {
    MyMatrix<Tint> eLatt = *set_remain.begin();
    std::vector<MyMatrix<Tint>> orbit = OrbitSublattice(ListGen, eLatt);
    for (auto &fLatt : orbit) {
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
      if (set_remain.count(fLatt) != 1) {
        std::cerr << "The list of sublattices is not invariant under the "
                     "group\n";
        throw TerminalException{1};
      }
#endif
      set_remain.erase(fLatt);
    }
    l_orbit.push_back({std::move(eLatt), orbit.size()});
  }
  return l_orbit;
}

// Enumerates the orbits of k-dimensional sublattices of determinant
// at most MaxDet under the automorphism group of the Gram matrix.
template <typename T, typename Tint, typename Tgroup>
std::vector<SublatticeOrbit<Tint>>
Rankin_k_level_orbits(MyMatrix<T> const &A, int const &k, T const &MaxDet,
                      std::ostream &os) {
  std::vector<MyMatrix<Tint>> l_latt = Rankin_k_level<T, Tint>(A, k, MaxDet, os);
  std::vector<MyMatrix<Tint>> ListGen =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(A, os);
#ifdef DEBUG_ENUMERATION_K_SPACE
  os << "ENUM_K_SPACE: Rankin_k_level_orbits |l_latt|=" << l_latt.size()
     << " |ListGen|=" << ListGen.size() << "\n";
#endif
  return OrbitSplittingSublattices(l_latt, ListGen);
}

template <typename T, typename Tint> struct ResultKRankinMinOrbits {
  T min;
  std::vector<SublatticeOrbit<Tint>> l_orbits;
};

// Computes the Rankin k-minimum and the orbits of the k-dimensional
// sublattices realizing it under the automorphism group of the
// Gram matrix. Only for exact arithmetic types.
template <typename T, typename Tint, typename Tgroup>
ResultKRankinMinOrbits<T, Tint>
Rankin_k_minimum_orbits(MyMatrix<T> const &A, int const &k, std::ostream &os) {
  T tol(0);
  ResultKRankinMin<T, Tint> result = Rankin_k_minimum<T, Tint>(A, k, tol, os);
  std::vector<MyMatrix<Tint>> ListGen =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(A, os);
#ifdef DEBUG_ENUMERATION_K_SPACE
  os << "ENUM_K_SPACE: Rankin_k_minimum_orbits min=" << result.min
     << " |l_spaces|=" << result.l_spaces.size()
     << " |ListGen|=" << ListGen.size() << "\n";
#endif
  std::vector<SublatticeOrbit<Tint>> l_orbits =
      OrbitSplittingSublattices(result.l_spaces, ListGen);
  return {std::move(result.min), std::move(l_orbits)};
}

// clang-format off
#endif  // SRC_LATT_ENUMERATION_K_SPACE_H_
// clang-format on
