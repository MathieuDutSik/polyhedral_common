// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_ENUMERATION_K_SPACE_H_
#define SRC_LATT_ENUMERATION_K_SPACE_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include "Shvec_exact.h"
#include "LatticeStabEquiCan.h"
#include <algorithm>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ENUMERATION_K_SPACE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ENUMERATION_K_SPACE
#endif

// The extensive checks re-run the full set-based enumeration in order
// to cross-validate the orbit-wise one, which multiplies the runtime
// by a large factor. They are therefore never enabled by the blanket
// SANITY_CHECK; define SANITY_CHECK_ENUMERATION_K_SPACE_EXTENSIVE
// explicitly to run them.

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
// at most MaxDet under the automorphism group of the Gram matrix,
// by generating the full set of sublattices and splitting it into
// orbits. Simple and reliable, but the full set can be much larger
// than the set of orbits; used as the reference implementation for
// checking the orbit-wise enumeration below.
template <typename T, typename Tint, typename Tgroup>
std::vector<SublatticeOrbit<Tint>>
Rankin_k_level_orbits_setbased(MyMatrix<T> const &A, int const &k,
                               T const &MaxDet, std::ostream &os) {
  std::vector<MyMatrix<Tint>> l_latt = Rankin_k_level<T, Tint>(A, k, MaxDet, os);
  std::vector<MyMatrix<Tint>> ListGen =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(A, os);
#ifdef DEBUG_ENUMERATION_K_SPACE
  os << "ENUM_K_SPACE: Rankin_k_level_orbits_setbased |l_latt|="
     << l_latt.size() << " |ListGen|=" << ListGen.size() << "\n";
#endif
  return OrbitSplittingSublattices(l_latt, ListGen);
}

//
// Orbit-wise enumeration.
// Instead of generating the full set of subspaces and then reducing by
// the group, we enumerate the orbits of the (k-1)-dimensional subspaces
// and for each representative X extend by the orbits, under the
// stabilizer of X, of the primitive short vectors of the projection of
// the lattice on the orthogonal complement of X. The candidates
// obtained are then reduced by equivalence. This mirrors the strategy
// of DoSpecificEnumeration in SublatticeEnumeration.g.
//
// The stabilizer of a subspace X and the equivalence of two subspaces
// are computed on a colored vector configuration: an invariant
// Z-spanning family W of L (color 0) together with the invariant
// family of the sublattice L cap X embedded into L (color 1), handled
// by the ListMat_Vdiag machinery. Any isometric color preserving
// permutation of the configuration lifts to an integral matrix (W
// spans L over Z) preserving the span of the color 1 vectors, that is
// stabilizing X; and conversely.
//

template <typename T, typename Tint> struct SublatticeStabEquiData {
  MyMatrix<T> GramMat;
  std::vector<MyMatrix<T>> ListMat;
  MyMatrix<Tint> W;
  MyMatrix<T> W_T;
};

template <typename T, typename Tint>
SublatticeStabEquiData<T, Tint>
GetSublatticeStabEquiData(MyMatrix<T> const &GramMat, std::ostream &os) {
  MyMatrix<Tint> W = ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, os);
  MyMatrix<T> W_T = UniversalMatrixConversion<T, Tint>(W);
  std::vector<MyMatrix<T>> ListMat{GramMat};
  return {GramMat, std::move(ListMat), std::move(W), std::move(W_T)};
}

template <typename T, typename Tint> struct SublatticeConfiguration {
  MyMatrix<Tint> X;
  T det;
  // The rows of W followed by the rows of invar(X); colors 0 and 1.
  MyMatrix<T> conf;
  std::vector<T> Vdiag;
  // The sorted norms of the rows of invar(X), a cheap invariant used
  // to avoid expensive equivalence tests.
  std::vector<T> invar_norms;
};

template <typename T, typename Tint>
SublatticeConfiguration<T, Tint>
BuildSublatticeConfiguration(SublatticeStabEquiData<T, Tint> const &data,
                             MyMatrix<Tint> const &X, std::ostream &os) {
  int n = data.GramMat.rows();
  MyMatrix<T> X_T = UniversalMatrixConversion<T, Tint>(X);
  MyMatrix<T> eGram = X_T * data.GramMat * X_T.transpose();
  T det = DeterminantMat(eGram);
  MyMatrix<Tint> SHV_X = ExtractInvariantVectorFamilyZbasis<T, Tint>(eGram, os);
  MyMatrix<Tint> invar = SHV_X * X;
  MyMatrix<T> invar_T = UniversalMatrixConversion<T, Tint>(invar);
  int n_w = data.W.rows();
  int n_i = invar.rows();
  MyMatrix<T> conf(n_w + n_i, n);
  std::vector<T> Vdiag(n_w + n_i);
  for (int i = 0; i < n_w; i++) {
    conf.row(i) = data.W_T.row(i);
    Vdiag[i] = T(0);
  }
  std::vector<T> invar_norms(n_i);
  for (int i = 0; i < n_i; i++) {
    conf.row(n_w + i) = invar_T.row(i);
    Vdiag[n_w + i] = T(1);
    MyVector<T> eV = GetMatrixRow(invar_T, i);
    invar_norms[i] = eV.dot(data.GramMat * eV);
  }
  std::sort(invar_norms.begin(), invar_norms.end());
  return {X, std::move(det), std::move(conf), std::move(Vdiag),
          std::move(invar_norms)};
}

// The generators of the stabilizer of the subspace of sc in the
// automorphism group of the Gram matrix.
template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
SublatticeStabilizerGenerators(SublatticeStabEquiData<T, Tint> const &data,
                               SublatticeConfiguration<T, Tint> const &sc,
                               std::ostream &os) {
  std::vector<MyMatrix<T>> LGen = GetIntAutomorphism_ListMat_Vdiag<T, Tgroup>(
      sc.conf, data.ListMat, sc.Vdiag, os);
  std::vector<MyMatrix<Tint>> LGenRet;
  for (auto &g_T : LGen) {
    MyMatrix<Tint> g = UniversalMatrixConversion<Tint, T>(g_T);
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
    MyMatrix<Tint> Ximg = sc.X * g;
    if (CanonicalizeSublatticeBasis(Ximg) != sc.X) {
      std::cerr << "SublatticeStabilizerGenerators: a generator does not "
                   "stabilize the subspace\n";
      throw TerminalException{1};
    }
#endif
    LGenRet.push_back(g);
  }
  return LGenRet;
}

// Returns, if it exists, an integral matrix P preserving the Gram
// matrix with X1 P = X2 as sublattices.
template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
SublatticeTestEquivalence(SublatticeStabEquiData<T, Tint> const &data,
                          SublatticeConfiguration<T, Tint> const &sc1,
                          SublatticeConfiguration<T, Tint> const &sc2,
                          std::ostream &os) {
  if (sc1.det != sc2.det) {
    return {};
  }
  if (sc1.invar_norms != sc2.invar_norms) {
    return {};
  }
  // TestIntEquivalence_ListMat_Vdiag(A, ..., B, ...) returns a matrix
  // mapping the rows of the configuration B onto the rows of the
  // configuration A, so the arguments are swapped to obtain X1 -> X2.
  std::optional<MyMatrix<T>> opt = TestIntEquivalence_ListMat_Vdiag<T, Tgroup>(
      sc2.conf, data.ListMat, sc2.Vdiag, sc1.conf, data.ListMat, sc1.Vdiag, os);
  if (!opt) {
    return {};
  }
  MyMatrix<Tint> P = UniversalMatrixConversion<Tint, T>(*opt);
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
  MyMatrix<Tint> Ximg = sc1.X * P;
  if (CanonicalizeSublatticeBasis(Ximg) != sc2.X) {
    std::cerr << "SublatticeTestEquivalence: the matrix does not map the "
                 "first subspace to the second\n";
    throw TerminalException{1};
  }
#endif
  return P;
}

// The order of the matrix group generated by LGen from its faithful
// permutation action on the rows of the invariant family W.
template <typename Tint, typename Tgroup>
typename Tgroup::Tint
MatrixGroupOrderOnFamily(std::vector<MyMatrix<Tint>> const &LGen,
                         MyMatrix<Tint> const &W) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n_row = W.rows();
  std::unordered_map<MyVector<Tint>, Tidx> map_row;
  for (int i = 0; i < n_row; i++) {
    MyVector<Tint> eV = GetMatrixRow(W, i);
    map_row[eV] = static_cast<Tidx>(i);
  }
  std::vector<Telt> ListPermGens;
  for (auto &g : LGen) {
    MyMatrix<Tint> g_tr = g.transpose();
    std::vector<Tidx> eList(n_row);
    for (int i = 0; i < n_row; i++) {
      MyVector<Tint> eV = GetMatrixRow(W, i);
      MyVector<Tint> fV = g_tr * eV;
      auto iter = map_row.find(fV);
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
      if (iter == map_row.end()) {
        std::cerr << "MatrixGroupOrderOnFamily: the family is not invariant "
                     "under the group\n";
        throw TerminalException{1};
      }
#endif
      eList[i] = iter->second;
    }
    ListPermGens.push_back(Telt(eList));
  }
  Tgroup grp(ListPermGens, n_row);
  return grp.size();
}

template <typename Tgroup>
size_t GroupIndexAsSizeT(typename Tgroup::Tint const &ord_big,
                         typename Tgroup::Tint const &ord_small) {
  typename Tgroup::Tint quot = ord_big / ord_small;
  std::stringstream ss;
  ss << quot;
  size_t val;
  std::stringstream(ss.str()) >> val;
  return val;
}

template <typename T, typename Tint> struct SubspaceProjection {
  MyMatrix<Tint> TheCompl;
  MyMatrix<T> RedGram;
  MyMatrix<T> Full_T;
  MyMatrix<T> FullInv_T;
};

// The projection of the lattice on the orthogonal complement of the
// span of X, together with the completion data needed for lifting
// vectors and transporting the stabilizer.
template <typename T, typename Tint>
SubspaceProjection<T, Tint>
GetSubspaceProjection(MyMatrix<T> const &GramMat, MyMatrix<Tint> const &X) {
  int n = GramMat.rows();
  MyMatrix<Tint> TheCompl = SubspaceCompletionInt(X, n);
  MyMatrix<T> TheCompl_T = UniversalMatrixConversion<T, Tint>(TheCompl);
  MyMatrix<T> TheProj = GetOrthogonalProjector(GramMat, X);
  MyMatrix<T> ProjCompl = TheCompl_T - TheCompl_T * TheProj.transpose();
  MyMatrix<T> RedGram = ProjCompl * GramMat * ProjCompl.transpose();
  MyMatrix<Tint> Full = Concatenate(X, TheCompl);
  MyMatrix<T> Full_T = UniversalMatrixConversion<T, Tint>(Full);
  MyMatrix<T> FullInv_T = Inverse(Full_T);
  return {std::move(TheCompl), std::move(RedGram), std::move(Full_T),
          std::move(FullInv_T)};
}

// The action induced on the projected lattice, in the coordinates of
// the projected completion basis, by a matrix g stabilizing the
// subspace: the bottom right block of the conjugated matrix.
template <typename T, typename Tint>
MyMatrix<Tint>
InducedProjectedAction(SubspaceProjection<T, Tint> const &sp, int const &dim_x,
                       MyMatrix<Tint> const &g) {
  int n = sp.Full_T.rows();
  int m = n - dim_x;
  MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
  MyMatrix<T> N_T = sp.Full_T * g_T * sp.FullInv_T;
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
  if (!IsIntegralMatrix(N_T)) {
    std::cerr << "InducedProjectedAction: the conjugated matrix should be "
                 "integral\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < dim_x; i++) {
    for (int j = 0; j < m; j++) {
      if (N_T(i, dim_x + j) != T(0)) {
        std::cerr << "InducedProjectedAction: the matrix does not stabilize "
                     "the subspace\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  MyMatrix<Tint> gp(m, m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      gp(i, j) = UniversalScalarConversion<Tint, T>(N_T(dim_x + i, dim_x + j));
    }
  }
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
  MyMatrix<T> gp_T = UniversalMatrixConversion<T, Tint>(gp);
  if (gp_T * sp.RedGram * gp_T.transpose() != sp.RedGram) {
    std::cerr << "InducedProjectedAction: the induced action does not "
                 "preserve the projected Gram matrix\n";
    throw TerminalException{1};
  }
#endif
  return gp;
}

// An upper bound for R^{1/k}, computed by exact bisection so that no
// floating point is involved. The bound is within a factor 2^{-20} of
// the true root, which is only a matter of speed, not of correctness.
template <typename T> T UpperKthRootBound(T const &R, int const &k) {
  if (R <= T(0)) {
    return T(0);
  }
  auto pow_ge = [&](T const &x) -> bool {
    T p(1);
    for (int i = 0; i < k; i++) {
      p *= x;
    }
    return p >= R;
  };
  T lo(0);
  T hi(1);
  while (!pow_ge(hi)) {
    lo = hi;
    hi *= T(2);
  }
  for (int iter = 0; iter < 20; iter++) {
    T mid = (lo + hi) / T(2);
    if (pow_ge(mid)) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return hi;
}

template <typename Tint> struct SublatticeRepInfo {
  MyMatrix<Tint> X;
  std::vector<MyMatrix<Tint>> stab_gens;
};

// The recursive kernel of the orbit-wise enumeration: the orbit
// representatives of the k-dimensional primitive sublattices of
// determinant at most MaxDet, together with generators of their
// stabilizers.
// Completeness of the recursion: a k-dimensional sublattice Y with
// det(Y) <= M contains a (k-1)-dimensional sublattice X (which can be
// saturated, only decreasing the determinant) with a bounded
// determinant, in two ways:
// * By the Minkowski and Hadamard inequalities, the span of vectors
//   realizing the first k-1 successive minima has
//   det(X) <= prod_{j<k} lambda_j(Y)^2
//   <= gamma_k^k det(Y) / lambda_k(Y)^2 <= gamma_k^k M / min(A).
// * By the Rankin inequality and the Rankin duality
//   gamma_{k,k-1} = gamma_{k,1} = gamma_k, we have
//   det(X)^k <= gamma_k^k det(Y)^{k-1} <= gamma_k^k M^{k-1}.
// The minimum of the two bounds is used; the Rankin bound makes the
// bounds shrink along the recursion instead of growing.
template <typename T, typename Tint, typename Tgroup>
std::vector<SublatticeRepInfo<Tint>>
Rankin_k_level_repinfo(SublatticeStabEquiData<T, Tint> const &data,
                       std::vector<MyMatrix<Tint>> const &l_autgen,
                       T const &minA, int const &k, T const &MaxDet,
                       std::ostream &os) {
  MyMatrix<T> const &A = data.GramMat;
  auto get_reps = [&]() -> std::vector<MyMatrix<Tint>> {
    if (k == 1) {
      // The full orbit of vectors has to be enumerated and reduced by
      // the group; this cannot be avoided.
      std::vector<MyVector<Tint>> short_vectors =
          computeLevel_GramMat<T, Tint>(A, MaxDet, os);
      std::vector<MyMatrix<Tint>> l_latt;
      for (auto &eV : short_vectors) {
        if (IsVectorPrimitive(eV)) {
          MyMatrix<Tint> M = MatrixFromVector(eV);
          l_latt.push_back(CanonicalizeSublatticeBasis(M));
        }
      }
      std::vector<SublatticeOrbit<Tint>> l_orb =
          OrbitSplittingSublattices(l_latt, l_autgen);
      std::vector<MyMatrix<Tint>> reps;
      for (auto &e_orb : l_orb) {
        reps.push_back(e_orb.representative);
      }
      return reps;
    }
    T HermPow = GetUpperBoundHermitePower<T>(k);
    T MprevA = HermPow * MaxDet / minA;
    T Mpow(1);
    for (int i = 0; i < k - 1; i++) {
      Mpow *= MaxDet;
    }
    T MprevB = UpperKthRootBound<T>(HermPow * Mpow, k);
    T Mprev = (MprevA < MprevB) ? MprevA : MprevB;
    std::vector<SublatticeRepInfo<Tint>> l_prev =
        Rankin_k_level_repinfo<T, Tint, Tgroup>(data, l_autgen, minA, k - 1,
                                                Mprev, os);
#ifdef DEBUG_ENUMERATION_K_SPACE
    os << "ENUM_K_SPACE: repinfo k=" << k << " MaxDet=" << MaxDet
       << " Mprev=" << Mprev << " |l_prev|=" << l_prev.size() << "\n";
#endif
    std::map<std::pair<T, std::vector<T>>, std::vector<size_t>> map_cand;
    std::vector<SublatticeConfiguration<T, Tint>> accepted;
    std::vector<MyMatrix<Tint>> reps;
    auto f_insert = [&](MyMatrix<Tint> const &Y) -> void {
      SublatticeConfiguration<T, Tint> sc =
          BuildSublatticeConfiguration(data, Y, os);
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE
      if (sc.det > MaxDet) {
        std::cerr << "Rankin_k_level_repinfo: a candidate exceeds MaxDet\n";
        throw TerminalException{1};
      }
#endif
      std::pair<T, std::vector<T>> key{sc.det, sc.invar_norms};
      std::vector<size_t> &l_idx = map_cand[key];
      for (auto &idx : l_idx) {
        std::optional<MyMatrix<Tint>> opt =
            SublatticeTestEquivalence<T, Tint, Tgroup>(data, sc, accepted[idx],
                                                       os);
        if (opt) {
          return;
        }
      }
      l_idx.push_back(accepted.size());
      reps.push_back(Y);
      accepted.push_back(std::move(sc));
    };
    for (auto &rep : l_prev) {
      SubspaceProjection<T, Tint> sp = GetSubspaceProjection(A, rep.X);
      MyMatrix<T> X_T = UniversalMatrixConversion<T, Tint>(rep.X);
      T detX = DeterminantMat(MyMatrix<T>(X_T * A * X_T.transpose()));
      T bound = MaxDet / detX;
      std::vector<MyMatrix<Tint>> l_gp;
      for (auto &g : rep.stab_gens) {
        l_gp.push_back(InducedProjectedAction(sp, k - 1, g));
      }
      std::vector<MyVector<Tint>> l_w =
          computeLevel_GramMat<T, Tint>(sp.RedGram, bound, os);
      std::vector<MyMatrix<Tint>> l_wmat;
      for (auto &w : l_w) {
        if (IsVectorPrimitive(w)) {
          MyMatrix<Tint> M = MatrixFromVector(w);
          l_wmat.push_back(CanonicalizeSublatticeBasis(M));
        }
      }
      std::vector<SublatticeOrbit<Tint>> l_worb =
          OrbitSplittingSublattices(l_wmat, l_gp);
#ifdef DEBUG_ENUMERATION_K_SPACE
      os << "ENUM_K_SPACE: repinfo k=" << k << " detX=" << detX
         << " |l_wmat|=" << l_wmat.size() << " |l_worb|=" << l_worb.size()
         << "\n";
#endif
      for (auto &e_worb : l_worb) {
        MyVector<Tint> w = GetMatrixRow(e_worb.representative, 0);
        MyVector<Tint> wC = sp.TheCompl.transpose() * w;
        MyMatrix<Tint> Ymat = ConcatenateMatVec(rep.X, wC);
        MyMatrix<Tint> Y = CanonicalizeSublatticeBasis(Ymat);
        f_insert(Y);
      }
    }
    return reps;
  };
  std::vector<MyMatrix<Tint>> reps = get_reps();
  std::vector<SublatticeRepInfo<Tint>> result;
  for (auto &X : reps) {
    SublatticeConfiguration<T, Tint> sc =
        BuildSublatticeConfiguration(data, X, os);
    std::vector<MyMatrix<Tint>> stab_gens =
        SublatticeStabilizerGenerators<T, Tint, Tgroup>(data, sc, os);
    result.push_back({X, std::move(stab_gens)});
  }
  return result;
}

// Enumerates the orbits of k-dimensional primitive sublattices of
// determinant at most MaxDet under the automorphism group of the Gram
// matrix, working orbit-wise throughout. The orbit sizes are obtained
// as index of the stabilizer in the automorphism group.
template <typename T, typename Tint, typename Tgroup>
std::vector<SublatticeOrbit<Tint>>
Rankin_k_level_orbits(MyMatrix<T> const &A, int const &k, T const &MaxDet,
                      std::ostream &os) {
  using Tgint = typename Tgroup::Tint;
  SublatticeStabEquiData<T, Tint> data =
      GetSublatticeStabEquiData<T, Tint>(A, os);
  std::vector<MyMatrix<Tint>> l_autgen =
      ArithmeticAutomorphismGroupMultiple_inner<T, Tint, Tgroup>(data.ListMat,
                                                                 data.W, os);
  Tgint aut_order = MatrixGroupOrderOnFamily<Tint, Tgroup>(l_autgen, data.W);
  T minA = T_ShortestVector<T, Tint>(A, os).min;
  std::vector<SublatticeRepInfo<Tint>> l_rep =
      Rankin_k_level_repinfo<T, Tint, Tgroup>(data, l_autgen, minA, k, MaxDet,
                                              os);
  std::vector<SublatticeOrbit<Tint>> l_orbit;
  for (auto &rep : l_rep) {
    Tgint stab_order =
        MatrixGroupOrderOnFamily<Tint, Tgroup>(rep.stab_gens, data.W);
    size_t orbit_size = GroupIndexAsSizeT<Tgroup>(aut_order, stab_order);
    l_orbit.push_back({rep.X, orbit_size});
  }
#ifdef SANITY_CHECK_ENUMERATION_K_SPACE_EXTENSIVE
  // Full cross-validation against the set-based enumeration. This can
  // be much slower than the enumeration being checked, so announce it.
  std::cerr << "ENUM_K_SPACE: entering the extensive sanity check: "
               "set-based cross-validation of the orbit-wise enumeration "
               "for k=" << k << " MaxDet=" << MaxDet << ". This is slow.\n";
  std::vector<SublatticeOrbit<Tint>> l_orbit_ref =
      Rankin_k_level_orbits_setbased<T, Tint, Tgroup>(A, k, MaxDet, os);
  auto get_signature = [&](std::vector<SublatticeOrbit<Tint>> const &l_orb)
      -> std::vector<std::pair<T, size_t>> {
    std::vector<std::pair<T, size_t>> l_sign;
    for (auto &e_orb : l_orb) {
      MyMatrix<T> Y_T =
          UniversalMatrixConversion<T, Tint>(e_orb.representative);
      T det = DeterminantMat(MyMatrix<T>(Y_T * A * Y_T.transpose()));
      l_sign.push_back({det, e_orb.orbit_size});
    }
    std::sort(l_sign.begin(), l_sign.end());
    return l_sign;
  };
  if (get_signature(l_orbit) != get_signature(l_orbit_ref)) {
    std::cerr << "Rankin_k_level_orbits: the orbit-wise enumeration does not "
                 "match the set-based one\n";
    std::cerr << "k=" << k << " MaxDet=" << MaxDet
              << " |l_orbit|=" << l_orbit.size()
              << " |l_orbit_ref|=" << l_orbit_ref.size() << "\n";
    throw TerminalException{1};
  }
#endif
  return l_orbit;
}

template <typename T, typename Tint> struct ResultKRankinMinOrbits {
  T min;
  std::vector<SublatticeOrbit<Tint>> l_orbits;
};

// Computes the Rankin k-minimum and the orbits of the k-dimensional
// sublattices realizing it under the automorphism group of the
// Gram matrix, working orbit-wise: all the orbits below the upper
// bound for the minimum are enumerated and the ones of minimal
// determinant are retained. Only for exact arithmetic types.
template <typename T, typename Tint, typename Tgroup>
ResultKRankinMinOrbits<T, Tint>
Rankin_k_minimum_orbits(MyMatrix<T> const &A, int const &k, std::ostream &os) {
  T MaxDet = UpperBoundRankinMinimalDeterminant<T, Tint>(A, k, os);
  std::vector<SublatticeOrbit<Tint>> l_orbit =
      Rankin_k_level_orbits<T, Tint, Tgroup>(A, k, MaxDet, os);
  if (l_orbit.empty()) {
    std::cerr << "Rankin_k_minimum_orbits: the upper bound should be "
                 "attained by some sublattice\n";
    throw TerminalException{1};
  }
  auto get_det = [&](MyMatrix<Tint> const &Y) -> T {
    MyMatrix<T> Y_T = UniversalMatrixConversion<T, Tint>(Y);
    return DeterminantMat(MyMatrix<T>(Y_T * A * Y_T.transpose()));
  };
  std::vector<T> l_det;
  T min = get_det(l_orbit[0].representative);
  l_det.push_back(min);
  for (size_t u = 1; u < l_orbit.size(); u++) {
    T det = get_det(l_orbit[u].representative);
    if (det < min) {
      min = det;
    }
    l_det.push_back(det);
  }
  std::vector<SublatticeOrbit<Tint>> l_orbits_min;
  for (size_t u = 0; u < l_orbit.size(); u++) {
    if (l_det[u] == min) {
      l_orbits_min.push_back(l_orbit[u]);
    }
  }
#ifdef DEBUG_ENUMERATION_K_SPACE
  os << "ENUM_K_SPACE: Rankin_k_minimum_orbits min=" << min
     << " |l_orbit|=" << l_orbit.size()
     << " |l_orbits_min|=" << l_orbits_min.size() << "\n";
#endif
  return {std::move(min), std::move(l_orbits_min)};
}

// clang-format off
#endif  // SRC_LATT_ENUMERATION_K_SPACE_H_
// clang-format on
