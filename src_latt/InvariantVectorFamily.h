// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_INVARIANTVECTORFAMILY_H_
#define SRC_LATT_INVARIANTVECTORFAMILY_H_

// clang-format off
#include "Shvec_exact.h"
#include "Positivity.h"
#include "Timings.h"
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <utility>
#include <map>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_INVARIANT_VECTOR_FAMILY
#endif

#ifdef DEBUG
#define DEBUG_INVARIANT_VECTOR_FAMILY
#endif

#ifdef DISABLE_DEBUG_INVARIANT_VECTOR_FAMILY
#undef DEBUG_INVARIANT_VECTOR_FAMILY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_INVARIANT_VECTOR_FAMILY
#endif

/*
  We are considering here the enumeration of configurations of vectors for
  positive definite quadratic forms. There are many possible optimization
  and choices.

  The basis assumption of the code are the following:
  --- If enumerating vectors at norm N, then the cost for vectors of norms
  smaller than N is very small compared to the cost for the vectors of norm
  exactly N. Therefore, it is ok, to redo many small enumeration until we get
  the configuration of vector that we want.
  --- The increment are usually monotonous like norm 2,3,4. So, going with
  the increment given by GCD is probably a good heuristic.

 */

template <typename T, typename Tint>
T GetMaxNorm(MyMatrix<T> const &eMat, std::ostream &os) {
  LLLreduction<T, Tint> recLLL = LLLreducedBasis<T, Tint>(eMat, os);
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(recLLL.Pmat);
  MyMatrix<T> eMatRed = Pmat_T * eMat * TransposedMat(Pmat_T);
  return MaximumDiagonal(eMatRed);
}

template <typename T> T GetSmallestIncrement(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  T eGcd = eMat(0, 0);
  // No early exit at eGcd = 1: over the rationals the gcd can decrease
  // below 1 (gcd(1, 2/3) = 1/3), so 1 is not a stopping point. A rational
  // Gram matrix does occur, e.g. the realizability matrices of the
  // short-vector configurations.
  for (int i = 1; i < n; i++) {
    eGcd = GcdPair(eGcd, eMat(i, i));
  }
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      T val = 2 * eMat(i, j);
      eGcd = GcdPair(eGcd, val);
    }
  }
  return eGcd;
}

template <typename T, typename Tint>
std::set<T> GetSetNormConsider(MyMatrix<T> const &eMat, std::ostream &os) {
  T incr = GetSmallestIncrement(eMat);
  T MaxNorm = GetMaxNorm<T, Tint>(eMat, os);
  std::set<T> AllowedNorms;
  T norm = incr;
  while (true) {
    if (norm > MaxNorm)
      break;
    AllowedNorms.insert(norm);
    norm += incr;
  }
  return AllowedNorms;
}

/*
  The rank ought to be self explanatory.
  The index is the index in the lattice obtaine by the saturation
  of the spanned lattce
 */
template <typename Tint> struct FundInvariantVectorFamily {
  int rank;
  Tint index;
};

template <typename Tint>
FundInvariantVectorFamily<Tint> TrivFundamentalInvariant() {
  return {0, 1};
}

template <typename Tint>
bool IsCompleteSystem(FundInvariantVectorFamily<Tint> const &fi, int n) {
  if (fi.rank != n) {
    return false;
  }
  return fi.index == 1;
}

template <typename T> bool is_antipodal(MyMatrix<T> const &SHV) {
  int n_row = SHV.rows();
  std::unordered_set<MyVector<T>> set;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = GetMatrixRow(SHV, i_row);
    set.insert(V);
  }
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = -GetMatrixRow(SHV, i_row);
    if (!set.contains(V)) {
      return false;
    }
  }
  return true;
}

template <typename T> MyMatrix<T> matrix_duplication(MyMatrix<T> const &SHV) {
  int dim = SHV.cols();
  int nbSHV = SHV.rows();
  MyMatrix<T> SHV_T(2 * nbSHV, dim);
  for (int i_row = 0; i_row < nbSHV; i_row++) {
    for (int i = 0; i < dim; i++) {
      T val = SHV(i_row, i);
      SHV_T(2 * i_row, i) = val;
      SHV_T(2 * i_row + 1, i) = -val;
    }
  }
  return SHV_T;
}

template <typename T> bool has_duplication(MyMatrix<T> const &SHV) {
  int n_row = SHV.rows();
  std::unordered_set<MyVector<T>> set;
  for (int i_row = 0; i_row < n_row; i_row++) {
    MyVector<T> V = GetMatrixRow(SHV, i_row);
    if (!set.insert(V).second) {
      return true;
    }
  }
  return false;
}

template <typename T> void check_antipodality_mymatrix(MyMatrix<T> const &SHV) {
  if (!is_antipodal(SHV)) {
    std::cerr << "TSPACE: The family SHV is not antipodal\n";
    throw TerminalException{1};
  }
}

template <typename Tint> bool IsFullDimZbasis(MyMatrix<Tint> const &M, [[maybe_unused]] std::ostream &os) {
  int n = M.cols();
  int rnk = RankMat(M);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: IsFullDimZbasis rnk=" << rnk << " n=" << n << "\n";
#endif
  if (rnk < n) {
    return false;
  }
  Tint indx = Int_IndexLattice(M);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: IsFullDimZbasis indx=" << indx << "\n";
#endif
  if (T_NormGen(indx) == 1) {
    return true;
  }
  return false;
}

template <typename Tint> bool IsFullDim(MyMatrix<Tint> const &M, [[maybe_unused]] std::ostream &os) {
  int n = M.cols();
  int rnk = RankMat(M);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: IsFullDimZbasis rnk=" << rnk << " n=" << n << "\n";
#endif
  if (rnk < n) {
    return false;
  }
  return true;
}

template <typename Tint>
FundInvariantVectorFamily<Tint>
ComputeFundamentalInvariant(MyMatrix<Tint> const &M,
                            [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: ComputeFundamentalInvariant, beginning\n";
  WriteMatrix(os, M);
#endif
  MyVector<Tint> eVect = SmithNormalFormInvariant(M);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: ComputeFundamentalInvariant, We have eVect\n";
#endif
  int dim = eVect.size();
  int rank = 0;
  Tint index = 1;
  for (int i = 0; i < dim; i++) {
    Tint val = eVect(i);
    if (val != 0) {
      rank++;
      index *= val;
    }
  }
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  int n = M.cols();
  if (rank != RankMat(M)) {
    std::cerr << "Something is inconsistent here\n";
    throw TerminalException{1};
  }
  if (rank == n) {
    Tint index_B = T_abs(Int_IndexLattice(M));
    if (index != index_B) {
      std::cerr << "index=" << index << " but index_B=" << index_B << "\n";
      throw TerminalException{1};
    }
  }
#endif
  return {rank, index};
}

template <typename Tint>
FundInvariantVectorFamily<Tint>
ComputeRankInvariant(MyMatrix<Tint> const &M,
                     [[maybe_unused]] std::ostream &os) {
  int rank = RankMat(M);
  Tint index(1);
  return {rank, index};
}

template <typename Tint>
bool operator>(FundInvariantVectorFamily<Tint> const &x,
               FundInvariantVectorFamily<Tint> const &y) {
  if (x.rank != y.rank) {
    return x.rank > y.rank;
  }
  if (x.index != y.index) {
    return x.index < y.index;
  }
  return false;
}

template <typename Tint>
bool operator<(FundInvariantVectorFamily<Tint> const &x,
               FundInvariantVectorFamily<Tint> const &y) {
  if (x.rank != y.rank) {
    return x.rank < y.rank;
  }
  if (x.index != y.index) {
    return x.index > y.index;
  }
  return false;
}

// Core search loop. Returns the family before antipodal duplication (one
// representative per +/-v pair); ExtractInvariantVectorFamily below
// duplicates it to restore the historical (fully duplicated) public
// return convention. Callers that only need one representative per pair
// (e.g. to feed conversion_and_duplication themselves later) should call
// this directly instead of duplicating then immediately halving again.
template <typename T, typename Tint, typename Ffinal, typename Finvariant>
MyMatrix<Tint> ExtractInvariantVectorFamilyHalf(MyMatrix<T> const &GramMat,
                                                Ffinal f_final,
                                                Finvariant f_invariant,
                                                std::ostream &os) {
  int n = GramMat.rows();
  T incr = GetSmallestIncrement(GramMat);
  T MaxNorm = GetMaxNorm<T, Tint>(GramMat, os);
  T norm = incr;
  MyMatrix<Tint> SHVret(0, n);
  FundInvariantVectorFamily<Tint> fi_ret = TrivFundamentalInvariant<Tint>();
  CVPSolver<T, Tint> solver(GramMat, os);
  while (true) {
    if (norm > MaxNorm) {
      std::cerr << "Failed to find a relevant vector configuration\n";
      throw TerminalException{1};
    }
    std::vector<MyVector<Tint>> ListVect = solver.fixed_norm_vectors(norm);
    // The shell is written straight into the concatenation: the intermediate
    // matrix of the new vectors was only ever fed to Concatenate.
    int n_prev = SHVret.rows();
    int n_new = ListVect.size();
    MyMatrix<Tint> SHVtest(n_prev + n_new, n);
    for (int i_row = 0; i_row < n_prev; i_row++) {
      for (int i_col = 0; i_col < n; i_col++) {
        SHVtest(i_row, i_col) = SHVret(i_row, i_col);
      }
    }
    for (int i_row = 0; i_row < n_new; i_row++) {
      for (int i_col = 0; i_col < n; i_col++) {
        SHVtest(n_prev + i_row, i_col) = ListVect[i_row](i_col);
      }
    }
    FundInvariantVectorFamily<Tint> fi_test = f_invariant(SHVtest);
    if (fi_ret < fi_test) {
      SHVret = SHVtest;
      fi_ret = fi_test;
      if (f_final(SHVret, fi_ret)) {
        return SHVret;
      }
    }
    norm += incr;
  }
}

template <typename T, typename Tint, typename Ffinal, typename Finvariant>
MyMatrix<Tint> ExtractInvariantVectorFamily(MyMatrix<T> const &GramMat,
                                            Ffinal f_final,
                                            Finvariant f_invariant,
                                            std::ostream &os) {
  MyMatrix<Tint> SHVhalf =
      ExtractInvariantVectorFamilyHalf<T, Tint, Ffinal, Finvariant>(
          GramMat, f_final, f_invariant, os);
  return matrix_duplication(SHVhalf);
}

// The half family, one vector per antipodal pair, before the duplication of
// ExtractInvariantVectorFamilyFullRank above; for callers that want the pairs.
template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyFullRankHalf(MyMatrix<T> const &eMat,
                                                        std::ostream &os) {
  int n = eMat.rows();
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M,
                     FundInvariantVectorFamily<Tint> const &fi) -> bool {
    return fi.rank == n;
  };
  auto f_invariant = [&](MyMatrix<Tint> const &M)
      -> FundInvariantVectorFamily<Tint> {
    return ComputeRankInvariant(M, os);
  };
  return ExtractInvariantVectorFamilyHalf<T, Tint, decltype(f_final),
                                          decltype(f_invariant)>(
      eMat, f_final, f_invariant, os);
}

// Representatives of the first K orbits of nonzero integer vectors under the
// group generated by ListGen (acting by v -> g v), each paired with its orbit
// size, ordered by increasing GramMat-norm v^T GramMat v. Vectors are enumerated
// shell by shell with the norm increasing by GetSmallestIncrement(GramMat), so
// there is no bound parameter. (Antipodes -v and v fall in the same orbit
// exactly when -Id is in the group, which holds e.g. for an arithmetic
// automorphism group.)
template <typename T, typename Tint>
std::vector<std::pair<MyVector<Tint>, size_t>>
get_k_short_orbit_vectors(MyMatrix<T> const &GramMat,
                          std::vector<MyMatrix<Tint>> const &ListGen, size_t K,
                          std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  T incr = GetSmallestIncrement(GramMat);
  T norm = incr;
  std::vector<std::pair<MyVector<Tint>, size_t>> reps;
  while (reps.size() < K) {
    std::vector<MyVector<Tint>> ListVect = solver.fixed_norm_vectors(norm);
    // The orbit of v under v -> g v stays within this norm shell (the ListGen
    // preserve the GramMat), so the dedup set only needs the current shell.
    // Vectors are sign-canonicalized, so antipodes -v and v count as one (the
    // orbit size is then the number of distinct lines, i.e. of distinct v v^T
    // directions).
    std::unordered_set<MyVector<Tint>> assigned;
    for (auto &v : ListVect) {
      if (reps.size() >= K) {
        break;
      }
      MyVector<Tint> cv = SignCanonicalizeVector(v);
      // insert().second is false when cv was already assigned to an earlier
      // orbit of this shell.
      if (!assigned.insert(cv).second) {
        continue;
      }
      std::vector<MyVector<Tint>> todo{cv};
      size_t orbit_size = 1;
      while (!todo.empty()) {
        MyVector<Tint> w = todo.back();
        todo.pop_back();
        for (auto &g : ListGen) {
          MyVector<Tint> gw = g * w;
          MyVector<Tint> cgw = SignCanonicalizeVector(gw);
          if (assigned.insert(cgw).second) {
            todo.push_back(cgw);
            orbit_size++;
          }
        }
      }
      reps.push_back({v, orbit_size});
    }
    norm += incr;
  }
  return reps;
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyFullRank(MyMatrix<T> const &eMat,
                                                    std::ostream &os) {
  int n = eMat.rows();
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M, FundInvariantVectorFamily<Tint> const& fi) -> bool {
    return fi.rank == n;
  };
  auto f_invariant=[&](MyMatrix<Tint> const &M) -> FundInvariantVectorFamily<Tint> {
    return ComputeRankInvariant(M, os);
  };
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_final), decltype(f_invariant)>(eMat, f_final, f_invariant, os);
}

// Half version (one representative per +/-v pair, no antipodal
// duplication) of ExtractInvariantVectorFamilyZbasis below, for callers
// that would otherwise immediately undo the duplication themselves.
template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasisHalf(MyMatrix<T> const &eMat,
                                                      std::ostream &os) {
  int n = eMat.cols();
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M, FundInvariantVectorFamily<Tint> const& fi) -> bool {
    return IsCompleteSystem(fi, n);
  };
  auto f_invariant=[&](MyMatrix<Tint> const &M) -> FundInvariantVectorFamily<Tint> {
    return ComputeFundamentalInvariant(M, os);
  };
  return ExtractInvariantVectorFamilyHalf<T, Tint, decltype(f_final), decltype(f_invariant)>(eMat, f_final, f_invariant, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasis(MyMatrix<T> const &eMat,
                                                  std::ostream &os) {
  return matrix_duplication(ExtractInvariantVectorFamilyZbasisHalf<T, Tint>(eMat, os));
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantBreakingVectorFamily(
    MyMatrix<T> const &eMat, std::vector<MyMatrix<Tint>> const &ListMatr,
    std::ostream &os) {
  auto f_check = [&](MyMatrix<Tint> const &M) -> bool {
    int n_row = M.rows();
    std::unordered_set<MyVector<Tint>> set;
    for (int i = 0; i < n_row; i++) {
      MyVector<Tint> V = GetMatrixRow(M, i);
      set.insert(V);
    }
    for (auto &eMatr : ListMatr) {
      MyMatrix<Tint> Mprod = M * eMatr;
      for (int i = 0; i < n_row; i++) {
        MyVector<Tint> V = GetMatrixRow(Mprod, i);
        if (!set.contains(V))
          return true;
      }
    }
    return false;
  };
  auto f_invariant = [&](MyMatrix<Tint> const &M) -> FundInvariantVectorFamily<Tint> {
    int artificial_rank = 0;
    if (f_check(M)) {
      artificial_rank = 1;
    }
    Tint index(1);
    return {artificial_rank, index};
  };
  auto f_final = [&]([[maybe_unused]] MyMatrix<Tint> const &M, FundInvariantVectorFamily<Tint> const& fi) -> bool {
    return fi.rank == 1;
  };
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_final), decltype(f_invariant)>(eMat, f_final, f_invariant, os);
}

template <typename Tint> bool CheckCentralSymmetry(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  std::unordered_map<MyVector<Tint>, int> map;
  for (int i = 0; i < nbRow; i++) {
    MyVector<Tint> V = GetMatrixRow(M, i);
    if (!IsZeroVector(V)) {
      MyVector<Tint> Vcan = SignCanonicalizeVector(V);
      map[Vcan]++;
    }
  }
  for (auto &[vcan, multiplicity] : map) {
    if (multiplicity != 2)
      return false;
  }
  return true;
}

template <typename T, typename Tint>
MyMatrix<Tint> ComputeVoronoiRelevantVector(MyMatrix<T> const &GramMat,
                                            std::ostream &os) {
  int n = GramMat.rows();
  std::vector<MyVector<Tint>> ListVect;
  BlockCppIterator blk(n, 2);
  CVPSolver<T, Tint> solver(GramMat, os);
  for (auto &eVect : blk) {
    int sum = 0;
    for (auto &eVal : eVect) {
      sum += eVal;
    }
    if (sum > 0) {
      MyVector<T> eV(n);
      for (int u = 0; u < n; u++) {
        T val = UniversalScalarConversion<T, int>(eVect[u]);
        eV(u) = val / T(2);
      }
      resultCVP<T, Tint> result = solver.nearest_vectors(eV);
      if (result.ListVect.rows() == 2) {
        MyVector<Tint> Vins(n);
        for (int u = 0; u < n; u++) {
          Vins(u) = result.ListVect(0, u) - result.ListVect(1, u);
        }
        ListVect.push_back(Vins);
        ListVect.push_back(-Vins);
      }
    }
  }
  return MatrixFromVectorFamily(ListVect);
}

template <typename T, typename Tint>
MyMatrix<Tint> FilterByNorm(MyMatrix<T> const &GramMat,
                            MyMatrix<Tint> const &ListVect, std::ostream &os) {
  int n = GramMat.rows();
  std::map<T, std::vector<MyVector<Tint>>> map;
  std::vector<T> LineMat = GetLineVector(GramMat);
  int n_vect = ListVect.rows();
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    MyVector<Tint> V = GetMatrixRow(ListVect, i_vect);
    MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
    T norm = EvaluateLineVector(LineMat, V_T);
    map[norm].push_back(V);
  }
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: FilterByNorm, map built\n";
#endif
  MyMatrix<Tint> SHV_ret(0, n);
  FundInvariantVectorFamily<Tint> fi_ret = TrivFundamentalInvariant<Tint>();
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  size_t pos = 0;
#endif
  for (auto &[norm, vectors] : map) {
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: FilterByNorm, pos=" << pos
       << " |vectors|=" << vectors.size() << "\n";
#endif
    MyMatrix<Tint> BlkMat = MatrixFromVectorFamily(vectors);
    MyMatrix<Tint> SHV_new = Concatenate(SHV_ret, BlkMat);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: FilterByNorm, We have SHV_new\n";
#endif
    FundInvariantVectorFamily<Tint> fi_new =
        ComputeFundamentalInvariant(SHV_new, os);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: FilterByNorm, We have fi_new\n";
#endif
    if (fi_new > fi_ret) {
      SHV_ret = SHV_new;
      fi_ret = fi_new;
      if (IsCompleteSystem(fi_ret, n)) {
        return SHV_ret;
      }
    }
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    pos += 1;
#endif
  }
  std::cerr
      << "Failed to terminate and find a correct family for the filtration\n";
  std::cerr
      << "This indicates that the original family was not a spanning one\n";
  throw TerminalException{1};
}

/*
  The characteristic vector set V_cv of Section 2.2 of the paper
  "A canonical form for positive definite matrices" (Dutour Sikiric,
  Haensch, Voight, van Woerden).

  A characteristic vector set function (Definition 1.2.1) is a map A -> V(A)
  from positive definite n x n matrices to finite subsets of Z^n such that
  --- V(A) generates Z^n as a group.
  --- U^{-1} V(A) = V(U^T A U) for every U in GL_n(Z).
  The second property is what makes it usable for canonicalization: the
  canonical form of A can be computed from V(A) alone.

  The construction of V_cv is inductive.  Write Min(A) for the set of
  vectors realizing the minimum of A and cvd(A,v) = min_{x in Z^n} A[x - v]
  for the squared distance from v to the lattice, and
      CV(A, v) = { x in Z^n : A[x - v] = cvd(A, v) }
  for the set of closest vectors.  Both are canonical: CV(U^T A U, U^{-1} v)
  = U^{-1} CV(A, v), which is what makes the whole construction work.

  Step 1: the well-rounded case (2.2.5).  If A is well-rounded, that is if
  Min(A) spans R^n, let L_min be the lattice generated by Min(A) and B a
  Z-basis of it.  Then
      V_wr-cv(A) = Min(A) union
                   union_{c in Z^n / L_min} (c - B CV(B^T A B, B^{-1} c))
  The deep-hole-like vectors c - B x pick up the part of Z^n that Min(A)
  misses, so V_wr-cv(A) generates Z^n even when Min(A) does not.

  Step 2: the general case (2.2.6).  Let L_1 = satspan(Min(A)) = QMin(A) cap
  Z^n with Z-basis B_1 (r x n), and A_1 = B_1 A B_1^T.  By construction A_1
  is well-rounded.  Let proj be the orthogonal projection onto the A-
  orthogonal complement of L_1, L_2 = proj(Z^n) with basis B_2
  ((n-r) x n, rational), and A_2 = B_2 A B_2^T.  Then
      V_cv(A) = B_1 V_wr-cv(A_1) union
                union_{v in B_2 V_cv(A_2)} CV(A, v)
  with V_cv of the empty matrix set to the empty set.  Since Z^n / L_1 is
  isomorphic to L_2 by proj, the closest vectors to the lifts of a
  generating set of L_2 complete B_1 V_wr-cv(A_1) to a generating set of
  Z^n.

  Both functions take a `spanning` flag. With spanning = true the family
  generates Z^n, as Definition 1.2.1 asks. With spanning = false the coset
  unions of (2.2.5) are dropped at every level of the recursion and the
  family is only of full rank, generating a finite index subgroup of Z^n.
  The second is usually much smaller and is the better choice whenever the
  consumer can work from a full rank family, since those coset unions are
  the expensive part of the construction: one closest-vector computation per
  coset, and there are [Z^n : L_min] of them at each level. What remains is
  still canonical, the decision to stop at full rank being uniform rather
  than read off the data.

  Two deviations from the text of the paper, both needed for correctness:
  --- The coset c = 0 in (2.2.5) contributes 0 - B CV(A_1, 0) = {0}.  The
      zero vector is dropped, it carries no information and a vector family
      containing it is not usable downstream.
  --- CV(A, v) in (2.2.6) is not guaranteed to consist of vectors x with
      proj(x) = v: a lattice point in a neighbouring L_1-coset can be
      closer to v than every point of the coset of v.  When that happens
      the union misses the coset of v altogether and V_cv(A) does not
      generate Z^n, so (2.2.6) as written is not a characteristic vector
      set function.  A family of counterexamples: for k >= 3 let A(k, N) be
      the Gram matrix of e_1, ..., e_k, f with the e_i pairwise orthogonal
      of norm 4, <f, e_i> = 2 and <f, f> = N.  Then Min(A) = {+- e_i} spans
      only L_1 = <e_1, ..., e_k>, the covering radius squared of L_1 is k,
      and proj(f) has norm N - k.  For N - k < k the origin is strictly
      closer to v = proj(f) than every point of f + L_1, so CV(A, v) = {0}
      and the formula of the paper returns a family generating a subgroup
      of index 2; A(3,5), A(8,11) and A(8,15) are worked out in
      results/charvect_cv.md of the strongly perfect lattices notes, where
      the measurement is recorded.  We therefore add, for each v, the
      closest vectors *within the coset determined by v*, which is the set
      the generation argument actually needs.  This is canonical for the
      same reason CV is, the coset being read off from v.  Both sets are
      kept so that no information is lost.

  Everything is in the row convention: vectors are rows, a basis matrix B
  has the basis vectors as rows, and the Gram matrix of the sublattice
  spanned by the rows of B is B A B^T.  A row vector y in the B-coordinates
  corresponds to the ambient vector y B, that is B^T y in column form.
 */

// Min(A), with both signs, as rows of the returned matrix.
template <typename T, typename Tint>
MyMatrix<Tint> CharVectSet_MinimalVectors(MyMatrix<T> const &GramMat,
                                          std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  Tshortest<T, Tint> RecSHV = solver.shortest_vectors();
  // Already antipodal, unlike fixed_norm_vectors which returns one vector per
  // pair and is the one the shell based families have to duplicate.
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  if (!CheckCentralSymmetry(RecSHV.SHV)) {
    std::cerr << "IVF: shortest_vectors did not return an antipodal family\n";
    throw TerminalException{1};
  }
#endif
  return RecSHV.SHV;
}

// Insertion into a list of rows with duplicate removal, the sets built below
// overlap a lot (CV(A,v) for nearby v, Min(A) and the c = 0 coset, ...).
template <typename Tint> struct CharVectSet_Accumulator {
  std::set<MyVector<Tint>> TheSet;
  // When set, the sign of every vector is normalized on its first nonzero
  // coordinate before insertion, so the set holds one vector per antipodal
  // pair and never the pair itself. That is what the consumers want: the
  // weight matrix and the graph are quadratic in the number of vectors, so
  // building the family already halved is worth more than halving it later.
  bool antipodal_half = false;
  void insert(MyVector<Tint> const &V) {
    if (IsZeroVector(V)) {
      return;
    }
    if (antipodal_half) {
      int n = V.size();
      for (int i = 0; i < n; i++) {
        if (V(i) != 0) {
          if (V(i) < 0) {
            TheSet.insert(MyVector<Tint>(-V));
            return;
          }
          break;
        }
      }
    }
    TheSet.insert(V);
  }
  void insert_rows(MyMatrix<Tint> const &M) {
    for (int i_row = 0; i_row < M.rows(); i_row++) {
      insert(GetMatrixRow(M, i_row));
    }
  }
  MyMatrix<Tint> get_matrix(int const &n) const {
    MyMatrix<Tint> M(TheSet.size(), n);
    int pos = 0;
    for (auto &V : TheSet) {
      for (int i = 0; i < n; i++) {
        M(pos, i) = V(i);
      }
      pos++;
    }
    return M;
  }
};

/*
  V_wr-cv(A) of (2.2.5), for a well-rounded A.

  With spanning = false only Min(A) is returned. That is enough to be of full
  rank, A being well-rounded, and it drops the union over Z^n / L_min, which
  is the only expensive part: it costs one closest-vector computation per
  coset and there are [Z^n : L_min] of them. What is lost is that the family
  then generates L_min instead of Z^n.
 */
template <typename T, typename Tint>
MyMatrix<Tint>
CharacteristicVectorSetWellRoundedCV(MyMatrix<T> const &GramMat,
                                     bool const &spanning,
                                     bool const &antipodal_half,
                                     std::ostream &os) {
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  MicrosecondTime time;
#endif
  int n = GramMat.rows();
  MyMatrix<Tint> SHV = CharVectSet_MinimalVectors<T, Tint>(GramMat, os);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  if (RankMat(SHV) != n) {
    std::cerr << "IVF: CharacteristicVectorSetWellRoundedCV called on a "
              << "matrix that is not well-rounded\n";
    throw TerminalException{1};
  }
#endif
  CharVectSet_Accumulator<Tint> acc;
  acc.antipodal_half = antipodal_half;
  acc.insert_rows(SHV);
  if (!spanning) {
    // Min(A) alone, which is of full rank since A is well-rounded. The union
    // over Z^n / L_min below is what makes the family span Z^n and it is also
    // the whole cost of this function.
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
    os << "|IVF: CharacteristicVectorSetWellRoundedCV(fullrank)|=" << time
       << "\n";
#endif
    return acc.get_matrix(n);
  }
  // B is a Z-basis of the lattice L_min generated by Min(A).
  MyMatrix<Tint> B = GetZbasis(SHV);
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  MyMatrix<T> GramB = B_T * GramMat * B_T.transpose();
  // c = c_B B in row convention, that is c = B^T c_B in column form.
  MyMatrix<T> B_Ttr = B_T.transpose();
  MyMatrix<T> BinvTr = Inverse(B_Ttr);
  std::vector<MyVector<Tint>> ListClasses =
      ComputeTranslationClasses<Tint, Tint>(B);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: WellRoundedCV n=" << n << " |Min|=" << SHV.rows()
     << " |Z^n/L_min|=" << ListClasses.size() << "\n";
#endif
  CVPSolver<T, Tint> solver(GramB, os);
  for (auto &eClass : ListClasses) {
    MyVector<T> eClass_T = UniversalVectorConversion<T, Tint>(eClass);
    MyVector<T> eClassB = BinvTr * eClass_T;
    resultCVP<T, Tint> res = solver.nearest_vectors(eClassB);
    for (int i_row = 0; i_row < res.ListVect.rows(); i_row++) {
      MyVector<Tint> y = GetMatrixRow(res.ListVect, i_row);
      // The vector is c - y B.
      MyVector<Tint> eVect = eClass - B.transpose() * y;
      acc.insert(eVect);
    }
  }
  MyMatrix<Tint> RetMat = acc.get_matrix(n);
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  os << "|IVF: CharacteristicVectorSetWellRoundedCV|=" << time << "\n";
#endif
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: WellRoundedCV |V_wr-cv|=" << RetMat.rows() << "\n";
#endif
  return RetMat;
}

/*
  V_cv(A) of (2.2.6).  Recursive over the successive saturated spans of the
  minimal vectors.
 */
template <typename T, typename Tint>
MyMatrix<Tint> CharacteristicVectorSetCV(MyMatrix<T> const &GramMat,
                                         bool const &spanning,
                                         bool const &antipodal_half,
                                         std::ostream &os) {
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  MicrosecondTime time;
#endif
  int n = GramMat.rows();
  MyMatrix<Tint> SHV = CharVectSet_MinimalVectors<T, Tint>(GramMat, os);
  int r = RankMat(SHV);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: CV n=" << n << " |Min|=" << SHV.rows() << " r=" << r << "\n";
#endif
  if (r == n) {
    // A is well-rounded, the recursion stops.
    return CharacteristicVectorSetWellRoundedCV<T, Tint>(GramMat, spanning,
                                                         antipodal_half, os);
  }
  CharVectSet_Accumulator<Tint> acc;
  acc.antipodal_half = antipodal_half;
  // B_1 is a Z-basis of L_1 = satspan(Min(A)), r x n.
  MyMatrix<Tint> B1 = IntegralSpaceSaturation(GetZbasis(SHV));
  MyMatrix<T> B1_T = UniversalMatrixConversion<T, Tint>(B1);
  MyMatrix<T> Gram1 = B1_T * GramMat * B1_T.transpose();
  MyMatrix<Tint> Vwr = CharacteristicVectorSetWellRoundedCV<T, Tint>(
      Gram1, spanning, antipodal_half, os);
  // B_1 V_wr-cv(A_1), that is the rows y of Vwr sent to y B_1.
  acc.insert_rows(Vwr * B1);
  // C completes B_1 to a Z-basis of Z^n, and B_2 = proj(C) is a basis of
  // L_2 = proj(Z^n) since proj has kernel L_1 on Z^n, L_1 being saturated.
  MyMatrix<Tint> C = SubspaceCompletionInt(B1, n);
  MyMatrix<T> C_T = UniversalMatrixConversion<T, Tint>(C);
  MyMatrix<T> CrossGram = C_T * GramMat * B1_T.transpose();
  MyMatrix<T> Gram1inv = Inverse(Gram1);
  MyMatrix<T> Coeff = CrossGram * Gram1inv;
  MyMatrix<T> B2 = C_T - Coeff * B1_T;
  MyMatrix<T> Gram2 = B2 * GramMat * B2.transpose();
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  MyMatrix<T> Ortho = B2 * GramMat * B1_T.transpose();
  if (!IsZeroMatrix(Ortho)) {
    std::cerr << "IVF: the projection B2 is not A-orthogonal to L_1\n";
    throw TerminalException{1};
  }
#endif
  // The recursion is on a rational Gram matrix.  Since V_cv only depends on
  // argmin sets it is invariant under scaling of the Gram matrix, but the
  // rational entries are kept as they are, no rescaling is needed for
  // correctness.
  /*
    The recursion may return one vector per pair. That loses nothing: the
    closest vectors to -v are the negatives of those to v, so the two
    contribute the same set once the signs are normalized.
   */
  MyMatrix<Tint> V2 =
      CharacteristicVectorSetCV<T, Tint>(Gram2, spanning, antipodal_half, os);
  CVPSolver<T, Tint> solver(GramMat, os);
  CVPSolver<T, Tint> solver1(Gram1, os);
  // Projection of Z^n onto L_2 written in the B_2 coordinates, that is
  // x |-> x - (x A B_1^T) Gram1^{-1} B_1 followed by the B_2 coordinates.
  // In the basis (B_1, C) of Z^n the projection of B_1 is 0 and that of C
  // is B_2, so the B_2 coordinate of proj(x) is the C part of x.
  MyMatrix<Tint> FullBasis(n, n);
  for (int i = 0; i < r; i++) {
    FullBasis.row(i) = B1.row(i);
  }
  for (int i = 0; i < n - r; i++) {
    FullBasis.row(r + i) = C.row(i);
  }
  for (int i_row = 0; i_row < V2.rows(); i_row++) {
    MyVector<Tint> y = GetMatrixRow(V2, i_row);
    MyVector<T> y_T = UniversalVectorConversion<T, Tint>(y);
    // v = y B_2 in R^n.
    MyVector<T> v = B2.transpose() * y_T;
    resultCVP<T, Tint> res = solver.nearest_vectors(v);
    acc.insert_rows(res.ListVect);
    // The closest vectors of the L_1-coset of v itself.  See the header
    // comment: CV(A, v) alone can miss that coset.
    MyVector<Tint> x0(n);
    for (int i = 0; i < r; i++) {
      x0(i) = 0;
    }
    for (int i = 0; i < n - r; i++) {
      x0(r + i) = y(i);
    }
    MyVector<Tint> xLift = FullBasis.transpose() * x0;
    // The coset is xLift + L_1, so we look for the point of L_1 closest to
    // v - xLift, in the B_1 coordinates.
    MyVector<T> xLift_T = UniversalVectorConversion<T, Tint>(xLift);
    MyVector<T> resid = v - xLift_T;
    // The B_1 coordinates of the L_1 component of resid.
    MyVector<T> residB1 = Gram1inv * (B1_T * (GramMat * resid));
    resultCVP<T, Tint> res1 = solver1.nearest_vectors(residB1);
    for (int i_r = 0; i_r < res1.ListVect.rows(); i_r++) {
      MyVector<Tint> z = GetMatrixRow(res1.ListVect, i_r);
      acc.insert(xLift + B1.transpose() * z);
    }
  }
  MyMatrix<Tint> RetMat = acc.get_matrix(n);
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  os << "|IVF: CharacteristicVectorSetCV|=" << time << "\n";
#endif
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: CV |V_cv|=" << RetMat.rows() << "\n";
#endif
  return RetMat;
}

/*
  The invariant vector family a canonicalization works from, together with
  whether it generates Z^n.

  Which family is best is not decidable in advance. The shell based
  ExtractInvariantVectorFamily is cheap when the successive minima are close
  together and catastrophic when they are not, a whole shell having to be
  taken at once: 16942 vectors on TestData/SlowCanonic/slow_canonic_2, 16816
  of them in the single shell of norm 14. V_cv does not care how far apart
  the minima are but pays 2^k closest-vector computations when the minimal
  vectors span a sublattice of index 2^k. Neither dominates: on
  slow_canonic_1 the shells give 2358 vectors against 8516 for the spanning
  V_cv, on slow_canonic_2 it is 16942 against 238.

  Building a family costs a fraction of a percent of what is done with it,
  the weight matrix being quadratic in its size and the canonical labelling
  worse than quadratic, so the cheapest way to avoid the bad case of either
  is to build both and keep the smaller. Both sizes are invariants of the
  isometry class, so the choice is itself canonical and the canonical forms
  obtained from it remain comparable.
 */
template <typename Tint> struct CanonicVectorFamily {
  /*
    One vector per antipodal pair. The families are closed under v -> -v, and
    everything downstream is quadratic in the number of vectors, so they are
    built halved rather than built whole and halved afterwards. Consumers
    that need both signs, the ones computing the automorphism group of the
    configuration as a permutation group, call get_full.
   */
  MyMatrix<Tint> SHVhalf;
  // Measured by FamilySpansLattice, never inferred from how the family was
  // built: a family built without asking to span usually spans anyway, and
  // every family that does can take the cheap canonicalization. The half
  // family spans what the whole one does, -v being in the span of v.
  bool spans_lattice;
  /*
    The whole family, the chosen representatives first and their negatives
    after, in the same order. That is the indexing the generators lifted by
    AbsTrick_LiftGenerators use, [0, nbPair) for +v and [nbPair, 2 nbPair)
    for -v, so the two agree without a translation table.
   */
  MyMatrix<Tint> get_full() const {
    int nbPair = SHVhalf.rows();
    int n = SHVhalf.cols();
    MyMatrix<Tint> SHV(2 * nbPair, n);
    for (int iPair = 0; iPair < nbPair; iPair++) {
      for (int i = 0; i < n; i++) {
        SHV(iPair, i) = SHVhalf(iPair, i);
        SHV(nbPair + iPair, i) = -SHVhalf(iPair, i);
      }
    }
    return SHV;
  }
  int n_pair() const { return SHVhalf.rows(); }
};

/*
  Whether the rows of SHV generate Z^n rather than a finite index subgroup.
  A Z-basis of the span followed by its determinant, which costs nothing
  next to what the answer saves.
 */
template <typename Tint>
bool FamilySpansLattice(MyMatrix<Tint> const &SHV) {
  int n = SHV.cols();
  if (RankMat(SHV) != n) {
    return false;
  }
  MyMatrix<Tint> Basis = GetZbasis(SHV);
  Tint det = DeterminantMat(Basis);
  return T_abs(det) == Tint(1);
}

template <typename Tint>
CanonicVectorFamily<Tint> get_canonic_vector_family(MyMatrix<Tint> &&SHVhalf) {
  bool spans = FamilySpansLattice<Tint>(SHVhalf);
  return {std::move(SHVhalf), spans};
}

/*
  The smaller of the two families, both stopped at full rank. Spanning Z^n is
  not required of them: the canonicalization only has to place Z^n relative
  to the span of the family, under the group preserving it, and a family
  built to span can be far larger, 8516 against 324 vectors on
  slow_canonic_1.
 */
template <typename T, typename Tint>
CanonicVectorFamily<Tint> GetCanonicVectorFamily(MyMatrix<T> const &GramMat,
                                                 std::ostream &os) {
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  MicrosecondTime time;
#endif
  const bool antipodal_half = true;
  MyMatrix<Tint> SHV_shell =
      ExtractInvariantVectorFamilyFullRankHalf<T, Tint>(GramMat, os);
  MyMatrix<Tint> SHV_cv =
      CharacteristicVectorSetCV<T, Tint>(GramMat, false, antipodal_half, os);
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
  os << "IVF: family choice, shells " << SHV_shell.rows() << " against V_cv "
     << SHV_cv.rows() << " (pairs)\n";
#endif
  auto f_ret = [&]() -> CanonicVectorFamily<Tint> {
    if (SHV_cv.rows() <= SHV_shell.rows()) {
      return get_canonic_vector_family<Tint>(std::move(SHV_cv));
    }
    return get_canonic_vector_family<Tint>(std::move(SHV_shell));
  };
  CanonicVectorFamily<Tint> fam = f_ret();
#ifdef TIMINGS_INVARIANT_VECTOR_FAMILY
  os << "|IVF: GetCanonicVectorFamily|=" << time << "\n";
#endif
  return fam;
}

// clang-format off
#endif  // SRC_LATT_INVARIANTVECTORFAMILY_H_
// clang-format on
