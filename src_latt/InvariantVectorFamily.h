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
  std::vector<T> ListVal;
  int n = eMat.rows();
  T eGcd = eMat(0, 0);
  for (int i = 1; i < n; i++)
    eGcd = GcdPair(eGcd, eMat(i, i));
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) {
      T val = 2 * eMat(i, j);
      eGcd = GcdPair(eGcd, val);
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
    if (set.count(V) == 0) {
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
    if (set.count(V) == 1) {
      return true;
    }
    set.insert(V);
  }
  return false;
}

template <typename T> void check_antipodality_mymatrix(MyMatrix<T> const &SHV) {
  if (!is_antipodal(SHV)) {
    std::cerr << "TSPACE: The family SHV is not antipodal\n";
    throw TerminalException{1};
  }
}

template <typename Tint> bool IsFullDimZbasis(MyMatrix<Tint> const &M) {
  int n = M.cols();
  if (RankMat(M) < n) {
    return false;
  }
  Tint indx = Int_IndexLattice(M);
  if (T_NormGen(indx) == 1) {
    return true;
  }
  return false;
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

template <typename T, typename Tint, typename Fcorrect>
MyMatrix<Tint> ExtractInvariantVectorFamily(MyMatrix<T> const &eMat,
                                            Fcorrect f_correct,
                                            std::ostream &os) {
  int n = eMat.rows();
  T incr = GetSmallestIncrement(eMat);
  T MaxNorm = GetMaxNorm<T, Tint>(eMat, os);
  T norm = incr;
  MyMatrix<Tint> SHVret(0, n);
  while (true) {
    if (norm > MaxNorm) {
      std::cerr << "Failed to find a relevant vector configuration\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> SHV_f = T_ShortVector_fixed<T, Tint>(eMat, norm, os);
    SHVret = Concatenate(SHVret, SHV_f);
    if (f_correct(SHVret))
      return matrix_duplication(SHVret);
    norm += incr;
  }
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyFullRank(MyMatrix<T> const &eMat,
                                                    std::ostream &os) {
  int n = eMat.rows();
  auto f_correct = [&](MyMatrix<Tint> const &M) -> bool {
    if (RankMat(M) == n)
      return true;
    return false;
  };
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_correct)>(
      eMat, f_correct, os);
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasis(MyMatrix<T> const &eMat,
                                                  std::ostream &os) {
  auto f_correct = [&](MyMatrix<Tint> const &M) -> bool {
    return IsFullDimZbasis(M);
  };
  MyMatrix<Tint> SHV =
      ExtractInvariantVectorFamily<T, Tint, decltype(f_correct)>(eMat,
                                                                 f_correct, os);
#ifdef SANITY_CHECK_INVARIANT_VECTOR_FAMILY
  check_antipodality_mymatrix(SHV);
#endif
  return SHV;
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantBreakingVectorFamily(
    MyMatrix<T> const &eMat, std::vector<MyMatrix<Tint>> const &ListMatr,
    std::ostream &os) {
  auto f_correct = [&](MyMatrix<Tint> const &M) -> bool {
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
        if (set.count(V) == 0)
          return true;
      }
    }
    return false;
  };
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_correct)>(
      eMat, f_correct, os);
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
  for (auto &kv : map) {
    if (kv.second != 2)
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
        eV(u) = val / 2;
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
  for (auto &kv : map) {
#ifdef DEBUG_INVARIANT_VECTOR_FAMILY
    os << "IVF: FilterByNorm, pos=" << pos
       << " |kv.second|=" << kv.second.size() << "\n";
#endif
    MyMatrix<Tint> BlkMat = MatrixFromVectorFamily(kv.second);
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
  See the paper "A canonical form for positive definite matrices"
  It has to be coded
 */
/*
template<typename T, typename Tint>
MyMatrix<Tint> WellRoundedInvariantVectorFamily(MyMatrix<Tint> const& GramMat) {

}
*/

// clang-format off
#endif  // SRC_LATT_INVARIANTVECTORFAMILY_H_
// clang-format on
