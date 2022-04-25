#ifndef SRC_LATT_INVARIANTVECTORFAMILY_H_
#define SRC_LATT_INVARIANTVECTORFAMILY_H_

#include "ShortestUniversal.h"
#include "Temp_Positivity.h"
#include <set>
#include <vector>

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
MyMatrix<Tint> EnumerateVectorsFixedNorm(MyMatrix<T> const &eMat,
                                         T const &norm) {
#ifdef TIMINGS
  std::cerr << "Begining of ExtractInvariantVectorFamily\n";
  SingletonTime time1;
#endif
  LLLreduction<T, Tint> recLLL = LLLreducedBasis<T, Tint>(eMat);
  MyMatrix<Tint> P_T = recLLL.Pmat.transpose();
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|LLL|=" << ms(time1, time2) << "\n";
#endif
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(recLLL.Pmat);
  //
  //  std::cerr << "eMat=\n";
  //  WriteMatrix(std::cerr, eMat);
  MyMatrix<T> eMatRed = Pmat_T * eMat * TransposedMat(Pmat_T);
  MyMatrix<Tint> SHVall = T_ShortVector_fixed<T, Tint>(eMatRed, norm);
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|T_ShortVector|=" << ms(time2, time3) << "\n";
#endif
  MyMatrix<Tint> SHV1_a = SHVall * recLLL.Pmat;
  MyMatrix<Tint> SHV1_b = -SHV1_a;
  return Concatenate(SHV1_a, SHV1_b);
}

template <typename T, typename Tint> T GetMaxNorm(MyMatrix<T> const &eMat) {
  LLLreduction<T, Tint> recLLL = LLLreducedBasis<T, Tint>(eMat);
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
std::set<T> GetSetNormConsider(MyMatrix<T> const &eMat) {
  T incr = GetSmallestIncrement(eMat);
  T MaxNorm = GetMaxNorm<T, Tint>(eMat);
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

template <typename T, typename Tint, typename Fcorrect>
MyMatrix<Tint> ExtractInvariantVectorFamily(MyMatrix<T> const &eMat,
                                            Fcorrect f_correct) {
  int n = eMat.rows();
  T incr = GetSmallestIncrement(eMat);
  T MaxNorm = GetMaxNorm<T, Tint>(eMat);
  T norm = incr;
  MyMatrix<Tint> SHVret(0, n);
  while (true) {
    if (norm > MaxNorm) {
      std::cerr << "Failed to find a relevant vector configuration\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> SHV_f = EnumerateVectorsFixedNorm<T, Tint>(eMat, norm);
    SHVret = Concatenate(SHVret, SHV_f);
    if (f_correct(SHVret))
      return SHVret;
    norm += incr;
  }
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyFullRank(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  auto f_correct = [&](MyMatrix<Tint> const &M) -> bool {
    if (RankMat(M) == n)
      return true;
    return false;
  };
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_correct)>(eMat,
                                                                    f_correct);
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasis(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  auto f_correct = [&](MyMatrix<Tint> const &M) -> bool {
    if (RankMat(M) < n)
      return false;
    //    std::cerr << "Before Int_IndexLattice computation\n";
    Tint indx = Int_IndexLattice(M);
    //    std::cerr << "indx=" << indx << "\n";
    if (T_Norm(indx) == 1)
      return true;
    return false;
  };
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_correct)>(eMat,
                                                                    f_correct);
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantBreakingVectorFamily(
    MyMatrix<T> const &eMat, std::vector<MyMatrix<Tint>> const &ListMatr) {
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
  return ExtractInvariantVectorFamily<T, Tint, decltype(f_correct)>(eMat,
                                                                    f_correct);
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

#endif //  SRC_LATT_INVARIANTVECTORFAMILY_H_
