#ifndef SRC_LATT_INVARIANTVECTORFAMILY_H_
#define SRC_LATT_INVARIANTVECTORFAMILY_H_

#include "ShortestUniversal.h"
#include "Temp_Positivity.h"
#include <set>
#include <vector>

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamily(
    MyMatrix<T> const &eMat,
    std::function<bool(MyMatrix<Tint> const &)> const &fCorrect) {
  T MaxNorm = MaximumDiagonal(eMat);
#ifdef TIMINGS
  std::cerr << "Begining of ExtractInvariantVectorFamily\n";
  SingletonTime time1;
#endif
  MyMatrix<Tint> SHVall = T_ShortVector<T, Tint>(eMat, MaxNorm);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|T_ShortVector|=" << ms(time1,time2) << "\n";
#endif
  //  std::cerr << "MaxNorm = " << MaxNorm << " eMat =\n";
  //  WriteMatrix(std::cerr, eMat);
  //  std::cerr << "|SHVall|=" << SHVall.rows() << "\n";
  //  WriteMatrix(std::cerr, SHVall);
  std::set<T> SetNorm;
  int nbSHV = SHVall.rows();
  std::vector<T> ListNorm(nbSHV);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<Tint> eRow = GetMatrixRow(SHVall, iSHV);
    T eNorm = EvaluationQuadForm(eMat, eRow);
    ListNorm[iSHV] = eNorm;
    if (eNorm > 0)
      SetNorm.insert(eNorm);
  }
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|Preparation|=" << ms(time2,time3) << "\n";
#endif
  std::vector<MyVector<Tint>> ListVect;
  for (auto const &eNorm : SetNorm) {
    //    std::cerr << "eNorm=" << eNorm << "\n";
    for (int iSHV = 0; iSHV < nbSHV; iSHV++)
      if (ListNorm[iSHV] == eNorm) {
        MyVector<Tint> eRow = GetMatrixRow(SHVall, iSHV);
        ListVect.push_back(eRow);
        ListVect.push_back(-eRow);
      }
    MyMatrix<Tint> SHVret = MatrixFromVectorFamily(ListVect);
    if (fCorrect(SHVret))
      return SHVret;
  }
  std::cerr
      << "ExtractInvariantVectorFamily : We should never reach that stage\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyFullRank(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  std::function<bool(MyMatrix<Tint> const &)> fCorrect =
      [&](MyMatrix<Tint> const &M) -> bool {
    if (RankMat(M) == n)
      return true;
    return false;
  };
  return ExtractInvariantVectorFamily(eMat, fCorrect);
}

template <typename T, typename Tint>
MyMatrix<Tint>
ExtractInvariantVectorFamilyZbasis_Kernel(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  std::function<bool(MyMatrix<Tint> const &)> fCorrect =
      [&](MyMatrix<Tint> const &M) -> bool {
    if (RankMat(M) < n)
      return false;
    //    std::cerr << "Before Int_IndexLattice computation\n";
    Tint indx = Int_IndexLattice(M);
    //    std::cerr << "indx=" << indx << "\n";
    if (T_Norm(indx) == 1)
      return true;
    return false;
  };
  return ExtractInvariantVectorFamily(eMat, fCorrect);
}

template <typename Tint> bool CheckCentralSymmetry(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  int n = M.cols();
  int res = nbRow % 2;
  if (res == 1)
    return false;
  int nbPair = nbRow / 2;
  for (int iPair = 0; iPair < nbPair; iPair++) {
    for (int i = 0; i < n; i++) {
      Tint eSum = M(2 * iPair, i) + M(2 * iPair + 1, i);
      if (eSum != 0)
        return false;
    }
  }
  return true;
}

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamilyZbasis(MyMatrix<T> const &eMat) {
#ifdef TIMINGS
  std::cerr << "Begining of ExtractInvariantVectorFamilyZbasis\n";
  SingletonTime time1;
#endif
  LLLreduction<T, Tint> recLLL = LLLreducedBasis<T, Tint>(eMat);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|LLL|=" << ms(time1,time2) << "\n";
#endif
  //  std::cerr << "recLLL.GramMatRed=\n";
  //  WriteMatrix(std::cerr, recLLL.GramMatRed);
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(recLLL.Pmat);
  //
  //  std::cerr << "eMat=\n";
  //  WriteMatrix(std::cerr, eMat);
  MyMatrix<T> eMatRed = Pmat_T * eMat * TransposedMat(Pmat_T);
  //  std::cerr << "eMatRed=\n";
  //  WriteMatrix(std::cerr, eMatRed);
  //
  MyMatrix<Tint> SHVred =
      ExtractInvariantVectorFamilyZbasis_Kernel<T, Tint>(eMatRed);
  MyMatrix<Tint> SHVret = SHVred * recLLL.Pmat;
#ifdef TIMINGS
  SingletonTime time3;
  std::cerr << "|ExtractInvariantVectorFamilyZbasis_Kernel|=" << ms(time2,time3) << "\n";
#endif
  return SHVret;
}

#endif //  SRC_LATT_INVARIANTVECTORFAMILY_H_
