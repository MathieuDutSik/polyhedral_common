#ifndef DEFINE_INVARIANT_VECTOR_FAMILY_H
#define DEFINE_INVARIANT_VECTOR_FAMILY_H

#include "ShortestUniversal.h"
#include "Temp_Positivity.h"

template <typename T, typename Tint>
MyMatrix<Tint> ExtractInvariantVectorFamily(
    MyMatrix<T> const &eMat,
    std::function<bool(MyMatrix<Tint> const &)> const &fCorrect) {
  T MaxNorm = MaximumDiagonal(eMat);
#ifdef TIMINGS
  std::cerr << "Begining of ExtractInvariantVectorFamily\n";
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
#endif
  MyMatrix<Tint> SHVall = T_ShortVector<T, Tint>(eMat, MaxNorm);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
  std::cerr
      << "Time(T_ShortVector) = "
      << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count()
      << "\n";
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
  std::chrono::time_point<std::chrono::system_clock> time3 =
      std::chrono::system_clock::now();
  std::cerr
      << "Time(Preparation) = "
      << std::chrono::duration_cast<std::chrono::seconds>(time3 - time2).count()
      << "\n";
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
  std::chrono::time_point<std::chrono::system_clock> time1 =
      std::chrono::system_clock::now();
#endif
  LLLreduction<T, Tint> recLLL = LLLreducedBasis<T, Tint>(eMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 =
      std::chrono::system_clock::now();
  std::cerr
      << "Time(LLL) = "
      << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count()
      << "\n";
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
  std::chrono::time_point<std::chrono::system_clock> time3 =
      std::chrono::system_clock::now();
  std::cerr
      << "|ExtractInvariantVectorFamilyZbasis_Kernel| = "
      << std::chrono::duration_cast<std::chrono::seconds>(time3 - time2).count()
      << "\n";
#endif
  return SHVret;
}

#endif
