#ifndef STRICT_POSITIVITY_INCLUDE
#define STRICT_POSITIVITY_INCLUDE

#include "Copositivity.h"
#include "Temp_PerfectForm.h"
#include "Temp_Tspace_General.h"

template <typename T, typename Tint> struct TestStrictPositivity {
  bool result;
  MyMatrix<Tint> RealizingFamily;
  MyVector<T> ListCoeff;
  MyMatrix<T> CertificateNonStrictlyPositive;
};

template <typename T, typename Tint>
PosRelRes<T> SearchForExistenceStrictPositiveRelation(MyMatrix<Tint> const &SHV,
                                                      MyMatrix<T> const &eMat) {
  MyVector<T> eMatVect = SymmetricMatrixToVector(eMat);
  int dimSymm = eMatVect.size();
  int nbBlock = SHV.rows();
  MyMatrix<T> ConeClassical = GetNakedPerfectConeClassical<T, Tint>(SHV);
  MyMatrix<T> TestMat(nbBlock + 1, dimSymm);
  for (int iBlock = 0; iBlock < nbBlock; iBlock++)
    TestMat.row(iBlock) = ConeClassical.row(iBlock);
  MyVector<T> V = -eMatVect;
  AssignMatrixRow(TestMat, nbBlock, V);
  return SearchPositiveRelationSimple(TestMat);
}

template <typename T, typename Tint>
TestStrictPositivity<T, Tint>
TestingAttemptStrictPositivity(MyMatrix<T> const &eMat,
                               MyMatrix<Tint> const &InitialBasis) {
  if (!IsSymmetricMatrix(eMat)) {
    std::cerr << "The matrix should be symmetric\n";
    throw TerminalException{1};
  }
  int n = eMat.rows();
  LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
  int dimLinSpa = LinSpa.ListMat.size();
  std::cerr << "dimLinSpa=" << dimLinSpa << "\n";
  MyVector<T> eMatExpr = LINSPA_GetVectorOfMatrixExpression(LinSpa, eMat);
  std::cerr << "eMatExpr=";
  WriteVector(std::cerr, eMatExpr);
  //
  std::function<bool(MyMatrix<T>)> IsAdmissible =
      [&](MyMatrix<T> const &eMatI) -> bool {
    CopositivityEnumResult<Tint> CopoRes;
    T MaxNorm = 1;
    RequestCopositivity<T> CopoReq{MaxNorm, true};
    std::cerr << "Case 1 eMatI=\n";
    WriteMatrix(std::cerr, eMatI);
    CopoRes =
        EnumerateCopositiveShortVector<T, Tint>(eMatI, InitialBasis, CopoReq);
    return CopoRes.test;
  };
  std::function<Tshortest<T, Tint>(MyMatrix<T>)> ShortestFunction =
      [&](MyMatrix<T> const &eMatI) -> Tshortest<T, Tint> {
    std::cerr << "Case 2 eMatI=\n";
    WriteMatrix(std::cerr, eMatI);
    return T_CopositiveShortestVector<T, Tint>(eMatI, InitialBasis);
  };
  RecShort<T, Tint> eRecShort{IsAdmissible, ShortestFunction};
  MyMatrix<T> SearchMatrix = AnLattice<T>(n) / T(2);
  int nbIter = 0;
  while (true) {
    nbIter++;
    //    std::cerr << "Before CopositiveShortestVector nbIter=" << nbIter <<
    //    "\n";
    Tshortest<T, Tint> RecSHV =
        T_CopositiveShortestVector<T, Tint>(SearchMatrix, InitialBasis);
    //    std::cerr << "Before GetNakedPerfectCone nbIter=" << nbIter << "\n";
    NakedPerfect<T, Tint> eNaked =
        GetNakedPerfectCone(LinSpa, SearchMatrix, RecSHV);
    int nbBlock = eNaked.ListBlock.size();

    T ScalMat = MatrixScalarProduct(SearchMatrix, eMat);
    //    std::cerr << "ScalMat=" << ScalMat << "\n";
    if (ScalMat <= 0) {
      std::cerr << "We obtained proof of non strict CP property\n";
      return {false, {}, {}, SearchMatrix};
    }
    //
    std::cerr << "Before computing ConeClassical\n";
    MyMatrix<T> ConeClassical =
        GetNakedPerfectConeClassical<T, Tint>(eNaked.SHVred);
    std::cerr << "Before SearchPositiveRelationSimple RankMat(ConeClassical)="
              << RankMat(ConeClassical) << " nbIter=" << nbIter << "\n";
    PosRelRes<T> eRel =
        SearchForExistenceStrictPositiveRelation(eNaked.SHVred, eMat);
    std::cerr << "We have PosRelRes eTestExist=" << eRel.eTestExist << "\n";
    if (eRel.eTestExist) {
      std::cerr << "TheRelat=";
      WriteVector(std::cerr, eRel.TheRelat);
      std::vector<int> ListIdx;
      for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
        if (eRel.TheRelat(iBlock) != 0)
          ListIdx.push_back(iBlock);
      }
      int nbIdx = ListIdx.size();
      MyMatrix<Tint> RealizingFamily(nbIdx, n);
      MyVector<T> eVectRet(nbIdx);
      MyMatrix<T> TotalSum = ZeroMatrix<T>(n, n);
      for (int iIdx = 0; iIdx < nbIdx; iIdx++) {
        int iBlock = ListIdx[iIdx];
        int iSHV = eNaked.ListBlock[iBlock][0];
        T eVal = eRel.TheRelat(iBlock) / eRel.TheRelat(nbBlock);
        RealizingFamily.row(iIdx) = eNaked.SHV.row(iSHV);
        eVectRet(iIdx) = eVal;
        for (int u = 0; u < n; u++)
          for (int v = 0; v < n; v++)
            TotalSum(u, v) +=
                eVal * RealizingFamily(iIdx, u) * RealizingFamily(iIdx, v);
      }
      std::cerr << "nbIter=" << nbIter << "\n";
      std::cerr << "TotalSum=\n";
      WriteMatrix(std::cerr, TotalSum);
      return {true, RealizingFamily, eVectRet, {}};
    }
    // Now trying to do the flipping
    std::cerr << "Before FindViolatedFace nbIter=" << nbIter << "\n";
    MyVector<T> eMatVect = SymmetricMatrixToVector(eMat);
    Face eFace = FindViolatedFace(ConeClassical, eMatVect);
    std::cerr << "eFace.count=" << eFace.count()
              << " eFace.size=" << eFace.size() << "\n";
    std::cerr << "Before FindFacetInequality nbIter=" << nbIter << "\n";
    std::cerr << "RankMat(ConeClassical)=" << RankMat(ConeClassical) << "\n";
    int dimSymm = ConeClassical.cols();
    MyVector<T> eFacet = FindFacetInequality(ConeClassical, eFace);
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      T eSum = 0;
      for (int u = 0; u < dimSymm; u++)
        eSum += eFacet(u) * ConeClassical(iBlock, u);
      std::cerr << "iBlock=" << iBlock << " eSum=" << eSum << "\n";
    }
    std::cerr << "Before MatDir creation nbIter=" << nbIter << "\n";
    //
    MyVector<T> Wvect = GetSymmetricMatrixWeightVector<T>(n);
    MyVector<T> Vexpand(dimSymm);
    for (int i = 0; i < dimSymm; i++)
      Vexpand(i) = eFacet(i) / Wvect(i);
    MyMatrix<T> eMatDir = VectorToSymmetricMatrix(Vexpand, n);
    T ScalDir = MatrixScalarProduct(eMatDir, eMat);
    std::cerr << "ScalDir=" << ScalDir << "\n";

    //    MyMatrix<T> eMatDir=LINSPA_GetMatrixInTspace(LinSpa, eFacet);
    std::cerr << "Before KernelFlipping nbIter=" << nbIter << "\n";
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> ePair =
        Kernel_Flipping_Perfect<T, Tint>(eRecShort, SearchMatrix, eMatDir);
    std::cerr << "NewMat=\n";
    WriteMatrix(std::cerr, ePair.first);
    std::cerr << "Before SearchMatrix assignation nbIter=" << nbIter << "\n";
    SearchMatrix = ePair.first;
  }
}

template <typename T, typename Tint>
void WriteStrictPositivityResult(
    std::ostream &os, TestStrictPositivity<T, Tint> const &StrictPos) {
  os << "result=" << StrictPos.result << "\n";
  if (StrictPos.result) {
    int n = StrictPos.RealizingFamily.cols();
    int nbBlock = StrictPos.RealizingFamily.rows();
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      T eVal = StrictPos.ListCoeff(iBlock);
      os << "iB=" << iBlock << " val=" << eVal << " V=";
      for (int i = 0; i < n; i++)
        os << StrictPos.RealizingFamily(iBlock, i) << " ";
      os << "\n";
    }
  } else {
    os << "Certificate of non strict completely positive matrix:\n";
    os << "Following matrix has negative scalar product with eMat and is "
          "copositive:\n";
    WriteMatrix(os, StrictPos.CertificateNonStrictlyPositive);
  }
}

#endif
