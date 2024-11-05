// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_COPOS_STRICTPOSITIVITY_H_
#define SRC_COPOS_STRICTPOSITIVITY_H_

// clang-format off
#include "Copositivity.h"
#include "Temp_PerfectForm.h"
#include "Tspace_General.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_STRICT_POSITIVITY
#endif

template <typename T, typename Tint> struct TestStrictPositivity {
  bool result;
  MyMatrix<Tint> RealizingFamily;
  MyVector<T> ListCoeff;
  MyMatrix<T> CertificateNonStrictlyPositive;
};

template <typename T, typename Tint>
PosRelRes<T> SearchForExistenceStrictPositiveRelation(MyMatrix<Tint> const &SHV,
                                                      MyMatrix<T> const &eMat,
                                                      std::ostream &os) {
  MyVector<T> eMatVect = SymmetricMatrixToVector(eMat);
  int dimSymm = eMatVect.size();
  int nbBlock = SHV.rows();
  MyMatrix<T> ConeClassical = GetNakedPerfectConeClassical<T, Tint>(SHV);
  MyMatrix<T> TestMat(nbBlock + 1, dimSymm);
  for (int iBlock = 0; iBlock < nbBlock; iBlock++)
    TestMat.row(iBlock) = ConeClassical.row(iBlock);
  MyVector<T> V = -eMatVect;
  AssignMatrixRow(TestMat, nbBlock, V);
  return SearchPositiveRelationSimple_Direct(TestMat, os);
}

template <typename T, typename Tint>
TestStrictPositivity<T, Tint>
TestingAttemptStrictPositivity(MyMatrix<T> const &eMat,
                               MyMatrix<Tint> const &InitialBasis,
                               std::ostream &os) {
  if (!IsSymmetricMatrix(eMat)) {
    std::cerr << "The matrix should be symmetric\n";
    throw TerminalException{1};
  }
  int n = eMat.rows();
  LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
  MyVector<T> eMatExpr = LINSPA_GetVectorOfMatrixExpression(LinSpa, eMat);
#ifdef DEBUG_STRICT_POSITIVITY
  int dimLinSpa = LinSpa.ListMat.size();
  std::cerr << "dimLinSpa=" << dimLinSpa << "\n";
  std::cerr << "eMatExpr=";
  WriteVectorNoDim(std::cerr, eMatExpr);
#endif
  //
  std::function<bool(MyMatrix<T>)> IsAdmissible =
      [&](MyMatrix<T> const &eMatI) -> bool {
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "IsAdmissible eMatI=\n";
    WriteMatrix(std::cerr, eMatI);
#endif
    CopositivityTestResult<Tint> result =
        TestCopositivity<T, Tint>(eMatI, InitialBasis, os);
    return result.test;
  };
  std::function<Tshortest<T, Tint>(MyMatrix<T>)> ShortestFunction =
      [&](MyMatrix<T> const &eMatI) -> Tshortest<T, Tint> {
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "ShortestFunction eMatI=\n";
    WriteMatrix(std::cerr, eMatI);
#endif
    return CopositiveShortestVector<T, Tint>(eMatI, InitialBasis, os);
  };
  RecShort<T, Tint> eRecShort{IsAdmissible, ShortestFunction};
  MyMatrix<T> SearchMatrix = AnLattice<T>(n) / T(2);
#ifdef DEBUG_STRICT_POSITIVITY
  int nbIter = 0;
#endif
  while (true) {
#ifdef DEBUG_STRICT_POSITIVITY
    nbIter++;
#endif
    Tshortest<T, Tint> RecSHV =
        CopositiveShortestVector<T, Tint>(SearchMatrix, InitialBasis, os);
    NakedPerfect<T, Tint> eNaked =
        GetNakedPerfectCone(LinSpa, SearchMatrix, RecSHV);
    int nbBlock = eNaked.ListBlock.size();

    T ScalMat = MatrixScalarProduct(SearchMatrix, eMat);
    if (ScalMat <= 0) {
#ifdef DEBUG_STRICT_POSITIVITY
      std::cerr << "We obtained proof of non strict CP property\n";
#endif
      return {false, {}, {}, SearchMatrix};
    }
    //
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "Before computing ConeClassical\n";
#endif
    MyMatrix<T> ConeClassical =
        GetNakedPerfectConeClassical<T, Tint>(eNaked.SHVred);
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "Before SearchPositiveRelationSimple RankMat(ConeClassical)="
              << RankMat(ConeClassical) << " nbIter=" << nbIter << "\n";
#endif
    PosRelRes<T> eRel =
        SearchForExistenceStrictPositiveRelation(eNaked.SHVred, eMat, os);
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "We have PosRelRes eTestExist=" << eRel.eTestExist << "\n";
#endif
    if (eRel.eTestExist) {
      MyVector<T> const &V = *eRel.TheRelat;
#ifdef DEBUG_STRICT_POSITIVITY
      std::cerr << "TheRelat=";
      WriteVectorNoDim(std::cerr, V);
#endif
      std::vector<int> ListIdx;
      for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
        if (V(iBlock) != 0)
          ListIdx.push_back(iBlock);
      }
      int nbIdx = ListIdx.size();
      MyMatrix<Tint> RealizingFamily(nbIdx, n);
      MyVector<T> eVectRet(nbIdx);
      MyMatrix<T> TotalSum = ZeroMatrix<T>(n, n);
      for (int iIdx = 0; iIdx < nbIdx; iIdx++) {
        int iBlock = ListIdx[iIdx];
        int iSHV = eNaked.ListBlock[iBlock][0];
        T eVal = V(iBlock) / V(nbBlock);
        RealizingFamily.row(iIdx) = eNaked.SHV.row(iSHV);
        eVectRet(iIdx) = eVal;
        for (int u = 0; u < n; u++)
          for (int v = 0; v < n; v++)
            TotalSum(u, v) +=
                eVal * RealizingFamily(iIdx, u) * RealizingFamily(iIdx, v);
      }
#ifdef DEBUG_STRICT_POSITIVITY
      std::cerr << "nbIter=" << nbIter << "\n";
      std::cerr << "TotalSum=\n";
      WriteMatrix(std::cerr, TotalSum);
#endif
      return {true, RealizingFamily, eVectRet, {}};
    }
    // Now trying to do the flipping
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "Before FindViolatedFace nbIter=" << nbIter << "\n";
#endif
    MyVector<T> eMatVect = SymmetricMatrixToVector(eMat);
    Face eFace = FindViolatedFace(ConeClassical, eMatVect, os);
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "eFace.count=" << eFace.count()
              << " eFace.size=" << eFace.size() << "\n";
    std::cerr << "Before FindFacetInequality nbIter=" << nbIter << "\n";
    std::cerr << "RankMat(ConeClassical)=" << RankMat(ConeClassical) << "\n";
#endif
    int dimSymm = ConeClassical.cols();
    MyVector<T> eFacet = FindFacetInequality(ConeClassical, eFace);
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      T eSum = 0;
      for (int u = 0; u < dimSymm; u++)
        eSum += eFacet(u) * ConeClassical(iBlock, u);
    }
    //
    MyVector<T> Wvect = GetSymmetricMatrixWeightVector<T>(n);
    MyVector<T> Vexpand(dimSymm);
    for (int i = 0; i < dimSymm; i++)
      Vexpand(i) = eFacet(i) / Wvect(i);
    MyMatrix<T> eMatDir = VectorToSymmetricMatrix(Vexpand, n);
    T ScalDir = MatrixScalarProduct(eMatDir, eMat);

    //    MyMatrix<T> eMatDir=LINSPA_GetMatrixInTspace(LinSpa, eFacet);
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "Before KernelFlipping nbIter=" << nbIter << "\n";
#endif
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> ePair =
        Kernel_Flipping_Perfect<T, Tint>(eRecShort, SearchMatrix, eMatDir);
#ifdef DEBUG_STRICT_POSITIVITY
    std::cerr << "NewMat=\n";
    WriteMatrix(std::cerr, ePair.first);
    std::cerr << "Before SearchMatrix assignation nbIter=" << nbIter << "\n";
#endif
    SearchMatrix = ePair.first;
  }
}

template <typename T, typename Tint>
void WriteStrictPositivityResult(
    std::ostream &os, std::string const &OutFormat,
    TestStrictPositivity<T, Tint> const &StrictPos) {
  if (OutFormat == "classic") {
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
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(result:=" << GAP_logical(StrictPos.result);
    if (StrictPos.result) {
      int nbBlock = StrictPos.RealizingFamily.rows();
      os << ", RealizingFamily:=[";
      for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
        if (iBlock > 0)
          os << ",\n";
        T eVal = StrictPos.ListCoeff(iBlock);
        MyVector<Tint> V = GetMatrixRow(StrictPos.RealizingFamily, iBlock);
        os << "rec(val:=" << eVal << ", V:=" << StringVectorGAP(V) << ")";
      }
      os << "]";
    } else {
      os << ", Certificate:="
         << StringMatrixGAP(StrictPos.CertificateNonStrictlyPositive);
    }
    os << ");\n";
    return;
  }
  if (OutFormat == "PYTHON") {
    os << "{\"result\":" << PYTHON_logical(StrictPos.result);
    if (StrictPos.result) {
      int nbBlock = StrictPos.RealizingFamily.rows();
      os << ", \"RealizingFamily\":[";
      for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
        if (iBlock > 0)
          os << ",";
        T eVal = StrictPos.ListCoeff(iBlock);
        MyVector<Tint> V = GetMatrixRow(StrictPos.RealizingFamily, iBlock);
        os << "{\"val\":" << eVal << ", \"V\":" << StringVectorPYTHON(V) << "}";
      }
      os << "]";
    } else {
      os << ", \"Certificate\":"
         << StringMatrixPYTHON(StrictPos.CertificateNonStrictlyPositive);
    }
    os << "}\n";
    return;
  }
  std::cerr << "WriteStrictPositivityResult: Failed to find a matching entry "
               "for output\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_COPOS_STRICTPOSITIVITY_H_
// clang-format on
