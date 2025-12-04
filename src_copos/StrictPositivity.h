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

#ifdef DISABLE_DEBUG_STRICT_POSITIVITY
#undef DEBUG_STRICT_POSITIVITY
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
#ifdef SANITY_CHECK_STRICT_POSITIVITY
  if (!IsSymmetricMatrix(eMat)) {
    std::cerr << "STR: The matrix should be symmetric\n";
    throw TerminalException{1};
  }
#endif
  int n = eMat.rows();
  LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
  MyVector<T> eMatExpr = LINSPA_GetVectorOfMatrixExpression(LinSpa, eMat);
#ifdef DEBUG_STRICT_POSITIVITY
  int dimLinSpa = LinSpa.ListMat.size();
  os << "STR: dimLinSpa=" << dimLinSpa << "\n";
  os << "STR: eMatExpr=";
  WriteVectorNoDim(os, eMatExpr);
#endif
  //
  auto f_admissible=[&](MyMatrix<T> const &eMatI) -> bool {
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: IsAdmissible eMatI=\n";
    WriteMatrix(os, eMatI);
#endif
    CopositivityTestResult<Tint> result =
        TestStrictCopositivity<T, Tint>(eMatI, InitialBasis, os);
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: IsAdmissible result=" << result.test << "\n";
#endif
    return result.test;
  };
  auto f_shortest=[&](MyMatrix<T> const &eMatI) -> Tshortest<T, Tint> {
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: ShortestFunction eMatI=\n";
    WriteMatrix(os, eMatI);
#endif
    return CopositiveShortestVector<T,Tint>(eMatI, InitialBasis, os);
  };
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
      GetNakedPerfectCone(LinSpa, SearchMatrix, RecSHV, os);
    int nbBlock = eNaked.ListBlock.size();

    T ScalMat = MatrixScalarProduct(SearchMatrix, eMat);
    if (ScalMat <= 0) {
#ifdef DEBUG_STRICT_POSITIVITY
      os << "STR: We obtained proof of non strict CP property\n";
#endif
      return {false, {}, {}, SearchMatrix};
    }
    //
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: Before computing ConeClassical\n";
#endif
    MyMatrix<T> ConeClassical =
        GetNakedPerfectConeClassical<T, Tint>(eNaked.SHVred);
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: Before SearchPositiveRelationSimple RankMat(ConeClassical)="
       << RankMat(ConeClassical) << " nbIter=" << nbIter << "\n";
#endif
    PosRelRes<T> eRel =
        SearchForExistenceStrictPositiveRelation(eNaked.SHVred, eMat, os);
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: We have PosRelRes eTestExist=" << eRel.eTestExist << "\n";
#endif
    if (eRel.eTestExist) {
      MyVector<T> const &V = *eRel.TheRelat;
#ifdef DEBUG_STRICT_POSITIVITY
      os << "STR: TheRelat=";
      WriteVectorNoDim(os, V);
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
        for (int u = 0; u < n; u++) {
          for (int v = 0; v < n; v++) {
            TotalSum(u, v) +=
                eVal * RealizingFamily(iIdx, u) * RealizingFamily(iIdx, v);
          }
        }
      }
#ifdef DEBUG_STRICT_POSITIVITY
      os << "STR: nbIter=" << nbIter << "\n";
      os << "STR: TotalSum=\n";
      WriteMatrix(os, TotalSum);
#endif
      return {true, RealizingFamily, eVectRet, {}};
    }
    // Now trying to do the flipping
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: Before FindViolatedFace nbIter=" << nbIter << "\n";
#endif
    MyVector<T> eMatVect = SymmetricMatrixToVector(eMat);
    Face eFace = FindViolatedFace(ConeClassical, eMatVect, os);
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: eFace.count=" << eFace.count()
       << " eFace.size=" << eFace.size() << "\n";
    os << "STR: Before FindFacetInequality nbIter=" << nbIter << "\n";
    os << "STR: RankMat(ConeClassical)=" << RankMat(ConeClassical) << "\n";
#endif
    int dimSymm = ConeClassical.cols();
    MyVector<T> eFacet = FindFacetInequality(ConeClassical, eFace);
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      T eSum(0);
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
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: Before KernelFlipping nbIter=" << nbIter << "\n";
#endif
    std::pair<MyMatrix<T>, Tshortest<T, Tint>> ePair =
      Kernel_Flipping_Perfect<T, Tint,decltype(f_admissible),decltype(f_shortest)>(f_admissible, f_shortest, SearchMatrix, eMatDir, os);
#ifdef DEBUG_STRICT_POSITIVITY
    os << "STR: Before SearchMatrix assignation nbIter=" << nbIter << "\n";
    os << "STR: NewMat=\n";
    WriteMatrix(os, ePair.first);
    os << "STR: Before SearchMatrix assignation nbIter=" << nbIter << "\n";
#endif
    SearchMatrix = ePair.first;
  }
}

template <typename T, typename Tint>
void WriteStrictPositivityResult(
    std::ostream &os_out, std::string const &OutFormat,
    TestStrictPositivity<T, Tint> const &StrictPos) {
  if (OutFormat == "classic") {
    os_out << "result=" << StrictPos.result << "\n";
    if (StrictPos.result) {
      int n = StrictPos.RealizingFamily.cols();
      int nbBlock = StrictPos.RealizingFamily.rows();
      for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
        T eVal = StrictPos.ListCoeff(iBlock);
        os_out << "iB=" << iBlock << " val=" << eVal << " V=";
        for (int i = 0; i < n; i++)
          os_out << StrictPos.RealizingFamily(iBlock, i) << " ";
        os_out << "\n";
      }
    } else {
      os_out << "Certificate of non strict completely positive matrix:\n";
      os_out << "Following matrix has negative scalar product with eMat and is "
            "copositive:\n";
      WriteMatrix(os_out, StrictPos.CertificateNonStrictlyPositive);
    }
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return rec(result:=" << GAP_logical(StrictPos.result);
    if (StrictPos.result) {
      int nbBlock = StrictPos.RealizingFamily.rows();
      os_out << ", RealizingFamily:=[";
      for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
        if (iBlock > 0)
          os_out << ",\n";
        T eVal = StrictPos.ListCoeff(iBlock);
        MyVector<Tint> V = GetMatrixRow(StrictPos.RealizingFamily, iBlock);
        os_out << "rec(val:=" << eVal << ", V:=" << StringVectorGAP(V) << ")";
      }
      os_out << "]";
    } else {
      os_out << ", Certificate:="
         << StringMatrixGAP(StrictPos.CertificateNonStrictlyPositive);
    }
    os_out << ");\n";
    return;
  }
  if (OutFormat == "PYTHON") {
    os_out << "{\"result\":" << PYTHON_logical(StrictPos.result);
    if (StrictPos.result) {
      int nbBlock = StrictPos.RealizingFamily.rows();
      os_out << ", \"RealizingFamily\":[";
      for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
        if (iBlock > 0)
          os_out << ",";
        T eVal = StrictPos.ListCoeff(iBlock);
        MyVector<Tint> V = GetMatrixRow(StrictPos.RealizingFamily, iBlock);
        os_out << "{\"val\":" << eVal << ", \"V\":" << StringVectorPYTHON(V) << "}";
      }
      os_out << "]";
    } else {
      os_out << ", \"Certificate\":"
         << StringMatrixPYTHON(StrictPos.CertificateNonStrictlyPositive);
    }
    os_out << "}\n";
    return;
  }
  std::cerr << "STR: WriteStrictPositivityResult: Failed to find a matching entry "
               "for output\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_COPOS_STRICTPOSITIVITY_H_
// clang-format on
