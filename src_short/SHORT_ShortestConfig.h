// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SHORT_SHORT_SHORTESTCONFIG_H_
#define SRC_SHORT_SHORT_SHORTESTCONFIG_H_

// clang-format off
#include "ShortestUniversal.h"
#include "COMB_Combinatorics.h"
#include "InvariantVectorFamily.h"
#include "LatticeDefinitions.h"
#include "MAT_MatrixInt.h"
#include "MAT_functions.h"
#include "MatrixGroup.h"
#include "POLY_LinearProgramming.h"
#include "POLY_LinearProgramming_GLPK.h"
#include "POLY_PolytopeInt.h"
#include "Parallel_Classes_Types.h"
#include "Positivity.h"
#include "Tspace_General.h"
#include <algorithm>
#include <map>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SHORTEST_CONFIG
#endif

#ifdef TIMINGS
#define TIMINGS_SHORTEST_CONFIG
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SHORTEST_CONFIG
#endif

template <typename Tint>
MyMatrix<Tint> SHORT_CleanAntipodality(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  std::unordered_set<MyVector<Tint>> setVect;
  for (int i = 0; i < nbRow; i++) {
    MyVector<Tint> eRow = GetMatrixRow(M, i);
    int iColFound = -1;
    for (int iCol = 0; iCol < nbCol; iCol++)
      if (iColFound == -1)
        if (eRow(iCol) != 0)
          iColFound = iCol;
    if (iColFound == -1) {
      std::cerr << "Bug in SHORT_CleanAntipodality\n";
      std::cerr << "The vector is 0 which is not allowed\n";
      throw TerminalException{1};
    }
    if (eRow(iColFound) < 0)
      eRow = -eRow;
    setVect.insert(eRow);
  }
  std::vector<MyVector<Tint>> ListVect;
  for (auto &eV : setVect)
    ListVect.push_back(eV);
  return MatrixFromVectorFamily(ListVect);
}

template <typename T, typename Tint> struct ReplyRealizability {
  int eCase;
  bool reply;
  bool replyCone;
  MyMatrix<T> eMat;
  MyMatrix<Tint> SHV;
  MyMatrix<Tint> SHVclean;
};

template <typename T, typename Tint>
MyVector<T> SHORT_GetIneq_Tspace(std::vector<MyMatrix<T>> const &ListMat,
                                 MyVector<Tint> const &eVect) {
  int nbMat = ListMat.size();
  int dim = nbMat + 1;
  MyVector<T> eIneq(dim);
  eIneq(0) = -1;
  for (int iMat = 0; iMat < nbMat; iMat++) {
    T scal = EvaluationQuadForm<T, Tint>(ListMat[iMat], eVect);
    eIneq(1 + iMat) = scal;
  }
  return eIneq;
}

template <typename T, typename Tint>
ReplyRealizability<T, Tint> SHORT_TestRealizabilityShortestFamilyEquivariant(
    std::vector<MyVector<Tint>> const &ListVect,
    std::vector<MyMatrix<T>> const &ListMat, bool const &NoExtension,
    std::string const &TheMethod, std::ostream &os) {
  ReplyRealizability<T, Tint> eRes;
  int n = ListVect[0].size();
  int nbVect = ListVect.size();
  int nbMat = ListMat.size();
  std::vector<MyVector<Tint>> ListVectTotUnred;
  for (auto &eVect : ListVect) {
    ListVectTotUnred.push_back(eVect);
    ListVectTotUnred.push_back(-eVect);
  }
  std::vector<MyVector<Tint>> ListVectTot = VectorAsSet(ListVectTotUnred);
  MyMatrix<T> MatrixValues(nbVect, nbMat);
  for (int iVect = 0; iVect < nbVect; iVect++)
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyVector<Tint> eVect = ListVect[iVect];
      MyMatrix<T> eMat = ListMat[iMat];
      T eVal = EvaluationQuadForm(eMat, eVect);
      MatrixValues(iVect, iMat) = eVal;
    }
  MyMatrix<T> MatrixValuesDiff(nbVect - 1, nbMat);
  for (int iVect = 0; iVect < nbVect - 1; iVect++)
    for (int iMat = 0; iMat < nbMat; iMat++)
      MatrixValuesDiff(iVect, iMat) =
          MatrixValues(iVect + 1, iMat) - MatrixValues(0, iMat);
  MyMatrix<T> NSP = NullspaceTrMat(MatrixValuesDiff);
  int dimSpa = NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    MyVector<T> eRow = GetMatrixRow(NSP, iDim);
    TheBasis[iDim] = GetMatrixFromBasis(ListMat, eRow);
  }
  //
  // Forming the vector family
  //
  std::unordered_set<MyVector<Tint>> ListForbiddenVector;
  for (auto &eVect : ListVectTot)
    ListForbiddenVector.insert(eVect);
  ListForbiddenVector.insert(ZeroVector<Tint>(n));
  std::unordered_set<MyVector<Tint>> TheFamilyVect;
  for (auto &eVect : ListVect) {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < 2; j++) {
        int shift = -1 + 2 * j;
        MyVector<Tint> rVect = eVect;
        rVect[i] += shift;
        auto iter = ListForbiddenVector.find(rVect);
        if (iter == ListForbiddenVector.end())
          TheFamilyVect.insert(rVect);
      }
  }
#ifdef DEBUG_SHORTEST_CONFIG
  std::cerr << "|TheFamilyVect|=" << TheFamilyVect.size() << "\n";
#endif
  auto GetListIneq = [&]() -> MyMatrix<T> {
    int siz = TheFamilyVect.size();
    MyMatrix<T> RetMat(siz, dimSpa + 1);
    int idx = 0;
    for (auto &eVect : TheFamilyVect) {
      MyVector<T> vVect = SHORT_GetIneq_Tspace<T, Tint>(TheBasis, eVect);
      for (int i = 0; i <= dimSpa; i++)
        RetMat(idx, i) = vVect(i);
      idx++;
    }
    return RetMat;
  };
  MyVector<T> ToBeMinimized(dimSpa + 1);
  ToBeMinimized(0) = 0;
  bool IsZero = true;
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    T scal = EvaluationQuadForm(TheBasis[iDim], ListVect[0]);
    if (scal != 0)
      IsZero = false;
    ToBeMinimized(iDim + 1) = scal;
  }
  if (IsZero) {
    eRes.eCase = 1;
    eRes.reply = false;
    eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
    std::cerr << "RETURN case 1\n";
#endif
    return eRes;
  }
  MyVector<T> ZerVect =
      SHORT_GetIneq_Tspace<T, Tint>(TheBasis, ZeroVector<Tint>(n));
  //
  int nbIter = 0;
  T eOptimal, eOptimalPrev;
  while (true) {
    if (nbIter > 2) {
      if (eOptimalPrev > eOptimal) {
        std::cerr << "Optimal values should be increasing\n";
        throw TerminalException{1};
      }
      eOptimalPrev = eOptimal;
    }
    nbIter++;
#ifdef DEBUG_SHORTEST_CONFIG
    std::cerr << "nbIter=" << nbIter << "\n";
#endif
    MyMatrix<T> ListIneq = GetListIneq();
    for (auto &eVect : ListVectTot) {
      auto iter = TheFamilyVect.find(eVect);
      if (iter != TheFamilyVect.end()) {
        std::cerr << "We find one of ListVectTot in TheFamilyVect\n";
        throw TerminalException{1};
      }
    }
    int sizFamVect = TheFamilyVect.size();
    for (int i = 0; i < sizFamVect; i++) {
      MyVector<T> rVect = GetMatrixRow(ListIneq, i);
      if (rVect == ZerVect) {
        eRes.eCase = 2;
        eRes.reply = false;
        eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
        std::cerr << "RETURN case 2\n";
#endif
        return eRes;
      }
    }
    MyMatrix<T> SetIneq = SortUnicizeMatrix(ListIneq);
    auto GetLpSolution = [&](MyMatrix<T> const &MatIneq,
                             MyVector<T> const &TheMinimized) -> LpSolution<T> {
      if (TheMethod == "cdd")
        return CDD_LinearProgramming(MatIneq, TheMinimized, os);
      if (TheMethod == "glpk_secure")
        return GLPK_LinearProgramming_Secure(MatIneq, TheMinimized, os);
      std::cerr << "We have TheMethod = " << TheMethod << "\n";
      throw TerminalException{1};
    };
    LpSolution<T> eSol = GetLpSolution(SetIneq, ToBeMinimized);
    if (!eSol.PrimalDefined && eSol.DualDefined) {
#ifdef DEBUG_SHORTEST_CONFIG
      std::cerr << "DualDefined but not primal defined\n";
#endif
      int nbIneqSet = SetIneq.size();
      MyVector<T> SumIneq = ZeroVector<T>(1 + nbIneqSet);
      for (int i = 0; i < nbIneqSet; i++) {
        MyMatrix<T> eRow = GetMatrixRow(SetIneq, i);
        SumIneq += eSol.DualSolution(i) * eRow;
      }
      MyVector<T> SumIneqRed = ZeroVector<T>(nbIneqSet);
      for (int i = 0; i < nbIneqSet; i++)
        SumIneqRed(i) = SumIneq(i + 1);
      if (SumIneq(0) < 0 && SumIneqRed == ZeroVector<T>(nbIneqSet)) {
        eRes.eCase = 3;
        eRes.reply = false;
        eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
        std::cerr << "RETURN case 3\n";
#endif
        return eRes;
      }
      std::cerr << "It seems we have a big problem here. Please correct\n";
      throw TerminalException{1};
    } else if (eSol.PrimalDefined && !eSol.DualDefined) {
#ifdef DEBUG_SHORTEST_CONFIG
      std::cerr << "PrimalDefined but not dual defined\n";
#endif
      // This is the primal_direction case
      MyVector<Tint> eVect = ListVect[0];
      MyVector<Tint> eVectB = 2 * eVect;
#ifdef DEBUG_SHORTEST_CONFIG
      std::cerr << "Inserting from PrimalDefined but not dual defined eVectB=";
      WriteVectorNoDim(std::cerr, eVectB);
#endif
      TheFamilyVect.insert(eVectB);
    } else {
#ifdef DEBUG_SHORTEST_CONFIG
      std::cerr << "We have optimal value\n";
#endif
#ifdef SANITY_CHECK_SHORTEST_CONFIG
      if (!eSol.PrimalDefined || !eSol.DualDefined) {
        std::cerr << "We have a real problem to solve. Please debug\n";
        throw TerminalException{1};
      }
#endif
      if (nbIter == 1)
        eOptimalPrev = eOptimal;
      eOptimal = eSol.OptimalValue;
#ifdef DEBUG_SHORTEST_CONFIG
      std::cerr << "eOptimal=" << eOptimal << " eOptimalPrev=" << eOptimalPrev
                << "\n";
#endif
      if (eOptimal > 1) {
        eRes.eCase = 4;
        eRes.reply = false;
        eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
        std::cerr << "RETURN case 4\n";
#endif
        return eRes;
      }
      if (eOptimal >= 1 && NoExtension) {
        eRes.eCase = 5;
        eRes.reply = false;
        eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
        std::cerr << "RETURN case 5\n";
#endif
        return eRes;
      }
      MyVector<T> eVectEmb = eSol.DirectSolution;
      MyVector<T> rVect = GetDirectSolutionExt(eSol);
      MyMatrix<T> eMatSec = GetMatrixFromBasis(TheBasis, eVectEmb);
      //
      // Some checks
      //
      int idx = 0;
      for (auto &eVect : TheFamilyVect) {
        MyVector<T> eRow = GetMatrixRow(ListIneq, idx);
        if (eRow != ZerVect) {
          T eNorm = EvaluationQuadForm(eMatSec, eVect);
          if (eNorm < 1) {
            std::cerr
                << "This should not happen if the matrix was well defined\n";
            throw TerminalException{1};
          }
          T eScal1 = ScalarProduct(rVect, eRow) + 1;
          if (eScal1 != eNorm) {
            std::cerr << "We should have eScal1 != eNorm\n";
            throw TerminalException{1};
          }
        }
        idx++;
      }
      bool IsFirst = true;
      T singleNorm;
      for (auto &eVect : ListVect) {
        T eNorm = EvaluationQuadForm(eMatSec, eVect);
        if (IsFirst) {
          singleNorm = eNorm;
          IsFirst = false;
        } else {
          if (eNorm != singleNorm) {
            std::cerr << "The norms are not consistent\n";
            throw TerminalException{1};
          }
        }
      }
      T eOne = 1;
      T CritNorm = std::min(eOne, singleNorm);
      if (CritNorm != eOptimal) {
        std::cerr << "CritNorm=" << CritNorm << " eOptimal=" << eOptimal
                  << "\n";
        std::cerr << "The value of eOptimal has not been well set\n";
        throw TerminalException{1};
      }
      bool testPosDef = IsPositiveDefinite(eMatSec);
#ifdef DEBUG_SHORTEST_CONFIG
      std::cerr << "testPosDef=" << testPosDef << "\n";
#endif
      if (testPosDef) {
#ifdef DEBUG_SHORTEST_CONFIG
        std::cerr << "Before computation of T_ShortestVector\n";
        WriteMatrix(std::cerr, eMatSec);
#endif
        Tshortest<T, Tint> RecSHV = T_ShortestVector<T, Tint>(eMatSec, os);
#ifdef DEBUG_SHORTEST_CONFIG
        std::cerr << " After computation of T_ShortestVector\n";
#endif
        int nbRow = RecSHV.SHV.rows();
        std::vector<MyVector<Tint>> SHV(nbRow);
        for (int iRow = 0; iRow < nbRow; iRow++) {
          MyVector<Tint> eVect = GetMatrixRow(RecSHV.SHV, iRow);
          SHV[iRow] = eVect;
        }
        std::vector<MyVector<Tint>> SetSHV = VectorAsSet(SHV);
        bool testEqua = SetSHV == ListVectTot;
        if (testEqua) {
          eRes.eCase = 6;
          eRes.reply = true;
          eRes.replyCone = true;
          eRes.SHV = RecSHV.SHV;
          eRes.SHVclean = SHORT_CleanAntipodality(RecSHV.SHV);
          eRes.eMat = eMatSec;
#ifdef DEBUG_SHORTEST_CONFIG
          std::cerr << "RETURN case 6\n";
#endif
          return eRes;
        } else {
          std::vector<MyVector<Tint>> SHVdiff;
          for (auto &eVect : SetSHV) {
            int pos = PositionVect(ListVectTot, eVect);
            if (pos == -1)
              SHVdiff.push_back(eVect);
          }
          std::vector<MyVector<Tint>> DiffNew;
          for (auto &eVect : SHVdiff) {
            auto iter = TheFamilyVect.find(eVect);
            if (iter == TheFamilyVect.end())
              DiffNew.push_back(eVect);
          }
          if (DiffNew.size() > 0) {
            for (auto &eVect : DiffNew) {
#ifdef DEBUG_SHORTEST_CONFIG
              std::cerr << "Inserting from DiffNew eVect=";
              WriteVectorNoDim(std::cerr, eVect);
#endif
              TheFamilyVect.insert(eVect);
            }
          } else {
            if (IsSubset(SetSHV, ListVectTot)) {
              eRes.eCase = 7;
              eRes.reply = false;
              eRes.replyCone = true;
              eRes.SHV = RecSHV.SHV;
              eRes.SHVclean = SHORT_CleanAntipodality(RecSHV.SHV);
              eRes.eMat = eMatSec;
#ifdef DEBUG_SHORTEST_CONFIG
              std::cerr << "RETURN case 7\n";
#endif
              return eRes;
            } else {
              eRes.eCase = 8;
              eRes.reply = false;
              eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
              std::cerr << "RETURN case 8\n";
#endif
              return eRes;
            }
          }
        }
      } else {
        if (RankMat(eMatSec) < n) {
          MyVector<Tint> eVect3 =
              GetShortVectorDegenerate<T, Tint>(eMatSec, CritNorm, os);
          if (PositionVect(ListVectTot, eVect3) != -1)
            eVect3 *= 2;
          if (TheFamilyVect.find(eVect3) != TheFamilyVect.end()) {
            std::cerr << "eMatSec=\n";
            WriteMatrix(std::cerr, eMatSec);
            std::cerr << "eVect3=";
            WriteVectorNoDim(std::cerr, eVect3);
            std::cerr << "We have a clear error here\n";
            throw TerminalException{1};
          }
#ifdef DEBUG_SHORTEST_CONFIG
          std::cerr << "Inserting from GetShortVectorDegenerate eVect3=";
          WriteVectorNoDim(std::cerr, eVect3);
#endif
          TheFamilyVect.insert(eVect3);
        } else {
          bool StrictIneq = true;
          MyVector<Tint> eVect = GetShortIntegralVector<T, Tint>(
              eMatSec, CritNorm, StrictIneq, os);
          if (TheFamilyVect.find(eVect) != TheFamilyVect.end()) {
            std::cerr << "We have a clear error here\n";
            throw TerminalException{1};
          }
#ifdef DEBUG_SHORTEST_CONFIG
          std::cerr << "Inserting from GetShortVector eVect=";
          WriteVectorNoDim(std::cerr, eVect);
#endif
          TheFamilyVect.insert(eVect);
        }
      }
    }
  }
}

template <typename T, typename Tint> struct ShortIso {
  MyMatrix<T> GramMat;
  MyMatrix<Tint> SHVdisc;
};

template <typename T, typename Tint>
ShortIso<T, Tint> SHORT_GetInformation(MyMatrix<Tint> const &M,
                                       std::ostream &os) {
  int n = M.cols();
  int nbVect = M.rows();
#ifdef DEBUG_SHORTEST_CONFIG
  std::cerr << "nbVect=" << nbVect << "\n";
#endif
  MyMatrix<T> Amat = ZeroMatrix<T>(n, n);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    MyVector<Tint> eVect_Tint = GetMatrixRow(M, iVect);
    MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect_Tint);
    MyMatrix<T> pMat = RankOneMatrix(eVect_T);
    Amat += pMat;
  }
  MyMatrix<T> TheGramMat = Inverse(Amat);
  MyMatrix<Tint> SHV =
      ExtractInvariantVectorFamilyZbasis<T, Tint>(TheGramMat, os);
  MyMatrix<Tint> Mc1 = 2 * M;
  MyMatrix<Tint> Mc2 = -2 * M;
  MyMatrix<Tint> Mcopy = Concatenate(Mc1, Mc2);
  MyMatrix<Tint> SHVret = Concatenate(SHV, Mcopy);
  return {std::move(TheGramMat), std::move(SHVret)};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>> SHORT_TestEquivalence(MyMatrix<Tint> const &M1,
                                                    MyMatrix<Tint> const &M2,
                                                    std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx_value = int16_t;
  ShortIso<T, Tint> eRec1 = SHORT_GetInformation<T, Tint>(M1, os);
  ShortIso<T, Tint> eRec2 = SHORT_GetInformation<T, Tint>(M2, os);
  WeightMatrix<true, T, Tidx_value> WMat1 =
      T_TranslateToMatrix_QM_SHV<T, Tint, Tidx_value>(eRec1.GramMat,
                                                      eRec1.SHVdisc, os);
  WeightMatrix<true, T, Tidx_value> WMat2 =
      T_TranslateToMatrix_QM_SHV<T, Tint, Tidx_value>(eRec2.GramMat,
                                                      eRec2.SHVdisc, os);
  std::optional<Telt> eResEquiv =
      TestEquivalenceWeightMatrix<T, Telt>(WMat1, WMat2, os);
  if (!eResEquiv)
    return {};
  MyMatrix<T> SHV1_T = UniversalMatrixConversion<T, Tint>(eRec1.SHVdisc);
  MyMatrix<T> SHV2_T = UniversalMatrixConversion<T, Tint>(eRec2.SHVdisc);
  MyMatrix<T> MatEquiv_T = FindTransformation(SHV1_T, SHV2_T, *eResEquiv);
  if (!IsIntegralMatrix(MatEquiv_T)) {
    std::cerr << "Error, the matrix is not integral\n";
    throw TerminalException{1};
  }
  return UniversalMatrixConversion<Tint, T>(MatEquiv_T);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>> SHORT_GetStabilizer(MyMatrix<Tint> const &M,
                                                std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tgr = GraphListAdj;
  using Tidx_value = int16_t;
  ShortIso<T, Tint> eRec1 = SHORT_GetInformation<T, Tint>(M, os);
  WeightMatrix<true, T, Tidx_value> WMat =
      T_TranslateToMatrix_QM_SHV<T, Tint, Tidx_value>(eRec1.GramMat,
                                                      eRec1.SHVdisc, os);
  Tgroup GRP = GetStabilizerWeightMatrix<T, Tgr, Tgroup, Tidx_value>(WMat, os);
#ifdef DEBUG_SHORTEST_CONFIG
  os << "|GRP| = " << GRP.size() << "\n";
#endif
  MyMatrix<Tint> Mneg = -M;
  MyMatrix<Tint> Mtot = Concatenate(M, Mneg);
  std::vector<MyMatrix<Tint>> ListMatrGen;
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(eRec1.SHVdisc);
  std::vector<Telt> LGen = GRP.GeneratorsOfGroup();
  for (auto const &eGen : LGen) {
    MyMatrix<T> MatEquiv_T = FindTransformation(SHV_T, SHV_T, eGen);
    if (!IsIntegralMatrix(MatEquiv_T)) {
      std::cerr << "Problem in SHORT_GetStabilizer\n";
      std::cerr << "Error, the matrix is not integral\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> MatEquiv_i = UniversalMatrixConversion<Tint, T>(MatEquiv_T);
    ListMatrGen.push_back(MatEquiv_i);
  }
#ifdef DEBUG_SHORTEST_CONFIG
  std::cerr << "Exiting SHORT_GetStabilizer\n";
#endif
  return ListMatrGen;
}

template <typename T>
int GetPositionAntipodal(std::vector<MyVector<T>> const &ListVect,
                         MyVector<T> const &eVect) {
  MyVector<T> rVect1 = -eVect;
  int nbVect = ListVect.size();
  for (int iVect = 0; iVect < nbVect; iVect++) {
    if (ListVect[iVect] == eVect)
      return iVect;
    if (ListVect[iVect] == rVect1)
      return iVect;
  }
  return -1;
}

int KissingNumberUpperBound(int const &n) {
  std::vector<int> ListVal = {-1,  2,   6,   12,   24,   45,   78,   135, 240,
                              366, 554, 870, 1357, 2069, 3183, 4866, 7355};
  if (n <= 16) {
    return ListVal[n];
  }
  throw TerminalException{1};
}

template <typename Tint> struct SHVreduced {
  MyMatrix<Tint> SHVred;
  MyMatrix<Tint> ReductionMatrix;
};

template <typename Tint>
SHVreduced<Tint> SHORT_GetLLLreduction_Kernel(MyMatrix<Tint> const &eSHV) {
  int n = eSHV.rows();
  int nbVect = eSHV.cols();
  using Tfield = typename overlying_field<Tint>::field_type;
  auto GetGram = [&](MyMatrix<Tint> const &uSHV) -> MyMatrix<Tfield> {
    MyMatrix<Tfield> TheGram = ZeroMatrix<Tfield>(n, n);
    for (int iVect = 0; iVect < nbVect; iVect++)
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          TheGram(i, j) += uSHV(iVect, i) * uSHV(iVect, j);
    return TheGram;
  };
  MyMatrix<Tfield> TheGram = GetGram(eSHV);
  LLLreduction<Tfield, Tint> res = LLLreducedBasis<Tfield, Tint>(TheGram);
  MyMatrix<Tfield> TheRemainder = res.GramMatRed;
  MyMatrix<Tint> TheTrans = res.Pmat;
  MyMatrix<Tint> Pmat = TransposedMat(TheTrans);
  MyMatrix<Tint> eSHVred = eSHV * Pmat;
  if (!TestEqualityMatrix(GetGram(eSHVred), TheRemainder)) {
    std::cerr << "Matrix error somewhere\n";
    throw TerminalException{1};
  }
  return {std::move(eSHVred), std::move(Pmat)};
}

template <typename T> std::vector<MyMatrix<T>> StandardSymmetricBasis(int n) {
  std::vector<MyMatrix<T>> ListMat;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
      eMat(i, j) = 1;
      eMat(j, i) = 1;
      ListMat.push_back(eMat);
    }
  return ListMat;
}

template <typename T, typename Tint, typename Tgroup>
ReplyRealizability<T, Tint>
SHORT_TestRealizabilityShortestFamily(MyMatrix<Tint> const &Minput,
                                      std::string const &TheMethod,
                                      std::ostream &os) {
  SHVreduced<Tint> RecLLL = SHORT_GetLLLreduction_Kernel(Minput);
  MyMatrix<Tint> M = RecLLL.SHVred;
  int n = M.cols();
  std::vector<MyMatrix<Tint>> ListMatrGen =
      SHORT_GetStabilizer<T, Tint, Tgroup>(M, os);
  std::vector<MyMatrix<T>> StdBasis = StandardSymmetricBasis<T>(n);
  std::vector<MyMatrix<T>> ListGen_T;
  for (auto &eGen : ListMatrGen)
    ListGen_T.push_back(UniversalMatrixConversion<T, Tint>(eGen));
#ifdef DEBUG_SHORTEST_CONFIG
  os << "Before BasisInvariantForm\n";
#endif
  std::vector<MyMatrix<T>> ListMat = BasisInvariantForm(n, ListGen_T, os);
#ifdef DEBUG_SHORTEST_CONFIG
  os << " After BasisInvariantForm\n";
#endif
  int nbVect = M.rows();
  std::vector<MyVector<T>> ListRankOne;
  std::vector<MyVector<Tint>> ListVectWork;
  auto FuncInsertVector = [&StdBasis, &ListVectWork,
                           &ListRankOne](MyVector<Tint> const &eRow) -> void {
    if (GetPositionAntipodal(ListVectWork, eRow) == -1) {
      ListVectWork.push_back(eRow);
      ListRankOne.push_back(SHORT_GetIneq_Tspace<T, Tint>(StdBasis, eRow));
    }
  };
#ifdef DEBUG_SHORTEST_CONFIG
  os << "nbVect=" << nbVect << "\n";
#endif
  for (int iRow = 0; iRow < nbVect; iRow++) {
    MyVector<Tint> eRow = GetMatrixRow(M, iRow);
    FuncInsertVector(eRow);
  }
  int InitialSize = ListVectWork.size();
  MyMatrix<T> ListRankOne_mat = MatrixFromVectorFamily(ListRankOne);
  auto IsNeededInsert =
      [&StdBasis, &ListRankOne_mat](MyVector<Tint> const &eVect) -> bool {
    MyVector<T> VectImg = SHORT_GetIneq_Tspace<T, Tint>(StdBasis, eVect);
    std::optional<MyVector<T>> opt = SolutionMat(ListRankOne_mat, VectImg);
    return opt.has_value();
  };
  int TheRank = RankMat(ListRankOne_mat);
  bool NoExtension = false;
  if (PositionVect({n, n + 1, n + 2}, TheRank) != -1)
    NoExtension = true;
  while (true) {
    ReplyRealizability<T, Tint> RecTest;
    if (InitialSize > KissingNumberUpperBound(n)) {
      RecTest.eCase = 9;
      RecTest.reply = false;
      RecTest.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
      std::cerr << "RETURN case 9\n";
#endif
      return RecTest;
    }
    RecTest = SHORT_TestRealizabilityShortestFamilyEquivariant<T, Tint>(
        ListVectWork, ListMat, NoExtension, TheMethod, os);
    if (RecTest.reply) {
      bool replyRet = static_cast<int>(ListVectWork.size()) == InitialSize;
      std::vector<MyVector<Tint>> ListVectComplete;
      for (auto &eVect : ListVectWork) {
        ListVectComplete.push_back(eVect);
        ListVectComplete.push_back(-eVect);
      }
      MyMatrix<Tint> SHVret = MatrixFromVectorFamily(ListVectComplete);
      MyMatrix<Tint> TheTrans = RecLLL.ReductionMatrix;
      MyMatrix<Tint> InvTrans = Inverse(TheTrans);
      MyMatrix<T> TheTrans_T = UniversalMatrixConversion<T, Tint>(TheTrans);
      MyMatrix<T> RetMat = TheTrans_T * RecTest.eMat * TheTrans_T.transpose();
      //
      RecTest.reply = replyRet;
      RecTest.replyCone = true;
      RecTest.eMat = RetMat;
      RecTest.SHV = SHVret * InvTrans;
      RecTest.SHVclean = MatrixFromVectorFamily(ListVectWork) * InvTrans;
      return RecTest;
    }
    if (!RecTest.replyCone) {
      RecTest.reply = false;
      RecTest.replyCone = false;
      return RecTest;
    }
    int SizPrev = ListVectWork.size();
    int nbSHVclean = RecTest.SHVclean.rows();
    for (int iRow = 0; iRow < nbSHVclean; iRow++) {
      MyVector<Tint> eVect = GetMatrixRow(RecTest.SHVclean, iRow);
      if (IsNeededInsert(eVect))
        FuncInsertVector(eVect);
    }
    int SizAfter = ListVectWork.size();
#ifdef DEBUG_SHORTEST_CONFIG
    std::cerr << "SizPrev=" << SizPrev << " SizAfter=" << SizAfter << "\n";
#endif
    if (SizAfter == SizPrev) {
      RecTest.reply = false;
      RecTest.replyCone = false;
      return RecTest;
    }
  }
}

template <typename T, typename Tint> struct SHVshortest {
  MyMatrix<Tint> SHV;
};

template <typename T, typename Tint>
std::istream &operator>>(std::istream &is, SHVshortest<T, Tint> &obj) {
  MyMatrix<Tint> M;
  M = ReadMatrix<Tint>(is);
  obj = {M};
  return is;
}

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, SHVshortest<T, Tint> const &obj) {
  WriteMatrix(os, obj.SHV);
  return os;
}

template <typename T, typename Tint> struct SHVinvariant {
  size_t eInv;
};

template <typename T, typename Tint>
std::istream &operator>>(std::istream &is, SHVinvariant<T, Tint> &obj) {
  SHVinvariant<T, Tint> eInv;
  is >> eInv;
  obj = {eInv};
  return is;
}

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, SHVinvariant<T, Tint> const &obj) {
  os << obj.eInv;
  return os;
}

template <typename T, typename Tint>
bool operator==(SHVinvariant<T, Tint> const &x,
                SHVinvariant<T, Tint> const &y) {
  return x.eInv == y.eInv;
}

template <typename T, typename Tint>
bool operator<(SHVinvariant<T, Tint> const &x, SHVinvariant<T, Tint> const &y) {
  return x.eInv < y.eInv;
}

template <typename T, typename Tint>
SHVinvariant<T, Tint> SHORT_Invariant(MyMatrix<Tint> const &eSpann,
                                      std::ostream &os) {
  SHVshortest<T, Tint> eEnt{eSpann};
  ShortIso<T, Tint> eShIso = SHORT_GetInformation<T, Tint>(eSpann, os);
  size_t seed = 146;
  size_t eInvGV =
      GetInvariantGramShortest(eShIso.GramMat, eShIso.SHVdisc, seed, os);
  return {eInvGV};
}

template <typename T, typename Tint> struct equiv_info<SHVshortest<T, Tint>> {
  typedef MyMatrix<Tint> equiv_type;
};

template <typename T, typename Tint>
struct invariant_info<SHVshortest<T, Tint>> {
  typedef SHVinvariant<T, Tint> invariant_type;
};

template <typename Tint>
Tint SHORT_GetMaximumDeterminant(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  IteratorBinomial<mpz_class> eIter(nbRow, nbCol);
  std::vector<Face> ListFace = eIter.ListAllFace();
  bool IsFirst = true;
  Tint eRet = -2;
  for (auto &eFace : ListFace) {
    MyMatrix<Tint> Mred = SelectRow(M, eFace);
    Tint eDet = DeterminantMat(Mred);
    Tint eDetA = T_abs(eDet);
    if (IsFirst) {
      eRet = eDetA;
      IsFirst = false;
    } else {
      if (eDetA > eRet)
        eRet = eDetA;
    }
  }
  return eRet;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
SHORT_SpannSimplicial(MyMatrix<Tint> const &M,
                      std::vector<MyMatrix<Tint>> const &ListSHVinp,
                      std::string const &TheMethod, std::ostream &os) {
  Tint eMaxDet = SHORT_GetMaximumDeterminant(M);
  int n = M.cols();
  int nbVect = M.rows();
  std::vector<MyMatrix<Tint>> ListMatrGen =
      SHORT_GetStabilizer<T, Tint, Tgroup>(M, os);
  //
  // Building the set of inequalities
  //
  IteratorBinomial<mpz_class> eIter(nbVect, n - 1);
  std::vector<Face> ListAllFace = eIter.ListAllFace();
  std::vector<MyVector<T>> ListIneq;
  MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
  for (auto const &eFace : ListAllFace) {
    MyMatrix<T> Mred = SelectRow(M_T, eFace);
    if (RankMat(Mred) == n - 1) {
      MyVector<T> eVect(n);
      for (int i = 0; i < n; i++) {
        MyMatrix<T> eLine(n, 1);
        for (int j = 0; j < n; j++)
          eLine(j, 0) = 0;
        eLine(i, 0) = 1;
        MyMatrix<T> Mdet = Concatenate(Mred, eLine);
        T eDet = DeterminantMat(Mdet);
        eVect(i) = eDet;
      }
      MyVector<T> eIneq(n + 1), fIneq(n + 1);
      eIneq(0) = eMaxDet;
      fIneq(0) = eMaxDet;
      for (int i = 0; i < n; i++) {
        eIneq(i + 1) = eVect(i);
        fIneq(i + 1) = -eVect(i);
      }
      ListIneq.push_back(eIneq);
      ListIneq.push_back(fIneq);
    }
  }
  MyMatrix<T> FAC = MatrixFromVectorFamily(ListIneq);
  MyMatrix<T> EXT = cdd::DualDescription(FAC, os);
  std::vector<MyVector<Tint>> ListPt =
      GetListIntegralPoint<T, Tint>(FAC, EXT, os);
  //
  // Breaking into orbits
  //
  std::vector<MyMatrix<Tint>> ListGen;
  for (auto &eGen : ListMatrGen)
    ListGen.push_back(TransposedMat(eGen));
  std::function<MyVector<Tint>(MyVector<Tint> const &, MyMatrix<Tint> const &)>
      TheAct = [](MyVector<Tint> const &x,
                  MyMatrix<Tint> const &M) -> MyVector<Tint> { return M * x; };
  std::vector<MyVector<Tint>> ListRepr =
      OrbitSplittingGeneralized(ListPt, ListGen, TheAct);
  //
  // Building the set of inequalities
  //
  auto IsPresent = [&](MyMatrix<Tint> const &P) -> bool {
    for (auto &P2 : ListSHVinp) {
      std::optional<MyMatrix<Tint>> eResEquiv =
          SHORT_TestEquivalence<T, Tint, Tgroup>(P, P2, os);
      if (eResEquiv)
        return true;
    }
    return false;
  };
  auto PassFacetIsoCheck = [&](MyMatrix<Tint> const &U) -> bool {
    int len = U.rows();
    for (int i = 0; i < len; i++) {
      Face eFace(len);
      for (int j = 0; j < len; j++)
        eFace[j] = 1;
      eFace[i] = 0;
      MyMatrix<Tint> Usel = SelectRow(U, eFace);
      if (!IsPresent(Usel))
        return false;
    }
    return true;
  };
  std::vector<MyMatrix<Tint>> ListSpann;
  auto FuncInsert = [&](MyMatrix<Tint> const &Mnew) -> void {
    if (!PassFacetIsoCheck(Mnew))
      return;
    for (auto &P2 : ListSpann) {
      std::optional<MyMatrix<Tint>> eResEquiv =
          SHORT_TestEquivalence<T, Tint, Tgroup>(Mnew, P2, os);
      if (eResEquiv)
        return;
    }
    ReplyRealizability<T, Tint> eTestRes =
        SHORT_TestRealizabilityShortestFamily<T, Tint, Tgroup>(Mnew, TheMethod,
                                                               os);
    if (eTestRes.reply && eTestRes.replyCone)
      ListSpann.push_back(Mnew);
  };
  for (auto &ePt : ListRepr) {
    MyMatrix<Tint> eLine(n, 1);
    Tint eNorm = 0;
    for (int i = 0; i < n; i++) {
      Tint eVal = ePt(i);
      eNorm += eVal * eVal;
      eLine(i) = eVal;
    }
    if (eNorm > 0) {
      MyMatrix<Tint> eCand = Concatenate(M, eLine);
      FuncInsert(eCand);
    }
  }
  return ListSpann;
}

template <typename Tint>
void WriteListConfigurationShortestVector(
    std::string const &eFile, std::vector<MyMatrix<Tint>> const &ListConf) {
  std::ofstream os(eFile);
  int nbConf = ListConf.size();
  os << nbConf << "\n";
  for (int i = 0; i < nbConf; i++)
    WriteMatrix(os, ListConf[i]);
}

template <typename Tint>
std::vector<MyMatrix<Tint>>
ReadListConfigurationShortestVector(std::string const &eFile) {
  std::ifstream is(eFile);
  int nbConf;
  is >> nbConf;
  std::vector<MyMatrix<Tint>> ListConf(nbConf);
  for (int i = 0; i < nbConf; i++)
    ListConf[i] = ReadMatrix<Tint>(is);
  return ListConf;
}

// Canonicalize the vector under the action of the symmetric group
// and inversion (group size 2^n n!)
//
MyVector<int> CyclicCanonicalization_SymN(MyVector<int> const &fCand,
                                          int const &d) {
  int n = fCand.size();
  std::vector<int> RetCand;
  for (int i = 0; i < n; i++) {
    int aVal = fCand(i);
    int bVal;
    if (aVal < 0) {
      bVal = (-aVal) % d;
    } else {
      bVal = aVal % d;
    }
    int cVal;
    if (2 * bVal > d) {
      cVal = d - bVal;
    } else {
      cVal = bVal;
    }
    RetCand.push_back(cVal);
  }
  int len = fCand.size();
  std::vector<int> eList = SortingPerm<int, int>(RetCand);
  MyVector<int> nCand(len);
  for (int i = 0; i < len; i++) {
    int j = eList[i];
    nCand(i) = RetCand[j];
  }
  return nCand;
}

// Canicalize the vector under symn action + the multiplication by
// non-zero element of Fd (field of the prime d)
MyVector<int> CyclicCanonicalization_SymN_fact(MyVector<int> const &V,
                                               int const &d) {
  MyVector<int> Vcan = CyclicCanonicalization_SymN(V, d);
  for (int mult = 2; mult <= d - 1; mult++) {
    MyVector<int> Vm = mult * V;
    MyVector<int> VmCan = CyclicCanonicalization_SymN(Vm, d);
    if (VmCan < Vcan)
      Vcan = VmCan;
  }
  return Vcan;
}

std::vector<std::vector<int>> SHORT_GetCandidateCyclic_Optimized(int const &n,
                                                                 int const &d) {
  int res = d % 2;
  int MaxVal = (d - res) / 2;
  auto Canonicalization =
      [&d](std::vector<int> const &fCand) -> std::vector<int> {
    std::vector<int> RetCand;
    for (auto &eVal : fCand) {
      int nVal;
      if (2 * eVal > d) {
        nVal = d - eVal;
      } else {
        nVal = eVal;
      }
      RetCand.push_back(nVal);
    }
    int len = fCand.size();
    std::vector<int> eList = SortingPerm<int, int>(RetCand);
    std::vector<int> nCand(len);
    for (int i = 0; i < len; i++) {
      int j = eList[i];
      nCand[i] = RetCand[j];
    }
    return nCand;
  };
  auto IsMinimal = [&d,
                    &Canonicalization](std::vector<int> const &eCand) -> bool {
    for (int mult = 2; mult <= d - 1; mult++) {
      std::vector<int> eProd;
      for (auto &x : eCand) {
        int y = x * mult;
        int res = y % d;
        eProd.push_back(res);
      }
      std::vector<int> NewCand = Canonicalization(eProd);
      if (NewCand < eCand)
        return false;
    }
    return true;
  };
  std::vector<std::vector<int>> ListCand = {{1}};
  for (int iDim = 1; iDim < n; iDim++) {
    std::vector<std::vector<int>> NewListCand;
    for (auto &eCand : ListCand) {
      int LastVal = eCand[iDim - 1];
      for (int i = LastVal; i <= MaxVal; i++) {
        std::vector<int> NewCand = ConcatenateVect(eCand, {i});
        if (IsMinimal(NewCand)) {
          NewListCand.push_back(NewCand);
        }
      }
    }
    ListCand = NewListCand;
#ifdef DEBUG_SHORTEST_CONFIG
    std::cerr << "iDim=" << iDim << " |ListCand|=" << ListCand.size() << "\n";
#endif
  }
  return ListCand;
}

struct PrimeListAllowed {
  int p;
  bool DoWeListFeasible;
  std::vector<MyVector<int>> ListCases;
};

template <typename T>
bool IsMatchingListOfPrimes(std::vector<PrimeListAllowed> const &ListPrime,
                            MyMatrix<T> const &M) {
  std::vector<MyVector<int>> ListClasses = ComputeTranslationClasses<T, int>(M);
  int n = M.rows();
  MyMatrix<T> eInv = Inverse(M);
  auto GetOrder = [&](MyVector<int> const &eV) -> int {
    if (IsZeroVector(eV))
      return -1;
    MyVector<T> V2(n);
    for (int i = 0; i < n; i++) {
      T sum = 0;
      for (int j = 0; j < n; j++)
        sum += eV(j) * eInv(j, i);
      V2(i) = sum;
    }
    int ord = 1;
    while (true) {
      bool IsProductCorr = true;
      for (int i = 0; i < n; i++) {
        T eProd = ord * V2(i);
        if (!IsInteger(eProd))
          IsProductCorr = false;
      }
      if (IsProductCorr)
        return ord;
      ord++;
    }
    return -1;
  };
  for (auto &eV : ListClasses) {
    int ord = GetOrder(eV);
    for (auto &eCaseP : ListPrime) {
      if (ord == eCaseP.p) {
        auto IsCorrectClass = [&](MyVector<int> const &W) -> bool {
          for (auto &W2 : eCaseP.ListCases) {
            if (W2 == W) {
              if (eCaseP.DoWeListFeasible)
                return true;
              else
                return false;
            }
          }
          if (eCaseP.DoWeListFeasible)
            return false;
          else
            return true;
        };
        MyVector<int> Vtest(n);
        for (int i = 0; i < n; i++) {
          T sum = 0;
          for (int j = 0; j < n; j++)
            sum += eV(j) * eInv(j, i);
          sum *= ord;
          if (!IsInteger(sum)) {
            std::cerr << "The sum should be integral\n";
            throw TerminalException{1};
          }
          int sum_i = UniversalScalarConversion<int, T>(sum);
          Vtest(i) = sum_i;
        }
        MyVector<int> VtestCan = CyclicCanonicalization_SymN_fact(Vtest, ord);
        if (!IsCorrectClass(VtestCan))
          return false;
      }
    }
  }
  return true;
}

template <typename T, typename Tint, typename Tgroup>
std::pair<std::vector<MyMatrix<Tint>>, std::vector<int>>
SHORT_ReduceByIsomorphism(std::vector<MyMatrix<Tint>> const &ListSHV,
                          std::ostream &os) {
  std::map<SHVinvariant<T, Tint>, std::vector<int>> TheMap;
  std::vector<MyMatrix<Tint>> ListRet;
  int siz = 0;
  //  SHVinvariant<T,Tint> eInv;
  auto FuncInsert = [&](MyMatrix<Tint> const &eSpann) -> int {
    SHVinvariant<T, Tint> eInv = SHORT_Invariant<T, Tint>(eSpann, os);
    auto search = TheMap.find(eInv);
    if (search == TheMap.end()) {
      TheMap[eInv] = {siz};
      ListRet.push_back(eSpann);
      int pos = siz;
      siz++;
      return pos;
    }
    for (auto &iSpann : TheMap[eInv]) {
      MyMatrix<Tint> fSpann = ListRet[iSpann];
      std::optional<MyMatrix<Tint>> RecTest =
          SHORT_TestEquivalence<T, Tint, Tgroup>(eSpann, fSpann, os);
      if (RecTest)
        return iSpann;
    }
    TheMap[eInv].push_back(siz);
    ListRet.push_back(eSpann);
    int pos = siz;
    siz++;
    return pos;
  };
  int nbConf = ListSHV.size();
  int pos = 0;
  std::vector<int> ListIdx(nbConf, -1);
  for (auto &eSHV : ListSHV) {
    int idx = FuncInsert(eSHV);
    ListIdx[pos] = idx;
    pos++;
  }
  return {std::move(ListRet), std::move(ListIdx)};
}

// clang-format off
#endif  // SRC_SHORT_SHORT_SHORTESTCONFIG_H_
// clang-format on
