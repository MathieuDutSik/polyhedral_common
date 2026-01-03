// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SHORT_SHORT_SHORTESTCONFIG_H_
#define SRC_SHORT_SHORT_SHORTESTCONFIG_H_

// clang-format off
#include "SHORT_ShortestConfigStabEquiCanInv.h"
#include "Shvec_exact.h"
#include "COMB_Combinatorics.h"
#include "InvariantVectorFamily.h"
#include "LatticeDefinitions.h"
#include "MAT_MatrixInt.h"
#include "MAT_functions.h"
#include "MatrixGroup.h"
#include "POLY_LinearProgramming.h"
#ifdef USE_GLPK
#include "POLY_LinearProgramming_GLPK.h"
#endif
#include "POLY_PolytopeInt.h"
#include "Positivity.h"
#include "Tspace_InvariantForm.h"
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
    auto get_icol=[&]() -> int {
      for (int iCol = 0; iCol < nbCol; iCol++)
        if (eRow(iCol) != 0)
          return iCol;
      return -1;
    };
    int iColFound = get_icol();
#ifdef SANITY_CHECK_SHORTEST_CONFIG
    if (iColFound == -1) {
      std::cerr << "SHORT: Bug in SHORT_CleanAntipodality\n";
      std::cerr << "SHORT: The vector is 0 which is not allowed\n";
      throw TerminalException{1};
    }
#endif
    if (eRow(iColFound) < 0) {
      eRow = -eRow;
    }
    setVect.insert(eRow);
  }
  std::vector<MyVector<Tint>> ListVect;
  for (auto &eV : setVect)
    ListVect.push_back(eV);
  return MatrixFromVectorFamily(ListVect);
}

template <typename T, typename Tint> struct ReplyRealizability {
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
T get_single_norm(MyMatrix<T> const& eMat, std::vector<MyVector<Tint>> const& ListVect) {
  bool IsFirst = true;
  T singleNorm;
  for (auto &eVect : ListVect) {
    T eNorm = EvaluationQuadForm(eMat, eVect);
    if (IsFirst) {
      singleNorm = eNorm;
      IsFirst = false;
    } else {
      if (eNorm != singleNorm) {
        std::cerr << "SHORT: The norms are not consistent\n";
        throw TerminalException{1};
      }
    }
  }
  return singleNorm;
}


/*
  It is an initial set of vector used for linear programming purposes.
  So, for the shortest vector but also for other searches.
  It is supposed to be full-dimensional in T-space, even if ListVect is empty.
 */
template <typename Tint>
std::unordered_set<MyVector<Tint>> get_initial_vector_test_s(int const& n, std::vector<MyVector<Tint>> const &ListVect, [[maybe_unused]] std::ostream &os) {
  auto f_can=[&](MyVector<Tint> const& V) -> MyVector<Tint> {
    for (int i=0; i<n; i++) {
      if (V(i) != 0) {
        if (V(i) > 0) {
          return V;
        }
        if (V(i) < 0) {
          return -V;
        }
      }
    }
    // Zero vector, nothing to be done!
    return V;
  };
  std::unordered_set<MyVector<Tint>> ListForbiddenVector;
  ListForbiddenVector.insert(ZeroVector<Tint>(n));
  for (auto & eV: ListVect) {
    MyVector<Tint> eVcan = f_can(eV);
    ListForbiddenVector.insert(eVcan);
  }
  std::vector<MyVector<Tint>> TestVect1;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2; j++) {
      int shift = -1 + 2 * j;
      MyVector<Tint> tVect1 = ZeroVector<Tint>(n);
      tVect1[i] += shift;
      TestVect1.push_back(tVect1);
    }
  }
  std::vector<MyVector<Tint>> TestVect2;
  for (size_t i_v = 0; i_v < TestVect1.size(); i_v++) {
    MyVector<Tint> tVect1 = TestVect1[i_v];
    for (size_t j_v = i_v+1; j_v < TestVect1.size(); j_v++) {
      MyVector<Tint> tVect2 = tVect1 + TestVect1[j_v];
      TestVect2.push_back(tVect2);
    }
    TestVect2.push_back(tVect1);
  }
  std::unordered_set<MyVector<Tint>> TestVect3;
  for (auto & eVect: ListVect) {
    for (auto &tVect2: TestVect2) {
      MyVector<Tint> tVect3_A = eVect + tVect2;
      MyVector<Tint> tVect3 = f_can(tVect3_A);
      if (!IsZeroVector(tVect3) && ListForbiddenVector.count(tVect3) == 0) {
        TestVect3.insert(tVect3);
      }
    }
  }
  return TestVect3;
}

template <typename Tint>
std::vector<MyVector<Tint>> get_initial_vector_test_v(int const& n, std::vector<MyVector<Tint>> const &ListVect, std::ostream &os) {
  std::unordered_set<MyVector<Tint>> vectors_s = get_initial_vector_test_s(n, ListVect, os);
  std::vector<MyVector<Tint>> vectors_v;
  for (auto & eV: vectors_s) {
    vectors_v.push_back(eV);
  }
  return vectors_v;
}

template <typename T, typename Tint>
std::vector<MyMatrix<T>> get_spec_shv_basis(std::vector<MyMatrix<T>> const& ListMat, std::vector<MyVector<Tint>> const &ListVect) {
  int nbVect = ListVect.size();
  int nbMat = ListMat.size();
  MyMatrix<T> MatrixValuesDiff(nbVect - 1, nbMat);
  for (int iMat = 0; iMat < nbMat; iMat++) {
    MyMatrix<T> const& eMat = ListMat[iMat];
    MyVector<Tint> const& V1 = ListVect[0];
    T val1 = EvaluationQuadForm(eMat, V1);
    for (int iVect = 0; iVect < nbVect - 1; iVect++) {
      MyVector<Tint> const& V2 = ListVect[iVect+1];
      T val2 = EvaluationQuadForm(eMat, V2);
      MatrixValuesDiff(iVect, iMat) = val2 - val1;
    }
  }
  MyMatrix<T> NSP = NullspaceTrMat(MatrixValuesDiff);
  int dimSpa = NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    MyVector<T> eRow = GetMatrixRow(NSP, iDim);
    TheBasis[iDim] = GetMatrixFromBasis(ListMat, eRow);
  }
  return TheBasis;
}



template <typename T, typename Tint>
ReplyRealizability<T, Tint> SHORT_TestRealizabilityShortestFamilyEquivariant(
    std::vector<MyVector<Tint>> const &ListVect,
    std::vector<MyMatrix<T>> const &ListMat, bool const &NoExtension,
    std::ostream &os) {
  ReplyRealizability<T, Tint> eRes;
#ifdef SANITY_CHECK_SHORTEST_CONFIG
  if (ListVect.size() == 0) {
    std::cerr << "SHORT: ListVect should not be empty\n";
    throw TerminalException{1};
  }
#endif
  int n = ListVect[0].size();
#ifdef DEBUG_SHORTEST_CONFIG
  os << "SHORT: SHORT_TestRealizabilityShortestFamilyEquivariant, n=" << n << " nbVect=" << ListVect.size() << " nbMat=" << ListMat.size() << "\n";
#endif
  std::unordered_set<MyVector<Tint>> SetVectTot;
  for (auto &eVect : ListVect) {
    SetVectTot.insert(eVect);
    SetVectTot.insert(-eVect);
  }
  std::vector<MyMatrix<T>> TheBasis = get_spec_shv_basis(ListMat, ListVect);
  int dimSpa = TheBasis.size();
  //
  // Forming the vector family. Start with something and extend later on.
  //
  std::unordered_set<MyVector<Tint>> TheFamilyVect = get_initial_vector_test_s(n, ListVect, os);
#ifdef DEBUG_SHORTEST_CONFIG
  os << "SHORT: |TheFamilyVect|=" << TheFamilyVect.size() << "\n";
#endif
  auto GetListIneq = [&]() -> MyMatrix<T> {
    int siz = TheFamilyVect.size();
    MyMatrix<T> RetMat(siz, dimSpa + 1);
    int idx = 0;
    for (auto &eVect : TheFamilyVect) {
      MyVector<T> vVect = SHORT_GetIneq_Tspace<T, Tint>(TheBasis, eVect);
      AssignMatrixRow(RetMat, idx, vVect);
      idx++;
    }
    return RetMat;
  };
  MyVector<T> ToBeMinimized(dimSpa + 1);
  ToBeMinimized(0) = 0;
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    T scal = EvaluationQuadForm(TheBasis[iDim], ListVect[0]);
    ToBeMinimized(iDim + 1) = scal;
  }
  if (IsZeroVector(ToBeMinimized)) {
    eRes.reply = false;
    eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
    os << "SHORT: RETURN case 1\n";
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
#ifdef SANITY_CHECK_SHORTEST_CONFIG
      if (eOptimalPrev > eOptimal) {
        std::cerr << "SHORT: Optimal values should be increasing\n";
        throw TerminalException{1};
      }
#endif
      eOptimalPrev = eOptimal;
    }
    nbIter++;
#ifdef DEBUG_SHORTEST_CONFIG
    os << "SHORT: nbIter=" << nbIter << "\n";
#endif
    MyMatrix<T> ListIneq = GetListIneq();
#ifdef SANITY_CHECK_SHORTEST_CONFIG
    for (auto &eVect : SetVectTot) {
      if (TheFamilyVect.count(eVect) == 1) {
        std::cerr << "SHORT: We find one of SetVectTot in TheFamilyVect\n";
        throw TerminalException{1};
      }
    }
#endif
    int sizFamVect = TheFamilyVect.size();
    for (int i = 0; i < sizFamVect; i++) {
      MyVector<T> rVect = GetMatrixRow(ListIneq, i);
      if (rVect == ZerVect) {
        eRes.reply = false;
        eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
        os << "SHORT: RETURN case 2\n";
#endif
        return eRes;
      }
    }
    MyMatrix<T> SetIneq = SortUnicizeMatrix(ListIneq);
    LpSolution<T> eSol = CDD_LinearProgramming(SetIneq, ToBeMinimized, os);
#ifdef DEBUG_SHORTEST_CONFIG
    os << "SHORT: DualDefined=" << eSol.DualDefined << " PrimalDefined=" << eSol.PrimalDefined << "\n";
#endif
    if (!eSol.PrimalDefined && eSol.DualDefined) {
#ifdef DEBUG_SHORTEST_CONFIG
      os << "SHORT: DualDefined but not primal defined\n";
#endif
      int nbIneqSet = SetIneq.size();
      MyVector<T> SumIneq = ZeroVector<T>(1 + nbIneqSet);
      for (int i = 0; i < nbIneqSet; i++) {
        MyMatrix<T> eRow = GetMatrixRow(SetIneq, i);
        SumIneq += eSol.DualSolution(i) * eRow;
      }
      MyVector<T> SumIneqRed = ZeroVector<T>(nbIneqSet);
      for (int i = 0; i < nbIneqSet; i++) {
        SumIneqRed(i) = SumIneq(i + 1);
      }
      if (SumIneq(0) < 0 && SumIneqRed == ZeroVector<T>(nbIneqSet)) {
        eRes.reply = false;
        eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
        os << "SHORT: RETURN case 3\n";
#endif
        return eRes;
      }
      std::cerr << "SHORT: It seems we have a big problem here. Please correct\n";
      throw TerminalException{1};
    } else if (eSol.PrimalDefined && !eSol.DualDefined) {
#ifdef DEBUG_SHORTEST_CONFIG
      os << "SHORT: PrimalDefined but not dual defined\n";
#endif
      // This is the primal_direction case
      MyVector<Tint> eVect = ListVect[0];
      MyVector<Tint> eVectB = 2 * eVect;
#ifdef DEBUG_SHORTEST_CONFIG
      os << "SHORT: Inserting from PrimalDefined but not dual defined eVectB=";
      WriteVectorNoDim(os, eVectB);
#endif
      TheFamilyVect.insert(eVectB);
    } else {
#ifdef DEBUG_SHORTEST_CONFIG
      os << "SHORT: We have optimal value\n";
#endif
#ifdef SANITY_CHECK_SHORTEST_CONFIG
      if (!eSol.PrimalDefined || !eSol.DualDefined) {
        std::cerr << "SHORT: We have a real problem to solve. Please debug\n";
        throw TerminalException{1};
      }
#endif
      if (nbIter == 1) {
        eOptimalPrev = eOptimal;
      }
      eOptimal = eSol.OptimalValue;
#ifdef DEBUG_SHORTEST_CONFIG
      os << "SHORT: eOptimal=" << eOptimal << " eOptimalPrev=" << eOptimalPrev << "\n";
#endif
      if (eOptimal > 1) {
        eRes.reply = false;
        eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
        os << "SHORT: RETURN case 4\n";
#endif
        return eRes;
      }
      if (eOptimal >= 1 && NoExtension) {
        eRes.reply = false;
        eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
        os << "SHORT: RETURN case 5\n";
#endif
        return eRes;
      }
      MyVector<T> eVectEmb = eSol.DirectSolution;
      MyVector<T> rVect = GetDirectSolutionExt(eSol);
      MyMatrix<T> eMatSec = GetMatrixFromBasis(TheBasis, eVectEmb);
#ifdef DEBUG_SHORTEST_CONFIG
      os << "SHORT: eMatSec=\n";
      WriteMatrix(os, eMatSec);
#endif
      //
      // Some checks
      //
#ifdef SANITY_CHECK_SHORTEST_CONFIG
      int idx = 0;
      for (auto &eVect : TheFamilyVect) {
        MyVector<T> eRow = GetMatrixRow(ListIneq, idx);
        if (eRow != ZerVect) {
          T eNorm = EvaluationQuadForm(eMatSec, eVect);
          if (eNorm < 1) {
            std::cerr
                << "SHORT: This should not happen if the matrix was well defined\n";
            throw TerminalException{1};
          }
          T eScal1 = ScalarProduct(rVect, eRow) + 1;
          if (eScal1 != eNorm) {
            std::cerr << "SHORT: We should have eScal1 != eNorm\n";
            throw TerminalException{1};
          }
        }
        idx++;
      }
#endif
      T singleNorm = get_single_norm(eMatSec, ListVect);
      T CritNorm = std::min(T(1), singleNorm);
#ifdef SANITY_CHECK_SHORTEST_CONFIG
      if (CritNorm != eOptimal) {
        std::cerr << "SHORT: CritNorm=" << CritNorm << " eOptimal=" << eOptimal
                  << "\n";
        std::cerr << "SHORT: The value of eOptimal has not been well set\n";
        throw TerminalException{1};
      }
#endif
      bool testPosDef = IsPositiveDefinite(eMatSec, os);
#ifdef DEBUG_SHORTEST_CONFIG
      os << "SHORT: testPosDef=" << testPosDef << "\n";
#endif
      if (testPosDef) {
#ifdef DEBUG_SHORTEST_CONFIG
        os << "SHORT: Before computation of T_ShortestVector\n";
        WriteMatrix(os, eMatSec);
#endif
        Tshortest<T, Tint> rec_shv = T_ShortestVector<T, Tint>(eMatSec, os);
#ifdef DEBUG_SHORTEST_CONFIG
        os << "SHORT: After computation of T_ShortestVector\n";
#endif
        int nbRow = rec_shv.SHV.rows();
        std::unordered_set<MyVector<Tint>> SetSHV;
        for (int iRow = 0; iRow < nbRow; iRow++) {
          MyVector<Tint> eVect = GetMatrixRow(rec_shv.SHV, iRow);
          SetSHV.insert(eVect);
        }
#ifdef DEBUG_SHORTEST_CONFIG
        os << "SHORT: SetSHV=\n";
        WriteMatrix(os, MatrixFromUnorderedSetFamily(SetSHV));
        os << "SHORT: SetVectTot=\n";
        WriteMatrix(os, MatrixFromUnorderedSetFamily(SetVectTot));
#endif
        bool testEqua = SetSHV == SetVectTot;
#ifdef DEBUG_SHORTEST_CONFIG
        os << "SHORT: testEqua=" << testEqua << "\n";
#endif
        if (testEqua) {
          eRes.reply = true;
          eRes.replyCone = true;
          eRes.SHV = rec_shv.SHV;
          eRes.SHVclean = SHORT_CleanAntipodality(rec_shv.SHV);
          eRes.eMat = eMatSec;
#ifdef DEBUG_SHORTEST_CONFIG
          os << "SHORT: RETURN case 6\n";
#endif
          return eRes;
        } else {
          std::vector<MyVector<Tint>> DiffNew;
          for (auto &eVect : SetSHV) {
            if (SetVectTot.count(eVect) == 0 && TheFamilyVect.count(eVect) == 0) {
              DiffNew.push_back(eVect);
            }
          }
          if (DiffNew.size() > 0) {
            for (auto &eVect : DiffNew) {
#ifdef DEBUG_SHORTEST_CONFIG
              os << "SHORT: Inserting from DiffNew eVect=";
              WriteVectorNoDim(os, eVect);
#endif
              TheFamilyVect.insert(eVect);
            }
          } else {
#ifdef DEBUG_SHORTEST_CONFIG
            MyMatrix<Tint> M1 = MatrixFromUnorderedSetFamily(SetSHV);
            os << "SHORT: M1=\n";
            WriteMatrix(os, M1);
            MyMatrix<Tint> M2 = MatrixFromUnorderedSetFamily(SetVectTot);
            os << "SHORT: M2=\n";
            WriteMatrix(os, M2);
#endif
            if (IsSubsetUnorderedSet(SetSHV, SetVectTot)) {
              eRes.reply = false;
              eRes.replyCone = true;
              eRes.SHV = rec_shv.SHV;
              eRes.SHVclean = SHORT_CleanAntipodality(rec_shv.SHV);
              eRes.eMat = eMatSec;
#ifdef DEBUG_SHORTEST_CONFIG
              os << "SHORT: RETURN case 7\n";
#endif
              return eRes;
            } else {
              eRes.reply = false;
              eRes.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
              os << "SHORT: RETURN case 8\n";
#endif
              return eRes;
            }
          }
        }
      } else {
        auto f_insert=[&](MyVector<Tint> const& eVect) -> void {
#ifdef SANITY_CHECK_SHORTEST_CONFIG
          if (TheFamilyVect.count(eVect) == 1) {
            std::cerr << "SHORT: We have a clear error here\n";
            throw TerminalException{1};
          }
#endif
#ifdef DEBUG_SHORTEST_CONFIG
          os << "SHORT: Inserting from GetShortVector eVect=";
          WriteVectorNoDim(os, eVect);
#endif
          TheFamilyVect.insert(eVect);
        };
        if (RankMat(eMatSec) < n) {
          std::vector<MyVector<Tint>> list_v3 =
              GetShortVectorDegenerate<T, Tint>(eMatSec, os);
          for (auto& eVect3 : list_v3) {
            if (SetVectTot.count(eVect3) == 1) { // This can happen.
              eVect3 *= 2;
            }
            f_insert(eVect3);
          }
        } else {
          bool StrictIneq = true;
          MyVector<Tint> eVect = GetShortIntegralVector<T, Tint>(
              eMatSec, CritNorm, StrictIneq, os);
          f_insert(eVect);
        }
      }
    }
  }
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
  std::cerr
      << "SHORT: The dimension is too large for use to give the kissing number upper bound\n";
  throw TerminalException{1};
}

template <typename Tint> struct SHVreduced {
  MyMatrix<Tint> SHVred;
  MyMatrix<Tint> ReductionMatrix;
};

template <typename Tint>
SHVreduced<Tint> SHORT_GetLLLreduction_Kernel(MyMatrix<Tint> const &eSHV,
                                              std::ostream &os) {
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
  LLLreduction<Tfield, Tint> res = LLLreducedBasis<Tfield, Tint>(TheGram, os);
  MyMatrix<Tfield> TheRemainder = res.GramMatRed;
  MyMatrix<Tint> TheTrans = res.Pmat;
  MyMatrix<Tint> Pmat = TransposedMat(TheTrans);
  MyMatrix<Tint> eSHVred = eSHV * Pmat;
#ifdef SANITY_CHECK_SHORTEST_CONFIG
  if (!TestEqualityMatrix(GetGram(eSHVred), TheRemainder)) {
    std::cerr << "SHORT: Matrix error somewhere\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(eSHVred), std::move(Pmat)};
}

template <typename T> std::vector<MyMatrix<T>> StandardSymmetricBasis(int n) {
  std::vector<MyMatrix<T>> ListMat;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
      eMat(i, j) = 1;
      eMat(j, i) = 1;
      ListMat.push_back(eMat);
    }
  }
  return ListMat;
}

template <typename T, typename Tint, typename Tgroup>
ReplyRealizability<T, Tint>
SHORT_TestRealizabilityShortestFamily(MyMatrix<Tint> const &Minput,
                                      std::ostream &os) {
  SHVreduced<Tint> RecLLL = SHORT_GetLLLreduction_Kernel(Minput, os);
  MyMatrix<Tint> M = RecLLL.SHVred;
  int n = M.cols();
  std::vector<MyMatrix<Tint>> ListMatrGen =
      SHORT_GetStabilizer<T, Tint, Tgroup>(M, os);
  std::vector<MyMatrix<T>> StdBasis = StandardSymmetricBasis<T>(n);
  std::vector<MyMatrix<T>> ListGen_T = UniversalStdVectorMatrixConversion<T,Tint>(ListMatrGen);
#ifdef DEBUG_SHORTEST_CONFIG
  os << "SHORT: Before BasisInvariantForm\n";
#endif
  std::vector<MyMatrix<T>> ListMat = BasisInvariantForm(n, ListGen_T, os);
#ifdef DEBUG_SHORTEST_CONFIG
  os << "SHORT: After BasisInvariantForm\n";
#endif
  int nbVect = M.rows();
  std::vector<MyVector<T>> ListRankOne;
  std::vector<MyVector<Tint>> ListVectWork;
  auto FuncInsertVector = [&](MyVector<Tint> const &eRow) -> void {
    if (GetPositionAntipodal(ListVectWork, eRow) == -1) {
      ListVectWork.push_back(eRow);
      ListRankOne.push_back(SHORT_GetIneq_Tspace<T, Tint>(StdBasis, eRow));
    }
  };
#ifdef DEBUG_SHORTEST_CONFIG
  os << "SHORT: nbVect=" << nbVect << "\n";
#endif
  for (int iRow = 0; iRow < nbVect; iRow++) {
    MyVector<Tint> eRow = GetMatrixRow(M, iRow);
    FuncInsertVector(eRow);
  }
  int InitialSize = ListVectWork.size();
  MyMatrix<T> ListRankOne_mat = MatrixFromVectorFamily(ListRankOne);
  auto IsNeededInsert = [&](MyVector<Tint> const &eVect) -> bool {
    MyVector<T> VectImg = SHORT_GetIneq_Tspace<T, Tint>(StdBasis, eVect);
    std::optional<MyVector<T>> opt = SolutionMat(ListRankOne_mat, VectImg);
    return opt.has_value();
  };
  int TheRank = RankMat(ListRankOne_mat);
  bool NoExtension = false;
  if (PositionVect({n, n + 1, n + 2}, TheRank) != -1) {
    NoExtension = true;
  }
  while (true) {
    ReplyRealizability<T, Tint> RecTest;
    if (InitialSize > KissingNumberUpperBound(n)) {
      RecTest.reply = false;
      RecTest.replyCone = false;
#ifdef DEBUG_SHORTEST_CONFIG
      std::cerr << "SHORT: RETURN case 9\n";
#endif
      return RecTest;
    }
    RecTest = SHORT_TestRealizabilityShortestFamilyEquivariant<T, Tint>(
        ListVectWork, ListMat, NoExtension, os);
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
    std::cerr << "SHORT: SizPrev=" << SizPrev << " SizAfter=" << SizAfter << "\n";
#endif
    if (SizAfter == SizPrev) {
      RecTest.reply = false;
      RecTest.replyCone = false;
      return RecTest;
    }
  }
}

template <typename T, typename Tint>
size_t SHORT_Invariant(MyMatrix<Tint> const &eSpann, std::ostream &os) {
  ShortIso<T, Tint> eShIso = SHORT_GetInformation<T, Tint>(eSpann, os);
  size_t seed = 146;
  return GetInvariantGramShortest(eShIso.GramMat, eShIso.SHVdisc, seed, os);
}

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
      if (eDetA > eRet) {
        eRet = eDetA;
      }
    }
  }
  return eRet;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
SHORT_SpannSimplicial(MyMatrix<Tint> const &M,
                      std::vector<MyMatrix<Tint>> const &ListSHVinp,
                      std::ostream &os) {
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
      if (eResEquiv) {
        return;
      }
    }
    ReplyRealizability<T, Tint> eTestRes =
        SHORT_TestRealizabilityShortestFamily<T, Tint, Tgroup>(Mnew,
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

// clang-format off
#endif  // SRC_SHORT_SHORT_SHORTESTCONFIG_H_
// clang-format on
