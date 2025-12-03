// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_LORENTZIAN_PERFECT_H_
#define SRC_LORENTZIAN_LORENTZIAN_PERFECT_H_

// clang-format off
#include "boost_serialization.h"
#include "FundamentalDelaunay.h"
#include "PolytopeEquiStabInt.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_AdjacencyScheme.h"
#include "SystemNamelist.h"
#include "fractions.h"
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <utility>
// clang-format on

#ifdef DEBUG
#define DEBUG_LORENTZIAN_PERFECT
#define DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LORENTZIAN_PERFECT
#define SANITY_CHECK_LORENTZIAN_FIND_POSITIVE_VECTORS
#endif

#ifdef TIMINGS
#define TIMINGS_LORENTZIAN_PERFECT
#endif

#ifdef GAPDEBUGOUT
#define GAPDEBUGOUT_LORENTZIAN_FIND_POSITIVE_VECTORS
#endif

static const int LORENTZIAN_PERFECT_OPTION_ISOTROP = 23;
static const int LORENTZIAN_PERFECT_OPTION_TOTAL = 47;

std::string GetNatureOption(int const &TheOption) {
  if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
    return "isotropic";
  }
  if (TheOption == LORENTZIAN_PERFECT_OPTION_TOTAL) {
    return "total";
  }
  std::cerr << "LORPERF: GetNatureOption failure with TheOption=" << TheOption << "\n";
  throw TerminalException{1};
}

template<typename T>
void check_correctness_lorentzian_perfect(MyMatrix<T> const& G, std::ostream& os) {
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(G, os);
  if (DiagInfo.nbZero != 0) {
    std::cerr << "LORPERF: matrix has non-zero kernel";
    throw TerminalException{1};
  }
  int nbPlus = DiagInfo.nbPlus;
  if (nbPlus != 1) {
    std::cerr << "LORPERF: The signature should be n_plus=1, n_minus=n\n";
    throw TerminalException{1};
  }
}


/*
  We search for the solution such that 0 < x * L * eVect <= MaxScal
  and x * L * x >= 0.
  with the following constraints:
    * if OnlyShortest = true, we select the found vector with
      x * L * eVect being of minimal values.
    * if TheOption is set to ISOTROP, then we limit ourseleves to the
      vectors such that x * L * x = 0.
  --
  How do we solve that problem.
  We iterate over the possible scalar products x * L * eVect
  which we write as
  x = eSol + z U
  That gets us
  norm = x * L * x
       = (eSol + z U) L (U^T z^T + eSol^T)
       = eSol L eSol^T + 2 z U L eSol^T + z U L U^T z^T
          with G = U L U^T
       = eSol L eSol^T + 2 z U L eSol^T + z G z^T
          with w = eSol L U^T G^{-1}
       = eSol L eSol^T + 2 z G w^T + z G z^T
       = eSol L eSol^T - w G w^T + (z + w) G (z + w)^T
       = C + (z + w) G (z + w)^T
          with C = eSol L eSol^T - w G w^T

  The middle point is going to be
    mid = eVect * beta
    with mid * L * eVect = eVal * TheRec.gcd
    That gets us   beta eNorm = eVal * TheRec.gcd
    beta = eVal TheRec.gcd / eNorm
    mid * L * mid = beta * beta * eNorm
    -------------------------------------------------
    The problem we have to solve is to enumerate the point x s.t.
    x.v <= C
    Q[x] >= 0 for Q a form of signature (1,n).
 */
template <typename T, typename Tint>
std::vector<MyVector<Tint>> LORENTZ_FindPositiveVectorsKernel(
    MyMatrix<T> const &LorMat, MyVector<T> const &eVect, T const &MaxScal,
    int const &TheOption, bool const &OnlyShortest, std::ostream &os) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
  MicrosecondTime time;
#endif
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: beginning\n";
  os << "LORPERF: LORENTZ_FindPositiveVectors: OnlyShortest=" << OnlyShortest
     << "\n";
  os << "LORPERF: LORENTZ_FindPositiveVectors: TheOption=" << TheOption << "\n";
  os << "LORPERF: LORENTZ_FindPositiveVectors: det(LorMat)="
     << DeterminantMat(LorMat) << "\n";
  os << "LORPERF: LORENTZ_FindPositiveVectors: LorMat=\n";
  WriteMatrix(os, LorMat);
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(LorMat, os);
  os << "LORPERF: LORENTZ_FindPositiveVectors: nbPlus=" << DiagInfo.nbPlus
     << " nbZero=" << DiagInfo.nbZero << " nbMinus=" << DiagInfo.nbMinus
     << "\n";
#endif
#ifdef GAPDEBUGOUT_LORENTZIAN_FIND_POSITIVE_VECTORS
  auto f_save = [&](std::vector<MyVector<Tint>> const &TotalList) -> void {
    std::string FileSave =
        FindAvailableFileFromPrefix("LORENTZ_FindPositiveVectors");
    MyMatrix<Tint> TotalListB = MatrixFromVectorFamily(TotalList);
    std::ofstream osf(FileSave);
    osf << "return rec(LorMat:=" << StringMatrixGAP(LorMat)
        << ", eVect:=" << StringVectorGAP(eVect) << ", MaxScal:=" << MaxScal
        << ", TheOption:=\"" << GetNatureOption(TheOption) << "\""
        << ", OnlyShortest:=" << GAP_logical(OnlyShortest)
        << ", TotalList:=" << StringMatrixGAP(TotalListB) << ");\n";
  };
#endif
  T eNorm = EvaluationQuadForm(LorMat, eVect);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: eNorm=" << eNorm << "\n";
#endif
#ifdef SANITY_CHECK_LORENTZIAN_FIND_POSITIVE_VECTORS
  if (MaxScal <= 0) {
    std::cerr << "LORPERF: MaxScal=" << MaxScal << " but we should have MaxScal > 0\n";
    throw TerminalException{1};
  }
  if (TheOption != LORENTZIAN_PERFECT_OPTION_ISOTROP &&
      TheOption != LORENTZIAN_PERFECT_OPTION_TOTAL) {
    std::cerr << "LORPERF: TheOption=" << TheOption << "\n";
    throw TerminalException{1};
  }
  if (eNorm <= 0) {
    std::cerr << "LORPERF: eNorm=" << eNorm << " but we should have eNorm > 0\n";
    throw TerminalException{1};
  }
  if (!IsIntegralVector(eVect)) {
    std::cerr << "LORPERF: eVect=\n";
    WriteVector(std::cerr, eVect);
    std::cerr << "LORPERF: eVect should be integral\n";
    throw TerminalException{1};
  }
  if (!IsIntegralMatrix(LorMat)) {
    std::cerr << "LORPERF: LorMat=\n";
    WriteMatrix(std::cerr, LorMat);
    std::cerr << "LORPERF: LorMat should be integral\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 1\n";
#endif
  MyVector<Tint> eVect_tint = UniversalVectorConversion<Tint, T>(eVect);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 2\n";
#endif
  MyVector<T> eVect_LorMat = LorMat * eVect;
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 3        eVect="
     << StringVector(eVect) << "\n";
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 3 eVect_LorMat="
     << StringVector(eVect_LorMat) << "\n";
#endif
  MyVector<Tint> eVect_LorMat_tint =
      UniversalVectorConversion<Tint, T>(eVect_LorMat);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 5 eVect_LorMat_tint="
     << StringVector(eVect_LorMat_tint) << "\n";
#endif
  MyMatrix<T> UbasisPre_T = NullspaceIntVect(eVect_LorMat);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 6\n";
#endif
  MyMatrix<T> Ubasis_T = SublatticeBasisReduction(UbasisPre_T, os);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 7\n";
#endif
  MyMatrix<Tint> Ubasis = UniversalMatrixConversion<Tint, T>(Ubasis_T);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 8\n";
#endif
  MyMatrix<T> GramMat = -Ubasis_T * LorMat * Ubasis_T.transpose();
#ifdef SANITY_CHECK_LORENTZIAN_FIND_POSITIVE_VECTORS
  if (!IsPositiveDefinite(GramMat, os)) {
    std::cerr << "LORPERF: GramMat=\n";
    WriteMatrix(std::cerr, GramMat);
    std::cerr << "LORPERF: GramMat should be positive definite\n";
    throw TerminalException{1};
  }
#endif
  CVPSolver<T, Tint> solver(GramMat, os);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 9\n";
#endif
  GCD_dot<Tint> TheRec_pre = ComputeGcdDot(eVect_LorMat_tint);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 10\n";
  os << "LORPERF: LORENTZ_FindPositiveVectors: TheRec_pre.V="
     << StringVector(TheRec_pre.V) << " gcd=" << TheRec_pre.gcd << "\n";
#endif
  GCD_dot<Tint> TheRec = PositivityNormalizeGcdDot(TheRec_pre);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 11 eNorm=" << eNorm << "\n";
  os << "LORPERF: LORENTZ_FindPositiveVectors: TheRec.V="
     << StringVector(TheRec.V) << " gcd=" << TheRec.gcd << "\n";
#endif
  std::vector<MyVector<Tint>> TotalListSol;
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
  os << "LORPERF: LORENTZ_FindPositiveVectors: step 12\n";
  int iter_findpos = 0;
#endif
  T alpha_basic = TheRec.gcd / eNorm;
  T const& beta_basic = alpha_basic;
  MyVector<Tint> const& eBasSol_basic = TheRec.V;
  MyVector<T> eBasSol_T_basic = UniversalVectorConversion<T, Tint>(eBasSol_basic);
  T scal1_basic = eBasSol_T_basic.dot(eVect_LorMat);
  T scal2_basic = alpha_basic * eVect.dot(eVect_LorMat);
  MyVector<T> eTrans_basic = alpha_basic * eVect - eBasSol_T_basic;
  std::optional<MyVector<T>> opt_basic = SolutionMat(Ubasis_T, eTrans_basic);
  MyVector<T> eSol_basic = unfold_opt(opt_basic, "Getting eSol_basic");
  T eSquareDist_basic = beta_basic * beta_basic * eNorm;
  Tint eVal = 1;
  while (true) {
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
    T scal_expe1 = eVal * TheRec.gcd;
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 1 eVal=" << eVal
       << " scal_expe=" << scal_expe1 << " MaxScal=" << MaxScal
       << " iter_findpos=" << iter_findpos << "\n";
    iter_findpos += 1;
#endif
    MyVector<Tint> eBasSol = eVal * TheRec.V; // A solution of
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 2\n";
#endif
    MyVector<T> eBasSol_T = UniversalVectorConversion<T, Tint>(eBasSol);
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 3\n";
#endif
    T alpha = eVal * alpha_basic;
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 4 alpha=" << alpha
       << "\n";
    T scal1 = eVal * scal1_basic;
    T scal2 = eVal * scal2_basic;
    os << "LORPERF: LORENTZ_FindPositiveVectors: scal1=" << scal1
       << " scal2=" << scal2 << "\n";
    if (scal1 != scal2) {
      std::cerr << "LORPERF: We should have scal1 = scal\n";
      throw TerminalException{1};
    }
    MyVector<T> eTrans = eVal * eTrans_basic;
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 5 eTrans="
       << StringVector(eTrans) << "\n";
#endif
    MyVector<T> eSol = eVal * eSol_basic;
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 7 eSol="
       << StringVector(eSol) << "\n";
#endif
    T beta = eVal * alpha_basic;
    T eSquareDist = eVal * eVal * eSquareDist_basic;
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 8 eSquareDist=" << eSquareDist
       << " det(GramMat)=" << DeterminantMat(GramMat) << "\n";
#endif
    auto iife_iele = [&]() -> std::vector<MyVector<Tint>> {
      if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
        return solver.fixed_dist_vectors(eSol, eSquareDist);
      } else {
        return solver.at_most_dist_vectors(eSol, eSquareDist);
      }
    };
    std::vector<MyVector<Tint>> LVect = iife_iele();
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 9 |LVect|="
       << LVect.size() << "\n";
#endif
    for (auto &eSolA : LVect) {
      MyVector<Tint> eSolC = eBasSol + Ubasis.transpose() * eSolA;
#ifdef SANITY_CHECK_LORENTZIAN_FIND_POSITIVE_VECTORS
      MyVector<T> eSolC_T = UniversalVectorConversion<T, Tint>(eSolC);
      MyVector<T> eSolA_T = UniversalVectorConversion<T, Tint>(eSolA);
      T scal = eSolC_T.dot(eVect_LorMat);
      T scal_expe2 = eVal * TheRec.gcd;
      if (scal != scal_expe2) {
        std::cerr << "LORPERF: scal=" << scal << " scal_expe2=" << scal_expe2 << "\n";
        throw TerminalException{1};
      }
      if (scal > MaxScal) {
        std::cerr << "LORPERF: scal=" << scal << " MaxScal=" << MaxScal << "\n";
        throw TerminalException{1};
      }
      MyVector<T> eDiffA = eSolA_T - eSol;
      T normA = EvaluationQuadForm(GramMat, eDiffA);
      T normC = EvaluationQuadForm(LorMat, eSolC);
      if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
        if (normA != eSquareDist) {
          std::cerr << "LORPERF: normA=" << normA << " eSquareDist=" << eSquareDist
                    << " but should be equal\n";
          throw TerminalException{1};
        }
        if (normC != 0) {
          std::cerr << "LORPERF: normC=" << normC << " but it should be isotrop\n";
          throw TerminalException{1};
        }
      } else {
        if (normA > eSquareDist) {
          std::cerr << "LORPERF: normA=" << normA << " eSquareDist=" << eSquareDist
                    << " but we should have normA <= eSquareDist\n";
          throw TerminalException{1};
        }
        if (normC < 0) {
          std::cerr << "LORPERF: normC=" << normC
                    << " but it should be of zero or positive norm\n";
          throw TerminalException{1};
        }
      }
#endif
      TotalListSol.emplace_back(std::move(eSolC));
    }
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
    os << "LORPERF: LORENTZ_FindPositiveVectors: while step 10 |TotalListSol|="
       << TotalListSol.size() << "\n";
#endif
    if (OnlyShortest && TotalListSol.size()) {
#ifdef GAPDEBUGOUT_LORENTZIAN_FIND_POSITIVE_VECTORS
      f_save(TotalListSol);
#endif
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: LORENTZ_FindPositiveVectors 1, eVal=" << eVal << " |TotalListSol|=" << TotalListSol.size() << "\n";
#endif
#ifdef TIMINGS_LORENTZIAN_PERFECT
      os << "|LORPERF: LORENTZ_FindPositiveVectors 1|=" << time << "\n";
#endif
      return TotalListSol;
    }
    eVal += 1;
    if (MaxScal > 0) {
      T scal = eVal * TheRec.gcd;
      if (scal > MaxScal) {
#ifdef DEBUG_LORENTZIAN_FIND_POSITIVE_VECTORS
        os << "LORPERF: LORENTZ_FindPositiveVectors: doing a break scal="
           << scal << " MaxScal=" << MaxScal << "\n";
#endif
        break;
      }
    }
  }
#ifdef GAPDEBUGOUT_LORENTZIAN_FIND_POSITIVE_VECTORS
  f_save(TotalListSol);
#endif
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: LORENTZ_FindPositiveVectors 2, eVal=" << eVal << " |TotalListSol|=" << TotalListSol.size() << "\n";
#endif
#ifdef TIMINGS_LORENTZIAN_PERFECT
  os << "|LORPERF: LORENTZ_FindPositiveVectors 2|=" << time << "\n";
#endif
  return TotalListSol;
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
LORENTZ_FindPositiveVectors(MyMatrix<T> const &LorMat, MyVector<T> const &eVect,
                            T const &MaxScal, int const &TheOption,
                            bool const &OnlyShortest, std::ostream &os) {
  FractionVector<T> eRec = RemoveFractionVectorPlusCoeff(eVect);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: LORENTZ_FindPositiveVectors, beginning, eRec.TheMult="
     << eRec.TheMult << "\n";
#endif
  MyVector<T> const &eVectNew = eRec.TheVect;
  T MaxScalNew = MaxScal * eRec.TheMult;
  return LORENTZ_FindPositiveVectorsKernel<T, Tint>(
      LorMat, eVectNew, MaxScalNew, TheOption, OnlyShortest, os);
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
LORENTZ_SearchInitialVector(MyMatrix<T> const &LorMat,
                            MyVector<T> const &PosVect, int const &TheOption,
                            std::ostream &os) {
  bool OnlyShortest = true;
  T MaxScal = EvaluationQuadForm(LorMat, PosVect);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: LORENTZ_SearchInitialVector: PosVect=" << StringVectorGAP(PosVect) << "\n";
  os << "LORPERF: LORENTZ_SearchInitialVector: MaxScal=" << MaxScal << "\n";
#endif
  // For the "total" it terminates at the first call.
  // For "isotropic", more iterations may be needed.
  while(true) {
    std::vector<MyVector<Tint>> LVect = LORENTZ_FindPositiveVectors<T, Tint>(LorMat, PosVect, MaxScal,
                                                                             TheOption, OnlyShortest, os);
    if (LVect.size() > 0) {
      return LVect;
    }
    MaxScal *= 2;
  }
}

template <typename T>
std::vector<MyVector<T>>
GetRationalIsotropyVectors(MyMatrix<T> const &LorMat22) {
  // The matrix is written as
  // | a b |
  // | b c |
  // The quadratic form is a x^2 + 2b x y + c y^2
  // We have Delta = 4 (b^2 - a c)
  // So { a (x / y)^2 + 2b (x/y) + c } y^2
  // This gets us
  // x1 / y1 = (-2 b + sqrt(Delta) ) / (2a)
  //         = (-b + sqrt(DeltaB)) / a
  // x2 / y2 = (-2 b - sqrt(Delta) ) / (2a)
  //         = (-b - sqrt(DeltaB) ) / a
  T a = LorMat22(0, 0);
  T b = LorMat22(0, 1);
  T c = LorMat22(1, 1);
  T DeltaB = b * b - a * c;
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
  if (DeltaB <= 0) {
    std::cerr << "LORPERF: If the discriminant is negative, then we cannot have "
                 "signature (1,1)";
    throw TerminalException{1};
  }
#endif
  std::optional<T> opt = UniversalSquareRoot(DeltaB);
  if (!opt) {
    return {};
  }
  T TheSqrt = *opt;
  MyVector<T> v1(2);
  MyVector<T> v2(2);
  if (a != 0) {
    T root1 = (-b + TheSqrt) / a;
    T root2 = (-b - TheSqrt) / a;
    v1(0) = root1;
    v1(1) = 1;
    v2(0) = root2;
    v2(1) = 1;
  } else {
    v1(0) = 1;
    v1(1) = 0;
    // Remaining equation is 2 b x + c y = 0
    v2(0) = -c;
    v2(1) = 2 * b;
  }
  return {v1, v2};
}

template <typename T> MyVector<T> GetReducedVector(MyVector<T> const &V) {
  int n = V.size() - 1;
  MyVector<T> Vret(n);
  for (int i = 0; i < n; i++) {
    Vret(i) = V(i + 1);
  }
  return Vret;
}

template <typename T>
MyVector<T> ConcatenateScalarVector(T const &scal, MyVector<T> const &V) {
  int n = V.size();
  MyVector<T> Vret(n + 1);
  Vret(0) = scal;
  for (int i = 0; i < n; i++) {
    Vret(i + 1) = V(i);
  }
  return Vret;
}

// We compute an upper bound TheMult on the possible values of x such that
// eNSP = eNSPbas + x * eNSPdir is admissible as a facet vector
//
// We also check for the isotropy situation that could happen and report it
// on the output. That is if the best x_upp gives eNSP[1] = 0 and eNSP{[2..n+1]}
// is isotropic than in the following call to the LORENTZ_FindPositiveVectors
// the following happens:
// -- MaxScal = 0 (because eNSP[1] = 0)
// -- eVect * LorMat * fVect = 0 (because 0 <= val <= 0)
// -- Thus fVect is in the orthogonal of an isotropic vector and is of norm >=
// 0.
// -- By the geometry this gets us to fVect a multiple of eVect
// That scenario is not acceptable for finding perfect domain.
template <typename T>
std::optional<T> GetUpperBound(MyMatrix<T> const &LorMat,
                               MyVector<T> const &eNSPbas,
                               MyVector<T> const &eNSPdir,
                               [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetUpperBound, step 1\n";
#endif
  MyMatrix<T> LorMatInv = Inverse(LorMat);
  T eCstBas = eNSPbas(0);
  T eCstDir = eNSPdir(0);
  MyVector<T> eBas = GetReducedVector(eNSPbas);
  MyVector<T> eDir = GetReducedVector(eNSPdir);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetUpperBound, step 2\n";
#endif
  // For an acceptable vector w of length n+1 in return we must have w[0] < 0.
  // Since we have w = u + TheMult * v we have a potential upper bound
  // on TheMult, but only if v[0] > 0
  std::vector<T> ListUpperBound;
  std::optional<T> UpperBound_constant;
  if (eCstDir > 0) {
    UpperBound_constant = -eCstBas / eCstDir;
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
    if (*UpperBound_constant <= 0) {
      std::cerr << "LORPERF: The upper bound from constant is not correct\n";
      throw TerminalException{1};
    }
#endif
    ListUpperBound.push_back(*UpperBound_constant);
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetUpperBound, step 3\n";
#endif
  //
  // Get raw upper bound
  //
  T iShift = 1;
#ifdef DEBUG_LORENTZIAN_PERFECT
  MyVector<T> eBas_V = LorMatInv * eBas;
  MyVector<T> eDir_V = LorMatInv * eDir;
  T normBas = EvaluationQuadForm<T, T>(LorMat, eBas_V);
  T normDir = EvaluationQuadForm<T, T>(LorMat, eDir_V);
  T scalBasDir = ScalarProductQuadForm(LorMat, eBas_V, eDir_V);
  os << "LORPERF: GetUpperBound, normBas=" << normBas << " normDir=" << normDir << " scalBasDir=" << scalBasDir << "\n";
#endif
  while (true) {
    MyVector<T> eV = eBas + iShift * eDir;
    MyVector<T> eVect = LorMatInv * eV;
    T eNorm = EvaluationQuadForm<T, T>(LorMat, eVect);
    if (eNorm < 0) {
      ListUpperBound.push_back(iShift);
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: UpperBound iShift=" << iShift << "\n";
#endif
      break;
    }
    iShift *= 2;
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetUpperBound, step 4\n";
#endif
  //
  // More subttle upper bound coming from isotropy computation
  //
  MyVector<T> eVectBas = LorMatInv * eBas;
  MyVector<T> eVectDir = LorMatInv * eDir;
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
  auto f_norm = [&](T const &val) -> T {
    MyVector<T> eVect = eVectBas + val * eVectDir;
    return EvaluationQuadForm<T, T>(LorMat, eVect);
  };
#endif
  std::vector<MyVector<T>> TheBasis{eVectBas, eVectDir};
  MyMatrix<T> TheBasis_mat = MatrixFromVectorFamily(TheBasis);
  MyMatrix<T> LorMat22 = TheBasis_mat * LorMat * TheBasis_mat.transpose();
  std::vector<MyVector<T>> ListIso = GetRationalIsotropyVectors(LorMat22);
  std::optional<T> UpperBound_isotropic_opt;
  T UpperBound_isotropic;
  std::vector<T> ListUpperBound_Iso;
  for (auto &eIso : ListIso) {
    if (eIso(0) != 0) {
      T fact = eIso(1) / eIso(0);
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
      if (f_norm(fact) != 0) {
        std::cerr << "LORPERF: eVect should be isotropic\n";
        throw TerminalException{1};
      }
#endif
      if (fact > 0) {
#ifdef DEBUG_LORENTZIAN_PERFECT
        os << "LORPERF: UpperBound fact=" << fact << "\n";
#endif
        ListUpperBound_Iso.push_back(fact);
      }
    }
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetUpperBound, step 5\n";
#endif
  if (ListUpperBound_Iso.size() > 0) {
    UpperBound_isotropic_opt = VectorMin(ListUpperBound_Iso);
    UpperBound_isotropic = *UpperBound_isotropic_opt;
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: UpperBound UpperBound_isotropic=" << UpperBound_isotropic << "\n";
#endif
    ListUpperBound.push_back(UpperBound_isotropic);
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetUpperBound, step 6\n";
#endif
  if (UpperBound_constant.has_value() && UpperBound_isotropic_opt.has_value()) {
    if (*UpperBound_constant == *UpperBound_isotropic_opt) {
      return {};
    }
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetUpperBound, step 7\n";
#endif
  // Need to see if better upper bounds are possible, but this is a secondary
  // question
  T BestUpper = VectorMin(ListUpperBound);
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
  if (f_norm(BestUpper) > 0) {
    std::cerr << "LORPERF: We should have eNorm <= 0\n";
    throw TerminalException{1};
  }
#endif
  return BestUpper;
}

template <typename T, typename Tint> struct ResultFlipping {
  std::vector<MyVector<Tint>> ListTotal;
  MyVector<T> eNSPtest;
  MyVector<T> eVectTest;
  T MaxScal;
};

// Given a Critical set a vector eNSPbas and a direction eNSPdir
// we are looking for a lambda > 0 such that eNSP = eNSPbas + lambda * eNSPdir
// such that the list of vectors satisfying eNSP * v = 0 defines a bigger
// set than CritSet.
//
// eNSPbas must be of positive norm. eNSPdir must be of negative norm.
//
template <typename T, typename Tint>
ResultFlipping<T, Tint> LORENTZ_Kernel_Flipping(
    MyMatrix<T> const &LorMat, std::vector<MyVector<Tint>> const &CritSet,
    MyVector<T> const &eNSPbas, MyVector<T> const &eNSPdir,
    int const &TheOption, std::ostream &os) {
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: Beginning of LORENTZ_Kernel_Flipping\n";
#endif
  MyMatrix<T> LorMatInv = Inverse(LorMat);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: We have LorMatInv. LorMat=\n";
  WriteMatrix(os, LorMat);
#endif
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
  std::vector<MyVector<T>> M1{eNSPbas, eNSPdir};
  MyMatrix<T> M2 = MatrixFromVectorFamily(M1);
  if (RankMat(M2) != 2) {
    std::cerr
        << "The vector eNSPbas and eNSPdir should be linearly independent\n";
    throw TerminalException{1};
  }
  for (auto &x : CritSet) {
    MyVector<Tint> xext = ConcatenateScalarVector(Tint(1), x);
    MyVector<T> xext_T = UniversalVectorConversion<T, Tint>(xext);
    T scal_bas = xext_T.dot(eNSPbas);
    T scal_dir = xext_T.dot(eNSPdir);
    if (scal_bas != 0) {
      std::cerr << "LORPERF: Kernel_Flipping: scal_bas=" << scal_bas
                << " scal_dir=" << scal_dir << "\n";
      std::cerr << "LORPERF: eNSPbas should have scalar product 0 with all entries in "
                   "CritSet\n";
      throw TerminalException{1};
    }
    if (scal_dir != 0) {
      std::cerr << "LORPERF: eNSPdir should have scalar product 0 with all entries in "
                   "CritSet\n";
      throw TerminalException{1};
    }
    T norm = EvaluationQuadForm<T, Tint>(LorMat, x);
    if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
      if (norm != 0) {
        std::cerr << "LORPERF: CritSet contains some non-isotrop vectors\n";
        throw TerminalException{1};
      }
    } else {
      if (norm < 0) {
        std::cerr << "LORPERF: CritSet contains some a vector of negative norm\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  bool OnlyShortest = true;
  T TheLowerBound = 0;
  std::optional<T> TheUpperBound_opt = GetUpperBound(LorMat, eNSPbas, eNSPdir, os);
  if (!TheUpperBound_opt) {
    std::cerr << "LORPERF: We have TheUpperBound that is none and that is a problem\n";
    throw TerminalException{1};
  }
  T TheUpperBound = *TheUpperBound_opt;
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: TheUpperBound=" << TheUpperBound << "\n";
#endif
#ifdef DEBUG_LORENTZIAN_PERFECT
  int n_iter = 0;
#endif
  MyVector<T> eVectBas = LorMatInv * GetReducedVector(eNSPbas);
  MyVector<T> eVectDir = LorMatInv * GetReducedVector(eNSPdir);
  T eNormBas = EvaluationQuadForm<T, T>(LorMat, eVectBas);
  T eNormDir = EvaluationQuadForm<T, T>(LorMat, eVectDir);
  std::vector<MyVector<Tint>> ListTotal;
  while (true) {
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: Kernel_Flipping: TheLowerBound=" << TheLowerBound
       << " TheUpperBound=" << TheUpperBound << " n_iter=" << n_iter << "\n";
#endif
    T TheMidVal = get_mid_val(TheLowerBound, TheUpperBound);
    MyVector<T> eNSPtest = eNSPbas + TheMidVal * eNSPdir;
    MyVector<T> eVectTest = LorMatInv * GetReducedVector(eNSPtest);
    T eNormTest = EvaluationQuadForm<T, T>(LorMat, eVectTest);
    MyVector<T> CritSet0_T = UniversalVectorConversion<T, Tint>(CritSet[0]);
    T MaxScal = ScalarProductQuadForm(LorMat, CritSet0_T, eVectTest);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: eNormTest=" << eNormTest << " MaxScal=" << MaxScal << "\n";
#endif
    if (eNormTest <= 0 || MaxScal <= 0) {
      TheUpperBound = TheMidVal;
    } else {
      ListTotal = LORENTZ_FindPositiveVectors<T, Tint>(
          LorMat, eVectTest, MaxScal, TheOption, OnlyShortest, os);
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
      if (IsSubset(CritSet, ListTotal) && CritSet.size() > ListTotal.size()) {
        std::cerr << "LORPERF: |ListTotal|=" << ListTotal.size()
                  << " |CritSet|=" << CritSet.size() << "\n";
        std::cerr << "LORPERF: LorMat=" << StringMatrixGAP(LorMat) << "\n";
        std::cerr << "LORPERF: eVectTest=" << StringVectorGAP(eVectTest)
                  << "\n";
        std::cerr << "LORPERF: MaxScal=" << MaxScal
                  << " TheOption=" << TheOption
                  << " OnlyShortest=" << OnlyShortest << "\n";
        std::cerr << "LORPERF: |CritSet|=" << CritSet.size() << "\n";
        if (CritSet.size() > 0) {
          WriteMatrix(std::cerr, MatrixFromVectorFamily(CritSet));
        }
        std::cerr << "LORPERF: |ListTotal|=" << ListTotal.size() << "\n";
        if (ListTotal.size() > 0) {
          WriteMatrix(std::cerr, MatrixFromVectorFamily(ListTotal));
        }
        std::cerr << "LORPERF: Bug: if included, it should be equal\n";
        throw TerminalException{1};
      }
#endif
      if (IsEqualSet(ListTotal, CritSet)) {
        TheLowerBound = TheMidVal;
      } else {
        if (IsSubset(ListTotal, CritSet)) {
#ifdef DEBUG_LORENTZIAN_PERFECT
          os << "LORPERF: Kernel_Flipping: EXIT 1 |ListTotal|="
             << ListTotal.size() << " MaxScal=" << MaxScal << "\n";
          os << "LORPERF: Kernel_Flipping: NSPtest="
             << StringVectorGAP(eNSPtest) << "\n";
#endif
          return {ListTotal, eNSPtest, eVectTest, MaxScal};
        } else {
          break;
        }
      }
    }
#ifdef DEBUG_LORENTZIAN_PERFECT
    n_iter += 1;
#endif
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: Kernel_Flipping: Going to the second scheme\n";
#endif
  while (true) {
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
    if (ListTotal.size() == 0) {
      std::cerr << "LORPERF: ListTotal is empty, so that is a bug\n";
      throw TerminalException{1};
    }
#endif
    MyVector<Tint> eVect = ListTotal[0];
    MyVector<Tint> eVert = ConcatenateScalarVector(Tint(1), eVect);
    MyVector<T> eVert_T = UniversalVectorConversion<T, Tint>(eVert);
    std::vector<MyVector<Tint>> ListIsoTest = CritSet;
    ListIsoTest.push_back(eVect);
    T aShift = -(eNSPbas.dot(eVert_T)) / (eNSPdir.dot(eVert_T));
    MyVector<T> eNSPtest = eNSPbas + aShift * eNSPdir;
    MyVector<T> eVectTest = LorMatInv * GetReducedVector(eNSPtest);
    MyVector<T> CritSet0_T = UniversalVectorConversion<T, Tint>(CritSet[0]);
    T MaxScal = ScalarProductQuadForm(LorMat, CritSet0_T, eVectTest);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: Kernel_Flipping: MaxScal=" << MaxScal
       << " n_iter=" << n_iter << "\n";
#endif
    ListTotal = LORENTZ_FindPositiveVectors<T, Tint>(
        LorMat, eVectTest, MaxScal, TheOption, OnlyShortest, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: Kernel_Flipping: |ListTotal|=" << ListTotal.size() << "\n";
    n_iter += 1;
#endif
    if (IsSubset(ListTotal, CritSet)) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: EXIT 2 |ListTotal|=" << ListTotal.size()
         << " MaxScal=" << MaxScal << "\n";
      os << "LORPERF: NSPtest=" << StringVectorGAP(eNSPtest) << "\n";
#endif
      return {ListTotal, eNSPtest, eVectTest, MaxScal};
    }
  }
}

template <typename T, typename Tint>
MyVector<T> LORENTZ_GetOneOutsideRay(MyMatrix<T> const &LorMat,
                                     MyMatrix<T> const &SpannBasis,
                                     std::vector<MyVector<Tint>> const &TheSet,
                                     MyVector<T> const &eNSPbas,
                                     std::ostream &os) {
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOneOutsideRay: begin\n";
#endif
  MyVector<T> TheSet0 = UniversalVectorConversion<T, Tint>(TheSet[0]);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOneOutsideRay: We have TheSet0\n";
#endif
  MyMatrix<T> TheMat = SpannBasis * LorMat * SpannBasis.transpose();
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
  DiagSymMat<T> dsm = DiagonalizeSymmetricMatrix(TheMat, os);
  if (dsm.nbMinus == 0) {
    std::cerr << "LORPERF: We should have a negative in the entry\n";
    throw TerminalException{1};
  }
#endif
  int n_dim = TheMat.rows();
  MyMatrix<T> ePerturb = IdentityMat<T>(n_dim);
#ifdef DEBUG_LORENTZIAN_PERFECT
  size_t iter = 0;
  os << "LORPERF: GetOneOutsideRay: Before while loop\n";
#endif
  while (true) {
    MyMatrix<T> TheMatPerturb = -ePerturb * TheMat * ePerturb.transpose();
#ifdef DEBUG_LORENTZIAN_PERFECT
    iter += 1;
    os << "LORPERF: GetOneOutsideRay: iter=" << iter << "\n";
#endif
    MyMatrix<T> uVect =
        GetIntegralPositiveVector_allmeth<T, T>(TheMatPerturb, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOneOutsideRay: We have uVect\n";
#endif
    MyMatrix<T> SpannMatrixPert = ePerturb * SpannBasis;
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOneOutsideRay: We have SpannMatrixPert\n";
#endif
    MyVector<T> RetVect = SpannMatrixPert.transpose() * uVect;
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
    T TheNorm = EvaluationQuadForm<T, T>(LorMat, RetVect);
    if (TheNorm == 0) {
      std::cerr << "LORPERF: GetOneOutsideRay: TheNorm=" << TheNorm << "\n";
      std::cerr << "LORPERF: The vector should be outside of the cone and so have "
                   "negative norm\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> tVect = LorMat * RetVect;
    T eScal = tVect.dot(TheSet0);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOneOutsideRay: tVect=" << tVect << " eScal=" << eScal << "\n";
#endif
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
    for (auto &eVect : TheSet) {
      MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect);
      T fScal = tVect.dot(eVect_T);
      if (eScal != fScal) {
        std::cerr << "LORPERF: eScal=" << eScal << " fScal=" << fScal << "\n";
        std::cerr << "LORPERF: The scalar products are incorrect\n";
        throw TerminalException{1};
      }
    }
#endif
    MyVector<T> eNSPdir = ConcatenateScalarVector(T(-eScal), tVect);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOneOutsideRay: We have eNSPdir\n";
#endif
    std::optional<T> TheUpperBound_opt =
      GetUpperBound(LorMat, eNSPbas, eNSPdir, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOneOutsideRay: We have TheUpperBound_opt\n";
#endif
    if (TheUpperBound_opt.has_value()) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: GetOneOutsideRay: iter=" << iter << " eNSPdir=" << StringVector(eNSPdir) << "\n";
#endif
      return eNSPdir;
    }
    ePerturb = ePerturb * GetRandomMatrixPerturbation<T>(n_dim);
  }
}

template <typename T, typename Tint>
MyMatrix<T> GetFullExpanded(std::vector<MyVector<Tint>> const &CritSet) {
  int dim = CritSet[0].size();
  int n_vect = CritSet.size();
  MyMatrix<T> ListVectExt(n_vect, dim + 1);
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    ListVectExt(i_vect, 0) = 1;
    for (int i = 0; i < dim; i++) {
      ListVectExt(i_vect, i + 1) =
          UniversalScalarConversion<T, Tint>(CritSet[i_vect](i));
    }
  }
  return ListVectExt;
}

template <typename T, typename Tint>
void LORENTZ_CheckCorrectnessVectorFamily(
    MyMatrix<T> const &LorMat, std::vector<MyVector<Tint>> const &CritSet) {
  int n = LorMat.cols();
  MyMatrix<T> ListVectExt = GetFullExpanded<T, Tint>(CritSet);
  int rnk = RankMat(ListVectExt);
  if (rnk != n) {
    std::cerr << "LORPERF: We have rnk=" << rnk << " n=" << n
              << " they should be equal\n";
    throw TerminalException{1};
  }
}

template <typename T, typename Tint> struct LorentzianPerfectEntry {
  std::vector<MyVector<Tint>> ListTotal;
  MyVector<T> eNSPtest;
  MyVector<T> eVectTest;
};

template <typename T, typename Tint>
MyMatrix<Tint> LORENTZ_GetEXT(LorentzianPerfectEntry<T, Tint> const &eRec) {
  return MatrixFromVectorFamily(eRec.ListTotal);
}

template <typename T, typename Tint>
LorentzianPerfectEntry<T, Tint> LORENTZ_GetOnePerfect(MyMatrix<T> const &LorMat,
                                                      int const &TheOption,
                                                      std::ostream &os) {
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOnePerfect: beginning\n";
  os << "LORPERF: LorMat=\n";
  WriteMatrix(os, LorMat);
#endif
  int n = LorMat.rows();
  MyMatrix<T> LorMatInv = Inverse(LorMat);
  MyVector<Tint> CentralVect =
      INDEFINITE_GetShortPositiveVector<T, Tint>(LorMat, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOnePerfect: We have CentralVect\n";
#endif
  MyVector<T> CentralVect_T = UniversalVectorConversion<T, Tint>(CentralVect);
  std::vector<MyVector<Tint>> CritSet = LORENTZ_SearchInitialVector<T, Tint>(
      LorMat, CentralVect_T, TheOption, os);
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
  if (CritSet.size() == 0) {
    std::cerr << "LORPERF: CritSet is empty which ruins everything\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOnePerfect: We have CritSet |CritSet|=" << CritSet.size() << "\n";
#endif
  MyVector<T> CritSet0_T = UniversalVectorConversion<T, Tint>(CritSet[0]);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOnePerfect: We have CritSet0_T\n";
#endif
  MyVector<T> LorMat_Central = LorMat * CentralVect_T;
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOnePerfect: We have LorMat_Central\n";
#endif
  T eScal = LorMat_Central.dot(CritSet0_T);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOnePerfect: We have eScal\n";
#endif
  MyVector<T> eNSPbas = ConcatenateScalarVector(T(-eScal), LorMat_Central);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: GetOnePerfect: Before loop\n";
#endif
  while (true) {
    int rnk = 0;
    if (CritSet.size() > 0) {
      rnk = RankMat(MatrixFromVectorFamily(CritSet));
    }
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOnePerfect n=" << n << " rnk=" << rnk
       << " |CritSet|=" << CritSet.size() << "\n";
#endif
    if (rnk == n) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      LORENTZ_CheckCorrectnessVectorFamily(LorMat, CritSet);
#endif
      return {CritSet, eNSPbas, CentralVect_T};
    }
    MyMatrix<T> EXT = GetFullExpanded<T, Tint>(CritSet);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOnePerfect: We have |EXT|=" << EXT.rows() << " / "
       << EXT.cols() << "\n";
#endif
    MyMatrix<T> NSP = NullspaceIntTrMat(EXT);
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
    if (NSP.rows() == 0) {
      std::cerr << "LORPERF: NSP should be non-empty\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<T> ListDir = DropColumn(NSP, 0) * LorMatInv;
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOnePerfect: |ListDir|=" << ListDir.rows() << " / " << ListDir.cols() << "\n";
    os << "LORPERF: GetOnePerfect: ListDir=\n";
    WriteMatrix(std::cerr, ListDir);
#endif
    MyMatrix<T> eNSPdir = LORENTZ_GetOneOutsideRay<T, Tint>(
        LorMat, ListDir, CritSet, eNSPbas, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOnePerfect: Before LORENTZ_Kernel_Flipping\n";
#endif
    ResultFlipping<T, Tint> eRecB = LORENTZ_Kernel_Flipping<T, Tint>(
        LorMat, CritSet, eNSPbas, eNSPdir, TheOption, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: GetOnePerfect: After LORENTZ_Kernel_Flipping\n";
#endif
    CritSet = eRecB.ListTotal;
    eNSPbas = eRecB.eNSPtest;
  }
}

template <typename T, typename Tint>
LorentzianPerfectEntry<T, Tint>
LORENTZ_DoFlipping(MyMatrix<T> const &LorMat,
                   std::vector<MyVector<Tint>> const &ListIso, Face eInc,
                   int const &TheOption, std::ostream &os) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
  MicrosecondTime time;
#endif
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: DoFlipping: beginning\n";
  os << "LORPERF: DoFlipping: |ListIso|=" << ListIso.size() << "\n";
  os << "LORPERF: LorMat=\n";
  WriteMatrix(os, LorMat);
  os << "LORPERF: ListIso=\n";
  WriteMatrix(os, MatrixFromVectorFamily(ListIso));
#endif
  size_t n_vect = eInc.size();
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: DoFlipping: |eInc|=" << eInc.size() << " / " << eInc.count()
     << "\n";
#endif
  MyMatrix<T> EXT = GetFullExpanded<T, Tint>(ListIso);
  auto get_eVert = [&]() -> size_t {
    for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
      if (eInc[i_vect] == 0) {
        return i_vect;
      }
    }
    std::cerr << "LORPERF: Failed to find a matching entry in get_eVert\n";
    throw TerminalException{1};
  };
  size_t eVert = get_eVert();
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: DoFlipping: eVert=" << eVert << "\n";
#endif
  std::vector<MyVector<T>> ListIsoSel;
  for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
    if (eInc[i_vect] == 1) {
      MyVector<T> eVect = UniversalVectorConversion<T, Tint>(ListIso[i_vect]);
      ListIsoSel.push_back(eVect);
    }
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: DoFlipping: |ListIsoSel|=" << ListIsoSel.size() << "\n";
#endif
  MyMatrix<T> MatrIsoSel = MatrixFromVectorFamily(ListIsoSel);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: DoFlipping: |MatrIsoSel|=" << MatrIsoSel.rows() << " / "
     << MatrIsoSel.cols() << "\n";
#endif
  MyMatrix<T> NSP = NullspaceTrMat(MatrIsoSel);
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
  if (NSP.rows() != 1) {
    std::cerr << "LORPERF: NSP should have size 1\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> TheDir = GetMatrixRow(NSP, 0);
  MyVector<T> eIso = UniversalVectorConversion<T, Tint>(ListIso[eVert]);
  auto get_eNSPdir = [&]() -> MyVector<T> {
    T eScal = TheDir.dot(eIso);
    if (eScal < 0) {
      MyVector<T> TheDirNeg = -TheDir;
      return ConcatenateScalarVector(T(0), TheDirNeg);
    } else {
      return ConcatenateScalarVector(T(0), TheDir);
    }
  };
  MyVector<T> eNSPdir = get_eNSPdir();
  MyMatrix<T> NSPb = NullspaceTrMat(EXT);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: DoFlipping: |NSPb|=" << NSPb.rows() << " / " << NSPb.cols()
     << "\n";
#endif
  MyVector<T> NSPb_0 = GetMatrixRow(NSPb, 0);
  MyVector<T> eVectB = GetReducedVector(NSPb_0);
  auto iife_get_eNSPbas = [&]() -> MyVector<T> {
    if (eVectB.dot(eIso) > 0) {
      return NSPb_0;
    } else {
      return -NSPb_0;
    }
  };
  MyVector<T> eNSPbas = iife_get_eNSPbas();
  std::vector<MyVector<Tint>> CritSet;
  for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
    if (eInc[i_vect] == 1) {
      CritSet.push_back(ListIso[i_vect]);
    }
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: DoFlipping: Before LORENTZ_Kernel_Flipping |CritSet|="
     << CritSet.size() << "\n";
#endif
  std::vector<MyVector<Tint>> TheFlip =
      LORENTZ_Kernel_Flipping<T, Tint>(LorMat, CritSet, eNSPbas, eNSPdir,
                                       TheOption, os)
          .ListTotal;
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: DoFlipping: After LORENTZ_Kernel_Flipping |TheFlip|="
     << TheFlip.size() << "\n";
#endif
#ifdef DEBUG_LORENTZIAN_PERFECT
  LORENTZ_CheckCorrectnessVectorFamily(LorMat, TheFlip);
#endif
  // No return of additional info so far.
#ifdef TIMINGS_LORENTZIAN_PERFECT
  os << "|LORPERF: LORENTZ_DoFlipping|=" << time << "\n";
#endif
  return {TheFlip, {}, {}};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LORENTZ_TestEquivalence(MyMatrix<T> const &LorMat1, MyMatrix<T> const &LorMat2,
                        MyMatrix<Tint> const &eFamEXT1,
                        MyMatrix<Tint> const &eFamEXT2, std::ostream &os) {
  return LinPolytopeIntegral_Isomorphism_GramMat<T, Tint, Tgroup>(
      eFamEXT1, LorMat1, eFamEXT2, LorMat2, os);
}

template <typename Tint, typename Tgroup> struct ResultStabilizer {
  std::vector<MyMatrix<Tint>> ListGen;
  Tgroup GRPperm;
};

template <typename T, typename Tint, typename Tgroup>
ResultStabilizer<Tint, Tgroup>
LORENTZ_ComputeStabilizer(MyMatrix<T> const &LorMat,
                          MyMatrix<Tint> const &eFamEXT, std::ostream &os) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
  MicrosecondTime time;
#endif
  MyMatrix<T> eFamEXT_T = UniversalMatrixConversion<T, Tint>(eFamEXT);
  Tgroup GRPisom =
      LinPolytope_Automorphism_GramMat<T, Tgroup>(eFamEXT_T, LorMat, os);
  Tgroup GRPperm =
      LinPolytopeIntegral_Stabilizer_Method8(eFamEXT_T, GRPisom, os);
  std::vector<MyMatrix<Tint>> ListGen;
  for (auto &eGen : GRPperm.SmallGeneratingSet()) {
    MyMatrix<Tint> eGenMatr =
        RepresentVertexPermutation(eFamEXT, eFamEXT, eGen);
    ListGen.push_back(eGenMatr);
  }
#ifdef TIMINGS_LORENTZIAN_PERFECT
  os << "|LORPERF: LORENTZ_ComputeStabilizer|=" << time << "\n";
#endif
  return {ListGen, GRPperm};
}

template <typename T, typename Tint>
size_t ComputeInvariantPerfectForm(size_t seed, MyMatrix<T> const &LorMat,
                                   MyMatrix<Tint> const &EXT,
                                   [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
  MicrosecondTime time;
#endif
  int n = LorMat.rows();
  int nbRow = EXT.rows();
  MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(EXT);
  std::map<T, size_t> ListDiagNorm;
  std::map<T, size_t> ListOffDiagNorm;
  MyVector<T> V(n);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int i = 0; i < n; i++) {
      T sum(0);
      for (int j = 0; j < n; j++) {
        sum += LorMat(i, j) * EXT_T(iRow, j);
      }
      V(i) = sum;
    }
    T scal = 0;
    for (int i = 0; i < n; i++) {
      scal += V(i) * EXT_T(iRow, i);
    }
    ListDiagNorm[scal] += 1;
    for (int jRow = iRow + 1; jRow < nbRow; jRow++) {
      T scal(0);
      for (int i = 0; i < n; i++) {
        scal += V(i) * EXT_T(jRow, i);
      }
      ListOffDiagNorm[scal] += 1;
    }
  }
  size_t hash = ComputeHashTwoMap(seed, ListDiagNorm, ListOffDiagNorm);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: ComputeInvariantPerfectForm |EXT|=" << EXT.rows() << " |ListDiagNorm|=" << ListDiagNorm.size() << " |ListOffDiagNorm|=" << ListOffDiagNorm.size() << " hash=" << hash << "\n";
#endif
#ifdef TIMINGS_LORENTZIAN_PERFECT
  os << "|LORPERF: ComputeInvariantPerfectForm|=" << time << "\n";
#endif
  return hash;
}

template <typename T, typename Tint, typename Tgroup>
struct DataPerfectLorentzian {
  int n;
  MyMatrix<T> LorMat;
  int TheOption;
  RecordDualDescOperation<T, Tgroup> rddo;
};

template <typename Tint, typename Tgroup> struct PerfLorentzian_Obj {
  MyMatrix<Tint> EXT;
  Tgroup GRP;
};

template <typename Tint, typename Tgroup>
void WriteEntryGAP(std::ostream &os_out,
                   PerfLorentzian_Obj<Tint, Tgroup> const &ent) {
  os_out << "rec(EXT:=";
  WriteMatrixGAP(os_out, ent.EXT);
  os_out << ", GRP:=" << ent.GRP.GapString() << ")";
}

template <typename T, typename Tint, typename Tgroup>
void WriteDetailedEntryGAP(std::ostream &os_out,
                           [[maybe_unused]] DataPerfectLorentzian<T, Tint, Tgroup> const& data,
                           PerfLorentzian_Obj<Tint, Tgroup> const &ent, [[maybe_unused]] std::ostream& os) {
  WriteEntryGAP(os_out, ent);
}

template <typename Tint, typename Tgroup>
void WriteEntryPYTHON(std::ostream &os_out,
                      PerfLorentzian_Obj<Tint, Tgroup> const &ent) {
  os_out << "{\"EXT\":";
  WriteMatrixPYTHON(os_out, ent.EXT);
  os_out << ", \"GRP\":" << ent.GRP.PythonString() << "}";
}

namespace boost::serialization {
template <class Archive, typename Tint, typename Tgroup>
inline void serialize(Archive &ar, PerfLorentzian_Obj<Tint, Tgroup> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("EXT", eRec.EXT);
  ar &make_nvp("GRP", eRec.GRP);
}
} // namespace boost::serialization

template <typename Tint> struct PerfLorentzian_AdjI {
  Face eInc;
  MyMatrix<Tint> EXT;
};

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, PerfLorentzian_AdjI<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("EXT", eRec.EXT);
}
} // namespace boost::serialization

template <typename Tint> struct PerfLorentzian_AdjO {
  Face eInc;
  MyMatrix<Tint> eBigMat;
};

template <typename Tint>
void WriteEntryGAP(std::ostream &os_out, PerfLorentzian_AdjO<Tint> const &ent) {
  os_out << "rec(eInc:=";
  WriteFaceGAP(os_out, ent.eInc);
  os_out << ", eBigMat:=" << StringMatrixGAP(ent.eBigMat) << ")";
}

template <typename Tint>
void WriteEntryPYTHON(std::ostream &os_out, PerfLorentzian_AdjO<Tint> const &ent) {
  os_out << "{\"eInc\":";
  WriteFacePYTHON(os_out, ent.eInc);
  os_out << ", \"eBigMat\":" << StringMatrixPYTHON(ent.eBigMat) << "}";
}

namespace boost::serialization {
template <class Archive, typename Tint>
inline void serialize(Archive &ar, PerfLorentzian_AdjO<Tint> &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("eInc", eRec.eInc);
  ar &make_nvp("eBigMat", eRec.eBigMat);
}
} // namespace boost::serialization

template <typename T, typename Tint, typename Tgroup>
struct DataPerfectLorentzianFunc {
  DataPerfectLorentzian<T, Tint, Tgroup> data;
  using Tobj = PerfLorentzian_Obj<Tint, Tgroup>;
  using TadjI = PerfLorentzian_AdjI<Tint>;
  using TadjO = PerfLorentzian_AdjO<Tint>;
  std::ostream &get_os() { return data.rddo.os; }
  Tobj f_init() {
#ifdef TIMINGS_LORENTZIAN_PERFECT
    MicrosecondTime time;
#endif
    std::ostream &os = get_os();
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_init, before LORENTZ_GetOnePerfect\n";
#endif
    LorentzianPerfectEntry<T, Tint> eRec =
        LORENTZ_GetOnePerfect<T, Tint>(data.LorMat, data.TheOption, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_init, after LORENTZ_GetOnePerfect\n";
#endif
    MyMatrix<Tint> EXT = LORENTZ_GetEXT(eRec);
    Tobj x{std::move(EXT), {}};
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_init, returning x\n";
#endif
#ifdef TIMINGS_LORENTZIAN_PERFECT
    os << "|LORPERF: LorPerfect(f_init)|=" << time << "\n";
#endif
    return x;
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return ComputeInvariantPerfectForm(seed, data.LorMat, x.EXT, data.rddo.os);
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
    MicrosecondTime time;
#endif
    std::ostream &os = get_os();
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_repr: Before LORENTZ_TestEquivalence\n";
#endif
    std::optional<MyMatrix<Tint>> opt =
        LORENTZ_TestEquivalence<T, Tint, Tgroup>(data.LorMat, data.LorMat,
                                                 x.EXT, y.EXT, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_repr: After LORENTZ_TestEquivalence\n";
#endif
    if (!opt) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: f_repr: returning not isomorphic\n";
#endif
#ifdef TIMINGS_LORENTZIAN_PERFECT
      os << "|LORPERF: LorPerfect(f_repr,None)|=" << time << "\n";
#endif
      return {};
    }
    MyMatrix<Tint> const &eBigMat = *opt;
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_repr: |eBigMat|=" << eBigMat.rows() << " / "
       << eBigMat.cols() << "\n";
#endif
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
    MyMatrix<T> eBigMat_T = UniversalMatrixConversion<T, Tint>(eBigMat);
    MyMatrix<T> prod = eBigMat_T * data.LorMat * eBigMat_T.transpose();
    if (prod != data.LorMat) {
      std::cerr << "LORPERF: eBigMat should preserve the quadratic form\n";
      throw TerminalException{1};
    }
#endif
    TadjO ret{y.eInc, eBigMat};
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_repr: Before returning ret\n";
#endif
#ifdef TIMINGS_LORENTZIAN_PERFECT
    os << "|LORPERF: LorPerfect(f_repr,Some)|=" << time << "\n";
#endif
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    MyMatrix<Tint> EXT = x.EXT;
    Tobj x_ret{EXT, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.n);
    TadjO ret{x.eInc, eBigMat};
    return {x_ret, ret};
  }
  std::vector<TadjI> f_adj(Tobj &x) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
    MicrosecondTime time;
#endif
    std::ostream &os = get_os();
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_adj: beginning\n";
#endif
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(x.EXT);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_adj: |EXT_T|=" << EXT_T.rows() << " / " << EXT_T.cols()
       << "\n";
#endif
    ResultStabilizer<Tint, Tgroup> res_stab =
        LORENTZ_ComputeStabilizer<T, Tint, Tgroup>(data.LorMat, x.EXT, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_adj: |GRPperm|=" << res_stab.GRPperm.size() << "\n";
#endif
    vectface TheOutput =
        DualDescriptionRecordFullDim(EXT_T, res_stab.GRPperm, data.rddo);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: f_adj: |TheOutput|=" << TheOutput.size() << "\n";
#endif
    x.GRP = res_stab.GRPperm;
    std::vector<MyVector<Tint>> ListIso;
    for (int i_row = 0; i_row < x.EXT.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(x.EXT, i_row);
      ListIso.emplace_back(std::move(V));
    }
    std::vector<TadjI> ListAdj;
    for (auto &eInc : TheOutput) {
      LorentzianPerfectEntry<T, Tint> eRecFlip = LORENTZ_DoFlipping<T, Tint>(
          data.LorMat, ListIso, eInc, data.TheOption, os);
      MyMatrix<Tint> EXTflip = LORENTZ_GetEXT(eRecFlip);
      TadjI eAdj{eInc, EXTflip};
      ListAdj.push_back(eAdj);
    }
#ifdef TIMINGS_LORENTZIAN_PERFECT
    os << "|LORPERF: LorPerfect(f_adj)|=" << time << "\n";
#endif
    return ListAdj;
  }
  Tobj f_adji_obj(TadjI const &x) {
    MyMatrix<Tint> EXT = x.EXT;
    Tobj x_ret{EXT, {}};
    return x_ret;
  }
};

FullNamelist NAMELIST_GetStandard_COMPUTE_PERFECT_LORENTZIAN() {
  std::map<std::string, SingleBlock> ListBlock;
  // SYSTEM
  ListBlock["SYSTEM"] = SINGLEBLOCK_Get_System();
  // DATA
  std::map<std::string, std::string> ListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["LorMatFile"] = "unset.gram";
  ListStringValues1["Option"] = "unset";
  ListStringValues1["FileDualDescription"] = "unset";
  SingleBlock BlockDATA;
  BlockDATA.setListStringValues(ListStringValues1);
  ListBlock["DATA"] = BlockDATA;
  // Merging all data
  return FullNamelist(ListBlock);
}

template <typename Tint, typename Telt, typename Tout>
std::vector<MyMatrix<Tint>>
LORENTZ_ExtractGeneratorsFromObjList(std::vector<Tout> const &l_obj) {
  std::unordered_set<MyMatrix<Tint>> s_gen;
  auto f_insert_gen = [&](MyMatrix<Tint> const &M) -> void {
    if (!IsIdentity(M)) {
      s_gen.insert(M);
    }
  };
  for (auto &ent : l_obj) {
    // Generators coming from the equivalences
    for (auto &eAdj : ent.ListAdj) {
      f_insert_gen(eAdj.x.eBigMat);
    }
    // Generators comoing from the stabilizers
    for (auto &ePermGen : ent.x.GRP.SmallGeneratingSet()) {
      MyMatrix<Tint> eMatrGen = RepresentVertexPermutation<Tint, Telt>(
          ent.x.EXT, ent.x.EXT, ePermGen);
      f_insert_gen(eMatrGen);
    }
  }
  std::vector<MyMatrix<Tint>> l_gen;
  for (auto &eGen : s_gen) {
    l_gen.push_back(eGen);
  }
  return l_gen;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
LORENTZ_GetGeneratorsAutom_Reduced(MyMatrix<T> const &LorMat, std::ostream &os) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
  MicrosecondTime time;
#endif
  int n = LorMat.rows();
  int dimEXT = n + 1;
  using TintGroup = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);

  int TheOption = LORENTZIAN_PERFECT_OPTION_TOTAL;

  DataPerfectLorentzian<T, Tint, Tgroup> data{n, LorMat, TheOption,
                                              std::move(rddo)};
  using Tdata = DataPerfectLorentzianFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  //
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  int max_runtime_second = 0;
  std::vector<Tout> l_obj =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);

  std::vector<MyMatrix<Tint>> l_gen = 
      LORENTZ_ExtractGeneratorsFromObjList<Tint, Telt>(l_obj);
#ifdef TIMINGS_LORENTZIAN_PERFECT
  os << "|LORPERF: LORENTZ_GetGeneratorsAutom|=" << time << "\n";
#endif
  return l_gen;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
LORENTZ_GetGeneratorsAutom_Kernel(MyMatrix<T> const &LorMat, std::ostream &os) {
  ResultReduction<T, Tint> res =
    IndefiniteReduction<T,Tint>(LorMat, os);
  // We have res.B * LorMat * res.B^T = res.Mred
  std::vector<MyMatrix<Tint>> ListGen =
    LORENTZ_GetGeneratorsAutom_Reduced<T,Tint,Tgroup>(res.Mred, os);
  MyMatrix<Tint> Binv = Inverse(res.B);
  //
  std::vector<MyMatrix<Tint>> ListGenRet;
  for (auto & eGen : ListGen) {
    MyMatrix<Tint> eGenRet = Binv * eGen * res.B;
    ListGenRet.push_back(eGenRet);
  }
  return ListGenRet;
}




template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
LORENTZ_GetGeneratorsAutom(MyMatrix<T> const &LorMat, std::ostream &os) {
  int dim = LorMat.rows();
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: LORENTZ_GetGeneratorsAutom, beginning\n";
#endif
  check_correctness_lorentzian_perfect(LorMat, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: LORENTZ_GetGeneratorsAutom, after correctness checks\n";
#endif
  bool test = has_isotropic_factorization(LorMat);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: LORENTZ_GetGeneratorsAutom, test=" << test << "\n";
#endif
  if (test) {
    if (dim == 1) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: LORENTZ_GetGeneratorsAutom, before OneDimIsotropic_AutomGenerator\n";
#endif
      return OneDimIsotropic_AutomGenerator<T,Tint>(LorMat);
    }
    if (dim == 2) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: LORENTZ_GetGeneratorsAutom, before TwoDimIsotropic_AutomGenerator\n";
#endif
      return TwoDimIsotropic_AutomGenerator<T,Tint>(LorMat);
    }
  }
  return LORENTZ_GetGeneratorsAutom_Kernel<T,Tint,Tgroup>(LorMat, os);
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyVector<Tint>>
LORENTZ_GetOrbitRepresentative_Kernel(MyMatrix<T> const &LorMat, T const &X,
                               std::ostream &os) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
  MicrosecondTime time;
#endif
  int n = LorMat.rows();
  int dimEXT = n + 1;
  using TintGroup = typename Tgroup::Tint;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);

  int TheOption = LORENTZIAN_PERFECT_OPTION_TOTAL;

  DataPerfectLorentzian<T, Tint, Tgroup> data{n, LorMat, TheOption,
                                              std::move(rddo)};
  using Tdata = DataPerfectLorentzianFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  //
  auto f_incorrect = [&]([[maybe_unused]] Tobj const &x) -> bool {
    return false;
  };
  int max_runtime_second = 0;
  std::vector<Tout> l_obj =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);
  if (X != 0) {
    std::cerr << "LORPERF: Some code needs to be written for one sign\n";
    std::cerr
        << "LORPERF: For the other sign, that seems to require different methods\n";
    throw TerminalException{1};
  }
  size_t miss_val = std::numeric_limits<size_t>::max();
  std::vector<MyVector<Tint>> ListVect;
  if (X == 0) {
    struct SingleDataVector {
      Face f;
      MyVector<Tint> V;
    };
    struct DataVector {
      std::unordered_map<MyVector<Tint>, size_t> map_vert;
      std::vector<MyVector<Tint>> l_vert;
      std::vector<size_t> belonging;
      std::vector<SingleDataVector> l_data;
    };
    std::vector<DataVector> l_orbit;
    size_t n_entries = 0;
    std::vector<size_t> ListSize;
    std::vector<size_t> ListPos;
    std::vector<std::pair<size_t, size_t>> l_pair;
    std::unordered_map<std::pair<size_t, size_t>, size_t> map_pair;
    size_t iOrbCone = 0;
    size_t iOrbFull = 0;
    for (auto &eObj : l_obj) {
      vectface vf = DecomposeOrbitPoint_Full(eObj.x.GRP);
      size_t n_vert = eObj.x.EXT.rows();
      std::vector<size_t> belonging(n_vert, miss_val);
      std::unordered_map<MyVector<Tint>, size_t> map_vert;
      std::vector<MyVector<Tint>> l_vert;
      std::vector<SingleDataVector> l_data;
      for (size_t u = 0; u < n_vert; u++) {
        MyVector<Tint> V = GetMatrixRow(eObj.x.EXT, u);
        map_vert[V] = u + 1;
        l_vert.push_back(V);
      }
      size_t iOrb = 0;
      for (auto &eFace : vf) {
        boost::dynamic_bitset<>::size_type aRow = eFace.find_first();
        MyVector<Tint> V = GetMatrixRow(eObj.x.EXT, aRow);
        T sum = EvaluationQuadForm(LorMat, V);
        if (sum == 0) {
          SingleDataVector sdv{eFace, V};
          l_data.push_back(sdv);
          for (auto &pos : FaceToVector<size_t>(eFace)) {
            belonging[pos] = iOrb;
          }
          iOrb += 1;
        }
      }
      size_t eSize = l_data.size();
      for (size_t u = 0; u < eSize; u++) {
        std::pair<size_t, size_t> pair{iOrbCone, u};
        l_pair.push_back(pair);
        map_pair[pair] = iOrbFull;
        iOrbFull += 1;
      }
      ListSize.push_back(eSize);
      ListPos.push_back(n_entries);
      n_entries += eSize;
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: iOrbCone=" << iOrbCone
         << " belonging=" << StringStdVectorGAP(belonging) << "\n";
#endif
      DataVector dv{map_vert, l_vert, belonging, l_data};
      l_orbit.push_back(dv);
      iOrbCone += 1;
    }
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: We have build |l_orbit|=" << l_orbit.size()
       << " n_entries=" << n_entries << "\n";
#endif
    using Tgr = GraphBitset;
    Tgr eG(n_entries);
    for (size_t iOrbCone = 0; iOrbCone < l_orbit.size(); iOrbCone++) {
      std::vector<MyVector<Tint>> const &l_vert = l_orbit[iOrbCone].l_vert;
      size_t n_adj = l_obj[iOrbCone].ListAdj.size();
      for (size_t i_adj = 0; i_adj < n_adj; i_adj++) {
        auto &eAdj = l_obj[iOrbCone].ListAdj[i_adj];
        int iOrbConeAdj = eAdj.iOrb;
        TadjO &adj = eAdj.x;
        Face const &eInc = adj.eInc;
        MyMatrix<Tint> const &eBigMat = adj.eBigMat;
        MyMatrix<Tint> eBigMatInv = Inverse(eBigMat);
#ifdef DEBUG_LORENTZIAN_PERFECT
        os << "LORPERF: iOrbCone=" << iOrbCone << " i_adj=" << i_adj
           << " |eInc|=" << eInc.size() << " / " << eInc.count() << "\n";
#endif
        for (auto &eIdx : FaceToVector<size_t>(eInc)) {
          size_t iOrb = l_orbit[iOrbCone].belonging[eIdx];
          if (iOrb != miss_val) {
            std::pair<size_t, size_t> pair{iOrbCone, iOrb};
            size_t iOrbFull = map_pair.at(pair);
#ifdef DEBUG_LORENTZIAN_PERFECT_DISABLE
            os << "LORPERF: LorMat=\n";
            WriteMatrix(os, LorMat);
            os << "LORPERF: eBigMatInv=\n";
            WriteMatrix(os, eBigMatInv);
            os << "LORPERF: eIdx=" << eIdx << " l_vert[eIdx]=\n";
            WriteVector(os, l_vert[eIdx]);
#endif
            MyVector<Tint> eV = eBigMatInv.transpose() * l_vert[eIdx];
#ifdef DEBUG_LORENTZIAN_PERFECT_DISABLE
            os << "LORPERF: eV=\n";
            WriteVector(os, eV);
#endif
            size_t pos = l_orbit[iOrbConeAdj].map_vert[eV];
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
            if (pos == 0) {
              std::cerr
                  << "The vector is missing and that is now what we want\n";
              throw TerminalException{1};
            }
#endif
            size_t uVert = pos - 1;
            size_t iOrbAdj = l_orbit[iOrbConeAdj].belonging[uVert];
#ifdef DEBUG_LORENTZIAN_PERFECT_DISABLE
            os << "LORPERF iOrbCone=" << iOrbCone << " iOrb=" << iOrb
               << " iOrbConeAdj=" << iOrbConeAdj << " iOrbAdj=" << iOrbAdj
               << "\n";
#endif
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
            if (iOrbAdj == miss_val) {
              size_t n_vert = l_orbit[iOrbCone].belonging.size();
              size_t n_vert_adj = l_orbit[iOrbConeAdj].belonging.size();
              std::cerr << "LORPERF: |l_orbit|=" << l_orbit.size() << "\n";
              std::cerr << "LORPERF: iOrbCone=" << iOrbCone
                        << " iOrbConeAdj=" << iOrbConeAdj << "\n";
              std::cerr << "LORPERF: n_vert=" << n_vert << " n_vert_adj=" << n_vert_adj
                        << "\n";
              std::cerr << "LORPERF: Belonging(iOrbCone)="
                        << StringStdVectorGAP(l_orbit[iOrbCone].belonging)
                        << "\n";
              std::cerr << "LORPERF: Belonging(iOrbConeAdj)="
                        << StringStdVectorGAP(l_orbit[iOrbConeAdj].belonging)
                        << "\n";
              std::cerr << "LORPERF: eBigMat=\n";
              WriteMatrix(std::cerr, eBigMat);
              std::cerr
                  << "LORPERF: uVert should belong to an orbit of isotrop vertices\n";
              throw TerminalException{1};
            }
#endif
            std::pair<size_t, size_t> pairAdj{iOrbConeAdj, iOrbAdj};
            size_t iOrbFullAdj = map_pair.at(pairAdj);
#ifdef DEBUG_LORENTZIAN_PERFECT
            os << "LORPERF: Adjacency iOrbFull=" << iOrbFull
               << " iOrbFullAdj=" << iOrbFullAdj << "\n";
#endif
            if (iOrbFull != iOrbFullAdj) {
              eG.AddAdjacent(iOrbFull, iOrbFullAdj);
              eG.AddAdjacent(iOrbFullAdj, iOrbFull);
            }
          }
        }
      }
    }
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: We have eG\n";
#endif
    std::vector<std::vector<size_t>> ListConn = ConnectedComponents_set(eG);
    for (auto &eConn : ListConn) {
      size_t iOrbFull = eConn[0];
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: eConn =";
      for (auto &jOrbFull : eConn) {
        std::pair<size_t, size_t> const &pair = l_pair[jOrbFull];
        os << " (" << pair.first << "," << pair.second << ")";
      }
      os << "\n";
#endif
      std::pair<size_t, size_t> const &pair = l_pair[iOrbFull];
      size_t iOrbCone = pair.first;
      size_t iOrb = pair.second;
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: iOrbFull=" << iOrbFull << " iOrbCone=" << iOrbCone
         << " iOrb=" << iOrb << "\n";
#endif
      SingleDataVector const &sdv = l_orbit[iOrbCone].l_data[iOrb];
      MyVector<Tint> const &V = sdv.V;
      ListVect.push_back(V);
    }
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: |ListVect|=" << ListVect.size() << "\n";
#endif
  }
  if (X > 0) {
    // Need to write the copositive stuff
  }
#ifdef TIMINGS_LORENTZIAN_PERFECT
  os << "|LORPERF: LORENTZ_GetOrbitRepresentative|=" << time << "\n";
#endif
  return ListVect;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyVector<Tint>>
LORENTZ_GetOrbitRepresentative(MyMatrix<T> const &LorMat, T const &X,
                               std::ostream &os) {
  int dim = LorMat.rows();
  check_correctness_lorentzian_perfect(LorMat, os);
  if (has_isotropic_factorization(LorMat)) {
    if (dim == 1) {
      std::cerr << "We need to write down the code for this case\n";
      throw TerminalException{1};
    }
    if (dim == 2) {
      return TwoDimIsotropic_OrbitRepresentative<T,Tint>(LorMat, X);
    }
  }
  return LORENTZ_GetOrbitRepresentative_Kernel<T,Tint,Tgroup>(LorMat, X, os);
}

template <typename T>
size_t INDEF_FORM_Invariant_NonDeg(MyMatrix<T> const &SymMat, size_t seed,
                                   std::ostream &os) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
  MicrosecondTime time;
#endif
  int n = SymMat.rows();
  T det = DeterminantMat<T>(SymMat);
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(SymMat, os);
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  size_t hash_ret = seed;
  size_t hash_det = std::hash<T>()(det);
  combine_hash(hash_ret, hash_det);
  size_t hash_n = n;
  combine_hash(hash_ret, hash_n);
  size_t hash_minus = DiagInfo.nbMinus;
  combine_hash(hash_ret, hash_minus);
#ifdef TIMINGS_LORENTZIAN_PERFECT
  os << "|LORPERF: INDEF_FORM_Invariant_NonDeg|=" << time << "\n";
#endif
  return hash_ret;
}

template<typename Tint>
struct ResultEquivLorentzian {
  MyMatrix<Tint> equiv;
  std::vector<MyMatrix<Tint>> list_gen;
};


template <typename T, typename Tint, typename Tgroup>
std::optional<ResultEquivLorentzian<Tint>>
LORENTZ_TestEquivalenceMatrices_Reduced(MyMatrix<T> const &LorMat1,
                                        MyMatrix<T> const &LorMat2, std::ostream &os) {
#ifdef TIMINGS_LORENTZIAN_PERFECT
  MicrosecondTime time;
#endif
  size_t seed = 1678;
  size_t inv1 = INDEF_FORM_Invariant_NonDeg(LorMat1, seed, os);
  size_t inv2 = INDEF_FORM_Invariant_NonDeg(LorMat2, seed, os);
  if (inv1 != inv2) {
    return {};
  }
  int n = LorMat1.rows();
  int dimEXT = n + 1;
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  int TheOption = LORENTZIAN_PERFECT_OPTION_TOTAL;
  //
  LorentzianPerfectEntry<T, Tint> eRec2 =
      LORENTZ_GetOnePerfect<T, Tint>(LorMat2, TheOption, os);
  MyMatrix<Tint> EXT2 = LORENTZ_GetEXT(eRec2);
  //
  DataPerfectLorentzian<T, Tint, Tgroup> data{n, LorMat1, TheOption,
                                              std::move(rddo)};
  using Tdata = DataPerfectLorentzianFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_Serial<Tobj, TadjO>;
  std::optional<MyMatrix<Tint>> opt;
  auto f_incorrect = [&](Tobj const &x) -> bool {
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: TestEquivalenceMatrices: Before LORENTZ_TestEquivalence\n";
#endif
    std::optional<MyMatrix<Tint>> opt_res =
        LORENTZ_TestEquivalence<T, Tint, Tgroup>(LorMat1, LorMat2, x.EXT, EXT2,
                                                 os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: TestEquivalenceMatrices: After LORENTZ_TestEquivalence\n";
#endif
    if (opt_res) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "LORPERF: TestEquivalenceMatrices: Matching opt_res\n";
#endif
      MyMatrix<Tint> eBigMat = Inverse(*opt_res);
      opt = eBigMat;
#ifdef DEBUG_LORENTZIAN_PERFECT
      MyMatrix<T> eBigMat_T = UniversalMatrixConversion<T, Tint>(eBigMat);
      MyMatrix<T> eProd = eBigMat_T * LorMat1 * eBigMat_T.transpose();
      if (eProd != LorMat2) {
        std::cerr << "LORPERF: It is not actually an equivalence\n";
        throw TerminalException{1};
      }
      os << "LORPERF: We found a matching entry. Exiting the enumeration\n";
#endif
      return true;
    }
    return false;
  };
  int max_runtime_second = 0;
  std::vector<Tout> l_obj =
    EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
      data_func, f_incorrect, max_runtime_second);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: After EnumerateAndStore_Serial function call\n";
#endif
#ifdef TIMINGS_LORENTZIAN_PERFECT
  os << "|LORPERF: LORENTZ_TestEquivalenceMatrices|=" << time << "\n";
#endif
  std::vector<MyMatrix<Tint>> list_gen =
      LORENTZ_ExtractGeneratorsFromObjList<Tint, Telt>(l_obj);
  if (opt) {
    MyMatrix<Tint> const& equiv = *opt;
    ResultEquivLorentzian<Tint> res{equiv, list_gen};
    return res;
  }
  return {};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LORENTZ_TestEquivalenceMatrices_Kernel(MyMatrix<T> const &LorMat1,
                                       MyMatrix<T> const &LorMat2, std::ostream &os) {
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: Before IndefiniteReduction, LorMat1=\n";
  WriteMatrix(os, LorMat1);
  os << "LORPERF: LorMat2=\n";
  WriteMatrix(os, LorMat2);
#endif
  ResultReduction<T, Tint> res1 =
    IndefiniteReduction<T,Tint>(LorMat1, os);
  ResultReduction<T, Tint> res2 =
    IndefiniteReduction<T,Tint>(LorMat2, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: After IndefiniteReduction, res1.Mred=\n";
  WriteMatrix(os, res1.Mred);
  os << "LORPERF: res2.Mred=\n";
  WriteMatrix(os, res2.Mred);
#endif
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
  auto f_test=[&](MyMatrix<Tint> const& eEquivRet) -> void {
    MyMatrix<T> eEquivRet_T = UniversalMatrixConversion<T,Tint>(eEquivRet);
    MyMatrix<T> eProd = eEquivRet_T * LorMat1 * eEquivRet_T.transpose();
    if (eProd != LorMat2) {
      std::cerr << "LORPERF: LorMat1 should map to LorMat2\n";
      throw TerminalException{1};
    }
  };
#endif
  // We have res1.B * LorMat1 * res1.B^T = res1.Mred
  //     and res2.B * LorMat2 * res2.B^T = res2.Mred
  std::optional<ResultEquivLorentzian<Tint>> opt = LORENTZ_TestEquivalenceMatrices_Reduced<T,Tint,Tgroup>(res1.Mred, res2.Mred, os);
  if (opt) {
    MyMatrix<Tint> B2_inv = Inverse(res2.B);
    MyMatrix<Tint> B1_inv = Inverse(res1.B);
    ResultEquivLorentzian<Tint> const& res = *opt;
    // Using the generators that we found. This is not the complete list
    // but they are bona fide generators and so can be used to get a smaller
    // equivalence.
    std::vector<MyMatrix<Tint>> list_gen1;
    for (auto & eGen: res.list_gen) {
      MyMatrix<Tint> eGen_B = B1_inv * eGen * res1.B;
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
      MyMatrix<T> eGen_BT = UniversalMatrixConversion<T,Tint>(eGen_B);
      MyMatrix<T> eProd = eGen_BT * LorMat1 * eGen_BT.transpose();
      if (eProd != LorMat1) {
        std::cerr << "LORPERF: eGen_B should preserve the LorMat1\n";
        throw TerminalException{1};
      }
#endif
      list_gen1.push_back(eGen_B);
    }
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: We have |list_gen1|=" << list_gen1.size() << "\n";
#endif
    std::vector<MyMatrix<Tint>> list_gen_red1 = ExhaustiveReductionComplexityGroupMatrix(list_gen1, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: We have |list_gen_red1|=" << list_gen_red1.size() << "\n";
#endif
    // Building an equivalence
    MyMatrix<Tint> const& eEquiv = res.equiv;
    // eEquiv * res1.Mred * eEquiv^T = res2.Mred
    MyMatrix<Tint> eEquivRet = B2_inv * eEquiv * res1.B;
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: We have eEquivRet=\n";
    WriteMatrix(os, eEquivRet);
#endif
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
    f_test(eEquivRet);
#endif
    MyMatrix<Tint> eEquivRet_red = ExhaustiveMatrixRightCosetSimplification(eEquivRet, list_gen_red1);
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORPERF: We have eEquivRet_red=\n";
    WriteMatrix(os, eEquivRet_red);
#endif
#ifdef SANITY_CHECK_LORENTZIAN_PERFECT
    f_test(eEquivRet_red);
#endif
    return eEquivRet_red;
  }
  return {};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LORENTZ_TestEquivalenceMatrices(MyMatrix<T> const &LorMat1,
                                MyMatrix<T> const &LorMat2, std::ostream &os) {
  int dim = LorMat1.rows();
  check_correctness_lorentzian_perfect(LorMat1, os);
  check_correctness_lorentzian_perfect(LorMat2, os);
  bool test1 = has_isotropic_factorization(LorMat1);
  bool test2 = has_isotropic_factorization(LorMat2);
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: LORENTZ_TestEquivalenceMatrices, test1=" << test1 << " test2=" << test2 << "\n";
#endif
  if (test1 != test2) {
    return {};
  }
#ifdef DEBUG_LORENTZIAN_PERFECT
  os << "LORPERF: LORENTZ_TestEquivalenceMatrices, before the calls with dim=" << dim << "\n";
#endif
  if (test1) {
    if (dim == 1) {
      return OneDimIsotropic_TestEquivalence<T,Tint>(LorMat1, LorMat2);
    }
    if (dim == 2) {
      return TwoDimIsotropic_TestEquivalence<T,Tint>(LorMat1, LorMat2);
    }
  }
  return LORENTZ_TestEquivalenceMatrices_Kernel<T,Tint,Tgroup>(LorMat1, LorMat2, os);
}

// clang-format off
#endif  // SRC_LORENTZIAN_LORENTZIAN_PERFECT_H_
// clang-format on
