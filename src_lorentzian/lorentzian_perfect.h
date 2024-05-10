// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_PERFECT_H_
#define SRC_LORENTZIAN_PERFECT_H_

// clang-format off
#include "boost_serialization.h"
#include "FundamentalDelaunay.h"
#include "POLY_RecursiveDualDesc.h"
#include "POLY_AdjacencyScheme.h"
#include "fractions.h"
#include <map>
#include <string>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LORENTZIAN_PERFECT
#endif

static const int LORENTZIAN_PERFECT_OPTION_ISOTROP = 23;
static const int LORENTZIAN_PERFECT_OPTION_TOTAL = 47;

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
LORENTZ_FindPositiveVectors(MyMatrix<T> const &LorMat, MyVector<T> const &eVect,
                            T const &MaxScal, int const &TheOption,
                            bool const &OnlyShortest, std::ostream &os) {
  int n = LorMat.rows();
  T eNorm = EvaluationQuadForm(LorMat, eVect);
#ifdef DEBUG_LORENTZIAN_PERFECT
  if (MaxScal <= 0) {
    std::cerr << "MaxScal=" << MaxScal << "\n";
    throw TerminalException{1};
  }
  if (TheOption != LORENTZIAN_PERFECT_OPTION_ISOTROP &&
      TheOption != LORENTZIAN_PERFECT_OPTION_TOTAL) {
    std::cerr << "TheOption=" << TheOption << "\n";
    throw TerminalException{1};
  }
  if (eNorm <= 0) {
    std::cerr << "eNorm=" << eNorm << "\n";
    std::cerr << "Wrong norm of vector, will not work\n";
    throw TerminalException{1};
  }
  if (!IsIntegralVector(eVect)) {
    std::cerr << "eVect should be integral\n";
    throw TerminalException{1};
  }
  if (!IsIntegralMatrix(LorMat)) {
    std::cerr << "LorMat should be integral\n";
    throw TerminalException{1};
  }
#endif
  MyVector<Tint> eVect_tint = UniversalVectorConversion<Tint, T>(eVect);
  MyVector<T> eVect_LorMat = LorMat * eVect;
  FractionVector<T> eRec = RemoveFractionVectorPlusCoeff(eVect_LorMat);
  MyVector<Tint> eVect_LorMat_tint =
      UniversalVectorConversion<Tint, T>(eRec.TheVect);
  MyMatrix<Tint> eVect_LorMat_tint_M(1, n);
  for (int u = 0; u < n; u++) {
    eVect_LorMat_tint_M(0, u) = eVect_LorMat_tint(u);
  }
  MyMatrix<Tint> Ubasis = NullspaceIntMat(eVect_LorMat_tint_M);
  MyMatrix<T> Ubasis_T = UniversalMatrixConversion<T, Tint>(Ubasis);
  MyMatrix<T> GramMat = -Ubasis_T * LorMat * Ubasis_T.transpose();
  CVPSolver<T, Tint> solver(GramMat, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
  if (!IsPositiveDefinite(GramMat)) {
    std::cerr << "GramMat should be positive definite\n";
    throw TerminalException{1};
  }
#endif
  GCD_dot<Tint> TheRec_pre = ComputeGcdDot(eVect_LorMat_tint);
  GCD_dot<Tint> TheRec = PositivityNormalizeGcdDot(TheRec_pre);
  std::vector<MyVector<Tint>> TotalListSol;
  Tint eVal = 1;
  while (true) {
    MyVector<Tint> eBasSol = eVal * TheRec.V; // A solution of
    MyVector<T> eBasSol_T = UniversalVectorConversion<T, Tint>(eBasSol);
    T alpha = eVal * TheRec.gcd / eNorm;
    MyVector<T> eTrans = eBasSol_T - alpha * eVect;
    std::optional<MyVector<T>> opt = SolutionMat(Ubasis_T, eTrans);
    MyVector<T> eSol = unfold_opt(opt, "Getting eSol");
    T eSquareDist = alpha * alpha * eNorm;
    auto iele = [&]() -> std::vector<MyVector<Tint>> {
      if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
        return solver.FixedNormVectors(eSol, eSquareDist);
      } else {
        return solver.AtMostNormVectors(eSol, eSquareDist);
      }
    };
    for (auto &eSolA : iele()) {
      MyVector<Tint> eSolC = eBasSol + Ubasis.transpose() * eSolA;
      TotalListSol.emplace_back(std::move(eSolC));
    }
    if (OnlyShortest && TotalListSol.size()) {
      return TotalListSol;
    }
    eVal += 1;
    if (MaxScal > 0) {
      T scal = TheRec.gcd * eVal;
      if (scal > MaxScal) {
        break;
      }
    }
  }
  return TotalListSol;
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
LORENTZ_SearchInitialVector(MyMatrix<T> const &LorMat,
                            MyVector<T> const &PosVect, int const &TheOption,
                            std::ostream &os) {
  bool OnlyShortest = true;
  T MaxScal = 0;
  return LORENTZ_FindPositiveVectors<T, Tint>(LorMat, PosVect, MaxScal,
                                              TheOption, OnlyShortest, os);
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
  if (DeltaB <= 0) {
    std::cerr << "If the discriminant is negative, then we cannot have "
                 "signature (1,1)";
    throw TerminalException{1};
  }
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
                               MyVector<T> const &eNSPdir) {
  MyMatrix<T> LorMatInv = Inverse(LorMat);
  T eCstBas = eNSPbas(0);
  T eCstDir = eNSPdir(0);
  MyVector<T> eBas = GetReducedVector(eNSPbas);
  MyVector<T> eDir = GetReducedVector(eNSPdir);
  // For an acceptable vector w of length n+1 in return we must have w[0] < 0.
  // Since we have w = u + TheMult * v we have a potential upper bound
  // on TheMult, but only if v[0] > 0
  std::vector<T> ListUpperBound;
  std::optional<T> UpperBound_constant;
  if (eCstDir > 0) {
    UpperBound_constant = -eCstBas / eCstDir;
    if (*UpperBound_constant <= 0) {
      std::cerr << "The upper bound from constant is not correct\n";
      throw TerminalException{1};
    }
    ListUpperBound.push_back(*UpperBound_constant);
  }
  //
  // Get raw upper bound
  //
  T iShift = 1;
  while (true) {
    MyVector<T> eV = eBas + iShift * eDir;
    MyVector<T> eVect = LorMatInv * eV;
    T eNorm = EvaluationQuadForm<T, T>(LorMat, eVect);
    if (eNorm < 0) {
      ListUpperBound.push_back(iShift);
      break;
    }
    iShift *= 2;
  }
  //
  // More subttle upper bound coming from isotropy computation
  //
  MyVector<T> eVectBas = LorMatInv * eBas;
  MyVector<T> eVectDir = LorMatInv * eDir;
#ifdef DEBUG_LORENTZIAN_PERFECT
  auto f_norm = [&](T const &val) -> T {
    MyVector<T> eVect = eVectBas + val * eVectDir;
    return EvaluationQuadForm<T, T>(LorMat, eVect);
  };
#endif
  std::vector<MyVector<T>> TheBasis{eVectBas, eVectDir};
  MyMatrix<T> TheBasis_mat = MatrixFromVectorFamily(TheBasis);
  MyMatrix<T> LorMat22 = TheBasis_mat * LorMat * TheBasis_mat.transpose();
  std::vector<MyVector<T>> ListIso = GetRationalIsotropyVectors(LorMat22);
  std::optional<T> UpperBound_isotropic;
  std::vector<T> ListUpperBound_Iso;
  for (auto &eIso : ListIso) {
    if (eIso(0) != 0) {
      T fact = eIso(1) / eIso(0);
#ifdef DEBUG_LORENTZIAN_PERFECT
      if (f_norm(fact) != 0) {
        std::cerr << "eVect should be isotropic\n";
        throw TerminalException{1};
      }
#endif
      if (fact > 0) {
        ListUpperBound_Iso.push_back(fact);
      }
    }
  }
  if (ListUpperBound_Iso.size() > 0) {
    UpperBound_isotropic = VectorMin(ListUpperBound_Iso);
    ListUpperBound.push_back(*UpperBound_isotropic);
  }
  if (UpperBound_constant.has_value() && UpperBound_isotropic.has_value()) {
    if (*UpperBound_constant == *UpperBound_isotropic) {
      return {};
    }
  }
  // Need to see if better upper bounds are possible, but this is a secondary
  // question
  T BestUpper = VectorMin(ListUpperBound);
#ifdef DEBUG_LORENTZIAN_PERFECT
  if (f_norm(BestUpper) > 0) {
    std::cerr << "We should have eNorm <= 0\n";
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
    int const &TheOption, [[maybe_unused]] std::ostream &os) {
  MyMatrix<T> LorMatInv = Inverse(LorMat);
#ifdef DEBUG_LORENTZIAN_PERFECT
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
    if (xext_T.dot(eNSPbas) != 0) {
      std::cerr << "eNSPbas should have scalar product 0 with all entries in "
                   "CritSet\n";
      throw TerminalException{1};
    }
    if (xext_T.dot(eNSPdir) != 0) {
      std::cerr << "eNSPdir should have scalar product 0 with all entries in "
                   "CritSet\n";
      throw TerminalException{1};
    }
    if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
      T norm = EvaluationQuadForm<T, Tint>(LorMat, x);
      if (norm != 0) {
        std::cerr << "CritSet contains some non-isotrop vectors\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  bool OnlyShortest = true;
  T TheLowerBound = 0;
  std::optional<T> TheUpperBound_opt = GetUpperBound(LorMat, eNSPbas, eNSPdir);
  if (!TheUpperBound_opt) {
    std::cerr << "We have TheUpperBound that is none and that is a problem\n";
    throw TerminalException{1};
  }
  T TheUpperBound = *TheUpperBound_opt;
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
    os << "TheLowerBound=" << TheLowerBound
       << " TheUpperBound=" << TheUpperBound << " n_iter=" << n_iter << "\n";
#endif
    T TheMidVal = get_mid_val(TheLowerBound, TheUpperBound);
    MyVector<T> eNSPtest = eNSPbas + TheMidVal * eNSPdir;
    MyVector<T> eVectTest = LorMatInv * GetReducedVector(eNSPtest);
    T eNormTest = EvaluationQuadForm<T, T>(LorMat, eVectTest);
    MyVector<T> CritSet0_T = UniversalVectorConversion<T, Tint>(CritSet[0]);
    T MaxScal = ScalarProductQuadForm(LorMat, CritSet0_T, eVectTest);
    if (eNormTest <= 0 && MaxScal <= 0) {
      TheUpperBound = TheMidVal;
    } else {
      ListTotal = LORENTZ_FindPositiveVectors<T, Tint>(
          LorMat, eVectTest, MaxScal, TheOption, OnlyShortest, os);
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "|ListTotal|=" << ListTotal.size() << "\n";
      if (IsSubset(CritSet, ListTotal) && CritSet.size() > ListTotal.size()) {
        std::cerr << "Bug: if included, it should be equal\n";
        throw TerminalException{1};
      }
#endif
      if (IsEqualSet(ListTotal, CritSet)) {
        TheLowerBound = TheMidVal;
      } else {
        if (IsSubset(ListTotal, CritSet)) {
#ifdef DEBUG_LORENTZIAN_PERFECT
          os << "EXIT 1 |ListTotal|=" << ListTotal.size()
             << " MaxScal=" << MaxScal << "\n";
          os << "  NSPtest=" << StringVectorGAP(eNSPtest) << "\n";
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
  os << "Going to the second scheme\n";
#endif
  while (true) {
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
    std::vector<MyVector<Tint>> ListTotal =
        LORENTZ_FindPositiveVectors<T, Tint>(LorMat, eVectTest, MaxScal,
                                             TheOption, OnlyShortest, os);
    if (IsSubset(ListTotal, CritSet)) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      os << "EXIT 2 |ListTotal|=" << ListTotal.size() << " MaxScal=" << MaxScal
         << "\n";
      os << "  NSPtest=" << StringVectorGAP(eNSPtest) << "\n";
#endif
      return {ListTotal, eNSPtest, eVectTest, MaxScal};
    }
  }
}

template <typename T, typename Tint>
MyVector<T> GetOneOutsideRay(MyMatrix<T> const &LorMat,
                             MyMatrix<T> const &SpannBasis,
                             std::vector<MyVector<Tint>> const &TheSet,
                             MyVector<T> const &eNSPbas, std::ostream &os) {
  MyVector<T> TheSet0 = UniversalVectorConversion<T, Tint>(TheSet[0]);
  MyMatrix<T> TheMat = SpannBasis * LorMat * SpannBasis.transpose();
#ifdef DEBUG_LORENTZIAN_PERFECT
  DiagSymMat<T> dsm = DiagonalizeSymmetricMatrix(TheMat);
  if (dsm.nbMinus == 0) {
    std::cerr << "We should have a negative in the entry\n";
    throw TerminalException{1};
  }
#endif
  int n_dim = TheMat.rows();
  MyMatrix<T> ePerturb = IdentityMat<T>(n_dim);
  while (true) {
    MyMatrix<T> TheMatPerturb = -ePerturb * TheMat * ePerturb.transpose();
    MyMatrix<T> uVect =
        GetIntegralPositiveVector_allmeth<T, T>(TheMatPerturb, os);
    MyMatrix<T> SpannMatrixPert = ePerturb * SpannBasis;
    MyVector<T> RetVect = SpannMatrixPert.transpose() * uVect;
#ifdef DEBUG_LORENTZIAN_PERFECT
    T TheNorm = EvaluationQuadForm<T, T>(LorMat, RetVect);
    if (TheNorm == 0) {
      std::cerr << "The vector should be outside of the cone and so have "
                   "negative norm\n";
      throw TerminalException{1};
    }
#endif
    MyVector<T> tVect = LorMat * RetVect;
    T eScal = tVect.dot(TheSet0);
#ifdef DEBUG_LORENTZIAN_PERFECT
    for (auto &eVect : TheSet) {
      MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect);
      T fScal = tVect.dot(eVect_T);
      if (eScal != fScal) {
        std::cerr << "The scalar products are incorrect\n";
        throw TerminalException{1};
      }
    }
#endif
    MyVector<T> eNSPdir = ConcatenateScalarVector(T(-eScal), tVect);
    std::optional<T> TheUpperBound_opt =
        GetUpperBound(LorMat, eNSPbas, eNSPdir);
    if (TheUpperBound_opt.has_value()) {
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
    ListVectExt(i_vect, 0);
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
    std::cerr << "We have rnk=" << rnk << " n=" << n
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
  int n = LorMat.rows();
  MyMatrix<T> LorMatInv = Inverse(LorMat);
  MyVector<Tint> CentralVect =
      INDEFINITE_GetShortPositiveVector<T, Tint>(LorMat, os);
  MyVector<T> CentralVect_T = UniversalVectorConversion<T, Tint>(CentralVect);
  std::vector<MyVector<Tint>> CritSet = LORENTZ_SearchInitialVector<T, Tint>(
      LorMat, CentralVect_T, TheOption, os);
  MyVector<T> CritSet0_T = UniversalVectorConversion<T, Tint>(CritSet[0]);
  MyVector<T> LorMat_Central = LorMat * CentralVect_T;
  T eScal = LorMat_Central.dot(CritSet0_T);
  MyVector<T> eNSPbas = ConcatenateScalarVector(T(-eScal), LorMat_Central);
  while (true) {
    int rnk = 0;
    if (CritSet.size() > 0) {
      rnk = RankMat(MatrixFromVectorFamily(CritSet));
    }
#ifdef DEBUG_LORENTZIAN_PERFECT
    os << "LORENTZ_GetOnePerfect rnk=" << rnk << " |CritSet|=" << CritSet.size()
       << "\n";
#endif
    if (rnk == n) {
#ifdef DEBUG_LORENTZIAN_PERFECT
      LORENTZ_CheckCorrectnessVectorFamily(LorMat, CritSet);
#endif
      return {CritSet, eNSPbas, CentralVect_T};
    }
    MyMatrix<T> EXT = GetFullExpanded<T, Tint>(CritSet);
    MyMatrix<T> NSP = NullspaceTrMat(EXT);
#ifdef DEBUG_LORENTZIAN_PERFECT
    if (NSP.rows() == 0) {
      std::cerr << "NSP should be non-empty\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<T> ListDir = DropColumn(NSP, 0) * LorMatInv;
    MyMatrix<T> eNSPdir =
        GetOneOutsideRay<T, Tint>(LorMat, ListDir, CritSet, eNSPbas, os);
    ResultFlipping<T, Tint> eRecB = LORENTZ_Kernel_Flipping<T, Tint>(
        LorMat, CritSet, eNSPbas, eNSPdir, TheOption, os);
    CritSet = eRecB.ListTotal;
    eNSPbas = eRecB.eNSPtest;
  }
}

template <typename T, typename Tint>
LorentzianPerfectEntry<T, Tint>
LORENTZ_DoFlipping(MyMatrix<T> const &LorMat,
                   std::vector<MyVector<Tint>> const &ListIso, Face eInc,
                   int const &TheOption, std::ostream &os) {
  size_t n_vect = eInc.size();
  MyMatrix<T> EXT = GetFullExpanded<T, Tint>(ListIso);
  auto get_eVert = [&]() -> size_t {
    for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
      if (eInc[i_vect] == 0) {
        return i_vect;
      }
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  size_t eVert = get_eVert();
  std::vector<MyVector<T>> ListIsoSel;
  for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
    if (eInc[i_vect] == 1) {
      MyVector<T> eVect = UniversalVectorConversion<T, Tint>(ListIso[i_vect]);
      ListIsoSel.push_back(eVect);
    }
  }
  MyMatrix<T> MatrIsoSel = MatrixFromVectorFamily(ListIsoSel);
  MyMatrix<T> NSP = NullspaceTrMat(MatrIsoSel);
#ifdef DEBUG_LORENTZIAN_PERFECT
  if (NSP.rows() != 1) {
    std::cerr << "NSP should have size 1\n";
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
  MyVector<T> NSPb_0 = GetMatrixRow(NSPb, 0);
  MyVector<T> eVectB = GetReducedVector(NSPb_0);
  auto get_eNSPbas = [&]() -> MyVector<T> {
    if (eVectB.dot(eIso) > 0) {
      return NSPb_0;
    } else {
      return -NSPb_0;
    }
  };
  MyVector<T> eNSPbas = get_eNSPbas();
  std::vector<MyVector<Tint>> CritSet;
  for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
    if (eInc[i_vect] == 1) {
      CritSet.push_back(ListIso[i_vect]);
    }
  }
  std::vector<MyVector<Tint>> TheFlip =
      LORENTZ_Kernel_Flipping<T, Tint>(LorMat, CritSet, eNSPbas, eNSPdir,
                                       TheOption, os)
          .ListTotal;
#ifdef DEBUG_LORENTZIAN_PERFECT
  LORENTZ_CheckCorrectnessVectorFamily(LorMat, TheFlip);
#endif
  // No return of additional info so far.
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
  return {ListGen, GRPperm};
}

template <typename T, typename Tint>
size_t ComputeInvariantPerfectForm(size_t seed, MyMatrix<T> const &LorMat,
                                   MyMatrix<Tint> const &EXT,
                                   [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_DELAUNAY_ENUMERATION
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
      T sum = 0;
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
      T scal = 0;
      for (int i = 0; i < n; i++) {
        scal += V(i) * EXT_T(jRow, i);
      }
      ListOffDiagNorm[scal] += 1;
    }
  }
  size_t hash = ComputeHashTwoMap(seed, ListDiagNorm, ListOffDiagNorm);
#ifdef TIMINGS_DELAUNAY_ENUMERATION
  os << "|ComputeInvariantDelaunay|=" << time << "\n";
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
    LorentzianPerfectEntry<T, Tint> eRec = LORENTZ_GetOnePerfect<T, Tint>(
        data.LorMat, data.TheOption, data.rddo.os);
    MyMatrix<Tint> EXT = LORENTZ_GetEXT(eRec);
    Tobj x{std::move(EXT), {}};
    return x;
  }
  size_t f_hash(size_t const &seed, Tobj const &x) {
    return ComputeInvariantPerfectForm(seed, data.LorMat, x.EXT, data.rddo.os);
  }
  std::optional<TadjO> f_repr(Tobj const &x, TadjI const &y) {
    std::optional<MyMatrix<Tint>> opt =
        LORENTZ_TestEquivalence<T, Tint, Tgroup>(data.LorMat, data.LorMat,
                                                 x.EXT, y.EXT, data.rddo.os);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const &eBigMat = *opt;
    TadjO ret{y.eInc, eBigMat};
    return ret;
  }
  std::pair<Tobj, TadjO> f_spann(TadjI const &x) {
    MyMatrix<Tint> EXT = x.EXT;
    Tobj x_ret{EXT, {}};
    MyMatrix<Tint> eBigMat = IdentityMat<Tint>(data.n + 1);
    TadjO ret{x.eInc, eBigMat};
    return {x_ret, ret};
  }
  std::vector<TadjI> f_adj(Tobj &x) {
    MyMatrix<T> EXT_T = UniversalMatrixConversion<T, Tint>(x.EXT);
    ResultStabilizer<Tint, Tgroup> res_stab =
        LORENTZ_ComputeStabilizer<T, Tint, Tgroup>(data.LorMat, x.EXT,
                                                   data.rddo.os);
    vectface TheOutput =
        DualDescriptionRecordFullDim(EXT_T, res_stab.GRPperm, data.rddo);
    x.GRP = res_stab.GRPperm;
    std::vector<MyVector<Tint>> ListIso;
    for (int i_row = 0; i_row < x.EXT.rows(); i_row++) {
      MyVector<Tint> V = GetMatrixRow(x.EXT, i_row);
      ListIso.push_back(V);
    }
    std::vector<TadjI> ListAdj;
    for (auto &eInc : TheOutput) {
      LorentzianPerfectEntry<T, Tint> eRecFlip = LORENTZ_DoFlipping<T, Tint>(
          data.LorMat, ListIso, eInc, data.TheOption, data.rddo.os);
      MyMatrix<Tint> EXTflip = LORENTZ_GetEXT(eRecFlip);
      TadjI eAdj{eInc, EXTflip};
      ListAdj.push_back(eAdj);
    }
    return ListAdj;
  }
  Tobj f_adji_obj(TadjI const &x) {
    MyMatrix<Tint> EXT = x.EXT;
    Tobj x_ret{EXT, {}};
    return x_ret;
  };
};

FullNamelist NAMELIST_GetStandard_COMPUTE_PERFECT_LORENTZIAN() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string>> ListListStringValues1;
  ListStringValues1["arithmetic_T"] = "gmp_rational";
  ListStringValues1["arithmetic_Tint"] = "gmp_integer";
  ListStringValues1["LorMatFile"] = "unset.gram";
  ListStringValues1["Option"] = "unset";
  ListStringValues1["OutFormat"] = "nothing";
  ListStringValues1["OutFile"] = "unset.out";
  ListStringValues1["FileDualDescription"] = "unset";
  ListIntValues1["max_runtime_second"] = 0;
  ListBoolValues1["ApplyStdUnitbuf"] = false;
  SingleBlock BlockDATA;
  BlockDATA.ListIntValues = ListIntValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListDoubleValues = ListDoubleValues1;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListListStringValues = ListListStringValues1;
  ListBlock["DATA"] = BlockDATA;
  // STORAGE
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<std::string>> ListListStringValues2;
  ListBoolValues2["Saving"] = false;
  ListStringValues2["Prefix"] = "/irrelevant/";
  SingleBlock BlockSTORAGE;
  BlockSTORAGE.ListIntValues = ListIntValues2;
  BlockSTORAGE.ListBoolValues = ListBoolValues2;
  BlockSTORAGE.ListDoubleValues = ListDoubleValues2;
  BlockSTORAGE.ListStringValues = ListStringValues2;
  BlockSTORAGE.ListListStringValues = ListListStringValues2;
  ListBlock["STORAGE"] = BlockSTORAGE;
  // Merging all data
  return {ListBlock, "undefined"};
}

template <typename T, typename Tint, typename Tgroup>
void ComputePerfectLorentzian(boost::mpi::communicator &comm,
                              FullNamelist const &eFull) {
  std::unique_ptr<std::ofstream> os_ptr = get_mpi_log_stream(comm, eFull);
  std::ostream &os = *os_ptr;
  SingleBlock BlockDATA = eFull.ListBlock.at("DATA");
  SingleBlock BlockSTORAGE = eFull.ListBlock.at("STORAGE");
  //
  bool STORAGE_Saving = BlockSTORAGE.ListBoolValues.at("Saving");
  std::string STORAGE_Prefix = BlockSTORAGE.ListStringValues.at("Prefix");
  CreateDirectory(STORAGE_Prefix);
  //
  int max_runtime_second = BlockDATA.ListIntValues.at("max_runtime_second");
  std::cerr << "max_runtime_second=" << max_runtime_second << "\n";
  std::string LorMatFile = BlockDATA.ListStringValues.at("LorMatFile");
  MyMatrix<T> LorMat = ReadMatrixFile<T>(LorMatFile);
  //
  std::string TheOption_str = BlockDATA.ListStringValues.at("TheOption");
  auto get_option = [&]() -> int {
    if (TheOption_str == "isotropic") {
      return LORENTZIAN_PERFECT_OPTION_ISOTROP;
    }
    if (TheOption_str == "total") {
      return LORENTZIAN_PERFECT_OPTION_TOTAL;
    }
    std::cerr << "Failed to find a matching entry for TheOption\n";
    throw TerminalException{1};
  };
  int TheOption = get_option();
  //
  std::string OutFormat = BlockDATA.ListStringValues.at("OutFormat");
  std::string OutFile = BlockDATA.ListStringValues.at("OutFile");
  std::cerr << "OutFormat=" << OutFormat << " OutFile=" << OutFile << "\n";
  //
  int n = LorMat.rows();
  int dimEXT = n + 1;
  using TintGroup = typename Tgroup::Tint;
  std::string FileDualDesc =
      BlockDATA.ListStringValues.at("FileDualDescription");
  PolyHeuristicSerial<TintGroup> AllArr =
      Read_AllStandardHeuristicSerial_File<TintGroup>(FileDualDesc, dimEXT, os);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, os);
  //
  DataPerfectLorentzian<T, Tint, Tgroup> data{n, LorMat, TheOption,
                                              std::move(rddo)};
  using Tdata = DataPerfectLorentzianFunc<T, Tint, Tgroup>;
  Tdata data_func{std::move(data)};
  using Tobj = typename Tdata::Tobj;
  using TadjO = typename Tdata::TadjO;
  using Tout = DatabaseEntry_MPI<Tobj, TadjO>;
  //
  std::pair<bool, std::vector<Tout>> pair = EnumerateAndStore_MPI<Tdata>(
      comm, data_func, STORAGE_Prefix, STORAGE_Saving, max_runtime_second);
#ifdef DEBUG_DELAUNAY_ENUMERATION
  os << "PERF_LOR: We now have IsFinished=" << pair.first << "\n";
  os << "PERF_LOR: We now have |ListPerfect|=" << pair.second.size() << "\n";
#endif
  //
  if (pair.first) {
    WriteFamilyObjects<Tobj, TadjO>(comm, OutFormat, OutFile, pair.second, os);
  }
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
LORENTZ_GetGeneratorsAutom(MyMatrix<T> const &LorMat, std::ostream &os) {
  int n = LorMat.rows();
  int dimEXT = n + 1;
  using TintGroup = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<Tint>(dimEXT, os);
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
  std::optional<std::vector<Tout>> opt =
      EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
          data_func, f_incorrect, max_runtime_second);
  std::vector<Tout> l_obj = unfold_opt(opt, "Enumeration unexpectedly failed");

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

template <typename T>
size_t INDEF_FORM_Invariant_NonDeg(MyMatrix<T> const &SymMat, size_t seed,
                                   [[maybe_unused]] std::ostream &os) {
  int n = SymMat.rows();
  T det = DeterminantMat<T>(SymMat);
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(SymMat);
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
  return hash_ret;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
LORENTZ_TestEquivalenceMatrices(MyMatrix<T> const &LorMat1,
                                MyMatrix<T> const &LorMat2, std::ostream &os) {
  size_t seed = 1678;
  size_t inv1 = INDEF_FORM_Invariant_NonDeg(LorMat1, seed, os);
  size_t inv2 = INDEF_FORM_Invariant_NonDeg(LorMat2, seed, os);
  if (inv1 != inv2) {
    return {};
  }
  int n = LorMat1.rows();
  int dimEXT = n + 1;
  using TintGroup = typename Tgroup::Tint;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<Tint>(dimEXT, os);
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
  std::optional<MyMatrix<Tint>> opt;
  auto f_incorrect = [&](Tobj const &x) -> bool {
    std::optional<MyMatrix<Tint>> opt_res =
        LORENTZ_TestEquivalence<T, Tint, Tgroup>(LorMat1, LorMat2, x.EXT, EXT2,
                                                 os);
    if (!opt_res) {
      MyMatrix<Tint> const &eBigMat = *opt_res;
      opt = eBigMat;
#ifdef DEBUG_LORENTZIAN_PERFECT
      MyMatrix<T> eBigMat_T = UniversalMatrixConversion<T, Tint>(eBigMat);
      MyMatrix<T> eProd = eBigMat_T * LorMat1 * eBigMat_T.transpose();
      if (eProd != LorMat2) {
        std::cerr << "It is not actually an equivalence\n";
        throw TerminalException{1};
      }
#endif
      return true;
      ;
    }
    return false;
  };
  int max_runtime_second = 0;
  (void)EnumerateAndStore_Serial<Tdata, decltype(f_incorrect)>(
      data_func, f_incorrect, max_runtime_second);
  return opt;
}

// clang-format off
#endif  // SRC_LORENTZIAN_PERFECT_FUND_H_
// clang-format on
