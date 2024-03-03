// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_PERFECT_H_
#define SRC_LORENTZIAN_PERFECT_H_

#ifdef DEBUG
# define DEBUG_LORENTZIAN_PERFECT_FUND
#endif


static const int LORENTZIAN_PERFECT_OPTION_ISOTROP = 23;
static const int LORENTZIAN_PERFECT_OPTION_TOTAL = 47;


template<typename T, typename Tint>
std::vector<MyVector<Tint>> LORENTZ_FindPositiveVectors(MyMatrix<T> const& LorMat, MyVector<T> const& eVect, T const& MaxScal, int const& TheOption, bool const& OnlyShortest, std::ostream & os) {
  int n = LorMat.rows();
  T eNorm = EvaluationQuadForm(LorMat, eVect);
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
  if (MaxScal <= 0) {
    std::cerr << "MaxScal=" << MaxScal << "\n";
    throw TerminalException{1};
  }
  if (TheOption != LORENTZIAN_PERFECT_OPTION_ISOTROP && TheOption != LORENTZIAN_PERFECT_OPTION_TOTAL) {
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
  MyVector<T> eVect_LorMat = LorMat * eVect;
  FractionVector<T> eRec = RemoveFractionVectorPlusCoeff(eVect_LorMat);
  MyVector<Tint> eVect_LorMat_tint = UniversalVectorConversion<Tint,T>(eRec.TheVect);
  MyMatrix<Tint> eVect_LorMat_tint_M(1,n);
  for (int u=0; u<n; u++) {
    eVect_LorMat_tint_M(0,u) = eVect_LorMat_tint(u);
  }
  MyMatrix<Tint> Ubasis = NullspaceIntMat(eVectM_LorMat_tint_M);
  MyMatrix<T> Ubasis_T = UniversalMatrixConversion<T,Tint>(Ubasis);
  MyMatrix<T> GramMat = - Ubasis_T * LorMat * Ubasis_T.transpose();
  CVPSolver<T,Tint> solver(GramMat, os);
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
  if (!IsPositiveDefiniteSymmetricMatrix(GramMat)) {
    std::cerr << "GramMat should be positive definite\n";
    throw TerminalException{1};
  }
#endif
  GCD_dot<Tint> TheRec_pre = ComputeGcdDot(eVect_LorMat_tint);
  GCD_dot<Tint> TheRec = PositivityNormalizeGcdDot(TheRec_pre);
  std::vector<MyVector<Tint>> TotalListSol;
  Tint eVal = 1;
  while(true) {
    MyVector<Tint> eBasSol = eVal * TheRec.V; // A solution of
    Tint alpha = eVal * TheRec.gcd / eNorm;
    MyVector<Tint> eTrans = eBasSol - alpha * eVect;
    std::optional<MyVector<Tint>> opt = SolutionMat(Ubasis, eTrans);
    MyVector<Tint> eSol = unfold_opt(opt);
    T eSquareDist = alpha * alpha * eNorm;
    std::vector<MyVector<Tint>> ListSolA = 
    
    if (OnlyShortest 
    
    eVal += 1;
  }
  return TotalListSol;
}


LORENTZ_SearchInitialVector

LORENTZ_GetShortPositiveVector


template<typename T, typename Tint>
MyMatrix<Tint> LORENTZ_GetOnePerfect(MyMatrix<Tint> const& LorMat) {
  

}




// clang-format off
#endif  // SRC_LORENTZIAN_PERFECT_H_
// clang-format on
