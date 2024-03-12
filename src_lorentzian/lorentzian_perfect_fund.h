// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_PERFECT_H_
#define SRC_LORENTZIAN_PERFECT_H_

#ifdef DEBUG
#define DEBUG_LORENTZIAN_PERFECT_FUND
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
    auto iele=[&]() -> std::vector<MyVector<Tint>> {
      if (TheOption == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
        return solver.FixedNormVectors(eSol, eSquareDist);
      } else {
        return solver.AtMostNormVectors(eSol, eSquareDist);
      }
    };
    for (auto & eSolA : iele()) {
      MyVector<Tint> eSolC = eBasSol + Ubasis.transpose() * eSolA;
      TotalListSol.emplace_back(std::move(eSolC));
    }
    if (OnlyShortest && TotalListSol.size()) {
      return ListSol;
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

template<typename T, typename Tint>
std::vector<MyVector<Tint>> LORENTZ_SearchInitialVector(MyMatrix<T> const& LorMat, MyVector<T> const& PosVect, int const& TheOption, std::ostream & os) {
  bool OnlyShortest = true;
  T MaxScal = 0;
  return LORENTZ_FindPositiveVectors(LorMat, PosVect, MaxScal, TheOption, OnlyShortest, os);
}


template<typename T>
std::vector<MyVector<T>> GetRationalIsotropyVectors(MyMatrix<T> const& LorMat22) {
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
  T a = LorMat22(0,0);
  T b = LorMat22(0,1);
  T c = LorMat22(1,1);
  T DeltaB = b*b - a*c;
  if (DeltaB <= 0) {
    std::cerr << "If the discriminant is negative, then we cannot have signature (1,1)";
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

template<typename T>
MyVector<T> GetReducedVector(MyVector<T> const& V) {
  int n = V.size() - 1;
  MyVector<T> Vret(n);
  for (int i=0; i<n; i++) {
    Vret(i) = V(i+1);
  }
  return Vret;
}

template<typename T>
MyVector<T> ConcatenateScalarVector(T const& scal, MyVector<T> const& V) {
  int n = V.size();
  MyVector<T> Vret(n + 1);
  Vret(0) = scal;
  for (int i=0; i<n; i++) {
    Vret(i+1) = V(i);
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
// -- Thus fVect is in the orthogonal of an isotropic vector and is of norm >= 0.
// -- By the geometry this gets us to fVect a multiple of eVect
// That scenario is not acceptable for finding perfect domain.
template<typename T>
std::optional<T> GetUpperBound(MyMatrix<T> const& LorMat, MyVector<T> const& eNSPbas, MyVector<T> const& eNSPdir) {
  int n = LorMat.rows();
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
  while(true) {
    MyVector<T> eV = eBas + iShift * eDir;
    MyVector<T> eVect = LorMatInv * eV;
    T eNorm = EvaluationQuadForm<T,T>(LorMat, eVect);
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
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
  auto f_norm=[&](T const& val) -> T {
    MyVector<T> eVect = eVectBas + fact * eVectDir;
    return EvaluationQuadForm<T,T>(LorMat, eVect);
  };
#endif
  std::vector<MyVector<T>> TheBasis{eVectBas, eVectDir};
  MyMatrix<T> TheBasis_mat = MatrixFromVectorFamily(TheBasis);
  MyMatrix<T> LorMat22 = TheBasis_mat * LorMat * TheBasis_mat.transpose();
  std::vector<MyVector<T>> ListIso = GetRationalIsotropyVectors(LorMat22);
  std::optional<T> UpperBound_isotropic;
  std::vector<T> ListUpperBound_Iso;
  for (auto & eIso : ListIso) {
    if (eIso(0) != 0) {
      T fact = eIso(1) / eIso(0);
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
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
  // Need to see if better upper bounds are possible, but this is a secondary question
  T BestUpper = VectorMin(ListUpperBound);
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
  if (f_norm(BestUpper) > 0) {
    std::cerr << "We should have eNorm <= 0\n";
    throw TerminalException{1};
  }
#endif
  return BestUpper;
}

template<typename T, typename Tint>
struct ResultFlipping {
  std::vector<MyVector<Tint>> ListTotal;
  MyVector<T> eNSPtest;
  MyVector<T> eVectTest;
  T MaxScal;
}

// Given a Critical set a vector eNSPbas and a direction eNSPdir
// we are looking for a lambda > 0 such that eNSP = eNSPbas + lambda * eNSPdir
// such that the list of vectors satisfying eNSP * v = 0 defines a bigger
// set than CritSet.
//
// eNSPbas must be of positive norm. eNSPdir must be of negative norm.
//
template<typename T, typename Tint>
ResultFlipping<T,Tint> LORENTZ_Kernel_Flipping(MyMatrix<T> const& LorMat, std::vector<MyVector<Tint>> const& CritSet, MyVector<T> const& eNSPbas, MyVector<T> const& eNSPdir, int option, [[maybe_unused]] std::ostream & os) {
  int n = LorMat.rows();
  MyMatrix<T> LorMatInv = Inverse(LorMat);
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
  std::vector<MyVector<T>> M1{eNSPbas, eNSPdir};
  MyMatrix<T> M2 = MatrixFromVectorFamily(M1);
  if (RankMat(M2) != 2) {
    std::cerr << "The vector eNSPbas and eNSPdir should be linearly independent\n";
    throw TerminalException{1};
  }
  for (auto & x : CritSet) {
    MyVector<Tint> xext = ConcatenateScalerVector(1, x);
    MyVector<T> xext_T = UniversalVectorConversion<T,Tint>(xext);
    if (xext_T.dot(eNSPbas) != 0) {
      std::cerr << "eNSPbas should have scalar product 0 with all entries in CritSet\n";
      throw TerminalException{1};
    }
    if (xext_T.dot(eNSPdir) != 0) {
      std::cerr << "eNSPdir should have scalar product 0 with all entries in CritSet\n";
      throw TerminalException{1};
    }
    if (option == LORENTZIAN_PERFECT_OPTION_ISOTROP) {
      T norm = EvaluationQuadForm<T,Tint>(LorMat, x);
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
  int n_iter = 0;
  MyVector<T> eVectBas = LorMatInv * GetReducedVector(eNSPbas);
  MyVector<T> eVectDir = LorMatInv * GetReducedVector(eNSPdir);
  T eNormBas = EvaluationQuadForm<T,T>(LorMat, eVectBas);
  T eNormDir = EvaluationQuadForm<T,T>(LorMat, eVectDir);
  int n_iter = 0;
  while(true) {
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
    os << "TheLowerBound=" << TheLowerBound << " TheUpperBound=" << TheUpperBound << " n_iter=" << n_iter << "\n";
#endif
    T TheMidVal = GetMidVal(TheLowerBound, TheUpperBound);
    MyVector<T> eNSPtest = eNSPbas + TheMidVal * eNSPdir;
    MyVector<T> eVectTest = LorMatInv * GetReducedVector(eNSPtest);
    T eNormTest = EvaluationQuadForm<T,T>(LorMat, eVectTest);
    T MaxScal = ScalarProductQuadForm(LorMat, CritSet[0], eVectTest);
    if (eNormTest <= 0 && MaxScal <= 0) {
      TheUpperBound = MidVal;
    } else {
      std::vector<MyVector<Tint>> ListTotal = LORENTZ_FindPositiveVectors(LorMat, eVectTest, MaxScal, TheOption, OnlyShortest);
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
      os << "|ListTotal|=" << ListTotal.size() << "\n";
      if (IsSubset(CritSet, ListTotal) && CritSet.size() > ListTotal.size()) {
        std::cerr << "Bug: if included, it should be equal\n";
        throw TerminalException{1};
      }
#endif
      if (Set(ListTotal) == Set(CritSet)) {
        TheLowerBound = MidVal;
      } else {
        if (IsSubset(ListTotal, CritSet)) {
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
          os << "EXIT 1 |ListTotal|=" << ListTotal.size() << " MaxScal=" << MaxScal << "\n";
          os << "  NSPtest=" << StringVectorGAP(eNSPtest) << "\n";
#endif
          return {ListTotal, eNSPtest, eVectTest, MaxScal};
        } else {
          break;
        }
      }
    }
    n_iter += 1;
  }
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
  os << "Going to the second scheme\n";
#endif
  while (true) {
    MyVector<Tint> eVect = ListTotal[0];
    MyVector<Tint> eVert = ConcatenateScalarVector(1, eVect);
    std::vector<MyVector<Tint>> ListIsoTest = CritSet;
    ListIsoTest.push_back(eVect);
    T aShift = - (eNSPbas.dot(eVert)) / (eNSPdir.dot(eVert));
    MyVector<Tint> eNSPtest = eNSPbas + aShift * eNSPdir;
    MyVector<T> eVectTest = LorMatInv * GetReducedVector(eNSPtest);
    T MaxScal = ScalarProductQuadForm(LorMat, CritSet[0], eVectTest);
    std::vector<MyVector<Tint>> ListTotal = LORENTZ_FindPositiveVectors(LorMat, eVectTest, MaxScal, TheOption, OnlyShortest);
    if (IsSubset(ListTotal, CritSet)) {
#ifdef DEBUG_LORENTZIAN_PERFECT_FUND
      os << "EXIT 2 |ListTotal|=" << ListTotal.size() << " MaxScal=" << MaxScal << "\n";
      os << "  NSPtest=" << StringVectorGAP(eNSPtest) << "\n";
#endif
      return {ListTotal, eNSPtest, eVectTest, MaxScal};
    }
  }

}








template<typename T, typename Tint>
MyVector<Tint> LORENTZ_GetShortPositiveVector(MyVector<T> const& LorMat) {
  int n = LorMat.rows();
  auto L1_norm=[&](

}


template<typename T, typename Tint>
MyMatrix<Tint> LORENTZ_GetOnePerfect(MyMatrix<Tint> const& LorMat) {
  

}




// clang-format off
#endif  // SRC_LORENTZIAN_PERFECT_H_
// clang-format on
