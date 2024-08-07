// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
#define SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_

// clang-format off
#include "ApproximateModels.h"
#include "lorentzian_perfect.h"
#include "lorentzian_linalg.h"
#include "EquiStabMemoization.h"
#include "MatrixGroup.h"
// clang-format on

// This is a reimplementation of the GAP code and implements the algorithm for indefinite forms.
// 
//
// The result of the paper were published in
// Mathieu Dutour Sikirić, Klaus Hulek, Moduli of polarized Enriques surfaces -- computational aspects,
// Journal of the London Mathematical Society (2023)
// preprint at https://arxiv.org/abs/2302.01679

#ifdef DEBUG
#define DEBUG_INDEFINITE_COMBINED_ALGORITHMS
#endif

static const int INDEFINITE_FORM_PLANE = 32;
static const int INDEFINITE_FORM_FLAG = 92;


template<typename T, typename Tint>
struct INDEF_FORM_GetVectorStructure {
public:
  T eNorm;
  MyMatrix<T> Qmat;
  MyVector<Tint> v;
  MyVector<T> v_T;
  MyMatrix<T> Pmat;
  MyMatrix<T> PmatInv;
  MyMatrix<Tint> NSP;
  MyMatrix<T> NSP_T;
  MyMatrix<T> GramMatRed;
  INDEF_FORM_GetVectorStructure(MyMatrix<T> const& _Qmat, MyVector<Tint> const& _v) : Qmat(_Qmat), v(_v) {
    int n = Qmat.rows();
    eNorm = EvaluationQuadForm<T,Tint>(Qmat, v);
    v_T = UniversalVectorConversion<T,Tint>(v);
    MyVector<T> eProd = Qmat * v_T;
    NSP_T = NullspaceIntVect(eProd);
    NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
    GramMatRed = NSP_T * Qmat * NSP_T.transpose();
    if (eNorm != 0) {
      Pmat = ZeroMatrix<T>(n,n);
      for (int i=0; i<n-1; i++) {
        for (int j=0; j<n; j++) {
          Pmat(i, j) = NSP(i, j);
        }
      }
      for (int i=0; i<n; i++) {
        Pmat(n-1,i) = v_T(i);
      }
      PmatInv = Inverse(Pmat);
    }
  }
  MyMatrix<T> MapOrthogonalSublatticeEndomorphism(MyMatrix<Tint> const& eEndoRed) {
    MyMatrix<T> eEndoRed_T = UniversalMatrixConversion<T,Tint>(eEndoRed);
    if (eNorm != 0) {
      MyMatrix<T> TheBigMat = ExpandMatrix(eEndoRed_T);
      MyMatrix<T> RetMat = PmatInv * TheBigMat * Pmat;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      MyVector<T> prodV = RetMat.transpose() * v_T;
      if (prodV != v_T) {
        std::cerr << "RetMat is not preserving the vector v\n";
        throw TerminalException{1};
      }
#endif
      return RetMat;
    } else {
      MyMatrix<T> Subspace1 = eEndoRed_T * NSP_T;
      MyMatrix<T> const& Subspace2 = NSP_T;
      std::optional<MyMatrix<T>> opt = LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(Qmat, Subspace1, Qmat, Subspace2);
      MyMatrix<T> RetMat = unfold_opt(opt, "opt should be something because NSP.rows = RankMat(NSP)");
      MyVector<T> vImg = RetMat.transpose() * v_T;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      if (vImg != v_T && vImg != -v_T) {
        std::cerr << "RetMat should map v to v or -v\n";
        throw TerminalException{1};
      }
#endif
      if (vImg == v_T) {
        return RetMat;
      } else {
        return -RetMat;
      }
    }
  }
  std::vector<MyMatrix<T>> MapOrthogonalSublatticeGroup(std::vector<MyMatrix<Tint>> const& GRPmatr) {
    std::vector<MyMatrix<T>> NewListGen;
    for (auto & eGen : GRPmatr) {
      NewListGen.push_back(MapOrthogonalSublatticeEndomorphism(eGen));
    }
    return NewListGen;
  }
};


template<typename T, typename Tint>
struct INDEF_FORM_GetRec_IsotropicKplane {
public:
  MyMatrix<T> Qmat;
  MyMatrix<Tint> Plane;
  MyMatrix<T> Plane_T;
  MyMatrix<Tint> PlaneExpr;
  MyMatrix<T> PlaneExpr_T;
  int dimSpace;
  int dim;
  MyMatrix<Tint> NSP;
  MyMatrix<T> NSP_T;
  MyMatrix<T> GramMatRed;
  MyMatrix<T> QmatRed;
  MyMatrix<T> FullBasis_T;
  MyMatrix<Tint> FullBasis;
  MyMatrix<Tint> FullBasisInv;
  int the_dim;
  int dimCompl;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
  void check_generator(MyMatrix<Tint> const& eEndoRed, MyMatrix<T> const& RetMat) {
    MyMatrix<T> eEndoRed_T = UniversalMatrixConversion<T,Tint>(eEndoRed);
    MyMatrix<T> PlaneImg = Plane_T * RetMat;
    MyMatrix<T> TransRed(dim, dim);
    for (int u=0; u<dim; u++) {
      MyVector<T> eV = GetMatrixRow(PlaneImg, u);
      std::optional<MyVector<T>> opt = SolutionMat(Plane_T, eV);
      MyVector<T> fV = unfold_opt(opt, "Get the vector fV");
      AssignMatrixRow(TransRed, u, fV);
    }
    T eDet = DeterminantMat(TransRed);
    if (T_abs(eDet) != 1) {
      std::cerr << "TransRed should have absolute determinant 1\n";
      throw TerminalException{1};
    }
    if (!TestEqualitySpaces(PlaneImg, Plane_T)) {
      std::cerr << "Plane should be invariant (isotropic case)\n";
      throw TerminalException{1};
    }
    int dimNSP = NSP_T.rows();
    MyMatrix<T> RetMat_red(dimNSP, dimNSP);
    for (int u=0; u<dimNSP; u++) {
      MyVector<T> eV = GetMatrixRow(NSP_T, u);
      MyVector<T> fV = RetMat.transpose() * eV;
      std::optional<MyVector<T>> opt = SolutionMat(NSP_T, fV);
      MyVector<T> gV = unfold_opt(opt, "getting gV");
      AssignMatrixRow(RetMat_red, u, gV);
    }
    if (RetMat_red != eEndoRed_T) {
      std::cerr << "RetMat_red restricted to the nullspace is not the original eEndoRed\n";
      throw TerminalException{1};
    }
  }
#endif
  INDEF_FORM_GetRec_IsotropicKplane(MyMatrix<T> const& _Qmat, MyMatrix<Tint> const& _Plane) : Qmat(_Qmat), Plane(_Plane), dimSpace(Qmat.rows()), dim(Plane.rows()) {
    Plane_T = UniversalMatrixConversion<T,Tint>(Plane);
    MyMatrix<T> eProd = Plane_T * Qmat;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> WitnessIsotropy = Plane_T * Qmat * Plane_T.transpose();
    if (!IsZeroMatrix(WitnessIsotropy)) {
      std::cerr << "The matrix Plane does not define a totally isotropic space\n";
      throw TerminalException{1};
    }
#endif
    NSP_T = NullspaceIntTrMat(eProd);
    NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
    GramMatRed = NSP_T * Qmat * NSP_T.transpose();
    the_dim = NSP_T.rows();
    PlaneExpr_T = MyMatrix<T>(dim, the_dim);
    for (int u=0; u<dim; u++) {
      MyVector<T> eV = GetMatrixRow(Plane_T, u);
      std::optional<MyVector<T>> opt = SolutionIntMat(NSP_T, eV);
      MyVector<T> fV = unfold_opt(opt, "getting fV");
      AssignMatrixRow(PlaneExpr_T, u, fV);
    }
    PlaneExpr = UniversalMatrixConversion<Tint,T>(PlaneExpr_T);
    MyMatrix<T> TheCompl = SubspaceCompletionInt(PlaneExpr_T, the_dim);
    dimCompl = TheCompl.rows();
    FullBasis_T = Concatenate(TheCompl, PlaneExpr_T);
    FullBasis = UniversalMatrixConversion<Tint,T>(FullBasis_T);
    FullBasisInv = Inverse(FullBasis);
    QmatRed = TheCompl * GramMatRed * TheCompl.transpose();
  }
  // We map automorphism group of a sublattice to the full group.
  // The difficulty os that the rational kernel is of strictly positive Q-rank
  // if dim(Plane) > 1.
  std::vector<MyMatrix<T>> MapOrthogonalSublatticeGroup(std::vector<MyMatrix<Tint>> const& GRPmatr) {
    MyMatrix<T> const& Subspace1 = NSP_T;
    std::vector<MyMatrix<T>> ListGenTotal;
    std::vector<T> ListD{1};
    for (auto & eEndoRed : GRPmatr) {
      MyMatrix<T> eEndoRed_T = UniversalMatrixConversion<T,Tint>(eEndoRed);
      MyMatrix<T> Subspace2 = eEndoRed_T * NSP_T;
      LORENTZ_ExtendOrthogonalIsotropicIsomorphism<T> TheRec(Qmat, Subspace1, Qmat, Subspace2);
      MyMatrix<T> RetMat = TheRec.get_one_transformation();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      check_generator(eEndoRed, RetMat);
#endif
      ListGenTotal.push_back(RetMat);
      T eDen = GetDenominatorMatrix(RetMat);
      ListD.push_back(eDen);
    }
    T TheDen = LCMlist(ListD);
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism<T> TheRec(Qmat, Subspace1, Qmat, Subspace1);
    std::vector<MyMatrix<T>> TheKer = TheRec.get_kernel_generating_set(TheDen);
    MyMatrix<Tint> eEndoRed = IdentityMat<Tint>(NSP.rows());
    for (auto & RetMat : TheKer) {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      check_generator(eEndoRed, RetMat);
#endif
      ListGenTotal.push_back(RetMat);
    }
    return ListGenTotal;
  }
};

template<typename T>
std::vector<MyMatrix<T>> GetAutomorphismOfFlag(int const& n) {
  std::vector<MyMatrix<T>> LGen;
  for (int i=0; i<n; i++) {
    MyMatrix<T> TheMat = IdentityMat<T>(n);
    TheMat(i,i) = -1;
    LGen.push_back(TheMat);
  }
  for (int i=0; i<n; i++) {
    for (int j=0; j<i; j++) {
      MyMatrix<T> TheMat = IdentityMat<T>(n);
      TheMat(i, j) = 1;
      LGen.push_back(TheMat);
    }
  }
  return LGen;
}

template<typename T>
std::vector<MyMatrix<T>> ExtendIsometryGroup_Triangular(std::vector<MyMatrix<T>> const& GRPmatr, int const& p, int const& n) {
  std::vector<MyMatrix<T>> ListGens;
  for (auto & eGen : GRPmatr) {
    MyMatrix<T> NewMat = IdentityMat<T>(n);
    for (int i=0; i<p; i++) {
      for (int j=0; j<p; j++) {
        NewMat(i, j) = eGen(i, j);
      }
    }
    ListGens.push_back(NewMat);
  }
  std::vector<MyMatrix<T>> SubNPgroup = GetAutomorphismOfFlag<T>(n-p);
  for (auto & eGen : SubNPgroup) {
    MyMatrix<T> NewMat = IdentityMat<T>(n);
    for (int i=0; i<n-p; i++) {
      for (int j=0; j<n-p; j++) {
        NewMat(i+p, j+p) = eGen(i, j);
      }
    }
    ListGens.push_back(NewMat);
  }
  for (int i=0; i<p; i++) {
    for (int idx=p; idx<n; idx++) {
      MyMatrix<T> NewMat = IdentityMat<T>(n);
      NewMat(i, idx) = 1;
      ListGens.push_back(NewMat);
    }
  }
  return ListGens;
}

/*
  As obtained by GeneratorsOfGroup(GeneralLinearGroup(6,Integers))
  in GAP.
 */
template<typename T>
std::vector<MyMatrix<T>> GeneralLinearGroup(int const& n) {
  std::vector<MyMatrix<T>> ListGens;
  if (n > 1) {
    MyMatrix<T> mat1 = ZeroMatrix<T>(n,n);
    for (int i=0; i<n; i++) {
      int iNext = 0;
      if (i<n-1) {
        iNext = i+1;
      }
      mat1(i,iNext) = 1;
    }
    ListGens.push_back(mat1);
  }
  //
  if (n > 2) {
    MyMatrix<T> mat2 = ZeroMatrix<T>(n,n);
    mat2(1,0) = 1;
    mat2(0,1) = 1;
    for (int i=2; i<n; i++) {
      mat2(i,i) = 1;
    }
    ListGens.push_back(mat2);
  }
  //
  MyMatrix<T> mat3 = IdentityMat<T>(n);
  mat3(0,0) = -1;
  ListGens.push_back(mat3);
  //
  if (n > 1) {
    MyMatrix<T> mat4 = IdentityMat<T>(n);
    mat4(0,1) = 1;
    ListGens.push_back(mat4);
  }
  //
  return ListGens;
}


/*
  Isometry group defined on a p dimensional space for a quadratic form Qp.
  We extend the quadratic form to dimension n with
  Qn = | Qp 0 |
       | 0  0 |
  If the original matrices satisfy Pp Qp Pp^T = Qp
  the the extended matrices must satisfy Pn Qn Pn^T = Qn
  So, in block formulation
  Pn = | A B |
       | C D |
  and so
  Pn Qn Pn^T = | A B |     | Qp 0 |     | A^T C^T |
               | C D |  x  | 0  0 |  x  | B^T D^T |
             = | A Qp 0 |     | A^T C^T |
               | C Qp 0 |  x  | B^T D^T |
             = | A Qp A^T  A Qp C^T |
               | C Qp A^T  C Qp C^T |
             = | Qp 0 |
               | 0  0 |
  And so we get A Qp A^T = Qp , C Qp A^T = 0 , C Qp C^T = 0
  The equation A Qp A^T forces A to be an isometry of Qp.
  The equation C Qp A^T forces C to be 0 from which the rest follows.
*/
template<typename T>
std::vector<MyMatrix<T>> ExtendIsometryGroup(std::vector<MyMatrix<T>> const& GRPmatr, int const& p, int const& n) {
  std::vector<MyMatrix<T>> ListGens;
  for (auto & eGen : GRPmatr) {
    MyMatrix<T> NewMat = IdentityMat<T>(n);
    for (int i=0; i<p; i++) {
      for (int j=0; j<p; j++) {
        NewMat(i,j) = eGen(i,j);
      }
    }
    ListGens.push_back(eGen);
  }
  if (n > p) {
    for (auto & eGen : GeneralLinearGroup<T>(n-p)) {
      MyMatrix<T> NewMat = IdentityMat<T>(n);
      for (int i=0; i<n-p; i++) {
        for (int j=0; j<n-p; j++) {
          NewMat(i+p, j+p) = eGen(i,j);
        }
      }
      ListGens.push_back(NewMat);
    }
    for (int i=0; i<p; i++) {
      MyMatrix<T> NewMat = IdentityMat<T>(n);
      NewMat(i, p) = 1;
      ListGens.push_back(NewMat);
    }
  }
  return ListGens;
}


template<typename T, typename Tint, typename Tgroup>
std::vector<MyVector<Tint>> INDEF_FORM_GetOrbitRepresentative_PosNeg(MyMatrix<T> const& Q, T const& X, std::ostream & os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Tgroup::Tidx;
  auto get_orbit_representatives=[&](MyMatrix<T> const& Qin) -> std::vector<MyVector<Tint>> {
    MyMatrix<Tint> SHV = EnumerateVectorsFixedNorm<T, Tint>(Qin, X, os);
    size_t len = SHV.rows();
    std::vector<MyVector<Tint>> ListVect;
    std::unordered_map<MyVector<Tint>, size_t> MapVect;
    for (int i=0; i<SHV.rows(); i++) {
      MyVector<Tint> eV = GetMatrixRow(SHV, i);
      MapVect[eV] = i+1;
      ListVect.push_back(eV);
    }
    std::vector<MyMatrix<Tint>> LGen = ArithmeticAutomorphismGroup<T,Tint>(Qin, os);
    std::vector<Telt> ListPerm;
    for (auto & eGen : LGen) {
      std::vector<Tidx> eList(len);
      for (size_t i=0; i<len; i++) {
        MyVector<Tint> eV = ListVect[i];
        MyVector<Tint> eVimg = eGen.transpose() * eV;
        size_t pos = MapVect[eVimg];
        if (pos == 0) {
          std::cerr << "The image of the vector is not where it should be\n";
          throw TerminalException{1};
        }
        eList[i] = pos - 1;
      }
      Telt ePerm(eList);
      ListPerm.push_back(ePerm);
    }
    Tgroup eG(ListPerm, len);
    std::vector<size_t> LPos = DecomposeOrbitPoint_FullRepr(eG);
    std::vector<MyVector<Tint>> ListSol;
    for (auto & ePos : LPos) {
      ListSol.push_back(ListVect[ePos]);
    }
    return ListSol;
  };
  DiagSymMat<T> DSM = DiagonalizeNonDegenerateSymmetricMatrix(Q);
  if (DSM.nbPlus == 0 && DSM.nbZero == 0) {
    MyMatrix<T> Qneg = -Q;
    return get_orbit_representatives(Qneg);
  }
  if (DSM.nbMinus == 0 && DSM.nbZero == 0) {
    return get_orbit_representatives(Q);
  }
  std::cerr << "Failed to find a matching entry for PosNeg\n";
  throw TerminalException{1};
}

template<typename T, typename Tint>
std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup_PosNeg(MyMatrix<T> const& Q, std::ostream & os) {
  DiagSymMat<T> DSM = DiagonalizeNonDegenerateSymmetricMatrix(Q);
  if (DSM.nbPlus == 0 && DSM.nbZero == 0) {
    MyMatrix<T> Qneg = -Q;
    return ArithmeticAutomorphismGroup<T,Tint>(Qneg, os);
  }
  if (DSM.nbMinus == 0 && DSM.nbZero == 0) {
    return ArithmeticAutomorphismGroup<T,Tint>(Q, os);
  }
  std::cerr << "Failed to find a matching entry for PosNeg\n";
  throw TerminalException{1};
}

template<typename T, typename Tint>
std::optional<MyMatrix<Tint>> INDEF_FORM_TestEquivalence_PosNeg(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, std::ostream & os) {
  DiagSymMat<T> DSM1 = DiagonalizeNonDegenerateSymmetricMatrix(Q1);
  DiagSymMat<T> DSM2 = DiagonalizeNonDegenerateSymmetricMatrix(Q2);
  if (DSM1.nbPlus != DSM2.nbPlus) {
    return {};
  }
  if (DSM1.nbMinus != DSM2.nbMinus) {
    return {};
  }
  if (DSM1.nbPlus == 0 && DSM1.nbZero == 0) {
    MyMatrix<T> Qneg1 = -Q1;
    MyMatrix<T> Qneg2 = -Q2;
    return ArithmeticEquivalence<T,Tint>(Qneg1, Qneg2, os);
  }
  if (DSM1.nbMinus == 0 && DSM1.nbZero == 0) {
    return ArithmeticEquivalence<T,Tint>(Q1, Q2, os);
  }
  std::cerr << "Failed to find a matching entry for PosNeg\n";
  throw TerminalException{1};
}





template<typename T, typename Tint, typename Tgroup>
struct IndefiniteCombinedAlgo {
private:
  std::ostream& os;
  DatabaseResultEquiStab<MyMatrix<T>, MyMatrix<Tint>> database;
private:
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup_Kernel(MyMatrix<T> const& Qmat) {
    AttackScheme<T> eBlock = INDEF_FORM_GetAttackScheme(Qmat);
    if (eBlock.h == 0) {
      return INDEF_FORM_AutomorphismGroup_PosNeg<T,Tint>(Qmat, os);
    }
    if (eBlock.h == 1) {
      LORENTZ_GetGeneratorsAutom<T, Tint, Tgroup>(Qmat, os);
    }
    std::vector<MyMatrix<Tint>> ListGenerators;
    auto f_insert=[&](MyMatrix<Tint> const& eGen) -> void {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      MyMatrix<T> eGen_T = UniversalMatrixConversion<T,Tint>(eGen);
      MyMatrix<T> prod = eGen_T * Qmat * eGen_T.transpose();
      if (prod != Qmat) {
        std::cerr << "eGen is not preserving Qmat\n";
        throw TerminalException{1};
      }
#endif
      ListGenerators.push_back(eGen);
    };
    ApproximateModel<T,Tint> approx = INDEF_FORM_GetApproximateModel<T,Tint,Tgroup>(Qmat, os);
    FirstNorm<T,Tint> first_norm = GetFirstNorm(approx);
    T const& X = first_norm.X;
    MyVector<Tint> const& v1 = first_norm.eVect;
    for (auto & eGen : approx.GetApproximateGroup()) {
      f_insert(eGen);
    }
    for (auto & eGen : INDEF_FORM_StabilizerVector(Qmat, v1)) {
      f_insert(eGen);
    }
    for (auto & v2 : approx.GetCoveringOrbitRepresentatives(X)) {
      std::optional<MyMatrix<Tint>> opt = INDEF_FORM_EquivalenceVector(Qmat, Qmat, v1, v2);
      if (opt) {
        f_insert(*opt);
      }
    }
    return ListGenerators;
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_TestEquivalence_Kernel(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2) {
    AttackScheme<T> eBlock1 = INDEF_FORM_GetAttackScheme(Qmat1);
    AttackScheme<T> eBlock2 = INDEF_FORM_GetAttackScheme(Qmat2);
    if (eBlock1.h == 0) {
      return INDEF_FORM_TestEquivalence_PosNeg<T,Tint>(Qmat1, Qmat2, os);
    }
    if (eBlock1.h == 1) {
      return LORENTZ_TestEquivalenceMatrices<T, Tint, Tgroup>(eBlock1.mat, eBlock2.mat, os);
    }
    ApproximateModel<T,Tint> approx1 = INDEF_FORM_GetApproximateModel<T,Tint,Tgroup>(Qmat1, os);
    FirstNorm<T,Tint> first_norm1 = GetFirstNorm(approx1);
    T const& X = first_norm1.X;
    MyVector<Tint> const& v1 = first_norm1.eVect;
    ApproximateModel<T,Tint> approx2 = INDEF_FORM_GetApproximateModel<T,Tint,Tgroup>(Qmat2, os);
    std::vector<MyVector<Tint>> ListCand2 = approx2.GetCoveringOrbitRepresentatives(X);
    for (auto & v2 : ListCand2) {
      std::optional<MyMatrix<Tint>> opt = INDEF_FORM_EquivalenceVector(Qmat1, Qmat2, v1, v2);
      if (opt) {
        return opt;
      }
    }
    return {};
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup_FullDim(MyMatrix<T> const& Qmat) {
    ResultReduction<T, Tint> ResRed =
      ComputeReductionIndefinitePermSign<T, Tint>(Qmat, os);
    MyMatrix<T> const& QmatRed = ResRed.Mred;
    MyMatrix<Tint> const& B = ResRed.B;
    MyMatrix<Tint> Binv = Inverse(B);
    auto get_stab=[&]() -> std::vector<MyMatrix<Tint>> {
      std::optional<std::vector<MyMatrix<Tint>>> opt = database.attempt_stabilizer(QmatRed);
      if (opt) {
        return *opt;
      } else {
        std::vector<MyMatrix<Tint>> ListGen = INDEF_FORM_AutomorphismGroup_Kernel(QmatRed);
        database.insert_stab(QmatRed, ListGen);
        return ListGen;
      }
    };
    std::vector<MyMatrix<Tint>> LGenFinal;
    for (auto & eGen : get_stab()) {
      MyMatrix<Tint> NewGen = Binv * eGen * B;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      MyMatrix<T> NewGen_T = UniversalMatrixConversion<T,Tint>(NewGen);
      MyMatrix<T> eProd = NewGen_T * Qmat * NewGen_T.transpose();
      if (eProd != Qmat) {
        std::cerr << "The matrix is not an equivalence for Qmat\n";
        throw TerminalException{1};
      }
#endif
      LGenFinal.push_back(NewGen);
    }
    return LGenFinal;
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_TestEquivalence_FullDim(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2) {
    ResultReduction<T, Tint> ResRed1 =
      ComputeReductionIndefinitePermSign<T, Tint>(Qmat1, os);
    ResultReduction<T, Tint> ResRed2 =
      ComputeReductionIndefinitePermSign<T, Tint>(Qmat2, os);
    MyMatrix<T> const& QmatRed1 = ResRed1.Mred;
    MyMatrix<T> const& QmatRed2 = ResRed2.Mred;
    auto get_equi=[&]() -> std::optional<MyMatrix<Tint>> {
      std::optional<std::optional<MyMatrix<Tint>>> opt = database.attempt_equiv(QmatRed1, QmatRed2);
      if (opt) {
        return *opt;
      } else {
        std::optional<MyMatrix<Tint>> optB = INDEF_FORM_TestEquivalence_Kernel(QmatRed1, QmatRed2);
        database.insert_equi(QmatRed1, QmatRed2, optB);
        return optB;
      }
    };
    std::optional<MyMatrix<Tint>> opt = get_equi();
    if (opt) {
      MyMatrix<Tint> const& eEquiv = *opt;
      MyMatrix<Tint> NewEquiv = Inverse(ResRed2.B) * eEquiv * ResRed1.B;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      MyMatrix<T> NewEquiv_T = UniversalMatrixConversion<T,Tint>(NewEquiv);
      MyMatrix<T> eProd = NewEquiv_T * Qmat1 * NewEquiv_T.transpose();
      if (eProd != Qmat2) {
        std::cerr << "The matrix is not an equivalence for Qmat1 / Qmat2\n";
        throw TerminalException{1};
      }
#endif
      return NewEquiv;
    } else {
      return {};
    }
  }
  //
  // Specific functions f_stab_*, f_equiv_*
  //
  std::vector<MyMatrix<Tint>> f_stab_plane(INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec) {
    return INDEF_FORM_AutomorphismGroup(eRec.GramMatRed);
  }
  std::vector<MyMatrix<Tint>> f_get_list_spaces(MyMatrix<Tint> const& ListVect) {
    int n_vect=ListVect.rows();
    int dim = ListVect.cols();
    std::vector<MyMatrix<Tint>> ListSpaces;
    for (int len=1; len<=n_vect; len++) {
      MyMatrix<Tint> eSpace(len, dim);
      for (int u=0; u<len; u++) {
        for (int iCol=0; iCol<dim; iCol++) {
          eSpace(u, iCol) = ListVect(u, iCol);
        }
      }
      ListSpaces.push_back(eSpace);
    }
    return ListSpaces;
  }
  std::vector<MyMatrix<Tint>> f_stab_flag(INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec) {
    std::vector<MyMatrix<Tint>> GRPred = INDEF_FORM_AutomorphismGroup(eRec.QmatRed);
    std::vector<MyMatrix<Tint>> GRPfull = ExtendIsometryGroup_Triangular(GRPred, eRec.dimCompl, eRec.the_dim);
    std::vector<MyMatrix<Tint>> ListGenTot;
    for (auto & eGen : GRPfull) {
      MyMatrix<Tint> eGenB = eRec.FullBasisInv * eGen * eRec.FullBasis;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      MyMatrix<T> eGenB_T = UniversalMatrixConversion<T,Tint>(eGenB);
      MyMatrix<T> eProd = eGenB_T * eRec.GramMatRed * eGenB_T.transpose();
      if (eProd != eRec.GramMatRed) {
        std::cerr << "eGenB should preserve eRec.GramMatRed\n";
        throw TerminalException{1};
      }
      std::vector<MyMatrix<Tint>> ListSpaces = f_get_list_spaces(eRec.PlaneExpr);
      for (auto & eSpace : ListSpaces) {
        MyMatrix<Tint> eSpaceImg = eSpace * eGenB;
        if (!TestEqualitySpaces(eSpaceImg, eSpace)) {
          std::cerr << "The space eSpace is not correctly preserved\n";
          throw TerminalException{1};
        }
      }
#endif
      ListGenTot.push_back(eGenB);
    }
    return ListGenTot;
  }
  std::vector<MyMatrix<Tint>> f_stab(INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec, int choice) {
    if (choice == INDEFINITE_FORM_PLANE) {
      return f_stab_plane(eRec);
    }
    if (choice == INDEFINITE_FORM_FLAG) {
      return f_stab_flag(eRec);
    }
    std::cerr << "Failed to have a valid input choice in f_stab\n";
    throw TerminalException{1};
  }
  std::optional<MyMatrix<Tint>> f_equiv_plane(INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec1, INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec2) {
    return INDEF_FORM_TestEquivalence(eRec1.GramMatRed, eRec2.GramMatRed);
  }
  std::optional<MyMatrix<Tint>> f_equiv_flag(INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec1, INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec2) {
    std::optional<MyMatrix<Tint>> opt = INDEF_FORM_TestEquivalence(eRec1.QmatRed, eRec2.QmatRed);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const& test = *opt;
    MyMatrix<Tint> TheEquivTest = IdentityMat<Tint>(eRec2.the_dim);
    int p = eRec1.QmatRed.rows();
    for (int i=0; i<p; i++) {
      for (int j=0; j<p; j++) {
        TheEquivTest(i,j) = test(i,j);
      }
    }
    MyMatrix<Tint> TheEquiv = eRec2.FullBasisInv * TheEquivTest * eRec1.FullBasis;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> TheEquiv_T = UniversalMatrixConversion<T,Tint>(TheEquiv);
    MyMatrix<T> eProd = TheEquiv_T * eRec1.GramMatRed * TheEquiv_T.transpose();
    if (eProd != eRec2.GramMatRed) {
      std::cerr << "This is not an equivalence\n";
    }
    std::vector<MyMatrix<Tint>> ListSpaces1 = f_get_list_spaces(eRec1.PlaneExpr);
    std::vector<MyMatrix<Tint>> ListSpaces2 = f_get_list_spaces(eRec2.PlaneExpr);
    MyMatrix<Tint> TheEquivInv = Inverse(TheEquiv);
    for (size_t iSpace=0; iSpace<ListSpaces1.size(); iSpace++) {
      MyMatrix<Tint> eSpace1_img = ListSpaces1[iSpace] * TheEquivInv;
      MyMatrix<Tint> const& eSpace2 = ListSpaces2[iSpace];
      if (!TestEqualitySpaces(eSpace1_img, eSpace2)) {
        std::cerr << "The space are not correctly mapped\n";
        throw TerminalException{1};
      }
    }
#endif
    return TheEquiv;
  }
  std::optional<MyMatrix<Tint>> f_equiv(INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec1, INDEF_FORM_GetRec_IsotropicKplane<T,Tint> const& eRec2, int const& choice) {
    if (choice == INDEFINITE_FORM_PLANE) {
      return f_equiv_plane(eRec1, eRec2);
    }
    if (choice == INDEFINITE_FORM_FLAG) {
      return f_equiv_flag(eRec1, eRec2);
    }
    std::cerr << "Failed to have a valid input choice in f_equiv\n";
    throw TerminalException{1};
  }
  //
  // The function using choice integer input
  //
  size_t INDEF_FORM_Invariant_IsotropicKstuff_Kernel(MyMatrix<T> const& Qmat, MyMatrix<Tint> const& Plane, int const& choice) {
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec(Qmat, Plane);
    std::vector<MyMatrix<Tint>> GRP1 = f_stab(eRec, choice);
    std::vector<MyMatrix<T>> GRP2 = eRec.MapOrthogonalSublatticeGroup(GRP1);
    size_t eInvRed = INDEF_FORM_Invariant_IsotropicKplane_Raw(Qmat, Plane);
    size_t GRP_inv = GetRationalInvariant(GRP2);
    return eInvRed + GRP_inv;
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2, int const& choice) {
    if (INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat1, Plane1, choice) != INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat2, Plane2, choice)) {
      return {};
    }
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec1(Qmat1, Plane1);
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec2(Qmat2, Plane2);
    std::optional<MyMatrix<Tint>> opt = f_equiv(eRec1, eRec2, choice);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const& test = *opt;
    MyMatrix<Tint> Subspace1 = Inverse(test) * eRec2.NSP;
    MyMatrix<Tint> Subspace2 = eRec1.NSP;
    MyMatrix<T> Subspace1_T = UniversalMatrixConversion<T,Tint>(Subspace1);
    MyMatrix<T> Subspace2_T = UniversalMatrixConversion<T,Tint>(Subspace2);
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism<T> TheRec(Qmat1, Subspace1_T, Qmat2, Subspace2_T);
    MyMatrix<T> EquivRat = TheRec.get_one_transformation();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> eProd = EquivRat * Qmat1 * EquivRat.transpose();
    if (eProd != Qmat2) {
      std::cerr << "The matrix EquivRat is not mapping Qmat1 to Qmat2\n";
      throw TerminalException{1};
    }
#endif
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> EquivRatInv = Inverse(EquivRat);
    MyMatrix<T> Plane1_T = UniversalMatrixConversion<T,Tint>(Plane1);
    MyMatrix<T> Plane2_T = UniversalMatrixConversion<T,Tint>(Plane2);
    MyMatrix<T> Plane1_img = Plane1_T * EquivRatInv;
    if (!TestEqualitySpaces(Plane1_img, Plane2_T)) {
      std::cerr << "Plane1 and Plane2 should be mapped\n";
      throw TerminalException{1};
    }
#endif
    if (IsIntegralMatrix(EquivRat)) {
      // This is an Ansatz. If the matrix is integral, no more work to be done.
      MyMatrix<Tint> EquivRat_tint = UniversalMatrixConversion<Tint,T>(EquivRat);
      return EquivRat_tint;
    }
    std::vector<MyMatrix<Tint>> GRP1_A = f_stab(eRec1, choice);
    std::vector<MyMatrix<T>> GRP1_B = eRec1.MapOrthogonalSublatticeGroup(GRP1_A);
    // Find a g1 in GRP1_B such that EquivRat * g1 in GL(n,Z)
    std::optional<MyMatrix<Tint>> optB = MatrixIntegral_Equivalence_Bis_General<T,Tint,Tgroup>(GRP1_B, EquivRat, os);
    if (!optB) {
      return {};
    }
    MyMatrix<Tint> const& TheRet = *optB;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> TheRet_T = UniversalMatrixConversion<T,Tint>(TheRet);
    MyMatrix<T> eProdB = TheRet_T * Qmat1 * TheRet_T.transpose();
    if (eProdB != Qmat2) {
      std::cerr << "The redution matrix did not work as we expected it. Please debug\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> TheRetInv = Inverse(TheRet);
    MyMatrix<Tint> Plane1_imgB = Plane1 * TheRetInv;
    if (!TestEqualitySpaces(Plane1_imgB, Plane2)) {
      std::cerr << "Plane1 and Plane2 should be mapped\n";
      throw TerminalException{1};
    }
#endif
    return TheRet;
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(MyMatrix<T> const& Qmat, MyMatrix<Tint> const& Plane, int const& choice) {
    int n = Qmat.rows();
    if (RankMat(Qmat) != n) {
      std::cerr << "Right now INDEF_FORM_StabilizerVector requires Qmat to be full dimensional\n";
      throw TerminalException{1};
    }
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec(Qmat, Plane);
    std::vector<MyMatrix<Tint>> GRP1 = f_stab(eRec, choice);
    std::vector<MyMatrix<T>> GRP2 = eRec.MapOrthogonalSublatticeGroup(GRP1);
    return MatrixIntegral_Stabilizer_General<T,Tint,Tgroup>(n, GRP2, os);
  }
  std::vector<MyMatrix<T>> INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(MyMatrix<T> const& Qmat, MyMatrix<Tint> const& ePlane, int const& choice) {
    // We have two groups:
    // -- The group stabilizing ePlane
    // -- The group stabilizing ePlane^{perp} and its mapping to the full group.
    // The group stabilizing ePlane^{perp} can be injected
    int n = Qmat.rows();
    if (RankMat(Qmat) != n) {
      std::cerr << "Right now INDEF_FORM_StabilizerVector requires Qmat to be full dimensional\n";
      throw TerminalException{1};
    }
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec(Qmat, ePlane);
    std::vector<MyMatrix<Tint>> GRP1 = f_stab(eRec, choice);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    for (auto & eGen : GRP1) {
      MyMatrix<T> eGen_T = UniversalMatrixConversion<T,Tint>(eGen);
      MyMatrix<T> eProd = eGen_T * eRec.GramMatRed * eGen_T.transpose();
      if (eProd != eRec.GramMatRed) {
        std::cerr << "GRP1 does not preserve eRec.GramMatRed\n";
        throw TerminalException{1};
      }
    }
#endif
    std::vector<MyMatrix<T>> GRP2 = eRec.MapOrthogonalSublatticeGroup(GRP1);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    for (auto & eGen : GRP2) {
      MyMatrix<T> eProd = eGen * Qmat * eGen.transpose();
      if (eProd != Qmat) {
        std::cerr << "GRP2 does not preserve Qmat\n";
        throw TerminalException{1};
      }
      MyMatrix<T> ePlane_T = UniversalMatrixConversion<T,Tint>(ePlane);
      MyMatrix<T> ePlaneImg = ePlane_T * eGen;
      int k = ePlane.rows();
      for (int u=0; u<k; u++) {
        MyVector<T> eLine = GetMatrixRow(ePlaneImg, u);
        std::optional<MyVector<T>> opt = SolutionIntMat(ePlane_T, eLine);
        if (!opt) {
          std::cerr << "ePlane should be left invariant by eGen\n";
          throw TerminalException{1};
        }
      }
    }
#endif
    std::vector<MyMatrix<T>> ListRightCoset = MatrixIntegral_RightCosets_General<T,Tgroup>(n, GRP2, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    for (auto & eCos : ListRightCoset) {
      MyMatrix<T> eProd = eCos * Qmat * eCos.transpose();
      if (eProd != Qmat) {
        std::cerr << "eCos is not preserving Qmat\n";
        throw TerminalException{1};
      }
    }
#endif
    return ListRightCoset;
  }
  size_t fiso_inv(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane, int const& choice) {
    if (choice == INDEFINITE_FORM_PLANE) {
      return INDEF_FORM_Invariant_IsotropicKplane(Q, Plane);
    }
    if (choice == INDEFINITE_FORM_FLAG) {
      return INDEF_FORM_Invariant_IsotropicKflag(Q, Plane);
    }
    std::cerr << "Failed to have a valid input choice in f_inv\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> fiso_coset(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane, int const& choice) {
    return INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(Q, Plane, choice);
  }
  std::optional<MyMatrix<Tint>> fiso_equiv(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2, int const& choice) {
    return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(Q1, Q2, Plane1, Plane2, choice);
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel(MyMatrix<T> const& Qmat, int k, int const& choice) {
    T eNorm(0);
    std::vector<MyMatrix<Tint>> ListOrbit;
    for (auto & eVect : INDEF_FORM_GetOrbitRepresentative(Qmat, eNorm)) {
      MyMatrix<Tint> ePlane = MatrixFromVector(eVect);
      ListOrbit.push_back(ePlane);
    }
    for (int iK=2; iK<=k; iK++) {
      struct PlaneInv {
        MyMatrix<Tint> ePlane;
        size_t eInv;
      };
      auto get_planeinv=[&](MyMatrix<Tint> const& Plane) -> PlaneInv {
        size_t eInv = fiso_inv(Qmat, Plane, choice);
        return {Plane, eInv};
      };
      std::vector<PlaneInv> ListRecReprKplane;
      auto fInsert=[&](PlaneInv const& fRecReprKplane) -> void {
        for (auto & eRecReprKplane : ListRecReprKplane) {
          if (eRecReprKplane.eInv == fRecReprKplane.eInv) {
            std::optional<MyMatrix<Tint>> opt = fiso_equiv(Qmat, Qmat, eRecReprKplane.ePlane, fRecReprKplane.ePlane, choice);
            if (opt) {
              return;
            }
          }
        }
        ListRecReprKplane.push_back(fRecReprKplane);
      };
      auto SpanRepresentatives=[&](MyMatrix<Tint> const& ePlane) {
        // Some possible improvement. Use the double cosets
        // The double coset consists in splitting an orbit x G as
        // y1 H \cup ..... yN H
        // Or in other words G = \cup_i Stab(x) y_i H
        //
        // In that case
        // G = group of rational transformation preserving L = S^{perp}
        //     and acting integrally on it.
        // Stab(x) = Integral transformations preserving L and some vector x in L.
        // H = group of integral transformations preserving the big lattice
        //     and preserving L.
        //
        // If G = H the only one entry to treat.
        // What that mean is that when we go down the chain of subgroup by breaking down
        // the orbit split. Can this happen in the same way over all the orbits?
        // The full of the lattice group is G(Qmat) and the full orbit is x G(Qmat)
        // We want to write the code as
        MyMatrix<T> ePlane_T = UniversalMatrixConversion<T,Tint>(ePlane);
        MyMatrix<T> ePlaneQ = ePlane_T * Qmat;
        MyMatrix<T> NSP_T = NullspaceIntTrMat(ePlaneQ);
        MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
        int dimNSP = NSP.rows();
        MyMatrix<Tint> ePlane_expr(k, dimNSP);
        for (int u=0; u<k; u++) {
          MyVector<Tint> eV = GetMatrixRow(ePlane, u);
          std::optional<MyVector<Tint>> opt = SolutionIntMat(NSP, eV);
          MyVector<Tint> eSol = unfold_opt(opt, "eSol should not be fail");
          AssignMatrixRow(ePlane_expr, u, eSol);
        }
        MyMatrix<Tint> ComplBasisInNSP = SubspaceCompletionInt(ePlane_expr, dimNSP);
        MyMatrix<Tint> NSP_sub = ComplBasisInNSP * NSP;
        MyMatrix<T> NSP_sub_T = UniversalMatrixConversion<T,Tint>(NSP_sub);
        MyMatrix<T> QmatRed = NSP_sub_T * Qmat * NSP_sub_T.transpose();
        std::vector<MyVector<Tint>> ListOrbitF = INDEF_FORM_GetOrbitRepresentative(QmatRed, eNorm);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: |ListOrbitF|=" <<  ListOrbitF.size() << "\n";
#endif
        std::vector<MyMatrix<T>> ListRightCosets = fiso_coset(Qmat, ePlane, choice);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
        os << "COMB: |ListRightCosets|=" << ListRightCosets.size() << "\n";
#endif
        for (auto & eVect : ListOrbitF) {
          MyVector<Tint> eVectB = NSP_sub.transpose() * eVect;
          MyVector<T> eVectB_T = UniversalVectorConversion<T,Tint>(eVectB);
          for (auto & eCos : ListRightCosets) {
            MyVector<T> eVectC_T = eCos.transpose() * eVectB_T;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
            if (EvaluationQuadForm(Qmat, eVectC_T) != 0) {
              std::cerr << "eVectC is not isotropic\n";
              throw TerminalException{1};
            }
            if (!IsIntegralVector(eVectC_T)) {
              std::cerr << "eVectC_T should be integral\n";
              throw TerminalException{1};
            }
#endif
            MyVector<Tint> eVectC = UniversalVectorConversion<Tint,T>(eVectC_T);
            MyMatrix<Tint> ePlaneB = ConcatenateMatVec(ePlane, eVectC);
            PlaneInv eRecReprKplane = get_planeinv(ePlaneB);
            fInsert(eRecReprKplane);
          }
        }
      };
      for (auto & eRepr : ListOrbit) {
        SpanRepresentatives(eRepr);
      }
      ListOrbit.clear();
      for (auto & eRec : ListRecReprKplane) {
        ListOrbit.push_back(eRec.ePlane);
      }
    }
    return ListOrbit;
  }
public:
  IndefiniteCombinedAlgo(std::ostream& _os) : os(_os) {
  }
  // Now the specific implementations
  size_t INDEF_FORM_Invariant_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
    return INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Q, Plane, INDEFINITE_FORM_PLANE);
  }
  size_t INDEF_FORM_Invariant_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
    return INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Q, Plane, INDEFINITE_FORM_FLAG);
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKplane(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2) {
    return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(Qmat1, Qmat2, Plane1, Plane2, INDEFINITE_FORM_PLANE);
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKflag(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2) {
    return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(Qmat1, Qmat2, Plane1, Plane2, INDEFINITE_FORM_FLAG);
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
    return INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(Q, Plane, INDEFINITE_FORM_PLANE);
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
    return INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(Q, Plane, INDEFINITE_FORM_FLAG);
  }
  std::vector<MyMatrix<T>> INDEF_FORM_RightCosets_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
    return INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(Q, Plane, INDEFINITE_FORM_PLANE);
  }
  std::vector<MyMatrix<T>> INDEF_FORM_RightCosets_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
    return INDEF_FORM_RightCosets_IsotropicKstuff_Kernel(Q, Plane, INDEFINITE_FORM_FLAG);
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_GetOrbit_IsotropicKplane(MyMatrix<T> const& Q, int k) {
    return INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel(Q, k, INDEFINITE_FORM_PLANE);
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_GetOrbit_IsotropicKflag(MyMatrix<T> const& Q, int k) {
    return INDEF_FORM_GetOrbit_IsotropicKstuff_Kernel(Q, k, INDEFINITE_FORM_FLAG);
  }
  std::vector<MyVector<Tint>> INDEF_FORM_GetOrbitRepresentative(MyMatrix<T> const& Q, T const& X) {
    AttackScheme<T> eBlock = INDEF_FORM_GetAttackScheme(Q);
    if (eBlock.h == 0) {
      return INDEF_FORM_GetOrbitRepresentative_PosNeg<T,Tint,Tgroup>(Q, X, os);
    }
    ApproximateModel<T,Tint> approx = INDEF_FORM_GetApproximateModel<T,Tint,Tgroup>(Q, os);
    std::vector<MyVector<Tint>> ListRepr;
    auto f_insert=[&](MyVector<Tint> fRepr) -> void {
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      T fNorm = EvaluationQuadForm<T,Tint>(Q, fRepr);
      if (fNorm != X) {
        std::cerr << "The norm is inconsistent\n";
        throw TerminalException{1};
      }
#endif
      for (auto & eRepr : ListRepr) {
        std::optional<MyMatrix<Tint>> opt = INDEF_FORM_EquivalenceVector(Q, Q, eRepr, fRepr);
        if (opt) {
          return;
        }
      }
      ListRepr.push_back(fRepr);
    };
    std::vector<MyVector<Tint>> ListCand = approx.GetCoveringOrbitRepresentatives(X);
    for (auto & eCand : ListCand) {
      f_insert(eCand);
    }
    return ListRepr;
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_StabilizerVector(MyMatrix<T> const& Qmat, MyVector<Tint> const& v) {
    if (RankMat(Qmat) != Qmat.rows()) {
      std::cerr << "Right now the matrix Qmat should be full dimensional\n";
      throw TerminalException{1};
    }
    int n = Qmat.rows();
    INDEF_FORM_GetVectorStructure<T,Tint> eRec(Qmat, v);
    std::vector<MyMatrix<Tint>> GRP1 = INDEF_FORM_AutomorphismGroup(eRec.GramMatRed);
    std::vector<MyMatrix<T>> GRP2_T = eRec.MapOrthogonalSublatticeGroup(GRP1);
    std::vector<MyMatrix<Tint>> ListMat = MatrixIntegral_Stabilizer_General<T,Tint,Tgroup>(n, GRP2_T, os);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    for (auto & eMat : ListMat) {
      MyMatrix<T> eMat_T = UniversalMatrixConversion<T,Tint>(eMat);
      MyMatrix<T> eProd = eMat_T * Qmat * eMat_T.transpose();
      if (eProd != Qmat) {
        std::cerr << "The matrix eMat does not preserves Qmat\n";
        throw TerminalException{1};
      }
      MyVector<Tint> v_eMat = eMat.transpose() * v;
      if (v_eMat != v) {
        std::cerr << "The matrix eMat does not preserves v\n";
        throw TerminalException{1};
      }
    }
#endif
    return ListMat;
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_EquivalenceVector(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyVector<Tint> const& v1, MyVector<Tint> const& v2) {
    if (INDEF_FORM_InvariantVector(Q1, v1) != INDEF_FORM_InvariantVector(Q2, v2)) {
      return {};
    }
    if (RankMat(Q1) != Q1.rows()) {
      std::cerr << "We need Q1 to be a square matrix\n";
      throw TerminalException{1};
    }
    T eNorm = EvaluationQuadForm<T,Tint>(Q1, v1);
    INDEF_FORM_GetVectorStructure<T,Tint> eRec1(Q1, v1);
    INDEF_FORM_GetVectorStructure<T,Tint> eRec2(Q2, v2);
    std::optional<MyMatrix<Tint>> opt = INDEF_FORM_TestEquivalence(eRec1.GramMatRed, eRec2.GramMatRed);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const& test = *opt;
    MyMatrix<T> test_T = UniversalMatrixConversion<T,Tint>(test);
    auto get_equiv_rat=[&]() -> MyMatrix<T> {
      if (eNorm != 0) {
        return eRec2.PmatInv * ExpandMatrix(test_T) * eRec1.Pmat;
      } else {
        MyMatrix<T> Subspace1 = Inverse(test_T) * eRec2.NSP_T;
        MyMatrix<T> Subspace2 = eRec1.NSP_T;
        std::optional<MyMatrix<T>> opt = LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(Q1, Subspace1, Q2, Subspace2);
        MyMatrix<T> EquivRat = unfold_opt(opt, "failed to find EquivRat");
        MyMatrix<T> EquivRatInv = Inverse(EquivRat);
        MyVector<T> v1_inv_equivrat = EquivRatInv.transpose() * eRec1.v_T;
        if (v1_inv_equivrat == eRec2.v_T) {
          return EquivRat;
        } else {
          return -EquivRat;
        }
      }
    };
    MyMatrix<T> EquivRat = get_equiv_rat();
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> EquivRatInv = Inverse(EquivRat);
    MyVector<T> v1_T = UniversalVectorConversion<T,Tint>(v1);
    MyVector<T> v2_T = UniversalVectorConversion<T,Tint>(v2);
    MyVector<T> v_EquivRatInv = EquivRatInv.transpose() * v1_T;
    if (v_EquivRatInv != v2_T) {
      std::cerr << "The vector v1 is not mapped correctly\n";
      throw TerminalException{1};
    }
#endif
    std::vector<MyMatrix<Tint>> GRP1_A = INDEF_FORM_AutomorphismGroup(eRec1.GramMatRed);
    std::vector<MyMatrix<T>> GRP1_B = eRec1.MapOrthogonalSublatticeGroup(GRP1_A);
    // Find a g1 in GRP1_B such that EquivRat * g1 in GL(n,Z)
    std::optional<MyMatrix<Tint>> optB = MatrixIntegral_Equivalence_Bis_General<T,Tint,Tgroup>(GRP1_B, EquivRat, os);
    if (!optB) {
      return {};
    }
    MyMatrix<Tint> const& TheRet = *optB;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> TheRet_T = UniversalMatrixConversion<T,Tint>(TheRet);
    MyMatrix<T> eProd = TheRet_T * Q1 * TheRet_T.transpose();
    if (eProd != Q2) {
      std::cerr << "The redution matrix did not work as we expected it. Please debug\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> TheRetInv = Inverse(TheRet);
    MyVector<Tint> v_TheRetInv = TheRetInv.transpose() * v1;
    if (v_TheRetInv != v2) {
      std::cerr << "The vector v1 is not mapped correctly\n";
      throw TerminalException{1};
    }
#endif
    return TheRet;
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup(MyMatrix<T> const& Q) {
    int n = Q.rows();
    MyMatrix<T> NSP_T = NullspaceIntMat(Q);
    MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
    MyMatrix<Tint> TheCompl = SubspaceCompletionInt(NSP, n);
    MyMatrix<Tint> FullBasis = Concatenate(TheCompl, NSP);
    MyMatrix<Tint> FullBasisInv = Inverse(FullBasis);
    int p = TheCompl.rows();
    MyMatrix<T> TheCompl_T = UniversalMatrixConversion<T,Tint>(TheCompl);
    MyMatrix<T> QmatRed = TheCompl_T * Q * TheCompl_T.transpose();
    std::vector<MyMatrix<Tint>> GRPred = INDEF_FORM_AutomorphismGroup_FullDim(QmatRed);
    std::vector<MyMatrix<Tint>> GRPfull = ExtendIsometryGroup(GRPred, p, n);
    std::vector<MyMatrix<Tint>> ListGenTot;
    for (auto & eGen : GRPfull) {
      MyMatrix<Tint> eGenB = FullBasisInv * eGen * FullBasis;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      MyMatrix<T> eGenB_T = UniversalMatrixConversion<T,Tint>(eGenB);
      MyMatrix<T> eProd = eGenB_T * Q * eGenB_T.transpose();
      if (eProd != Q) {
        std::cerr << "eGenB should preserve Qmat\n";
        throw TerminalException{1};
      }
#endif
      ListGenTot.push_back(eGenB);
    }
    return ListGenTot;
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_TestEquivalence(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2) {
    int n = Q1.rows();
    MyMatrix<T> NSP1_T = NullspaceIntMat(Q1);
    MyMatrix<Tint> NSP1 = UniversalMatrixConversion<Tint,T>(NSP1_T);
    MyMatrix<Tint> TheCompl1 = SubspaceCompletionInt(NSP1, n);
    MyMatrix<T> TheCompl1_T = UniversalMatrixConversion<T,Tint>(TheCompl1);
    MyMatrix<Tint> FullBasis1 = Concatenate(TheCompl1, NSP1);
    MyMatrix<T> QmatRed1 = TheCompl1_T * Q1 * TheCompl1_T.transpose();
    MyMatrix<T> NSP2_T = NullspaceIntMat(Q2);
    MyMatrix<Tint> NSP2 = UniversalMatrixConversion<Tint,T>(NSP2_T);
    MyMatrix<Tint> TheCompl2 = SubspaceCompletionInt(NSP2, n);
    MyMatrix<T> TheCompl2_T = UniversalMatrixConversion<T,Tint>(TheCompl2);
    MyMatrix<Tint> FullBasis2 = Concatenate(TheCompl2, NSP2);
    MyMatrix<T> QmatRed2 = TheCompl2_T * Q2 * TheCompl2_T.transpose();
    int p = TheCompl1.rows();
    std::optional<MyMatrix<Tint>> opt = INDEF_FORM_TestEquivalence_FullDim(QmatRed1, QmatRed2);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const& test = *opt;
    MyMatrix<Tint> TheEquivTest = IdentityMat<Tint>(n);
    for (int i=0; i<p; i++) {
      for (int j=0; j<p; j++) {
        TheEquivTest(i,j) = test(i,j);
      }
    }
    MyMatrix<Tint> TheEquiv = Inverse(FullBasis2) * TheEquivTest * FullBasis1;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> TheEquiv_T = UniversalMatrixConversion<T,Tint>(TheEquiv);
    MyMatrix<T> eProd = TheEquiv_T * Q1 * TheEquiv_T.transpose();
    if (eProd != Q2) {
      std::cerr << "TheEquiv is not mapping Q1 to Q2\n";
      throw TerminalException{1};
    }
#endif
    return TheEquiv;
  }
};

// clang-format off
#endif  // SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
// clang-format on
