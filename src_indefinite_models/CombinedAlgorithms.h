// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
#define SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_

// clang-format off
#include "ApproximateModels.h"
#include "lorentzian_perfect.h"
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
  MyMatrix<T> GramMatRed;
  INDEF_FORM_GetVectorStructure(MyMatrix<T> const& _Qmat, MyVector<Tint> const& _v) : Qmat(_Qmat), v(_v) {
    eNorm = EvaluationQuadForm<T,Tint>(Qmat, v);
    MyVector<T> v_T = UniversalVectorConversion<T,Tint>(v);
    MyVector<T> eProd = Qmat * v_T;
    MyMatrix<T> NSP_T = NullspaceIntVect(eProd);
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
    if (eNorm != 0) {
      MyMatrix<T> eEndoRed_T = UniversalMatrixConversion<T,Tint>(eEndoRed);
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
      MyMatrix<Tint> Subspace1 = eEndoRed * NSP;
      MyMatrix<Tint> const& Subspace2 = NSP;
      MyMatrix<T> RetMat = LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(Qmat, Subspace1, Qmat, Subspace2);
      MyVector<T> vImg = RetMat.transpose() * v_T;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      if (vImg != v_T && vImg != -v) {
        std::cerr << "RetMat should map v to v or -v\n";
        throw TerminalException{1};
      }
#endif
      if (vImg == v) {
        return RetMat;
      } else {
        return -RetMat;
      }
    }
  }
  std::vector<MyMatrix<T>> MapOrthogonalSublatticeGroup(std::vector<MyMatrix<Tint>> const& GRPmatr) {
    std::vector<T> NewListGen;
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
  MyMatrix<T> Plane_T;
  MyMatrix<T> PlaneExpr_T;
  int dimSpace;
  int dim;
  MyMatrix<T> NSP_T;
  MyMatrix<T> GramMatRed;
  MyMatrix<T> QmatRed;
  MyMatrix<Tint> FullBasis;
  MyMatrix<Tint> FullBasisInv;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
  void check_generator(MyMatrix<Tint> const& eEndoRed, MyMatrix<T> const& RetMat) {
    MyMatrix<T> eEndoRed_T = UniversalMatrixConversion<T,Tint>(eEndoRed);
    MyMatrix<T> PlaneImg = Plane * RetMat;
    MyMatrix<T> TransRed(dim, dim);
    for (int u=0; u<dim; u++) {
      MyVector<T> eV = GetMatrixRow(PlaneImg, u);
      std::optional<MyVector<T>> opt = SolutionMat(Plane, eV);
      MyVector<T> fV = unfold_opt(opt, "Get the vector fV");
      AssignMatrixRow(TransRed, u, fV);
    }
    T eDet = DeterminantMat(TransRed);
    if (T_abs(eDet) != 1) {
      std::cerr << "TransRed should have absolute determinant 1\n";
      throw TerminalException{1};
    }
    if (!TestEqualitySpace(PlaneImg, Plane)) {
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
    MyMatrix<T> eProd = Plane * Qmat;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> WitnessIsotropy = Plane_T * Qmat * Plane_T.transpose();
    if (!IsZeroMatrix(WitnessIsotropy)) {
      std::cerr << "The matrix Plane does not define a totally isotropic space\n";
      throw TerminalException{1};
    }
#endif
    NSP_T = NullspaceIntTrMat(eProd);
    GramMatRed = NSP * Qmat * MSP.transpose();
    Plane_T = UniversalMatrixConversion<T,Tint>(Plane);

    int dimNSP = eRec.NSP.rows();
    PlaneExpr_T = MyMatrix<T>(dim, dimNSP);
    for (int u=0; u<k; u++) {
      MyMatrix<T> eV = GetMatrixRow(ePlane_T, u);
      std::optional<MyVector<T>> opt = SolutionIntMat(NSP_T, eV);
      MyVector<T> fV = unfold_opt(opt, "getting fV");
      AssignMatrixRow(PlaneExpr_T, u, fV);
    }
    int the_dim = eRec.NSP_T.rows();
    MyMatrix<T> TheCompl = SubspaceCompletionInt(PlaneExpr_T, the_dim);
    int dimCompl = TheCompl.rows();
    FullBasis = Concatenation(TheCompl, PlaneExpr_T);
    FullBasisInv = Inverse(FullBasis);
    QmatRed = TheCompl * eRec.GramMatRed * TheCompl.transpose();
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
      MyMatrix<T> Subspace2 = eEndoRed_T * NSP;
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
}

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



template<typename T, typename Tint, typename Tgroup>
struct IndefiniteCombinedAlgo {
private:
  std::ostream& os;
  DatabaseResultEquiStab<MyMatrix<T>, MyMatrix<Tint>> database;
private:
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup_Kernel(MyMatrix<T> const& Qmat) {
    AttackScheme<T> eBlock = INDEF_FORM_GetAttackScheme(Qmat);
    if (eBlock.h == 0) {
      return INDEF_FORM_AutomorphismGroup_PosNeg(Qmat);
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
    ApproxModel<T,Tint> approx = INDEF_FORM_GetApproximateModel<T,Tint,Tgroup>(Qmat);
    FirstNorm<T,Tint> first_norm = GetFirstNorm(approx);
    T const& X = first_norm.X;
    MyVector<Tint> const& v1 = first_norm.eVect;
    for (auto & eGen : approx.GetApproximateGroup()) {
      f_insert(eGen);
    }
    for (auto & eGen : approx.INDEF_FORM_StabilizerVector(Qmat, v1)) {
      f_insert(eGen);
    }
    for (auto & v2 : approx.GetCoveringOrbitRepresentatives(X)) {
      std::optional<MyVector<Tint>> opt = INDEF_FORM_EquivalenceVector(Qmat, Qmat, v1, v2);
      if (opt) {
        f_insert(*opt);
      }
    }
    return ListGenerators;
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_TestEquivalence_Kernel(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2) {
    AttackScheme<T> eBlock1 = INDEF_FORM_GetAttackScheme(Qmat1);
    AttackScheme<T> eBlock2 = INDEF_FORM_GetAttackScheme(Qmat2);
    if (Block1.h == 0) {
      return INDEF_FORM_TestEquivalence_PosNeg(Qmat1, Qmat2);
    }
    if (Block1.h == 1) {
      return LORENTZ_TestEquivalenceMatrices<T, Tint, Tgroup>(LorMat1, LorMat2, os);
    }
    ApproximateModel<T,Tint> approx1 = INDEF_FORM_GetApproximateModel<T,Tint,Tgroup>(Qmat1);
    FirstNorm<T,Tint> first_norm1 = GetFirstNorm(approx1);
    T const& X = first_norm1.X;
    MyVector<Tint> const& v1 = first_norm1.eVect;
    ApproximateModel<T,Tint> approx2 = INDEF_FORM_GetApproximateModel<T,Tint,Tgroup>(Qmat2);
    std::vector<MyVector<Tint>> ListCand2 = Approx2.GetCoveringOrbitRepresentatives(X);
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
        std::vector<MyMatrix<Tint>> LGen = INDEF_FORM_AutomorphismGroup_Kernel(QmatRed);
        database.insert_stab(eMat, ListGen);
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
        database.insert_equi(eMat1, eMat2, optB);
        return optB;
      }
    };
    std::optional<MyMatrix<Tint>> opt = get_equi();
    if (opt) {
      MyMatrix<Tint> const& eEquiv = *opt;
      MyMatrix<Tint> NewEquiv = Inverse(RecRed2.B) * eEquiv * RecRed1.B;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      MyMatrix<T> NewGen_T = UniversalMatrixConversion<T,Tint>(NewGen);
      MyMatrix<T> eProd = NewGen_T * Qmat1 * NewGen_T.transpose();
      if (eProd != Qmat2) {
        std::cerr << "The matrix is not an equivalence for Qmat1 / Qmat2\n";
        throw TerminalException{1};
      }
#endif
      return NewEquiv;
    } else {
      return {}:
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
    std::vector<MyMatrix<Tint>> GRPfull = ExtendIsometryGroup_Triangular(GRPred, eRec.dimCompl, the_dim);
    std::vector<MyMatrix<Tint>> ListGenTot;
    for (auto & eGen : GRPfull) {
      MyMatrix<T> eGenB = eRec.FullBasisInv * eGen * eRec.FullBasis;
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
      MyMatrix<T> eProd = eGenB * eRec.GramMatRed * eGenB.transpose();
      if (eProd != eRec.GramMatRed) {
        std::cerr << "eGenB should preserve eRec.GramMatRed\n";
        throw TerminalException{1};
      }
      std::vector<MyMatrix<Tint>> ListSpaces = f_get_list_spaces(eRec.PlaneExpr);
      for (auto & eSpace : ListSpaces) {
        MyMatrix<Tint> eSpaceImg = eSpace * eGenB;
        if (!TestEqualitySpace(eSpaceImg, eSpace)) {
          std::cerr << "The space eSpace is not correctly preserved\n";
          throw TerminalException{1};
        }
      }
#endif
      ListGenTot.push_back(eGenB);
    }
    return ListGenTot;
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
    MyMatrix<Tint> TheEquivTest = IdentityMat<Tint>(the_dim);
    int p = eRec1.QmatRed1.rows();
    for (int i=0; i<p; i++) {
      for (int j=0; j<p; j++) {
        TheEquivTest(i,j) = test(i,j);
      }
    }
    MyMatrix<Tint> TheEquiv = eRec2.FullBasisInv * TheEquivTest * FullBasis1;
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
      MyMatrix<Tint> eSpace1_img = ListSpaces1[i] * TheEquivInv;
      MyMatrix<Tint> const& eSpace2 = ListSpaces2[i];
      if (!TestEqualitySpace(eSpace1_img, eSpace2)) {
        std::cerr << "The space are not correctly mapped\n";
        throw TerminalException{1};
      }
    }
#endif
    return TheEquiv;
  }
  //
  // The function using template arguments
  //
  template<typename Fstab>
  size_t INDEF_FORM_Invariant_IsotropicKstuff_Kernel(MyMatrix<T> const& Qmat, MyMatrix<Tint> const& Plane, Fstab f_stab) {
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec(Qmat, Plane);
    std::vector<MyMatrix<T>> GRP1 = f_stab(eRec);
    std::vector<MyMatrix<T>> GRP2 = eRec.MapOrthogonalSublatticeGroup(GRP1);
    size_t eInvRed = INDEF_FORM_Invariant_IsotropicKplane_Raw(Qmat, Plane);
    size_t GRP_inv = GetRationalInvariant(GRP2);
    return eInvRed + GRP_inv;
  }
  template<typename Fequiv, typename Fstab>
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2, Fequiv f_equiv, Fstab f_stab) {
    if (INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat1, Plane1, f_stab) != INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat2, Plane2, f_stab)) {
      return {};
    }
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec1(Qmat1, Plane1);
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec2(Qmat2, Plane2);
    std::optional<MyVector<Tint>> opt = f_equiv(eRec1, eRec2);
    if (!opt) {
      return {};
    }
    MyMatrix<Tint> const& test = *opt;
    MyMatrix<Tint> Subspace1 = Inverse(test) * eRec2.NSP;
    MyMatrix<Tint> Subspace2 = eRec1.NSP;
    MyMatrix<T> Subspace1_T = UniversalMatrixConversion<T,Tint>(Subspace1);
    MyMatrix<T> Subspace2_T = UniversalMatrixConversion<T,Tint>(Subspace2);
    LORENTZ_ExtendOrthogonalIsotropicIsomorphism<T,> TheRec(Qmat1, Subspace1, Qmat2, Subspace2);
    MyMatrix<T> EquivRat = TheRec.get_one_transformation();
    MyMatrix<T> EquivRatInv = Inverse(EquivRat);
#ifdef DEBUG_INDEFINITE_COMBINED_ALGORITHMS
    MyMatrix<T> Plane1_img = Plane1 * EquivRatInv;
    if (!TestEqualitySpace(Plane1_img, Plane2)) {
      std::cerr << "Plane1 and Plane2 should be mapped\n";
      throw TerminalException{1};
    }
#endif
    auto fTest=[&](MyMatrix<Tint> const& eTest) -> bool {
      MyMatrix<T> eProd = Plane1 * Inverse(eTest);
      return TestEqualitySpace(eProd, Plane2);
    };
    return Kernel_Equivalence_Qmat(Qmat1, Qmat2, EquivRat, eRec1, fTest, f_stab);
  }
  template<typename Fstab>
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel(MyMatrix<T> comnst& Qmat, MyMatrix<Tint> const& Plane, Fstab f_stab) {
    int n Qmat.rows();
    if (RankMat(Qmat) != n) {
      std::cerr << "Right now INDEF_FORM_StabilizerVector requires Qmat to be full dimensional\n";
      throw TerminalException{1};
    }
    INDEF_FORM_GetRec_IsotropicKplane<T,Tint> eRec(Qmat, Plane);
    std::vector<MyMatrix<Tint>> GRP1 = f_stab(eRec);
    std::vector<MyMatrix<T>> GRP2 = eRec.MapOrthogonalSublatticeGroup(GRP1);
    return MatrixIntegral_Stabilizer(n, GRP2);
  }

public:
  // Now the specific implementations
  InvariantIsotropic<T,Tint> INDEF_FORM_Invariant_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
    return INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat, Plane, f_stab_plane);
  }
  InvariantIsotropic<T,Tint> INDEF_FORM_Invariant_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
    return INDEF_FORM_Invariant_IsotropicKstuff_Kernel(Qmat, Plane, f_stab_flag);
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKflag(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2) {
    return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(Qmat1, Qmat2, Plane1, Plane2, f_equiv_flag, f_stab_flag);
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKplane(MyMatrix<T> const& Qmat1, MyMatrix<T> const& Qmat2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2) {
    return INDEF_FORM_Equivalence_IsotropicKstuff_Kernel(Qmat1, Qmat2, Plane1, Plane2, f_equiv_plane, f_stab_plane);
  }



  
  std::vector<MyVector<Tint>> INDEF_FORM_GetOrbitRepresentative(MyMatrix<T> const& Qmat) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_StabilizerVector(MyMatrix<T> const& Qmat) {
    if (RankMat(Qmat) != Qmat.rows()) {
      std::cerr << "Right now the matrix Qmat should be full dimensional\n";
      throw TerminalExpression{1};
    }
    INDEF_FORM_GetVectorStructure<T,Tint> eRec(Qmat, v);
    std::vector<MyMatrix<Tint>> GRP1 = INDEF_FORM_AutomorphismGroup(eRec.GramMatRed);
    std::vector<MyMatrix<T>> GRP2_T = eRec.MapOrthogonalSublatticeGroup(GRP1);
    std::vector<MyMatrix<T>> ListMat_T = MatrixIntegral_Stabilizer(GRP2_T);
    std::vector<MyMatrix<Tint>> ListMat;
    for (auto & eMat_T : ListMat_T) {
      MyMatrix<Tint> eMat = UniversalMatrixConversion<Tint,T>(eMat_T);
      ListMat.push_back(eMat);
    }
    return ListMat;
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_EquivalenceVector(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyVector<Tint> const& v1, MyVector<Tint> const& v2) {
    if (INDEF_FORM_InvariantVector(Q1, v1) != INDEF_FORM_InvariantVector(Q2, v2)) {
      return {};
    }
    
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_StabilizerVector(MyMatrix<T> const& Q, MyVector<Tint> const& v) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup(MyMatrix<T> const& Q) {
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_TestEquivalence(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2) {
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKplane(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_RightCosets_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_GetOrbit_IsotropicKplane(MyMatrix<T> const& Q) {
  }
  std::optional<MyMatrix<Tint>> INDEF_FORM_Equivalence_IsotropicKflag(MyMatrix<T> const& Q1, MyMatrix<T> const& Q2, MyMatrix<Tint> const& Plane1, MyMatrix<Tint> const& Plane2) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_Stabilizer_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_RightCosets_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  std::vector<MyMatrix<Tint>> INDEF_FORM_GetOrbit_IsotropicKflag(MyMatrix<T> const& Q) {
  }

};

// clang-format off
#endif  // SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
// clang-format on
