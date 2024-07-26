// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
#define SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_

// clang-format off
#include "ApproximateModels.h"
#include "lorentzian_perfect.h"
#include "MatrixGroup.h"
// clang-format on

// This is a reimplementation of the GAP code and implements the 
//
// The result of the paper were published in
// Mathieu Dutour Sikirić, Klaus Hulek, Moduli of polarized Enriques surfaces -- computational aspects,
// Journal of the London Mathematical Society (2023)
// preprint at https://arxiv.org/abs/2302.01679

#ifdef DEBUG
#define DEBUG_INDEFINITE_COMBINED_ALGORITHMS
#endif


// Used for the memoization process for isomorphism
template<typename T, typename Tint>
struct ResultIsomorphism {
  MyMatrix<T> Q1;
  MyMatrix<T> Q2;
  std::optional<MyMatrix<Tint>> result;
};

template<typename T, typename Tint>
struct ResultStabilizer {
  MyMatrix<T> Q;
  std::vector<MyMatrix<Tint>> ListGens;
};

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
  int dimSpace;
  int dim;
  MyMatrix<T> NSP_T;
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
  std::vector<ResultIsomorphism<T,Tint>> ListResultIsomorphism;
  std::vector<ResultStabilizer<T,Tint>> ListResultStabilizer;

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
  std::vector<MyMatrix<Tint>> INDEF_FORM_AutomorphismGroup_Memoize(MyMatrix<T> const& Qmat) {
  }

  template<typename Fstab>
  size_t INDEF_FORM_Invariant_IsotropicKstuff_Kernel(MyMatrix<T> const& Qmat, MyMatrix<Tint> const& Plane, Fstab f_stab) {
  eRec:=INDEF_FORM_GetRec_IsotropicKplane(Qmat, Plane);
    std::vector<MyMatrix<T>> GRP1 = f_stab
    
        GRP1:=f_stab(eRec);
        GRP2:=eRec.MapOrthogonalSublatticeGroup(GRP1);
        eInvRed:=INDEF_FORM_Invariant_IsotropicKplane_Raw(Qmat, Plane);
        return [GetRationalInvariant(GRP2), eInvRed];
    end;

  
  InvariantIsotropic<T,Tint> INDEF_FORM_Invariant_IsotropicKplane(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }
  InvariantIsotropic<T,Tint> INDEF_FORM_Invariant_IsotropicKflag(MyMatrix<T> const& Q, MyMatrix<Tint> const& Plane) {
  }

public:
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
