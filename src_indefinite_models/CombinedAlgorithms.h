// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
#define SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_

// clang-format off
#include "ApproximateModels.h"
#include "lorentzian_perfect.h"
#include "MatrixGroup.h"
// clang-format on

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
struct InvariantIsotropic {
};

template<typename T>
T GetRationalInvariant(std::vector<MyMatrix<T>> const& ListGen) {
  std::set<T> set_den;
  for (auto & eGen : ListGen) {
    if (!IsIntegralMatrix(eGen)) {
      int n = eGen.rows();
      for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
          T eDen = GetDenominator(eGen(i,j));
          set_den.insert(eDen);
        }
      }
    }
  }
  std::set<T> primes;
  for (auto & eDen : set_den) {
    std::map<T, size_t> l_primes = FactorsIntMap(eDen);
    for (auto & kv: l_primes) {
      primes.insert(kv.first);
    }
  }
  T prod(1);
  for (auto & p: primes) {
    prod *= p;
  }
  return prod;
}

template<typename T>
MyMatrix<T> ExpandMatrix(MyMatrix<T> const& M) {
  int n = M.rows();
  MyMatrix<T> TheBigMat = IdentityMat<T>(n+1);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      TheBigMat(i,j) = M(i,j);
    }
  }
  return TheBigMat;
}


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
