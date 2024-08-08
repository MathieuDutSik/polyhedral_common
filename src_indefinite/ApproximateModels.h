// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_
#define SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_

// clang-format off
#include "Isotropic.h"
#include "GRP_GroupFct.h"
#include "MAT_MatrixInt.h"
#include "LatticeStabEquiCan.h"
#include "IndefiniteFormFundamental.h"
#include "MatrixGroup.h"
#include <memory>
// clang-format on

#ifdef DEBUG
#define DEBUG_APPROXIMATE_MODELS
#endif


template<typename T, typename Tint>
struct ApproximateModel {
  std::function<std::vector<MyMatrix<Tint>>(std::ostream&)> GetApproximateGroup;
  std::function<void(std::vector<MyMatrix<Tint>>,std::ostream&)> SetListClassesOrbitwise;
  std::function<std::vector<MyVector<Tint>>(T,std::ostream&)> GetCoveringOrbitRepresentatives;
  std::function<std::optional<MyVector<Tint>>(T,std::ostream&)> GetOneOrbitRepresentative;
};

// The Eichler transvection are defined for an even lattice (so integral and even norms of vectors.
// For such a lattice L* is defined as Q^{-1} Z^n
// The formula describing them is
// E_{f,x}(y) = y + (y,x) f - (x,x) / 2 (y,f) f - (y,f) x
// The Eichler transvection preserve the lattice L because (x,x) / 2 in Z, (y,x) in Z and (y,f) in Z.
// Also if y in L* since the above is still true, we have that E_{f,x} preserve L* and the image of
// the transformation in L*/L is the identity.
//
// Test of composition property:
// We have
// E_{f,y}(z) = z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y
// E_{f,x}E_{f,y}(z) =
//   (z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y)
//   + (z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y, x) f
//   - (x,x) / 2 (z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y, f) f
//   - (z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y, f) x
// = z + (z,y) f - (y,y) / 2 (z,f) f - (z,f) y
//   + (z - (z,f) y, x) f
//   - (x,x) / 2 (z,f) f
//   - (z,f) x
// = z + (z,x + y) f - (z,f) (x+y) - (y,y)/2 (z,f) f - (x,x)/2 (z,f) f
//     -(z,f) (y,x) f
// = z + (z,x + y) f - (z,f) (x+y) - (x+y,x+y)/2 (z,f) f
template<typename T, typename Tint>
MyMatrix<Tint> INDEF_FORM_Eichler_Transvection(MyMatrix<T> const& Qmat, MyVector<Tint> const& f, MyVector<Tint> const& x) {
  if (!INDEF_FORM_IsEven(Qmat)) {
    std::cerr << "The lattice Qmat should be even in order to define the Eichler transvection\n";
    throw TerminalException{1};
  }
  int n = Qmat.rows();
  MyVector<T> x_T = UniversalVectorConversion<T,Tint>(x);
  MyVector<T> f_T = UniversalVectorConversion<T,Tint>(f);
#ifdef DEBUG_APPROXIMATE_MODELS
  T fNorm = EvaluationQuadForm(Qmat, f);
  MyVector<T> eProd = Qmat * x_T;
  T scal = f_T.dot(eProd);
  if (fNorm != 0 || scal != 0) {
    std::cerr << "eNorm or scal are inconsistent\n";
    throw TerminalException{1};
  }
#endif
  T xNorm = EvaluationQuadForm(Qmat, x);
  MyMatrix<Tint> RetMat(n,n);
  for (int u=0; u<n; u++) {
    MyVector<Tint> eImg = ZeroVector<Tint>(n);
    MyVector<Tint> y = ZeroVector<Tint>(n);
    y(u) = 1;
    eImg += y;
    //
    T scal1 = ScalarProductQuadForm(Qmat, y, x);
    Tint scal1_tint = UniversalScalarConversion<Tint,T>(scal1);
    eImg += scal1_tint * f;
    T scal2 = ScalarProductQuadForm(Qmat, y, f);
    T scal3 = (xNorm/2) * scal2;
    Tint scal2_tint = UniversalScalarConversion<Tint,T>(scal2);
    Tint scal3_tint = UniversalScalarConversion<Tint,T>(scal3);
    eImg -= scal3_tint * f;
    eImg -= scal2_tint * x;
    AssignMatrixRow(RetMat, u, eImg);
  }
#ifdef DEBUG_APPROXIMATE_MODELS
  MyMatrix<T> RetMat_T = UniversalMatrixConversion<T,Tint>(RetMat);
  MyMatrix<T> prod = RetMat_T * Qmat * RetMat_T.transpose();
  if (prod != Qmat) {
    std::cerr << "RetMat should prerve the Qmat\n";
    throw TerminalException{1};
  }
#endif
  return RetMat;
}


template<typename T, typename Tint>
std::optional<MyVector<Tint>> INDEF_FindIsotropic(MyMatrix<T> const& M, std::ostream& os) {
  int rnk = RankMat(M);
  int n = M.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: INDEF_FindIsotropic, M=\n";
  WriteMatrix(os, M);
  os << "MODEL: rnk=" << rnk << " n=" << n << "\n";
#endif
  if (rnk < n) {
    MyMatrix<T> NSP_T = NullspaceIntMat(M);
    MyVector<T> eV_T = GetMatrixRow(NSP_T, 0);
    MyVector<Tint> eV = UniversalVectorConversion<Tint,T>(eV_T);
    return eV;
  }
  std::optional<MyVector<T>> opt = FindIsotropic(M, os);
  if (opt) {
    MyVector<T> const& eV = *opt;
    MyVector<T> fV = RemoveFractionVector(eV);
    MyVector<Tint> gV = UniversalVectorConversion<Tint,T>(fV);
    return gV;
  } else {
    return {};
  }
}


template<typename Tint>
std::vector<MyMatrix<Tint>> GeneratorsSL2Z() {
  MyMatrix<Tint> eGenS = ZeroMatrix<Tint>(2,2);
  eGenS(0,1) = -1;
  eGenS(1,0) = 1;
  MyMatrix<Tint> eGenT = ZeroMatrix<Tint>(2,2);
  eGenT(0,0) = 1;
  eGenT(0,1) = 1;
  eGenT(1,1) = 1;
  return {eGenS, eGenT};
}


template<typename T, typename Tint>
struct InternalEichler {
  MyMatrix<T> Qmat;
  MyMatrix<T> Gmat;
  std::vector<MyVector<Tint>> ListClasses;
  std::vector<MyMatrix<Tint>> GRPstart;
};


template<typename T>
std::vector<T> GetSquareDivisors(T const& X) {
  std::map<T, size_t> map = FactorsIntMap(T_abs(X));
  std::vector<std::vector<T>> ListListProd;
  int64_t two(2);
  for (auto & kv : map) {
    T const& p = kv.first;
    int64_t mult = kv.second;
    int64_t mult_div = QuoInt(mult, two);
    std::vector<T> ListProd;
    T val(1);
    for (int64_t u=0; u<=mult_div; u++) {
      ListProd.push_back(val);
      val *= p;
    }
    ListListProd.push_back(ListProd);
  }
  std::vector<T> ListSquareDiv{ T(1) };
  for (auto & ListProd : ListListProd) {
    std::vector<T> V;
    for (auto & x : ListSquareDiv) {
      for (auto & y : ListProd) {
        T z = x * y;
        V.push_back(z);
      }
    }
    ListSquareDiv = V;
  }
  return ListSquareDiv;
}


/*
  Based on paragraph 10 of Eichler book
  Quadratische Formen und orthogonale gruppen
  Also used the Scattone thesis as basis.

  Note: There are other works based on Eichler work that are analogous:
  ---Wall CTC, On The Orthogonal Groups of Unimodular Quadratic Forms
  ---James D.G. Indefinite Quadratic Forms of Determinant pm2p
  They also require the same hypothesis of having two hyperplanes.
*/
template<typename T, typename Tint, typename Tgroup>
ApproximateModel<T,Tint> INDEF_FORM_EichlerCriterion_TwoHyperplanesEven(MyMatrix<T> const& Qmat) {
  int n = Qmat.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
  std::vector<int> LIdx{1,0,3,2};
  for (int iRow=0; iRow<4; iRow++) {
    int iColCrit = LIdx[iRow];
    for (int iCol=0; iCol<4; iCol++) {
      if (iCol == iColCrit) {
        if (Qmat(iRow,iCol) != 1) {
          std::cerr << "We do not have the structure of 2U, A\n";
          throw TerminalException{1};
        }
      } else {
        if (Qmat(iRow,iCol) != 0) {
          std::cerr << "We do not have the structure of 2U, B\n";
          throw TerminalException{1};
        }
      }
    }
  }
  for (int iRow=0; iRow<4; iRow++) {
    for (int iCol=4; iCol<n; iCol++) {
      if (Qmat(iRow,iCol) != 0) {
        std::cerr << "We do not have the structure of 2U, C\n";
        throw TerminalException{1};
      }
    }
  }
  if (!INDEF_FORM_IsEven(Qmat)) {
    std::cerr << "MODIND: The lattice is not even\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> Gmat(n-4, n-4);
  for (int i=4; i<n; i++) {
    for (int j=4; j<n; j++) {
      Gmat(i-4, j-4) = Qmat(i,j);
    }
  }
  // ComputeClasses
  std::vector<MyVector<Tint>> ListClasses = ComputeTranslationClasses<T,Tint>(Gmat);
  std::vector<MyMatrix<Tint>> GRPstart;
  InternalEichler<T,Tint> ie{Qmat, Gmat, ListClasses, GRPstart};
  std::shared_ptr<InternalEichler<T,Tint>> shr_ptr = std::make_shared<InternalEichler<T,Tint>>(ie);
  std::function<void(std::vector<MyMatrix<Tint>> const&, std::ostream&)> SetListClassesOrbitwise = [=](std::vector<MyMatrix<Tint>> const& GRPmatr, [[maybe_unused]] std::ostream& os) -> void {
    shr_ptr->GRPstart = GRPmatr;
    MyMatrix<T> const& Qmat = shr_ptr->Qmat;
    int n = Qmat.rows();
    std::vector<MyVector<Tint>> ListClassesExt;
    for (auto & eV : ListClasses) {
      MyVector<Tint> NewV = ZeroVector<Tint>(n);
      for (int i=0; i<n-4; i++) {
        NewV(i+4) = eV(i);
      }
      ListClassesExt.emplace_back(std::move(NewV));
    }
    auto GetPosition=[&](MyVector<Tint> const& eVect) -> size_t {
      for (size_t i=0; i<ListClassesExt.size(); i++) {
        MyVector<Tint> eDiff = ListClassesExt[i] - eVect;
        MyVector<T> eDiff_T = UniversalVectorConversion<T,Tint>(eDiff);
        std::optional<MyVector<T>> opt = SolutionIntMat(Qmat, eDiff_T);
        if (opt) {
          return i;
        }
      }
      std::cerr << "Failed to find a matching coset in GetPosition\n";
      throw TerminalException{1};
    };
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    std::vector<Telt> ListPermGens;
    for (auto & eMatrGen : GRPmatr) {
      std::vector<Tidx> eList;
      for (auto & eClassExt : ListClassesExt) {
        MyVector<Tint> x_eM = eMatrGen.transpose() * eClassExt;
        eList.push_back(GetPosition(x_eM));
      }
      Telt ePerm(eList);
      ListPermGens.push_back(ePerm);
    }
    int n_classes = ListClassesExt.size();
    Tgroup GRPperm(ListPermGens, n_classes);
    std::vector<size_t> ListRepr = DecomposeOrbitPoint_FullRepr<Tgroup>(GRPperm);
    std::vector<MyVector<Tint>> ListClasses;
    for (auto & iRepr : ListRepr) {
      MyVector<Tint> eV = shr_ptr->ListClasses[iRepr];
      ListClasses.push_back(eV);
    }
    shr_ptr->ListClasses = ListClasses;
  };
  std::function<std::vector<MyMatrix<Tint>>(std::ostream&)> GetApproximateGroup = [=]([[maybe_unused]] std::ostream& os) -> std::vector<MyMatrix<Tint>> {
    MyMatrix<T> Qmat = shr_ptr->Qmat;
    int n = shr_ptr->Qmat.rows();
    // The quadratic form 2 x1 x2 + 2 x3 x4
    // We formulate it as (1/2) det U(x1,x2,x3,x4)
    // with U(x1,x2,x3,x4) = [ [x1 , -x4] , [x3 , x2] ]
    // If we multiply by a matrix in SL(2,Z) on the left and right
    auto GetLeftMultiplication = [&](MyMatrix<Tint> const& P) -> MyMatrix<Tint> {
      // We compute the product [[a, b],[c,d]] U(x1,x2,x3,x4)
      // This gives us
      // [[a x1 + b x3 , -a x4 + b x2] , [ c x1 + d x3 , -c x4 + d x2]]
      // So Y = XA,
      // Y1 =  a x1 + b x3
      // Y2 = -c x4 + d x2
      // Y3 =  c x1 + d x3
      // Y4 =  a x4 - b x2
      Tint a = P(0,0);
      Tint b = P(0,1);
      Tint c = P(1,0);
      Tint d = P(1,1);
      MyMatrix<Tint> RetMat = ZeroMatrix<Tint>(4,4);
      RetMat(0,0) = a;
      RetMat(0,2) = c;
      RetMat(1,1) = d;
      RetMat(1,3) = -b;
      RetMat(2,0) = b;
      RetMat(2,2) = d;
      RetMat(3,1) = -c;
      RetMat(3,3) = a;
      return RetMat;
    };
    auto GetRightMultiplication = [&](MyMatrix<Tint> const& P) -> MyMatrix<Tint> {
      // We compute the product [ [x1 , -x4] , [x3 , x2] ]      [[a, b],[c,d]]
      // This gives us
      // [[a x1 - c x4 , b x1 - d x4] , [a x3 + c x2 , b x3 + d x2]]
      // So Y = XA,
      // Y1 = a x1 - c x4
      // Y2 = b x3 + d x2
      // Y3 = a x3 + c x2
      // Y4 = -b x1 + d x4
      Tint a = P(0,0);
      Tint b = P(0,1);
      Tint c = P(1,0);
      Tint d = P(1,1);
      MyMatrix<Tint> RetMat = ZeroMatrix<Tint>(4,4);
      RetMat(0,0) = a;
      RetMat(0,3) = -b;
      RetMat(1,1) = d;
      RetMat(1,2) = c;
      RetMat(2,1) = b;
      RetMat(2,2) = a;
      RetMat(3,0) = -c;
      RetMat(3,3) = d;
      return RetMat;
    };
    auto ExpandGenerator=[&](MyMatrix<Tint> const& eGen) -> MyMatrix<Tint> {
      MyMatrix<Tint> BigMat = IdentityMat<Tint>(n);
      for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
          BigMat(i,j) = eGen(i,j);
        }
      }
      return BigMat;
    };
    std::vector<MyMatrix<Tint>> ListGenerators;
    auto FuncInsert=[&](MyMatrix<Tint> const& TheMat) -> void {
#ifdef DEBUG_APPROXIMATE_MODELS
      MyMatrix<T> TheMat_T = UniversalMatrixConversion<T,Tint>(TheMat);
      MyMatrix<T> eProd = TheMat_T * Qmat * TheMat_T.transpose();
      if (eProd != Qmat) {
        std::cerr << "The matrix is not preserving the quadratic form\n";
        throw TerminalException{1};
      }
#endif
      ListGenerators.push_back(TheMat);
    };
    for (auto & eGen : shr_ptr->GRPstart) {
      FuncInsert(eGen);
    }
    for (auto & eGen : GeneratorsSL2Z<Tint>()) {
      FuncInsert(ExpandGenerator(GetLeftMultiplication(eGen)));
      FuncInsert(ExpandGenerator(GetRightMultiplication(eGen)));
    }
    // Now looking at the isotropic vectors, generating the Eichler transvections
    for (int i=0; i<4; i++) {
      MyVector<Tint> eVect = ZeroVector<Tint>(n);
      eVect(i) = 1;
      MyVector<T> eVect_T = UniversalVectorConversion<T,Tint>(eVect);
      MyVector<T> eProd = Qmat * eVect_T;
      MyMatrix<T> eProd_M(1,n);
      for (int i=0; i<n; i++) {
        eProd_M(0,i) = eProd(i);
      }
      MyMatrix<T> BasisOrth = NullspaceIntMat(eProd_M);
      for (int i=0; i<BasisOrth.rows(); i++) {
        MyVector<T> eOrth_T = GetMatrixRow(BasisOrth, i);
        MyVector<Tint> eOrth = UniversalVectorConversion<Tint,T>(eOrth_T);
        FuncInsert(INDEF_FORM_Eichler_Transvection(Qmat, eVect, eOrth));
      }
    }
    return ListGenerators;
  };
  // Notion of divisor
  // For an element y in a lattice L, define div(y) the integer such that
  // (L,y) = div(y) Z.
  // For each v in L* we can find smallest d such that w = dv in L.
  // And we then have div(w) = d
  std::function<std::vector<MyVector<Tint>>(T const&, std::ostream&)> EnumerateVectorOverDiscriminant = [=](T const& X,[[maybe_unused]] std::ostream& os) -> std::vector<MyVector<Tint>> {
    // For each vector v in M* / M
    // We want (x, v) divided by d.
    //
    // We should have e_I.x = d alpha_I
    // We should have (x,x) = X
    // For the first coordinates, we can write them as (0 , 0 , d , du)
    // By changing the value of u, we can shift the value of u which shifts the value
    // by 2 d^2, + or -.
    // The lattice is even so the 2 is not a problem.
    //
    // So, if we have some X, we want to find u\in Z and v in Z^{n-4} such that X = 2d^2 u + A[v]
    // with A the Gram matrix of M.
    // v must actually be of the form v = v0 + d A^{-1} h.
    // We are only interested on the value modulo 2d^2 of A on those vectors. So, it is a finite problem.
    // But it is a hard one as the search space might be very large.
    //
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of EnumerateVectorOverDiscriminant X=" << X << "\n";
#endif
    MyMatrix<T> const& Qmat = shr_ptr->Qmat;
    MyMatrix<T> const& Gmat = shr_ptr->Gmat;
    int n = Qmat.rows();
    T two(2);
    if (ResInt(X, two) == T(1)) {
      return {};
    }
    T Xdiv2 = X / 2;
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: EnumerateVectorOverDiscriminant Xdiv2=" << Xdiv2 << "\n";
#endif
    // The two hyperbolic planes are unimodular so for the discriminant, we only need
    // to consider the Gmat part.
    MyMatrix<T> Ginv = Inverse(Gmat);
    std::vector<MyVector<Tint>> ListSolution;
    for (auto & eClass1 : ListClasses) {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: EnumerateVectorOverDiscriminant eClass1=" << StringVectorGAP(eClass1) << "\n";
#endif
      // The vector eClass1 is defined modulo an element of Gmat Z^n
      // Thus eClass2 is defined modulo an element of Z^n. This is an elements of L*
      // Thus eClass3 is defined modulo an element of d Z^n
      // So, we can write v = eClass3 + d w
      // which gets us A[v] = A[eClass3] + 2d w^T Gmat eClass3 + d^2 A[w]
      // So, the problem is actually linear.
      MyVector<T> eClass1_T = UniversalVectorConversion<T,Tint>(eClass1);
      MyVector<T> eClass2 = Ginv * eClass1_T;
      T DivD_frac = RemoveFractionVectorPlusCoeff(eClass2).TheMult;
      T DivD = GetDenominator(DivD_frac);
      T DivD_sqr = DivD * DivD;
      MyVector<T> eClass3 = DivD * eClass2;
      T X_res1 = ResInt(Xdiv2, DivD);
      T X_res2 = ResInt(Xdiv2, DivD_sqr);
      T Aclass3_norm = EvaluationQuadForm(Gmat, eClass3) / 2;
      T Aclass3_res1 = ResInt(Aclass3_norm, DivD);
      T Aclass3_res2 = ResInt(Aclass3_norm, DivD_sqr);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: EnumerateVectorOverDiscriminant Aclass3_norm=" << Aclass3_norm << "\n";
#endif
      if (Aclass3_res1 == X_res1) { // if not we cannot find a solution
        // and so
        // (X - A[eClass3])/2 = 2d w^T Gmat eClass3 + ....
        T diff = (X_res2 - Aclass3_res2) / DivD;
        MyVector<T> eProd = Gmat * eClass3;
        GCD_dot<T> eRec = ComputeGcdDot(eProd);
        auto iife_quot=[&]() -> T {
          if (eRec.gcd == 0) {
            return T(0);
          } else {
            return diff / eRec.gcd;
          }
        };
        T quot = iife_quot();
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: EnumerateVectorOverDiscriminant quot=" << quot << "\n";
#endif
        if (IsInteger(quot)) {
          MyVector<T> w = quot * eRec.V;
          MyVector<T> eClass4 = eClass3 + DivD * w;
          T Aclass4_norm = EvaluationQuadForm(Gmat, eClass4) / 2;
          T Aclass4_res2 = ResInt(Aclass4_norm, DivD_sqr);
#ifdef DEBUG_APPROXIMATE_MODELS
          if (Aclass4_res2 != X_res2) {
            std::cerr << "A bug to resolve\n";
            throw TerminalException{1};
          }
#endif
          T u = (Xdiv2 - Aclass4_norm) / DivD_sqr;
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: EnumerateVectorOverDiscriminant u=" << u << "\n";
#endif
          MyVector<T> eSolution_T = ZeroVector<T>(n);
          eSolution_T(2) = DivD;
          eSolution_T(3) = DivD * u;
          for (int u=0; u<n-4; u++) {
            eSolution_T(u+4) = eClass4(u);
          }
#ifdef DEBUG_APPROXIMATE_MODELS
          T eNorm = EvaluationQuadForm(Qmat, eSolution_T);
          if (eNorm != X) {
            std::cerr << "eNorm / X is inconsistent\n";
            throw TerminalException{1};
          }
#endif
          MyVector<Tint> eSolution = UniversalVectorConversion<Tint,T>(eSolution_T);
          ListSolution.push_back(eSolution);
        }
      }
    }
    return ListSolution;
  };
  std::function<std::vector<MyVector<Tint>>(T const&,std::ostream&)> GetCoveringOrbitRepresentatives = [=](T const& X, [[maybe_unused]] std::ostream& os) -> std::vector<MyVector<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of GetCoveringOrbitRepresentatives 1: X=" << X << "\n";
#endif
    if (X == 0) {
      return EnumerateVectorOverDiscriminant(0, os);
    }
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Before GetSquareDivisor\n";
#endif
    std::vector<T> ListDiv = GetSquareDivisors(X);
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: After GetSquareDivisor\n";
#endif
    std::vector<MyVector<Tint>> ListSolution;
    for (auto & eDiv_T : ListDiv) {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: eDiv_T=" << eDiv_T << "\n";
#endif
      Tint eDiv = UniversalScalarConversion<Tint,T>(eDiv_T);
      T Xtarget = X / (eDiv_T * eDiv_T);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Xtarget=" << Xtarget << "\n";
#endif
      for (auto & eSol : EnumerateVectorOverDiscriminant(Xtarget, os)) {
        MyVector<Tint> eV = eDiv * eSol;
        ListSolution.push_back(eV);
      }
    }
    return ListSolution;
  };
  std::function<std::optional<MyVector<Tint>>(T const&,std::ostream&)> GetOneOrbitRepresentative = [=](T const& X, [[maybe_unused]] std::ostream& os) -> std::optional<MyVector<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of GetOneOrbitRepresentative 1: X=" << X << "\n";
#endif
    T two(2);
    T res = ResInt(X, two);
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: GetOneOrbitRepresentative res=" << res << "\n";
#endif
    if (res == T(1)) {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: No solution\n";
#endif
      return {};
    }
    T Xdiv2 = X / 2;
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Xdiv2=" << Xdiv2 << "\n";
#endif
    MyMatrix<T> const& Qmat = shr_ptr->Qmat;
    int n = Qmat.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: n=" << n << "\n";
#endif
    MyVector<Tint> eSol = ZeroVector<Tint>(n);
    eSol(0) = 1;
    eSol(1) = Xdiv2;
    return eSol;
  };
  return {GetApproximateGroup, SetListClassesOrbitwise, GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative};
}


/*
  We compute some isometries of the lattice coming from the blocks.
  The only constraint is that the matrices being returned are isometries.
  
 */
template<typename T, typename Tint>
std::vector<MyMatrix<Tint>> GetEasyIsometries(MyMatrix<T> const& Qmat, std::ostream& os) {
  std::vector<MyMatrix<Tint>> ListGenerators;
  auto f_insert=[&](MyMatrix<Tint> const& P) -> void {
#ifdef DEBUG_APPROXIMATE_MODELS
    MyMatrix<T> P_T = UniversalMatrixConversion<T,Tint>(P);
    MyMatrix<T> prod = P_T * Qmat * P_T.transpose();
    if (prod != Qmat) {
      std::cerr << "MODEL: P does not preserve the quadratic form\n";
      throw TerminalException{1};
    }
#endif
    ListGenerators.push_back(P);
  };
  size_t n = Qmat.rows();
  ListGenerators.push_back(IdentityMat<Tint>(n));
  //
  GraphBitset eG(n);
  for (size_t i=0; i<n; i++) {
    for (size_t j=i+1; j<n; j++) {
      if (Qmat(i, j) != 0) {
        eG.AddAdjacent(i, j);
        eG.AddAdjacent(j, i);
      }
    }
  }
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(eG);
  std::vector<MyMatrix<T>> ListQ;
  std::vector<size_t> ListPosIdx;
  for (size_t iConn=0; iConn<LConn.size(); iConn++) {
    std::vector<size_t> const& eConn = LConn[iConn];
    size_t dim = eConn.size();
    MyMatrix<T> eQ(dim, dim);
    for (size_t i=0; i<dim; i++) {
      for (size_t j=i+1; j<n; j++) {
        eQ(i,j) = Qmat(eConn[i], eConn[j]);
      }
    }
    ListQ.push_back(eQ);
    if (IsPositiveDefinite(eQ)) {
      ListPosIdx.push_back(iConn);
    }
  }
  for (auto & iConn : ListPosIdx) {
    MyMatrix<T> const& eQ = ListQ[iConn];
    std::vector<size_t> const& eConn = LConn[iConn];
    size_t dim = eConn.size();
    std::vector<MyMatrix<Tint>> LGen = ArithmeticAutomorphismGroup<T,Tint>(eQ, os);
    for (auto & eGen : LGen) {
      MyMatrix<Tint> eBigGen = IdentityMat<Tint>(n);
      for (size_t i=0; i<dim; i++) {
        for (size_t j=0; j<dim; j++) {
          eBigGen(eConn[i], eConn[j]) = eGen(i,j);
        }
      }
      f_insert(eBigGen);
    }
  }
  /*
    A = (A1 0)
        (0 A2)
    now P A1 P^T = A2
    ----
    W = (0 R)
        (S 0)
    We have
    W A W^T = (0 R) (A1 0) (0   S^T)
              (S 0) (0 A2) (R^T   0)
            = (0    R A2) (0   S^T)
              (S A1    0) (R^T   0)
            = (R A2 R^T    0    )
              (0        S A1 S^T)
    S = P
    R = P^{-1}
   */
  size_t nbConnRed = ListPosIdx.size();
  for (size_t iConnRed=0; iConnRed<nbConnRed; iConnRed++) {
    size_t iConn = ListPosIdx[iConnRed];
    MyMatrix<T> const& eQ = ListQ[iConn];
    std::vector<size_t> const& eConn = LConn[iConn];
    for (size_t jConnRed=iConnRed+1; jConnRed<nbConnRed; jConnRed++) {
      size_t jConn = ListPosIdx[jConnRed];
      MyMatrix<T> const& fQ = ListQ[jConn];
      std::vector<size_t> const& fConn = LConn[jConn];
      size_t dim = eConn.size();
      std::optional<MyMatrix<Tint>> opt = ArithmeticEquivalence<T,Tint>(eQ, fQ, os);
      if (opt) {
        MyMatrix<Tint> const& P = *opt;
        MyMatrix<Tint> Pinv = Inverse(P);
        MyMatrix<Tint> BigP = IdentityMat<Tint>(n);
        for (size_t i=0; i<dim; i++) {
          for (size_t j=0; j<dim; j++) {
            BigP(eConn[i], eConn[j]) = 0;
            BigP(fConn[i], fConn[j]) = 0;
            BigP(eConn[i], fConn[j]) = P(i,j);
            BigP(fConn[i], eConn[j]) = Pinv(i,j);
          }
        }
        f_insert(BigP);
      }
    }
  }
  return ListGenerators;
}

/*
  Ideally, given an isotropic vector v, and c, we want to find a vector w such that
  v.w = c and w is isotrop.
  ----
  Compute V = v^{perp}. We want to find the vectors in w in V with phi(w) = c for some linear
  form v. We can work in the quotient W = V / v and the function phi lifts to W.
  If the initial space has signature (p,q), then the resulting space W has signature (p-1,q-1).
 */
template<typename T, typename Tint>
std::pair<MyVector<Tint>, MyVector<Tint>> GetHyperbolicPlane(MyMatrix<T> const& Qmat, std::ostream& os) {
  int n = Qmat.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: Beginning of GetHyperbolicPlane, Qmat=\n";
  WriteMatrix(os, Qmat);
#endif
  std::set<MyVector<Tint>> SetVect;
  MyMatrix<Tint> ThePerturb = IdentityMat<Tint>(n);

  int n_iter = 100;
  auto f_insert_vect=[&]() -> MyVector<Tint> {
    int iter = 0;
    while(true) {
      MyMatrix<Tint> ePerturb = GetRandomMatrixPerturbation<Tint>(n);
      ThePerturb = ThePerturb * ePerturb;
      MyMatrix<T> ThePerturb_T = UniversalMatrixConversion<T,Tint>(ThePerturb);
      MyMatrix<T> M = ThePerturb_T * Qmat * ThePerturb_T.transpose();
      std::optional<MyVector<Tint>> opt = INDEF_FindIsotropic<T,Tint>(M, os);
      MyVector<Tint> eVect = unfold_opt(opt, "Failed to find an isotropic vector");
      MyVector<Tint> NewVect = ThePerturb.transpose() * eVect;
#ifdef DEBUG_APPROXIMATE_MODELS
      T sum = EvaluationQuadForm<T,Tint>(Qmat, NewVect);
      if (sum != 0) {
        std::cerr << "MODEL: NewVect is not an isotropic vector\n";
        throw TerminalException{1};
      }
      os << "MODEL: GetHyperbolicPlane NewVect=" << StringVectorGAP(NewVect) << "\n";
#endif
      if (SetVect.count(NewVect) == 0 || iter == n_iter) {
        SetVect.insert(NewVect);
        return NewVect;
      }
      iter += 1;
    }
  };

  auto has_non_orthogonal=[&](MyVector<Tint> const& NewVect) -> bool {
    MyVector<T> NewVect_T = UniversalVectorConversion<T,Tint>(NewVect);
    MyVector<T> prod = Qmat * NewVect_T;
    for (auto & eVect : SetVect) {
      MyVector<T> eVect_T = UniversalVectorConversion<T,Tint>(eVect);
      T scal = eVect_T.dot(prod);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: has_non_orthogonal, scal=" << scal << "\n";
#endif
      if (scal != 0) {
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: has_non_orthogonal, eVect=" << StringVectorGAP(eVect) << " NewVect=" << StringVectorGAP(NewVect) << "\n";
#endif
        return true;
      }
    }
    return false;
  };
  auto f_get_minimal=[&](MyVector<Tint> const& NewVect) -> std::pair<T, MyVector<Tint>> {
    T min_scal;
    bool is_first = true;
    MyVector<Tint> BestVect;
    //
    MyVector<T> NewVect_T = UniversalVectorConversion<T,Tint>(NewVect);
    MyVector<T> prod = Qmat * NewVect_T;
    for (auto & eVect : SetVect) {
      MyVector<Tint> rVect = eVect;
      MyVector<T> eVect_T = UniversalVectorConversion<T,Tint>(eVect);
      T scal = eVect_T.dot(prod);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: f_get_minimal, scal=" << scal << "\n";
#endif
      if (scal < 0) {
        scal = -scal;
        rVect = - eVect;
      }
      if (scal != 0) {
        if (!is_first) {
          if (scal < min_scal) {
            min_scal = scal;
            BestVect = rVect;
          }
        } else {
          min_scal = scal;
          BestVect = rVect;
          is_first = false;
        }
      }
    }
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: f_get_minimal, min_scal=" << min_scal << "\n";
    os << "MODEL: f_get_minimal, BestVect=" << StringVectorGAP(BestVect) << " NewVect=" << StringVectorGAP(NewVect) << "\n";
#endif
    return {min_scal, BestVect};
  };

  // Now iterating
  bool is_first = true;
  std::pair<MyVector<Tint>, MyVector<Tint>> pairA;
  T min_scal;
  int n_at_level = 0;
  int max_level = 10;
  while(true) {
    MyVector<Tint> NewVect = f_insert_vect();
    if (has_non_orthogonal(NewVect)) {
      std::pair<T, MyVector<Tint>> pairB = f_get_minimal(NewVect);
      if (is_first) {
        min_scal = pairB.first;
        pairA = {NewVect, pairB.second};
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: f_get_minimal, first, min_scal=" << min_scal << "\n";
        os << "MODEL: f_get_minimal, 1: pairA=(" << StringVectorGAP(pairA.first) << " / " << StringVectorGAP(pairA.second) << ")\n";
#endif
        n_at_level = 0;
        is_first = false;
      } else {
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: f_get_minimal, pairB.first=" << pairB.first << " min_scal=" << min_scal << "\n";
#endif
        if (pairB.first < min_scal) {
          min_scal = pairB.first;
          pairA = {NewVect, pairB.second};
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: f_get_minimal, now min_scal=" << min_scal << "\n";
          os << "MODEL: f_get_minimal, 2: pairA=(" << StringVectorGAP(pairA.first) << " / " << StringVectorGAP(pairA.second) << ")\n";
#endif
          n_at_level = 0;
        } else {
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: f_get_minimal, n_at_level=" << n_at_level << "\n";
#endif
          if (pairB.first == min_scal) {
            n_at_level += 1;
          }
          if (n_at_level == max_level) {
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: f_get_minimal, 3: pairA=(" << StringVectorGAP(pairA.first) << " / " << StringVectorGAP(pairA.second) << ")\n";
#endif
            return pairA;
          }
        }
      }
    }
  }
}

template<typename T, typename Tint>
MyMatrix<Tint> GetEichlerHyperplaneBasis(MyMatrix<T> const& Qmat, std::ostream& os) {
  MyMatrix<Tint> Basis1 = MatrixFromPairVector(GetHyperbolicPlane<T,Tint>(Qmat, os));
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: GetEichlerHyperplaneBasis, Basis1=\n";
  WriteMatrix(os, Basis1);
#endif
  MyMatrix<T> Basis1_T = UniversalMatrixConversion<T,Tint>(Basis1);
  MyMatrix<T> prod1 = Basis1_T * Qmat;
  MyMatrix<T> NSP_T = NullspaceIntTrMat(prod1);
  MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
  MyMatrix<T> Qmat2 = NSP_T * Qmat * NSP_T.transpose();
  MyMatrix<Tint> Basis2 = MatrixFromPairVector(GetHyperbolicPlane<T,Tint>(Qmat2, os));
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: GetEichlerHyperplaneBasis, Basis2=\n";
  WriteMatrix(os, Basis2);
#endif
  MyMatrix<Tint> Basis2_NSP = Basis2 * NSP;
  MyMatrix<Tint> HyperBasis = Concatenate(Basis1, Basis2_NSP);
  MyMatrix<T> HyperBasis_T = UniversalMatrixConversion<T,Tint>(HyperBasis);
  MyMatrix<T> prod2 = HyperBasis_T * Qmat;
  MyMatrix<T> NSP2_T = NullspaceIntTrMat(prod2);
  MyMatrix<Tint> NSP2 = UniversalMatrixConversion<Tint,T>(NSP2_T);
  MyMatrix<Tint> FullBasis = Concatenate(HyperBasis, NSP2);
#ifdef DEBUG_APPROXIMATE_MODELS
  Tint det = DeterminantMat(FullBasis);
  if (T_abs(det) != 1) {
    std::cerr << "MODEL: FullBasis should be of determinant 1\n";
    throw TerminalException{1};
  }
#endif
  return FullBasis;
}



template<typename T, typename Tint>
struct InternalApproxCaseA {
  MyMatrix<T> Qmat;
  MyMatrix<Tint> FullBasis;
  ApproximateModel<T,Tint> approx;
};


template<typename T, typename Tint>
struct InternalApproxCaseB {
  MyMatrix<T> Qmat;
  MyMatrix<Tint> FullBasis;
  MyMatrix<Tint> eEmbed;
  std::vector<MyMatrix<Tint>> ListCoset;
  std::vector<MyMatrix<Tint>> ListGenerators;
  ApproximateModel<T,Tint> approx;
};


template<typename T>
MyMatrix<T> GetTwoPlanes() {
  MyMatrix<T> M = ZeroMatrix<T>(4,4);
  std::vector<int> V{1,0,3,2};
  for (int u=0; u<4; u++) {
    int idxW = V[u];
    M(u,idxW) = 1;
  }
  return M;
}

template<typename T>
MyMatrix<T> GetSubBlock11(MyMatrix<T> const& M) {
  MyMatrix<T> Mret(4,4);
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      Mret(i,j) = M(i,j);
    }
  }
  return Mret;
}

template<typename T>
MyMatrix<T> GetSubBlock22(MyMatrix<T> const& M) {
  int n = M.rows();
  int dim = n - 4;
  MyMatrix<T> Mret(dim,dim);
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      Mret(i,j) = M(i+4,j+4);
    }
  }
  return Mret;
}

template<typename T>
MyMatrix<T> GetSubBlock12(MyMatrix<T> const& M) {
  int n = M.rows();
  int dim = n - 4;
  MyMatrix<T> Mret(4,dim);
  for (int i=0; i<4; i++) {
    for (int j=0; j<dim; j++) {
      Mret(i,j) = M(i,j+4);
    }
  }
  return Mret;
}

template<typename T>
MyMatrix<T> GetTwoEmbedding(MyMatrix<T> const& eBlock) {
  MyMatrix<T> eEmbed = IdentityMat<T>(4);
  T h1 = eBlock(0,1);
  T h2 = eBlock(2,3);
  eEmbed(0,0) = h1;
  eEmbed(2,2) = h2;
#ifdef DEBUG_APPROXIMATE_MODELS
  MyMatrix<T> TwoPlanes = GetTwoPlanes<T>();
  MyMatrix<T> eProd = eEmbed * TwoPlanes * eEmbed.transpose();
  if (eProd != eBlock) {
    std::cerr << "eEmbed does not provide an embedding\n";
    throw TerminalException{1};
  }
#endif
  return eEmbed;
}

template<typename T>
MyMatrix<T> AssembleTwoDiagBlock(MyMatrix<T> const& M1, MyMatrix<T> const& M2) {
  int n1 = M1.rows();
  int n2 = M2.rows();
  MyMatrix<T> Mret = ZeroMatrix<T>(n1+n2, n1+n2);
  for (int i=0; i<n1; i++) {
    for (int j=0; j<n1; j++) {
      Mret(i,j) = M1(i,j);
    }
  }
  for (int i=0; i<n2; i++) {
    for (int j=0; j<n2; j++) {
      Mret(i + n1,j + n1) = M1(i,j);
    }
  }
  return Mret;
}



template<typename T, typename Tint, typename Tgroup>
ApproximateModel<T,Tint> INDEF_FORM_GetApproximateModel(MyMatrix<T> const& Qmat, std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  int n = Qmat.rows();
  MyMatrix<Tint> FullBasis = GetEichlerHyperplaneBasis<T,Tint>(Qmat, os);
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: INDEF_FORM_GetApproximateModel, FullBasis=\n";
  WriteMatrix(os, FullBasis);
#endif
  MyMatrix<T> FullBasis_T = UniversalMatrixConversion<T,Tint>(FullBasis);
  MyMatrix<T> QmatRed = FullBasis_T * Qmat * FullBasis_T.transpose();
  MyMatrix<T> Block11 = GetSubBlock11(QmatRed);
  MyMatrix<T> TwoPlanes = GetTwoPlanes<T>();

  if (TwoPlanes == Block11) {
    InternalApproxCaseA<T,Tint> casea{Qmat, FullBasis, INDEF_FORM_EichlerCriterion_TwoHyperplanesEven<T,Tint,Tgroup>(QmatRed)};
    std::shared_ptr<InternalApproxCaseA<T,Tint>> shr_ptr_a = std::make_shared<InternalApproxCaseA<T,Tint>>(casea);
    std::vector<MyMatrix<Tint>> GRPeasy = GetEasyIsometries<T,Tint>(QmatRed, os);
    shr_ptr_a->approx.SetListClassesOrbitwise(GRPeasy, os);
    std::function<void(std::vector<MyMatrix<Tint>>,std::ostream&)> SetListClassesOrbitwise = [=]([[maybe_unused]] std::vector<MyMatrix<Tint>> const& GRPmatr, [[maybe_unused]] std::ostream& os) -> void {
      std::cerr << "That function should not be called\n";
      throw TerminalException{1};
    };
    std::function<std::vector<MyMatrix<Tint>>(std::ostream&)> GetApproximateGroup = [=]([[maybe_unused]] std::ostream& os) -> std::vector<MyMatrix<Tint>> {
      std::vector<MyMatrix<Tint>> ListGenerators;
      MyMatrix<Tint> const& FullBasis = shr_ptr_a->FullBasis;
      MyMatrix<Tint> FullBasisInv = Inverse(FullBasis);
      for (auto & eGen: shr_ptr_a->approx.GetApproximateGroup(os)) {
        MyMatrix<Tint> NewGen = FullBasisInv * eGen * FullBasis;
#ifdef DEBUG_APPROXIMATE_MODELS
        MyMatrix<T> const& Qmat = shr_ptr_a->Qmat;
        MyMatrix<T> NewGen_T = UniversalMatrixConversion<T,Tint>(NewGen);
        MyMatrix<T> prod = NewGen_T * Qmat * NewGen_T.transpose();
        if (prod != Qmat) {
          std::cerr << "NewGen is not preserving the matrix Qmat\n";
          throw TerminalException{1};
        }
#endif
        ListGenerators.push_back(NewGen);
      }
      return ListGenerators;
    };
    std::function<std::vector<MyVector<Tint>>(T,std::ostream&)> GetCoveringOrbitRepresentatives = [=](T const& X, [[maybe_unused]] std::ostream& os) -> std::vector<MyVector<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Beginning of GetCoveringOrbitRepresentatives 2: X=" << X << "\n";
#endif
      std::vector<MyVector<Tint>> ListSolution;
      MyMatrix<Tint> const& FullBasis = shr_ptr_a->FullBasis;
      for (auto & eVect : shr_ptr_a->approx.GetCoveringOrbitRepresentatives(X, os)) {
        MyVector<Tint> xNew = FullBasis.transpose() * eVect;
#ifdef DEBUG_APPROXIMATE_MODELS
        MyMatrix<T> const& Qmat = shr_ptr_a->Qmat;
        T eNorm = EvaluationQuadForm<T,Tint>(Qmat, xNew);
        if (eNorm != X) {
          std::cerr << "xNew is not of the right norm\n";
          throw TerminalException{1};
        }
#endif
        ListSolution.push_back(xNew);
      }
      return ListSolution;
    };
    std::function<std::optional<MyVector<Tint>>(T const&,std::ostream&)> GetOneOrbitRepresentative = [=](T const& X, [[maybe_unused]] std::ostream& os) -> std::optional<MyVector<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Beginning of GetOneOrbitRepresentative 2: X=" << X << "\n";
#endif
      MyMatrix<Tint> const& FullBasis = shr_ptr_a->FullBasis;
      std::optional<MyVector<Tint>> opt = shr_ptr_a->approx.GetOneOrbitRepresentative(X, os);
      if (opt) {
        MyVector<Tint> const& V = *opt;
        MyVector<Tint> xNew = FullBasis.transpose() * V;
#ifdef DEBUG_APPROXIMATE_MODELS
        MyMatrix<T> const& Qmat = shr_ptr_a->Qmat;
        T eNorm = EvaluationQuadForm<T,Tint>(Qmat, xNew);
        if (eNorm != X) {
          std::cerr << "xNew is not of the right norm\n";
          throw TerminalException{1};
        }
#endif
        return xNew;
      } else {
        return {};
      }
    };
    return {GetApproximateGroup, SetListClassesOrbitwise, GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative};
  }
  MyMatrix<T> Block12 = GetSubBlock12(QmatRed);
  if (IsZeroMatrix(Block12)) {
    MyMatrix<T> eEmbed11 = GetTwoEmbedding(Block11);
    MyMatrix<T> Block22 = GetSubBlock22(QmatRed);
    MyMatrix<T> QmatExt = AssembleTwoDiagBlock(TwoPlanes, Block22);
    MyMatrix<T> eEmbed_T = AssembleTwoDiagBlock(eEmbed11, IdentityMat<T>(n-4));
    MyMatrix<Tint> eEmbed = UniversalMatrixConversion<Tint,T>(eEmbed_T);
#ifdef DEBUG_APPROXIMATE_MODELS
    MyMatrix<T> prod = eEmbed_T * QmatExt * eEmbed_T.transpose();
    if (prod != QmatRed) {
      std::cerr << "eEmbed is not an embedding\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<T> eProdEmbed = Inverse(FullBasis_T) * eEmbed_T;
    MyMatrix<T> eProdEmbedInv = Inverse(eProdEmbed);
    ApproximateModel<T,Tint> approx = INDEF_FORM_EichlerCriterion_TwoHyperplanesEven<T,Tint,Tgroup>(QmatExt);
    std::vector<MyMatrix<Tint>> ListGen = approx.GetApproximateGroup(os);
    std::vector<MyMatrix<T>> ListGen_T = UniversalStdVectorMatrixConversion<T,Tint>(ListGen);
    int dim_ext = QmatExt.rows();
    GeneralMatrixGroupHelper<T,Telt> helper{dim_ext};
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: INDEF_FORM_GetApproximateModel, Before LinearSpace_Stabilizer_RightCoset\n";
#endif
    Stab_RightCoset<T> stab_right_coset = LinearSpace_Stabilizer_RightCoset<T,Tgroup,GeneralMatrixGroupHelper<T,Telt>>(ListGen_T, helper, eEmbed_T, os);
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: INDEF_FORM_GetApproximateModel, After LinearSpace_Stabilizer_RightCoset\n";
#endif
    std::vector<MyMatrix<Tint>> ListCoset = stab_right_coset.coset_desc.template expand<Tint>();
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: We have ListCoset\n";
#endif
    std::vector<MyMatrix<Tint>> ListGenerators;
    for (auto & eGen_T : stab_right_coset.list_gen) {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: eProdEmbed=\n";
      WriteMatrix(os, eProdEmbed);
      os << "MODEL: eGen_T=\n";
      WriteMatrix(os, eGen_T);
      os << "MODEL: eProdEmbedInv=\n";
      WriteMatrix(os, eProdEmbedInv);
#endif
      MyMatrix<T> fGen_T = eProdEmbed * eGen_T * eProdEmbedInv;
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: fGen_T=\n";
      WriteMatrix(os, fGen_T);
      if (!IsIntegralMatrix(fGen_T)) {
        std::cerr << "fGen should be integral\n";
        throw TerminalException{1};
      }
#endif
      MyMatrix<Tint> fGen = UniversalMatrixConversion<Tint,T>(fGen_T);
      ListGenerators.push_back(fGen);
    }
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: We have ListGenerators\n";
#endif

    InternalApproxCaseB<T,Tint> caseb{Qmat, FullBasis, eEmbed, ListCoset, ListGenerators, approx};
    std::shared_ptr<InternalApproxCaseB<T,Tint>> shr_ptr_b = std::make_shared<InternalApproxCaseB<T,Tint>>(caseb);
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: We have shr_ptr_b\n";
#endif


    std::function<void(std::vector<MyMatrix<Tint>>,std::ostream&)> SetListClassesOrbitwise = [=]([[maybe_unused]] std::vector<MyMatrix<Tint>> const& GRPmatr, [[maybe_unused]] std::ostream& os) -> void {
      std::cerr << "That function should not be called\n";
      throw TerminalException{1};
    };
    std::function<std::vector<MyMatrix<Tint>>(std::ostream&)> GetApproximateGroup = [=]([[maybe_unused]] std::ostream& os) -> std::vector<MyMatrix<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Beginning of GetApproximateGroup\n";
#endif
      return shr_ptr_b->ListGenerators;
    };
    std::function<std::vector<MyVector<Tint>>(T,std::ostream&)> GetCoveringOrbitRepresentatives = [=](T const& X, [[maybe_unused]] std::ostream& os) -> std::vector<MyVector<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Beginning of GetCoveringOrbitRepresentatives 3: X=" << X << "\n";
#endif
      std::vector<MyVector<Tint>> ListRepr;
      for (auto & eRepr : shr_ptr_b->approx.GetCoveringOrbitRepresentatives(X, os)) {
        for (auto & eCos : shr_ptr_b->ListCoset) {
          MyVector<Tint> fRepr = eCos.transpose() * eRepr;
          std::optional<MyVector<Tint>> opt = SolutionIntMat(shr_ptr_b->eEmbed, fRepr);
          if (opt) {
            MyVector<Tint> const& eSol = *opt;
            MyVector<Tint> fSol = shr_ptr_b->FullBasis.transpose() * eSol;
#ifdef DEBUG_APPROXIMATE_MODELS
            MyMatrix<T> const& Qmat = shr_ptr_b->Qmat;
            T eNorm = EvaluationQuadForm<T,Tint>(Qmat, fSol);
            if (eNorm != X) {
              std::cerr << "fSol is not of the right norm\n";
              throw TerminalException{1};
            }
#endif
            ListRepr.push_back(fSol);
          }
        }
      }
      return ListRepr;
    };
    std::function<std::optional<MyVector<Tint>>(T const&,std::ostream&)> GetOneOrbitRepresentative = [=](T const& X, [[maybe_unused]] std::ostream& os) -> std::optional<MyVector<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Beginning of GetOneOrbitRepresentative 3: X=" << X << "\n";
#endif
      std::optional<MyVector<Tint>> opt = shr_ptr_b->approx.GetOneOrbitRepresentative(X, os);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: We have opt\n";
#endif
      if (!opt) {
        return {};
      }
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Before iteration over |ListCoset|=" << ListCoset.size() << "\n";
#endif
      MyVector<Tint> const& eRepr = *opt;
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: eRepr=" << StringVectorGAP(eRepr) << "\n";
#endif
      for (auto & eCos : shr_ptr_b->ListCoset) {
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: eCos=" << StringMatrixGAP(eCos) << "\n";
#endif
        MyVector<Tint> fRepr = eCos.transpose() * eRepr;
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: fRepr=" << StringVectorGAP(fRepr) << "\n";
#endif
        std::optional<MyVector<Tint>> optB = SolutionIntMat(shr_ptr_b->eEmbed, fRepr);
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: We have optB\n";
#endif
        if (optB) {
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: Find a matching solution\n";
#endif
          MyVector<Tint> const& eSol = *optB;
          MyVector<Tint> fSol = shr_ptr_b->FullBasis.transpose() * eSol;
#ifdef DEBUG_APPROXIMATE_MODELS
          MyMatrix<T> const& Qmat = shr_ptr_b->Qmat;
          T eNorm = EvaluationQuadForm<T,Tint>(Qmat, fSol);
          if (eNorm != X) {
            std::cerr << "fSol is not of the right norm\n";
            throw TerminalException{1};
          }
#endif
          return fSol;
        }
      }
      // We wanted to avoid that but sometimes we cannot
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Nothing found, now enumerating the possibilities\n";
#endif
      std::vector<MyVector<Tint>> ListRepr = GetCoveringOrbitRepresentatives(X, os);
      if (ListRepr.size() > 0) {
        return ListRepr[0];
      }
      return {};
    };
    return {GetApproximateGroup, SetListClassesOrbitwise, GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative};
  }
  std::cerr << "Failed to find a relevant embedding\n";
  throw TerminalException{1};
}

template<typename T, typename Tint>
struct FirstNorm {
  T X;
  MyVector<Tint> eVect;
};

template<typename T, typename Tint>
FirstNorm<T,Tint> GetFirstNorm(ApproximateModel<T,Tint> const& approx, std::ostream& os) {
  T X = 1;
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: Beginning of GetFirstNorm\n";
#endif
  while(true) {
    std::optional<MyVector<Tint>> opt = approx.GetOneOrbitRepresentative(X, os);
    if (opt) {
      MyVector<Tint> const& V = *opt;
      return {X, V};
    }
    X += 1;
  }
}



// clang-format off
#endif  // SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_
// clang-format on
