// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_
#define SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_

// clang-format off
#include "Isotropic.h"
#include "GRP_GroupFct.h"
#include "MAT_MatrixInt.h"
#include <memory>
// clang-format on

#ifdef DEBUG
#define DEBUG_APPROXIMATE_MODELS
#endif


template<typename T>
struct ApproximateModel {
  std::function<std::vector<MyMatrix<T>>(void)> GetApproximateGroup;
  std::function<void(std::vector<MyMatrix<T>>)> SetListClassesOrbitwise;
  std::function<std::vector<MyVector<T>>(T)> GetCoveringOrbitRepresentatives;
  std::function<std::optional<MyVector<T>>(T)> GetOneOrbitRepresentative;
};

template<typename T>
bool INDEF_FORM_IsEven(MyMatrix<T> const& Qmat) {
  if (!IsIntegralMatrix(Qmat)) {
    return false;
  }
  T two(2);
  for (int u=0; u<Qmat.rows(); u++) {
    T res = ResInt(Qmat(u,u), two);
    if (res != 0) {
      return false;
    }
  }
  return true;
}

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
template<typename T>
MyMatrix<T> INDEF_FORM_Eichler_Transvection(MyMatrix<T> const& Qmat, MyVector<T> const& f, MyVector<T> const& x) {
  if (!INDEF_FORM_IsEven(Qmat)) {
    std::cerr << "The lattice Qmat should be even in order to define the Eichler transvection\n";
    throw TerminalException{1};
  }
  int n = Qmat.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
  T fNorm = EvaluationQuadForm(Qmat, f);
  MyVector<T> eProd = Qmat * x;
  T scal = f.dot(eProd);
  if (fNorm != 0 || scal != 0) {
    std::cerr << "eNorm or scal are inconsistent\n";
    throw TerminalException{1};
  }
#endif
  T xNorm = EvaluationQuadForm(Qmat, x);
  MyMatrix<T> RetMat(n,n);
  for (int u=0; u<n; u++) {
    MyVector<T> eImg = ZeroVector<T>(n);
    MyVector<T> y = ZeroVector<T>(n);
    y(u) = 1;
    eImg += y;
    //
    T scal_yx = ScalarProductQuadForm(Qmat, y, x);
    eImg += scal_yx * f;
    T scal_yf = ScalarProductQuadForm(Qmat, y, f);
    eImg -= (xNorm/2) * scal_yf * f;
    eImg -= scal_yf * x;
    AssignMatrixRow(RetMat, u, eImg);
  }
  return RetMat;
}



template<typename T>
std::vector<MyMatrix<T>> GeneratorsSL2Z() {
  MyMatrix<T> eGenS = ZeroMatrix<T>(2,2);
  eGenS(0,1) = -1;
  eGenS(1,0) = 1;
  MyMatrix<T> eGenT = ZeroMatrix<T>(2,2);
  eGenT(0,0) = 1;
  eGenT(0,1) = 1;
  eGenT(1,1) = 1;
  return {eGenS, eGenT};
}


template<typename T>
struct InternalEichler {
  MyMatrix<T> Qmat;
  MyMatrix<T> Gmat;
  std::vector<MyVector<T>> ListClasses;
  std::vector<MyMatrix<T>> GRPstart;
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
template<typename T, typename Tgroup>
ApproximateModel<T> INDEF_FORM_EichlerCriterion_TwoHyperplanesEven(MyMatrix<T> const& Qmat) {
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
  std::vector<MyVector<T>> ListClasses = ComputeTranslationClasses<T,T>(Gmat);
  std::vector<MyMatrix<T>> GRPstart;
  InternalEichler<T> ie{Qmat, Gmat, ListClasses, GRPstart};
  std::shared_ptr<InternalEichler<T>> shr_ptr = std::make_shared<InternalEichler<T>>(ie);
  std::function<void(std::vector<MyMatrix<T>>)> SetListClassesOrbitwise = [=](std::vector<MyMatrix<T>> const& GRPmatr) -> void {
    shr_ptr->GRPstart = GRPmatr;
    MyMatrix<T> const& Qmat = shr_ptr->Qmat;
    int n = Qmat.rows();
    std::vector<MyVector<T>> ListClassesExt;
    for (auto & eV : ListClasses) {
      MyVector<T> NewV = ZeroVector<T>(n);
      for (int i=0; i<n-4; i++) {
        NewV(i+4) = eV(i);
      }
      ListClassesExt.emplace_back(std::move(NewV));
    }
    auto GetPosition=[&](MyVector<T> const& eVect) -> size_t {
      for (size_t i=0; i<ListClassesExt.size(); i++) {
        MyVector<T> eDiff = ListClassesExt[i] - eVect;
        std::optional<MyVector<T>> opt = SolutionIntMat(Qmat, eDiff);
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
        MyVector<T> x_eM = eMatrGen.transpose() * eClassExt;
        eList.push_back(GetPosition(x_eM));
      }
      Telt ePerm(eList);
      ListPermGens.push_back(ePerm);
    }
    int n_classes = ListClassesExt.size();
    Tgroup GRPperm(ListPermGens, n_classes);
    std::vector<size_t> ListRepr = DecomposeOrbitPoint_FullRepr<Tgroup>(GRPperm);
    std::vector<MyVector<T>> ListClasses;
    for (auto & iRepr : ListRepr) {
      MyVector<T> eV = shr_ptr->ListClasses[iRepr];
      ListClasses.push_back(eV);
    }
    shr_ptr->ListClasses = ListClasses;
  };
  std::function<std::vector<MyMatrix<T>>(void)> GetApproximateGroup = [=]() -> std::vector<MyMatrix<T>> {
    MyMatrix<T> Qmat = shr_ptr->Qmat;
    int n = shr_ptr->Qmat.rows();
    // The quadratic form 2 x1 x2 + 2 x3 x4
    // We formulate it as (1/2) det U(x1,x2,x3,x4)
    // with U(x1,x2,x3,x4) = [ [x1 , -x4] , [x3 , x2] ]
    // If we multiply by a matrix in SL(2,Z) on the left and right
    auto GetLeftMultiplication = [&](MyMatrix<T> const& P) -> MyMatrix<T> {
      // We compute the product [[a, b],[c,d]] U(x1,x2,x3,x4)
      // This gives us
      // [[a x1 + b x3 , -a x4 + b x2] , [ c x1 + d x3 , -c x4 + d x2]]
      // So Y = XA,
      // Y1 =  a x1 + b x3
      // Y2 = -c x4 + d x2
      // Y3 =  c x1 + d x3
      // Y4 =  a x4 - b x2
      T a = P(0,0);
      T b = P(0,1);
      T c = P(1,0);
      T d = P(1,1);
      MyMatrix<T> RetMat = ZeroMatrix<T>(4,4);
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
    auto GetRightMultiplication = [&](MyMatrix<T> const& P) -> MyMatrix<T> {
      // We compute the product [ [x1 , -x4] , [x3 , x2] ]      [[a, b],[c,d]]
      // This gives us
      // [[a x1 - c x4 , b x1 - d x4] , [a x3 + c x2 , b x3 + d x2]]
      // So Y = XA,
      // Y1 = a x1 - c x4
      // Y2 = b x3 + d x2
      // Y3 = a x3 + c x2
      // Y4 = -b x1 + d x4
      T a = P(0,0);
      T b = P(0,1);
      T c = P(1,0);
      T d = P(1,1);
      MyMatrix<T> RetMat = ZeroMatrix<T>(4,4);
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
    auto ExpandGenerator=[&](MyMatrix<T> const& eGen) -> MyMatrix<T> {
      MyMatrix<T> BigMat = IdentityMat<T>(n);
      for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
          BigMat(i,j) = eGen(i,j);
        }
      }
      return BigMat;
    };
    std::vector<MyMatrix<T>> ListGenerators;
    auto FuncInsert=[&](MyMatrix<T> const& TheMat) -> void {
#ifdef DEBUG_APPROXIMATE_MODELS
      if (!IsIntegralMatrix(TheMat)) {
        std::cerr << "The matrix TheMat is not integral\n";
        throw TerminalException{1};
      }
      MyMatrix<T> eProd = TheMat * Qmat * TheMat.transpose();
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
    for (auto & eGen : GeneratorsSL2Z<T>()) {
      FuncInsert(ExpandGenerator(GetLeftMultiplication(eGen)));
      FuncInsert(ExpandGenerator(GetRightMultiplication(eGen)));
    }
    // Now looking at the isotropic vectors, generating the Eichler transvections
    for (int i=0; i<4; i++) {
      MyVector<T> eVect = ZeroVector<T>(n);
      eVect(i) = 1;
      MyVector<T> eProd = Qmat * eVect;
      MyMatrix<T> eProd_M(1,n);
      for (int i=0; i<n; i++) {
        eProd_M(0,i) = eProd(i);
      }
      MyMatrix<T> BasisOrth = NullspaceIntMat(eProd_M);
      for (int i=0; i<BasisOrth.rows(); i++) {
        MyVector<T> eOrth = GetMatrixRow(BasisOrth, i);
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
  std::function<std::vector<MyVector<T>>(T)> EnumerateVectorOverDiscriminant = [=](T const& X) -> std::vector<MyVector<T>> {
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
    MyMatrix<T> const& Qmat = shr_ptr->Qmat;
    MyMatrix<T> const& Gmat = shr_ptr->Gmat;
    int n = Qmat.rows();
    T two(2);
    if (ResInt(X, two) == T(1)) {
      return {};
    }
    T Xdiv2 = X / 2;
    // The two hyperbolic planes are unimodular so for the discriminant, we only need
    // to consider the Gmat part.
    MyMatrix<T> Ginv = Inverse(Gmat);
    std::vector<MyVector<T>> ListSolution;
    for (auto & eClass1 : ListClasses) {
      // The vector eClass1 is defined modulo an element of Gmat Z^n
      // Thus eClass2 is defined modulo an element of Z^n. This is an elements of L*
      // Thus eClass3 is defined modulo an element of d Z^n
      // So, we can write v = eClass3 + d w
      // which gets us A[v] = A[eClass3] + 2d w^T Gmat eClass3 + d^2 A[w]
      // So, the problem is actually linear.
      MyVector<T> eClass2 = Ginv * eClass1;
      T DivD_frac = RemoveFractionVectorPlusCoeff(eClass2).TheMult;
      T DivD = GetDenominator(DivD_frac);
      T DivD_sqr = DivD * DivD;
      MyVector<T> eClass3 = DivD * eClass2;
      T X_res1 = ResInt(Xdiv2, DivD);
      T X_res2 = ResInt(Xdiv2, DivD_sqr);
      T Aclass3_norm = EvaluationQuadForm(Gmat, eClass3) / 2;
      T Aclass3_res1 = ResInt(Aclass3_norm, DivD);
      T Aclass3_res2 = ResInt(Aclass3_norm, DivD_sqr);
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
          MyVector<T> eSolution = ZeroVector<T>(n);
          eSolution(2) = DivD;
          eSolution(3) = DivD * u;
          for (int u=0; u<n-4; u++) {
            eSolution(u+4) = eClass4(u);
          }
#ifdef DEBUG_APPROXIMATE_MODELS
          T eNorm = EvaluationQuadForm(Qmat, eSolution);
          if (eNorm != X) {
            std::cerr << "eNorm / X is inconsistent\n";
            throw TerminalException{1};
          }
#endif
          ListSolution.push_back(eSolution);
        }
      }
    }
    return ListSolution;
  };
  std::function<std::vector<MyVector<T>>(T)> GetCoveringOrbitRepresentatives = [=](T const& X) -> std::vector<MyVector<T>> {
    if (X == 0) {
      return EnumerateVectorOverDiscriminant(0);
    }
    std::vector<T> ListDiv = GetSquareDivisors(X);
    std::vector<MyVector<T>> ListSolution;
    for (auto & eDiv : ListDiv) {
      T Xtarget = X / (eDiv * eDiv);
      for (auto & eSol : EnumerateVectorOverDiscriminant(Xtarget)) {
        MyVector<T> eV = eDiv * eSol;
        ListSolution.push_back(eV);
      }
    }
    return ListSolution;
  };
  std::function<std::optional<MyVector<T>>(T const&)> GetOneOrbitRepresentative = [=](T const& X) -> std::optional<MyVector<T>> {
    T two(2);
    if (ResInt(X, two) == T(1)) {
      return {};
    }
    T Xdiv2 = X / 2;
    MyMatrix<T> const& Qmat = shr_ptr->Gmat;
    int n = Qmat.rows();
    MyVector<T> eSol = ZeroVector<T>(n);
    eSol(0) = 1;
    eSol(1) = Xdiv2;
    return eSol;
  };
  return {GetApproximateGroup, SetListClassesOrbitwise, GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative};
}



// clang-format off
#endif  // SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_
// clang-format on
