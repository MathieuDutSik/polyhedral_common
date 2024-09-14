// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_
#define SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_

// clang-format off
#include "Isotropic.h"
#include "GRP_GroupFct.h"
#include "MAT_MatrixInt.h"
#include "PositiveNegative.h"
#include "LatticeStabEquiCan.h"
#include "IndefiniteFormFundamental.h"
#include "MatrixGroup.h"
#include <memory>
#include <map>
#include <set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_APPROXIMATE_MODELS
#endif

#ifdef TIMINGS
#define TIMINGS_APPROXIMATE_MODELS
#endif

// This is the construction of the model from the Eichler criterion.
// This is used for the combined algorithms.
//
// The reference that we use is
// Scattone F. On the compactification of moduli spaces for algebraic K3
// surfaces Specific references: (1) Section 3.7 (page 41): the construction of
// the Eichler transvection
//     That gives us the approximate group.
// (2) Proposition 3.7.3 should give us the list of possible vectors.
//     That gives us the covering orbit representatives.

template <typename T, typename Tint> struct ApproximateModel {
  std::function<std::vector<MyMatrix<Tint>>(std::ostream &)>
      GetApproximateGroup;
  std::function<void(std::vector<MyMatrix<Tint>>, std::ostream &)>
      SetListClassesOrbitwise;
  std::function<std::vector<MyVector<Tint>>(T, std::ostream &)>
      GetCoveringOrbitRepresentatives;
  std::function<std::optional<MyVector<Tint>>(T, std::ostream &)>
      GetOneOrbitRepresentative;
};

// The Eichler transvection are defined for an even lattice (so integral and
// even norms of vectors. For such a lattice L* is defined as Q^{-1} Z^n The
// formula describing them is E_{f,x}(y) = y + (y,x) f - (x,x) / 2 (y,f) f -
// (y,f) x The Eichler transvection preserve the lattice L because (x,x) / 2 in
// Z, (y,x) in Z and (y,f) in Z. Also if y in L* since the above is still true,
// we have that E_{f,x} preserve L* and the image of the transformation in L*/L
// is the identity.
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
template <typename T, typename Tint>
MyMatrix<Tint> INDEF_FORM_Eichler_Transvection(MyMatrix<T> const &Qmat,
                                               MyVector<Tint> const &f,
                                               MyVector<Tint> const &x) {
  if (!INDEF_FORM_IsEven(Qmat)) {
    std::cerr << "The lattice Qmat should be even in order to define the "
                 "Eichler transvection\n";
    throw TerminalException{1};
  }
  int n = Qmat.rows();
  MyVector<T> x_T = UniversalVectorConversion<T, Tint>(x);
  MyVector<T> f_T = UniversalVectorConversion<T, Tint>(f);
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
  MyMatrix<Tint> RetMat(n, n);
  for (int u = 0; u < n; u++) {
    MyVector<Tint> eImg = ZeroVector<Tint>(n);
    MyVector<Tint> y = ZeroVector<Tint>(n);
    y(u) = 1;
    eImg += y;
    //
    T scal1 = ScalarProductQuadForm(Qmat, y, x);
    Tint scal1_tint = UniversalScalarConversion<Tint, T>(scal1);
    eImg += scal1_tint * f;
    T scal2 = ScalarProductQuadForm(Qmat, y, f);
    T scal3 = (xNorm / 2) * scal2;
    Tint scal2_tint = UniversalScalarConversion<Tint, T>(scal2);
    Tint scal3_tint = UniversalScalarConversion<Tint, T>(scal3);
    eImg -= scal3_tint * f;
    eImg -= scal2_tint * x;
    AssignMatrixRow(RetMat, u, eImg);
  }
#ifdef DEBUG_APPROXIMATE_MODELS
  MyMatrix<T> RetMat_T = UniversalMatrixConversion<T, Tint>(RetMat);
  MyMatrix<T> prod = RetMat_T * Qmat * RetMat_T.transpose();
  if (prod != Qmat) {
    std::cerr << "RetMat should prerve the Qmat\n";
    throw TerminalException{1};
  }
#endif
  return RetMat;
}

template <typename T, typename Tint>
std::optional<MyVector<Tint>> INDEF_FindIsotropic(MyMatrix<T> const &M,
                                                  std::ostream &os) {
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
    MyVector<Tint> eV = UniversalVectorConversion<Tint, T>(eV_T);
#ifdef DEBUG_APPROXIMATE_MODELS
    T sum = EvaluationQuadForm<T, Tint>(M, eV);
    if (sum != 0) {
      std::cerr << "MODEL: INDEF_FindIsotropic, error in the isotrop by rank\n";
      throw TerminalException{1};
    }
#endif
    return eV;
  }
  std::optional<MyVector<T>> opt = FindIsotropic(M, os);
  if (opt) {
    MyVector<T> const &eV = *opt;
    MyVector<T> fV = RemoveFractionVector(eV);
    MyVector<Tint> gV = UniversalVectorConversion<Tint, T>(fV);
#ifdef DEBUG_APPROXIMATE_MODELS
    T sum = EvaluationQuadForm<T, Tint>(M, gV);
    if (sum != 0) {
      std::cerr << "MODEL: INDEF_FindIsotropic, error in the isotrop by "
                   "FindIsotropic\n";
      throw TerminalException{1};
    }
#endif
    return gV;
  } else {
    return {};
  }
}

template <typename Tint> std::vector<MyMatrix<Tint>> GeneratorsSL2Z() {
  MyMatrix<Tint> eGenS = ZeroMatrix<Tint>(2, 2);
  eGenS(0, 1) = -1;
  eGenS(1, 0) = 1;
  MyMatrix<Tint> eGenT = ZeroMatrix<Tint>(2, 2);
  eGenT(0, 0) = 1;
  eGenT(0, 1) = 1;
  eGenT(1, 1) = 1;
  return {eGenS, eGenT};
}

template <typename T, typename Tint> struct InternalEichler {
  MyMatrix<T> Qmat;
  MyMatrix<T> Gmat;
  std::vector<MyVector<Tint>> ListClasses;
  std::vector<MyMatrix<Tint>> GRPstart;
};

template <typename T> std::vector<T> GetSquareDivisors(T const &X) {
  std::map<T, size_t> map = FactorsIntMap(T_abs(X));
  std::vector<std::vector<T>> ListListProd;
  int64_t two(2);
  for (auto &kv : map) {
    T const &p = kv.first;
    int64_t mult = kv.second;
    int64_t mult_div = QuoInt(mult, two);
    std::vector<T> ListProd;
    T val(1);
    for (int64_t u = 0; u <= mult_div; u++) {
      ListProd.push_back(val);
      val *= p;
    }
    ListListProd.push_back(ListProd);
  }
  std::vector<T> ListSquareDiv{T(1)};
  for (auto &ListProd : ListListProd) {
    std::vector<T> V;
    for (auto &x : ListSquareDiv) {
      for (auto &y : ListProd) {
        T z = x * y;
        V.push_back(z);
      }
    }
    ListSquareDiv = V;
  }
  return ListSquareDiv;
}

template <typename T> bool is_eichler_canonical(MyMatrix<T> const &Qmat) {
  int n = Qmat.rows();
  if (n < 4) {
    return false;
  }
  std::vector<int> LIdx{1, 0, 3, 2};
  for (int iRow = 0; iRow < 4; iRow++) {
    int iColCrit = LIdx[iRow];
    for (int iCol = 0; iCol < 4; iCol++) {
      if (iCol == iColCrit) {
        if (Qmat(iRow, iCol) != 1) {
          return false;
        }
      } else {
        if (Qmat(iRow, iCol) != 0) {
          return false;
        }
      }
    }
  }
  for (int iRow = 0; iRow < 4; iRow++) {
    for (int iCol = 4; iCol < n; iCol++) {
      if (Qmat(iRow, iCol) != 0) {
        return false;
      }
    }
  }
  if (!INDEF_FORM_IsEven(Qmat)) {
    return false;
  }
  return true;
}

/*
  Dualities and action of the automorphism group on the disciminant.
  The Gram matrix G realize a scalar product <x, y> = x G y.
  ---
  The dual is L^*.
  If L is integral (that is G is integral) then L is a subset of L^*.
  The quotient L^* / L is a finite group has an order of the size det(G).
  The representatives are v in Z^n quotiented by the relation
  x equiv y <=> x - y = z G
  ---
  The automorphism group acts like
  P G P^T = G.
  x P^T - y P^T = z G P^T
                = z' G with z' = z P^{-1}
  That would mean that the action is the direct one.

 */

/*
  Based on paragraph 10 of Eichler book
  Quadratische Formen und orthogonale gruppen
  Also used the Scattone thesis as basis.

  Note: There are other works based on Eichler work that are analogous:
  ---Wall CTC, On The Orthogonal Groups of Unimodular Quadratic Forms
  ---James D.G. Indefinite Quadratic Forms of Determinant pm2p
  They also require the same hypothesis of having two hyperplanes.
*/
template <typename T, typename Tint, typename Tgroup>
ApproximateModel<T, Tint>
INDEF_FORM_EichlerCriterion_TwoHyperplanesEven(MyMatrix<T> const &Qmat) {
  int n = Qmat.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
  if (!is_eichler_canonical(Qmat)) {
    std::cerr << "The matrix should be Eichler canonical form\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> Gmat(n - 4, n - 4);
  for (int i = 4; i < n; i++) {
    for (int j = 4; j < n; j++) {
      Gmat(i - 4, j - 4) = Qmat(i, j);
    }
  }
  // ComputeClasses
  std::vector<MyVector<Tint>> ListClasses =
      ComputeTranslationClasses<T, Tint>(Gmat);
  std::vector<MyMatrix<Tint>> GRPstart;
  InternalEichler<T, Tint> ie{Qmat, Gmat, ListClasses, GRPstart};
  std::shared_ptr<InternalEichler<T, Tint>> shr_ptr =
      std::make_shared<InternalEichler<T, Tint>>(ie);
  std::function<void(std::vector<MyMatrix<Tint>> const &, std::ostream &)>
      SetListClassesOrbitwise = [=](std::vector<MyMatrix<Tint>> const &GRPmatr,
                                    [[maybe_unused]] std::ostream &os) -> void {
#ifdef TIMINGS_APPROXIMATE_MODELS
    MicrosecondTime time;
#endif
    shr_ptr->GRPstart = GRPmatr;
    MyMatrix<T> const &Qmat = shr_ptr->Qmat;
    int n = Qmat.rows();
    int n_classes = ListClasses.size();
    std::vector<MyVector<Tint>> ListClassesExt;
    for (auto &eV : ListClasses) {
      MyVector<Tint> NewV = ZeroVector<Tint>(n);
      for (int i = 0; i < n - 4; i++) {
        NewV(i + 4) = eV(i);
      }
      ListClassesExt.emplace_back(std::move(NewV));
    }
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: SetListClassesOrbitwise, |ListClassesExt|=" << n_classes
       << "\n";
#endif
    auto GetPosition = [&](MyVector<Tint> const &eVect) -> size_t {
      for (size_t i = 0; i < ListClassesExt.size(); i++) {
        MyVector<Tint> eDiff = ListClassesExt[i] - eVect;
        MyVector<T> eDiff_T = UniversalVectorConversion<T, Tint>(eDiff);
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
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: SetListClassesOrbitwise, |GRPmatr|=" << GRPmatr.size()
       << "\n";
#endif
    for (auto &eMatrGen : GRPmatr) {
      std::vector<Tidx> eList;
#ifdef DEBUG_APPROXIMATE_MODELS_DISABLE
      std::vector<int> status(n_classes, 0);
#endif
      for (auto &eClassExt : ListClassesExt) {
#ifdef DEBUG_APPROXIMATE_MODELS
        MyMatrix<T> eMatrGen_T = UniversalMatrixConversion<T, Tint>(eMatrGen);
        MyMatrix<T> prod = eMatrGen_T * Qmat * eMatrGen_T.transpose();
        if (prod != Qmat) {
          std::cerr << "The Matrix is not preserving the quadratic form\n";
          throw TerminalException{1};
        }
#endif
        MyVector<Tint> x_eM = eMatrGen * eClassExt;
        Tidx pos = GetPosition(x_eM);
#ifdef DEBUG_APPROXIMATE_MODELS_DISABLE
        os << "MODEL: SetListClassesOrbitwise, pos=" << static_cast<int>(pos)
           << "\n";
        int &val = status[pos];
        if (val == 1) {
          std::cerr << "The value has already been attained\n";
          throw TerminalException{1};
        }
        val = 1;
#endif
        eList.push_back(pos);
      }
      Telt ePerm(eList);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: SetListClassesOrbitwise, ePerm=" << ePerm << "\n";
#endif
      ListPermGens.push_back(ePerm);
    }
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: SetListClassesOrbitwise, we have ListPermGen\n";
#endif
    Tgroup GRPperm(ListPermGens, n_classes);
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: SetListClassesOrbitwise, |GRPperm|=" << GRPperm.size()
       << "\n";
#endif
    std::vector<size_t> ListRepr =
        DecomposeOrbitPoint_FullRepr<Tgroup>(GRPperm);
    std::vector<MyVector<Tint>> ListClasses;
    for (auto &iRepr : ListRepr) {
      MyVector<Tint> eV = shr_ptr->ListClasses[iRepr];
      ListClasses.push_back(eV);
    }
    shr_ptr->ListClasses = ListClasses;
#ifdef TIMINGS_APPROXIMATE_MODELS
    os << "MODEL: |SetListClassesOrbitwise|=" << time << "\n";
#endif
  };
  std::function<std::vector<MyMatrix<Tint>>(std::ostream &)>
      GetApproximateGroup = [=]([[maybe_unused]] std::ostream &os)
      -> std::vector<MyMatrix<Tint>> {
#ifdef TIMINGS_APPROXIMATE_MODELS
    MicrosecondTime time;
#endif
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of GetApproximateGroup : 1\n";
#endif
    MyMatrix<T> const &Qmat = shr_ptr->Qmat;
    int n = Qmat.rows();
    // The quadratic form 2 x1 x2 + 2 x3 x4
    // We formulate it as (1/2) det U(x1,x2,x3,x4)
    // with U(x1,x2,x3,x4) = [ [x1 , -x4] , [x3 , x2] ]
    // If we multiply by a matrix in SL(2,Z) on the left and right
    auto GetLeftMultiplication =
        [&](MyMatrix<Tint> const &P) -> MyMatrix<Tint> {
      // We compute the product [[a, b],[c,d]] U(x1,x2,x3,x4)
      // This gives us
      // [[a x1 + b x3 , -a x4 + b x2] , [ c x1 + d x3 , -c x4 + d x2]]
      // So Y = XA,
      // Y1 =  a x1 + b x3
      // Y2 = -c x4 + d x2
      // Y3 =  c x1 + d x3
      // Y4 =  a x4 - b x2
      Tint a = P(0, 0);
      Tint b = P(0, 1);
      Tint c = P(1, 0);
      Tint d = P(1, 1);
      MyMatrix<Tint> RetMat = ZeroMatrix<Tint>(4, 4);
      RetMat(0, 0) = a;
      RetMat(0, 2) = c;
      RetMat(1, 1) = d;
      RetMat(1, 3) = -b;
      RetMat(2, 0) = b;
      RetMat(2, 2) = d;
      RetMat(3, 1) = -c;
      RetMat(3, 3) = a;
      return RetMat;
    };
    auto GetRightMultiplication =
        [&](MyMatrix<Tint> const &P) -> MyMatrix<Tint> {
      // We compute the product [ [x1 , -x4] , [x3 , x2] ]      [[a, b],[c,d]]
      // This gives us
      // [[a x1 - c x4 , b x1 - d x4] , [a x3 + c x2 , b x3 + d x2]]
      // So Y = XA,
      // Y1 = a x1 - c x4
      // Y2 = b x3 + d x2
      // Y3 = a x3 + c x2
      // Y4 = -b x1 + d x4
      Tint a = P(0, 0);
      Tint b = P(0, 1);
      Tint c = P(1, 0);
      Tint d = P(1, 1);
      MyMatrix<Tint> RetMat = ZeroMatrix<Tint>(4, 4);
      RetMat(0, 0) = a;
      RetMat(0, 3) = -b;
      RetMat(1, 1) = d;
      RetMat(1, 2) = c;
      RetMat(2, 1) = b;
      RetMat(2, 2) = a;
      RetMat(3, 0) = -c;
      RetMat(3, 3) = d;
      return RetMat;
    };
    auto ExpandGenerator = [&](MyMatrix<Tint> const &eGen) -> MyMatrix<Tint> {
      MyMatrix<Tint> BigMat = IdentityMat<Tint>(n);
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          BigMat(i, j) = eGen(i, j);
        }
      }
      return BigMat;
    };
    std::vector<MyMatrix<Tint>> ListGenerators;
    auto FuncInsert = [&](MyMatrix<Tint> const &TheMat) -> void {
#ifdef DEBUG_APPROXIMATE_MODELS
      MyMatrix<T> TheMat_T = UniversalMatrixConversion<T, Tint>(TheMat);
      MyMatrix<T> eProd = TheMat_T * Qmat * TheMat_T.transpose();
      if (eProd != Qmat) {
        std::cerr << "The matrix is not preserving the quadratic form\n";
        throw TerminalException{1};
      }
#endif
      ListGenerators.push_back(TheMat);
    };
    for (auto &eGen : shr_ptr->GRPstart) {
      FuncInsert(eGen);
    }
    for (auto &eGen : GeneratorsSL2Z<Tint>()) {
      FuncInsert(ExpandGenerator(GetLeftMultiplication(eGen)));
      FuncInsert(ExpandGenerator(GetRightMultiplication(eGen)));
    }
    // Now looking at the isotropic vectors, generating the Eichler
    // transvections
    for (int i = 0; i < 4; i++) {
      MyVector<Tint> eVect = ZeroVector<Tint>(n);
      eVect(i) = 1;
      MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect);
      MyVector<T> eProd = Qmat * eVect_T;
      MyMatrix<T> eProd_M(1, n);
      for (int i = 0; i < n; i++) {
        eProd_M(0, i) = eProd(i);
      }
      MyMatrix<T> BasisOrth = NullspaceIntMat(eProd_M);
      for (int i = 0; i < BasisOrth.rows(); i++) {
        MyVector<T> eOrth_T = GetMatrixRow(BasisOrth, i);
        MyVector<Tint> eOrth = UniversalVectorConversion<Tint, T>(eOrth_T);
        FuncInsert(INDEF_FORM_Eichler_Transvection(Qmat, eVect, eOrth));
      }
    }
#ifdef TIMINGS_APPROXIMATE_MODELS
    os << "MODEL: |GetApproximateGroup|=" << time << "\n";
#endif
    return ListGenerators;
  };
  // Relevant notations (see [S, v]):
  // (1) G_L = L^* / L
  // (2) div(v) is the integer d such that (v, L) = d Z
  // (3) v^* is the class of v / div(v) in G_L
  //
  // The theorem says that among primitive vectors for a fixed norm
  // we have that v equivalent to w under the resticted isometry group
  // if and only if v^* = w^*.
  //
  // Therefore, the algorithm should proceed in the following way:
  // (1) Iterating over the elements of G_L
  // (2) Identifying the Div

  // Notion of divisor
  // For an element y in a lattice L, define div(y) the integer such that
  // (L,y) = div(y) Z.
  // For each v in L* we can find smallest d such that w = dv in L.
  // And we then have div(w) = d
  std::function<std::vector<MyVector<Tint>>(T const &, std::ostream &)>
      EnumerateVectorOverDiscriminant =
          [=](T const &X, [[maybe_unused]] std::ostream &os)
      -> std::vector<MyVector<Tint>> {
#ifdef TIMINGS_APPROXIMATE_MODELS
    MicrosecondTime time;
#endif
    // For each vector v in M* / M
    // We want (x, v) divided by d.
    //
    // We should have e_I.x = d alpha_I
    // We should have (x,x) = X
    // For the first coordinates, we can write them as (0 , 0 , d , du)
    // By changing the value of u, we can shift the value of u which shifts the
    // value by 2 d^2, + or -. The lattice is even so the 2 is not a problem.
    //
    // So, if we have some X, we want to find u\in Z and v in Z^{n-4} such that
    // X = 2d^2 u + A[v] with A the Gram matrix of M. v must actually be of the
    // form v = v0 + d A^{-1} h. We are only interested on the value modulo 2d^2
    // of A on those vectors. So, it is a finite problem. But it is a hard one
    // as the search space might be very large.
    //
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of EnumerateVectorOverDiscriminant X=" << X << "\n";
#endif
    MyMatrix<T> const &Qmat = shr_ptr->Qmat;
    MyMatrix<T> const &Gmat = shr_ptr->Gmat;
    std::vector<MyVector<Tint>> const &ListClasses = shr_ptr->ListClasses;
    int n = Qmat.rows();
    T two(2);
    if (ResInt(X, two) == T(1)) {
      return {};
    }
    T Xdiv2 = X / 2;
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: EnumerateVectorOverDiscriminant Xdiv2=" << Xdiv2 << "\n";
#endif
    // The two hyperbolic planes are unimodular so for the discriminant, we only
    // need to consider the Gmat part.
    std::vector<MyVector<Tint>> ListSolution;
    if (n > 4) {
      MyMatrix<T> Ginv = Inverse(Gmat);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: EnumerateVectorOverDiscriminant |ListClasses|="
         << ListClasses.size() << "\n";
#endif
      for (auto &eClass1 : ListClasses) {
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: ----------------------------------------\n";
        os << "MODEL: EnumerateVectorOverDiscriminant eClass1="
           << StringVectorGAP(eClass1) << "\n";
#endif
        // The vector eClass1 is defined modulo an element of Gmat Z^n
        // Thus eClass2 is defined modulo an element of Z^n. This is an elements
        // of L* Then d = div(v). Thus eClass3 is defined modulo an element of d
        // Z^n So, we can write v = eClass3 + d w which gets us X = A[v] =
        // A[eClass3] + 2d w^T Gmat eClass3 + d^2 A[w] So, the problem is
        // actually linear.
        MyVector<T> eClass1_T = UniversalVectorConversion<T, Tint>(eClass1);
        MyVector<T> eClass2 = Ginv * eClass1_T;
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: EnumerateVectorOverDiscriminant eClass2="
           << StringVectorGAP(eClass2) << "\n";
#endif
        // DivV_frac is the fraction such that DivV_frac eClass2 is a primitive
        // vector.
        T DivV_frac = RemoveFractionVectorPlusCoeff(eClass2).TheMult;
        T DivV = GetNumerator(DivV_frac);
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: EnumerateVectorOverDiscriminant DivV_frac=" << DivV_frac
           << " DivV=" << DivV << "\n";
#endif
        T DivV_sqr = DivV * DivV;
        MyVector<T> eClass3 = DivV * eClass2;
        T X_res1 = ResInt(Xdiv2, DivV);
        T X_res2 = ResInt(Xdiv2, DivV_sqr);
        T Aclass3_norm = EvaluationQuadForm(Gmat, eClass3) / 2;
        T Aclass3_res1 = ResInt(Aclass3_norm, DivV);
        T Aclass3_res2 = ResInt(Aclass3_norm, DivV_sqr);
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: EnumerateVectorOverDiscriminant Aclass3_norm="
           << Aclass3_norm << "\n";
        os << "MODEL: EnumerateVectorOverDiscriminant Aclass3_res1="
           << Aclass3_res1 << " X_res1=" << X_res1 << "\n";
#endif
        if (Aclass3_res1 == X_res1) { // if not we cannot find a solution
          // and so
          // (X - A[eClass3])/2 = d w^T Gmat eClass3 + ....
          T diff = (X_res2 - Aclass3_res2) / DivV;
          MyVector<T> eProd = Gmat * eClass3;
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: EnumerateVectorOverDiscriminant diff=" << diff << "\n";
          os << "MODEL: EnumerateVectorOverDiscriminant eProd="
             << StringVectorGAP(eProd) << "\n";
#endif
          GCD_dot<T> eRec = ComputeGcdDot(eProd);
          auto iife_quot = [&]() -> T {
            if (eRec.gcd == 0) {
              return T(0);
            } else {
              return diff / eRec.gcd;
            }
          };
          T quot = iife_quot();
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: EnumerateVectorOverDiscriminant quot=" << quot << "\n";
          os << "MODEL: EnumerateVectorOverDiscriminant eRec.V="
             << StringVectorGAP(eRec.V) << "\n";
#endif
          if (IsInteger(quot)) { // If not then the equation has no solution
            MyVector<T> w = quot * eRec.V;
            MyVector<T> eClass4 = eClass3 + DivV * w;
            T Aclass4_norm = EvaluationQuadForm(Gmat, eClass4) / 2;
            T Aclass4_res2 = ResInt(Aclass4_norm, DivV_sqr);
#ifdef DEBUG_APPROXIMATE_MODELS
            os << "MODEL: EnumerateVectorOverDiscriminant w="
               << StringVectorGAP(w) << "\n";
            os << "MODEL: EnumerateVectorOverDiscriminant eClass4="
               << StringVectorGAP(eClass4) << "\n";
            os << "MODEL: EnumerateVectorOverDiscriminant Aclass4_res2="
               << Aclass4_res2 << " X_res2=" << X_res2 << "\n";
            if (Aclass4_res2 != X_res2) {
              std::cerr << "A bug to resolve\n";
              throw TerminalException{1};
            }
#endif
            T u = (Xdiv2 - Aclass4_norm) / DivV_sqr;
#ifdef DEBUG_APPROXIMATE_MODELS
            os << "MODEL: EnumerateVectorOverDiscriminant u=" << u << "\n";
#endif
            MyVector<T> eSolution_T = ZeroVector<T>(n);
            eSolution_T(2) = DivV;
            eSolution_T(3) = DivV * u;
            for (int u = 0; u < n - 4; u++) {
              eSolution_T(u + 4) = eClass4(u);
            }
#ifdef DEBUG_APPROXIMATE_MODELS
            T eNorm = EvaluationQuadForm(Qmat, eSolution_T);
            if (eNorm != X) {
              std::cerr << "eNorm / X is inconsistent\n";
              throw TerminalException{1};
            }
            os << "MODEL: EnumerateVectorOverDiscriminant DivV=" << DivV
               << " u=" << u << "\n";
            os << "MODEL: EnumerateVectorOverDiscriminant eSolution_T="
               << StringVectorGAP(eSolution_T) << "\n";
#endif
            MyVector<Tint> eSolution =
                UniversalVectorConversion<Tint, T>(eSolution_T);
            ListSolution.push_back(eSolution);
          }
        }
      }
    } else {
      MyVector<T> eSolution_T = ZeroVector<T>(4);
      eSolution_T(2) = 1;
      eSolution_T(3) = Xdiv2;
      MyVector<Tint> eSolution =
          UniversalVectorConversion<Tint, T>(eSolution_T);
      ListSolution.push_back(eSolution);
    }
#ifdef TIMINGS_APPROXIMATE_MODELS
    os << "MODEL: |EnumerateVectorOverDiscriminant|=" << time << "\n";
#endif
    return ListSolution;
  };
  std::function<std::vector<MyVector<Tint>>(T const &, std::ostream &)>
      GetCoveringOrbitRepresentatives =
          [=](T const &X, [[maybe_unused]] std::ostream &os)
      -> std::vector<MyVector<Tint>> {
#ifdef TIMINGS_APPROXIMATE_MODELS
    MicrosecondTime time;
#endif
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of GetCoveringOrbitRepresentatives 1: X=" << X
       << "\n";
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
    for (auto &eDiv_T : ListDiv) {
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: eDiv_T=" << eDiv_T << "\n";
#endif
      Tint eDiv = UniversalScalarConversion<Tint, T>(eDiv_T);
      T Xtarget = X / (eDiv_T * eDiv_T);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Xtarget=" << Xtarget << "\n";
#endif
      for (auto &eSol : EnumerateVectorOverDiscriminant(Xtarget, os)) {
        MyVector<Tint> eV = eDiv * eSol;
        ListSolution.push_back(eV);
      }
    }
#ifdef TIMINGS_APPROXIMATE_MODELS
    os << "MODEL: |GetCoveringOrbitRepresentatives|=" << time << "\n";
#endif
    return ListSolution;
  };
  std::function<std::optional<MyVector<Tint>>(T const &, std::ostream &)>
      GetOneOrbitRepresentative =
          [=](T const &X, [[maybe_unused]] std::ostream &os)
      -> std::optional<MyVector<Tint>> {
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
    MyMatrix<T> const &Qmat = shr_ptr->Qmat;
    int n = Qmat.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: n=" << n << "\n";
#endif
    MyVector<Tint> eSol = ZeroVector<Tint>(n);
    eSol(0) = 1;
    eSol(1) = Xdiv2;
    return eSol;
  };
  return {GetApproximateGroup, SetListClassesOrbitwise,
          GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative};
}

/*
  We compute some isometries of the lattice coming from the blocks.
  The only constraint is that the matrices being returned are isometries.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>> GetEasyIsometries(MyMatrix<T> const &Qmat,
                                              std::ostream &os) {
#ifdef TIMINGS_APPROXIMATE_MODELS
  MicrosecondTime time;
#endif
  std::vector<MyMatrix<Tint>> ListGenerators;
  auto f_insert = [&](MyMatrix<Tint> const &P) -> void {
#ifdef DEBUG_APPROXIMATE_MODELS
    MyMatrix<T> P_T = UniversalMatrixConversion<T, Tint>(P);
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
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if (Qmat(i, j) != 0) {
        eG.AddAdjacent(i, j);
        eG.AddAdjacent(j, i);
      }
    }
  }
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(eG);
  std::vector<MyMatrix<T>> ListQ;
  std::vector<size_t> ListPosIdx;
  for (size_t iConn = 0; iConn < LConn.size(); iConn++) {
    std::vector<size_t> const &eConn = LConn[iConn];
    size_t dim = eConn.size();
    MyMatrix<T> eQ(dim, dim);
    for (size_t i = 0; i < dim; i++) {
      for (size_t j = 0; j < dim; j++) {
        eQ(i, j) = Qmat(eConn[i], eConn[j]);
      }
    }
    ListQ.push_back(eQ);
    bool test = INDEF_FORM_IsPosNeg(eQ);
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: iConn=" << iConn << " dim=" << dim << " test=" << test
       << "\n";
    os << "MODEL: eQ=\n";
    WriteMatrix(os, eQ);
#endif
    if (test) {
      ListPosIdx.push_back(iConn);
    }
  }
  for (auto &iConn : ListPosIdx) {
    MyMatrix<T> const &eQ = ListQ[iConn];
    std::vector<size_t> const &eConn = LConn[iConn];
    size_t dim = eConn.size();
    std::vector<MyMatrix<Tint>> LGen =
        INDEF_FORM_AutomorphismGroup_PosNeg<T, Tint>(eQ, os);
    for (auto &eGen : LGen) {
      MyMatrix<Tint> eBigGen = IdentityMat<Tint>(n);
      for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
          eBigGen(eConn[i], eConn[j]) = eGen(i, j);
        }
      }
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: Before f_insert, case 1\n";
#endif
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
  for (size_t iConnRed = 0; iConnRed < nbConnRed; iConnRed++) {
    size_t iConn = ListPosIdx[iConnRed];
    MyMatrix<T> const &eQ = ListQ[iConn];
    std::vector<size_t> const &eConn = LConn[iConn];
    for (size_t jConnRed = iConnRed + 1; jConnRed < nbConnRed; jConnRed++) {
      size_t jConn = ListPosIdx[jConnRed];
      MyMatrix<T> const &fQ = ListQ[jConn];
      std::vector<size_t> const &fConn = LConn[jConn];
      size_t dim = eConn.size();
      std::optional<MyMatrix<Tint>> opt =
          INDEF_FORM_TestEquivalence_PosNeg<T, Tint>(eQ, fQ, os);
      if (opt) {
        MyMatrix<Tint> const &P = *opt;
        MyMatrix<Tint> Pinv = Inverse(P);
        MyMatrix<Tint> BigP = IdentityMat<Tint>(n);
        for (size_t i = 0; i < dim; i++) {
          for (size_t j = 0; j < dim; j++) {
            BigP(eConn[i], eConn[j]) = 0;
            BigP(fConn[i], fConn[j]) = 0;
            BigP(eConn[i], fConn[j]) = Pinv(i, j);
            BigP(fConn[i], eConn[j]) = P(i, j);
          }
        }
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: Before f_insert, case 2\n";
#endif
        f_insert(BigP);
      }
    }
  }
#ifdef TIMINGS_APPROXIMATE_MODELS
  os << "MODEL: |GetEasyIsometries|=" << time << "\n";
#endif
  return ListGenerators;
}

template <typename T, typename Tint> struct ResultHyperbolicPlane {
  MyMatrix<Tint> Basis;
  T scal;
};

template <typename T, typename Tint>
ResultHyperbolicPlane<T, Tint> get_result_hyperbolic_plane(
    MyMatrix<T> const &M,
    std::pair<MyVector<Tint>, MyVector<Tint>> const &pair) {
  MyMatrix<Tint> Basis = MatrixFromPairVector(pair);
  MyVector<T> V1 = UniversalVectorConversion<T, Tint>(pair.first);
  MyVector<T> V2 = UniversalVectorConversion<T, Tint>(pair.second);
  MyVector<T> prod = M * V1;
  T scal = prod.dot(V2);
  return {Basis, scal};
}

/*
  Ideally, given an isotropic vector v, and c, we want to find a vector w such
  that v.w = c and w is isotrop.
  ----
  Compute V = v^{perp}. We want to find the vectors in w in V with phi(w) = c
  for some linear form v. We can work in the quotient W = V / v and the function
  phi lifts to W. If the initial space has signature (p,q), then the resulting
  space W has signature (p-1,q-1).
 */
template <typename T, typename Tint>
ResultHyperbolicPlane<T, Tint> GetHyperbolicPlane(MyMatrix<T> const &Qmat,
                                                  std::ostream &os) {
#ifdef TIMINGS_APPROXIMATE_MODELS
  MicrosecondTime time;
#endif
  int n = Qmat.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: Beginning of GetHyperbolicPlane, Qmat=\n";
  WriteMatrix(os, Qmat);
#endif
  std::set<MyVector<Tint>> SetVect;
  MyMatrix<Tint> ThePerturb = IdentityMat<Tint>(n);

  int n_iter = 100;
  auto f_insert_vect = [&]() -> MyVector<Tint> {
    int iter = 0;
    while (true) {
      MyMatrix<Tint> ePerturb = GetRandomMatrixPerturbation<Tint>(n);
      ThePerturb = ThePerturb * ePerturb;
      MyMatrix<T> ThePerturb_T = UniversalMatrixConversion<T, Tint>(ThePerturb);
      MyMatrix<T> M = ThePerturb_T * Qmat * ThePerturb_T.transpose();
      std::optional<MyVector<Tint>> opt = INDEF_FindIsotropic<T, Tint>(M, os);
      MyVector<Tint> eVect =
          unfold_opt(opt, "Failed to find an isotropic vector");
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: GetHyperbolicPlane eVect=" << StringVectorGAP(eVect)
         << "\n";
      T sum1 = EvaluationQuadForm<T, Tint>(M, eVect);
      if (sum1 != 0) {
        std::cerr << "MODEL: eVect is not an isotropic vector for M\n";
        throw TerminalException{1};
      }
#endif
      MyVector<Tint> NewVect = ThePerturb.transpose() * eVect;
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: GetHyperbolicPlane NewVect=" << StringVectorGAP(NewVect)
         << "\n";
      T sum2 = EvaluationQuadForm<T, Tint>(Qmat, NewVect);
      if (sum2 != 0) {
        std::cerr << "MODEL: NewVect is not an isotropic vector for Qmat\n";
        throw TerminalException{1};
      }
#endif
      if (SetVect.count(NewVect) == 0 || iter == n_iter) {
        SetVect.insert(NewVect);
        return NewVect;
      }
      iter += 1;
    }
  };

  auto has_non_orthogonal = [&](MyVector<Tint> const &NewVect) -> bool {
    MyVector<T> NewVect_T = UniversalVectorConversion<T, Tint>(NewVect);
    MyVector<T> prod = Qmat * NewVect_T;
    for (auto &eVect : SetVect) {
      MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect);
      T scal = eVect_T.dot(prod);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: has_non_orthogonal, scal=" << scal << "\n";
#endif
      if (scal != 0) {
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: has_non_orthogonal, eVect=" << StringVectorGAP(eVect)
           << " NewVect=" << StringVectorGAP(NewVect) << "\n";
#endif
        return true;
      }
    }
    return false;
  };
  auto f_get_minimal =
      [&](MyVector<Tint> const &NewVect) -> std::pair<T, MyVector<Tint>> {
    T min_scal;
    bool is_first = true;
    MyVector<Tint> BestVect;
    //
    MyVector<T> NewVect_T = UniversalVectorConversion<T, Tint>(NewVect);
    MyVector<T> prod = Qmat * NewVect_T;
    for (auto &eVect : SetVect) {
      MyVector<Tint> rVect = eVect;
      MyVector<T> eVect_T = UniversalVectorConversion<T, Tint>(eVect);
      T scal = eVect_T.dot(prod);
#ifdef DEBUG_APPROXIMATE_MODELS
      os << "MODEL: f_get_minimal, scal=" << scal << "\n";
#endif
      if (scal < 0) {
        scal = -scal;
        rVect = -eVect;
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
    os << "MODEL: f_get_minimal, BestVect=" << StringVectorGAP(BestVect)
       << " NewVect=" << StringVectorGAP(NewVect) << "\n";
#endif
    return {min_scal, BestVect};
  };

  // Now iterating
  bool is_first = true;
  std::pair<MyVector<Tint>, MyVector<Tint>> pairA;
  T min_scal;
  int n_at_level = 0;
  int max_level = 10;
  int n_no_change = 0;
  int max_no_change = 20;
  while (true) {
    MyVector<Tint> NewVect = f_insert_vect();
    if (has_non_orthogonal(NewVect)) {
      std::pair<T, MyVector<Tint>> pairB = f_get_minimal(NewVect);
      if (is_first) {
        min_scal = pairB.first;
        pairA = {NewVect, pairB.second};
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: f_get_minimal, first, min_scal=" << min_scal << "\n";
        os << "MODEL: f_get_minimal, 1: pairA=(" << StringVectorGAP(pairA.first)
           << " / " << StringVectorGAP(pairA.second) << ")\n";
#endif
        n_at_level = 0;
        is_first = false;
      } else {
#ifdef DEBUG_APPROXIMATE_MODELS
        os << "MODEL: f_get_minimal, pairB.first=" << pairB.first
           << " min_scal=" << min_scal << "\n";
#endif
        if (pairB.first < min_scal) {
          min_scal = pairB.first;
          pairA = {NewVect, pairB.second};
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: f_get_minimal, now min_scal=" << min_scal << "\n";
          os << "MODEL: f_get_minimal, 2: pairA=("
             << StringVectorGAP(pairA.first) << " / "
             << StringVectorGAP(pairA.second) << ")\n";
#endif
          n_at_level = 0;
          n_no_change = 0;
        } else {
#ifdef DEBUG_APPROXIMATE_MODELS
          os << "MODEL: f_get_minimal, n_at_level=" << n_at_level << "\n";
#endif
          if (pairB.first == min_scal) {
            n_at_level += 1;
          } else {
            n_no_change += 1;
          }
          bool terminate = false;
          if (n_no_change == max_no_change) {
#ifdef DEBUG_APPROXIMATE_MODELS
            os << "MODEL: f_get_minimal, 3: pairA=("
               << StringVectorGAP(pairA.first) << " / "
               << StringVectorGAP(pairA.second) << ")\n";
#endif
            terminate = true;
          }
          if (n_at_level == max_level) {
#ifdef DEBUG_APPROXIMATE_MODELS
            os << "MODEL: f_get_minimal, 4: pairA=("
               << StringVectorGAP(pairA.first) << " / "
               << StringVectorGAP(pairA.second) << ")\n";
#endif
            terminate = true;
          }
          if (terminate) {
#ifdef TIMINGS_APPROXIMATE_MODELS
            os << "MODEL: |GetHyperbolicPlane|=" << time << "\n";
#endif
            return get_result_hyperbolic_plane(Qmat, pairA);
          }
        }
      }
    }
  }
}

/*
  What are the strategies for getting decompositions?
  The objective is to get all of this things working.
  ---The embedding for Eichler actually gets us an embedding
  of the lattice for the quadratic form (aU + bU + L) into U + U + L.
  The index of that lattice into U + U + L is a b.
  ---What if the hyperbolic planes aU + bU + (aU + bU)^{perp}
  actually define a sublattice of the lattice K.
  ---We can probably rescale and that would work for us, though maybe at a cost.

     ----

  The function should return the following:
  ---A matrix P and a factor C such that
            P Q P^T = C (U + U + L)
     and L even integral.
  ---P does not have to be invertible but it has to be.

  P1 = FullBasis
  P1 Q P1^T = aU + bU + L
  Define Embed = (1/a , 1 , 1/b , 1, 1^{n-4})
  P2 = Embed P1
  P2 Q P2^T = U + U + L

 */
template <typename T, typename Tint> struct EichlerReduction {
  MyMatrix<Tint> Embed;
  MyMatrix<T> Embed_T;
  MyMatrix<T> EmbedInv_T;
  T scal;
  MyMatrix<T> QmatEichler;
};

template <typename T, typename Tint>
EichlerReduction<T, Tint> GetEichlerHyperplaneBasis(MyMatrix<T> const &Qmat,
                                                    std::ostream &os) {
#ifdef TIMINGS_APPROXIMATE_MODELS
  MicrosecondTime time;
#endif
  int n = Qmat.rows();
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: Beginning of GetEichlerHyperplaneBasis Qmat=\n";
  WriteMatrix(os, Qmat);
#endif
  ResultHyperbolicPlane<T, Tint> hyp1 = GetHyperbolicPlane<T, Tint>(Qmat, os);
  MyMatrix<Tint> const &Basis1 = hyp1.Basis;
  T a1 = hyp1.scal;
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: GetEichlerHyperplaneBasis, Basis1=\n";
  WriteMatrix(os, Basis1);
#endif
  MyMatrix<T> Basis1_T = UniversalMatrixConversion<T, Tint>(Basis1);
  MyMatrix<T> prod1 = Basis1_T * Qmat;
  MyMatrix<T> NSP_T = NullspaceIntTrMat(prod1);
  MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint, T>(NSP_T);
  MyMatrix<T> Qmat2 = NSP_T * Qmat * NSP_T.transpose();
  ResultHyperbolicPlane<T, Tint> hyp2 = GetHyperbolicPlane<T, Tint>(Qmat2, os);
  MyMatrix<Tint> const &Basis2 = hyp2.Basis;
  T a2 = hyp2.scal;
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: GetEichlerHyperplaneBasis, a1=" << a1 << " a2=" << a2 << "\n";
  os << "MODEL: GetEichlerHyperplaneBasis, Basis2=\n";
  WriteMatrix(os, Basis2);
#endif
  MyMatrix<Tint> Basis2_NSP = Basis2 * NSP;
  MyMatrix<Tint> HyperBasis = Concatenate(Basis1, Basis2_NSP);
  MyMatrix<T> HyperBasis_T = UniversalMatrixConversion<T, Tint>(HyperBasis);
  MyMatrix<T> prod2 = HyperBasis_T * Qmat;
  MyMatrix<T> NSP2_T = NullspaceIntTrMat(prod2);
  MyMatrix<Tint> NSP2 = UniversalMatrixConversion<Tint, T>(NSP2_T);
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: GetEichlerHyperplaneBasis, NSP2=\n";
  WriteMatrix(os, NSP2);
#endif
  MyMatrix<Tint> FullBasis = Concatenate(HyperBasis, NSP2);
  MyMatrix<T> FullBasis_T = UniversalMatrixConversion<T, Tint>(FullBasis);
  MyMatrix<T> EmbedSimpInv = IdentityMat<T>(n);
  EmbedSimpInv(0, 0) = 1 / a1;
  EmbedSimpInv(2, 2) = 1 / a2;
  MyMatrix<T> preEmbedInv_T = EmbedSimpInv * FullBasis_T;
  MyMatrix<T> QmatEichler = preEmbedInv_T * Qmat * preEmbedInv_T.transpose();
#ifdef DEBUG_APPROXIMATE_MODELS
  if (!is_eichler_canonical(QmatEichler)) {
    std::cerr << "The matrix should be reduced according to the Eichler "
                 "canonical form\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> preEmbed_T = Inverse(preEmbedInv_T);
  FractionMatrix<T> frac = RemoveFractionMatrixPlusCoeff(preEmbed_T);
  MyMatrix<T> Embed_T = frac.TheMat;
  MyMatrix<T> EmbedInv_T = Inverse(Embed_T);
  MyMatrix<Tint> Embed = UniversalMatrixConversion<Tint, T>(Embed_T);
  T scal = frac.TheMult * frac.TheMult;
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: GetEichlerHyperplaneBasis Embed=\n";
  WriteMatrix(os, Embed);
  os << "MODEL: GetEichlerHyperplaneBasis scal=" << scal << "\n";
  os << "MODEL: GetEichlerHyperplaneBasis Qmat=\n";
  WriteMatrix(os, Qmat);
  os << "MODEL: GetEichlerHyperplaneBasis QmatEichler=\n";
  WriteMatrix(os, QmatEichler);
#endif
#ifdef TIMINGS_APPROXIMATE_MODELS
  os << "MODEL: |GetEichlerHyperplaneBasis|=" << time << "\n";
#endif
  return {Embed, Embed_T, EmbedInv_T, scal, QmatEichler};
}

template <typename T, typename Tint> struct InternalApproxCase {
  MyMatrix<T> Qmat;
  EichlerReduction<T, Tint> er;
  std::vector<MyMatrix<Tint>> ListCoset;
  std::vector<MyMatrix<Tint>> ListGenerators;
  ApproximateModel<T, Tint> approx;
};

template <typename T, typename Tint, typename Tgroup>
ApproximateModel<T, Tint>
INDEF_FORM_GetApproximateModel(MyMatrix<T> const &Qmat, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  int n = Qmat.rows();
  EichlerReduction<T, Tint> er = GetEichlerHyperplaneBasis<T, Tint>(Qmat, os);
  ApproximateModel<T, Tint> approx =
      INDEF_FORM_EichlerCriterion_TwoHyperplanesEven<T, Tint, Tgroup>(
          er.QmatEichler);
  std::vector<MyMatrix<Tint>> ListGen = approx.GetApproximateGroup(os);
  std::vector<MyMatrix<T>> ListGen_T =
      UniversalStdVectorMatrixConversion<T, Tint>(ListGen);
  GeneralMatrixGroupHelper<T, Telt> helper{n};
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: INDEF_FORM_GetApproximateModel, Before "
        "LinearSpace_Stabilizer_RightCoset\n";
#endif
  Stab_RightCoset<T> stab_right_coset =
      LinearSpace_Stabilizer_RightCoset<T, Tgroup,
                                        GeneralMatrixGroupHelper<T, Telt>>(
          ListGen_T, helper, er.Embed_T, os);
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: INDEF_FORM_GetApproximateModel, After "
        "LinearSpace_Stabilizer_RightCoset\n";
#endif
  std::vector<MyMatrix<Tint>> ListCoset =
      stab_right_coset.coset_desc.template expand<Tint>();
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: We have ListCoset\n";
#endif
  std::vector<MyMatrix<Tint>> ListGenerators;
  for (auto &eGen_T : stab_right_coset.list_gen) {
    MyMatrix<T> fGen_T = er.Embed_T * eGen_T * er.EmbedInv_T;
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: fGen_T=\n";
    WriteMatrix(os, fGen_T);
    if (!IsIntegralMatrix(fGen_T)) {
      std::cerr << "fGen should be integral\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<Tint> fGen = UniversalMatrixConversion<Tint, T>(fGen_T);
    ListGenerators.push_back(fGen);
  }
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: We have ListGenerators\n";
#endif
  std::vector<MyMatrix<Tint>> GRPeasy =
      GetEasyIsometries<T, Tint>(er.QmatEichler, os);
  approx.SetListClassesOrbitwise(GRPeasy, os);
  InternalApproxCase<T, Tint> iac{Qmat, er, ListCoset, ListGenerators, approx};
  std::shared_ptr<InternalApproxCase<T, Tint>> shr_ptr =
      std::make_shared<InternalApproxCase<T, Tint>>(iac);
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: We have shr_ptr\n";
#endif

  std::function<std::vector<MyVector<Tint>>(MyVector<Tint>, T)>
      get_vector_representatives =
          [=](MyVector<Tint> const &eRepr,
              T const &X) -> std::vector<MyVector<Tint>> {
    std::vector<MyVector<Tint>> ListV;
    for (auto &eCos : shr_ptr->ListCoset) {
      MyVector<Tint> fRepr = eCos.transpose() * eRepr;
      std::optional<MyVector<Tint>> opt =
          SolutionIntMat(shr_ptr->er.Embed, fRepr);
      if (opt) {
        MyVector<Tint> const &eSol = *opt;
#ifdef DEBUG_APPROXIMATE_MODELS
        MyMatrix<T> const &Qmat = shr_ptr->Qmat;
        T eNorm = EvaluationQuadForm<T, Tint>(Qmat, eSol);
        if (eNorm != X) {
          std::cerr << "fSol is not of the right norm\n";
          throw TerminalException{1};
        }
#endif
        ListV.push_back(eSol);
      }
    }
    return ListV;
  };

  std::function<void(std::vector<MyMatrix<Tint>>, std::ostream &)>
      SetListClassesOrbitwise =
          [=]([[maybe_unused]] std::vector<MyMatrix<Tint>> const &GRPmatr,
              [[maybe_unused]] std::ostream &os) -> void {
    std::cerr << "That function should not be called\n";
    throw TerminalException{1};
  };
  std::function<std::vector<MyMatrix<Tint>>(std::ostream &)>
      GetApproximateGroup = [=]([[maybe_unused]] std::ostream &os)
      -> std::vector<MyMatrix<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of GetApproximateGroup : 3\n";
#endif
    return shr_ptr->ListGenerators;
  };
  std::function<std::vector<MyVector<Tint>>(T, std::ostream &)>
      GetCoveringOrbitRepresentatives =
          [=](T const &X, [[maybe_unused]] std::ostream &os)
      -> std::vector<MyVector<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of GetCoveringOrbitRepresentatives 3: X=" << X
       << "\n";
#endif
    T Xscal = X * shr_ptr->er.scal;
    std::vector<MyVector<Tint>> ListRepr;
    for (auto &eRepr :
         shr_ptr->approx.GetCoveringOrbitRepresentatives(Xscal, os)) {
      for (auto &fRepr : get_vector_representatives(eRepr, X)) {
        ListRepr.push_back(fRepr);
      }
    }
    return ListRepr;
  };
  std::function<std::optional<MyVector<Tint>>(T const &, std::ostream &)>
      GetOneOrbitRepresentative =
          [=](T const &X, [[maybe_unused]] std::ostream &os)
      -> std::optional<MyVector<Tint>> {
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Beginning of GetOneOrbitRepresentative 3: X=" << X << "\n";
#endif
    T Xscal = X * shr_ptr->er.scal;
    std::optional<MyVector<Tint>> opt =
        shr_ptr->approx.GetOneOrbitRepresentative(Xscal, os);
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: We have opt\n";
#endif
    if (!opt) {
      return {};
    }
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Before iteration over |ListCoset|=" << ListCoset.size()
       << "\n";
#endif
    MyVector<Tint> const &eRepr = *opt;
    std::vector<MyVector<Tint>> ListV = get_vector_representatives(eRepr, X);
    if (ListV.size() > 0) {
      return ListV[0];
    }
    // We wanted to avoid that but sometimes we cannot
#ifdef DEBUG_APPROXIMATE_MODELS
    os << "MODEL: Nothing found, now enumerating the possibilities\n";
#endif
    std::vector<MyVector<Tint>> ListRepr =
        GetCoveringOrbitRepresentatives(X, os);
    if (ListRepr.size() > 0) {
      return ListRepr[0];
    }
    return {};
  };
  return {GetApproximateGroup, SetListClassesOrbitwise,
          GetCoveringOrbitRepresentatives, GetOneOrbitRepresentative};
}

template <typename T, typename Tint> struct FirstNorm {
  T X;
  MyVector<Tint> eVect;
};

template <typename T, typename Tint>
FirstNorm<T, Tint> GetFirstNorm(ApproximateModel<T, Tint> const &approx,
                                std::ostream &os) {
  T X = 1;
#ifdef DEBUG_APPROXIMATE_MODELS
  os << "MODEL: Beginning of GetFirstNorm\n";
#endif
  while (true) {
    std::optional<MyVector<Tint>> opt = approx.GetOneOrbitRepresentative(X, os);
    if (opt) {
      MyVector<Tint> const &V = *opt;
      return {X, V};
    }
    X += 1;
  }
}

// clang-format off
#endif  // SRC_INDEFINITE_MODELS_APPROXIMATEMODELS_H_
// clang-format on
