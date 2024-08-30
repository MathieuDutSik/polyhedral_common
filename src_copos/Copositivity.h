// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_COPOS_COPOSITIVITY_H_
#define SRC_COPOS_COPOSITIVITY_H_

// clang-format off
#include "LatticeDefinitions.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "POLY_LinearProgramming.h"
#include "Positivity.h"
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_COPOSITIVITY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_COPOSITIVITY
#endif

template <typename T>
bool TestCopositivityByPositivityCoeff(MyMatrix<T> const &eSymmMatB) {
  int n = eSymmMatB.rows();
  for (int i = 0; i < n; i++) {
    T eNorm = eSymmMatB(i, i);
    if (eNorm <= 0)
      return false;
  }
  for (int i = 0; i < n - 1; i++)
    for (int j = i + 1; j < n; j++) {
      T eScal = eSymmMatB(i, j);
      if (eScal < 0)
        return false;
    }
  return true;
}

template <typename T> struct CopositivityInfoReduction {
  MyMatrix<T> TheBasisReduced;
  MyMatrix<T> eSymmMatReduced;
  MyMatrix<T> Pmat;
  std::vector<int> IdxVector;
  std::vector<int> IdxColumn;
};

// #define DEBUG_UNDER_POS_COND

// Given C >0 a >=0 and b>0 find the maximum
// value of x>=0 such that a x + b x^2 <= C
// reduced form is a2 x + x^2 <= C2
// polynomial is   x^2 + a2 x - C2 = 0
//               a x^2 +  b x + c  = 0
template <typename T> int FindLargest(T const &a, T const &b, T const &C) {
  T C2 = C / b;
  T a2 = a / b;
  double a2_doubl = UniversalScalarConversion<double, T>(a2);
  double C2_doubl = UniversalScalarConversion<double, T>(C2);
#ifdef SANITY_CHECK_COPOSITIVITY
  if (b <= 0 || C <= 0) {
    std::cerr << "b should be strictly positive. b=" << b << "\n";
    std::cerr << "C should be strictly positive. C=" << C << "\n";
    throw TerminalException{1};
  }
#endif
  double delta = a2_doubl * a2_doubl + 4 * C2_doubl;
  double x1 = 0.5 * (-a2_doubl + sqrt(delta));
  int eReturn = UniversalScalarConversion<int, double>(x1);
  auto f = [&](int const &x) -> bool {
    T eDiff = C2 - a2 * x - x * x;
    if (eDiff >= 0)
      return true;
    return false;
  };
  while (true) {
    bool test1 = f(eReturn);
    bool test2 = f(eReturn + 1);
    if (test1 && !test2)
      break;
    if (!test1)
      eReturn--;
    if (test2)
      eReturn++;
  }
  return eReturn;
}

// This function is the kernel of our methodology
// We have a symmetric matrix eSymmMat and
// a family of vectors in TheBasis.
//
// Assumptions
// ---TheBasis is a set of vectors (v1, ..., vn)
// ---We have vi in the positive orthant R_{\geq 0}^n
// ---(v1, ...., vn) spans a simplicial cone. We have no constraint on det(v1,
// ..., vn). Just be non-zero
// ---For all 1<= i,j <= n  we have vi eSymmMat * vj >= 0
template <typename T, typename Tint, typename F>
void EnumerateShortVectorInCone_UnderPositivityCond_F(MyMatrix<T> const &eSymmMat,
                                                      MyMatrix<Tint> const &TheBasis,
                                                      T const &MaxNorm, F f) {
#ifdef SANITY_CHECK_COPOSITIVITY
  MyMatrix<T> TheBasis_T = UniversalMatrixConversion<T, Tint>(TheBasis);
  MyMatrix<T> tstSymmMatB = TheBasis_T * eSymmMat * TheBasis_T.transpose();
  bool test1 = TestCopositivityByPositivityCoeff(tstSymmMatB);
  if (!test1) {
    std::cerr << "Inconsistency in the positivity of coefficients\n";
    throw TerminalException{1};
  }
#endif
  BasisReduction<Tint> eRedBas = ComputeBasisReduction(TheBasis);
  MyMatrix<Tint> TheBasisReord = eRedBas.TheBasisReord;
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(eRedBas.Pmat);
  MyMatrix<T> Pinv = Pmat_T.inverse().eval();
  MyMatrix<Tint> Pinv_int = UniversalMatrixConversion<Tint, T>(Pinv);
  MyMatrix<Tint> Pinv_cgr = Pinv_int.transpose();
  MyMatrix<T> eSymmMatB = Pinv * eSymmMat * Pinv.transpose();
#ifdef DEBUG_COPOSITIVITY
  std::cerr << "eRedBas.Pmat=\n";
  WriteMatrix(std::cerr, eRedBas.Pmat);
  std::cerr << "eSymmMat=\n";
  WriteMatrix(std::cerr, eSymmMat);
  std::cerr << "TheBasis=\n";
  WriteMatrix(std::cerr, TheBasis);
  //
  std::cerr << "eSymmMatB=\n";
  WriteMatrix(std::cerr, eSymmMatB);
  std::cerr << "TheBasisReord=\n";
  WriteMatrix(std::cerr, TheBasisReord);
#endif
  MyMatrix<T> TheBasisReord_T =
      UniversalMatrixConversion<T, Tint>(TheBasisReord);
  MyMatrix<T> eSymmMatC =
      TheBasisReord_T * eSymmMatB * TheBasisReord_T.transpose();
#ifdef SANITY_CHECK_COPOSITIVITY
  bool test2 = TestCopositivityByPositivityCoeff(eSymmMatC);
  if (!test2) {
    std::cerr << "Inconsistency in the positivity of coefficients\n";
    std::cerr << "and a matrix error as well (We should have test1=test2)\n";
    throw TerminalException{1};
  }
#endif
  int n = eSymmMat.cols();
  std::vector<T> AII(n);
  for (int i = 0; i < n; i++)
    AII[i] = TheBasisReord(i, i);
#ifdef DEBUG_COPOSITIVITY
  std::cerr << "AII=";
  for (int i = 0; i < n; i++)
    std::cerr << " " << AII[i];
  std::cerr << "\n";
#endif
  //
  // Now the variables whose value will change over the iteration
  //
  // x = sum_i lambda_i v_i
  // with lambda_i = h_i + z_i/a_ii  with 0 <= z_i <= z_max(i)
  //
  // We define the partial quadratic forms
  // q_i = q[sum_{j=i}^{n-1} lambda_j v_j].
  // this expression depends only on the X[j] for j>=i.
  // q_0 is the full quadratic form.
  // q_i(Z_i, Z_{i+1}, Z_{n-1}) = q_ii Z_i^2 + q_{i+1}(Z_{i+1}, Z_{i+2}, ...,
  // Z_{n-1})
  //        + W_{i} Z_i
  //

  // Partial norms of vectors
  std::vector<T> ListNorm(n);
  // The sought after vector
  MyVector<Tint> X(n);
  // Xmax values
  std::vector<Tint> Zmax(n);
  // X values over which we iterate
  std::vector<Tint> Z(n, 0);
  // List of H: smallest values of lambda
  std::vector<T> ListH(n);
  // values of the Lambda_i:
  std::vector<T> Lambda(n);
  // PartialSum is defined as
  // PartialSum[i][j]
  MyMatrix<T> PartialSum = ZeroMatrix<T>(n + 1, n + 1);
  auto ComputeHXZmaxNorm = [&](int const &nupdt) -> void {
  // First we can compute the ListH just from the X
#ifdef DEBUG_COPOSITIVITY
    std::cerr << "ComputeHXZmaxNorm : Before computing ListH, Lambda and X\n";
#endif
    for (int i = nupdt - 1; i >= 0; i--) {
      T eNumber = 0;
      for (int j = i + 1; j < n; j++)
        eNumber += TheBasisReord(j, i) * Lambda[j];
      T eFr = FractionalPart(-eNumber) / AII[i];
      ListH[i] = eFr;
      T eZ = Z[i];
      Lambda[i] = eZ / AII[i] + ListH[i];
      T eX_q = eNumber + Lambda[i] * AII[i];
#ifdef SANITY_CHECK_COPOSITIVITY
      mpz_class eDen = eX_q.get_den();
      if (eDen > 1) {
        std::cerr << "eDen should be equal to 1\n";
        std::cerr << "eZ=" << eZ << "\n";
        std::cerr << "eX_q=" << eX_q << "\n";
        throw TerminalException{1};
      }
#endif
      X(i) = UniversalScalarConversion<Tint, T>(eX_q);
    }
#ifdef DEBUG_COPOSITIVITY
    std::cerr << "X=";
    for (int i = 0; i < n; i++)
      std::cerr << " " << X(i);
    std::cerr << "\n";
    std::cerr << "Z=";
    for (int i = 0; i < n; i++)
      std::cerr << " " << Z[i];
    std::cerr << "\n";
    std::cerr << "Lambda=";
    for (int i = 0; i < n; i++)
      std::cerr << " " << Lambda[i];
    std::cerr << "\n";
    std::cerr << "ListH=";
    for (int i = 0; i < n; i++)
      std::cerr << " " << ListH[i];
    std::cerr << "\n";
    std::cerr << "Computing norms and Zmax\n";
    std::cerr << "n=" << n << "\n";
#endif
    for (int i = nupdt - 1; i >= 0; i--) {
#ifdef DEBUG_COPOSITIVITY
      std::cerr << "i=" << i << " / " << nupdt << "\n";
#endif
      T Qii = eSymmMatC(i, i);
      T eNorm = Qii * Lambda[i] * Lambda[i];
      T eW = 0;
#ifdef DEBUG_COPOSITIVITY
      std::cerr << "Norm compute step 1\n";
#endif
      for (int j = i + 1; j < n; j++)
        eW += 2 * eSymmMatC(i, j) * Lambda[j];
#ifdef DEBUG_COPOSITIVITY
      std::cerr << "Norm compute step 2\n";
#endif
      T NextNorm;
      if (i < n - 1) {
        NextNorm = ListNorm[i + 1];
      } else {
        NextNorm = 0;
      }
      eNorm += eW * Lambda[i] + NextNorm;
      ListNorm[i] = eNorm;
      T eH = ListH[i];
      T eVal0 = NextNorm + eW * eH + Qii * eH * eH;
      if (eVal0 >= MaxNorm) {
        Zmax[i] = 0;
      } else {
        T eVal1 = (eW + 2 * eH * Qii) / AII[i];
        T eVal2 = Qii / (AII[i] * AII[i]);
        T C = MaxNorm - eVal0;
        Zmax[i] = FindLargest(eVal1, eVal2, C);
      }
    }
#ifdef DEBUG_COPOSITIVITY
    std::cerr << "Zmax=";
    for (int i = 0; i < n; i++)
      std::cerr << " " << Zmax[i];
    std::cerr << "\n";
    std::cerr << "ListNorm=";
    for (int i = 0; i < n; i++)
      std::cerr << " " << ListNorm[i];
    std::cerr << "\n";
#endif
  };
  //
  // The initialization
  //
  ComputeHXZmaxNorm(n);
#ifdef DEBUG_COPOSITIVITY
  std::cerr << "After ComputeHXZmaxNorm\n";
#endif
  //
  // The tree search function
  //
  auto NextInTree = [&]() -> bool {
#ifdef DEBUG_COPOSITIVITY
    for (int i = 0; i < n; i++)
      std::cerr << "i=" << i << " Zmax=" << Zmax[i] << "\n";
#endif
    for (int i = 0; i < n; i++) {
      if (Z[i] < Zmax[i]) {
#ifdef DEBUG_COPOSITIVITY
        std::cerr << "Clear update here Z[i]=" << Z[i] << "\n";
#endif
        Z[i]++;
#ifdef DEBUG_COPOSITIVITY
        std::cerr << "i=" << i << " Z[i]=" << Z[i] << "\n";
#endif
        for (int j = 0; j < i; j++)
          Z[j] = 0;
        ComputeHXZmaxNorm(i + 1);
        return true;
      }
    }
    return false;
  };
  while (true) {
    if (ListNorm[0] <= MaxNorm) {
#ifdef DEBUG_COPOSITIVITY
      std::cerr << "Before computing NewX\n";
      std::cerr << "Pinv_int=\n";
      WriteMatrix(std::cerr, Pinv_int);
#endif
      MyVector<Tint> NewX = Pinv_cgr * X;
#ifdef SANITY_CHECK_COPOSITIVITY
      T eVal = EvaluationQuadForm(eSymmMat, NewX);
      if (eVal != ListNorm[0]) {
        std::cerr << "Norm inconsistency in the code\n";
        std::cerr << "eVal=" << eVal << " ListNorm[0]=" << ListNorm[0] << "\n";
        std::cerr << "NewX=";
        WriteVectorNoDim(std::cerr, NewX);
        std::cerr << "X=";
        WriteVectorNoDim(std::cerr, X);
        std::cerr << "\n";
        for (int i = 0; i < n; i++) {
          std::cerr << "i=" << i << " norm=" << ListNorm[i] << "\n";
        }
        throw TerminalException{1};
      }
      if (NewX.minCoeff() < 0) {
        std::cerr << "X=";
        WriteVectorNoDim(std::cerr, X);
        std::cerr << "Error. We found a negative vector\n";
        std::cerr << "eSymmMat=\n";
        WriteMatrix(std::cerr, eSymmMat);
        std::cerr << "TheBasis=\n";
        WriteMatrix(std::cerr, TheBasis);
        std::cerr << "MaxNorm=" << MaxNorm << "\n";
        throw TerminalException{1};
      }
#endif
      f(NewX, ListNorm[0]);
    }
#ifdef DEBUG_COPOSITIVITY
    std::cerr << "X=[";
    for (int i = 0; i < n; i++)
      std::cerr << " " << X(i);
    std::cerr << "]\n";
    std::cerr << "Z=[";
    for (int i = 0; i < n; i++)
      std::cerr << " " << Z[i];
    std::cerr << "]\n";
    std::cerr << "Lambda=[";
    for (int i = 0; i < n; i++)
      std::cerr << " " << Lambda[i];
    std::cerr << "]\n";
#endif
    bool test = NextInTree();
    if (!test)
      break;
  }
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
EnumerateShortVectorInCone_UnderPositivityCond(MyMatrix<T> const &eSymmMat,
                                               MyMatrix<Tint> const &TheBasis,
                                               T const &MaxNorm) {
  std::vector<MyVector<Tint>> ListVectRet;
  auto f=[&](MyVector<Tint> const& V, [[maybe_unused]] T const& norm) -> void {
    ListVectRet.push_back(V);
  };
  EnumerateShortVectorInCone_UnderPositivityCond_F(eSymmMat, TheBasis, MaxNorm, f);
  return ListVectRet;
}

template <typename Tint> struct SingleTestResult {
  bool test;
  std::string strNature;
  MyVector<Tint> eVectResult1;
};

template <typename Tint>
void WriteSingleTestResult(
    std::ostream &os, std::string const &OutFormat,
    std::pair<SingleTestResult<Tint>, size_t> const &eRes) {
  SingleTestResult<Tint> const &eResult = eRes.first;
  if (OutFormat == "classic") {
    os << "Total number of cones used in the computation nbCone=" << eRes.second
       << "\n";
    if (eResult.test) {
      os << "The matrix is indeed copositive\n";
    } else {
      os << "The matrix is not copositive\n";
      os << "Nature of violation=" << eResult.strNature << "\n";
      os << "V=";
      WriteVectorNoDim(os, eResult.eVectResult1);
    }
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(nbCone:=" << eRes.second
       << ", isCopositive:=" << GAP_logical(eResult.test);
    if (!eResult.test) {
      os << ", violation_nature:=\"" << eResult.strNature << "\"";
      os << ", V:=" << StringVectorGAP(eResult.eVectResult1);
    }
    os << ");\n";
    return;
  }
  std::cerr << "WriteSingleTestResult: Failed to find a matching entry\n";
  throw TerminalException{1};
}

template <typename Tint> struct CopositivityEnumResult {
  bool test;
  std::vector<MyVector<Tint>> TotalListVect;
  SingleTestResult<Tint> eResult;
};

template <typename T, typename Tint>
void WriteCopositivityEnumResult(std::ostream &os, std::string const &OutFormat,
                                 MyMatrix<T> const &eSymmMat,
                                 CopositivityEnumResult<Tint> const &CopoRes) {
  if (OutFormat == "classic") {
    if (CopoRes.test == false) {
      os << "Matrix is not Copositive\n";
      os << "Nature=" << CopoRes.eResult.strNature << "\n";
      os << "eVect1=";
      WriteVectorNoDim(os, CopoRes.eResult.eVectResult1);
    } else {
      int n = eSymmMat.rows();
      os << "Matrix is Copositive\n";
      int nbVect = CopoRes.TotalListVect.size();
      os << "nbVect=" << nbVect << "\n";
      for (int iVect = 0; iVect < nbVect; iVect++) {
        MyVector<Tint> eVect = CopoRes.TotalListVect[iVect];
        T eNorm = EvaluationQuadForm(eSymmMat, eVect);
        os << "iVect=" << iVect << " eVect=[";
        for (int i = 0; i < n; i++) {
          if (i > 0)
            os << ",";
          os << eVect[i];
        }
        os << "] norm=" << eNorm << "\n";
      }
    }
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(isCopositive:=" << GAP_logical(CopoRes.test);
    if (CopoRes.test) {
      os << ", nature:=\"" << CopoRes.eResult.strNature << "\"";
      os << ", eVect:=" << StringVectorGAP(CopoRes.eResult.eVectResult1);
    } else {
      int nbVect = CopoRes.TotalListVect.size();
      os << ", ListVect:=[";
      for (int iVect = 0; iVect < nbVect; iVect++) {
        if (iVect > 0)
          os << ",\n";
        MyVector<Tint> eVect = CopoRes.TotalListVect[iVect];
        T eNorm = EvaluationQuadForm(eSymmMat, eVect);
        os << "rec(vect:=" << StringVectorGAP(eVect) << ", norm:=" << eNorm
           << ")";
      }
      os << "]";
    }
    os << ");\n";
    return;
  }
  std::cerr << "WriteCopositivityEnumResult: Failed to find a matching entry\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
std::vector<MyMatrix<Tint>> PairDecomposition(MyMatrix<T> const &eSymmMat,
                                              MyMatrix<Tint> const &TheBasis) {
  int n = eSymmMat.rows();
  MyMatrix<T> TheBasis_T = UniversalMatrixConversion<T, Tint>(TheBasis);
  MyMatrix<T> eSymmMatB = TheBasis_T * eSymmMat * TheBasis_T.transpose();
  std::pair<Tint, Tint> PairXY_found{-1, -1};
  std::pair<int, int> PairIJ_found{-1, -1};
  // We form a matrix
  // | a b |
  // | b c |
  // and we look for an integer vector (x,y) with x >= 0 and y >= 0 such that
  // a x + b y >= 0 and a b + c y >= 0
  // In practice it means that when we have a negative value
  auto GetSplit = [](T const &a, T const &b,
                     T const &c) -> std::pair<Tint, Tint> {
    int eSum = 2;
    while (true) {
      for (int x = 1; x < eSum; x++) {
        int y = eSum - x;
        T val1 = a * x + b * y;
        T val2 = b * x + c * y;
        if (val1 >= 0 && val2 >= 0)
          return {Tint(x), Tint(y)};
      }
      eSum++;
    }
  };
  int OptMinimum = 1;
  bool IsFirst = true;
  Tint MinSumXY = 0;
  T MaxValQuot = 0;
  for (int i = 0; i < n - 1; i++)
    for (int j = i + 1; j < n; j++) {
      T a = eSymmMatB(i, i);
      T b = eSymmMatB(i, j);
      T c = eSymmMatB(j, j);
      if (b < 0) {
        std::pair<Tint, Tint> CandXY;
        Tint sumCoef = 0;
        T eQuot = 0;
        if (OptMinimum == 1) {
          CandXY = {1, 1};
          eQuot = (b * b) / (a * c);
        }
        if (OptMinimum == 2) {
          CandXY = GetSplit(a, b, c);
          sumCoef = CandXY.first + CandXY.second;
        }
        auto InsertValue = [&]() -> bool {
          if (IsFirst) {
            IsFirst = false;
            MaxValQuot = eQuot;
            MinSumXY = sumCoef;
            return true;
          }
          if (OptMinimum == 1 && eQuot > MaxValQuot) {
            MaxValQuot = eQuot;
            return true;
          }
          if (OptMinimum == 2 && sumCoef < MinSumXY) {
            MinSumXY = sumCoef;
            return true;
          }
          return false;
        };
        bool test = InsertValue();
        if (test) {
          PairXY_found = CandXY;
          PairIJ_found = {i, j};
        }
      }
    }
  int iFound = PairIJ_found.first;
  int jFound = PairIJ_found.second;
#ifdef SANITY_CHECK_COPOSITIVITY
  if (iFound == -1) {
    std::cerr << "iFound=" << iFound << " jFound=" << jFound
              << " x=" << PairXY_found.first << " y=" << PairXY_found.second
              << "\n";
    std::cerr << "eSymmMatB=\n";
    WriteMatrix(std::cerr, eSymmMatB);
    std::cerr << "TheBasis=\n";
    WriteMatrix(std::cerr, TheBasis);
    std::cerr << "Looks like the matrix is all non-negative. No need for "
                 "splitting. Most likely a bug\n";
    throw TerminalException{1};
  }
#endif
#ifdef PRINT_SPLIT_CONE
  std::cerr << "iFound=" << iFound << " jFound=" << jFound
            << " PairXY_found=" << PairXY_found.first << " / "
            << PairXY_found.second << "\n";
#endif

  MyVector<Tint> eVectMid(n);
  for (int i = 0; i < n; i++)
    eVectMid(i) = PairXY_found.first * TheBasis(iFound, i) +
                  PairXY_found.second * TheBasis(jFound, i);
  MyMatrix<Tint> TheBasis1(n, n);
  MyMatrix<Tint> TheBasis2(n, n);
  for (int i = 0; i < n; i++) {
    MyVector<Tint> eVect1;
    MyVector<Tint> eVect2;
    if (i == iFound) {
      eVect1 = eVectMid;
    } else {
      eVect1 = TheBasis.row(i);
    }
    //
    if (i == jFound) {
      eVect2 = eVectMid;
    } else {
      eVect2 = TheBasis.row(i);
    }
    TheBasis1.row(i) = eVect1;
    TheBasis2.row(i) = eVect2;
  }
  MyMatrix<T> TheBasis1_T = UniversalMatrixConversion<T, Tint>(TheBasis1);
  MyMatrix<T> TheBasis2_T = UniversalMatrixConversion<T, Tint>(TheBasis2);
  MyMatrix<T> eSymmMatB1 = TheBasis1_T * eSymmMat * TheBasis1_T.transpose();
  MyMatrix<T> eSymmMatB2 = TheBasis2_T * eSymmMat * TheBasis2_T.transpose();
#ifdef SANITY_CHECK_COPOSITIVITY
  if (OptMinimum == 2) {
    if (eSymmMatB1(iFound, jFound) < 0 || eSymmMatB2(iFound, jFound) < 0) {
      std::cerr << "We clearly have a bug in our computation\n";
      throw TerminalException{1};
    }
  }
#endif
  return {std::move(TheBasis1), std::move(TheBasis2)};
}

template <typename T, typename Tint>
SingleTestResult<Tint>
SingleTestStrictCopositivity(MyMatrix<T> const &eSymmMat,
                             MyMatrix<Tint> const &TheBasis,
                             MyMatrix<T> const &eSymmMatB) {
  int n = eSymmMat.rows();
  MyVector<Tint> eVectZero;
  for (int i = 0; i < n; i++)
    if (eSymmMatB(i, i) <= 0)
      return {false, "Not Positive diag", TheBasis.row(i)};
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) {
      T a = eSymmMatB(i, i);
      T b = eSymmMatB(i, j);
      T c = eSymmMatB(j, j);
      if (b < 0) {
        if (b * b == a * c) {
          MyVector<T> V(2);
          V(0) = -b;
          V(1) = a;
          MyVector<T> Vred = RemoveFractionVector(V);
          MyVector<Tint> Vint = UniversalVectorConversion<Tint, T>(Vred);
          MyVector<Tint> eVect =
              Vint(0) * TheBasis.row(i) + Vint(1) * TheBasis.row(j);
          return {false, "zero vector detected", std::move(eVect)};
        }
        if (b * b > a * c) {
          int eSum = 2;
          while (true) {
            for (int x = 1; x < eSum; x++) {
              int y = eSum - x;
              T qVal = a * x * x + 2 * b * x * y + c * y * y;
              if (qVal < 0) {
                MyVector<Tint> eVect =
                    x * TheBasis.row(i) + y * TheBasis.row(j);
                return {false, "Off diagonal violation", std::move(eVect)};
              }
            }
            eSum++;
          }
        }
      }
    }
  return {true, "on the surface ok", std::move(eVectZero)};
}

template <typename T, typename Tint>
SingleTestResult<Tint> SingleTestCopositivity(MyMatrix<T> const &eSymmMat,
                                              MyMatrix<Tint> const &TheBasis,
                                              MyMatrix<T> const &eSymmMatB) {
  int n = eSymmMat.rows();
  MyVector<Tint> eVectZero;
  for (int i = 0; i < n; i++)
    if (eSymmMatB(i, i) < 0)
      return {false, "Not Positive diag", TheBasis.row(i)};
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) {
      T a = eSymmMatB(i, i);
      T b = eSymmMatB(i, j);
      T c = eSymmMatB(j, j);
      if (b < 0 && b * b > a * c) {
        int eSum = 2;
        while (true) {
          for (int x = 1; x < eSum; x++) {
            int y = eSum - x;
            T qVal = a * x * x + 2 * b * x * y + c * y * y;
            if (qVal < 0) {
              MyVector<Tint> eVect = x * TheBasis.row(i) + y * TheBasis.row(j);
              return {false, "Off diagonal violation", std::move(eVect)};
            }
          }
          eSum++;
        }
      }
    }
  return {true, "on the surface ok", std::move(eVectZero)};
}

template <typename T, typename Tint>
SingleTestResult<Tint> SearchByZeroInKernel(MyMatrix<T> const &eSymmMat,
                                            std::ostream &os) {
  MyVector<Tint> eVectZero;
  MyMatrix<T> NSP = NullspaceMat(eSymmMat);
  int nbCol = NSP.cols();
  int nbNSP = NSP.rows();
  if (nbNSP == 0)
    return {true, "on the surface ok", eVectZero};
  MyMatrix<T> FAC(nbCol + 1, nbNSP + 1);
  for (int iCol = 0; iCol < nbCol; iCol++) {
    FAC(iCol, 0) = 0;
    for (int iNSP = 0; iNSP < nbNSP; iNSP++)
      FAC(iCol, iNSP + 1) = NSP(iNSP, iCol);
  }
  MyVector<T> ToBeMinimized(nbNSP + 1);
  ToBeMinimized(0) = 0;
  FAC(nbCol, 0) = -1;
  for (int iNSP = 0; iNSP < nbNSP; iNSP++) {
    T eSum = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      eSum += NSP(iNSP, iCol);
    ToBeMinimized(iNSP + 1) = eSum;
    FAC(nbCol, iNSP + 1) = eSum;
  }
  LpSolution<T> eSol = CDD_LinearProgramming(FAC, ToBeMinimized, os);
  if (!eSol.PrimalDefined)
    return {true, "on the surface ok", eVectZero};
  MyVector<T> eVect1(nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++) {
    T eSum = 0;
    for (int iNSP = 0; iNSP < nbNSP; iNSP++)
      eSum += eSol.DirectSolution(iNSP) * NSP(iNSP, iCol);
    eVect1(iCol) = eSum;
  }
  MyVector<T> eVect2 = RemoveFractionVector(eVect1);
  MyVector<Tint> eVect3 = UniversalVectorConversion<Tint, T>(eVect2);
  return {false, "Zero vector from kernel", std::move(eVect3)};
}

template <typename T, typename Tint>
CopositivityEnumResult<Tint>
KernelEnumerateShortVectorInCone(MyMatrix<T> const &eSymmMat,
                                 MyMatrix<Tint> const &TheBasis,
                                 T const& MaxNorm) {
  std::vector<MyVector<Tint>> TotalList;
  //
  // First clear up with all signs being positive case
  //
  MyMatrix<T> TheBasis_T = UniversalMatrixConversion<T, Tint>(TheBasis);
  MyMatrix<T> eSymmMatB = TheBasis_T * eSymmMat * TheBasis_T.transpose();
  bool test = TestCopositivityByPositivityCoeff(eSymmMatB);
  if (test) {
    TotalList = EnumerateShortVectorInCone_UnderPositivityCond<T, Tint>(
        eSymmMat, TheBasis, MaxNorm);
    return {true, std::move(TotalList), {}};
  }
  //
  // Second see if the vectors allow to conclude directly.
  //
  SingleTestResult<Tint> eResult =
      SingleTestStrictCopositivity(eSymmMat, TheBasis, eSymmMatB);
  if (!eResult.test)
    return {false, {}, eResult};
  //
  // Third, split the domain.
  //
  std::vector<MyMatrix<Tint>> ListBasis = PairDecomposition(eSymmMat, TheBasis);
  std::vector<MyMatrix<Tint>> ListBasisInt;
  for (auto &eBasis : ListBasis) {
    CopositivityEnumResult<Tint> fEnumResult =
        KernelEnumerateShortVectorInCone(eSymmMat, eBasis, MaxNorm);
    if (!fEnumResult.test)
      return {false, {}, fEnumResult.eResult};
    TotalList.insert(TotalList.end(), fEnumResult.TotalListVect.begin(),
                     fEnumResult.TotalListVect.end());
  }
  return {true, std::move(TotalList), {}};
}

// This is supposed to be a tree version and so faster than
// KernelEnumerateShortVectorInCone The other program is recursive and so has a
// problem with too large level.
//
// ---fSingleTest checks if the cone is adequate for conclusion.
//    If yes, then the conclusion is returned.
// ---fInsert checks if the cone is adequate and do the business (insert vector,
//       increase nbCone, etc.)
//    if it is not a boolean is returned and the cone is further subdivided.
template <typename T, typename Tint, typename Finsert, typename Ftest>
SingleTestResult<Tint> EnumerateCopositiveShortVector_Kernel(
    MyMatrix<T> const &eSymmMat, MyMatrix<Tint> const &InitialBasis,
    Finsert f_insert, Ftest f_test, [[maybe_unused]] std::ostream &os) {
  struct DataPair {
    int idx;
    MyMatrix<Tint> TheBasis0;
    MyMatrix<Tint> TheBasis1;
  };
  std::vector<DataPair> ListPair;
  auto SetValue = [&](int const &pos, DataPair const &ePair) -> void {
    int siz = ListPair.size();
    if (pos == siz) {
      ListPair.push_back(ePair);
    } else {
      ListPair[pos] = ePair;
    }
  };
  //
  SingleTestResult<Tint> eResult{true, "all ok", {}};
  //
  int position = -1;
  auto GetBasis = [&]() -> MyMatrix<Tint> {
    if (position == -1)
      return InitialBasis;
    int idx = ListPair[position].idx;
    if (idx == 0)
      return ListPair[position].TheBasis0;
    return ListPair[position].TheBasis1;
  };
  auto GoUpNextInTree = [&]() -> bool {
    while (true) {
      if (position == -1)
        return false;
      int idx = ListPair[position].idx;
      if (idx == 0) {
        ListPair[position].idx = 1;
        return true;
      }
      position--;
    }
  };
  auto GetDataPair =
      [&](std::vector<MyMatrix<Tint>> const &ListBasis) -> DataPair {
    return {0, ListBasis[0], ListBasis[1]};
  };
  auto NextInTree = [&]() -> bool {
    MyMatrix<Tint> TheBasis = GetBasis();
    MyMatrix<T> TheBasis_T = UniversalMatrixConversion<T, Tint>(TheBasis);
    MyMatrix<T> eSymmMatB = TheBasis_T * eSymmMat * TheBasis_T.transpose();
    //
    // If we can conclude the computation then we do so.
    //
    SingleTestResult<Tint> fResult = f_test(TheBasis, eSymmMatB);
    if (!fResult.test) {
      eResult = fResult;
      return false;
    }
    //
    // If we have vi A vj >= 0 then direct approach works.
    //
    bool testB = f_insert(TheBasis, eSymmMatB);
    if (testB) {
      return GoUpNextInTree();
    }
    //
    // Now doing the split
    //
    std::vector<MyMatrix<Tint>> ListBasis =
        PairDecomposition(eSymmMat, TheBasis);
    position++;
    SetValue(position, GetDataPair(ListBasis));
    return true;
  };
  while (true) {
    bool res = NextInTree();
    if (!res)
      break;
  }
  return eResult;
}

template <typename T, typename Tint>
Tshortest<T, Tint>
CopositiveShortestVector(MyMatrix<T> const &eSymmMat,
                         MyMatrix<Tint> const &InitialBasis,
                         std::ostream &os) {
#ifdef SANITY_CHECK_COPOSITIVITY
  SingleTestResult<Tint> kerResult =
      SearchByZeroInKernel<T, Tint>(eSymmMat, os);
  if (!kerResult.test) {
    std::cerr << "Inconsistency in the run. A bug to be solved\n";
    throw TerminalException{1};
  }
#endif
  int n = eSymmMat.rows();
  T MinNorm = MinimumDiagonal(eSymmMat);
  std::unordered_set<MyVector<Tint>> TotalListVect_set;
  auto f_insert = [&](MyMatrix<Tint> const &TheBasis,
                      MyMatrix<T> const &eSymmMatB) -> bool {
    bool test = TestCopositivityByPositivityCoeff(eSymmMatB);
    if (!test) {
      return false;
    }
    std::vector<MyVector<Tint>> TotalList =
        EnumerateShortVectorInCone_UnderPositivityCond(eSymmMat, TheBasis,
                                                       MinNorm);
    for (auto &eVect : TotalList) {
      T eNorm = EvaluationQuadForm(eSymmMat, eVect);
      if (eNorm > 0) {
        if (eNorm < MinNorm) {
          MinNorm = eNorm;
          TotalListVect_set.clear();
          TotalListVect_set.insert(eVect);
        } else {
          if (eNorm == MinNorm) {
            TotalListVect_set.insert(eVect);
          }
        }
      }
    }
    return true;
  };
  auto f_test = [&](MyMatrix<Tint> const &TheBasis,
                    MyMatrix<T> const &eSymmMatB) -> SingleTestResult<Tint> {
    return SingleTestStrictCopositivity(eSymmMat, TheBasis, eSymmMatB);
  };
#ifdef SANITY_CHECK_COPOSITIVITY
  SingleTestResult<Tint> eResult = EnumerateCopositiveShortVector_Kernel(
      eSymmMat, InitialBasis, f_insert, f_test, os);
  if (!eResult.test) {
    std::cerr << "We should not have non-copositivity in "
                 "CopositiveShortestVector\n";
    throw TerminalException{1};
  }
#endif
  int nbVect = TotalListVect_set.size();
  int iVect = 0;
  MyMatrix<Tint> SHV(nbVect, n);
  for (auto &eVect : TotalListVect_set) {
    AssignMatrixRow(SHV, iVect, eVect);
    iVect++;
  }
  return {MinNorm, std::move(SHV)};
}

template <typename T, typename Tint>
CopositivityEnumResult<Tint> EnumerateCopositiveShortVector_V2(
    MyMatrix<T> const &eSymmMat, MyMatrix<Tint> const &InitialBasis,
    T const& MaxNorm, std::ostream &os) {
  SingleTestResult<Tint> kerResult =
      SearchByZeroInKernel<T, Tint>(eSymmMat, os);
  if (!kerResult.test) {
    return {false, {}, kerResult};
  }
  std::unordered_set<MyVector<Tint>> TotalListVect_set;
  auto f_insert = [&](MyMatrix<Tint> const &TheBasis,
                      MyMatrix<T> const &eSymmMatB) -> bool {
    bool test = TestCopositivityByPositivityCoeff(eSymmMatB);
    if (!test)
      return false;
    std::vector<MyVector<Tint>> TotalList =
        EnumerateShortVectorInCone_UnderPositivityCond(eSymmMat, TheBasis,
                                                       MaxNorm);
    for (auto &eVect : TotalList)
      TotalListVect_set.insert(eVect);
    return true;
  };
  auto f_test = [&](MyMatrix<Tint> const &TheBasis,
                    MyMatrix<T> const &eSymmMatB) -> SingleTestResult<Tint> {
    return SingleTestStrictCopositivity(eSymmMat, TheBasis, eSymmMatB);
  };
  SingleTestResult<Tint> eResult = EnumerateCopositiveShortVector_Kernel(
      eSymmMat, InitialBasis, f_insert, f_test, os);
  if (!eResult.test) {
    return {false, {}, eResult};
  }
  std::vector<MyVector<Tint>> TotalListVect;
  for (auto &eVect : TotalListVect_set)
    TotalListVect.push_back(eVect);
  return {true, std::move(TotalListVect), eResult};
}

template<typename Tint>
struct ResultListCone {
  bool test;
  std::vector<MyMatrix<Tint>> ListCone;
};

template <typename T, typename Tint>
ResultListCone<Tint> EnumerateListConeCopositive(
    MyMatrix<T> const &eSymmMat, MyMatrix<Tint> const &InitialBasis,
    std::ostream &os) {
  SingleTestResult<Tint> kerResult =
      SearchByZeroInKernel<T, Tint>(eSymmMat, os);
  if (!kerResult.test) {
    return {false, {}};
  }
  std::vector<MyMatrix<Tint>> ListCone;
  auto f_insert = [&](MyMatrix<Tint> const &TheBasis,
                      MyMatrix<T> const &eSymmMatB) -> bool {
    bool test = TestCopositivityByPositivityCoeff(eSymmMatB);
    if (!test)
      return false;
    ListCone.push_back(TheBasis);
    return true;
  };
  auto f_test = [&](MyMatrix<Tint> const &TheBasis,
                    MyMatrix<T> const &eSymmMatB) -> SingleTestResult<Tint> {
    return SingleTestStrictCopositivity(eSymmMat, TheBasis, eSymmMatB);
  };
  SingleTestResult<Tint> eResult = EnumerateCopositiveShortVector_Kernel(
      eSymmMat, InitialBasis, f_insert, f_test, os);
  if (!eResult.test) {
    return {false, {}};
  }
  return {true, std::move(ListCone)};
}







template <typename T, typename Tint>
std::pair<SingleTestResult<Tint>, size_t>
TestCopositivity(MyMatrix<T> const &eSymmMat,
                 MyMatrix<Tint> const &InitialBasis, std::ostream &os) {
  int n = eSymmMat.rows();
  size_t nbCone = 0;
  auto f_insert = [&]([[maybe_unused]] MyMatrix<Tint> const &TheBasis,
                      MyMatrix<T> const &eSymmMatB) -> bool {
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++)
        if (eSymmMatB(i, j) < 0)
          return false;
    nbCone++;
    return true;
  };
  auto f_test = [&](MyMatrix<Tint> const &TheBasis,
                    MyMatrix<T> const &eSymmMatB) -> SingleTestResult<Tint> {
    return SingleTestCopositivity(eSymmMat, TheBasis, eSymmMatB);
  };
  SingleTestResult<Tint> eResult = EnumerateCopositiveShortVector_Kernel(
      eSymmMat, InitialBasis, f_insert, f_test, os);
  return {eResult, nbCone};
}

template <typename T, typename Tint>
CopositivityEnumResult<Tint> EnumerateCopositiveShortVector_V1(
    MyMatrix<T> const &eSymmMat, MyMatrix<Tint> const &InitialBasis,
    T const& MaxNorm, [[maybe_unused]] std::ostream &os) {
  CopositivityEnumResult<Tint> CopoRes =
      KernelEnumerateShortVectorInCone(eSymmMat, InitialBasis, MaxNorm);
  std::unordered_set<MyVector<Tint>> TheSet;
  // We have to do something non-clever to remove duplicates
  for (auto &eVect : CopoRes.TotalListVect) {
    TheSet.insert(eVect);
  }
  std::vector<MyVector<Tint>> TotalListRed;
  for (auto &eVect : TheSet) {
    TotalListRed.push_back(eVect);
  }
  return {CopoRes.test,
          std::move(TotalListRed),
          CopoRes.eResult};
}

template <typename T, typename Tint>
CopositivityEnumResult<Tint> EnumerateCopositiveShortVector(
    MyMatrix<T> const &eSymmMat, MyMatrix<Tint> const &InitialBasis,
    T const& MaxNorm, std::ostream &os) {
#ifdef SANITY_CHECK_COPOSITIVITY
  CopositivityEnumResult<Tint> res1 =
      EnumerateCopositiveShortVector_V1(eSymmMat, InitialBasis, MaxNorm, os);
  CopositivityEnumResult<Tint> res2 =
      EnumerateCopositiveShortVector_V2(eSymmMat, InitialBasis, MaxNorm, os);
  if (res1.test != res2.test) {
    std::cerr << "Inconsistency in computing res1.test and res2.test\n";
    throw TerminalException{1};
  }
  if (!res1.test)
    return res1;
  if (res1.TotalListVect != res2.TotalListVect) {
    std::cerr << "res1.TotalListVect, |res1.TotalListVect|="
              << res1.TotalListVect.size() << "\n";
    for (auto &eVect : res1.TotalListVect)
      WriteVectorNoDim(std::cerr, eVect);
    //
    std::cerr << "res2.TotalListVect, |res2.TotalListVect|="
              << res2.TotalListVect.size() << "\n";
    for (auto &eVect : res2.TotalListVect)
      WriteVectorNoDim(std::cerr, eVect);
    //
    throw TerminalException{1};
  }
  return res1;
#else
  return EnumerateCopositiveShortVector_V2<T, Tint>(eSymmMat, InitialBasis,
                                                    MaxNorm, os);
#endif
}

template <typename T, typename Tint>
Tshortest<T, Tint>
CopositiveShortestVector_V1(MyMatrix<T> const &eSymmMat,
                              MyMatrix<Tint> const &InitialBasis,
                              std::ostream &os) {
  int n = eSymmMat.rows();
  T MaxNorm = MinimumDiagonal(eSymmMat);
  CopositivityEnumResult<Tint> CopoRes =
      EnumerateCopositiveShortVector(eSymmMat, InitialBasis, MaxNorm, os);
#ifdef SANITY_CHECK_COPOSITIVITY
  if (!CopoRes.test) {
    std::cerr << "Inconsistency in result\n";
    throw TerminalException{1};
  }
#endif
  int nbVect = CopoRes.TotalListVect.size();
  MyMatrix<int> SHV(nbVect, n);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    MyVector<Tint> V = CopoRes.TotalListVect[iVect];
    AssignMatrixRow(SHV, iVect, V);
  }
  return SelectShortestVector<T, Tint>(eSymmMat, SHV);
}

// clang-format off
#endif  // SRC_COPOS_COPOSITIVITY_H_
// clang-format on
