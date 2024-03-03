// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TEMP_POSITIVITY_H_
#define SRC_LATT_TEMP_POSITIVITY_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixFLT.h"
#include "MAT_MatrixInt.h"
#include <vector>
// clang-format on

// This code does the following:
// * diagonalization of symmetric matrices
// * Finding a vector of positive norm in a vector.

#ifdef DEBUG
#define DEBUG_POSITIVITY
#endif


template <typename T> T MinimumDiagonal(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  T MinNorm = eMat(0, 0);
  for (int i = 1; i < n; i++) {
    T eVal = eMat(i, i);
    if (eVal < MinNorm)
      MinNorm = eVal;
  }
  return MinNorm;
}

template <typename T> T MaximumDiagonal(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  T MaxNorm = eMat(0, 0);
  for (int i = 1; i < n; i++) {
    T eVal = eMat(i, i);
    if (eVal > MaxNorm)
      MaxNorm = eVal;
  }
  return MaxNorm;
}

template <typename T> bool IsPositiveDefinite(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  for (int siz = 1; siz <= n; siz++) {
    MyMatrix<T> eMatRed(siz, siz);
    for (int i = 0; i < siz; i++)
      for (int j = 0; j < siz; j++)
        eMatRed(i, j) = eMat(i, j);
    T eDet = DeterminantMat(eMatRed);
    if (eDet <= 0)
      return false;
  }
  return true;
}

template <typename T> MyMatrix<T> AnLattice(int const &n) {
  MyMatrix<T> eMat(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      eMat(i, j) = 0;
  for (int i = 0; i < n; i++)
    eMat(i, i) = 2;
  for (int i = 1; i < n; i++) {
    eMat(i, i - 1) = -1;
    eMat(i - 1, i) = -1;
  }
  return eMat;
}

template <typename T>
MyVector<T> FindNonIsotropeVector(MyMatrix<T> const &SymMat) {
  int n = SymMat.rows();
  MyVector<T> eVect = ZeroVector<T>(n);
  for (int i = 0; i < n; i++)
    if (SymMat(i, i) != 0) {
      eVect(i) = 1;
      return eVect;
    }
  for (int i = 0; i < n - 1; i++)
    for (int j = i + 1; j < n; j++)
      if (SymMat(i, j) != 0) {
        eVect(i) = 1;
        eVect(j) = 1;
        return eVect;
      }
  std::cerr << "SymMat=\n";
  WriteMatrix(std::cerr, SymMat);
  std::cerr << "Clear error in FindNonIsotropeVector\n";
  throw TerminalException{1};
}

template <typename T> struct DiagSymMat {
  MyMatrix<T> Transform;
  MyMatrix<T> RedMat;
  int nbZero;
  int nbPlus;
  int nbMinus;
};

template <typename T>
DiagSymMat<T>
DiagonalizeNonDegenerateSymmetricMatrix(MyMatrix<T> const &SymMat) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
#ifdef DEBUG_POSITIVITY
  std::cerr << "Beginning of DiagonalizeNonDegenerateSymmetricMatrix\n";
  WriteMatrix(std::cerr, SymMat);
#endif
  int n = SymMat.rows();
  if (n == 0) {
    return {{}, {}, 0, 0, 0};
  }
  std::vector<MyVector<T>> ListVect;
  for (int i = 0; i < n; i++) {
    MyMatrix<T> BasisOrthogonal;
    if (i == 0) {
      BasisOrthogonal = IdentityMat<T>(n);
    } else {
#ifdef DEBUG_POSITIVITY
      std::cerr << "|ListVect|=" << ListVect.size() << "\n";
#endif
      MyMatrix<T> TheBasis = MatrixFromVectorFamily(ListVect);
      MyMatrix<T> eProd = SymMat * TheBasis.transpose();
      BasisOrthogonal = NullspaceMat(eProd);
    }
    MyMatrix<T> RedInBasis =
        BasisOrthogonal * SymMat * BasisOrthogonal.transpose();
    MyVector<T> V = FindNonIsotropeVector(RedInBasis);
    MyVector<T> Vins = (BasisOrthogonal.transpose()) * V;
    ListVect.push_back(Vins);
  }
  MyMatrix<T> TheBasis = MatrixFromVectorFamily(ListVect);
  MyMatrix<T> RedMat = TheBasis * SymMat * TheBasis.transpose();
#ifdef DEBUG_POSITIVITY
  std::cerr << "RedMat=\n";
  WriteMatrix(std::cerr, RedMat);
#endif
  int nbPlus = 0;
  int nbMinus = 0;
  int nbZero = 0;
  for (int i = 0; i < n; i++) {
    if (RedMat(i, i) > 0)
      nbPlus++;
    if (RedMat(i, i) < 0)
      nbMinus++;
  }
#ifdef DEBUG_POSITIVITY
  std::cerr << "nbPlus=" << nbPlus << " nbMinus=" << nbMinus << "\n";
#endif
  return {TheBasis, RedMat, nbZero, nbPlus, nbMinus};
}

template <typename T> struct NSPreduction {
  MyMatrix<T> RedMat;
  MyMatrix<T> Transform;
  MyMatrix<T> NonDegenerate;
};

template <typename T>
MyMatrix<T> SymmetricExtractSubMatrix(MyMatrix<T> const &M,
                                      std::vector<int> const &S) {
  int siz = S.size();
  MyMatrix<T> Mret(siz, siz);
  for (int i = 0; i < siz; i++)
    for (int j = 0; j < siz; j++) {
      int i2 = S[i];
      int j2 = S[j];
      Mret(i, j) = M(i2, j2);
    }
  return Mret;
}

template <typename T>
NSPreduction<T> NullspaceReduction(MyMatrix<T> const &SymMat) {
  int n = SymMat.rows();
  MyMatrix<T> PreTransfMat = NullspaceMat(SymMat);
  int DimKern = PreTransfMat.rows();
  MyMatrix<T> TransfMat = ZeroMatrix<T>(n, n);
  for (int iRow = 0; iRow < DimKern; iRow++)
    for (int i = 0; i < n; i++)
      TransfMat(iRow, i) = PreTransfMat(iRow, i);
  std::vector<int> eSet = ColumnReductionSet(PreTransfMat);
  std::vector<int> totSet(n);
  for (int i = 0; i < n; i++)
    totSet[i] = i;
  std::vector<int> diffSet = DifferenceVect(totSet, eSet);
  //  std::cerr << "|eSet|=" << eSet.size() << " |diffSet|=" << diffSet.size()
  //  << "\n";
  int idx = DimKern;
  for (auto &eCol : diffSet) {
    TransfMat(idx, eCol) = 1;
    idx++;
  }
  /*  std::cerr << "TransfMat=\n";
      WriteMatrix(std::cerr, TransfMat);*/
  MyMatrix<T> SymMat2 = TransfMat * SymMat * TransfMat.transpose();
  std::vector<int> gSet;
  for (int i = DimKern; i < n; i++)
    gSet.push_back(i);
  MyMatrix<T> MatExt = SymmetricExtractSubMatrix(SymMat2, gSet);
  return {SymMat2, TransfMat, MatExt};
}

template <typename T>
DiagSymMat<T> DiagonalizeSymmetricMatrix(MyMatrix<T> const &SymMat) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  //  std::cerr << "DiagonalizeSymmetricMatrix\n";
  //  WriteMatrix(std::cerr, SymMat);
  //  std::cerr << "DiagonalizeSymmetricMatrix, RankMat(SymMat)=" <<
  //  RankMat(SymMat) << "\n";
  int n1 = SymMat.rows();
  NSPreduction<T> NSP1 = NullspaceReduction(SymMat);
  MyMatrix<T> RMat1 = NSP1.Transform;
  MyMatrix<T> SymMat2 = NSP1.NonDegenerate;
  //  std::cerr << "DiagonalizeSymmetricMatrix, RankMat(SymMat2)=" <<
  //  RankMat(SymMat2) << "\n";
  DiagSymMat<T> NSP2 = DiagonalizeNonDegenerateSymmetricMatrix(SymMat2);
  int n2 = SymMat2.rows();
  MyMatrix<T> RMat2 = ZeroMatrix<T>(n1, n1);
  for (int i = 0; i < n2; i++)
    for (int j = 0; j < n2; j++)
      RMat2(i + n1 - n2, j + n1 - n2) = NSP2.Transform(i, j);
  for (int i = 0; i < n1 - n2; i++)
    RMat2(i, i) = 1;
  MyMatrix<T> RMat = RMat2 * RMat1;
  MyMatrix<T> RedMat = RMat * SymMat * RMat.transpose();
  int nbPlus = 0;
  int nbMinus = 0;
  int nbZero = 0;
  for (int i = 0; i < n1; i++) {
    if (RedMat(i, i) > 0)
      nbPlus++;
    if (RedMat(i, i) < 0)
      nbMinus++;
    if (RedMat(i, i) == 0)
      nbZero++;
  }
  return {RMat, RedMat, nbZero, nbPlus, nbMinus};
}

template <typename T>
MyVector<T> GetPositiveNormVector(MyMatrix<T> const &SymMat) {
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(SymMat);
  MyMatrix<T> const &Transform = DiagInfo.Transform;
  MyMatrix<T> const &RedMat = DiagInfo.RedMat;
  int n = SymMat.rows();
  for (int i = 0; i < n; i++) {
    if (RedMat(i, i) > 0) {
      MyVector<T> eVect = GetMatrixRow(Transform, i);
      T norm = eVect.dot(SymMat * eVect);
      if (norm != RedMat(i, i)) {
        std::cerr << "We have a consistency error\n";
        throw TerminalException{1};
      }
      return eVect;
    }
  }
  std::cerr << "Failed to find a negative norm vector\n";
  throw TerminalException{1};
}


template<typename T, typename Tint, typename Ttest>
MyVector<Tint> GetIntegralPositiveVector_family(std::vector<MyVector<Ttest>> const& ListVect, MyMatrix<T> const& M) {
  int scal = 1;
  MyVector<Tint> V_ret(n);
  while(true) {
    for (auto & eV : ListVect) {
      for (int i=0; i<n; ui++) {
        Ttest val_d = scal * eV(i);
        Tint val = UniversalNearestScalarInteger<Tint, Ttest>(val_d);
        V_ret(i) = val;
      }
      T eNorm = EvaluationQuadForm<T,Tint>(M, V_ret);
      if (eNorm > 0) {
        return V_ret;
      }
    }
#ifdef DEBUG_POSITIVITY
    os << "POS: GetPositiveNormVector_family |ListVect|=" << ListVect.size() << " scal=" << scal << "\n";
#endif
    scal += 1;
  }
}




template <typename T, typename Tint>
MyVector<Tint> GetIntegralPositiveVector_diag(MyMatrix<T> const &M) {
  int n = M.rows();
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(M);
  MyMatrix<T> const &Transform = DiagInfo.Transform;
  MyMatrix<T> const &RedMat = DiagInfo.RedMat;
  std::vector<MyVector<T>> ListVect;
  for (int i = 0; i < n; i++) {
    if (RedMat(i, i) > 0) {
      MyVector<T> eVect = GetMatrixRow(Transform, i);
      T sum = 0;
      for (int j=0; j<n; j++) {
        sum += T_abs(eVect(j));
      }
      MyVector<T> fVect = eVect / sum;
      ListVect.push_back(fVect);
    }
  }
  return GetIntegralPositiveVector_family<T,Tint,T>(ListVect, M);
}


template<typename T, typename Tint>
MyVector<Tint> GetIntegralPositiveVector_eigen(MyMatrix<T> const& M, [[maybe_unused]] std::ostream & os) {
  int n = M.rows();
  MyMatrix<double> M_double = UniversalMatrixConversion<double,T>(M);
  Eigen::SelfAdjointEigenSolver<MyMatrix<double>> eig(M_double);
  MyVector<T> ListEig = eig.eigenvalues();
  MyMatrix<T> ListVect = eig.eigenvectors();
  std::vector<MyVector<double>> ListEigVect;
  for (int i=0; i<n; i++) {
    if (ListEig(i) > 0) {
      MyVector<double> V(n);
      for (int j=0; j<n; j++) {
        V(j) = ListVect(i, j);
      }
      ListEigVect.push_back(V);
    }
  }
  return GetIntegralPositiveVector_family<T,Tint,double>(ListEigVect, M);
}

template <typename T>
std::vector<MyVector<T>> GetSetNegativeOrZeroVector(MyMatrix<T> const &SymMat) {
  DiagSymMat<T> eRecDiag = DiagonalizeSymmetricMatrix(SymMat);
  std::vector<MyVector<T>> TheSet;
  int n = SymMat.rows();
  for (int i = 0; i < n; i++)
    if (eRecDiag.RedMat(i, i) <= 0) {
      MyVector<T> eVect = ZeroVector<T>(n);
      eVect(i) = 1;
      MyVector<T> fVect = (eRecDiag.Transform.transpose()) * eVect;
      T eEval = EvaluationQuadForm(SymMat, fVect);
      if (eEval > 0) {
        std::cerr << "Big bad error\n";
        throw TerminalException{1};
      }
      T eMax = fVect.maxCoeff();
      MyVector<T> gVect = fVect / eMax;
      TheSet.push_back(gVect);
    }
  return TheSet;
}

template<typename T>
std::vector<T> GetLineVector(MyMatrix<T> const& M) {
  int n = M.rows();
  int dim = n * (n+1) / 2;
  std::vector<T> V(dim);
  int pos = 0;
  for (int i=0; i<n; i++) {
    for (int j=i; j<n; j++) {
      if (i == j) {
        V[pos] = M(i, i);
      } else {
        V[pos] = 2 * M(i, j);
      }
      pos++;
    }
  }
  return V;
}

template<typename T>
T EvaluateLineVector(std::vector<T> const& V_mat, MyVector<T> const& V) {
  T sum = 0;
  T pSum = 0;
  int n = V.size();
  size_t pos = 0;
  for (int i=0; i<n; i++) {
    pSum = 0;
    for (int j=i; j<n; j++) {
      pSum += V_mat[pos] * V[j];
      pos++;
    }
    sum += pSum * V[i];
  }
  return sum;
}

template<typename T>
T EvaluateLineVectorShift(std::vector<T> const& V_mat, MyVector<T> const& V, int shift) {
  T sum = 0;
  T pSum = 0;
  int n = V.size();
  size_t pos = 0;
  for (int i=shift; i<n; i++) {
    pSum = 0;
    for (int j=i; j<n; j++) {
      pSum += V_mat[pos] * V[j];
      pos++;
    }
    sum += pSum * V[i];
  }
  return sum;
}


// clang-format off
#endif  // SRC_LATT_TEMP_POSITIVITY_H_
// clang-format on
