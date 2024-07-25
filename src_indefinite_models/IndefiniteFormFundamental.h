// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_MODELS_FUNDAMENTALS_H_
#define SRC_INDEFINITE_MODELS_FUNDAMENTALS_H_

template<typename T>
struct AttackScheme {
  int h;
  MyMatrix<T> mat;
};

template<typename T>
AttackScheme<T> INDEF_FORM_GetAttackScheme(MyMatrix<T> const& Qmat) {
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(Qmat);
  int nbPlus=DiagInfo.nbPlus;
  int nbMinus=DiagInfo.nbMinus;
  if (nbMinus < nbPlus) {
    MyMatrix<T> RetMat = -Qmat;
    return {nbMinus, RetMat};
  } else {
    return {nbPlus, Qmat};
  }
}

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

template<typename T>
struct INDEF_InvariantQ {
  int n;
  T eDet;
  int nbPlus;
  int nbMinus;
  int nbZero;
  bool IsEven;
};

template<typename T>
INDEF_InvariantQ<T> INDEF_FORM_Invariant(MyMatrix<T> const& Qmat) {
  using Tring = typename underlying_ring<T>::ring_type;
  int n = Qmat.rows();
  MyMatrix<T> NSP = NullspaceIntMat(Qmat);
  MyMatrix<T> TheCompl = SubspaceCompletionInt(NSP, n);
  MyMatrix<T> GramRed = TheCompl * Qmat * TheCompl.transpose();
  T eDet = DeterminantMat(GramRed);
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(Qmat);
  int nbPlus = DiagInfo.nbPlus;
  int nbMinus = DiagInfo.nbMinus;
  int nbZero = DiagInfo.nbZero;
  bool IsEven = INDEF_FORM_IsEven(Qmat);
  return {n, eDet, nbPlus, nbMinus, nbZero, IsEven};
}







// clang-format off
#endif  // SRC_INDEFINITE_MODELS_COMBINEDALGORITHMS_H_
// clang-format on

