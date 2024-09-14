// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_POSITIVENEGATIVE_H_
#define SRC_INDEFINITE_POSITIVENEGATIVE_H_

// clang-format off
#include "LatticeStabEquiCan.h"
// clang-format on

// In some contexts, we need to deal with positive-definite or negative-definite
// matrices.

template <typename T, typename Tint, typename Tgroup>
std::vector<MyVector<Tint>>
INDEF_FORM_GetOrbitRepresentative_PosNeg(MyMatrix<T> const &Q, T const &X,
                                         std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Tgroup::Tidx;
  auto get_orbit_representatives =
      [&](MyMatrix<T> const &Qin) -> std::vector<MyVector<Tint>> {
    MyMatrix<Tint> SHV = EnumerateVectorsFixedNorm<T, Tint>(Qin, X, os);
    size_t len = SHV.rows();
    std::vector<MyVector<Tint>> ListVect;
    std::unordered_map<MyVector<Tint>, size_t> MapVect;
    for (int i = 0; i < SHV.rows(); i++) {
      MyVector<Tint> eV = GetMatrixRow(SHV, i);
      MapVect[eV] = i + 1;
      ListVect.push_back(eV);
    }
    std::vector<MyMatrix<Tint>> LGen =
        ArithmeticAutomorphismGroup<T, Tint>(Qin, os);
    std::vector<Telt> ListPerm;
    for (auto &eGen : LGen) {
      std::vector<Tidx> eList(len);
      for (size_t i = 0; i < len; i++) {
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
    for (auto &ePos : LPos) {
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

template <typename T, typename Tint>
std::vector<MyMatrix<Tint>>
INDEF_FORM_AutomorphismGroup_PosNeg(MyMatrix<T> const &Q, std::ostream &os) {
  DiagSymMat<T> DSM = DiagonalizeNonDegenerateSymmetricMatrix(Q);
  if (DSM.nbPlus == 0 && DSM.nbZero == 0) {
    MyMatrix<T> Qneg = -Q;
    return ArithmeticAutomorphismGroup<T, Tint>(Qneg, os);
  }
  if (DSM.nbMinus == 0 && DSM.nbZero == 0) {
    return ArithmeticAutomorphismGroup<T, Tint>(Q, os);
  }
  std::cerr << "Failed to find a matching entry for PosNeg\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
std::optional<MyMatrix<Tint>>
INDEF_FORM_TestEquivalence_PosNeg(MyMatrix<T> const &Q1, MyMatrix<T> const &Q2,
                                  std::ostream &os) {
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
    return ArithmeticEquivalence<T, Tint>(Qneg1, Qneg2, os);
  }
  if (DSM1.nbMinus == 0 && DSM1.nbZero == 0) {
    return ArithmeticEquivalence<T, Tint>(Q1, Q2, os);
  }
  std::cerr << "Failed to find a matching entry for PosNeg\n";
  throw TerminalException{1};
}

template <typename T> bool INDEF_FORM_IsPosNeg(MyMatrix<T> const &M) {
  DiagSymMat<T> DSM = DiagonalizeNonDegenerateSymmetricMatrix(M);
  if (DSM.nbZero == 0) {
    if (DSM.nbPlus == 0 || DSM.nbMinus == 0) {
      return true;
    }
  }
  return false;
}

// clang-format off
#endif  // SRC_INDEFINITE_POSITIVENEGATIVE_H_
// clang-format on
