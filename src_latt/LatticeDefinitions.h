// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEDEFINITIONS_H_
#define SRC_LATT_LATTICEDEFINITIONS_H_

// clang-format off
#include "Boost_bitset.h"
#include "COMB_Combinatorics.h"
#include "ClassicLLL.h"
#include "MAT_Matrix.h"
#include <utility>
#include <string>
#include <vector>
// clang-format on

template <typename T, typename Tint> struct Tshortest {
  T min;
  MyMatrix<Tint> SHV;
};

template <typename T, typename Tint>
Tshortest<T,Tint> shortest_get_half(Tshortest<T,Tint> const& rec_shv) {
  int n_row = rec_shv.SHV.rows();
  int n_col = rec_shv.SHV.cols();
  int n_row2 = n_row / 2;
  MyMatrix<Tint> SHV(n_row2, n_col);
  for (int i_row2=0; i_row2<n_row2; i_row2++) {
    for (int i_col=0; i_col<n_col; i_col++) {
      SHV(i_row2, i_col) = rec_shv.SHV(2*i_row2, i_col);
    }
  }
  return {rec_shv.min, std::move(SHV)};
}


template <typename T, typename Tint>
Tshortest<T, Tint> SelectShortestVector(MyMatrix<T> const &eMat,
                                        MyMatrix<Tint> const &SHV) {
  int n = eMat.rows();
  int nbRow = SHV.rows();
  std::vector<int> ListStatus(nbRow, 0);
  int nbShort = 0;
  T MinNorm = -1;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<Tint> eVect = SHV.row(iRow);
    T eNorm = EvaluationQuadForm(eMat, eVect);
    if (eNorm > 0) {
      if (nbShort == 0) {
        MinNorm = eNorm;
        nbShort = 1;
        ListStatus[iRow] = 1;
      } else {
        if (eNorm < MinNorm) {
          MinNorm = eNorm;
          nbShort = 1;
          for (int jRow = 0; jRow < iRow; jRow++)
            ListStatus[jRow] = 0;
          ListStatus[iRow] = 1;
        } else {
          if (eNorm == MinNorm) {
            nbShort++;
            ListStatus[iRow] = 1;
          }
        }
      }
    }
  }
  MyMatrix<Tint> TheSHV(nbShort, n);
  int iShort = 0;
  for (int iRow = 0; iRow < nbRow; iRow++)
    if (ListStatus[iRow] == 1) {
      for (int i = 0; i < n; i++) {
        TheSHV(iShort, i) = SHV(iRow, i);
      }
      iShort++;
    }
  //  std::cerr << "MinNorm=" << MinNorm << "\n";
  return {MinNorm, std::move(TheSHV)};
}

template <typename T, typename Tint> struct resultCVP {
  T TheNorm;
  MyMatrix<Tint> ListVect;
};

template <typename Tint>
void EnumerateOrbitPrimitiveVector(
    std::vector<MyMatrix<int>> const &ListMat, int const &pPrime,
    int const &dim,
    std::function<void(std::vector<int> const &, int const &)> const &FCT) {
  Tint expo(1);
  for (int i = 0; i < dim; i++) {
    expo *= pPrime;
  }
  Tint TotalNb = (expo - 1) / (pPrime - 1);
  Face StatusVect(TotalNb);
  for (Tint i = 0; i < TotalNb; i++) {
    StatusVect[i] = 1;
  }
  //
  // The modulo p Operations
  //
  auto Canonicalization = [&](int const &x) -> int {
    int res = x % pPrime;
    return res;
  };
  auto ProdModP = [&](int const &a, int const &b) -> int {
    return Canonicalization(a * b);
  };
  std::vector<int> ListInverse(pPrime);
  for (int i = 0; i < pPrime; i++) {
    int eInv = -400;
    for (int j = 0; j < pPrime; j++) {
      int res = ProdModP(i, j);
      if (res == 1)
        eInv = j;
    }
    ListInverse[i] = eInv;
  }
  auto PrimitivizeVector = [&](std::vector<int> &V) -> void {
    int ifound = -1;
    for (int i = 0; i < dim; i++) {
      if (V[i] != 0) {
        ifound = i;
      }
    }
    int eVal = V[ifound];
    int eInv = ListInverse[eVal];
    for (int i = 0; i <= ifound; i++) {
      V[i] = ProdModP(eInv, V[i]);
    }
  };
  auto ImagePrimitiveVector =
      [&](MyMatrix<int> const &M,
          std::vector<int> const &V) -> std::vector<int> {
    std::vector<int> Vret(dim);
    for (int i = 0; i < dim; i++) {
      int sum = 0;
      for (int j = 0; j < dim; j++)
        sum += M(j, i) * V[j];
      int res = Canonicalization(sum);
      Vret[i] = res;
    }
    PrimitivizeVector(Vret);
    return Vret;
  };
  auto ImageNumber = [&](MyMatrix<int> const &M, Tint const &valIn) -> Tint {
    std::vector<int> Vin = ConvertNumberToPrimitiveVector(valIn, pPrime, dim);
    int NewVal = ConvertPrimitiveVectorToNumber<Tint>(Vin, pPrime);
    if (NewVal != valIn) {
      std::cerr << "Inconsistency to be solved\n";
      throw TerminalException{1};
    }
    std::vector<int> Vout = ImagePrimitiveVector(M, Vin);
    return ConvertPrimitiveVectorToNumber<Tint>(Vout, pPrime);
  };
  //
  // Combinatorial enumeration code
  //
  std::vector<std::vector<int>> ListRepresentative;
  Tint cnt = TotalNb;
  StackStorage<Tint> ActivePoints;
  int ThePos = 0;
  auto GetPosition = [&]() -> int {
    while (true) {
      if (StatusVect[ThePos] == 1)
        return ThePos;
      ThePos++;
    }
  };
  int nbOrbit = 0;
  while (true) {
    if (cnt == 0)
      break;
    int OrbSize = 0;
    int ePos = GetPosition();
    std::vector<int> eRepr = ConvertNumberToPrimitiveVector(ePos, pPrime, dim);
    auto FuncInsert = [&](Tint const &val) -> void {
      if (StatusVect[val] == 1) {
        StatusVect[val] = 0;
        ActivePoints.push_back(val);
        cnt--;
        OrbSize++;
      }
    };
    FuncInsert(ePos);
    while (true) {
      if (ActivePoints.size() == 0)
        break;
      Tint ePt = ActivePoints.pop();
      for (auto &M : ListMat) {
        Tint fPt = ImageNumber(M, ePt);
        FuncInsert(fPt);
      }
    }
    FCT(eRepr, OrbSize);
    std::cerr << "Find an orbit of size " << OrbSize << "\n";
    nbOrbit++;
  }
  std::cerr << "nbOrbit=" << nbOrbit << "\n";
}

// clang-format off
#endif  // SRC_LATT_LATTICEDEFINITIONS_H_
// clang-format on
