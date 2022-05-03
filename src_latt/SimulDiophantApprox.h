#ifndef SRC_LATT_SIMULDIOPHANTAPPROX_H_
#define SRC_LATT_SIMULDIOPHANTAPPROX_H_

#include "LatticeDefinitions.h"
#include <string>
#include <utility>

template <typename Tint> struct DiophantResult {
  MyVector<Tint> Numerators;
  Tint Denominator;
};

template <typename T, typename Tint>
T ComputeDiophantinePenalty(MyVector<T> const &V,
                            DiophantResult<Tint> const &Res) {
  int n = V.size();
  T maxErr = 0;
  T Denom_T = UniversalScalarConversion<T, Tint>(Res.Denominator);
  for (int i = 0; i < n; i++) {
    T eNum_T = UniversalScalarConversion<T, Tint>(Res.Numerators(i));
    T eFrac = eNum_T / Denom_T;
    T diff = eFrac - V(i);
    T err = T_abs(diff);
    maxErr = T_max(maxErr, err);
  }
  return maxErr;
}

template <typename Tint>
std::string GapStringDiophantineApprox(DiophantResult<Tint> const &Res) {
  std::string strO = "[ ";
  int n = Res.Numerators.size();
  for (int i = 0; i < n; i++) {
    if (i > 0)
      strO += ", ";
    strO += std::to_string(Res.Numerators(i)) + "/" +
            std::to_string(Res.Denominator);
  }
  strO += " ]";
  return strO;
}

template <typename T, typename Tint>
DiophantResult<Tint> SimultaneousDiophantineApproximation(MyVector<T> const &V,
                                                          T const &epsilon) {
  //
  // Building the Gram matrix according to Lenstra, Lenstra, Lovasz, Factoring
  // Polynomials with Rational Coefficients Proposition 1.39, Page 525.
  //
  int n = V.size();
  int totdim = n * (n + 1) / 2;
  MyMatrix<T> GramMat = ZeroMatrix<T>(n + 1, n + 1);
  T eCoeff = 0;
  for (int i = 0; i < n; i++) {
    GramMat(i, i) = 1;
    T eV = -V(i);
    GramMat(n, i) = eV;
    GramMat(i, n) = eV;
    eCoeff += eV * eV;
  }
  auto IntPow = [](T eVal, int nb) -> T {
    T retV = 1;
    for (int idx = 0; idx < nb; idx++)
      retV *= eVal;
    return retV;
  };
  T eVal = 2;
  eCoeff += IntPow(epsilon, n) / IntPow(eVal, totdim);
  GramMat(n, n) = eCoeff;
  //
  // Computing the LLL reduction
  //
  LLLreduction<T, Tint> LLLinfo = LLLreducedBasis<T, Tint>(GramMat);
  //
  // Extrqcting the approimation
  //
  MyVector<Tint> Vapprox(n);
  for (int i = 0; i < n; i++)
    Vapprox(i) = LLLinfo.Pmat(0, i);
  Tint q = LLLinfo.Pmat(0, n);
  if (q > 0) {
    return {std::move(Vapprox), q};
  } else {
    return {std::move(-Vapprox), -q};
  }
}

// clang-format off
#endif  // SRC_LATT_SIMULDIOPHANTAPPROX_H_
// clang-format on
