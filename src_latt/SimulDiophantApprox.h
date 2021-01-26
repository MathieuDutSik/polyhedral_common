#ifndef INCLUDE_SIMULTANEOUS_DIOPHANTINE_APPROX
#define INCLUDE_SIMULTANEOUS_DIOPHANTINE_APPROX

#include "LatticeDefinition.h"


template<typename Tint>
struct DiopantResult {
  MyVector<Tint> Numerators;
  Tint Denominator;
};



template<typename T, typename Tint>
DiophantResult<Tint> SimultaneousDiophantineApproximation(MyVector<T> const& V, T const& epsilon)
{
  //
  // Building the Gram matrix according to Lenstra, Lenstra, Lovasz, Factoring Polynomials with Rational Coefficients
  // Proposition 1.39, Page 525.
  //
  int n = V.size();
  int totdim = n * (n+1) / 2;
  MyMatrix<T> GramMat = ZeroMatrix<T>(n+1, n+1);
  T eCoeff = 0;
  for (int i=0; i<n; i++) {
    GramMat(i,i) = 1;
    T eV = - V(i);
    GramMat(n,i) = eV;
    GramMat(i,n) = eV;
    eCoeff += eV * eV;
  }
  auto IntPow=[](T eVal, int nb) -> T {
    T retV = 1;
    for (int idx=0; idx<nb; idx++)
      retV *= eVal;
    return retV;
  };
  T eVal = 2;
  eCoeff += IntPow(epsilon, n) / IntPow(eVal, totdim);
  GramMat(n,n) = eCoeff;
  //
  // Computing the LLL reduction
  //
  LLLreduction<Tmat,Tint> LLLinfo = LLLreducedBasis<T, Tint>(GramMat);
  //
  // Extrqcting the approimation
  //
  MyVector<Tint> V(n);
  for (int i=0; i<n; i++)
    V(i) = LLLinfo.Pmat(0, i);
  Tint q = LLLinfo.Pmat(0, n);
  return {std::move(V), q};
}



#endif
