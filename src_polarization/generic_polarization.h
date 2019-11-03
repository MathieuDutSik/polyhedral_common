#ifndef INCLUDE_GENERIC_POLARIZATION
#define INCLUDE_GENERIC_POLARIZATION

#include "Temp_common.h"
#include "autodiff.h"
#include "COMB_Combinatorics.h"

// Maybe use MaxNorm?
dual2nd Compute_Polarization(int const& n, std::vector<int> const& Ranges, std::vector<dual2nd> const& cA, std::function<dual2nd(dual2nd const&)> const& f)
{
  std::vector<int> Ranges_Ext(n);
  for (int i=0; i<n; i++)
    Ranges_Ext[i] = 1 + 2*Ranges[i];
  int TotDim=n*(n+1)/2;
  std::vector<int> ListI(TotDim), ListJ(TotDim), ListFactor(TotDim);
  int idx=0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<=i; j++) {
      int eFactor = i == j ? 1 : 2;
      //      std::cerr << "i=" << i << " j=" << j << " eF=" << eFactor << "\n";
      ListI[idx] = i;
      ListJ[idx] = j;
      ListFactor[idx] = eFactor;
      idx++;
    }
  }
  //
  BlockIterationMultiple Blk(Ranges_Ext);
  dual2nd RetVal=0;

  while(true) {
    std::vector<int> eVect = Blk.GetVect();
    std::vector<int> xVect(n);
    for (int i=0; i<n; i++)
      xVect[i] = eVect[i] - Ranges[i];
    //
    dual2nd Aval=0;
    for (int idx=0; idx<TotDim; idx++) {
      double eFactor=double(ListFactor[idx]);
      int i=ListI[idx];
      int j=ListJ[idx];
      Aval += eFactor * (cA[i] - xVect[i]) * cA[n+idx] * (cA[j] - xVect[j]);
    }
    //    Aval += A[0];
    //    RetVal += f(Aval);
    RetVal += f(Aval);
    int test = Blk.IncrementShow();
    if (test == -1)
      break;
  }

  /*
  while(true) {
    std::vector<int> eVect = Blk.GetVect();
    std::vector<int> xVect(n);
    for (int i=0; i<n; i++)
      xVect[i] = eVect[i] - Ranges[i];
    //
    dual2nd Aval=0;
    for (int idx=0; idx<TotDim; idx++) {
      double eFactor=double(ListFactor[idx]);
      int i=ListI[idx];
      int j=ListJ[idx];
      //      Aval += eFactor * (c[i] - xVect[i]) * A[idx] * (c[j] - xVect[j]);
      Aval += eFactor * A[idx];
    }
    //    RetVal += f(Aval);
    double fact = 2;
    RetVal += fact * A[0];
    int test = Blk.IncrementShow();
    if (test == -1)
      break;
  }
*/
  return RetVal;
}





#endif
