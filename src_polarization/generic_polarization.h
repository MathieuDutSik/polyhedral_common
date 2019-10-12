#ifndef INCLUDE_GENERIC_POLARIZATION
#define INCLUDE_GENERIC_POLARIZATION

#include "Temp_common.h"
#include "autodiff.h"


dual2nd Compute_Polarization(int const& n, std::vector<int> const& Ranges, double const& MaxNorm, std::vector<dual2nd> const& c, std::vector<dual2nd> const& A, dualnd std::function<dualnd(dualnd const&)> const& f)
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
      ListI.push_back(i);
      ListJ.push_back(j);
      ListFactor.push_back(eFactor);
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
      int eFactor=ListFactor[idx];
      int i=ListI[idx];
      int j=ListJ[idx];
      Aval += eFactor * (c[i] - xVect[i]) * A[idx] * (c[j] - xVect[j]);
    }
    RetVal += f(Aval);
    int test = Blk.IncrementShow();
    if (test == -1)
      break;
  }
  return RetVal;
}

#endif
