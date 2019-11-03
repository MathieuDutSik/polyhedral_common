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
    RetVal += f(Aval);
    int test = Blk.IncrementShow();
    if (test == -1)
      break;
  }
  return RetVal;
}





// Computation using the Senergy.
dual2nd Senergy(int const& n, std::vector<int> const& Ranges, std::vector<dual2nd> const& A, std::function<dual2nd(dual2nd const&)> const& f)
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
    int sumAbs = 0;
    for (int i=0; i<n; i++) {
      int eVal = eVect[i] - Ranges[i];
      sumAbs += T_abs(eVal);
      xVect[i] = eVal;
    }
    if (sumAbs > 0) {
      //
      dual2nd Aval=0;
      for (int idx=0; idx<TotDim; idx++) {
        double eFactor=double(ListFactor[idx]);
        int i=ListI[idx];
        int j=ListJ[idx];
        Aval += eFactor * xVect[i] * A[idx] * xVect[j];
      }
      RetVal += f(Aval);
    }
    int test = Blk.IncrementShow();
    if (test == -1)
      break;
  }
  return RetVal;
}

/*
MyMatrix<double> NewtonIterationOptimization(MyMatrix<double> const& Ainit, std::function<dual2nd(dual2nd const&)> const& f)
{
  int n=Ainit.rows();
  std::vector<int> Ranges(n, 10);
  while(true) {
    
    Senergy(n, Ranges, std::vector<dual2nd> const& A, f);

    
  }
}
*/




#endif
