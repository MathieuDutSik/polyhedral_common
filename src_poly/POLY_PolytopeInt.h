#ifndef TEMP_POLYTOPE_INT_INCLUDE
#define TEMP_POLYTOPE_INT_INCLUDE

#include "COMB_Combinatorics.h"


template<typename T>
std::vector<MyVector<T>> GetListIntegralPoint(MyMatrix<T> const& FAC, MyMatrix<T> const& EXT)
{
  int nbVert=EXT.rows();
  int n=EXT.cols();
  int dim=n-1;
  MyMatrix<T> EXTnorm(nbVert,n);
  for (int iVert=0; iVert<nbVert; iVert++) {
    T eVal=EXT(iVert,0);
    T eValInv=1/eVal;
    for (int i=0; i<n; i++)
      EXTnorm(iVert,i)=eValInv*EXT(iVert,i);
  }
  std::vector<T> ListLow(dim);
  std::vector<T> ListUpp(dim);
  std::vector<int> ListSize(dim);
  for (int iDim=0; iDim<dim; iDim++) {
    MyVector<T> eCol(nbVert);
    for (int iVert=0; iVert<nbVert; iVert++)
      eCol(iVert)=EXTnorm(iVert,iDim+1);
    T eLow=Ceil(eCol.minCoeff());
    T eUpp=Floor(eCol.maxCoeff());
    ListLow[iDim]=eLow;
    ListUpp[iDim]=eUpp;
    T eVal_T=1 + eUpp - eLow;
    int eVal_I=UniversalTypeConversion<int,T>(eVal_T);
    ListSize[iDim]=eVal_I;
  }
  BlockIterationMultiple BlIter(ListSize);
  auto IsCorrect=[&](MyVector<T> const& eVect) -> bool {
    int nbFac=FAC.rows();
    for (int iFac=0; iFac<nbFac; iFac++) {
      T eScal=0;
      for (int i=0; i<n; i++)
	eScal += FAC(iFac,i) * eVect(i);
      if (eScal <0)
	return false;
    }
    return true;
  };
  std::vector<MyVector<T>> ListPoint;
  for (auto const& eVect : BlIter) {
    MyVector<T> ePoint(n);
    ePoint(0)=1;
    for (int iDim=0; iDim<dim; iDim++)
      ePoint(iDim+1)=ListLow[iDim] + eVect[iDim];
    bool test=IsCorrect(ePoint);
    if (test)
      ListPoint.push_back(ePoint);
  }
  return ListPoint;
}




#endif
