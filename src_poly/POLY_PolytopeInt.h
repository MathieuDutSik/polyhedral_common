#ifndef TEMP_POLYTOPE_INT_INCLUDE
#define TEMP_POLYTOPE_INT_INCLUDE

#include "COMB_Combinatorics.h"


template<typename T, typename Tint, typename Finsert>
void Kernel_GetListIntegralPoint(MyMatrix<T> const& FAC, MyMatrix<T> const& EXT, Finsert f_insert)
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
  std::vector<Tint> ListLow(dim);
  std::vector<Tint> ListUpp(dim);
  std::vector<int> ListSize(dim);
  for (int iDim=0; iDim<dim; iDim++) {
    MyVector<T> eCol(nbVert);
    for (int iVert=0; iVert<nbVert; iVert++)
      eCol(iVert) = EXTnorm(iVert,iDim+1);
    Tint eLow=Ceil(eCol.minCoeff());
    Tint eUpp=Floor(eCol.maxCoeff());
    ListLow[iDim] = eLow;
    ListUpp[iDim] = eUpp;
    Tint eSiz=1 + eUpp - eLow;
    ListSize[iDim]=eSiz;
  }
  int nbFac=FAC.rows();
  auto IsCorrect=[&](MyVector<T> const& eVect) -> bool {
    for (int iFac=0; iFac<nbFac; iFac++) {
      T eScal=0;
      for (int i=0; i<n; i++)
	eScal += FAC(iFac,i) * eVect(i);
      if (eScal <0)
	return false;
    }
    return true;
  };
  MyVector<Tint> ePoint(n);
  ePoint(0)=1;
  BlockIterationMultiple BlIter(ListSize);
  for (auto const& eVect : BlIter) {
    for (int iDim=0; iDim<dim; iDim++)
      ePoint(iDim+1) = ListLow[iDim] + eVect[iDim];
    bool test=IsCorrect(ePoint);
    if (test)
      f_insert(ePoint);
  }
}


template<typename T, typename Tint, typename Tint>
std::vector<MyVector<Tint>> GetListIntegralPoint(MyMatrix<T> const& FAC, MyMatrix<T> const& EXT)
{
  std::vector<MyVector<Tint>> ListPoint;
  auto f_insert=[&](const MyVector<Tint>& ePoint) -> void {
    ListPoint.push_back(ePoint);
  };
  Kernel_GetListIntegralPoint(FAC, EXT, f_insert);
  return ListPoint;
}




template<typename T, typename Tint, typename Finsert>
void Kernel_GetListIntegralPoint_LP(MyMatrix<T> const& FAC, Finsert f_insert)
{
  // Getting the bounds with the coordinates i<pos1 being fixed with values
  // set by eVert. The values of the coordinate pos2 is then computed.
  auto get_bound=[&](const MyVector<T>& eVert, const size_t& pos1, const size_t& pos2) -> std::pair<Tint,Tint> {
  };
  while(true) {
  }
}









#endif
