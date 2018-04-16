#ifndef INCLUDE_COMBINATORICS_ELEM
#define INCLUDE_COMBINATORICS_ELEM

#include "Temp_common.h"
#include "MAT_Matrix.h"

std::vector<int> BinomialStdvect_First(int const& k)
{
  std::vector<int> eVect(k);
  for (int i=0; i<k; i++)
    eVect[i]=i;
  return eVect;
}


bool BinomialStdvect_Increment(int const&n, int const&k, std::vector<int> & Tvect)
{
  Tvect[0]++;
  int xy2=1;
  while ((xy2 < k) && (Tvect[xy2-1] >= Tvect[xy2])) {
    Tvect[xy2]++;
    xy2++;
  }
  if (xy2 != 1) {
    for (int xy1=0; xy1<xy2-1; xy1++)
      Tvect[xy1]=xy1;
  }
  if (Tvect[k-1] == n)
    return false;
  return true;
}



MyMatrix<int> BuildSet(int const& n, int const& Nval)
{
  int TotalNb=MyPow<int>(Nval, n);
  std::vector<int> eVect(n, 0);
  std::function<int(void)> fUpdate=[&]() -> int {
    for (int iVar=0; iVar<n; iVar++)
      if (eVect[iVar] < Nval-1) {
	eVect[iVar]++;
	for (int jVar=0; jVar<iVar; jVar++)
	  eVect[jVar]=0;
	return 0;
      }
    return -1;
  };
  MyMatrix<int> RetMatrix(TotalNb, n);
  int idx=0;
  while(true) {
    for (int i=0; i<n; i++)
      RetMatrix(idx,i)=eVect[i];
    idx++;
    int test=fUpdate();
    if (test == -1)
      break;
  }
  return RetMatrix;
}


int PositionBuildSet(int const& n, int const& Nval, MyVector<int> const& V)
{
  int pos=0;
  int eProd=1;
  for (int i=0; i<n; i++) {
    pos += eProd * V(i);
    eProd *= Nval;
  }
  return pos;
}






#endif


