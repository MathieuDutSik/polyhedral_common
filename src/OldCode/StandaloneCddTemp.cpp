#include <iostream>
#include "gmpxx.h"
#include "Temp_cdd.h"


using namespace std;


int main(int argc, char *argv[])
{
  int n;
  int i, j, k;
  int iFac, iCol, nbFac, nbCol;
  MyMatrix<mpq_class> TheEXT, TheFAC;
  mpq_class a, b, smallVal;
  mpq_class eVal;
  if (argc != 3)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneLrsTemp n k\n");
      fprintf(stderr, "and creates a Cyclic Polytope\n");
      return -1;
    }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &k);
  fprintf(stderr, "n=%d  k=%d\n", n, k);
  TMat_Allocate(&TheEXT, n, k+1);
  for (i=1; i<=n; i++)
    {
      a=1;
      b=i+1;
      for (j=0; j<=k; j++)
	{
	  TMat_Assign(&TheEXT, i-1,j, a);
	  eVal=a;
	  cout << " " << eVal;
	  a=a*b;
	}
      cout << "\n";
    }
  smallVal=0;
  cdd::DualDescription<mpq_class>(&TheEXT, &TheFAC, smallVal);
  fprintf(stderr, "We are finished\n");
  nbFac=TheFAC.nbRow;
  nbCol=TheFAC.nbCol;
  for (iFac=0; iFac<nbFac; iFac++)
    {
      cout << "iFac=" << iFac << " : ";
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  eVal=TMat_Get(&TheFAC, iFac,iCol);
	  cout << " " << eVal;
	  a=a*b;
	}
      cout << "\n";
    }
  TMat_Free(&TheFAC);
  TMat_Free(&TheEXT);
}
