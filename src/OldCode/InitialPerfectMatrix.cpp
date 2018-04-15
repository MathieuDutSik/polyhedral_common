#include <iostream>
#include <vector>
#include "gmpxx.h"
#include "Temp_PerfectForm.h"

using namespace std;


int main(int argc, char *argv[])
{
  int n;
  int i,j;
  LinSpaceMatrix<mpq_class> LinSpa;
  MyMatrix<mpq_class> *ThePerfMat;
  mpq_class eVal;
  ThePerfMat   =new MyMatrix<mpq_class>;
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "InitialPerfectMatrix n\n");
      return -1;
    }
  sscanf(argv[1], "%d", &n);
  ComputeCanonicalSpace(&LinSpa, n);
  GetOnePerfectForm(&LinSpa, ThePerfMat);
  MatrixSpace_Free(&LinSpa);
  for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)
	{
	  eVal=TMat_Get(ThePerfMat, i, j);
	  cout << " " << eVal;
	}
      cout << "\n";
    }
  TMat_Free(ThePerfMat);
  delete ThePerfMat;
}
