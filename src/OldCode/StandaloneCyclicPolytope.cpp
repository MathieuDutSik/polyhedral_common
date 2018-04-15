#include <iostream>
#include <vector>
#include "stdlib.h"
#include "MPQ_CDD_LRS.h"

using namespace std;

int main(int argc, char *argv[])
{
  int n;
  int i, j, k;
  vector<vector<mpq_class> > TheEXT, TheFAC;
  vector<mpq_class> eLine;
  mpq_class a, b;
  MyMatrix TheEXT_old, TheFAC_old;
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
  MATRIX_Allocate(&TheEXT_old, n, k+1);
  for (i=1; i<=n; i++)
    {
      a=1;
      b=i+1;
      for (j=0; j<=k; j++)
	{
	  fprintf(stderr, " ");
	  mpq_out_str(stderr, 10, a.get_mpq_t());
	  MATRIX_Assign(&TheEXT_old, i, j+1, a.get_mpq_t());
	  a=a*b;
	}
      fprintf(stderr, "\n");
    }
  Initialize_All();
  LRS_DualDescription(&TheEXT_old, &TheFAC_old);
  Free_All();
}
