#include <iostream>
#include <vector>

#include "stdlib.h"

using namespace std;

#include "gmpxx.h"
#include "Temp_lrslib.h"

int main(int argc, char *argv[])
{
  int n;
  int i, j, k;
  vector<vector<mpq_class> > TheEXT, TheFAC;
  vector<mpq_class> eLine;
  mpq_class a, b;
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
  for (i=1; i<=n; i++)
    {
      a=1;
      b=i+1;
      for (j=0; j<=k; j++)
	{
	  eLine.push_back(a);
	  fprintf(stderr, " ");
	  mpq_out_str(stderr, 10, a.get_mpq_t());
	  a=a*b;
	}
      fprintf(stderr, "\n");
      TheEXT.push_back(eLine);
      eLine.clear();
    }
  lrs::DualDescription_temp<mpq_class>(TheEXT, TheFAC);
  fprintf(stderr, "We are finished\n");
}
