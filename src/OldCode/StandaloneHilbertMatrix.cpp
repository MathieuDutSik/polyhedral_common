#include "mpreal.h"
#include "stdlib.h"

using namespace mpfr;

#include "Temp_Matrix.h"

int main(int argc, char *argv[])
{
  int n, theprec;
  int i,j;
  if (argc != 3) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "StandaloneHilbertMatrix n theprec\n";
    return -1;
  }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &theprec);
  mpreal::set_default_prec(theprec);
  MyMatrix<mpreal> eMat(n, n);
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      mpreal a1=1;
      mpreal aIJ=i+j+1;
      mpreal b=a1/aIJ;
      eMat(i,j)=b;
    }
  std::cout << "Printing eMat\n";
  WriteMatrix(std::cout, eMat);
  MyMatrix<mpreal> eMatInv=Inverse(eMat);
  std::cout << "Printing eMatInv\n";
  WriteMatrix(std::cout, eMatInv);
}
