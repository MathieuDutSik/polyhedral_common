#include "Temp_PolytopeFct.h"
#include "Temp_PolytopeEquiStab.h"

#include "gmpxx.h"
#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
  std::ifstream fs;
  MyMatrix<mpq_class> *TheEXT;
  vector<int> eInc;
  
  TheEXT = new MyMatrix<mpq_class>;
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneInitialVertex [DATAINP]\n");
      fprintf(stderr, "DATAINP: The input data of the form\n");
      fprintf(stderr, "n\n");
      fprintf(stderr, "A00  A01 ..... A0p\n");
      fprintf(stderr, ".     .          .\n");
      fprintf(stderr, ".     .          .\n");
      fprintf(stderr, "An0  A01 ..... Anp\n");
      return -1;
    }
  fprintf(stderr, "  input file:%s\n", argv[1]);
  fs.open(argv[1]);
  ReadPolytope(fs, TheEXT);
  fs.close();
  FindOneInitialVertex(TheEXT, eInc);

  TMat_Free(TheEXT);
  delete TheEXT;
  fprintf(stderr, "Completion of the program\n");
  
}
