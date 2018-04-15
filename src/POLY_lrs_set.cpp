#include "NumberTheory.h"
#include "POLY_lrslib.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_lrs [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "DATAIN : The polyhedral cone inequalities\n";
      return -1;
    }
    //
    std::ifstream is(argv[1]);
    MyMatrix<mpq_class> EXT=ReadMatrixLrsCdd<mpq_class>(is);
    //
    std::vector<Face> ListFace=lrs::DualDescription_temp_incd(EXT);
    std::cerr << "nbVert = " << ListFace.size() << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
