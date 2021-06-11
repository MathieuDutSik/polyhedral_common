#include "NumberTheory.h"
#include "Permutation.h"
#include "Group.h"
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinPolytope_TestAutomorphism [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN : The list of vertices (or inequalities for that matter)\n";
      return -1;
    }
    //
    using T=mpq_class;
    using Tint = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt,Tint>;
    //
    std::ifstream is(argv[1]);
    MyMatrix<T> EXT=ReadMatrix<T>(is);
    //
    const bool use_scheme1 = true;
    Tgroup GRP1 = LinPolytope_Automorphism<T,use_scheme1,Tgroup>(EXT);
    //
    const bool use_scheme2 = false;
    Tgroup GRP2 = LinPolytope_Automorphism<T,use_scheme2,Tgroup>(EXT);
    //
    bool test = GRP1 == GRP2;
    if (!test) {
      std::cerr << "The groups are different. It is the clear bug\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
