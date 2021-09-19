#include "NumberTheory.h"
#include "Permlib_specific.h"
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 4 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Isomorphism [EXT1] [EXT2] [OutEquiv]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Isomorphism [EXT1] [EXT2]\n";
      std::cerr << "\n";
      std::cerr << "OutEquiv : The equivalence information file (otherwise printed to screen)\n";
      return -1;
    }
    //
    using Tint = mpz_class;
    //    using Tidx = uint16_t;
    using Tidx = uint32_t;
    std::ifstream is1(argv[1]);
    MyMatrix<Tint> EXT1=ReadMatrix<Tint>(is1);
    std::ifstream is2(argv[2]);
    MyMatrix<Tint> EXT2=ReadMatrix<Tint>(is2);
    size_t nbCol=EXT1.cols();
    size_t nbRow=EXT1.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    const bool use_scheme = true;
    std::optional<std::vector<Tidx>> equiv = LinPolytope_Isomorphism<Tint,Tidx,use_scheme>(EXT1, EXT2);
    //
    auto print_info=[&](std::ostream& os) -> void {
      if (equiv) {
        const std::vector<Tidx>& V = *equiv;
        os << "return [";
        for (size_t iRow=0; iRow<nbRow; iRow++) {
          if (iRow > 0)
            os << ",";
          os << (V[iRow] + 1);
        }
        os << "];\n";
      } else {
        os << "return fail;\n";
      }
    };
    if (argc == 4) {
      std::ofstream os(argv[3]);
      print_info(os);
    }
    if (argc == 3) {
      print_info(std::cerr);
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
