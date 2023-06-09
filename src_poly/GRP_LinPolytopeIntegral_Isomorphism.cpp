// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Isomorphism [EXT1] [EXT2] [OutEquiv]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Isomorphism [EXT1] [EXT2]\n";
      std::cerr << "\n";
      std::cerr << "OutEquiv : The equivalence information file (otherwise "
                   "printed to screen)\n";
      return -1;
    }
    //
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    using Tidx_value = uint32_t;
    using Tgr = GraphBitset;
    //
    std::string FileExt1 = argv[1];
    std::string FileExt2 = argv[2];
    MyMatrix<Tint> EXT1 = ReadMatrixFile<Tint>(FileExt1);
    MyMatrix<Tint> EXT2 = ReadMatrixFile<Tint>(FileExt2);
    size_t nbCol = EXT1.cols();
    size_t nbRow = EXT1.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    //    const bool use_scheme = true;
    const bool use_scheme = false;
    std::optional<MyMatrix<Tint>> equiv =
        LinPolytopeIntegral_Isomorphism<Tint, Tidx, Tgroup, Tidx_value, Tgr,
                                        use_scheme>(EXT1, EXT2);
    //
    auto print_info = [&](std::ostream &os) -> void {
      if (equiv) {
        os << "return ";
        WriteMatrixGAP(os, *equiv);
        os << ";\n";
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
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytopeIntegral_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
