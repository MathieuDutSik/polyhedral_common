// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 5 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Isomorphism [EXT1] [EXT2] [OutFormat] "
                   "[OutEquiv]\n";
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
    std::string FileExt1 = argv[1];
    std::string FileExt2 = argv[2];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    MyMatrix<Tint> EXT1 = ReadMatrixFile<Tint>(FileExt1);
    MyMatrix<Tint> EXT2 = ReadMatrixFile<Tint>(FileExt2);
    size_t nbCol = EXT1.cols();
    size_t nbRow = EXT1.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    const bool use_scheme = true;
    std::optional<std::vector<Tidx>> equiv =
        LinPolytope_Isomorphism<Tint, Tidx, use_scheme>(EXT1, EXT2, std::cerr);
    //
    auto print_info = [&](std::ostream &os) -> void {
      if (OutFormat == "Oscar") {
        if (equiv) {
          const std::vector<Tidx> &V = *equiv;
          os << nbRow << "\n";
          for (size_t iRow = 0; iRow < nbRow; iRow++) {
            if (iRow > 0)
              os << " ";
            int eVal = V[iRow] + 1;
            os << eVal;
          }
        } else {
          os << "-1\n";
        }
        return;
      }
      if (OutFormat == "GAP") {
        if (equiv) {
          const std::vector<Tidx> &V = *equiv;
          os << "return [";
          for (size_t iRow = 0; iRow < nbRow; iRow++) {
            if (iRow > 0)
              os << ",";
            os << (V[iRow] + 1);
          }
          os << "];\n";
        } else {
          os << "return fail;\n";
        }
        return;
      }
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      print_info(std::cerr);
    } else {
      if (FileO == "stdout") {
        print_info(std::cout);
      } else {
        std::ofstream os(FileO);
        print_info(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
