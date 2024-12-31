// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "GRP_GroupFct.h"
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 7 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Isomorphism_GramMat [EXT1] [GramMat1] "
                   "[EXT2] [GramMat2] [OutFormat] [OutEquiv]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Isomorphism_GramMat [EXT1] [GramMat1] "
                   "[EXT2] [GramMat2]\n";
      std::cerr << "\n";
      std::cerr << "OutEquiv : The equivalence information file (otherwise "
                   "printed to screen)\n";
      return -1;
    }
    //
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Tint = boost::multiprecision::mpz_int;
#else
    using Tint = mpz_class;
#endif
    using Tidx = uint32_t;
    std::string FileExt1 = argv[1];
    std::string FileGram1 = argv[2];
    std::string FileExt2 = argv[3];
    std::string FileGram2 = argv[4];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      FileO = argv[6];
    }
    MyMatrix<Tint> EXT1 = ReadMatrixFile<Tint>(FileExt1);
    MyMatrix<Tint> GramMat1 = ReadMatrixFile<Tint>(FileGram1);
    MyMatrix<Tint> EXT2 = ReadMatrixFile<Tint>(FileExt2);
    MyMatrix<Tint> GramMat2 = ReadMatrixFile<Tint>(FileGram2);
    size_t nbRow = EXT1.rows();
    //
    std::optional<std::vector<Tidx>> equiv =
        LinPolytope_Isomorphism_GramMat<Tint, Tidx>(EXT1, GramMat1, EXT2,
                                                    GramMat2, std::cerr);
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
          os << "0\n";
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
    std::cerr << "Error in GRP_LinPolytope_Isomorphism_GramMat\n";
    exit(e.eVal);
  }
  runtime(time1);
}
