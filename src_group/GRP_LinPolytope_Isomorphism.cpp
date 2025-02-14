// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryQuadField.h"
#include "GRP_GroupFct.h"
#include "PolytopeEquiStab.h"
// clang-format on

template <typename T>
void process(std::string const &FileExt1, std::string const &FileExt2,
             std::string const &OutFormat, std::string const &FileO) {
  using Tidx = uint32_t;
  MyMatrix<T> EXT1 = ReadMatrixFile<T>(FileExt1);
  MyMatrix<T> EXT2 = ReadMatrixFile<T>(FileExt2);
  size_t nbRow = EXT1.rows();
  //
  std::optional<std::vector<Tidx>> equiv =
      LinPolytope_Isomorphism<T, Tidx>(EXT1, EXT2, std::cerr);
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
          int eVal = V[iRow] + 1;
          os << eVal;
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
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 6 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "GRP_LinPolytope_Isomorphism [arith] [EXT1] [EXT2] [OutFormat] "
             "[OutEquiv]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Isomorphism [arith] [EXT1] [EXT2]\n";
      std::cerr << "\n";
      std::cerr << "OutEquiv : The equivalence information file (otherwise "
                   "printed to screen)\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileExt1 = argv[2];
    std::string FileExt2 = argv[3];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[4];
      FileO = argv[5];
    }
    auto f = [&]() -> void {
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(FileExt1, FileExt2, OutFormat, FileO);
      }
      if (arith == "Qsqrt3") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 3>;
        return process<T>(FileExt1, FileExt2, OutFormat, FileO);
      }
      std::cerr << "Failed to find a matching arithmetic\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of GRP_LinPolytope_Isomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
