// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "LatticeStabEquiCan.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "LATT_Isomorphism [ListMat1] [ListMat2] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Isomorphism [ListMat1] [ListMat2]\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
#else
    using T = mpq_class;
    using Tint = mpz_class;
#endif
    //
    std::string FileListMat1 = argv[1];
    std::string FileListMat2 = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    std::vector<MyMatrix<T>> ListMat1 = ReadListMatrixFile<T>(FileListMat1);
    std::vector<MyMatrix<T>> ListMat2 = ReadListMatrixFile<T>(FileListMat2);

    std::optional<MyMatrix<Tint>> equiv =
        ArithmeticEquivalenceMultiple<T, Tint>(ListMat1, ListMat2, std::cerr);
    //
    auto prt = [&](std::ostream &os) -> void {
      if (OutFormat == "GAP") {
        if (equiv) {
          os << "return ";
          WriteMatrixGAP(os, *equiv);
          os << ";\n";
        } else {
          os << "return false;\n";
        }
        return;
      }
      if (OutFormat == "Oscar") {
        if (equiv) {
          WriteMatrix(os, *equiv);
        } else {
          os << "0 0\n";
        }
        return;
      }
      std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat
                << "\n";
      throw TerminalException{1};
    };
    if (OutFile == "stderr") {
      prt(std::cerr);
    } else {
      if (OutFile == "stdout") {
        prt(std::cout);
      } else {
        std::ofstream os(OutFile);
        prt(os);
      }
    }
    std::cerr << "Normal termination of LATT_Isomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Raised exception led to premature end of LATT_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
