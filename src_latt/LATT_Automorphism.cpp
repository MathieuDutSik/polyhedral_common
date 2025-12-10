// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Automorphism [ListMat] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Automorphism [ListMat]\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
#else
    using T = mpq_class;
    using Tint = mpz_class;
#endif
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    std::string FileListMat = argv[1];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      OutFile = argv[3];
    }
    std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(FileListMat);

    std::vector<MyMatrix<Tint>> ListGen =
        ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat,
                                                             std::cerr);
    //
    auto prt = [&](std::ostream &os) -> void {
      if (OutFormat == "GAP") {
        os << "return [";
        bool IsFirst = true;
        for (auto &eMat : ListGen) {
          if (!IsFirst)
            os << ",\n";
          IsFirst = false;
          WriteMatrixGAP(os, eMat);
        }
        os << "];\n";
        return;
      }
      if (OutFormat == "Oscar") {
        os << ListGen.size() << "\n";
        for (auto &eMat : ListGen) {
          WriteMatrix(os, eMat);
        }
        return;
      }
      std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat
                << "\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, prt);
    std::cerr << "Normal termination of LATT_Automorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
