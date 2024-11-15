// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "CombinedAlgorithms.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &File1, std::string const &File2,
             std::string const &OutFormat, std::ostream &os_out) {
  MyMatrix<T> Qmat1 = ReadMatrixFile<T>(File1);
  MyMatrix<T> Qmat2 = ReadMatrixFile<T>(File2);
  std::cerr << "We have Q\n";
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  std::optional<MyMatrix<Tint>> opt =
      comb.INDEF_FORM_TestEquivalence(Qmat1, Qmat2);
  if (OutFormat == "PYTHON") {
    if (opt) {
      WriteMatrixPYTHON(os_out, *opt);
    } else {
      os_out << "None";
    }
    os_out << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    if (opt) {
      WriteMatrixGAP(os_out, *opt);
    } else {
      os_out << "fail";
    }
    os_out << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "INDEF_FORM_TestEquivalence [arith] [File1] [File2]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_TestEquivalence [arith] [File1] [File2] "
                   "[OutFormat] [OutFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string File1 = argv[2];
    std::string File2 = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint, Tgroup>(File1, File2, OutFormat, os);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    if (OutFile == "stderr") {
      f(std::cerr);
    } else {
      if (OutFile == "stdout") {
        f(std::cout);
      } else {
        std::ofstream os(OutFile);
        f(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_TestEquivalence\n";
    exit(e.eVal);
  }
  runtime(time);
}
