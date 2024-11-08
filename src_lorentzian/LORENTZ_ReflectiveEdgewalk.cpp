// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "Group.h"
#include "Permutation.h"
#include "edgewalk.h"
// clang-format on


template<typename T, typename Tint>
void process(std::string const& MatFile, std::string const& OutFormat, std::ostream& os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> LorMat = ReadMatrixFile<T>(MatFile);
  //
  ResultEdgewalk<T,Tint> re = StandardEdgewalkAnalysis<T,Tint,Tgroup>(LorMat, std::cerr);
  bool ComputeAllSimpleRoots = true;
  PrintResultEdgewalk(LorMat, re, os_out, OutFormat, ComputeAllSimpleRoots, std::cerr);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "LORENTZ_ReflectiveEdgewalk [arith] [MatFile]\n";
      std::cerr << "     or\n";
      std::cerr << "LORENTZ_ReflectiveEdgewalk [arith] [MatFile] [OutFormat] [OutFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string MatFile = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f=[&](std::ostream& os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T,Tint>(MatFile, OutFormat, os_out);
      }
      std::cerr << "Failed to find matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    //
    if (OutFile == "stderr") {
      f(std::cerr);
    } else {
      if (OutFile == "stdout") {
        f(std::cout);
      } else {
        std::ofstream os_out(OutFile);
        f(os_out);
      }
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_FundDomain_AllcockEdgewalk\n";
    exit(e.eVal);
  }
  runtime(time);
}
