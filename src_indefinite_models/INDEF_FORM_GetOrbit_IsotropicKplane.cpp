// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "CombinedAlgorithms.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &MatFile, std::string const& KStr, std::string choice, std::string const &OutFormat,
             std::ostream &os_out) {
  MyMatrix<T> Qmat = ReadMatrixFile<T>(MatFile);
  int k = ParseScalar<int>(KStr);
  std::cerr << "We have Q\n";
  IndefiniteCombinedAlgo<T,Tint,Tgroup> comb(std::cerr);
  auto f_get=[&]() -> std::vector<MyMatrix<Tint>> {
    if (choice == "plane") {
      return comb.INDEF_FORM_GetOrbit_IsotropicKplane(Qmat, k);
    }
    if (choice == "flag") {
      return comb.INDEF_FORM_GetOrbit_IsotropicKflag(Qmat, k);
    }
    std::cerr << "No correct choice\n";
    throw TerminalException{1};
  };
  std::vector<MyMatrix<Tint>> l_planes = f_get();
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteListMatrixGAP(os_out, l_gen);
    os_out << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "INDEF_FORM_GetOrbit_IsotropicKplane_flag [arith] [MatFile] [k] [choice]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_GetOrbit_IsotropicKplane_flag [arith] [MatFile] [k] [choice] [OutFormat] [OutFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string MatFile = argv[2];
    std::string KStr = argv[3];
    std::string choice = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
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
        return process<T,Tint,Tgroup>(MatFile, KStr, choice, OutFormat, os);
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
    std::cerr << "Error in INDEF_FORM_GetOrbit_IsotropicKplane_flag\n";
    exit(e.eVal);
  }
  runtime(time1);
}
