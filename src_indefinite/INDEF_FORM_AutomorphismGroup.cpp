// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "CombinedAlgorithms.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &MatFile, std::string const &OutFormat,
             std::ostream &os_out) {
  MyMatrix<T> Qmat = ReadMatrixFile<T>(MatFile);
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  std::vector<MyMatrix<Tint>> l_gen = comb.INDEF_FORM_AutomorphismGroup(Qmat);
  for (auto & e_gen : l_gen) {
    MyMatrix<T> e_gen_T = UniversalMatrixConversion<T,Tint>(e_gen);
    MyMatrix<T> prod = e_gen_T * Qmat * e_gen_T.transpose();
    if (prod != Qmat) {
      std::cerr << "The matrix is not preserving the space\n";
      throw TerminalException{1};
    }
  }
  if (OutFormat == "CPP") {
    return WriteListMatrix(os_out, l_gen);
  }
  if (OutFormat == "PYTHON") {
    return WriteListMatrixPYTHON(os_out, l_gen);
  }
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
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "INDEF_FORM_AutomorphismGroup [arith] [MatFile]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_AutomorphismGroup [arith] [MatFile] [OutFormat] "
                   "[OutFile]\n";
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
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint, Tgroup>(MatFile, OutFormat, os);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of INDEF_FORM_AutomorphismGroup\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_AutomorphismGroup, runtime=" << time << "\n";
    exit(e.eVal);
  }
  runtime(time);
}
