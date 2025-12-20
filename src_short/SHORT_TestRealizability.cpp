// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "SHORT_ShortestConfig.h"
#include "rational.h"

template<typename T, typename Tint>
void test_realizability(std::string const& FileSHV, std::string const& OutFormat, std::string const& OutFile) {
  using Tidx = uint16_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;

  MyMatrix<Tint> SHV = ReadMatrixFile<Tint>(FileSHV);


  ReplyRealizability<T, Tint> eRes = SHORT_TestRealizabilityShortestFamily<T, Tint, Tgroup>(SHV, std::cerr);

  auto f_out=[&](std::ostream& osf) -> void {
    if (OutFormat == "TXT") {
      if (eRes.reply) {
        osf << "Realizable with following matrix\n";
        WriteMatrix(osf, eRes.eMat);
      } else {
        osf << "Non Realizable\n";
      }
      return;
    }
    if (OutFormat == "GAP") {
      osf << "return rec(realizable:=";
      if (eRes.reply) {
        osf << "true, matrix:=" << StringMatrixGAP(eRes.eMat);
      } else {
        osf << "false";
      }
      osf << ");\n";
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
    std::cerr << "Allowed OutFormat are GAP, TXT\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f_out);
}



int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_TestRealizability [arith] [FileIn]\n";
      std::cerr << "or\n";
      std::cerr << "SHORT_TestRealizability [arith] [FileIn] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "[arith]        : The arithmetic being used, only gmp is possible\n";
      std::cerr << "[FileIn]       : The input file of the system\n";
      std::cerr << "[OutFormat]    : The output format, TXT or GAP\n";
      std::cerr << "[OutFile]      : The output file of the program\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileSHV = argv[2];
    std::string OutFormat = argv[3];
    std::string OutFile = argv[4];
    //
    auto f_treat=[&]() -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return test_realizability<T,Tint>(FileSHV, OutFormat, OutFile);
      }
      std::cerr << "SHORT_TestRealizability failed to find mathching arothmetic\n";
      throw TerminalException{1};
    };
    f_treat();
    //
    //
    std::cerr << "Normal termination of SHORT_TestRealizability\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_TestRealizability\n";
    exit(e.eVal);
  }
  runtime(time1);
}
