// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "rational.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "SHORT_ShortestConfig.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_TestRealizability [FileIn] [TheMethod] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "[FileIn]    : The input file of the system\n";
      std::cerr << "[TheMethod] : Can be cdd or glpk_secure\n";
      std::cerr
          << "[FileOut]   : The output file of the program (GAP readable)\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = int;

    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, mpz_class>;
    //
    std::string FileMat = argv[1];
    std::string TheMethod = argv[2];
    std::string FileOut = argv[3];
    //
    std::ifstream is(FileMat);
    MyMatrix<Tint> SHV = ReadMatrix<Tint>(is);
    //
    ReplyRealizability<T, Tint> eRes =
        SHORT_TestRealizabilityShortestFamily<T, Tint, Tgroup>(SHV, TheMethod);
    //
    std::ofstream os(FileOut);
    if (eRes.reply) {
      os << "Realizable with following matrix\n";
      WriteMatrix(os, eRes.eMat);
    } else {
      os << "Non Realizable\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_TestRealizability\n";
    exit(e.eVal);
  }
  runtime(time1);
}
