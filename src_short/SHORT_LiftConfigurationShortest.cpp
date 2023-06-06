// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "rational.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "SHORT_ShortestConfig.h"
#include "SHORT_ShortestConfig_Parall.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LiftConfigurationShortest [NPROC] [DATAIN] [TheMethod] "
                   "[DATAOUT]\n";
      std::cerr << "\n";
      std::cerr << "NPROC     : The number of processor used\n";
      std::cerr << "DATAIN    : The list of configuration of shortest vector "
                   "as input\n";
      std::cerr << "TheMethod : can be cdd or glpk_secure\n";
      std::cerr << "DATAOUT   : The list of configuration of shortest vector "
                   "as output\n";
      return -1;
    }
    //
    using T = mpq_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = int;
    using Tgroup = permutalib::Group<Telt, mpz_class>;
    int NPROC;
    sscanf(argv[1], "%d", &NPROC);
    //
    std::string eFileIN(argv[2]);
    std::string TheMethod(argv[3]);
    std::vector<MyMatrix<Tint>> ListConfIn =
        ReadListConfigurationShortestVector<Tint>(eFileIN);
    //
    std::vector<MyMatrix<Tint>> ListConfOut =
        SHORT_SimplicialEnumeration<T, Tint, Tgroup>(ListConfIn, NPROC,
                                                     TheMethod);
    //
    std::string eFileOUT(argv[4]);
    WriteListConfigurationShortestVector(eFileOUT, ListConfOut);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_LiftConfigurationShortest\n";
    exit(e.eVal);
  }
  runtime(time1);
}
