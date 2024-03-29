// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "rational.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "SHORT_ShortestConfig.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_ReduceVectorFamilyGAP [FileIn] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "[FileIn]   : The input file of the system\n";
      std::cerr
          << "[FileOut]  : The output file of the program (GAP readable)\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = int;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, mpz_class>;
    //
    std::string FileIn = argv[1];
    std::vector<MyMatrix<Tint>> ListSHV =
        ReadListConfigurationShortestVector<Tint>(FileIn);
    std::pair<std::vector<MyMatrix<Tint>>, std::vector<Tint>> RecRet =
        SHORT_ReduceByIsomorphism<T, Tint, Tgroup>(ListSHV, std::cerr);
    //
    std::string FileOut = argv[2];
    std::ofstream os(FileOut);
    os << "return rec(ListReduced:=[";
    bool IsFirst = true;
    for (auto &eSHV : RecRet.first) {
      if (!IsFirst)
        os << ",\n";
      IsFirst = false;
      WriteMatrixGAP(os, eSHV);
    }
    os << "],\nVectPos:=[";
    int nbConfTot = RecRet.second.size();
    for (int iConf = 0; iConf < nbConfTot; iConf++) {
      if (iConf > 0)
        os << ",";
      os << RecRet.second[iConf];
    }
    os << "]);\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_ReduceVectorFamilyGAP\n";
    exit(e.eVal);
  }
  runtime(time1);
}
