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
    //
    std::string FileIn = argv[1];
    std::vector<MyMatrix<Tint>> ListSHV =
        ReadListConfigurationShortestVector<Tint>(FileIn);
    std::unordered_set<MyMatrix<Tint>> set;
    for (auto & SHV : ListSHV) {
      MyMatrix<Tint> SHV_can = SHORT_Canonicalize<T,Tint>(SHV, std::cerr);
      set.insert(SHV_can);
    }
    std::vector<MyMatrix<Tint>> v;
    for (auto & fSHV : set) {
      v.push_back(fSHV);
    }
    //
    std::string FileOut = argv[2];
    std::ofstream os(FileOut);
    os << "return rec(ListReduced:=[";
    bool IsFirst = true;
    for (auto &eSHV : v) {
      if (!IsFirst)
        os << ",\n";
      IsFirst = false;
      WriteMatrixGAP(os, eSHV);
    }
    os << "]);\n";
    std::cerr << "Normal termination of the SHORT_ReduceVectorFamilyGAP\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_ReduceVectorFamilyGAP\n";
    exit(e.eVal);
  }
  runtime(time1);
}
