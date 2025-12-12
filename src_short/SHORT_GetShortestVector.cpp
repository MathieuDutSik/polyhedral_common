// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "rational.h"
#include "SHORT_ShortestConfig.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_GetShortestVector [FileIn] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "[FileIn]   : The input file of the system\n";
      std::cerr
          << "[FileOut]  : The output file of the program (GAP readable)\n";
      return -1;
    }
    //
    using T = mpq_class;
    //    using Tint=mpz_class;
    using Tint = int;
    //    using Tgroup=TheGroupFormat<mpz_class>;
    std::ifstream is(argv[1]);
    MyMatrix<T> M = ReadMatrix<T>(is);
    //
    Tshortest<T, Tint> rec_shv = T_ShortestVector<T, Tint>(M, std::cerr);
    int nbRow = rec_shv.SHV.rows();
    std::cerr << "nbRow=" << nbRow << "\n";
    std::cerr << "Normal termination of the SHORT_GetShortestVector\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_GetShortestVector\n";
    exit(e.eVal);
  }
  runtime(time1);
}
