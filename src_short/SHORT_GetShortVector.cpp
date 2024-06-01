// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "Positivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_GetShortVector [FileIn] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "[FileIn]   : The input file of the system\n";
      std::cerr
          << "[FileOut]  : The output file of the program (GAP readable)\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::ifstream is(argv[1]);
    MyMatrix<T> M = ReadMatrix<T>(is);
    std::cerr << "We have |M|=" << M.rows() << " / " << M.cols() << "\n";
    //
    int StrictIneq_i;
    is >> StrictIneq_i;
    bool StrictIneq;
    if (StrictIneq_i == 1) {
      StrictIneq = true;
    } else {
      StrictIneq = false;
    }
    std::cerr << "StrictIneq=" << StrictIneq << "\n";
    //
    int NeedNonZero_i;
    is >> NeedNonZero_i;
    bool NeedNonZero;
    if (NeedNonZero_i == 1) {
      NeedNonZero = true;
    } else {
      NeedNonZero = false;
    }
    std::cerr << "NeedNonZero=" << NeedNonZero << "\n";
    //
    T CritNorm;
    //
    is >> CritNorm;
    std::cerr << "CritNorm=" << CritNorm << "\n";
    //
    //    std::cerr << "StrictIneq=" << StrictIneq << " NeedNonZero=" <<
    //    NeedNonZero << " CritNorm=" << CritNorm << "\n";
    //
    MyVector<Tint> eVect = GetShortIntegralVector<T, Tint>(M, CritNorm, StrictIneq, NeedNonZero, std::cerr);
    //
    std::string FileOut = argv[2];
    std::ofstream os(FileOut);
    os << "return [";
    int n = eVect.size();
    for (int i = 0; i < n; i++) {
      if (i > 0)
        os << ",";
      os << eVect(i);
    }
    os << "];\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_GetShortVector\n";
    exit(e.eVal);
  }
  runtime(time1);
}
