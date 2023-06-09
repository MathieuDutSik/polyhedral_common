// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "Temp_ShortVectorUndefinite.h"

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
    //
    std::ifstream is(argv[1]);
    MyMatrix<mpq_class> M = ReadMatrix<mpq_class>(is);
    //
    int StrictIneq_i;
    is >> StrictIneq_i;
    bool StrictIneq;
    if (StrictIneq_i == 1) {
      StrictIneq = true;
    } else {
      StrictIneq = false;
    }
    //
    int NeedNonZero_i;
    is >> NeedNonZero_i;
    bool NeedNonZero;
    if (NeedNonZero_i == 1) {
      NeedNonZero = true;
    } else {
      NeedNonZero = false;
    }
    //
    mpq_class CritNorm;
    //
    is >> CritNorm;
    //
    //    std::cerr << "StrictIneq=" << StrictIneq << " NeedNonZero=" <<
    //    NeedNonZero << " CritNorm=" << CritNorm << "\n";
    //
    MyVector<mpz_class> eVect =
        GetShortVector_unlimited_float<mpz_class, mpq_class>(
            M, CritNorm, StrictIneq, NeedNonZero);
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
