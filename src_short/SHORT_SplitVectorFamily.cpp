// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "SHORT_ShortestConfig.h"

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "SHORT_SplitVectorFamily [FileIn] [N] [mod] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "[FileIn]   : The input file of the system\n";
      std::cerr << "[N]        : The total number of blocks\n";
      std::cerr << "[mod]      : The residue in output\n";
      std::cerr << "[FileOut]  : The output file of the program\n";
      return -1;
    }
    using Tint = int;
    //
    std::string FileIn = argv[1];
    std::vector<MyMatrix<Tint>> ListSHV =
        ReadListConfigurationShortestVector<Tint>(FileIn);
    int nbShort = ListSHV.size();
    //
    int N;
    sscanf(argv[2], "%d", &N);
    //
    int mod;
    sscanf(argv[3], "%d", &mod);
    //
    int nbOut = 0;
    for (int iShort = 0; iShort < nbShort; iShort++) {
      int res = iShort % N;
      if (res == mod)
        nbOut++;
    }
    //
    std::string FileOut = argv[4];
    std::ofstream os(FileOut);
    os << nbOut << "\n";
    for (int iShort = 0; iShort < nbShort; iShort++) {
      int res = iShort % N;
      if (res == mod)
        WriteMatrix(os, ListSHV[iShort]);
    }
    std::cerr << "Normal termination of SHORT_SplitVectorFamily\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in SHORT_SplitVectorFamily\n";
    exit(e.eVal);
  }
  runtime(time1);
}
