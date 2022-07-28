// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

#include "POLY_Heuristics.h"

int main(int argc, char *argv[]) {
  try {
    //  using T = mpq_class;
    using T = T_uint64_t;
    // The chosenoptions
    FullNamelist eFull = StandardHeuristicDualDescriptionProgram_TS();
    if (argc == 2) {
      std::cerr << "Reading one file\n";
      std::string filename = argv[1];
      NAMELIST_ReadNamelistFile(filename, eFull);
    }
    std::cerr << "Working with:\n";
    NAMELIST_WriteNamelistFile(std::cerr, eFull);


    ThompsonSamplingHeuristic<T> TSH(std::cerr, eFull);

    size_t N = 10000;
    std::map<std::string,T> TheCand;
    for (size_t i=0; i<N; i++) {
      int incidence = 30 + (random() % 100);
      TheCand["incidence"] = incidence;
      std::string choice = TSH.get_eval(TheCand);
      TSH.pop();
      std::cerr << "i=" << i << " incidence=" << incidence << " choice=" << choice << "\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
}
