// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

#include "POLY_Heuristics.h"

int main() {
  //  using T = mpq_class;
  using T = T_uint64_t;

  FullNamelist eFull = StandardHeuristicDualDescriptionProgram_TS<T>();

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
}
