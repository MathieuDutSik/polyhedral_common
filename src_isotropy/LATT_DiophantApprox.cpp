// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "SimulDiophantApprox.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_DiophantApprox [FileIn] [epsilon]\n";
      std::cerr << "\n";
      std::cerr << "FileIn  : The vector file on on input\n";
      std::cerr << "epsilon : The epsilon value on input\n";
      return -1;
    }
    //    using T = mpq_class;
    using T = double;
    using Tint = mpz_class;
    //
    std::ifstream is(argv[1]);
    MyVector<T> V = ReadVector<T>(is);
    //
    std::string epsilon_str = argv[2];
    std::stringstream s(epsilon_str);
    T epsilon;
    s >> epsilon;
    //
    DiophantResult<Tint> eRes =
      SimultaneousDiophantineApproximation<T, Tint>(V, epsilon, std::cerr);
    //
    std::cerr << "Proposed solution = " << GapStringDiophantineApprox(eRes)
              << "\n";
    T penalty = ComputeDiophantinePenalty(V, eRes);
    std::cerr << "penalty=" << penalty << "\n";
    //
    std::cerr << "Normal termination of LATT_DiophantApprox\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_DiophantApprox\n";
    exit(e.eVal);
  }
  runtime(time);
}
