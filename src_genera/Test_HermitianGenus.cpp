// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "HermitianGenus.h"
#include <string>
#include <vector>
// clang-format on

/*
  The genera of positive definite hermitian lattices over Q(zeta_d) /
  Q(zeta_d)^+ (HermitianGenus.h), for comparison with hermitian_genera of
  Hecke.jl:
    count d b N : a line "d b N number", the genera of rank b, scales
                  between O_E and N O_E
 */

int main(int argc, char *argv[]) {
  try {
    if (argc != 5 || std::string(argv[1]) != "count") {
      std::cerr << "Test_HermitianGenus count d b N\n";
      return -1;
    }
    int64_t d = std::stoll(argv[2]);
    int b = std::stoi(argv[3]);
    int64_t N = std::stoll(argv[4]);
    std::vector<HermPrimeSymbols> ps = HermitianPrimeSymbols(d, b, N);
    std::cout << d << " " << b << " " << N << " "
              << HermitianGlobalGenera(ps).size() << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_HermitianGenus\n";
    exit(e.eVal);
  }
}
