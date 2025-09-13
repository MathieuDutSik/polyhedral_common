// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryQuadField.h"
#include "Temp_PerfectForm_Enum.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    int n;
    int i, j;
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "InitialPerfectMatrix n\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    sscanf(argv[1], "%d", &n);
    LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(n);
    MyMatrix<T> ThePerfMat = GetOnePerfectForm<T,Tint>(LinSpa, std::cerr).first;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
        std::cout << " " << ThePerfMat(i, j);
      std::cout << "\n";
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
