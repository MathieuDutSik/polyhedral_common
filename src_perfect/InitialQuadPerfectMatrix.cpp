// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryQuadField.h"
#include "Temp_PerfectForm_Enum.h"
// clang-format on

typedef QuadField<mpq_class, 2> field;

int main(int argc, char *argv[]) {
  try {
    int n;
    int i, j;
    field eVal;
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "InitialPerfectMatrix n\n";
      return -1;
    }
    using field = QuadField<mpq_class, 2>;
    sscanf(argv[1], "%d", &n);
    LinSpaceMatrix<field> LinSpa = ComputeCanonicalSpace<field>(n);
    MyMatrix<field> ThePerfMat = GetOnePerfectForm(LinSpa);
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
        std::cout << " " << ThePerfMat(i, j);
      std::cout << "\n";
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
