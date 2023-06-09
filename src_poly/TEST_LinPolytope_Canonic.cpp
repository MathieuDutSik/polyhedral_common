// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_LinPolytope_Canonic [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that "
                   "matter)\n";
      return -1;
    }
    //
    using Tint = mpz_class;
    const bool use_scheme = true;
    std::string eFile = argv[1];
    MyMatrix<Tint> EXT = ReadMatrixFile<Tint>(eFile);
    int nbCol = EXT.cols();
    int nbRow = EXT.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    MyMatrix<Tint> EXT_can = LinPolytope_CanonicForm<Tint, use_scheme>(EXT);
    std::cerr
        << "------------------------------------------------------------\n";
    //
    auto get_random_equivalent =
        [](MyMatrix<Tint> const &eMat) -> MyMatrix<Tint> {
      int nbRow = eMat.rows();
      int n = eMat.cols();
      std::vector<int> ePerm = RandomPermutation<int>(nbRow);
      MyMatrix<Tint> eUnimod = RandomUnimodularMatrix<Tint>(n);
      MyMatrix<Tint> eProd = eMat * eUnimod;
      MyMatrix<Tint> RetMat(nbRow, n);
      for (int iRow = 0; iRow < nbRow; iRow++) {
        int jRow = ePerm[iRow];
        RetMat.row(iRow) = eProd.row(jRow);
      }
      return RetMat;
    };
    //
    int n_iter = 10;
    for (int i_iter = 0; i_iter < n_iter; i_iter++) {
      std::cerr << "i_iter=" << i_iter << " / " << n_iter << "\n";
      MyMatrix<Tint> EXT2 = get_random_equivalent(EXT);
      MyMatrix<Tint> EXT2_can = LinPolytope_CanonicForm<Tint, use_scheme>(EXT2);
      std::cerr
          << "------------------------------------------------------------\n";
      if (!TestEqualityMatrix(EXT_can, EXT2_can)) {
        std::cerr << "Inconsistency in the canonical code\n";
        std::cerr << "EXT_Can=\n";
        WriteMatrix(std::cerr, EXT_can);
        std::cerr << "EXT2_Can=\n";
        WriteMatrix(std::cerr, EXT2_can);
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_LinPolytope_Canonic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
