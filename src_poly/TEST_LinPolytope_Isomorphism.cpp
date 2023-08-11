// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_LinPolytope_Canonic [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that "
                   "matter)\n";
      return -1;
    }
    //
    using Tint = mpz_class;
    using Tidx = uint16_t;
    const bool use_scheme = true;
    std::string eFile = argv[1];
    MyMatrix<Tint> EXT1 = ReadMatrixFile<Tint>(eFile);
    int nbCol = EXT1.cols();
    int nbRow = EXT1.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
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
      MyMatrix<Tint> EXT2 = get_random_equivalent(EXT1);
      std::optional<std::vector<Tidx>> equiv =
          LinPolytope_Isomorphism<Tint, Tidx, use_scheme>(EXT1, EXT2);
      if (!equiv) {
        std::cerr << "The isomorhism check return wrong results\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_LinPolytope_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
