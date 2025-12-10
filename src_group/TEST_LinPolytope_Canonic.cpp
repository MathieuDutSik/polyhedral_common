// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PolytopeEquiStab.h"
// clang-format on

template <typename T>
void test_canonic(MyMatrix<T> const &EXT, size_t threshold) {
  MyMatrix<T> EXT_can = LinPolytope_CanonicForm<T>(EXT, threshold, std::cerr);
  std::cerr << "-----------------------------------------------------\n";
  //
  auto get_random_equivalent = [](MyMatrix<T> const &eMat) -> MyMatrix<T> {
    int nbRow = eMat.rows();
    int n = eMat.cols();
    std::vector<int> ePerm = RandomPermutation<int>(nbRow);
    MyMatrix<T> eUnimod = RandomUnimodularMatrix<T>(n);
    MyMatrix<T> eProd = eMat * eUnimod;
    MyMatrix<T> RetMat(nbRow, n);
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
    MyMatrix<T> EXT2 = get_random_equivalent(EXT);
    MyMatrix<T> EXT2_can =
        LinPolytope_CanonicForm<T>(EXT2, threshold, std::cerr);
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
}

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
    std::string eFile = argv[1];
    MyMatrix<Tint> EXT = ReadMatrixFile<Tint>(eFile);

    test_canonic(EXT, THRESHOLD_USE_SUBSET_SCHEME_CANONIC);
    test_canonic(EXT, THRESHOLD_USE_SUBSET_SCHEME_TEST_CANONIC);
    //
    std::cerr << "Normal termination of TEST_LinPolytope_Canonic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_LinPolytope_Canonic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
