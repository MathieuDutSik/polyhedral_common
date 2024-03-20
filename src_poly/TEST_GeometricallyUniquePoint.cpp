// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "POLY_PolytopeInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "TEST_GeometricallyUniquePoint [FAC]\n";
      return -1;
    }
    using T = mpq_class;
    //
    //  std::cerr << "Reading input\n";
    //
    std::string FileFAC = argv[1];
    MyMatrix<T> FAC = ReadMatrixFile<T>(FileFAC);
    int n_row = FAC.rows();
    int n_col = FAC.cols();
    MyVector<T> eCent = GetGeometricallyUniqueInteriorPoint(FAC);
    for (int i=0; i<100; i++) {
      std::vector<int> ePerm = RandomPermutation<int>(n_row);
      MyMatrix<T> eUnitMod = RandomUnimodularMatrix<T>(n_col);
      MyMatrix<T> eUnitMod_cgr = CongrMap(eUnitMod);
      MyMatrix<T> FACprod = FAC * eUnitMod;
      MyMatrix<T> NewFAC(n_row, n_col);
      for (int i_row = 0; i_row < n_row; i_row++) {
        int j_row = ePerm[i_row];
        NewFAC.row(i_row) = FACprod.row(j_row);
      }
      MyVector<T> eCentProd = GetGeometricallyUniqueInteriorPoint(NewFAC);
      MyVector<T> eCentMap = eUnitMod_cgr.transpose() * eCent;
      if (eCentMap != eCentProd) {
        std::cerr << "Inconsistency in the transformation\n";
        throw TerminalException{1};
      }
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_GeometricallyUniquePoint\n";
    exit(e.eVal);
  }
  runtime(time);
}
