// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "POLY_PolytopeInt.h"
// clang-format on

template<typename T>
bool TestMatrixFile(std::string const& FileFAC) {
  MyMatrix<T> FAC = ReadMatrixFile<T>(FileFAC);
  int n_row = FAC.rows();
  int n_col = FAC.cols();
  MyVector<T> eCent = GetGeometricallyUniqueInteriorPoint(FAC, std::cerr);
  for (int i=0; i<20; i++) {
    std::vector<int> ePerm = RandomPermutation<int>(n_row);
    MyMatrix<T> eUnitMod = RandomUnimodularMatrix<T>(n_col);
    MyMatrix<T> eUnitMod_cgr = CongrMap(eUnitMod);
    MyMatrix<T> FACprod = FAC * eUnitMod;
    MyMatrix<T> NewFAC(n_row, n_col);
    for (int i_row = 0; i_row < n_row; i_row++) {
      int j_row = ePerm[i_row];
      NewFAC.row(i_row) = FACprod.row(j_row);
    }
    MyVector<T> eCentProd = GetGeometricallyUniqueInteriorPoint(NewFAC, std::cerr);
    MyVector<T> eCentMap = eUnitMod_cgr.transpose() * eCent;
    std::cerr << "eCentMap =" << StringVectorGAP(eCentMap) << "\n";
    std::cerr << "eCentProd=" << StringVectorGAP(eCentProd) << "\n";
    if (eCentMap != eCentProd) {
      return false;
    }
  }
  return true;
}


int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "TEST_GeometricallyUniquePoint [FAC]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_GeometricallyUniquePoint [FAC] [resultFile]\n";
      return -1;
    }
    using T = mpq_class;
    //
    //  std::cerr << "Reading input\n";
    //
    std::string FileFAC = argv[1];
    bool test = TestMatrixFile<T>(FileFAC);
    if (argc == 2) {
      if (!test) {
        std::cerr << "Inconsistent result found\n";
        throw TerminalException{1};
      }
    } else {
      std::string FileResult = argv[2];
      std::ofstream osF(FileResult);
      osF << "return " << GAP_logical(test) << ";\n";
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_GeometricallyUniquePoint\n";
    exit(e.eVal);
  }
  runtime(time);
}
