// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "lorentzian_linalg.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "COXDYN_GetFacetOneDomain [FileI] [FileO]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    //
    std::string FileI = argv[1];
    std::vector<MyVector<T>> l_vect;
    {
      MyMatrix<T> M = ReadMatrixFile<T>(FileI);
      int n_rows = M.rows();
      for (int i_row = 0; i_row < n_rows; i_row++) {
        MyVector<T> V = GetMatrixRow(M, i_row);
        l_vect.push_back(V);
      }
    }
    std::vector<MyVector<T>> l_vect_red = GetFacetOneDomain(l_vect);
    //
    std::string FileO = argv[2];
    MyMatrix<T> Mred = MatrixFromVectorFamily(l_vect_red);
    std::ofstream os(FileO);
    os << "return ";
    WriteMatrixGAP(os, Mred);
    os << ";\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in COXDYN_GetFacetOneDomain\n";
    exit(e.eVal);
  }
  runtime(time);
}
