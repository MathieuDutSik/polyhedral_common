// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "lorentzian_linalg.h"

int main(int argc, char *argv[]) {
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
      for (int i_row=0; i_row<n_rows; i_row++) {
        MyVector<T> V = GetMatrixRow(M, i_row);
        l_vect.push_back(V);
      }
    }
    std::vector<MyVector<T>> l_vect_red = GetFacetOneDomain(l_vect);
    //
    std::string FileO = argv[2];
    MyMatrix<T> Mred = MatrixFromVectorFamily(l_vect_red);
    WriteMatrixFile(FileO, Mred);
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
