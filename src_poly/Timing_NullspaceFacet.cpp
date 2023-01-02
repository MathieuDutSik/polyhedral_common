// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "MatrixLinbox.h"
#include "NumberTheory.h"
#include "POLY_LinearProgramming.h"
#include "POLY_PolytopeFct.h"
#include <stdio.h>
int main(int argc, char *argv[]) {
  using T = mpq_class;
  try {
    if (argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceComputationLinbox [inputMat] [n_iter]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> EXT = ReadMatrix<T>(INmat);
    int n_col = EXT.cols();
    int n_iter;
    sscanf(argv[2], "%d", &n_iter);
    std::cerr << "n_iter=" << n_iter << "\n";

    vectface ListFace = FindVertices(EXT, n_iter);
    std::cerr << "We have the facets\n";

    auto f = [&](int method) -> void {
      std::chrono::time_point<std::chrono::system_clock> time1 =
          std::chrono::system_clock::now();
      int sumrank = 0;
      for (auto &eFace : ListFace) {
        size_t nb = eFace.count();
        MyMatrix<T> TheProv(nb, n_col);
        int aRow = eFace.find_first();
        for (size_t iRow = 0; iRow < nb; iRow++) {
          TheProv.row(iRow) = EXT.row(aRow);
          aRow = eFace.find_next(aRow);
        }
        MyMatrix<T> NSP;
        if (method == 1)
          NSP = NullspaceTrMat(TheProv);
        if (method == 2)
          NSP = NullspaceTrMat_linbox(TheProv);
        if (method == 3) {
          Eigen::FullPivLU<MyMatrix<T>> lu(TheProv.transpose());
          NSP = lu.kernel();
        }
        sumrank += NSP.rows();
      }
      std::chrono::time_point<std::chrono::system_clock> time2 =
          std::chrono::system_clock::now();
      std::cerr << "method=" << method << " time="
                << std::chrono::duration_cast<std::chrono::microseconds>(time2 -
                                                                         time1)
                       .count()
                << " sumrank=" << sumrank << "\n";
    };
    f(1);
    f(2);
    f(3);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
