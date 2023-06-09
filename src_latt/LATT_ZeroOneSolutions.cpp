// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "SimulDiophantApprox.h"
// clang-format on

template <typename T, typename F>
void DoProcessing(MyMatrix<T> const &M, MyVector<T> const &B, F f) {
  MyMatrix<T> Mtr = TransposedMat(M);
  std::option<MyVector<T>> opt = SolutionIntMat(Mtr, B);
  if (!opt)
    return;
  MyMatrix<T> NSP = NullspaceIntMat(Mtr);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ZeroOneSolutions [FileIn] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "FileIn  : The input file in the same format as of Vedran "
                   "Krcadinac\n";
      std::cerr << "FileOut : The output fileepsilon value on input\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::string FileI = argv[1];
    std::string FileO = argv[1];
    std::ifstream is(FileI);
    std::ofstream os(FileO);

    int nbRow, nbCol, contVal;
    is >> nbRow;
    is >> nbCol;
    is >> contVal;
    MyMatrix<T> M(nbRow, nbCol);
    MyVector<T> B(nbRow);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      T val;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        is >> val;
        M(iRow, iCol) = val;
      }
      is >> val;
      B(iRow) = val;
    }

    auto f = [&](Face f) -> void {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        int val = f[iCol];
        os << val;
      }
      os << "\n";
    };
    DoProcessing(M, B, f);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ZeroOneSolutions\n";
    exit(e.eVal);
  }
  runtime(time);
}
