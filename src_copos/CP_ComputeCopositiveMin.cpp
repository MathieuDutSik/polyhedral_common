// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_ComputeCopositiveMin [DATASYMM]\n";
      std::cerr << "or\n";
      std::cerr << "CP_ComputeCopositiveMin [DATASYMM] [InitialBasois]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The input data of the strict copositive "
                   "symmetric matrix A\n";
      std::cerr << "Returns the copositive minimum min_{COP}(A) = min_{v in "
                   "Z^n_{>=0}} A[v] "
                   "of A and a list of all v in Z^n_{>= 0} such that A[v] = "
                   "min_{COP}(A) ";
      std::cerr << "\n";
      std::cerr << "If InitialBasis is not put in argument, then it is the "
                   "standard basis {e1, ...., en}\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    using T = mpq_class;
    using Tint = int;
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(argv[1]);
    std::cerr << "eSymmMat=\n";
    WriteMatrix(std::cerr, eSymmMat);
    //
    MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
    if (argc == 3)
      InitialBasis = ReadMatrixFile<Tint>(argv[2]);
    //
    Tshortest<T, Tint> eSh =
        T_CopositiveShortestVector<T, Tint>(eSymmMat, InitialBasis, std::cerr);
    std::cerr << "eMin=" << eSh.eMin << "\n";
    WriteMatrix(std::cerr, eSh.SHV);
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_ComputeCopositiveMin\n";
    exit(e.eVal);
  }
  runtime(time1);
}
