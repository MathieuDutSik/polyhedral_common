// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveKernelMin [DATASYMM] [DATABASIS] [MaxNorm]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The input data of the symmetric matrix\n";
      std::cerr << "DATABASIS: The input data of the matrix basis\n";
      std::cerr << "It returns the copositive minimum vectors assuming that "
                   "the basis give non-negative pairwise scalar products\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    using T = mpq_class;
    using Tint = mpz_class;
    //
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(argv[1]);
    //
    MyMatrix<Tint> TheBasis = ReadMatrixFile<Tint>(argv[2]);
    //
    int MaxNorm_i;
    sscanf(argv[3], "%d", &MaxNorm_i);
    T MaxNorm = MaxNorm_i;
    //
    std::vector<MyVector<Tint>> LVect =
        EnumerateShortVectorInCone_UnderPositivityCond<T, Tint>(
            eSymmMat, TheBasis, MaxNorm);
    std::cerr << "|LVect|=" << LVect.size() << "\n";
    for (auto &eVect : LVect) {
      std::cerr << "eVect=";
      WriteVector(std::cerr, eVect);
    }
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
