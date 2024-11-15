// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "IndefApproxCanonical.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "TestPermutationSignCanonic [FileI]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    using Tint = mpz_class;

    std::string FileI = argv[1];
    //
    MyMatrix<T> M = ReadMatrixFile<T>(argv[1]);
    if (!IsSymmetricMatrix(M)) {
      std::cerr << "The matrix M should be symmetric\n";
      throw TerminalException{1};
    }
    int nbRow = M.rows();
    //
    MyMatrix<T> Mcan =
        CanonicalizationPermutationSigns<T, Tint>(M, std::cerr).Mred;
    for (int iter = 0; iter < 50; iter++) {
      std::cerr << "iter=" << iter << "\n";
      std::vector<int> ePerm = RandomPermutation<int>(nbRow);
      MyMatrix<T> eP = ZeroMatrix<T>(nbRow, nbRow);
      for (int iRow = 0; iRow < nbRow; iRow++) {
        eP(iRow, ePerm[iRow]) = 2 * (random() % 2) - 1;
      }
      MyMatrix<T> Mnew = eP * M * eP.transpose();
      MyMatrix<T> Mnew_can =
          CanonicalizationPermutationSigns<T, Tint>(Mnew, std::cerr).Mred;
      if (Mcan != Mnew_can) {
        std::cerr << "The matrices are not equal\n";
        std::cerr << "Mcan=\n";
        WriteMatrix(std::cerr, Mcan);
        std::cerr << "Mnew_can=\n";
        WriteMatrix(std::cerr, Mnew_can);
        std::cerr << "This error is actually expected since if the matrix has "
                     "many entries of the\n";
        std::cerr << "Same absolute value (e.g. Hadamard matrix) then this "
                     "method does not work\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TestPermutationSignCanonic\n";
    exit(e.eVal);
  }
  runtime(time);
}
