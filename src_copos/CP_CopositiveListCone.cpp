// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveListCone [DATASYMM] [OUTFILE]\n";
      std::cerr << "\n";
      std::cerr
          << "DATASYMM: The input data of the copositive symmetric matrix\n";
      std::cerr << "It returns the list of cones that is used to test that the "
                   "matrix is copositive\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::cerr << "Reading input\n";
    //
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(argv[1]);
    //
    int n = eSymmMat.rows();
    std::cerr << "eSymmMat.rows=" << eSymmMat.rows()
              << " eSymmMat.cols=" << eSymmMat.cols() << "\n";
    MyMatrix<Tint> TheBasis = IdentityMat<Tint>(n);
    ResultListCone<Tint> res =
        EnumerateListConeCopositive<T, Tint>(eSymmMat, TheBasis, std::cerr);
    size_t nbCone = res.ListCone.size();
    std::cerr << "nbCone=" << nbCone << "\n";
    if (res.test == false) {
      std::cerr << "Matrix is not Copositive\n";
    } else {
      std::ofstream OUTfs(argv[3]);
      for (size_t iCone = 0; iCone < nbCone; iCone++) {
        OUTfs << "iCone=" << iCone << "/" << nbCone << " Basis\n";
        MyMatrix<Tint> const &eMat = res.ListCone[iCone];
        WriteMatrix(OUTfs, eMat);
      }
    }
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_CopositiveListCone\n";
    exit(e.eVal);
  }
  runtime(time1);
}
