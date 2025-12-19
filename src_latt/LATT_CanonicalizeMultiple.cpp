// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "LatticeStabEquiCan.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_canonicalizeMultiple opt [FileMatrices] OutFile\n";
      std::cerr << "\n";
      std::cerr << "If opt=1 then only the matrix is in output.\n";
      std::cerr << "If opt=2 then the basis, list of vectors and matrix is in "
                   "GAP formatted output\n";
      std::cerr << "FileGrams (input) : The list of matrices used on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      return -1;
    }
    using T = mpz_class;
    //    using T=long;
    // using T=mpq_class;

    // using Tint=long;
    using Tint = mpz_class;
    //
    int opt;
    sscanf(argv[1], "%d", &opt);
    //
    std::ifstream is(argv[2]);
    std::vector<MyMatrix<T>> ListMatrix;
    int nbMat;
    is >> nbMat;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<T> eMat = ReadMatrix<T>(is);
      ListMatrix.push_back(eMat);
    }
    MyMatrix<Tint> B = ComputeCanonicalFormMultiple<T, Tint>(ListMatrix, std::cerr);
    MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
    MyMatrix<T> Mat_red = B_T * ListMatrix[0] * B_T.transpose();
    //
    if (opt == 1) {
      std::ofstream os(argv[3]);
      WriteMatrix(os, Mat_red);
    }
    if (opt == 2) {
      std::ofstream os(argv[3]);
      os << "return rec(Basis:=";
      WriteMatrixGAP(os, B);
      os << ", eG:=";
      WriteMatrixGAP(os, Mat_red);
      os << ");\n";
    }
    std::cerr << "Normal termination of LATT_canonicalizeMultiple\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_canonicalizeMultiple\n";
    exit(e.eVal);
  }
  runtime(time);
}
