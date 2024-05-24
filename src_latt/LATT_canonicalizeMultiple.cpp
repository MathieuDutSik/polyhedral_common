// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
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
    Canonic_PosDef<T, Tint> RetF =
      ComputeCanonicalFormMultiple<T, Tint>(ListMatrix, std::cerr);
    //
    if (opt == 1) {
      std::ofstream os(argv[3]);
      WriteMatrix(os, RetF.Mat);
    }
    if (opt == 2) {
      std::ofstream os(argv[3]);
      os << "return rec(Basis:=";
      WriteMatrixGAP(os, RetF.Basis);
      os << ", SHV:=";
      WriteMatrixGAP(os, TransposedMat(RetF.SHV));
      os << ", eG:=";
      WriteMatrixGAP(os, RetF.Mat);
      os << ");\n";
    }
    std::cerr << "Normal termination of LATT_canonicalizeMultiple\n";
  } catch (TerminalException const &e) {
    std::cerr << "Raised exception led to premature end of "
                 "LATT_canonicalizeMultiple\n";
    exit(e.eVal);
  }
  runtime(time);
}
