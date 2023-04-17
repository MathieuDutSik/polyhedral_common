// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "MatrixCanonicalForm.h"
#include "NumberTheory.h"

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_InvariantZBasisShortVectorNoGroup  [GramI] OutFile\n";
      std::cerr << "or\n";
      std::cerr << "LATT_InvariantZBasisShortVectorNoGroup  [GramI]\n";
      std::cerr << "\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      return -1;
    }
    using T = mpz_class;
    //    using T=long;
    // using T=mpq_class;

    // using Tint=long;
    using Tint = mpz_class;
    //
    std::ifstream is(argv[1]);
    MyMatrix<T> eMat = ReadMatrix<T>(is);
    MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat);
    auto do_print = [&](std::ostream &os) -> void {
      os << "SHV=\n";
      WriteMatrix(os, SHV);
    };
    if (argc == 2) {
      do_print(std::cerr);
    } else {
      std::string FileO = argv[2];
      std::ofstream os(FileO);
      do_print(os);
    }
    std::cerr << "Normal termination of LATT_canonicalize\n";
  } catch (TerminalException const &e) {
    std::cerr << "Raised exception led to premature end of LATT_canonicalize\n";
    exit(e.eVal);
  }
  runtime(time1);
}
