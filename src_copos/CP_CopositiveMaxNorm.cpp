// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveMaxNorm [DATASYMM] [MaxNorm]\n";
      std::cerr << "or\n";
      std::cerr << "CP_CopositiveMaxNorm [DATASYMM] [MaxNorm] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The symmetric matrix on input\n";
      std::cerr << "MaxNorm: The maximm norm considered\n";
      std::cerr << "OutFormat: classic or GAP. Default is classic\n";
      std::cerr << "OutFile: if assigned output goes to OutFile. Otherwise to "
                   "std::cerr\n";
      std::cerr << "\n";
      std::cerr << "It returns the list of integer vectors v in Z^n_{>= 0} "
                   "such that A[v] <= MaxNorm\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::cerr << "Reading input\n";
    std::string FileI = argv[1];
    std::string strMaxNorm = argv[2];
    std::string OutFormat = "classic";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
    T MaxNorm = ParseScalar<T>(strMaxNorm);
    //
    auto process = [&](std::ostream &os) -> void {
      MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
      //
      CopositivityEnumResult<Tint> CopoRes =
          EnumerateCopositiveShortVector<T, Tint>(eSymmMat, InitialBasis,
                                                  MaxNorm, std::cerr);
      //
      WriteCopositivityEnumResult(os, OutFormat, eSymmMat, CopoRes);
    };
    if (OutFile == "stderr") {
      process(std::cerr);
    } else {
      if (OutFile == "stdout") {
        process(std::cout);
      } else {
        std::ofstream os(OutFile);
        process(os);
      }
    }
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_CopositiveMaxNorm\n";
    exit(e.eVal);
  }
  runtime(time1);
}
