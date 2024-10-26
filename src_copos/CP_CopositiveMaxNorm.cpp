// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

template<typename T, typename Tint>
void compute(std::string const& FileI, std::string const& strMaxNorm, std::string const& OutFormat, std::ostream& os) {
  MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
  T MaxNorm = ParseScalar<T>(strMaxNorm);
  //
  MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
  //
  CopositivityEnumResult<Tint> CopoRes =
    EnumerateCopositiveShortVector<T, Tint>(eSymmMat, InitialBasis,
                                            MaxNorm, std::cerr);
  //
  WriteCopositivityEnumResult(os, OutFormat, eSymmMat, CopoRes);
}


int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveMaxNorm [arith] [DATASYMM] [MaxNorm]\n";
      std::cerr << "or\n";
      std::cerr << "CP_CopositiveMaxNorm [arith] [DATASYMM] [MaxNorm] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: The arithmetic being used\n";
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
    //
    std::cerr << "Reading input\n";
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string strMaxNorm = argv[3];
    std::string OutFormat = "classic";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return compute<T,Tint>(FileI, strMaxNorm, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    //
    if (OutFile == "stderr") {
      f(std::cerr);
    } else {
      if (OutFile == "stdout") {
        f(std::cout);
      } else {
        std::ofstream os(OutFile);
        f(os);
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
