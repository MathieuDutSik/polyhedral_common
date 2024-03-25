// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc < 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveMin [DATASYMM] [MaxNorm]\n";
      std::cerr << "CP_CopositiveMin [DATASYMM] [MaxNorm] [OutFormat]\n";
      std::cerr
          << "CP_CopositiveMin [DATASYMM] [MaxNorm] [OutFormat] [OutFile]\n";
      std::cerr << "CP_CopositiveMin [DATASYMM] [MaxNorm] [OutFormat] "
                   "[OutFile] [InitialBasis]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The symmetric matrix on input\n";
      std::cerr << "MaxNorm: The maximm norm considered\n";
      std::cerr << "OutFormat: classic or GAP. Default is classic\n";
      std::cerr << "OutFile: if assigned output goes to OutFile. Otherwise to "
                   "std::cerr\n";
      std::cerr << "InitialBasis: By default, identity mat, otherwise, the "
                   "matrix in the input file\n";
      std::cerr << "\n";
      std::cerr << "It returns the list of integer vectors v in Z^n_{>= 0} "
                   "such that A[v] <= MaxNorm\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::cerr << "Reading input\n";
    //
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(argv[1]);
    std::cerr << "eSymmMat=\n";
    WriteMatrix(std::cerr, eSymmMat);
    //
    int MaxNorm_i;
    sscanf(argv[2], "%d", &MaxNorm_i);
    T MaxNorm = MaxNorm_i;
    //
    std::string OutFormat = "classic";
    if (argc >= 4)
      OutFormat = argv[3];
    //
    auto process = [&](std::ostream &os) -> void {
      MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
      if (argc >= 6)
        InitialBasis = ReadMatrixFile<Tint>(argv[5]);
      //
      RequestCopositivity<T> CopoReq{MaxNorm, false};
      CopositivityEnumResult<Tint> CopoRes =
          EnumerateCopositiveShortVector<T, Tint>(eSymmMat, InitialBasis,
                                                  CopoReq, std::cerr);
      //
      WriteCopositivityEnumResult(os, OutFormat, eSymmMat, CopoRes);
    };
    if (argc >= 5) {
      std::string file = argv[4];
      std::ofstream os(file);
      process(os);
    } else {
      process(std::cerr);
    }
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_CopositiveMin\n";
    exit(e.eVal);
  }
  runtime(time1);
}
