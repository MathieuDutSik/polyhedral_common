// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

template<typename T, typename Tint>
void compute(std::string const& FileI, std::string const& OutFormat, std::ostream& os) {
  std::cerr << "Reading input\n";
  MyMatrix<T> eSymmMat = ReadMatrixFile<T>(FileI);
  std::cerr << "eSymmMat=\n";
  WriteMatrix(std::cerr, eSymmMat);

  MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
  //
  Tshortest<T, Tint> eSh =
    CopositiveShortestVector<T, Tint>(eSymmMat, InitialBasis, std::cerr);

  if (OutFormat == "clear") {
    os << "eMin=" << eSh.min << "\n";
    WriteMatrix(os, eSh.SHV);
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(eMin:=" << eSh.min << ", SHV:=";
    WriteMatrixGAP(os, eSh.SHV);
    os << ");\n";
    return;
  }
  std::cerr << "No matching for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}


int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveMin [arith] [DATASYMM]\n";
      std::cerr << "    or\n";
      std::cerr << "CP_CopositiveMin [arith] [DATASYMM] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "arith: The chosen arithmetic\n";
      std::cerr << "DATASYMM: The input data of the strict copositive "
                   "symmetric matrix A\n";
      std::cerr << "Returns the copositive minimum min_{COP}(A) = min_{v in "
                   "Z^n_{>=0}} A[v] "
                   "of A and a list of all v in Z^n_{>= 0} such that A[v] = "
                   "min_{COP}(A) ";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "classic";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
    auto f_print=[&](std::ostream & os_out) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return compute<T,Tint>(FileI, OutFormat, os_out);
      }
      std::cerr << "No matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    //
    print_stderr_stdout_file(OutFile, f_print);
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CP_CopositiveMin\n";
    exit(e.eVal);
  }
  runtime(time1);
}
