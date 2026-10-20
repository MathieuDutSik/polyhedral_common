// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Shvec_sieve.h"
#include "SignatureSymmetric.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

/*
  The minimal vectors of a positive definite Gram matrix by the Gauss sieve
  of Shvec_sieve.h followed by the certifying enumeration of Shvec_exact.h.
  The result is exact whatever the sieve found; the sieve only decides the
  bound of the enumeration. The rows of SHV come in pairs v, -v.
 */

template <typename T, typename Tint>
void process(std::string const &FileGram, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileGram);
  if (!IsSymmetricMatrix(GramMat) || !IsPositiveDefinite(GramMat, std::cerr)) {
    std::cerr << "LATT_ShvecSieve: the matrix is not symmetric positive "
                 "definite\n";
    throw TerminalException{1};
  }
  SieveOptions<T> opts;
  Tshortest<T, Tint> rec = T_ShortestVectorSieve<T, Tint>(GramMat, opts,
                                                          std::cerr);
  if (OutFormat == "GAP") {
    os << "return rec(min:=" << rec.min << ", SHV:=";
    WriteMatrixGAP(os, rec.SHV);
    os << ");\n";
    return;
  }
  if (OutFormat == "CPP") {
    os << rec.min << "\n";
    WriteMatrix(os, rec.SHV);
    return;
  }
  std::cerr << "LATT_ShvecSieve: No matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string FileGram;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library(
        dispatch, [&]<typename T, typename Tint>() -> void {
          process<T, Tint>(FileGram, OutFormat, *os_ptr);
        });
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_ShvecSieve [arith] [GramFile] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_ShvecSieve [arith] [GramFile]\n";
      std::cerr << "\n";
      std::cerr << "GramFile  : a symmetric positive definite matrix\n";
      std::cerr << "OutFormat : GAP (default) or CPP\n";
      std::cerr << "OutFile   : the output file, or stderr (default) / "
                   "stdout\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    std::string arith = argv[1];
    FileGram = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of LATT_ShvecSieve\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_ShvecSieve\n";
    exit(e.eVal);
  }
  runtime(time);
}
