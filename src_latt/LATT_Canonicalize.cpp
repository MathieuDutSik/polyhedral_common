// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
#include "SignatureSymmetric.h"
// clang-format on

template <typename T, typename Tint>
void ComputeCanonical(std::string const &FileI, std::string const &method,
                      std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> eMat = ReadMatrixFile<T>(FileI);
  if (!IsSymmetricMatrix(eMat) || !IsPositiveDefinite(eMat, std::cerr)) {
    std::cerr << "LATT_Canonicalize: The input Gram matrix in " << FileI
              << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  auto get_basis = [&]() -> MyMatrix<Tint> {
    if (method == "zbasis") {
      return ComputeCanonicalForm<T, Tint>(eMat, std::cerr);
    }
    if (method == "fullrank") {
      using Tidx = uint32_t;
      using Telt = permutalib::SingleSidedPerm<Tidx>;
      using TintGroup = mpz_class;
      using Tgroup = permutalib::Group<Telt, TintGroup>;
      return ComputeCanonicalFormFullRank<T, Tint, Tgroup>(eMat, std::cerr);
    }
    std::cerr << "LATT_Canonicalize: The method " << method
              << " is not among the supported ones: zbasis, fullrank\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> B = get_basis();
  MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
  MyMatrix<T> eMat_red = B_T * eMat * B_T.transpose();
  if (OutFormat == "CPP") {
    WriteMatrix(os, eMat_red);
    return;
  }
  if (OutFormat == "PYTHON") {
    os << "{\"Basis\":" << StringMatrixPYTHON(B)
       << ", \"eG\":" << StringMatrixPYTHON(eMat_red) << "}\n";
    return;
  }
  if (OutFormat == "GAP") {
    os << "return ";
    WriteMatrixGAP(os, eMat_red);
    os << ";\n";
    return;
  }
  if (OutFormat == "GAP_full") {
    os << "return rec(Basis:=";
    WriteMatrixGAP(os, B);
    os << ", eG:=";
    WriteMatrixGAP(os, eMat_red);
    os << ");\n";
    return;
  }
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc < 3 || argc > 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Canonicalize [arith] [method] [GramI]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Canonicalize [arith] [method] [GramI] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "method: zbasis (default) uses a vector family spanning "
                   "Z^n, fullrank a full rank one\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string method = "zbasis";
    std::string FileI;
    std::string OutFormat = "CPP";
    std::string OutFile = "stderr";
    if (argc == 3) {
      FileI = argv[2];
    }
    if (argc == 4) {
      method = argv[2];
      FileI = argv[3];
    }
    if (argc == 5) {
      FileI = argv[2];
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    if (argc == 6) {
      method = argv[2];
      FileI = argv[3];
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return ComputeCanonical<T, Tint>(FileI, method, OutFormat, os);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return ComputeCanonical<T, Tint>(FileI, method, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return ComputeCanonical<T, Tint>(FileI, method, OutFormat, os);
      }
      std::cerr << "Failed to find a matching entry for arith\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of LATT_Canonicalize\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Canonicalize\n";
    exit(e.eVal);
  }
  runtime(time);
}
