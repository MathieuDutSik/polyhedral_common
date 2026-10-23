// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "LatticeStabEquiCan.h"
#include "SignatureSymmetric.h"
#include "LatticeRootDecomposition.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ComputeCanonical(std::string const &FileI, std::string const &OutFormat,
                      std::string const &method, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> eMat = ReadMatrixFile<T>(FileI);
  if (!IsSymmetricMatrix(eMat) || !IsPositiveDefinite(eMat, std::cerr)) {
    std::cerr << "LATT_Canonicalize: The input Gram matrix in " << FileI
              << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<Tint> B;
  if (method == "graph") {
    B = ComputeCanonicalForm<T, Tint, Tgroup>(eMat, std::cerr);
  } else if (method == "vinberg") {
    // The roots are the vectors of norm 2, which is the decomposition of an
    // integral form.
    if (!IsIntegralMatrix(eMat)) {
      std::cerr << "LATT_Canonicalize: the vinberg method needs an integral "
                   "Gram matrix\n";
      throw TerminalException{1};
    }
    B = ComputeCanonicalFormRootDecomposed<T, Tint, Tgroup>(eMat, std::cerr);
  } else {
    std::cerr << "LATT_Canonicalize: unknown method=" << method << "\n";
    std::cerr << "The methods are graph and vinberg\n";
    throw TerminalException{1};
  }
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
    if (argc != 3 && argc != 5 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Canonicalize [arith] [GramI] [OutFormat] [OutFile] "
                   "[method]\n";
      std::cerr << "\n";
      std::cerr << "GramI (input) : The gram matrix on input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      std::cerr << "method values:\n";
      std::cerr << "  graph   : from the characteristic vector family "
                   "(default)\n";
      std::cerr << "  vinberg : along the root decomposition, for an "
                   "integral Gram matrix\n";
      std::cerr << "The canonical forms of two methods are not comparable "
                   "with each other.\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "CPP";
    std::string OutFile = "stderr";
    std::string method = "graph";
    if (argc >= 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    if (argc == 6) {
      method = argv[5];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            ComputeCanonical<T>(FileI, OutFormat, method, os);
          });
      dispatch.run(arith);
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
