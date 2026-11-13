// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "Group.h"
#include "Permutation.h"
#include "LatticeMethods.h"
#include "SignatureSymmetric.h"
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
  LatticeMethodChoice emethod = ParseLatticeMethodChoice(method);
  MyMatrix<Tint> B =
      LatticeCanonicalBasis<T, Tint, Tgroup>(eMat, emethod, std::cerr);
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
  std::cerr << "LATT_Canonicalize: unknown OutFormat=" << OutFormat << "\n";
  std::cerr << "The OutFormat values are CPP, PYTHON, GAP and GAP_full\n";
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
      std::cerr << "arith values:\n";
      std::cerr << "  gmp         : mpq_class / mpz_class (default choice)\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  gmp_boost   : the boost bindings to the gmp types\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "  multi_boost : the boost multiprecision types\n";
#endif
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "  flint       : fmpq_class / fmpz_class\n";
#endif
#ifdef QUADRATIC_FIELDS
      std::cerr << "  Qsqrt2      : the form over the field Q(sqrt(2))\n";
#endif
#ifdef QUADRATIC_FIELDS
      std::cerr << "  Qsqrt5      : the form over the field Q(sqrt(5))\n";
#endif
      std::cerr << "  RealAlgebraic=FileDesc : the form over the real\n";
      std::cerr << "              algebraic field described in FileDesc\n";
      std::cerr << "              (the lattice is Z^n in every case)\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      std::cerr << "OutFormat values:\n";
      std::cerr << "  CPP      : the canonical Gram matrix (default)\n";
      std::cerr << "  PYTHON   : the basis and the canonical Gram matrix\n";
      std::cerr << "  GAP      : the canonical Gram matrix\n";
      std::cerr << "  GAP_full : the basis and the canonical Gram matrix\n";
      std::cerr << "method values (default graph):\n";
      PrintLatticeStabCanIsomMethods(std::cerr);
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
      auto f_arith = [&]<typename T>() -> void {
        ComputeCanonical<T>(FileI, OutFormat, method, os);
      };
      add_rational_arithmetics_by_library(dispatch, f_arith);
      add_algebraic_field_arithmetics(dispatch, f_arith);
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
