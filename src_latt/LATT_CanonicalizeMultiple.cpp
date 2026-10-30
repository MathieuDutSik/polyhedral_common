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
#include "LatticeOrbitSplitCanonical.h"
#include "SignatureSymmetric.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ComputeCanonicalMultiple(std::string const &FileI,
                              std::string const &OutFormat,
                              std::string const &method, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListMatrix = ReadListMatrixFile<T>(FileI);
  if (ListMatrix.empty()) {
    std::cerr << "LATT_CanonicalizeMultiple: The input matrix list in "
              << FileI << " is empty\n";
    throw TerminalException{1};
  }
  auto f_canonical = [&]() -> MyMatrix<Tint> {
    if (method == "graph") {
      return ComputeCanonicalFormMultiple<T, Tint, Tgroup>(ListMatrix,
                                                           std::cerr);
    }
    if (method == "plesken_souvignier") {
      // The engine needs the first matrix positive definite, the others
      // are arbitrary.
      if (!IsSymmetricMatrix(ListMatrix[0]) ||
          !IsPositiveDefinite(ListMatrix[0], std::cerr)) {
        std::cerr << "LATT_CanonicalizeMultiple: the first matrix has to be "
                     "symmetric positive definite\n";
        throw TerminalException{1};
      }
      std::optional<MyMatrix<Tint>> opt =
          OrbitSplitCanonicalizeListMat<T, Tint, Tgroup>(ListMatrix,
                                                         std::cerr);
      if (!opt) {
        std::cerr << "LATT_CanonicalizeMultiple: the orbit splitting could "
                     "not select a canonical prefix, no canonical form for "
                     "this input\n";
        throw TerminalException{1};
      }
      return *opt;
    }
    std::cerr << "LATT_CanonicalizeMultiple: unknown method=" << method
              << "\n";
    std::cerr << "The methods are graph and plesken_souvignier\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> B = f_canonical();
  // The canonical form is the whole list of reduced matrices: two
  // inequivalent configurations can agree on the first one.
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  std::vector<MyMatrix<T>> ListMatRed;
  for (auto &eMat : ListMatrix) {
    ListMatRed.push_back(B_T * eMat * B_T.transpose());
  }
  if (OutFormat == "CPP") {
    WriteListMatrix(os, ListMatRed);
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(Basis:=";
    WriteMatrixGAP(os, B);
    os << ", ListMatRed:=";
    WriteListMatrixGAP(os, ListMatRed);
    os << ");\n";
    return;
  }
  std::cerr << "LATT_CanonicalizeMultiple: unknown OutFormat=" << OutFormat
            << "\n";
  std::cerr << "The OutFormat values are CPP and GAP\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 5 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_CanonicalizeMultiple [arith] [FileMatrices] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_CanonicalizeMultiple [arith] [FileMatrices]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_CanonicalizeMultiple [arith] [FileMatrices] "
                   "[OutFormat] [OutFile] [method]\n";
      std::cerr << "\n";
      std::cerr << "FileMatrices (input) : The list of matrices used on "
                   "input\n";
      std::cerr << "OutFile: The filename of the data in output\n";
      std::cerr << "\n";
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
      std::cerr << "OutFormat values:\n";
      std::cerr << "  CPP : the list of reduced matrices B * M_i * B^T\n";
      std::cerr << "  GAP : the basis B and the list of reduced matrices\n";
      std::cerr << "method values:\n";
      std::cerr << "  graph              : from the characteristic vector "
                   "family (default)\n";
      std::cerr << "  plesken_souvignier : orbit splitting over the "
                   "automorphism group\n";
      std::cerr << "                       computed by the Plesken-Souvignier "
                   "backtrack;\n";
      std::cerr << "                       the first matrix positive "
                   "definite\n";
      std::cerr << "The canonical forms of two methods are not comparable "
                   "with each other.\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    std::string method = "graph";
    if (argc >= 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    if (argc == 6) {
      method = argv[5];
    }
    //
    // A field is required: the canonicalization of a family that does not
    // span Z^n goes through LinPolytopeIntegral_Canonicalization_Subspaces,
    // which divides. Integral input is unaffected, a rational type prints an
    // integer as an integer.
    auto f = [&](std::ostream &os_out) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            ComputeCanonicalMultiple<T>(FileI, OutFormat, method, os_out);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_CanonicalizeMultiple\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_CanonicalizeMultiple\n";
    exit(e.eVal);
  }
  runtime(time);
}
