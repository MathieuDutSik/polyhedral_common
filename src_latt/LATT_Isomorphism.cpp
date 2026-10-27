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
#include "LatticePleskenSouvignier.h"
#include "LatticeAutomorphismVinberg.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void ComputeIsomorphism(std::string const &FileListMat1,
                        std::string const &FileListMat2,
                        std::string const &OutFormat,
                        std::string const &method, std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  std::vector<MyMatrix<T>> ListMat1 = ReadListMatrixFile<T>(FileListMat1);
  std::vector<MyMatrix<T>> ListMat2 = ReadListMatrixFile<T>(FileListMat2);
  if (ListMat1.empty()) {
    std::cerr << "LATT_Isomorphism: The input matrix list in " << FileListMat1
              << " is empty\n";
    throw TerminalException{1};
  }
  if (ListMat2.empty()) {
    std::cerr << "LATT_Isomorphism: The input matrix list in " << FileListMat2
              << " is empty\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(ListMat1[0]) ||
      !IsPositiveDefinite(ListMat1[0], std::cerr)) {
    std::cerr << "LATT_Isomorphism: The first input Gram matrix in "
              << FileListMat1 << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  if (!IsSymmetricMatrix(ListMat2[0]) ||
      !IsPositiveDefinite(ListMat2[0], std::cerr)) {
    std::cerr << "LATT_Isomorphism: The first input Gram matrix in "
              << FileListMat2 << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  std::optional<MyMatrix<Tint>> equiv;
  if (method == "graph") {
    equiv = ArithmeticEquivalenceMultiple<T, Tint, Tgroup>(ListMat1, ListMat2,
                                                           std::cerr);
  } else if (method == "plesken_souvignier") {
    if (ListMat1.size() != ListMat2.size()) {
      std::cerr << "LATT_Isomorphism: the two lists have different lengths\n";
      throw TerminalException{1};
    }
    std::vector<MyMatrix<Tint>> ListGenAut2;
    equiv = PleskenSouvignierLatticeIsometry<T, Tint>(ListMat1, ListMat2,
                                                      ListGenAut2, std::cerr);
  } else if (method == "vinberg" || method == "weyl_chamber") {
    // The reflections in the roots of norm 2 are integral only for an
    // integral form, and the decomposition is the one of a single form.
    if (ListMat1.size() != 1 || ListMat2.size() != 1 ||
        !IsIntegralMatrix(ListMat1[0]) || !IsIntegralMatrix(ListMat2[0])) {
      std::cerr << "LATT_Isomorphism: the " << method << " method needs a "
                << "single integral Gram matrix on each side\n";
      throw TerminalException{1};
    }
    if (method == "vinberg") {
      equiv = ComputeIsometryVinberg<T, Tint, Tgroup>(ListMat1[0],
                                                      ListMat2[0], std::cerr);
    } else {
      equiv = ComputeIsometryWeylChamber<T, Tint, Tgroup>(
          ListMat1[0], ListMat2[0], std::cerr);
    }
  } else {
    std::cerr << "LATT_Isomorphism: unknown method=" << method << "\n";
    std::cerr << "The methods are graph, plesken_souvignier, vinberg and "
                 "weyl_chamber\n";
    throw TerminalException{1};
  }
  if (OutFormat == "GAP") {
    if (equiv) {
      os << "return ";
      WriteMatrixGAP(os, *equiv);
      os << ";\n";
    } else {
      os << "return false;\n";
    }
    return;
  }
  if (OutFormat == "Oscar") {
    if (equiv) {
      WriteMatrix(os, *equiv);
    } else {
      os << "0 0\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4 && argc != 6 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Isomorphism [arith] [ListMat1] [ListMat2] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Isomorphism [arith] [ListMat1] [ListMat2]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Isomorphism [arith] [ListMat1] [ListMat2] [OutFormat] "
                   "[OutFile] [method]\n";
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
      std::cerr << "  GAP   : the equivalence matrix or false, GAP readable\n";
      std::cerr << "  Oscar : the equivalence matrix in the Oscar format\n";
      std::cerr << "method values:\n";
      std::cerr << "  graph              : the isomorphism of vector families "
                   "as graphs\n";
      std::cerr << "                       (default)\n";
      std::cerr << "  plesken_souvignier : the backtrack of Plesken and "
                   "Souvignier\n";
      std::cerr << "  vinberg            : through the root system, for a "
                   "single integral\n";
      std::cerr << "                       Gram matrix on each side\n";
      std::cerr << "  weyl_chamber       : as vinberg, on the Weyl-chamber "
                   "families\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileListMat1 = argv[2];
    std::string FileListMat2 = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::string method = "graph";
    if (argc >= 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    if (argc == 7) {
      method = argv[6];
    }
    //
    auto prt = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            ComputeIsomorphism<T>(FileListMat1, FileListMat2, OutFormat,
                                  method, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of LATT_Isomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
