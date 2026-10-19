/* sv.c  simple driver for shvec                              */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#include "Shvec_exact.h"
#include "SignatureSymmetric.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &choice, std::string const &FileGram,
             std::string const &FileVect, std::string const &OutFormat,
             std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileGram);
  if (!IsSymmetricMatrix(GramMat) ||
      !IsPositiveDefinite(GramMat, std::cerr)) {
    std::cerr << "LATT_near: The input Gram matrix in " << FileGram
              << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  int n = GramMat.rows();
  auto get_vect = [&]() -> MyVector<T> {
    if (FileVect == "zero") {
      return ZeroVector<T>(n);
    }
    return ReadVectorFile<T>(FileVect);
  };
  MyVector<T> eV = get_vect();
  auto get_result = [&]() -> MyMatrix<Tint> {
    if (choice == "nearest") {
      resultCVP<T, Tint> res = NearestVectors<T, Tint>(GramMat, eV, std::cerr);
      return res.ListVect;
    }
    std::optional<std::string> opt_near = get_postfix(choice, "near=");
    if (opt_near) {
      T norm = ParseScalar<T>(*opt_near);
      MyVector<T> fV = -eV;
      std::vector<MyVector<Tint>> ListVect =
          FindAtMostDistVectors<T, Tint>(GramMat, fV, norm, std::cerr);
      return MatrixFromVectorFamilyDim(n, ListVect);
    }
    std::optional<std::string> opt_fixed = get_postfix(choice, "fixed=");
    if (opt_fixed) {
      T norm = ParseScalar<T>(*opt_fixed);
      MyVector<T> fV = -eV;
      std::vector<MyVector<Tint>> ListVect =
          FindFixedDistVectors<T, Tint>(GramMat, fV, norm, std::cerr);
      return MatrixFromVectorFamilyDim(n, ListVect);
    }
    std::cerr << "Failed to find a matching entry for choice=" << choice
              << "\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> result = get_result();
  auto write_result = [&]() -> void {
    if (OutFormat == "GAP") {
      os << "return ";
      WriteMatrixGAP(os, result);
      os << ";\n";
      return;
    }
    if (OutFormat == "Oscar") {
      WriteMatrix(os, result);
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  };
  write_result();
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_near arith choice [FileGram] [FileV] [OutFormat] "
                   "[FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_near arith choice [FileGram] [FileV]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic (see below)\n";
      std::cerr << "choice    : The choose option (see below)\n";
      std::cerr << "FileGram  : The list of inequalities\n";
      std::cerr << "FileV     : The vector for which we want to\n";
      std::cerr << "       compute the distance\n";
      std::cerr << "OutFormat : The format of output, GAP or Oscar\n";
      std::cerr << "FileOut   : The file of output (if present, otherwise "
                   "std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "gmp         : T=mpq_class, Tint=mpz_class\n";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "gmp_boost   : the boost bindings over gmp\n";
#endif
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "multi_boost : the boost multiprecision integers\n";
#endif
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "flint       : T=fmpq_class, Tint=fmpz_class\n";
#endif
      std::cerr << "\n";
      std::cerr << "        --- choice ---\n";
      std::cerr << "\n";
      std::cerr << "nearest   : The nearest points to that vector\n";
      std::cerr << "near=dist : The vector up to distance near from\n";
      std::cerr << "     that vector\n";
      std::cerr << "fixed=dist : The vector at an exact distance\n";
      std::cerr << "     from that vector\n";

      return -1;
    }
    //
    std::string arith = argv[1];
    std::string choice = argv[2];
    std::string FileGram = argv[3];
    std::string FileVect = argv[4];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      FileOut = argv[6];
    }

    auto call_SV = [&](std::ostream &os) -> void {
      ArithmeticDispatch dispatch;
      add_rational_arithmetics_by_library(
          dispatch, [&]<typename T>() -> void {
            process<T>(choice, FileGram, FileVect, OutFormat, os);
          });
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileOut, call_SV);
    std::cerr << "Normal termination of LATT_near\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_near\n";
    exit(e.eVal);
  }
  runtime(time);
}
