// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "NumberTheoryCommon.h"
#include "CombinedAlgorithms.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T,
          typename Tint = typename underlying_z_ring<T>::ring_type>
void process(std::string const &MatFile, std::string const &KStr,
             std::string choice, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  // The group order is counted in TintGroup, which is unrelated to the
  // arithmetic of the coordinates: fixing it to mpz_class keeps the
  // permutation group out of whichever arithmetic Tint happens to be.
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> Qmat = ReadMatrixFile<T>(MatFile);
  int k = ParseScalar<int>(KStr);
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  auto f_get = [&]() -> std::vector<MyMatrix<Tint>> {
    if (choice == "plane") {
      return comb.INDEF_FORM_GetOrbit_IsotropicKplane(Qmat, k);
    }
    if (choice == "flag") {
      return comb.INDEF_FORM_GetOrbit_IsotropicKflag(Qmat, k);
    }
    std::cerr << "No correct choice. choice=" << choice << "\n";
    throw TerminalException{1};
  };
  std::vector<MyMatrix<Tint>> l_planes = f_get();
  for (auto &e_plane : l_planes) {
    MyMatrix<T> e_plane_T = UniversalMatrixConversion<T, Tint>(e_plane);
    MyMatrix<T> prod = e_plane_T * Qmat * e_plane_T.transpose();
    if (!IsZeroMatrix(prod)) {
      std::cerr << "The plane is not an isotropic plane\n";
      throw TerminalException{1};
    }
  }
  if (OutFormat == "PYTHON") {
    return WriteListMatrixPYTHON(os_out, l_planes);
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteListMatrixGAP(os_out, l_planes);
    os_out << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 5 && argc != 7) {
      std::cerr << "INDEF_FORM_GetOrbit_IsotropicKplane [arith] [MatFile] "
                   "[k] [choice]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_GetOrbit_IsotropicKplane [arith] [MatFile] "
                   "[k] [choice] [OutFormat] [OutFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string MatFile = argv[2];
    std::string KStr = argv[3];
    std::string choice = argv[4];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return process<T>(MatFile, KStr, choice, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(MatFile, KStr, choice, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(MatFile, KStr, choice, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return process<T>(MatFile, KStr, choice, OutFormat, os);
      }
#endif
      std::cerr << "Failed to find matching type for arith\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Allowed values: gmp, flint (and, with ENABLE_BOOST_TYPES";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost";
#endif
      std::cerr << " and ";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << "multi_boost";
#endif
      std::cerr << ")\n";
#else
      std::cerr << "Allowed values: gmp (build with "
                << "ENABLE_FLINT_SUPPORT=1 for flint; with "
                << "ENABLE_BOOST_TYPES for the boost ones)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of INDEF_FORM_GetOrbit_IsotropicKplane\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_GetOrbit_IsotropicKplane runtime=" << time
              << "\n";
    exit(e.eVal);
  }
  runtime(time);
}
