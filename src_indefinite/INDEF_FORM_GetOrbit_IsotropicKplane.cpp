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
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &MatFile, std::string const &KStr,
             std::string choice, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
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
    std::string MatFile;
    std::string KStr;
    std::string choice;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      process<T>(MatFile, KStr, choice, OutFormat, *os_ptr);
    });
    if (argc != 5 && argc != 7) {
      std::cerr << "INDEF_FORM_GetOrbit_IsotropicKplane [arith] [MatFile] "
                   "[k] [choice]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_GetOrbit_IsotropicKplane [arith] [MatFile] "
                   "[k] [choice] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    MatFile = argv[2];
    KStr = argv[3];
    choice = argv[4];
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
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
