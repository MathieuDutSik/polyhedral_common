// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "CombinedAlgorithms.h"
#include "Group.h"
#include "Permutation.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &QFile, std::string const &PlaneFile,
             std::string const &choice, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  // The group order is counted in TintGroup, which is unrelated to the
  // arithmetic of the coordinates: fixing it to mpz_class keeps the
  // permutation group out of whichever arithmetic Tint happens to be.
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> Qmat = ReadMatrixFile<T>(QFile);
  MyMatrix<Tint> Plane = ReadMatrixFile<Tint>(PlaneFile);
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  auto f_get = [&]() -> std::vector<MyMatrix<Tint>> {
    if (choice == "plane") {
      return comb.INDEF_FORM_Stabilizer_IsotropicKplane(Qmat, Plane);
    }
    if (choice == "flag") {
      return comb.INDEF_FORM_Stabilizer_IsotropicKflag(Qmat, Plane);
    }
    std::cerr << "No correct choice. choice=" << choice << "\n";
    throw TerminalException{1};
  };
  std::vector<MyMatrix<Tint>> l_gen = f_get();
  if (OutFormat == "CPP") {
    return WriteListMatrix(os_out, l_gen);
  }
  if (OutFormat == "PYTHON") {
    return WriteListMatrixPYTHON(os_out, l_gen);
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteListMatrixGAP(os_out, l_gen);
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
    std::string QFile;
    std::string PlaneFile;
    std::string choice;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      process<T>(QFile, PlaneFile, choice, OutFormat, *os_ptr);
    });
    if (argc != 5 && argc != 7) {
      std::cerr << "INDEF_FORM_StabilizerIsotropicPlane [arith] [QFile] "
                   "[PlaneFile] [choice]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_StabilizerIsotropicPlane [arith] [QFile] "
                   "[PlaneFile] [choice] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    QFile = argv[2];
    PlaneFile = argv[3];
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
    std::cerr << "Normal termination of INDEF_FORM_StabilizerIsotropicPlane\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_StabilizerIsotropicPlane, runtime="
              << time << "\n";
    exit(e.eVal);
  }
  runtime(time);
}
