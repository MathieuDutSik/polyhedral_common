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
#include "edgewalk.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &MatFile, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  // The group order is counted in TintGroup, which is unrelated to the
  // arithmetic of the coordinates: fixing it to mpz_class keeps the
  // permutation group out of whichever arithmetic Tint happens to be.
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> LorMat = ReadMatrixFile<T>(MatFile);
  //
  ResultEdgewalk<T, Tint> re =
      StandardEdgewalkAnalysis<T, Tint, Tgroup>(LorMat, std::cerr);
  bool ComputeAllSimpleRoots = true;
  PrintResultEdgewalk(LorMat, re, os_out, OutFormat, ComputeAllSimpleRoots,
                      std::cerr);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string MatFile;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      process<T>(MatFile, OutFormat, *os_ptr);
    });
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "LORENTZ_ReflectiveEdgewalk [arith] [MatFile]\n";
      std::cerr << "     or\n";
      std::cerr << "LORENTZ_ReflectiveEdgewalk [arith] [MatFile] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    MatFile = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    auto f = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of LORENTZ_ReflectiveEdgewalk\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_ReflectiveEdgewalk\n";
    exit(e.eVal);
  }
  runtime(time);
}
