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
void process(std::string const &FileM, std::string const &FileV1,
             std::string const &FileV2, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  // The group order is counted in TintGroup, which is unrelated to the
  // arithmetic of the coordinates: fixing it to mpz_class keeps the
  // permutation group out of whichever arithmetic Tint happens to be.
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> Q = ReadMatrixFile<T>(FileM);
  MyMatrix<Tint> v1 = ReadVectorFile<Tint>(FileV1);
  MyMatrix<Tint> v2 = ReadVectorFile<Tint>(FileV2);
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  std::optional<MyMatrix<Tint>> opt =
      comb.INDEF_FORM_EquivalenceVector(Q, Q, v1, v2);
  if (OutFormat == "PYTHON") {
    if (opt) {
      WriteMatrixPYTHON(os_out, *opt);
    } else {
      os_out << "None";
    }
    os_out << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    if (opt) {
      WriteMatrixGAP(os_out, *opt);
    } else {
      os_out << "fail";
    }
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
    std::string FileM;
    std::string FileV1;
    std::string FileV2;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      process<T>(FileM, FileV1, FileV2, OutFormat, *os_ptr);
    });
    if (argc != 5 && argc != 7) {
      std::cerr << "INDEF_FORM_TestEquivalenceVector [arith] [FileM] [FileV1] "
                   "[FileV2]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_TestEquivalenceVector [arith] [FileM] [FileV1] "
                   "[FileV2] [OutFormat] [OutFile]\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    FileM = argv[2];
    FileV1 = argv[3];
    FileV2 = argv[4];
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of INDEF_FORM_TestEquivalence\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_TestEquivalenceVector, runtime=" << time
              << "\n";
    exit(e.eVal);
  }
  runtime(time);
}
