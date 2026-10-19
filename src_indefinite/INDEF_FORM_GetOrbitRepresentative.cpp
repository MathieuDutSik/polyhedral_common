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
void process(std::string const &MatFile, std::string const &XnormStr,
             std::string const &OutFormat, std::ostream &os_out) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  // The group order is counted in TintGroup, which is unrelated to the
  // arithmetic of the coordinates: fixing it to mpz_class keeps the
  // permutation group out of whichever arithmetic Tint happens to be.
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> Qmat = ReadMatrixFile<T>(MatFile);
  T Xnorm = ParseScalar<T>(XnormStr);
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  std::vector<MyVector<Tint>> l_repr =
      comb.INDEF_FORM_GetOrbitRepresentative(Qmat, Xnorm);
  for (auto &e_repr : l_repr) {
    T norm = EvaluationQuadForm(Qmat, e_repr);
    if (norm != Xnorm) {
      std::cerr << "Wrong norm. norm=" << norm << " Xnorm=" << Xnorm << "\n";
      throw TerminalException{1};
    }
  }
  if (OutFormat == "PYTHON") {
    if (l_repr.empty()) {
      os_out << "[]";
    } else {
      WriteMatrixPYTHON(os_out, MatrixFromVectorFamily(l_repr));
    }
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    if (l_repr.empty()) {
      os_out << "[]";
    } else {
      WriteMatrixGAP(os_out, MatrixFromVectorFamily(l_repr));
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
    std::string MatFile;
    std::string XnormStr;
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      process<T>(MatFile, XnormStr, OutFormat, *os_ptr);
    });
    if (argc != 4 && argc != 6) {
      std::cerr
          << "INDEF_FORM_GetOrbitRepresentative [arith] [MatFile] [Xnorm] \n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_GetOrbitRepresentative [arith] [MatFile] "
                   "[Xnorm] [OutFormat] [OutFile]\n";
      std::cerr << "     -----\n";
      std::cerr << "arith can be: gmp\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    MatFile = argv[2];
    XnormStr = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    std::cerr << "Normal termination of INDEF_FORM_GetOrbitRepresentative\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_GetOrbitRepresentative, runtime=" << time
              << "\n";
    exit(e.eVal);
  }
  runtime(time);
}
