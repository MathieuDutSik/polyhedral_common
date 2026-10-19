// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "DeterminantMinimization.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> Q1 = ReadMatrixFile<T>(FileI);
  //
  auto f_print = [&](MyMatrix<T> const &Qinp) -> void {
    T det1 = DeterminantMat(Qinp);
    ResultDetMin<T> res = DeterminantMinimization(Qinp, true, std::cerr);
    T det2 = DeterminantMat(res.Mred);
    if (OutFormat == "GAP") {
      os << "return rec(det1:=" << det1 << ", det2:=" << det2 << ", P:=";
      WriteMatrixGAP(os, res.P);
      os << ", Mred:=";
      WriteMatrixGAP(os, res.Mred);
      os << ");\n";
      return;
    }
    if (OutFormat == "CPP") {
      WriteMatrix(os, res.Mred);
      return;
    }
    std::cerr << "Failed to find a matching OutFormat\n";
    throw TerminalException{1};
  };
  MyMatrix<T> Qinv = Inverse(Q1);
  MyMatrix<T> Q2 = RemoveFractionMatrix(Qinv);
  f_print(Q1);
  f_print(Q2);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string FileI;
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics(dispatch, [&]<typename T>() -> void {
      process<T>(FileI, OutFormat, *os_ptr);
    });
    if (argc != 3 && argc != 5) {
      std::cerr << "LATT_AutomorphyReduction arithmetic [FileI] [OutFormat] "
                   "[FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_AutomorphyReduction arithmetic [FileI]\n";
      std::cerr << "or\n";
      std::cerr << "Determining if two matrices are equivalent and computing "
                   "the automorphism group\n";
      std::cerr << "are tough problems. Two tools may help:\n";
      std::cerr << "* The equivalence for the inverse may be easier\n";
      std::cerr << "* The determinant may be reduceable and the problem easier "
                   "to solve for the reduced form\n";
      std::cerr << "\n";
      std::cerr << "     ------- arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    FileI = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    //
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileO, f);
    std::cerr << "Normal termination of LATT_AutomorphyReduction\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_AutomorphyReduction\n";
    exit(e.eVal);
  }
  runtime(time);
}
