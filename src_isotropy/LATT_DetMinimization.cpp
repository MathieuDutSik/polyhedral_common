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
  MyMatrix<T> Q = ReadMatrixFile<T>(FileI);
  //
  ResultDetMin<T> res = DeterminantMinimization(Q, false, std::cerr);
  T det1 = DeterminantMat(Q);
  T det2 = DeterminantMat(res.Mred);
  if (OutFormat == "GAP") {
    os << "return rec(det1:=" << det1 << ", det2:=" << det2 << ", P:=";
    WriteMatrixGAP(os, res.P);
    os << ", Mred:=";
    WriteMatrixGAP(os, res.Mred);
    os << ");\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
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
      std::cerr
          << "LATT_DetMinimization arithmetic [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_DetMinimization arithmetic [FileI]\n";
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
    std::cerr << "Normal termination of LATT_DetMinimization\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_DetMinimization\n";
    exit(e.eVal);
  }
  runtime(time);
}
