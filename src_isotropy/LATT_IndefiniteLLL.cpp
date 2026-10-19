// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Indefinite_LLL.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> M = ReadMatrixFile<T>(FileI);

  ResultIndefiniteLLL<T, Tint> res = Indefinite_LLL<T, Tint>(M, std::cerr);
  if (OutFormat == "GAP") {
    os << "return rec(B:=" << StringMatrixGAP(res.B);
    os << ", Mred:=" << StringMatrixGAP(res.Mred);
    if (res.Xisotrop) {
      MyVector<T> const &Xisotrop = *res.Xisotrop;
      os << ", Xisotrop:=" << StringVectorGAP(Xisotrop);
    }
    os << ");\n";
    return;
  }
  if (OutFormat == "CPP") {
    os << "B_T=\n";
    WriteMatrix(os, res.B);
    os << "Mred=\n";
    WriteMatrix(os, res.Mred);
    if (res.Xisotrop) {
      MyVector<T> const &Xisotrop = *res.Xisotrop;
      WriteVector(os, Xisotrop);
    } else {
      os << "no isotrop vector found\n";
    }
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
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      process<T>(FileI, OutFormat, *os_ptr);
    });
    if (argc != 3 && argc != 5) {
      std::cerr
          << "LATT_IndefiniteLLL arithmetic [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_IndefiniteLLL arithmetic [FileI]\n";
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
    std::cerr << "Normal termination of LATT_IndefiniteLLL\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_IndefiniteLLL\n";
    exit(e.eVal);
  }
  runtime(time);
}
