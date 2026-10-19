// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Isotropic.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &OutFormat,
             std::ostream &os_out) {
  MyMatrix<T> Q = ReadMatrixFile<T>(FileI);
  //
  std::optional<MyVector<T>> opt = FindIsotropic(Q, std::cerr);
  if (OutFormat == "GAP") {
    if (opt) {
      MyVector<T> const &eV = *opt;
      T val = EvaluationQuadForm(Q, eV);
      if (val != 0) {
        std::cerr << "LATT_FindIsotropic: eV is not an isotropic vector\n";
        throw TerminalException{1};
      }
      os_out << "return rec(has_isotropic:=true, ";
      os_out << "V:=" << StringVectorGAP(eV) << ");\n";
    } else {
      os_out << "return rec(has_isotropic:=false);\n";
    }
    return;
  }
  if (OutFormat == "PYTHON") {
    if (opt) {
      MyVector<T> const &eV = *opt;
      WriteVectorPYTHON(os_out, eV);
    } else {
      os_out << "None";
    }
    os_out << "\n";
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
      std::cerr << "LATT_FindIsotropic [arith] [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_FindIsotropic [arith] [FileI]\n";
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
    std::cerr << "Normal termination of LATT_FindIsotropic\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindIsotropic\n";
    exit(e.eVal);
  }
  runtime(time);
}
