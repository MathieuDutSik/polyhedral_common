// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "PolytopeEquiStab.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tfield>
void process(std::string const &FileExt, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(FileExt);
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  int n_rows = EXT.rows();
  //
  MyMatrix<T> Qinv = GetQmatrix(EXTred, std::cerr);
  std::vector<MyMatrix<T>> ListMat = {Qinv};
  std::vector<T> Vdiag(n_rows, 0);
  //
  size_t seed = 15;
  size_t e_hash = GetInvariant_ListMat_Vdiag<T, Tfield>(seed, EXTred, ListMat,
                                                        Vdiag, std::cerr);
  //
  if (OutFormat == "GAP") {
    os << "return " << e_hash << ";\n";
    return;
  }
  if (OutFormat == "CPP") {
    os << e_hash << "\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string FileExt;
    std::string OutFormat = "CPP";
    std::string FileOut = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_arithmetics(dispatch, [&]<typename T>() -> void {
      using Tfield = T;
      process<T, Tfield>(FileExt, OutFormat, *os_ptr);
    });
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Invariant [Arith] [FileEXT] [OutFormat] "
                   "[FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Invariant [Arith] [FileEXT]\n";
      std::cerr << "\n";
      std::cerr << "FileEXT : The list of vertices\n";
      std::cerr << "FileOut : The file for the hash\n";
      std::cerr << "\n";
      std::cerr << "     ------- Arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    FileExt = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }

    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    std::cerr << "Normal termination of GRP_LinPolytope_Invariant\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Invariant\n";
    exit(e.eVal);
  }
  runtime(time);
}
