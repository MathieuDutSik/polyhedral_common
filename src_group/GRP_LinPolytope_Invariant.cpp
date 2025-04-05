// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "PolytopeEquiStab.h"
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
  size_t e_hash =
      GetInvariant_ListMat_Vdiag<T, Tfield>(EXTred, ListMat, Vdiag, std::cerr);
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
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Invariant [Arith] [FileEXT] [OutFormat] "
                   "[FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Invariant [Arith] [FileEXT]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices\n";
      std::cerr << "OutCan : The file for the hash\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileExt = argv[2];
    std::string OutFormat = "CPP";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }

    auto f = [&](std::ostream &os) -> void {
      if (arith == "rational") {
        using T = mpq_class;
        using Tfield = T;
        return process<T, Tfield>(FileExt, OutFormat, os);
      }
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        using Tfield = T;
        return process<T, Tfield>(FileExt, OutFormat, os);
      }
      std::cerr << "Failed to find a matching arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileOut, f);
    std::cerr << "Normal termination of GRP_LinPolytope_Invariant\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Invariant\n";
    exit(e.eVal);
  }
  runtime(time);
}
