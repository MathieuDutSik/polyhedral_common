// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Indefinite_LLL.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4 && argc != 2) {
      std::cerr
          << "LATT_SublatticeBasisReduction [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_SublatticeBasisReduction [FileI]\n";
      std::cerr << "or\n";
      std::cerr << "FileI     : The input file\n";
      std::cerr << "OutFormat : Possible values, GAP and Oscar\n";
      std::cerr << "            Default: GAP\n";
      std::cerr << "FileO     : The output file\n";
      std::cerr << "            stdout: write to std::cout\n";
      std::cerr << "            stderr: write to std::cerr\n";
      std::cerr << "            Other filename get written to files\n";
      std::cerr << "            Default: stderr\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    std::string FileI = argv[1];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
    //
    MyMatrix<T> M = ReadMatrixFile<T>(argv[1]);
    std::cerr << "We have M\n";
    //
    auto print_result = [&](std::ostream &os) -> void {
      MyMatrix<T> Mred = SublatticeBasisReduction(M);
      if (OutFormat == "GAP") {
        os << "return rec(Mred:=";
        WriteMatrixGAP(os, Mred);
        os << ");\n";
        return;
      }
      if (OutFormat == "Oscar") {
        WriteMatrix(os, Mred);
        return;
      }
      std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
                << "\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileO, print_result);
    std::cerr << "Normal termination of LATT_SublatticeBasisReduction\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_SublatticeBasisReduction\n";
    exit(e.eVal);
  }
  runtime(time);
}
