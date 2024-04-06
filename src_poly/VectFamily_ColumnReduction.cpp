// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "POLY_Fundamental.h"
#include "LatticeDefinitions.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using T = mpq_class;
    if (argc != 2 && argc != 4) {
      std::cerr << "VectFamily_ColumnReduction [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "VectFamily_ColumnReduction [FileI]\n";
      std::cerr << "\n";
      std::cerr << " ------- OutFormat --------\n";
      std::cerr << "Possible values for OutFormat\n";
      std::cerr << "GAP : for writing in GAP readable file\n";
      std::cerr << "CPP : for writing in CPP polyhedral readable file\n";
      std::cerr << "\n";
      std::cerr << " ------- FileO --------\n";
      std::cerr << "Possible values for FileO\n";
      std::cerr << "stderr : for writing to std::cerr\n";
      std::cerr << "stdout : for writing in std::cout\n";
      std::cerr << "otherwise written to the named file in output\n";
      throw TerminalException{1};
    }
    std::string FileInput = argv[1];
    MyMatrix<T> M = ReadMatrixFile<T>(FileInput);
    std::string OutFormat = "CPP";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
    MyMatrix<T> Mred = ColumnReduction(M);

    auto print_mat = [&](std::ostream &os) -> void {
      if (OutFormat == "GAP") {
        os << "return ";
        WriteMatrixGAP(os, Mred);
        os << ";\n";
        return;
      }
      if (OutFormat == "CPP") {
        return WriteMatrix(os, Mred);
      }
      std::cerr << "No matching format in print_mat\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      print_mat(std::cerr);
    } else {
      if (FileO == "stdout") {
        print_mat(std::cout);
      } else {
        std::ofstream os(FileO);
        print_mat(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in VectFamily_Reduction\n";
    exit(e.eVal);
  }
  runtime(time1);
}
