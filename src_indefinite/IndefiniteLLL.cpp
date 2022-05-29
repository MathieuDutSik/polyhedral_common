// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Indefinite_LLL.h"
#include "NumberTheory.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3 && argc != 2) {
      std::cerr << "IndefiniteLLL [FileI] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "IndefiniteLLL [FileI]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    using Tint = mpz_class;

    std::string FileI = argv[1];
    //
    MyMatrix<T> M = ReadMatrixFile<T>(argv[1]);
    std::cerr << "We have M\n";
    //
    auto print_result = [&](std::ostream &os) -> void {
      ResultIndefiniteLLL<T, Tint> ResLLL = Indefinite_LLL<T, Tint>(M);
      if (ResLLL.success) {
        std::cerr << "B_T=\n";
        WriteMatrix(std::cerr, ResLLL.B);
        std::cerr << "Mred=\n";
        WriteMatrix(std::cerr, ResLLL.Mred);

        os << "return rec(B:=";
        WriteMatrixGAP(os, ResLLL.B);
        os << ", Mred:=";
        WriteMatrixGAP(os, ResLLL.Mred);
        os << ");\n";
      } else {
        os << "return rec(Xisotrop:=";
        WriteVectorGAP(os, ResLLL.Xisotrop);
        os << ");\n";
      }
    };
    if (argc == 2) {
      print_result(std::cerr);
    } else {
      std::string FileO = argv[2];
      std::ofstream os(FileO);
      print_result(os);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
