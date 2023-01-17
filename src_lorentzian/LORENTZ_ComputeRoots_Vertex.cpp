// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "Permutation.h"

#include "edgewalk.h"
#include "vinberg.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 3 && argc != 4) {
      std::cerr << "Program is used as\n";
      std::cerr << "LORENTZ_ComputeRoots_Vertex [G] [Vertex]\n";
      std::cerr << "or\n";
      std::cerr << "LORENTZ_ComputeRoots_Vertex [G] [Vertex] [outfile]\n";
      throw TerminalException{1};
    }
    std::string FileLorMat = argv[1];
    std::string FileVertex = argv[2];
    using T = mpq_class;
    using Tint = mpz_class;
    //
    MyMatrix<T> G = ReadMatrixFile<T>(FileLorMat);
    //
    std::ifstream is2(FileVertex);
    MyVector<T> gen = ReadVector<T>(is2);
    //
    // This code is mainly for testing
    std::string OptionNorms = "all";
    std::vector<T> l_norms = get_initial_list_norms<T, Tint>(G, OptionNorms);
    //
    MyMatrix<Tint> MatRoot = get_simple_cone<T, Tint>(G, l_norms, gen);
    //
    auto prt = [&](std::ostream &os) -> void {
      os << "return " << StringMatrixGAP(MatRoot) << ";\n";
    };
    if (argc == 3) {
      prt(std::cerr);
    }
    if (argc == 4) {
      std::string file = argv[3];
      std::ofstream os(file);
      prt(os);
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LORENTZ_ComputeRoots_Vertex\n";
    exit(e.eVal);
  }
  runtime(time1);
}
