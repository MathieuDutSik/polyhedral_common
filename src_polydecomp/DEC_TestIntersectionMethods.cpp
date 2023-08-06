// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "Decompositions.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "DEC_TestIntersectionMethods [FileI]\n";
      std::cerr << "\n";
      std::cerr << "  ------ FileI -------\n";
      std::cerr << "\n";
      std::cerr << "FileI is the input file containing the cones\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    //
    std::string FileI = argv[1];
    std::ifstream is(FileI);
    //
    // The polyhedral cones.
    std::vector<MyMatrix<T>> l_FAC;
    //
    auto f_test=[](MyMatrix<T> const& FACtot) -> void {
      bool test1 = IsFullDimensional(FACtot);
      bool test2 = IsFullDimensionalNextGen(FACtot);
      if (test1 != test2) {
        std::cerr << "We have test1=" << test1 << " test2=" << test2 << "\n";
        std::cerr << "We have a bug to resolve\n";
        throw TerminalException{1};
      }
    };
    //
    size_t n_domain;
    is >> n_domain;
    for (size_t i = 0; i < n_domain; i++) {
      std::cerr << "i=" << i << " / " << n_domain << "\n";
      MyMatrix<T> FAC = ReadMatrix<T>(is);
      std::cerr << "We have read FAC, |FAC|=" << FAC.rows() << "/" << FAC.cols()
                << "\n";
      l_FAC.push_back(FAC);
    }
    //
    for (size_t i=0; i<n_domain; i++) {
      for (size_t j=i+1; j<n_domain; j++) {
        MyMatrix<T> FACtot = Concatenate(l_FAC[i], l_FAC[j]);
        f_test(FACtot);
      }
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
