// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_SamplingFacet.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << " [DATAIN] [DATAOUT]\n";
      return -1;
    }
    //
    std::ifstream is(argv[1]);
    using T = mpq_class;
    MyMatrix<T> EXT = ReadMatrix<T>(is);
    size_t n_rows = EXT.rows();
    //
    vectface ListFace = GetFullRankFacetSet(EXT);
    std::cerr << "We have ListFace\n";
    //
    std::ofstream os(argv[2]);
    os << "return [";
    bool IsFirst = true;
    for (auto &eFace : ListFace) {
      if (!IsFirst)
        os << ",\n";
      bool IsFirstB = true;
      os << "[";
      for (size_t i_row = 0; i_row < n_rows; i_row++) {
        if (eFace[i_row] == 1) {
          if (!IsFirstB)
            os << ",";
          os << (i_row + 1);
          IsFirstB = false;
        }
      }
      os << "]";
      IsFirst = false;
    }
    os << "];\n";
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
