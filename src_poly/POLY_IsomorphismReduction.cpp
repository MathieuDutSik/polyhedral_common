// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheory.h"
#include "POLY_lrslib.h"

int main(int argc, char *argv[]) {
  try {
    using T = mpq_class;
    //
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_IsomorphismReduction [EXT] [FAC] [GRP] [OUT]\n";
      std::cerr << "\n";
      std::cerr << "EXT : The vertices\n";
      std::cerr << "FAC : The facets\n";
      std::cerr << "GRP : The group\n";
      std::cerr << "OUT : The gap readable file of output\n";
      return -1;
    }
    //
    std::string FileEXT = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
    //
    std::string FileFAC = argv[2];
    MyMatrix<T> FAC = ReadMatrixFile<T>(FileFAC);
    //
    std::string FileGRP = argv[3];
    Tgroup GRP = ReadGRoupFile<Tgroup>(FileGRP);
    //
    std::string FileOUT = argv[4];
    //
    std::ifstream is(argv[1]);
    MyMatrix<T> EXT = ReadMatrixLrsCdd<T>(is);
    //
    vectface ListFace = lrs::DualDescription_temp_incd_reduction(EXT);
    std::cerr << "nbVert = " << ListFace.size() << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
