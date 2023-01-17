// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GRP_GroupFct.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 3 && argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinPolytope_Automorphism [EXTIN] [OutGroup]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_LinPolytope_Automorphism [EXTIN]\n";
      std::cerr << "\n";
      std::cerr
          << "EXTIN : The list of vertices (or inequalities for that matter)\n";
      std::cerr << "OutGroup : The automorphism group file\n";
      return -1;
    }
    //
    using T = mpq_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::ifstream is(argv[1]);
    MyMatrix<T> EXT = ReadMatrix<T>(is);
    int nbCol = EXT.cols();
    int nbRow = EXT.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    const bool use_scheme = true;
    Tgroup GRP = LinPolytope_Automorphism<T, use_scheme, Tgroup>(EXT);
    std::cerr << "|GRP|=" << GRP.size() << "\n";
    //
    if (argc == 3) {
      std::ofstream os(argv[2]);
      WriteGroup(os, GRP);
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
