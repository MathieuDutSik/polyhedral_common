// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_LinPolytope_Automorphism [EXTIN]\n";
      std::cerr << "\n";
      std::cerr
          << "EXTIN : The list of vertices (or inequalities for that matter)\n";
      return -1;
    }
    //
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::string eFile = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(eFile);
    //
    Tgroup GRP =
        LinPolytope_Automorphism<T, Tgroup>(EXT, std::cerr);
    //
    for (auto & eGen : GRP.GeneratorsOfGroup()) {
      std::optional<MyMatrix<T>> opt = FindTransformationGeneral(EXT, EXT, eGen);
      if (!opt) {
        std::cerr << "Failed to find a transformation\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_LinPolytope_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
