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
    const bool use_scheme1 = true;
    Tgroup GRP1 =
        LinPolytope_Automorphism<T, use_scheme1, Tgroup>(EXT, std::cerr);
    //
    const bool use_scheme2 = false;
    Tgroup GRP2 =
        LinPolytope_Automorphism<T, use_scheme2, Tgroup>(EXT, std::cerr);
    //
    bool test = GRP1 == GRP2;
    if (!test) {
      std::cerr << "test = False\n";
      std::cerr << "The groups are different. It is the clear bug\n";
      throw TerminalException{1};
    } else {
      std::cerr << "test = True\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_LinPolytope_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
