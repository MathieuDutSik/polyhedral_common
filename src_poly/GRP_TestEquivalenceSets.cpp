// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GRP_GroupFct.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TestEquivalenceSets [GRPfile] [PairFile]\n";
      std::cerr << "\n";
      std::cerr << "GRPfile    : The file containing the group\n";
      std::cerr << "PairFile   : The file containing the two pairs\n";
      return -1;
    }
    //
    using Tint = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::cerr << "Reading input\n";
    //
    std::ifstream is1(argv[1]);
    Tgroup GRP = ReadGroup<Tgroup>(is1);
    //
    std::ifstream is2(argv[2]);
    vectface ListFace = ReadListFace(is2);
    //
    std::optional<Telt> eReply =
        GRP.RepresentativeAction_OnSets(ListFace[0], ListFace[1]);
    //
    if (eReply)
      std::cerr << "result=true\n";
    else
      std::cerr << "result=false\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
