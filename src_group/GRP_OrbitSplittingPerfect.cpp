// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_DoubleCoset.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 6) {
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_OrbitSplittingPerfect [BigGRP] [SmaGRP] [ListFace] "
                   "[FileOut2] [FileOut3]\n";
      std::cerr << "\n";
      std::cerr << "BigGRP     : The file containing the big group\n";
      std::cerr << "SmaGRP     : The file containing the small group\n";
      std::cerr << "ListFace   : The file containing the list of faces\n";
      std::cerr << "FileOut2   : The second file of the output\n";
      std::cerr << "FileOut3   : The third file of the output\n";
      return -1;
    }
    //
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::cerr << "Reading input\n";
    //
    std::string FileGroupBig = argv[1];
    Tgroup BigGRP = ReadGroupFile<Tgroup>(FileGroupBig);
    //
    std::string FileGroupSmall = argv[2];
    Tgroup SmaGRP = ReadGroupFile<Tgroup>(FileGroupSmall);
    //
    std::ifstream is3(argv[3]);
    vectface ListFaceBig = ReadListFace(is3);
    std::cerr << "|ListFaceBig|=" << ListFaceBig.size() << "\n";
    //
    // Now the output
    //
    std::ofstream os2(argv[4]);
    std::ofstream os3(argv[5]);
    OrbitSplittingPerfectFacet(BigGRP, SmaGRP, ListFaceBig, os2, os3,
                               std::cerr);
    //
    std::cerr << "Normal termination of GRP_OrbitSplittingPerfect\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_OrbitSplittingPerfect\n";
    exit(e.eVal);
  }
  runtime(time1);
}
