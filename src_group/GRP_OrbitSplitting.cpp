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
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_OrbitSplitting [BigGRP] [SmaGRP] [ListFace]\n";
      std::cerr << "\n";
      std::cerr << "BigGRP     : The file containing the big group\n";
      std::cerr << "SmaGRP     : The file containing the small group\n";
      std::cerr << "ListFace   : The file containing the list of faces\n";
      return -1;
    }
    //
    using Tint = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::cerr << "Reading input\n";
    //
    std::string BigGRP_file = argv[1];
    Tgroup BigGRP = ReadGroupFile<Tgroup>(BigGRP_file);
    //
    std::string SmaGRP_file = argv[2];
    Tgroup SmaGRP = ReadGroupFile<Tgroup>(SmaGRP_file);
    //
    std::ifstream is3(argv[3]);
    vectface ListFaceBig = ReadListFace(is3);
    std::cerr << "|ListFaceBig|=" << ListFaceBig.size() << "\n";
    //
    FaceOrbitsizeGrpContainer ListFaceOrbitsizes(BigGRP,
                                                 std::move(ListFaceBig));
    vectface ListFaceSma =
        OrbitSplittingListOrbit(BigGRP, SmaGRP, ListFaceOrbitsizes, std::cerr);
    std::cerr << "|ListFaceSma|=" << ListFaceSma.size() << "\n";
    //
    std::cerr << "Normal termination of GRP_OrbitSplitting\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_OrbitSplitting\n";
    exit(e.eVal);
  }
  runtime(time1);
}
