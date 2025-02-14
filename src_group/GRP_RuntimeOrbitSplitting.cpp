// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_DoubleCoset.h"
#include "GRP_GroupFct.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    std::vector<std::string> ListMethod = {"repr",
                                           "canonic",
                                           "canonic_initial_triv",
                                           "exhaustive_std",
                                           "exhaustive_sparse",
                                           "exhaustive_robin",
                                           "exhaustive_hopscotch",
                                           "single_cosets",
                                           "double_cosets"};
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_RuntimeOrbitSplitting [FileDoubleCoset] [method]\n";
      std::cerr << "\n";
      std::cerr << "with FileDoubleCoset containing the BigGRP, the SmaGRP and "
                   "the list of orbits\n";
      std::cerr << "method can be all in which case all methods are used\n";
      std::cerr << "methods =";
      for (auto &method : ListMethod)
        std::cerr << " " << method;
      std::cerr << "\n";
      return -1;
    }
    //
    using Tint = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::cerr << "Reading input\n";
    //
    std::string eFile = argv[1];
    std::ifstream is(eFile);
    Tgroup BigGRP = ReadGroup<Tgroup>(is);
    std::cerr << "|BigGRP|=" << BigGRP.size() << "\n";
    Tgroup SmaGRP = ReadGroup<Tgroup>(is);
    std::cerr << "|SmaGRP|=" << SmaGRP.size() << "\n";
    vectface ListFaceBig = ReadListFace(is);
    //
    FaceOrbitsizeGrpContainer ListFaceOrbitsizes(BigGRP,
                                                 std::move(ListFaceBig));
    //
    std::string str_methods = argv[2];
    if (str_methods != "all") {
      ListMethod = STRING_Split(str_methods, ",");
    }
    std::cerr << "ListMethod =";
    for (auto &method : ListMethod)
      std::cerr << " " << method;
    std::cerr << "\n";
    //
    std::set<size_t> SetSizes;
    for (auto &method : ListMethod) {
      std::cerr << "-----------------------------\n";
      HumanTime time;
      vectface ListFaceSma = OrbitSplittingListOrbit_spec(
          BigGRP, SmaGRP, ListFaceOrbitsizes, method, std::cerr);
      std::cerr << "Result for method=" << method
                << " |ListFaceSma|=" << ListFaceSma.size() << " time=" << time
                << "\n";
      SetSizes.insert(ListFaceSma.size());
    }
    std::cerr << "-----------------------------\n";
    if (SetSizes.size() != 1) {
      std::cerr << "Incoherency in the orbit splitting computation. Methods "
                   "give different results\n";
      throw TerminalException{1};
    }
    //
    std::cerr << "Normal termination of GRP_RuntimeOrbitSplitting\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_RuntimeOrbitSplitting\n";
    exit(e.eVal);
  }
  runtime(time1);
}
