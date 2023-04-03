// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GRP_DoubleCoset.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_RuntimeOrbitSplitting [FileDoubleCoset]\n";
      std::cerr << "   or\n";
      std::cerr << "GRP_RuntimeOrbitSplitting [FileDoubleCoset] [method]\n";
      std::cerr << "\n";
      std::cerr << "with FileDoubleCoset containing the BigGRP, the SmaGRP and the list of orbits\n";
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
    std::vector<std::string> ListMethod={"repr", "canonic", "exhaustive", "single_cosets"};
    if (argc == 3) {
      std::string method = argv[2];
      ListMethod = {method};
    }
    //
    for (auto & method : ListMethod) {
      HumanTime time;
      vectface ListFaceSma =
        OrbitSplittingListOrbit_spec(BigGRP, SmaGRP, ListFaceBig, method, std::cerr);
      std::cerr << "Result for method=" << method << " |ListFaceSma|=" << ListFaceSma.size() << " time=" << time << "\n";
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_OrbitSplitting\n";
    exit(e.eVal);
  }
  runtime(time1);
}
