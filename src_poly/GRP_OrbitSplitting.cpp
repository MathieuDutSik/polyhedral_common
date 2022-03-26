#include "GRP_DoubleCoset.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TestEquivalenceSets [BigGRP] [SmaGRP] [ListFace]\n";
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
    std::ifstream is1(argv[1]);
    Tgroup BigGRP = ReadGroup<Tgroup>(is1);
    //
    std::ifstream is2(argv[2]);
    Tgroup SmaGRP = ReadGroup<Tgroup>(is2);
    //
    std::ifstream is3(argv[3]);
    vectface ListFaceBig = ReadListFace(is3);
    std::cerr << "|ListFaceBig|=" << ListFaceBig.size() << "\n";
    //
    vectface ListFaceSma =
        OrbitSplittingListOrbit(BigGRP, SmaGRP, ListFaceBig, std::cerr);
    std::cerr << "|ListFaceSma|=" << ListFaceSma.size() << "\n";
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
