#include "Permutation.h"
#include "Group.h"
#include "GRP_GroupFct.h"


int main(int argc, char *argv[])
{
  try {
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TestEquivalenceSets [BigGRP] [SmaGRP] [ListFace] [FileOut2] [FileOut3]\n";
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
    using Tgroup = permutalib::Group<Telt,Tint>;
    std::cerr << "Reading input\n";
    //
    std::ifstream is1(argv[1]);
    Tgroup BigGRP=ReadGroup<Tgroup>(is1);
    //
    std::ifstream is2(argv[2]);
    Tgroup SmaGRP=ReadGroup<Tgroup>(is2);
    //
    std::ifstream is3(argv[3]);
    vectface ListFaceBig=ReadListFace(is3);
    std::cerr << "|ListFaceBig|=" << ListFaceBig.size() << "\n";
    //
    // Now the output
    //
    std::ofstream os2(argv[4]);
    std::ofstream os3(argv[5]);
    OrbitSplittingPerfectFacet(BigGRP, SmaGRP, ListFaceBig, os2, os3);
    //
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
