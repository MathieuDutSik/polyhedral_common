#include "Permlib_specific.h"
#include "GRP_GroupFct.h"
int main(int argc, char *argv[])
{
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
    using Tgroup=TheGroupFormat<mpz_class>;
    std::cerr << "Reading input\n";
    //
    std::ifstream is1(argv[1]);
    Tgroup BigGRP=ReadGroup<Tgroup>(is1);
    //
    std::ifstream is2(argv[2]);
    Tgroup SmaGRP=ReadGroup<Tgroup>(is2);
    //
    std::ifstream is3(argv[3]);
    std::vector<Face> ListFaceBig=ReadListFace(is3);
    std::cerr << "|ListFaceBig|=" << ListFaceBig.size() << "\n";
    //
    std::vector<Face> ListFaceSma=OrbitSplittingListOrbit(BigGRP, SmaGRP, ListFaceBig, std::cerr);
    std::cerr << "|ListFaceSma|=" << ListFaceSma.size() << "\n";
    //
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
