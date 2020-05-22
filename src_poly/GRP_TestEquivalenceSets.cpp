#include "GRP_GroupFct.h"
int main(int argc, char *argv[])
{
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
    std::cerr << "Reading input\n";
    //
    std::ifstream is1(argv[1]);
    TheGroupFormat GRP=ReadGroup(is1);
    //
    std::ifstream is2(argv[2]);
    std::vector<Face> ListFace=ReadListFace(is2);
    //
    EquivTest<permlib::Permutation>  eReply=TestEquivalenceGeneralNaked(GRP, ListFace[0], ListFace[1], 1);
    //
    std::cerr << "result=" << eReply.TheReply << "\n";
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
