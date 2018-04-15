#include "GroupFct.h"

int main(int argc, char *argv[])
{
  if (argc != 2) {
    fprintf(stderr, "Number of argument is = %d\n", argc);
    fprintf(stderr, "This program is used as\n");
    fprintf(stderr, "PrintAllGroupElements [TheFile]\n");
    return -1;
  }
  std::ifstream GRPfs;
  GRPfs.open(argv[1]);
  TheGroupFormat GRP=ReadGroup(GRPfs);
  GRPfs.close();
  //
  int nbPoint=GRP.n;
  std::cerr << "nbPoint=" << nbPoint << "\n";
  std::cerr << "group=" << *(GRP.group) << "\n";

  IteratorGrp eIter=GetInitialIterator(GRP);
  int nbElt=0;
  std::set<permlib::Permutation> ListPerm;
  while(1) {
    permlib::Permutation ePerm=GetPermutation(eIter);
    std::cerr << "ePerm=" << ePerm << "\n";
    ListPerm.insert(ePerm);
    nbElt++;
    int res=IteratorIncrease(eIter);
    if (res == -1)
      break;
  }
  std::cerr << "|ListPerm|=" << ListPerm.size() << "\n";
  std::cerr << "nbElt=" << nbElt << "\n";
}
