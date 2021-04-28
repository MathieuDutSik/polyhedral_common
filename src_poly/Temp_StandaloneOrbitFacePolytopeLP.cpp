#include "Permlib_specific.h"
#include "Group.h"
#include "Permutation.h"
#include "NumberTheory.h"
#include "POLY_Kskeletton.h"

int main(int argc, char *argv[])
{
  try {
    int LevSearch;
    int nbLev;
    if (argc != 5)
      {
	std::cerr << "Number of argument is = " << argc << "\n";
	std::cerr << "This program is used as\n";
	std::cerr << "StandaloneGroupPolytope [DATAEXT] [DATAGRP] [LevSearch] [DATAOUT]\n";
	std::cerr << "DATAEXT: The input data of the polytope\n";
	std::cerr << "DATAGRP: The input data of the group\n";
	std::cerr << "LevSearch: The level of the search\n";
	std::cerr << "DATAOUT: The output of the enumeration\n";
	std::cerr << "n\n";
	return -1;
      }
    std::cerr << "Reading input\n";
    std::ifstream EXTfs(argv[1]);
    using T=mpq_class;
    MyMatrix<T> TheEXT=ReadMatrix<T>(EXTfs);

    std::ifstream GRPfs(argv[2]);
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt,Tint>;
    //    using Tgroup=TheGroupFormat<mpz_class>;
    Tgroup TheGRP=ReadGroup<Tgroup>(GRPfs);

    LevSearch=atoi(argv[3]);
    std::cerr << "LevSearch=" << LevSearch << "\n";

    std::cerr << "Step main 5\n";
    std::vector<vectface> ListListOrb=EnumerationFaces<T,Tgroup>(TheGRP, TheEXT, LevSearch);
    std::cerr << "Completion of the program\n";

    std::ofstream OUTfs(argv[4]);
    PrintListListOrb_IntGAP(OUTfs, ListListOrb);

    nbLev=ListListOrb.size();
    std::cerr << "nbLev=" << nbLev << "\n";
    for (int iLev=0; iLev<nbLev; iLev++) {
      int nbOrb=ListListOrb[iLev].size();
      std::cerr << "  iLev=" << iLev << " nbOrb=" << nbOrb << "\n";
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
