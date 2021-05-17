#include "Permutation.h"
#include "Group.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_AutPolytope [INfile]\n";
      std::cerr << "\n";
      std::cerr << "INfile    : The file containing the group\n";
      return -1;
    }
    using T = mpz_class;
    using Tidx = unsigned int;
    using Tgr = GraphListAdj;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt,Tint>;
    //
    std::cerr << "GRP_AutPolytope : Reading input\n";
    //
    std::ifstream is(argv[1]);
    MyMatrix<T> TheEXT = ReadMatrix<T>(is);
    int nbVert = TheEXT.rows();
    WeightMatrix<true, T> WMat = GetWeightMatrix(TheEXT);
    Tgr eGR=GetGraphFromWeightedMatrix<T,Tgr>(WMat);
    std::vector<std::vector<Tidx>> ListGen = GetGroupCanonicalizationVector_Kernel<Tgr,Tidx>(eGR, nbVert).second;
    Tgroup GRP = GetGroupListGen<Tgroup>(ListGen, nbVert);
    std::cerr << "|GRP|=" << GRP.size() << "\n";

  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
