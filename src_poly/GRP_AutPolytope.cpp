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
    using Tidx_value = int16_t;
    //
    std::cerr << "GRP_AutPolytope : Reading input\n";
    //
    std::ifstream is(argv[1]);
    MyMatrix<T> TheEXT = ReadMatrix<T>(is);
    int nbVert = TheEXT.rows();
    WeightMatrix<true, T, Tidx_value> WMat = GetWeightMatrix<T,Tidx_value>(TheEXT);
    std::cerr << "We have WMat\n";
    PrintWeightedMatrix(std::cerr, WMat);
    Tgr eGR=GetGraphFromWeightedMatrix<T,Tgr,Tidx_value>(WMat);
    std::cerr << "We have eGR\n";
    GRAPH_PrintOutput(std::cerr , eGR);
    std::vector<std::vector<Tidx>> ListGen = GetGroupCanonicalizationVector_Kernel<Tgr,Tidx>(eGR, nbVert).second;
    Tgroup GRP = GetGroupListGen<Tgroup>(ListGen, nbVert);
    std::cerr << "|GRP|=" << GRP.size() << "\n";

  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
