#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT [INfile] [OUTfile]\n";
      std::cerr << "\n";
      std::cerr << "INfile    : The file containing the group\n";
      std::cerr << "OUTfile   : The file containing the two pairs\n";
      return -1;
    }
    using T = mpz_class;
    //
    std::cerr << "Reading input\n";
    //
    std::ifstream is(argv[1]);
    int nbMat;
    is >> nbMat;
    std::vector<MyMatrix<T>> ListMat;
    for (int iMat=0; iMat<nbMat; iMat++) {
      MyMatrix<T> eMat = ReadMatrix<T>(is);
      ListMat.push_back(eMat);
    }
    MyMatrix<T> EXT = ReadMatrix<T>(is);
    int n_rows= EXT.rows();
    Face eSubset = ReadFace(is);
    //
    std::vector<std::vector<unsigned int>> ListGen = GetListGenAutomorphism_ListMat_Subset(EXT, ListMat, eSubset);
    //
    std::ofstream os(argv[2]);
    os << "return [";
    bool IsFirst=true;
    for (auto & eList : ListGen) {
      if (!IsFirst)
        os << ",\n";
      IsFirst=false;
      os << "[";
      for (int i=0; i<n_rows; i++) {
        if (i > 0)
          os << ",";
        os << (eList[i] + 1);
      }
      os << "]";
    }
    os << "];\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
