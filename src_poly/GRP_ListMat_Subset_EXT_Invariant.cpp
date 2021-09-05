#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_Invariant_ListMat_Subset_EXT [INfile] [OUTfile]\n";
      std::cerr << "\n";
      std::cerr << "INfile    : The file containing the group\n";
      std::cerr << "OUTfile   : The file containing the two pairs\n";
      return -1;
    }
    using T = mpz_class;
    //
    std::cerr << "GRP_ComputeAut_ListMat_Subset_EXT : Reading input\n";
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
    std::cerr << "n_rows=" << n_rows << "\n";
    Face eSubset = ReadFace(is);
    std::cerr << "|eSubset|=" << eSubset.size() << " / " << eSubset.count() << "\n";
    //
    size_t e_hash = GetInvariant_ListMat_Subset(EXT, ListMat, eSubset);
    //
    std::ofstream os(argv[2]);
    os << "return " << e_hash << ";\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
