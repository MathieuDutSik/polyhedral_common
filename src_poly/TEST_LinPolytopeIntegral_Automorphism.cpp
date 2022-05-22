#include "GRP_GroupFct.h"
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_LinPolytope_Automorphism [EXTIN]\n";
      std::cerr << "\n";
      std::cerr
          << "EXTIN : The list of vertices (or inequalities for that matter)\n";
      return -1;
    }
    //
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint32_t;
    using Tidx_value = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::ifstream is(argv[1]);
    MyMatrix<T> EXT = ReadMatrix<T>(is);
    MyMatrix<T> preWMat_EXT = ReadMatrix<int>(is);
    int nbRow = WMat_EXT.rows();
    int nbCol = WMat_EXT.cols();
    MyMatrix<std::vector<T>> WMat_EXT(nbRow,nbCol);
    for (int iRow=0; iRow<nbRow; iRow++) {
      for (int iCol=0; iCol<nbCol; iCol++) {
        std::vector<T> val{preWMat_EXT(iRow,iCol)};
        WMat_EXT(iRow,iCol) = val;
      }
    }
    //
    using pair_char =
        std::pair<MyMatrix<T>, WeightMatrix<true, std::vector<T>, Tidx_value>>;
    //
    pair_char ep{
        EXT,
        WeightedMatrixFromMyMatrix<true, std::vector<T>, Tidx_value>(WMat_EXT)};
    std::vector<MyMatrix<T>> LGen =
        LinPolytopeIntegralWMat_Automorphism<T, Tgroup, std::vector<T>,
                                             Tidx_value>(ep);
    std::cerr << "We have LGen\n";
    std::cerr << "|LGen|=" << LGen.size() << "\n";
    //    Tgroup GRP =
    //    LinPolytopeIntegralWMat_Automorphism<T,Tgroup,std::vector<T>,Tidx_value>(ep))
    //    std::cerr << GRP << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
