#include "NumberTheory.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_LinPolytope_Canonic [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that matter)\n";
      return -1;
    }
    //
    using Tint = mpz_class;
    using Tidx_value = uint16_t;
    using Tgr=GraphBitset;
    using Tidx = int;
    const bool use_scheme = true;
    //
    std::ifstream is(argv[1]);
    MyMatrix<Tint> EXT=ReadMatrix<Tint>(is);
    int nbCol=EXT.cols();
    int nbRow=EXT.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    auto get_canonicalized_wmat=[&](MyMatrix<Tint> const& EXT) -> WeightMatrix<true, Tint, Tidx_value> {
      WeightMatrix<true, T, Tidx_value> WMat = GetWeightMatrix(EXT);
      WMat.ReorderingSetWeight();
      std::vector<Tidx> ListIdx = GetGroupCanonicalizationVector_Kernel<T,Tgr,Tidx,Tidx_value>(WMat).first;
      WMat.RowColumnReordering(ListIdx);
      return WMat;
    };
    WeightMatrix<true, Tint, Tidx_value> WMat1 = get_canonicalized_wmat(EXT);
    std::cerr << "------------------------------------------------------------\n";
    //
    auto get_random_equivalent=[](MyMatrix<Tint> const& eMat) -> MyMatrix<Tint> {
      int nbRow = eMat.rows();
      int n = eMat.cols();
      std::vector<int> ePerm = RandomPermutation<int>(nbRow);
      MyMatrix<Tint> eUnimod = RandomUnimodularMatrix<Tint>(n);
      MyMatrix<Tint> eProd = eMat * eUnimod;
      MyMatrix<Tint> RetMat(nbRow, n);
      for (int iRow=0; iRow<nbRow; iRow++) {
        int jRow = ePerm[iRow];
        RetMat.row(iRow) = eProd.row(jRow);
      }
      return RetMat;
    };
    //
    int n_iter = 10;
    for (int i_iter=0; i_iter<n_iter; i_iter++) {
      std::cerr << "i_iter=" << i_iter << " / " << n_iter << "\n";
      MyMatrix<Tint> EXT2 = get_random_equivalent(EXT);
      WeightMatrix<true, Tint, Tidx_value> WMat2 = get_canonicalized_wmat(EXT2);
      std::cerr << "------------------------------------------------------------\n";
      if (WMat1 != WMat1) {
        std::cerr << "The reordering of the column matrix failed\n";
        std::cerr << "EXT_Can=\n";
        WriteMatrix(std::cerr, EXT_can);
        std::cerr << "EXT2_Can=\n";
        WriteMatrix(std::cerr, EXT2_can);
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
