#include "NumberTheory.h"
#include "Permutation.h"
#include "Group.h"
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
    using Tgr = GraphBitset;
    using Tidx = int;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt,Tint>;
    using Twmat = WeightMatrix<true, Tint, Tidx_value>;
    //
    std::ifstream is(argv[1]);
    MyMatrix<Tint> EXT=ReadMatrix<Tint>(is);
    int nbCol=EXT.cols();
    int nbRow=EXT.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    auto get_canonicalized_wmat=[](MyMatrix<Tint> const& EXT) -> std::pair<Twmat,Tgroup> {
      size_t n_row = EXT.rows();
      Twmat WMat = GetWeightMatrix<Tint,Tidx_value>(EXT);
      WMat.ReorderingSetWeight();
      std::pair<std::vector<Tidx>,std::vector<std::vector<Tidx>>> epair = GetGroupCanonicalizationVector_Kernel<Tint,Tgr,Tidx,Tidx_value>(WMat);
      const std::vector<Tidx>& ListIdx = epair.first;
      const std::vector<std::vector<Tidx>>& ListGen = epair.second;
      WMat.RowColumnReordering(ListIdx);
      std::vector<Telt> LGen;
      std::vector<Tidx> ListIdxRev(n_row);
      for (size_t i=0; i<n_row; i++)
        ListIdxRev[ListIdx[i]] = i;
      for (auto & eGen : ListGen) {
        std::vector<Tidx> V(n_row);
        for (size_t i1=0; i1<n_row; i1++) {
          size_t i2 = ListIdx[i1];
          size_t i3 = eGen[i2];
          size_t i4 = ListIdxRev[i3];
          V[i1] = i4;
        }
        Telt ePerm(V);
        LGen.push_back(ePerm);
      }
      Tgroup GRP(LGen, n_row);
      return {std::move(WMat), std::move(GRP)};
    };
    std::pair<Twmat,Tgroup> Pair1 = get_canonicalized_wmat(EXT);
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
      std::pair<Twmat,Tgroup> Pair2 = get_canonicalized_wmat(EXT2);
      std::cerr << "------------------------------------------------------------\n";
      if (Pair1.first != Pair2.first) {
        std::cerr << "The reordering of the column of the matrix failed\n";
        throw TerminalException{1};
      }
      if (Pair1.second != Pair2.second) {
        std::cerr << "The groups are not coherent\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
