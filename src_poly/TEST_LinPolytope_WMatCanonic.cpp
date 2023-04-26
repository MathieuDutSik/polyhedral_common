// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Group.h"
#include "NumberTheory.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_LinPolytope_Canonic [EXTIN]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices (or inequalities for that "
                   "matter)\n";
      return -1;
    }
    //
    using T = mpq_class;
    /*
      We need to have T a field since otherwise, the rescaling will rescale the
      weights and so the comparison will fail because of that.
     */
    static_assert(is_ring_field<T>::value, "Requires T to be a field");
    using Tidx_value = uint16_t;
    using Tgr = GraphBitset;
    using Tidx = int;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    using Twmat = WeightMatrix<true, T, Tidx_value>;
    //
    std::string eFile = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(eFile);
    int nbCol = EXT.cols();
    int nbRow = EXT.rows();
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    //
    auto get_canonicalized_wmat =
        [](MyMatrix<T> const &EXT) -> std::pair<Twmat, Tgroup> {
      size_t n_row = EXT.rows();
      Twmat WMat = GetWeightMatrix<T, Tidx_value>(EXT);
      WMat.ReorderingSetWeight();
      std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> epair =
          GetGroupCanonicalizationVector_Kernel<T, Tgr, Tidx, Tidx_value>(WMat);
      const std::vector<Tidx> &ListIdx = epair.first;
      const std::vector<std::vector<Tidx>> &ListGen = epair.second;
      WMat.RowColumnReordering(ListIdx);
      std::vector<Telt> LGen;
      std::vector<Tidx> ListIdxRev(n_row);
      for (size_t i = 0; i < n_row; i++)
        ListIdxRev[ListIdx[i]] = i;
      for (auto &eGen : ListGen) {
        std::vector<Tidx> V(n_row);
        for (size_t i1 = 0; i1 < n_row; i1++) {
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
    std::pair<Twmat, Tgroup> Pair1 = get_canonicalized_wmat(EXT);
    std::cerr
        << "------------------------------------------------------------\n";
    //
    auto get_random_equivalent = [](MyMatrix<T> const &eMat) -> MyMatrix<T> {
      int nbRow = eMat.rows();
      int n = eMat.cols();
      std::vector<int> ePerm = RandomPermutation<int>(nbRow);
      MyMatrix<T> eUnimod = RandomUnimodularMatrix<T>(n);
      MyMatrix<T> eProd = eMat * eUnimod;
      MyMatrix<T> RetMat(nbRow, n);
      for (int iRow = 0; iRow < nbRow; iRow++) {
        int jRow = ePerm[iRow];
        RetMat.row(iRow) = eProd.row(jRow);
      }
      return RetMat;
    };
    //
    int n_iter = 10;
    for (int i_iter = 0; i_iter < n_iter; i_iter++) {
      std::cerr << "i_iter=" << i_iter << " / " << n_iter << "\n";
      MyMatrix<T> EXT2 = get_random_equivalent(EXT);
      std::pair<Twmat, Tgroup> Pair2 = get_canonicalized_wmat(EXT2);
      std::cerr
          << "------------------------------------------------------------\n";
      if (Pair1.first != Pair2.first) {
        std::cerr << "The reordering of the column of the matrix failed\n";
        std::cerr << "Pair1.first=\n";
        PrintWeightedMatrix(std::cerr, Pair1.first);
        std::cerr << "Pair2.first=\n";
        PrintWeightedMatrix(std::cerr, Pair2.first);
        throw TerminalException{1};
      }
      if (Pair1.second != Pair2.second) {
        std::cerr << "The groups are not coherent\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_LinPolytope_WMatCanonic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
