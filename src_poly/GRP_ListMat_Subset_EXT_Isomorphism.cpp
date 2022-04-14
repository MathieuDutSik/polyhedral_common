#include "GRP_GroupFct.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism [INfile]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism [INfile] [OUTfile]\n";
      std::cerr << "\n";
      std::cerr << "INfile    : The file containing the group\n";
      std::cerr << "OUTfile   : The file containing the two pairs\n";
      return -1;
    }
    using T = mpz_class;
    using Tidx = unsigned int;
    //
    std::cerr << "Reading input\n";
    //
    std::ifstream is(argv[1]);
    int nbMat, len;
    is >> nbMat;
    std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism, step 1 nbMat=" << nbMat << "\n";
    //
    std::vector<MyMatrix<T>> ListMat1;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<T> eMat = ReadMatrix<T>(is);
      ListMat1.push_back(eMat);
    }
    std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism, step 2 |ListMat1|=" << ListMat1.size() << "\n";
    MyMatrix<T> EXT1 = ReadMatrix<T>(is);
    std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism, step 3 |EXT1|=" << EXT1.rows() << " / " << EXT1.cols() << "\n";
    int n_row1 = EXT1.rows();
    is >> len;
    if (len != n_row1) {
      std::cerr << "We have n_row1=" << n_row1 << " but len=" << len << "\n";
      throw TerminalException{1};
    }
    std::vector<T> Vdiag1(n_row1);
    for (int i = 0; i < n_row1; i++) {
      T val;
      is >> val;
      Vdiag1[i] = val;
    }
    std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism, step 4 |Vdiag1|=" << Vdiag1.size() << "\n";
    //
    is >> nbMat;
    std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism, step 4.1 nbMat=" << nbMat << "\n";
    std::vector<MyMatrix<T>> ListMat2;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<T> eMat = ReadMatrix<T>(is);
      std::cerr << "iMat=" << iMat << " |eMat|=" << eMat.rows() << " / " << eMat.cols() << "\n";
      ListMat2.push_back(eMat);
    }
    std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism, step 5 |ListMat2|=" << ListMat2.size() << "\n";
    MyMatrix<T> EXT2 = ReadMatrix<T>(is);
    std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism, step 6 |EXT2|=" << EXT2.rows() << " / " << EXT2.cols() << "\n";
    int n_row2 = EXT2.rows();
    is >> len;
    if (len != n_row2) {
      std::cerr << "We have n_row2=" << n_row2 << " but len=" << len << "\n";
      throw TerminalException{1};
    }
    std::vector<T> Vdiag2(n_row2);
    for (int i = 0; i < n_row2; i++) {
      T val;
      is >> val;
      Vdiag2[i] = val;
    }
    std::cerr << "GRP_ListMat_Subset_EXT_Isomorphism, step 7 |Vdiag2|=" << Vdiag2.size() << "\n";
    //
    std::cerr << "Input read, now testing for equivalence\n";
    const bool use_scheme = true;
    std::optional<std::vector<Tidx>> PairTest =
        TestEquivalence_ListMat_Vdiag<T, Tidx, use_scheme>(
            EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2);
    //
    auto prt=[&](std::ostream & os) -> void {
      if (!PairTest) {
        os << "return false;\n";
      } else {
        const std::vector<Tidx> &V = *PairTest;
        int n_rows = EXT1.rows();
        os << "return [";
        for (int i = 0; i < n_rows; i++) {
          if (i > 0)
            os << ",";
          os << (V[i] + 1);
        }
        os << "];\n";
      }
    };
    if (argc == 2) {
      prt(std::cerr);
    } else {
      std::ofstream os(argv[2]);
      prt(os);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
