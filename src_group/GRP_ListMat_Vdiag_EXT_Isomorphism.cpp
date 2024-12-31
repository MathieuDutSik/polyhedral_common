// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "GRP_GroupFct.h"
#include "PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "GRP_ListMat_Vdiag_EXT_Isomorphism [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_ListMat_Vdiag_EXT_Isomorphism [FileI]\n";
      std::cerr << "\n";
      std::cerr << "INfile    : The file containing the group\n";
      std::cerr << "OUTfile   : The file containing the two pairs\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpz_int;
    using Tfield = boost::multiprecision::mpq_rational;
#else
    using T = mpz_class;
    using Tfield = mpq_class;
#endif
    using Tidx = unsigned int;
    //
    std::cerr << "Reading input\n";
    //
    std::string FileI = argv[1];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
    std::ifstream is(FileI);
    int nbMat;
    is >> nbMat;
    //
    std::vector<MyMatrix<T>> ListMat1;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<T> eMat = ReadMatrix<T>(is);
      ListMat1.push_back(eMat);
    }
    MyMatrix<T> EXT1 = ReadMatrix<T>(is);
    for (auto &eMat1 : ListMat1) {
      if (!IsSymmetricMatrix(eMat1)) {
        std::cerr << "The matrix should be symmetric\n";
        throw TerminalException{1};
      }
      if (eMat1.cols() != EXT1.cols()) {
        std::cerr << "|eMat1|=" << eMat1.cols() << " |EXT1|=" << EXT1.cols()
                  << "\n";
        throw TerminalException{1};
      }
    }
    int n_row1 = EXT1.rows();
    int len;
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
    //
    is >> nbMat;
    std::vector<MyMatrix<T>> ListMat2;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<T> eMat = ReadMatrix<T>(is);
      ListMat2.push_back(eMat);
    }
    MyMatrix<T> EXT2 = ReadMatrix<T>(is);
    for (auto &eMat2 : ListMat2) {
      if (!IsSymmetricMatrix(eMat2)) {
        std::cerr << "The matrix should be symmetric\n";
        throw TerminalException{1};
      }
      if (eMat2.cols() != EXT2.cols()) {
        std::cerr << "|eMat2|=" << eMat2.cols() << " |EXT2|=" << EXT2.cols()
                  << "\n";
        throw TerminalException{1};
      }
    }
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
    //
    std::cerr << "Input read, now testing for equivalence\n";
    std::optional<std::vector<Tidx>> PairTest =
        TestEquivalence_ListMat_Vdiag<T, Tfield, Tidx>(
            EXT1, ListMat1, Vdiag1, EXT2, ListMat2, Vdiag2, std::cerr);
    //
    auto prt = [&](std::ostream &os) -> void {
      if (OutFormat == "GAP") {
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
        return;
      }
      if (OutFormat == "Oscar") {
        if (!PairTest) {
          os << "0\n";
        } else {
          const std::vector<Tidx> &V = *PairTest;
          int n_rows = EXT1.rows();
          os << n_rows << "\n";
          for (int i = 0; i < n_rows; i++) {
            int pos = V[i] + 1;
            os << " " << pos;
          }
        }
        return;
      }
      std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
                << "\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      prt(std::cerr);
    } else {
      if (FileO == "stdout") {
        prt(std::cout);
      } else {
        std::ofstream os(FileO);
        prt(os);
      }
    }
    std::cerr << "Normal termination of GRP_ListMat_Vdiag_EXT_Isomorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_ListMat_Vdiag_EXT_Isomorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
