// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "MatrixCanonicalForm.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Isomorphism [ListMat1] [ListMat2] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Isomorphism [ListMat1] [ListMat2]\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
#else
    using T = mpq_class;
    using Tint = mpz_class;
#endif
    using Tidx = uint32_t;
    //
    std::string FileListMat1 = argv[1];
    std::string FileListMat2 = argv[2];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    std::vector<MyMatrix<T>> ListMat1 = ReadListMatrixFile<T>(FileListMat1);
    std::vector<MyMatrix<T>> ListMat2 = ReadListMatrixFile<T>(FileListMat2);

    MyMatrix<Tint> SHV1 = ExtractInvariantVectorFamilyZbasis<T, Tint>(ListMat1[0]);
    MyMatrix<Tint> SHV2 = ExtractInvariantVectorFamilyZbasis<T, Tint>(ListMat2[0]);

    MyMatrix<T> SHV1_T = UniversalMatrixConversion<T,Tint>(SHV1);
    MyMatrix<T> SHV2_T = UniversalMatrixConversion<T,Tint>(SHV2);

    auto get_equiv=[&]() -> std::optional<MyMatrix<Tint>> {
      if (SHV1_T.rows() != SHV1.rows())
        return {};
      int n_rows = SHV1_T.rows();
      std::vector<T> Vdiag1(n_rows, 0);
      std::vector<T> Vdiag2(n_rows, 0);
      const bool use_scheme = true;
      std::optional<std::vector<Tidx>> opt = TestEquivalence_ListMat_Vdiag<T, Tidx, use_scheme>(SHV1_T, ListMat1, Vdiag1, SHV2_T, ListMat2, Vdiag2);
      if (!opt)
        return {};
      std::optional<MyMatrix<T>> optB = FindMatrixTransformationTest(SHV1_T, SHV2_T, *opt);
      if (!optB) {
        std::cerr << "We have a matrix bug\n";
        throw TerminalException{1};
      }
      MyMatrix<T> const& M_T = *optB;
      if (!IsIntegralMatrix(M_T)) {
        std::cerr << "Bug: The matrix should be integral\n";
        throw TerminalException{1};
      }
      MyMatrix<Tint> M = UniversalMatrixConversion<Tint,T>(M_T);
      return M;
    };
    std::optional<MyMatrix<Tint>> equiv = get_equiv();

    //
    auto prt = [&](std::ostream &os) -> void {
      if (OutFormat == "GAP") {
        if (equiv) {
          os << "return ";
          WriteMatrixGAP(os, *equiv);
          os << ";\n";
        } else {
          os << "return false;\n";
        }
      }
      if (OutFormat == "Oscar") {
        if (equiv) {
          WriteMatrix(os, *equiv);
        } else {
          os << "0\n";
        }
      }
      std::cerr << "Failed to find a matching type for OutFormat=" << OutFormat << "\n";
      throw TerminalException{1};
    };
    if (OutFile == "stderr") {
      prt(std::cerr);
    } else {
      if (OutFile == "stdout") {
        prt(std::cout);
      } else {
        std::ofstream os(OutFile);
        prt(os);
      }
    }
    std::cerr << "Normal termination of LATT_canonicalize\n";
  } catch (TerminalException const &e) {
    std::cerr << "Raised exception led to premature end of LATT_canonicalize\n";
    exit(e.eVal);
  }
  runtime(time);
}
