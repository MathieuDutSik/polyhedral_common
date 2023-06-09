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
    if (argc != 2 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_Isomorphism [ListMat] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "LATT_Isomorphism [ListMat]\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpq_rational;
    using Tfield = T;
    using Tint = boost::multiprecision::mpz_int;
#else
    using T = mpq_class;
    using Tfield = T;
    using Tint = mpz_class;
#endif
    using Tidx = uint32_t;
    //
    std::string FileListMat = argv[1];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      OutFile = argv[3];
    }
    std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(FileListMat);

    MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(ListMat[0]);

    MyMatrix<T> SHV_T = UniversalMatrixConversion<T,Tint>(SHV);
    int n_row = SHV_T.rows();
    std::vector<T> Vdiag(n_row,0);

    const bool use_scheme = true;
    std::vector<std::vector<Tidx>> ListGen =
      GetListGenAutomorphism_ListMat_Vdiag<T, Tfield, Tidx, use_scheme>(SHV_T, ListMat, Vdiag);

    std::vector<MyMatrix<Tint>> ListGenEquiv;
    for (auto & eList : ListGen) {
      std::optional<MyMatrix<T>> opt = FindMatrixTransformationTest(SHV_T, SHV_T, eList);
      if (!opt) {
        std::cerr << "Failed to find the matrix\n";
        throw TerminalException{1};
      }
      MyMatrix<T> const& M_T = *opt;
      if (!IsIntegralMatrix(M_T)) {
        std::cerr << "Bug: The matrix should be integral\n";
        throw TerminalException{1};
      }
      MyMatrix<Tint> M = UniversalMatrixConversion<Tint,T>(M_T);
      ListGenEquiv.push_back(M);
    }
    //
    auto prt = [&](std::ostream &os) -> void {
      if (OutFormat == "GAP") {
        os << "return [";
        bool IsFirst = true;
        for (auto & eMat : ListGenEquiv) {
          if (!IsFirst)
            os << ",\n";
          IsFirst = false;
          WriteMatrixGAP(os, eMat);
        }
        os << "];\n";
        return;
      }
      if (OutFormat == "Oscar") {
        os << ListGenEquiv.size() << "\n";
        for (auto & eMat : ListGenEquiv) {
          WriteMatrix(os, eMat);
        }
        return;
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
    std::cerr << "Normal termination of LATT_Automorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Raised exception led to premature end of LATT_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time);
}
