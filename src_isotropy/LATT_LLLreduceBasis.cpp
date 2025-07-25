// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "Indefinite_LLL.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4 && argc != 2) {
      std::cerr << "LATT_LLLreduceBasis [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_LLLreduceBasis [FileI]\n";
      std::cerr << "or\n";
      std::cerr << "FileI     : The input file\n";
      std::cerr << "OutFormat : Possible values, GAP and Oscar\n";
      std::cerr << "            Default: GAP\n";
      std::cerr << "FileO     : The output file\n";
      std::cerr << "            stdout: write to std::cout\n";
      std::cerr << "            stderr: write to std::cerr\n";
      std::cerr << "            Other filename get written to files\n";
      std::cerr << "            Default: stderr\n";
      throw TerminalException{1};
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
#else
    using T = mpq_class;
    using Tint = mpz_class;
#endif
    std::string FileI = argv[1];
    std::string OutFormat = "GAP";
    std::string FileO = "stderr";
    if (argc == 4) {
      OutFormat = argv[2];
      FileO = argv[3];
    }
    //
    MyMatrix<T> M = ReadMatrixFile<T>(argv[1]);
    std::cerr << "We have M\n";
    //
    auto print_result = [&](std::ostream &os) -> void {
      LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(M, std::cerr);
      MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(rec.Pmat);
      MyMatrix<T> M_Control = B_T * M * B_T.transpose();
      T det = T_abs(DeterminantMat(B_T));
      if (det != 1) {
        std::cerr << "det=" << det << "\n";
        std::cerr << "B_T should have determinant 1 or -1\n";
        throw TerminalException{1};
      }
      if (M_Control != rec.GramMatRed) {
        std::cerr << "M_Control is not what it should be\n";
        throw TerminalException{1};
      }
      if (OutFormat == "GAP") {
        os << "return rec(B:=";
        WriteMatrixGAP(os, rec.Pmat);
        os << ", Mred:=";
        WriteMatrixGAP(os, rec.GramMatRed);
        os << ");\n";
        return;
      }
      if (OutFormat == "Oscar") {
        WriteMatrix(os, rec.Pmat);
        WriteMatrix(os, rec.GramMatRed);
        return;
      }
      std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
                << "\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileO, print_result);
    std::cerr << "Normal termination of LATT_LLLreduceBasis\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_LLLreduceBasis\n";
    exit(e.eVal);
  }
  runtime(time);
}
