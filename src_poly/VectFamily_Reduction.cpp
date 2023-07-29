// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "POLY_PolytopeFct.h"
#include "LatticeDefinitions.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using T = mpq_class;
    if (argc != 3 && argc != 5) {
      std::cerr << "VectFamily_Reduction [FileI] [method] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "VectFamily_Reduction [FileI] [method]\n";
      std::cerr << "\n";
      std::cerr << " ------- method --------\n";
      std::cerr << "Possible values for method\n";
      std::cerr << "direct : apply the direct LLL method\n";
      std::cerr << "dual   : apply the dual LLL method\n";
      std::cerr << "\n";
      std::cerr << " ------- OutFormat --------\n";
      std::cerr << "Possible values for OutFormat\n";
      std::cerr << "GAP : for writing in GAP readable file\n";
      std::cerr << "CPP : for writing in CPP polyhedral readable file\n";
      std::cerr << "\n";
      std::cerr << " ------- FileO --------\n";
      std::cerr << "Possible values for FileO\n";
      std::cerr << "stderr : for writing to std::cerr\n";
      std::cerr << "stdout : for writing in std::cout\n";
      std::cerr << "otherwise written to the named file in output\n";
      throw TerminalException{1};
    }
    std::string FileInput = argv[1];
    MyMatrix<T> M = ReadMatrixFile<T>(FileInput);
    std::string method = argv[2];
    std::cerr << "method=" << method << "\n";
    std::string OutFormat = "CPP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    auto matrix_measure=[&](MyMatrix<T> const& Minp) -> void {
      T max_coeff = sqr_estimate_facet_coefficients(Minp);
      double max_coeff_d = UniversalScalarConversion<double,T>(max_coeff);
      double sqr_max = sqrt(max_coeff_d);
      double max_int64_d = std::numeric_limits<int64_t>::max();
      T l1_norm = L1_norm_mat(Minp);
      std::cerr << "max_coeff=" << max_coeff << "\n";
      std::cerr << "max_coeff_d=" << max_coeff_d << "\n";
      std::cerr << "sqr_max=" << sqr_max << "\n";
      std::cerr << "max_int64_d=" << max_int64_d << "\n";
      std::cerr << "L1(M)=" << l1_norm << "\n";
    };
    std::cerr << "Original complexity measures\n";
    matrix_measure(M);
    std::pair<MyMatrix<T>,MyMatrix<T>> pair = ReduceVectorFamily(M, method);
    std::cerr << "Output complexity measures\n";
    matrix_measure(pair.first);
    auto print_mat=[&](std::ostream & os) -> void {
      if (OutFormat == "GAP") {
        os << "return ";
        WriteMatrixGAP(os, pair.first);
        os << ";\n";
        return;
      }
      if (OutFormat == "CPP") {
        return WriteMatrix(os, pair.first);
      }
      std::cerr << "No matching format in print_mat\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      print_mat(std::cerr);
    } else {
      if (FileO == "stdout") {
        print_mat(std::cout);
      } else {
        std::ofstream os(FileO);
        print_mat(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in VectFamily_Reduction\n";
    exit(e.eVal);
  }
  runtime(time1);
}
