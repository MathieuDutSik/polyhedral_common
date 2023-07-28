// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "POLY_PolytopeFct.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using T = mpq_class;
    if (argc != 2) {
      std::cerr << "POLY_MaxCoefficients [M]\n";
      throw TerminalException{1};
    }
    std::string FileInput = argv[1];
    MyMatrix<T> M = ReadMatrixFile<T>(FileInput);
    T max_coeff = sqr_estimate_facet_coefficients(M);
    double max_coeff_d = UniversalScalarConversion<double,T>(max_coeff);
    double sqr_max = sqrt(max_coeff_d);
    double max_int64_d = std::numeric_limits<int64_t>::max();
    std::cerr << "max_coeff=" << max_coeff << "\n";
    std::cerr << "max_coeff_d=" << max_coeff_d << "\n";
    std::cerr << "sqr_max=" << sqr_max << "\n";
    std::cerr << "max_int64_d=" << max_int64_d << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_MaxCoefficients\n";
    exit(e.eVal);
  }
  runtime(time1);
}
