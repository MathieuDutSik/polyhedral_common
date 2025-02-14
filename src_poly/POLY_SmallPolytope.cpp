// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "POLY_LinearProgramming.h"
#include "POLY_lrslib.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "POLY_SmallPolytope [n] [k+] [k-]\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      return -1;
    }
    using T = mpq_class;
    //
    int n = ParseScalar<int>(argv[1]);
    int k_p = ParseScalar<int>(argv[2]);
    int k_m = ParseScalar<int>(argv[3]);
    //
    MyMatrix<T> EXT = ZeroMatrix<T>(n + 1, n);
    for (int i = 0; i < n; i++)
      EXT(i, i) = 1;
    for (int i = 0; i < k_p; i++)
      EXT(n, i) = 1;
    for (int i = 0; i < k_m; i++)
      EXT(n, i) = -1;
    MyMatrix<T> NSP = NullspaceMat(EXT);
    if (NSP.rows() != 1) {
      std::cerr << "The rank is incorrect\n";
      throw TerminalException{1};
    }
    vectface vf = lrs::DualDescription_incd(EXT);
    std::cerr << "|vf|=" << vf.size() << "\n";
    int n_p = 0;
    int n_m = 0;
    int n_z = 0;
    for (int i = 0; i <= n; i++) {
      T val = NSP(0, i);
      if (val > 0)
        n_p++;
      if (val < 0)
        n_m++;
      if (val == 0)
        n_z++;
    }
    int n_facet = n_p * n_m + n_z;
    std::cerr << "n_facet=" << n_facet << "\n";
    std::cerr << "Normal termination of POLY_SmallPolytopes\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SmallPolytopes\n";
    exit(e.eVal);
  }
  runtime(time1);
}
