// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Compare the jet-based moment derivative with the 2n+1-sample interpolation,
// and report the runtime of each. Emits one summary line on stdout:
//   RESULT <file> n=<n> match=<0|1> interp_ms=<..> jet_ms=<..>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "QuantizationDeformation.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on
#include <chrono>

using T = mpq_class;
using Tint = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, mpz_class>;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "TEST_moment_jet [Gram.file]\n";
    return -1;
  }
  std::string file = argv[1];
  MyMatrix<T> Q = ReadMatrixFile<T>(file);
  int n = Q.rows();
  // A rank-one direction B = v v^T with v = e_0 - e_1.
  MyVector<T> v = ZeroVector<T>(n);
  v(0) = 1;
  v(1) = -1;
  MyMatrix<T> B = v * v.transpose();

  auto now = []() { return std::chrono::steady_clock::now(); };
  auto ms = [](auto d) {
    return std::chrono::duration<double, std::milli>(d).count();
  };

  auto t0 = now();
  MyMatrix<T> DM_interp =
      compute_moment_derivative<T, Tint, Tgroup>(Q, B, std::cerr);
  auto t1 = now();
  MyMatrix<T> DM_jet =
      compute_moment_derivative_jet<T, Tint, Tgroup>(Q, B, std::cerr);
  auto t2 = now();

  bool match = (DM_interp == DM_jet);
  if (!match) {
    std::cerr << "MISMATCH on " << file << "\nDM_interp =\n"
              << DM_interp << "\nDM_jet =\n"
              << DM_jet << "\n";
  }
  std::cout << "RESULT " << file << " n=" << n << " match=" << (match ? 1 : 0)
            << " interp_ms=" << ms(t1 - t0) << " jet_ms=" << ms(t2 - t1)
            << "\n";
  return match ? 0 : 1;
}
