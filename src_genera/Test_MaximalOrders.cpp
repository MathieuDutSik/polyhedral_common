// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MaximalOrders.h"
#include <iostream>
#include <vector>
// clang-format on

/*
  The maximal orders of MaximalOrders.h on the orders Z[A] of the companion
  matrices A of monic irreducible polynomials f (a basis of the order: the
  powers of A): Round 2 (CommutativeMaximalOrder) and the search of
  p^{-1} Lambda / Lambda at each prime of the index (MaximalOrderOfSpanAt)
  give orders of the same discriminant, the known discriminant of the
  field. The polynomials: x^2 - 5 (Z[sqrt 5], index 2), x^2 + 3, x^2 + 16
  (Z[4 i], index 4), x^3 - 2 (maximal), x^3 + x^2 - 2 x + 8 (Dedekind's
  field, 2 a common index divisor), Phi_9 (maximal), x^4 - 10 x^2 + 1
  (Z[sqrt 2 + sqrt 3] in Q(sqrt 2, sqrt 3), of discriminant 8 12 24) and
  the minimal polynomial of 2 zeta_5 (Z[2 zeta_5], of index a power of
  2).
 */

using T = mpq_class;
using Tint = mpz_class;

// The companion matrix of f (coefficients by increasing degree, monic).
MyMatrix<T> CompanionMatrix(std::vector<int64_t> const &f) {
  int n = static_cast<int>(f.size()) - 1;
  MyMatrix<T> A = ZeroMatrix<T>(n, n);
  for (int i = 0; i + 1 < n; i++) {
    A(i + 1, i) = 1;
  }
  for (int i = 0; i < n; i++) {
    A(i, n - 1) = T(static_cast<long>(-f[i]));
  }
  return A;
}

// |det Tr(x_i x_j)|, the trace of the n x n matrices being Tr_{K/Q}.
T OrderDiscriminant(std::vector<MyMatrix<T>> const &O) {
  int d = static_cast<int>(O.size());
  MyMatrix<T> G(d, d);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      MyMatrix<T> P = O[i] * O[j];
      T tr(0);
      for (int u = 0; u < P.rows(); u++) {
        tr += P(u, u);
      }
      G(i, j) = tr;
    }
  }
  return T_abs(DeterminantMat(G));
}

int main() {
  try {
    struct Case {
      std::vector<int64_t> f;
      int64_t disc;
    };
    std::vector<Case> cases = {
        {{-5, 0, 1}, 5},
        {{3, 0, 1}, 3},
        {{16, 0, 1}, 4},
        {{-2, 0, 0, 1}, 108},
        {{8, -2, 1, 1}, 503},
        {{1, 0, 0, 1, 0, 0, 1}, 19683},
        {{1, 0, -10, 0, 1}, 2304},
        {{16, 8, 4, 2, 1}, 125},
    };
    size_t n_err = 0;
    for (auto const &c : cases) {
      MyMatrix<T> A = CompanionMatrix(c.f);
      int n = A.rows();
      std::vector<MyMatrix<T>> O{IdentityMat<T>(n)};
      for (int i = 1; i < n; i++) {
        O.push_back(O.back() * A);
      }
      std::vector<MyMatrix<T>> O2 = CommutativeMaximalOrder<T, Tint>(O, n);
      T d2 = OrderDiscriminant(O2);
      // the search, prime by prime over the index
      std::vector<MyMatrix<T>> Os = O;
      T d0 = OrderDiscriminant(O);
      Tint num = UniversalScalarConversion<Tint, T>(d0);
      for (auto const &[p, e] : FactorsIntMap(num)) {
        if (e < 2) {
          continue;
        }
        std::optional<MaximalOrderData<T, Tint>> opt =
            MaximalOrderOfSpanAt<T, Tint>(Os, p, 100000000);
        if (!opt) {
          std::cerr << "the search is beyond its bound\n";
          n_err++;
          continue;
        }
        Os = opt->basis;
      }
      T ds = OrderDiscriminant(Os);
      T expected = UniversalScalarConversion<T, int64_t>(c.disc);
      bool ok = (d2 == expected) && (ds == expected);
      std::cout << "f of degree " << n << ": disc Z[A] = " << d0
                << ", Round 2 " << d2 << ", search " << ds << ", field "
                << expected << (ok ? "" : "  FAIL") << "\n";
      if (!ok) {
        n_err++;
      }
    }
    std::cout << cases.size() << " orders, " << n_err << " errors\n";
    if (n_err > 0) {
      return 1;
    }
    std::cerr << "Normal termination of Test_MaximalOrders\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_MaximalOrders\n";
    exit(e.eVal);
  }
}
