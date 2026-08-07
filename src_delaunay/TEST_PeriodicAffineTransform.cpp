// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PeriodicStructures.h"
// clang-format on

/*
  Deterministic tests of the periodic affine transformation algebra against
  the corresponding (n+1) x (n+1) field matrices:
  --- composition and inversion agree with the matrix product / inversion,
  --- the denominators are stable under both operations,
  --- from_field / to_field round trip,
  --- the periodic point set reduction and coset lookup.
 */

using T = mpq_class;
using Tring = mpz_class;
using Ttrans = PeriodicAffineTransform<Tring>;
using Ttraits = transform_traits<Ttrans>;

// Small deterministic generator, enough for test data.
static uint64_t lcg_state = 123456789;
int rand_int(int bound) {
  lcg_state = 6364136223846793005ULL * lcg_state + 1442695040888963407ULL;
  return static_cast<int>((lcg_state >> 33) % bound);
}

MyMatrix<Tring> RandomUnimodular(int n, int n_oper) {
  MyMatrix<Tring> A = IdentityMat<Tring>(n);
  for (int iter = 0; iter < n_oper; iter++) {
    int i = rand_int(n);
    int j = rand_int(n);
    if (i == j) {
      for (int k = 0; k < n; k++) {
        A(i, k) = -A(i, k);
      }
    } else {
      Tring c(rand_int(5) - 2);
      for (int k = 0; k < n; k++) {
        A(i, k) += c * A(j, k);
      }
    }
  }
  return A;
}

Ttrans RandomTransform(int n, Tring const &d) {
  MyMatrix<Tring> A = RandomUnimodular(n, 8);
  MyVector<Tring> w(n);
  for (int i = 0; i < n; i++) {
    w(i) = Tring(rand_int(41) - 20);
  }
  return {std::move(A), std::move(w), d};
}

void check(bool test, std::string const &context) {
  if (!test) {
    std::cerr << "TEST_PeriodicAffineTransform failed: " << context << "\n";
    throw TerminalException{1};
  }
}

int main() {
  try {
    Tring d(12);
    for (int n = 2; n <= 6; n++) {
      for (int iter = 0; iter < 20; iter++) {
        Ttrans x = RandomTransform(n, d);
        Ttrans y = RandomTransform(n, d);
        MyMatrix<T> Mx = Ttraits::to_field<T>(x);
        MyMatrix<T> My = Ttraits::to_field<T>(y);
        // Composition matches the field matrix product, denominator kept.
        Ttrans xy = x * y;
        check(xy.d == d, "denominator stability of the composition");
        MyMatrix<T> Mxy = Mx * My;
        check(TestEqualityMatrix(Ttraits::to_field<T>(xy), Mxy),
              "composition against the field product");
        // Inversion matches the field inversion, denominator kept.
        Ttrans xinv = Inverse(x);
        check(xinv.d == d, "denominator stability of the inversion");
        check(TestEqualityMatrix(Ttraits::to_field<T>(xinv), Inverse(Mx)),
              "inversion against the field inversion");
        Ttrans id = IdentityPeriodicAffineTransform<Tring>(n, d);
        check(x * xinv == id, "x * Inverse(x) is the identity");
        check(xinv * x == id, "Inverse(x) * x is the identity");
        // Round trip through the field matrix.
        Ttrans x_back = Ttraits::from_field(Mx);
        check(x_back == x, "from_field(to_field(x)) equals x");
      }
    }
    // Periodic point set: cosets (0,0), (1/2,1/2) of the checkerboard.
    MyMatrix<T> Cosets(2, 2);
    Cosets(0, 0) = 0;
    Cosets(0, 1) = 0;
    Cosets(1, 0) = T(1) / T(2);
    Cosets(1, 1) = T(1) / T(2);
    PeriodicPointSet<Tring> pps =
        PeriodicPointSetFromRational<Tring>(Cosets);
    check(pps.N == 2, "denominator of the checkerboard");
    MyVector<Tring> u(2);
    u(0) = 5;
    u(1) = -3;
    std::optional<size_t> opt = GetCosetIndex(pps, u);
    check(opt.has_value() && *opt == 1, "coset lookup of (5,-3)/2");
    u(0) = 4;
    u(1) = -3;
    check(!GetCosetIndex(pps, u).has_value(),
          "coset lookup of a point outside the set");
    u(0) = -4;
    u(1) = 6;
    std::optional<size_t> opt2 = GetCosetIndex(pps, u);
    check(opt2.has_value() && *opt2 == 0, "coset lookup of (-4,6)/2");
    std::cerr << "Normal termination of TEST_PeriodicAffineTransform\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TEST_PeriodicAffineTransform\n";
    exit(e.eVal);
  }
}
