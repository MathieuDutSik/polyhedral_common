// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MaximalOrders.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

/*
  The maximal order of a commutative order given by a polynomial:

    GENUS_MaximalOrder c0,c1,...,cn

  prints "disc_order disc_maximal_order" for the order Z[x] / (f), f =
  c_0 + c_1 x + ... + x^n monic irreducible, and its maximal order
  (CommutativeMaximalOrder, MaximalOrders.h): the discriminants of their
  trace forms.
 */

// The integers of a comma separated list.
inline std::vector<int64_t> ParseIntList(std::string const &lst) {
  std::vector<int64_t> res;
  std::stringstream ss(lst);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (!tok.empty()) {
      res.push_back(std::stoll(tok));
    }
  }
  return res;
}

/*
  The order Z[x] / (f) of a monic irreducible f = c_0 + c_1 x + ... + x^n,
  as the powers of its companion matrix, and its maximal order
  (CommutativeMaximalOrder): "disc_order disc_maximal_order", the
  discriminants of their trace forms.
 */
template <typename T, typename Tint>
void PrintMaximalOrder(std::vector<int64_t> const &f) {
  int n = static_cast<int>(f.size()) - 1;
  MyMatrix<T> A = ZeroMatrix<T>(n, n);
  for (int i = 0; i + 1 < n; i++) {
    A(i + 1, i) = 1;
  }
  for (int i = 0; i < n; i++) {
    A(i, n - 1) = UniversalScalarConversion<T, int64_t>(-f[i]);
  }
  std::vector<MyMatrix<T>> O{IdentityMat<T>(n)};
  for (int i = 1; i < n; i++) {
    O.push_back(O.back() * A);
  }
  auto disc_of = [&](std::vector<MyMatrix<T>> const &Ob) -> T {
    MyMatrix<T> G(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        MyMatrix<T> P = Ob[i] * Ob[j];
        T tr(0);
        for (int u = 0; u < n; u++) {
          tr += P(u, u);
        }
        G(i, j) = tr;
      }
    }
    return T_abs(DeterminantMat(G));
  };
  std::cout << disc_of(O) << " "
            << disc_of(CommutativeMaximalOrder<T, Tint>(O, n)) << "\n";
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_MaximalOrder c0,c1,...,cn\n";
      return -1;
    }
    PrintMaximalOrder<mpq_class, mpz_class>(ParseIntList(argv[1]));
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_MaximalOrder\n";
    exit(e.eVal);
  }
}
