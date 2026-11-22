// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "HermitianGenus.h"
#include <string>
#include <vector>
// clang-format on

/*
  The genera of positive definite hermitian lattices over Q(zeta_d) /
  Q(zeta_d)^+ (HermitianGenus.h), for comparison with hermitian_genera of
  Hecke.jl:
    count d b N : a line "d b N number", the genera of rank b, scales
                  between O_E and N O_E
    decomposition : the primes of K_0 over l and their type in E
                  (HermPrimesAbove, from the Dirichlet characters) against
                  the explicit formulas below, for d <= 150 and l < 100
 */

// The explicit formulas: d = l^a m, e_E = phi(l^a), f = ord(l mod m), the
// primes of E over l ramified over K_0 when d = l^a, inert when -1 is a
// power of l modulo m, split otherwise.
HermPrimeData HermPrimesAboveFormula(int64_t d_in, int64_t l) {
  int64_t d = NormalizeCyclotomicIndex(d_in);
  int64_t la = HermPrimePowerPart(d, l);
  int64_t m = d / la;
  HermPrimeData res;
  res.l = l;
  res.dyadic = (l == 2);
  int64_t e_E = EulerPhi(la);
  int64_t f = MultiplicativeOrder(l, m);
  int64_t g_E = EulerPhi(m) / f;
  if (la > 1 && m == 1) {
    res.kind = 2;
    res.n_prime = g_E;
    res.e0 = e_E / 2;
  } else if (IsPowerModulo(-1, l, m)) {
    res.kind = 1;
    res.n_prime = g_E;
    res.e0 = e_E;
  } else {
    res.kind = 0;
    res.n_prime = g_E / 2;
    res.e0 = e_E;
  }
  return res;
}

int main(int argc, char *argv[]) {
  try {
    if (argc == 2 && std::string(argv[1]) == "decomposition") {
      size_t n_case = 0, n_err = 0;
      for (int64_t d = 3; d <= 150; d++) {
        if (d % 4 == 2) {
          continue;
        }
        for (int64_t l = 2; l < 100; l++) {
          if (!IsPrime(l)) {
            continue;
          }
          HermPrimeData a = HermPrimesAbove(d, l);
          HermPrimeData b = HermPrimesAboveFormula(d, l);
          n_case++;
          if (a.kind != b.kind || a.n_prime != b.n_prime || a.e0 != b.e0) {
            n_err++;
            std::cerr << "d=" << d << " l=" << l << ": kind " << a.kind
                      << " / " << b.kind << ", primes " << a.n_prime << " / "
                      << b.n_prime << ", e0 " << a.e0 << " / " << b.e0
                      << "\n";
          }
        }
      }
      std::cout << "cases=" << n_case << " errors=" << n_err << "\n";
      return n_err == 0 ? 0 : 1;
    }
    if (argc != 5 || std::string(argv[1]) != "count") {
      std::cerr << "Test_HermitianGenus count d b N\n";
      std::cerr << "Test_HermitianGenus decomposition\n";
      return -1;
    }
    int64_t d = std::stoll(argv[2]);
    int b = std::stoi(argv[3]);
    int64_t N = std::stoll(argv[4]);
    std::vector<HermPrimeSymbols> ps = HermitianPrimeSymbols(d, b, N);
    std::cout << d << " " << b << " " << N << " "
              << HermitianGlobalGenera(ps).size() << "\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_HermitianGenus\n";
    exit(e.eVal);
  }
}
