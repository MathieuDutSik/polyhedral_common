// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "FiniteField.h"
#include <iostream>
#include <vector>
// clang-format on

/*
  An exhaustive test of FiniteField.h on the fields F_p^f, p <= 13, f <= 2,
  and p <= 5, f = 3 (one irreducible modulus each): the field axioms on all the elements
  (pairs for the products, triples on a sample for associativity and
  distributivity), the inverses, Frobenius of order f, the traces in F_p
  and additive, the codes a bijection onto [0, q) with the operations on
  the codes (fast for f <= 2) and the tables equal to those on the
  polynomials, and the untrimmed inputs.
 */

int main() {
  try {
    size_t n_err = 0;
    size_t n_field = 0;
    auto check = [&](bool ok, char const *what, FiniteField const &F) {
      if (!ok) {
        n_err++;
        if (n_err < 20) {
          std::cerr << "FAIL " << what << " in F_" << F.p << "^" << F.f << "\n";
        }
      }
    };
    for (T_small_prime p : {2, 3, 5, 7, 11, 13}) {
      for (int f = 1; f <= 3; f++) {
        if (f == 3 && p > 5) {
          continue;
        }
        // the first monic irreducible of degree f, in the order of the codes
        PolyQ mu;
        uint64_t pf = 1;
        for (int i = 0; i < f; i++) {
          pf *= p;
        }
        for (uint64_t c = 0; c < pf; c++) {
          PolyQ g(f + 1, 0);
          uint64_t t = c;
          for (int i = 0; i < f; i++) {
            g[i] = t % p;
            t /= p;
          }
          g[f] = 1;
          if (IsIrreducibleModP(g, p)) {
            mu = g;
            break;
          }
        }
        FiniteField F = MakeFiniteField(p, mu);
        n_field++;
        check(F.q == pf && F.f == f, "q", F);
        std::vector<FqElt> elts;
        for (uint64_t c = 0; c < F.q; c++) {
          FqElt a = FqDecode(F, c);
          check(FqEncode(F, a) == c, "codes", F);
          check(a.empty() || a.back() != 0, "trimmed", F);
          elts.push_back(a);
        }
        FqTables tb = BuildFqTables(F);
        for (uint64_t i = 0; i < F.q; i++) {
          FqElt const &a = elts[i];
          check(FqAdd(F, a, FqNeg(F, a)).empty(), "negation", F);
          check(FqFrobenius(F, a, F.f) == a, "Frobenius order", F);
          check(FqPow(F, a, F.q) == a, "a^q", F);
          check(FqTrace(F, a) < p, "trace", F);
          // untrimmed input: the vector of length f + 2
          FqElt b = a;
          b.resize(F.f + 2, 0);
          check(FqMul(F, b, FqOne()) == a, "untrimmed", F);
          if (i != 0) {
            FqElt ai = FqInv(F, a);
            check(FqMul(F, a, ai) == FqOne(), "inverse", F);
            check(FqCodeInv(F, i) == FqEncode(F, ai), "code inverse", F);
            check(static_cast<uint64_t>(tb.inv[i]) == FqEncode(F, ai),
                  "table inverse", F);
          }
          for (uint64_t j = 0; j < F.q; j++) {
            FqElt const &c = elts[j];
            FqElt s = FqAdd(F, a, c), m = FqMul(F, a, c);
            check(s == FqAdd(F, c, a), "commutative sum", F);
            check(m == FqMul(F, c, a), "commutative product", F);
            check(FqSub(F, s, c) == a, "subtraction", F);
            check(FqCodeAdd(F, i, j) == FqEncode(F, s), "code sum", F);
            check(FqCodeMul(F, i, j) == FqEncode(F, m), "code product", F);
            check(static_cast<uint64_t>(tb.add[i * F.q + j]) == FqEncode(F, s),
                  "table sum", F);
            check(static_cast<uint64_t>(tb.mul[i * F.q + j]) == FqEncode(F, m),
                  "table product", F);
            check(FqTrace(F, s) == ResAdd(FqTrace(F, a), FqTrace(F, c), p),
                  "trace additive", F);
            // associativity and distributivity on a sample of the triples
            for (uint64_t k = 0; k < F.q; k += 1 + F.q / 7) {
              FqElt const &d = elts[k];
              check(FqMul(F, m, d) == FqMul(F, a, FqMul(F, c, d)),
                    "associative", F);
              check(FqMul(F, a, FqAdd(F, c, d)) ==
                        FqAdd(F, m, FqMul(F, a, d)),
                    "distributive", F);
            }
          }
        }
        // f = 2: the formulas of the codes a + p b of a + b omega, omega^2 =
        // s omega - n, of OrderLatticeNeighbor (s = -c_1, n = c_0)
        if (F.f == 2) {
          int64_t ell = p, s = (p - F.mu[1]) % p, nn = F.mu[0];
          for (uint64_t x = 0; x < F.q; x++) {
            for (uint64_t y = 0; y < F.q; y++) {
              int64_t a1 = x % ell, b1 = x / ell, a2 = y % ell, b2 = y / ell;
              int64_t bb = (b1 * b2) % ell;
              int64_t a = ((a1 * a2 - nn * bb) % ell + ell) % ell;
              int64_t b = ((a1 * b2 + a2 * b1 + s * bb) % ell + ell) % ell;
              check(FqCodeMul(F, x, y) == static_cast<uint64_t>(a + ell * b),
                    "order formula", F);
            }
          }
        }
      }
    }
    // a field too large for the codes
    FiniteField Fbig = MakeFiniteField(1000003, PolyQ{2, 0, 0, 0, 1});
    if (Fbig.q != 0) {
      n_err++;
      std::cerr << "FAIL the large field has codes\n";
    }
    std::cout << n_field << " fields, " << n_err << " errors\n";
    if (n_err > 0) {
      return 1;
    }
    std::cerr << "Normal termination of Test_FiniteField\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_FiniteField\n";
    exit(e.eVal);
  }
}
