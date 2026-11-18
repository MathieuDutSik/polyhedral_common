// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "IsometryGlueFactors.h"
#include <random>
#include <set>
#include <string>
#include <vector>
// clang-format on

/*
  Tests of the arithmetic of IsometryGlueFactors.h against independent
  computations:
  - CyclotomicPolynomialCoeffs: prod_{e | d} Phi_e = x^d - 1 and the degree
    phi(d), for d <= 60;
  - PolyResultant: |Res(Phi_a, Phi_b)| for a, b <= 40 against the closed
    form (Apostol): p^phi(b) when a / b is a power of the prime p (a > b),
    Phi_a(1) for b = 1, 1 otherwise; Res(x - 1, x + 1) = 2;
  - MinimalPolynomialCoeffs: on block diagonal matrices of companion
    matrices of cyclotomic polynomials, conjugated by a unimodular matrix,
    against the product of the distinct factors;
  - SublatticeDeterminants: against the definition, the D with s | D and
    D | d s for some s | S, by a search over D <= d S;
  - IsAdmissibleGlue: on the involutions of E8 (d+ = d- = 2^a, glue 2^a)
    and on cases that fail each condition;
  - EvenSublatticeBasis: on random odd lattices, a sublattice of index 2
    with even norms;
  - OddOverlatticeOrbits: on D_n, the even sublattice of I_n, with sigma =
    1 and no group, the odd overlattices of D(D_n) of order 2: the vector
    class always, the two spinor classes (of norm n / 4) when n / 4 is an
    odd integer: 3 for n = 4, 12, 1 for n = 5, 6, 8; all unimodular.
 */

using T = mpq_class;
using Tint = mpz_class;

std::vector<T> PolyPower(std::vector<T> const &p, int e) {
  std::vector<T> res{T(1)};
  for (int i = 0; i < e; i++) {
    res = PolyProductCoeffs(res, p);
  }
  return res;
}

T PolyEval(std::vector<T> const &p, T const &x) {
  T res(0);
  for (size_t i = p.size(); i-- > 0;) {
    res = res * x + p[i];
  }
  return res;
}

int EulerPhi(int n) {
  int res = 0;
  for (int k = 1; k <= n; k++) {
    if (std::gcd(k, n) == 1) {
      res++;
    }
  }
  return res;
}

// p when n = p^k, k >= 1, 0 otherwise
int PrimeOfPower(int n) {
  if (n < 2) {
    return 0;
  }
  int p = 2;
  while (n % p != 0) {
    p++;
  }
  while (n % p == 0) {
    n /= p;
  }
  return (n == 1) ? p : 0;
}

MyMatrix<T> Companion(std::vector<T> const &p) {
  int m = static_cast<int>(p.size()) - 1;
  MyMatrix<T> C = ZeroMatrix<T>(m, m);
  for (int i = 0; i + 1 < m; i++) {
    C(i + 1, i) = 1;
  }
  for (int i = 0; i < m; i++) {
    C(i, m - 1) = -p[i];
  }
  return C;
}

MyMatrix<T> BlockDiag(std::vector<MyMatrix<T>> const &blocks) {
  int n = 0;
  for (auto const &B : blocks) {
    n += B.rows();
  }
  MyMatrix<T> M = ZeroMatrix<T>(n, n);
  int pos = 0;
  for (auto const &B : blocks) {
    for (int i = 0; i < B.rows(); i++) {
      for (int j = 0; j < B.rows(); j++) {
        M(pos + i, pos + j) = B(i, j);
      }
    }
    pos += B.rows();
  }
  return M;
}

int main() {
  try {
    int n_error = 0;
    auto fail = [&](std::string const &what) -> void {
      std::cerr << "ERROR: " << what << "\n";
      n_error++;
    };
    // the cyclotomic polynomials
    for (int d = 1; d <= 60; d++) {
      std::vector<T> prod{T(1)};
      for (int e = 1; e <= d; e++) {
        if (d % e == 0) {
          prod = PolyProductCoeffs(prod, CyclotomicPolynomialCoeffs<T>(e));
        }
      }
      std::vector<T> target(d + 1, T(0));
      target[0] = -1;
      target[d] = 1;
      if (prod != target) {
        fail("prod Phi_e != x^d - 1 for d=" + std::to_string(d));
      }
      if (static_cast<int>(CyclotomicPolynomialCoeffs<T>(d).size()) - 1 !=
          EulerPhi(d)) {
        fail("deg Phi_d != phi(d) for d=" + std::to_string(d));
      }
    }
    // the resultants
    for (int a = 1; a <= 40; a++) {
      for (int b = 1; b < a; b++) {
        T res = T_abs(PolyResultant(CyclotomicPolynomialCoeffs<T>(a),
                                    CyclotomicPolynomialCoeffs<T>(b)));
        T expected(1);
        if (b == 1) {
          expected = T_abs(PolyEval(CyclotomicPolynomialCoeffs<T>(a), T(1)));
        } else if (a % b == 0 && PrimeOfPower(a / b) != 0) {
          expected = MyPow<T>(T(PrimeOfPower(a / b)), EulerPhi(b));
        }
        if (res != expected) {
          fail("Res(Phi_" + std::to_string(a) + ", Phi_" + std::to_string(b) +
               ")");
        }
        T res_sym = T_abs(PolyResultant(CyclotomicPolynomialCoeffs<T>(b),
                                        CyclotomicPolynomialCoeffs<T>(a)));
        if (res_sym != res) {
          fail("the resultant is not symmetric up to sign");
        }
      }
    }
    if (T_abs(PolyResultant(std::vector<T>{T(-1), T(1)},
                            std::vector<T>{T(1), T(1)})) != 2) {
      fail("Res(x - 1, x + 1) != 2");
    }
    // the minimal polynomials, of block diagonal matrices conjugated
    std::vector<std::vector<int>> cases{
        {1}, {2}, {1, 2}, {3, 3}, {1, 3, 3}, {4, 2, 1},
        {5, 10}, {6, 3, 2, 1}, {12, 4}, {8, 8, 1}, {7}, {9, 3}};
    MyMatrix<T> U = IdentityMat<T>(12);
    for (int i = 0; i + 1 < 12; i++) {
      U(i, i + 1) = (i % 3) - 1;
    }
    for (auto const &c : cases) {
      std::vector<MyMatrix<T>> blocks;
      std::set<int> distinct;
      for (int d : c) {
        blocks.push_back(Companion(CyclotomicPolynomialCoeffs<T>(d)));
        distinct.insert(d);
      }
      MyMatrix<T> M = BlockDiag(blocks);
      int n = M.rows();
      MyMatrix<T> Un = U.block(0, 0, n, n);
      MyMatrix<T> Mc = Un * M * Inverse(Un);
      std::vector<T> expected{T(1)};
      for (int d : distinct) {
        expected = PolyProductCoeffs(expected, CyclotomicPolynomialCoeffs<T>(d));
      }
      if (MinimalPolynomialCoeffs(Mc) != expected) {
        fail("the minimal polynomial of a block diagonal matrix");
      }
    }
    if (MinimalPolynomialCoeffs(MyMatrix<T>(-IdentityMat<T>(4))) !=
        std::vector<T>{T(1), T(1)}) {
      fail("the minimal polynomial of -Id");
    }
    // a power of a polynomial: the minimal polynomial of a Jordan block
    {
      MyMatrix<T> J = Companion(PolyPower({T(-1), T(1)}, 3));
      if (MinimalPolynomialCoeffs(J) != PolyPower({T(-1), T(1)}, 3)) {
        fail("the minimal polynomial of a unipotent companion matrix");
      }
    }
    // the determinants of the sublattices, against the definition
    for (int d : {1, 2, 3, 4, 6, 9, 12, 81}) {
      for (int res : {1, 2, 3, 4, 5, 6}) {
        for (int m = 0; m <= 3; m++) {
          std::vector<Tint> got =
              SublatticeDeterminants<Tint>(Tint(d), MyPow<Tint>(Tint(res), m));
          Tint S = MyPow<Tint>(Tint(res), m);
          std::vector<Tint> brute;
          for (Tint D = 1; D <= Tint(d) * S; D++) {
            bool ok = false;
            for (Tint s = 1; s <= S; s++) {
              if (S % s == 0 && D % s == 0 && (Tint(d) * s) % D == 0) {
                ok = true;
              }
            }
            if (ok) {
              brute.push_back(D);
            }
          }
          if (got != brute) {
            fail("SublatticeDeterminants d=" + std::to_string(d) +
                 " res=" + std::to_string(res) + " m=" + std::to_string(m));
          }
        }
      }
    }
    // the admissible glues
    for (int a = 0; a <= 4; a++) {
      T da = MyPow<T>(T(2), a);
      if (!IsAdmissibleGlue<T, Tint>(da, da, T(1), Tint(16))) {
        fail("an involution of E8 is not admissible");
      }
    }
    if (IsAdmissibleGlue<T, Tint>(T(2), T(8), T(1), Tint(16))) {
      fail("a glue of order 4 in D(X) of order 2 is admissible");
    }
    if (IsAdmissibleGlue<T, Tint>(T(4), T(4), T(1), Tint(2))) {
      fail("a glue beyond S is admissible");
    }
    if (IsAdmissibleGlue<T, Tint>(T(2), T(4), T(1), Tint(16))) {
      fail("a glue of non square order is admissible");
    }
    // the even sublattices
    std::mt19937 rng(17);
    for (int iter = 0; iter < 100; iter++) {
      int n = 1 + static_cast<int>(rng() % 6);
      MyMatrix<T> A = ZeroMatrix<T>(n, n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
          int x = static_cast<int>(rng() % 5) - 2;
          A(i, j) = x;
          A(j, i) = x;
        }
        A(i, i) = static_cast<int>(rng() % 7);
      }
      A(0, 0) = 1 + 2 * static_cast<int>(rng() % 3);
      MyMatrix<T> B = EvenSublatticeBasis(A);
      MyMatrix<T> A0 = B * A * B.transpose();
      if (T_abs(DeterminantMat(B)) != 2) {
        fail("the even sublattice is not of index 2");
      }
      for (int i = 0; i < n; i++) {
        if (!IsInteger(A0(i, i) / T(2))) {
          fail("the even sublattice has an odd norm");
        }
      }
    }
    // the odd overlattices of D_n
    for (auto const &[n, expected] :
         std::vector<std::pair<int, size_t>>{{4, 3}, {5, 1}, {6, 1}, {8, 1}, {12, 3}}) {
      MyMatrix<T> In = IdentityMat<T>(n);
      MyMatrix<T> B = EvenSublatticeBasis(In);
      MyMatrix<T> A0 = B * In * B.transpose();
      std::vector<OddOverlattice<T>> ov = OddOverlatticeOrbits<T, Tint>(
          A0, std::vector<MyMatrix<T>>{IdentityMat<T>(n)},
          std::vector<MyMatrix<Tint>>{});
      if (ov.size() != expected) {
        fail("the odd overlattices of D_" + std::to_string(n) + ": " +
             std::to_string(ov.size()));
      }
      for (auto const &o : ov) {
        MyMatrix<T> AL = o.basis * A0 * o.basis.transpose();
        bool is_odd = false;
        for (int i = 0; i < n; i++) {
          is_odd = is_odd || !IsInteger(AL(i, i) / T(2));
        }
        if (!IsIntegralMatrix(AL) || T_abs(DeterminantMat(AL)) != 1 ||
            !is_odd || o.orbit_size != 1) {
          fail("an odd overlattice of D_" + std::to_string(n));
        }
      }
    }
    std::cerr << "n_error=" << n_error << "\n";
    if (n_error > 0) {
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of Test_IsometryGlueFactors\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_IsometryGlueFactors\n";
    exit(e.eVal);
  }
}
