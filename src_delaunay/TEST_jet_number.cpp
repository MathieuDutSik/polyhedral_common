// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
#include "jet_number.h"
// clang-format on

// Standalone unit tests for the jet<T, N> numeric type. Build with
//   make CHOICE_COMPILATION= TEST_jet_number
// and, to exercise the expensive self-checks,
//   make CHOICE_COMPILATION=-DSANITY_CHECK_EXTENSIVE_JET_NUMBER TEST_jet_number

template <typename T, int N> using J = jet<T, N>;

int main() {
  using T = mpq_class;
  int n_error = 0;
  auto check = [&](bool ok, std::string const &name) {
    std::cerr << (ok ? "ok   " : "FAIL ") << name << "\n";
    if (!ok)
      n_error++;
  };

  // Arithmetic identities at order 3.
  {
    constexpr int Nd = 3;
    J<T, Nd> t = J<T, Nd>::var();
    J<T, Nd> a = J<T, Nd>(2) + J<T, Nd>(3) * t; // 2 + 3 t
    J<T, Nd> b = J<T, Nd>(1) - t;               // 1 - t
    // (a + b) = 3 + 2 t
    J<T, Nd> s = a + b;
    check(s[0] == T(3) && s[1] == T(2) && s[2] == T(0) && s[3] == T(0),
          "sum 2+3t plus 1-t");
    // a * b = (2 + 3t)(1 - t) = 2 + t - 3 t^2
    J<T, Nd> p = a * b;
    check(p[0] == T(2) && p[1] == T(1) && p[2] == T(-3) && p[3] == T(0),
          "product (2+3t)(1-t)");
    // 1/(1 - t) = 1 + t + t^2 + t^3 (+ O(t^4))
    J<T, Nd> inv = inverse(b);
    check(inv[0] == T(1) && inv[1] == T(1) && inv[2] == T(1) && inv[3] == T(1),
          "geometric series 1/(1-t)");
    // (a / b) * b == a exactly (truncated ring identity).
    J<T, Nd> q = a / b;
    J<T, Nd> back = q * b;
    check(back == a, "division round-trip (a/b)*b == a");
    // eval is a ring homomorphism for + (no truncation): eval(a+b) matches.
    T t0(1);
    t0 /= T(10);
    check((a + b).eval(t0) == a.eval(t0) + b.eval(t0), "eval homomorphism +");
  }

  // Leading-coefficient ordering / positivity.
  {
    constexpr int Nd = 2;
    J<T, Nd> t = J<T, Nd>::var();
    J<T, Nd> zero;
    // t > 0 (leading coeff of t is +1).
    check(t > zero, "t > 0 (infinitesimal positive)");
    check(!(t < zero), "not t < 0");
    // -t < 0.
    check((-t) < zero, "-t < 0");
    // (t - t^2) > 0: leading coeff +1.
    J<T, Nd> f = t - t * t;
    check(f > zero, "t - t^2 > 0");
    // (-t + 5 t^2) < 0: leading coeff -1 dominates near 0 despite +5 t^2.
    J<T, Nd> g = -t + J<T, Nd>(5) * t * t;
    check(g < zero, "-t + 5 t^2 < 0 (leading term wins)");
    check(g.sign() == -1, "sign(-t + 5 t^2) == -1");
    // A constant compares by its constant term.
    check(J<T, Nd>(3) > J<T, Nd>(2), "constant 3 > 2");
    // Equal jets.
    check(f == (t - t * t), "equality of identical jets");
    // is_nonnegative rule: first non-zero coeff >= 0; zero jet is >= 0.
    check(zero >= zero, "0 >= 0");
    check(f >= zero, "t - t^2 >= 0");
  }

  // Matrix-determinant lemma sanity: det(1 + t) style via a 2x2 by hand,
  // det([[1+t, t],[t, 1-t]]) = (1+t)(1-t) - t^2 = 1 - 2 t^2.
  {
    constexpr int Nd = 3;
    J<T, Nd> t = J<T, Nd>::var();
    J<T, Nd> m00 = J<T, Nd>(1) + t, m01 = t, m10 = t, m11 = J<T, Nd>(1) - t;
    J<T, Nd> det = m00 * m11 - m01 * m10;
    check(det[0] == T(1) && det[1] == T(0) && det[2] == T(-2) && det[3] == T(0),
          "2x2 jet determinant 1 - 2 t^2");
    // jet_deriv: d^2/dt^2 at 0 of (1 - 2 t^2) = -4.
    check(jet_deriv(det, 0) == T(1) && jet_deriv(det, 1) == T(0) &&
              jet_deriv(det, 2) == T(-4),
          "jet_deriv of 1 - 2 t^2 (0th, 1st, 2nd)");
  }

  // Feasibility of the numeric integral core over jets: MyMatrix<jet<T,N>> with
  // DeterminantMat / Inverse (the operations center_homog and direct_integral
  // rely on). Q + t H with Q = [[2,-1],[-1,2]], H = [[0,1],[1,0]].
  {
    constexpr int Nd = 2;
    using JT = jet<T, Nd>;
    JT t = JT::var();
    MyMatrix<JT> M(2, 2);
    M(0, 0) = JT(2);
    M(0, 1) = JT(-1) + t;
    M(1, 0) = JT(-1) + t;
    M(1, 1) = JT(2);
    // det(Q + tH) = (2)(2) - (-1+t)^2 = 4 - (1 - 2t + t^2) = 3 + 2t - t^2.
    JT det = DeterminantMat(M);
    check(det[0] == T(3) && det[1] == T(2) && det[2] == T(-1),
          "DeterminantMat over jets: 3 + 2t - t^2");
    // Inverse * M == I, checked at order N.
    MyMatrix<JT> Inv = Inverse(M);
    MyMatrix<JT> Prod = Inv * M;
    bool inv_ok = Prod(0, 0) == JT(1) && Prod(1, 1) == JT(1) &&
                  Prod(0, 1) == JT(0) && Prod(1, 0) == JT(0);
    check(inv_ok, "Inverse over jets: Inv * M == I");
  }

  if (n_error == 0)
    std::cerr << "All jet_number tests passed.\n";
  else
    std::cerr << n_error << " jet_number test(s) FAILED.\n";
  return n_error == 0 ? 0 : 1;
}
