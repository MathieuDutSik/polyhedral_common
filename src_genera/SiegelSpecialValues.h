// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_SIEGELSPECIALVALUES_H_
#define SRC_GENERA_SIEGELSPECIALVALUES_H_

// clang-format off
#include "NumberTheory.h"
#include "factorizations.h"
#include <cstdint>
#include <map>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_SIEGEL_SPECIAL_VALUES
#endif

/*
  The special values that Siegel's mass and representation formulas are
  made of, computed exactly.

  The values zeta(2k), L(k, chi_D) and Gamma(j/2) are transcendental, but
  in every formula they enter, the transcendental parts cancel and the
  result is rational. Rather than evaluating them numerically, or handing
  them to a symbolic algebra system to simplify, each is carried as a triple

      q * pi^a * sqrt(r),        q, a, r rational, r > 0,

  which is closed under products and inverses. A computation multiplies
  such triples and asserts at the end that a = 0 and that r is the square
  of a rational, which makes the cancellation a checked property rather
  than an assumption.

  Conventions: B_k are the Bernoulli numbers with B_1 = -1/2; chi_D is the
  Kronecker symbol (D / .) of a fundamental discriminant D, a primitive
  quadratic character of conductor |D|; D = 1 stands for the trivial
  character, whose L-function is zeta.
 */

//
// Integers: Kronecker symbol, square-free parts, fundamental discriminants
//

// A 64-bit integer as an element of T. int64_t is long long on some
// platforms, for which the GMP and FLINT classes have no constructor; long
// is 64 bits wide on every LP64 platform this code targets.
template <typename T> T ScalarFromInt64(int64_t x) {
  static_assert(sizeof(long) == sizeof(int64_t), "long must be 64 bits");
  return T(static_cast<long>(x));
}

// The Kronecker symbol (a / n) for n > 0.
inline int KroneckerSymbol(int64_t a, int64_t n) {
#ifdef SANITY_CHECK_SIEGEL_SPECIAL_VALUES
  if (n <= 0) {
    std::cerr << "SSV: KroneckerSymbol needs n > 0, got n=" << n << "\n";
    throw TerminalException{1};
  }
#endif
  int result = 1;
  while (n % 2 == 0) {
    n /= 2;
    if (a % 2 == 0) {
      return 0;
    }
    int64_t r = ((a % 8) + 8) % 8;
    if (r == 3 || r == 5) {
      result = -result;
    }
  }
  // n odd and positive: the Jacobi symbol
  int64_t aa = ((a % n) + n) % n;
  int64_t nn = n;
  while (aa != 0) {
    while (aa % 2 == 0) {
      aa /= 2;
      int64_t r = nn % 8;
      if (r == 3 || r == 5) {
        result = -result;
      }
    }
    std::swap(aa, nn);
    if (aa % 4 == 3 && nn % 4 == 3) {
      result = -result;
    }
    aa %= nn;
  }
  return (nn == 1) ? result : 0;
}

// (-1)^r x.
inline int64_t SignedPower(int r, int64_t x) { return (r % 2 == 0) ? x : -x; }

// The square-free part of a non-zero integer, with its sign.
inline int64_t SquareFreePart(int64_t x) {
  if (x == 0) {
    std::cerr << "SSV: SquareFreePart called on zero\n";
    throw TerminalException{1};
  }
  int64_t sign = (x < 0) ? -1 : 1;
  int64_t abs_x = (x < 0) ? -x : x;
  int64_t sqf = 1;
  if (abs_x > 1) {
    std::map<int64_t, size_t> fact = FactorsIntMap<int64_t>(abs_x);
    for (auto &kv : fact) {
      if (kv.second % 2 == 1) {
        sqf *= kv.first;
      }
    }
  }
  return sign * sqf;
}

// The discriminant of Q(sqrt(x)) for x non-zero; 1 when x is a square.
inline int64_t FundamentalDiscriminant(int64_t x) {
  int64_t d = SquareFreePart(x);
  if (d == 1) {
    return 1;
  }
  int64_t r = ((d % 4) + 4) % 4;
  return (r == 1) ? d : 4 * d;
}

// chi_D(a) for a fundamental discriminant D (D = 1: the trivial character).
inline int QuadraticCharacter(int64_t D, int64_t a) {
  if (D == 1) {
    return 1;
  }
  if (a > 0) {
    return KroneckerSymbol(D, a);
  }
  // chi_D(-1) is the sign of D
  int sign = (D > 0) ? 1 : -1;
  return sign * KroneckerSymbol(D, -a);
}

//
// Bernoulli numbers
//

template <typename T> T BinomialCoefficient(int n, int k) {
  if (k < 0 || k > n) {
    return T(0);
  }
  T res(1);
  for (int i = 1; i <= k; i++) {
    res *= T(n - k + i);
    res /= T(i);
  }
  return res;
}

// B_0, ..., B_kmax with B_1 = -1/2, from sum_{j=0}^{m} C(m+1, j) B_j = 0.
template <typename T> std::vector<T> BernoulliNumbers(int kmax) {
  std::vector<T> B(kmax + 1, T(0));
  B[0] = T(1);
  for (int m = 1; m <= kmax; m++) {
    T sum(0);
    for (int j = 0; j < m; j++) {
      sum += BinomialCoefficient<T>(m + 1, j) * B[j];
    }
    B[m] = -sum / T(m + 1);
  }
  return B;
}

/*
  The generalized Bernoulli number B_{k, chi_D} of the quadratic character of
  conductor f = |D|:

      B_{k,chi} = f^{k-1} sum_{a=1}^{f} chi(a) B_k(a/f)
                = sum_{j=0}^{k} C(k,j) B_j f^{j-1} S_{k-j},
      S_t = sum_{a=1}^{f} chi(a) a^t .

  For D = 1 this is B_k with the convention B_1 = +1/2 of the trivial
  character. The cost is O(f k), so the conductor is what bounds its use.
 */
template <typename T> T GeneralizedBernoulli(int k, int64_t D) {
  std::vector<T> B = BernoulliNumbers<T>(k);
  if (D == 1) {
    if (k == 1) {
      return T(1) / T(2);
    }
    return B[k];
  }
  int64_t f = (D < 0) ? -D : D;
  std::vector<T> S(k + 1, T(0));
  for (int64_t a = 1; a <= f; a++) {
    int c = KroneckerSymbol(D, a);
    if (c == 0) {
      continue;
    }
    T a_T = ScalarFromInt64<T>(a);
    T pw(1);
    for (int t = 0; t <= k; t++) {
      if (c == 1) {
        S[t] += pw;
      } else {
        S[t] -= pw;
      }
      pw *= a_T;
    }
  }
  T f_T = ScalarFromInt64<T>(f);
  T res(0);
  // f^{j-1} for j = 0 is 1/f
  T fpow = T(1) / f_T;
  for (int j = 0; j <= k; j++) {
    res += BinomialCoefficient<T>(k, j) * B[j] * fpow * S[k - j];
    fpow *= f_T;
  }
  return res;
}

//
// The q * pi^a * sqrt(r) representation
//

template <typename T> struct PiRadical {
  T q;
  T pi_exp;
  T rad;
};

template <typename T> PiRadical<T> PiRadicalRational(T const &q) {
  return {q, T(0), T(1)};
}

template <typename T>
PiRadical<T> operator*(PiRadical<T> const &x, PiRadical<T> const &y) {
  return {x.q * y.q, x.pi_exp + y.pi_exp, x.rad * y.rad};
}

template <typename T> PiRadical<T> PiRadicalInverse(PiRadical<T> const &x) {
  if (x.q == 0) {
    std::cerr << "SSV: inverse of a zero special value\n";
    throw TerminalException{1};
  }
  return {T(1) / x.q, -x.pi_exp, T(1) / x.rad};
}

template <typename T>
PiRadical<T> PiRadicalPower(PiRadical<T> const &x, int e) {
  PiRadical<T> res = PiRadicalRational<T>(T(1));
  PiRadical<T> base = (e >= 0) ? x : PiRadicalInverse(x);
  int abs_e = (e >= 0) ? e : -e;
  for (int i = 0; i < abs_e; i++) {
    res = res * base;
  }
  return res;
}

// x^(num/2) for a positive rational x, as q * sqrt(r).
template <typename T> PiRadical<T> RationalHalfPower(T const &x, int num) {
  int abs_num = (num >= 0) ? num : -num;
  T q(1);
  for (int i = 0; i < abs_num / 2; i++) {
    q *= x;
  }
  T rad = (abs_num % 2 == 1) ? x : T(1);
  PiRadical<T> res = {q, T(0), rad};
  return (num >= 0) ? res : PiRadicalInverse(res);
}

// The rational value of x, which must have pi exponent 0 and a radical that
// is the square of a rational.
template <typename T> T PiRadicalToRational(PiRadical<T> const &x) {
  if (x.pi_exp != 0) {
    std::cerr << "SSV: the value is not rational, pi exponent " << x.pi_exp
              << "\n";
    throw TerminalException{1};
  }
  std::optional<T> opt = UniversalSquareRoot(x.rad);
  if (!opt) {
    std::cerr << "SSV: the value is not rational, radical " << x.rad
              << " is not a square\n";
    throw TerminalException{1};
  }
  return x.q * *opt;
}

//
// Gamma, zeta, L
//

// Gamma(j/2) for j >= 1.
template <typename T> PiRadical<T> GammaHalfInteger(int j) {
  if (j < 1) {
    std::cerr << "SSV: GammaHalfInteger needs j >= 1, got j=" << j << "\n";
    throw TerminalException{1};
  }
  if (j % 2 == 0) {
    T fact(1);
    for (int i = 2; i < j / 2; i++) {
      fact *= T(i);
    }
    return PiRadicalRational<T>(fact);
  }
  // Gamma(j/2) = (j-2)!! / 2^{(j-1)/2} * sqrt(pi)
  T dfact(1);
  for (int i = j - 2; i >= 1; i -= 2) {
    dfact *= T(i);
  }
  T pow2(1);
  for (int i = 0; i < (j - 1) / 2; i++) {
    pow2 *= T(2);
  }
  return {dfact / pow2, T(1) / T(2), T(1)};
}

// zeta(k) for k >= 2 even: (-1)^{k/2+1} B_k (2 pi)^k / (2 k!).
template <typename T> PiRadical<T> ZetaEvenValue(int k) {
  if (k < 2 || k % 2 != 0) {
    std::cerr << "SSV: ZetaEvenValue needs an even k >= 2, got k=" << k
              << "\n";
    throw TerminalException{1};
  }
  std::vector<T> B = BernoulliNumbers<T>(k);
  T q = B[k];
  for (int i = 0; i < k; i++) {
    q *= T(2);
  }
  T fact(1);
  for (int i = 2; i <= k; i++) {
    fact *= T(i);
  }
  q /= T(2) * fact;
  if ((k / 2 + 1) % 2 == 1) {
    q = -q;
  }
  return {q, T(k), T(1)};
}

/*
  L(k, chi_D) for D = 1 (then k must be even, and this is zeta(k)) or D a
  fundamental discriminant with chi_D(-1) = (-1)^k, the parity at which the
  value is a rational multiple of pi^k sqrt(|D|):

      L(k, chi) = (-1)^{1 + (k - delta)/2} * sqrt(f)/2 * (2 pi / f)^k
                  * B_{k,chi} / k!,      delta = 0 if D > 0, 1 if D < 0.
 */
template <typename T> PiRadical<T> LQuadraticValue(int k, int64_t D) {
  if (D == 1) {
    return ZetaEvenValue<T>(k);
  }
  int delta = (D > 0) ? 0 : 1;
  if ((k - delta) % 2 != 0 || k < 1) {
    std::cerr << "SSV: L(k, chi_D) at k=" << k << " D=" << D
              << " is not at a critical point of the right parity\n";
    throw TerminalException{1};
  }
  int64_t f = (D < 0) ? -D : D;
  T f_T = ScalarFromInt64<T>(f);
  T q = GeneralizedBernoulli<T>(k, D) / T(2);
  T fact(1);
  for (int i = 2; i <= k; i++) {
    fact *= T(i);
  }
  q /= fact;
  for (int i = 0; i < k; i++) {
    q *= T(2);
    q /= f_T;
  }
  if ((1 + (k - delta) / 2) % 2 == 1) {
    q = -q;
  }
  return {q, T(k), f_T};
}

// clang-format off
#endif  // SRC_GENERA_SIEGELSPECIALVALUES_H_
// clang-format on
