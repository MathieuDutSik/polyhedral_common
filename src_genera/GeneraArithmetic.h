// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENERAARITHMETIC_H_
#define SRC_GENERA_GENERAARITHMETIC_H_

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryModular.h"
#include "ResidueArithmetic.h"
#include "Basic_functions.h"
#include "factorizations.h"
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>
// clang-format on

/*
  The arithmetic of integers and rationals used across src_genera: primes,
  valuations, residues, quadratic symbols and discriminants, binomials.
  Nothing here takes or returns a matrix.
 */

// A 64-bit integer as an element of T. int64_t is long long on some
// platforms, for which the GMP and FLINT classes have no constructor; long
// is 64 bits wide on every LP64 platform this code targets.
template <typename T> T ScalarFromInt64(int64_t x) {
  static_assert(sizeof(long) == sizeof(int64_t), "long must be 64 bits");
  return T(static_cast<long>(x));
}

// The rational a / b in the field type T.
template <typename T> T RationalFraction(int64_t a, int64_t b) {
  T x = UniversalScalarConversion<T, int64_t>(a);
  x /= UniversalScalarConversion<T, int64_t>(b);
  return x;
}

// v_p(x) for a nonzero integer x.
template <typename Tint> int IntegerValuation(Tint const &x, int64_t p) {
  return static_cast<int>(
      prime_valuation(x, UniversalScalarConversion<Tint, int64_t>(p)));
}

// v_p(x) for a nonzero rational x, of either sign.
template <typename T> int RationalValuation(T const &x, int64_t p) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  return rational_prime_valuation(x, UniversalScalarConversion<Tint, int64_t>(p));
}

inline int NextPrimeAfter(int q) {
  int cand = q + 1;
  while (true) {
    bool is_prime = cand >= 2;
    for (int d = 2; d * d <= cand; d++) {
      if (cand % d == 0) {
        is_prime = false;
        break;
      }
    }
    if (is_prime) {
      return cand;
    }
    cand++;
  }
}

template <typename Tint> int ResidueModP(Tint const &x, int p) {
  Tint p_T(p);
  Tint r = x % p_T;
  int r_int = UniversalScalarConversion<int, Tint>(r);
  return ((r_int % p) + p) % p;
}

// The odd prime factors of the positive integer N, without multiplicity.
template <typename Tint> std::vector<int64_t> OddPrimeDivisors(Tint const &N) {
  std::set<int64_t> primes;
  if (N > Tint(1)) {
    for (auto const &q : FactorsInt(N)) {
      int64_t q_i = UniversalScalarConversion<int64_t, Tint>(q);
      if (q_i != 2) {
        primes.insert(q_i);
      }
    }
  }
  return std::vector<int64_t>(primes.begin(), primes.end());
}

// The prime and the exponent of the prime power q (1 counts as 2^0).
inline std::pair<int64_t, int> PrimePowerDecomposition(int64_t q) {
  if (q < 1) {
    std::cerr << "GSYM: the scale " << q << " is not a positive integer\n";
    throw TerminalException{1};
  }
  if (q == 1) {
    return {2, 0};
  }
  int64_t p = 2;
  while (q % p != 0) {
    p++;
  }
  int k = 0;
  while (q % p == 0) {
    q /= p;
    k++;
  }
  if (q != 1) {
    std::cerr << "GSYM: a scale is not a prime power\n";
    throw TerminalException{1};
  }
  return {p, k};
}

inline int64_t LeastNonResidue(int64_t p) {
  int64_t u = 2;
  while (kronecker_symbol<int64_t>(u, p) != -1) {
    u++;
  }
  return u;
}

inline std::vector<int64_t> SmallPrimeDivisors(int64_t x) {
  std::vector<int64_t> res;
  for (int64_t q = 2; q * q <= x; q++) {
    if (x % q == 0) {
      res.push_back(q);
      while (x % q == 0) {
        x /= q;
      }
    }
  }
  if (x > 1) {
    res.push_back(x);
  }
  return res;
}

inline bool IsSquareFree(int64_t x) {
  if (x < 0) {
    x = -x;
  }
  for (int64_t q = 2; q * q <= x; q++) {
    if (x % (q * q) == 0) {
      return false;
    }
  }
  return true;
}

// A square root of the square a modulo the odd prime p (Tonelli-Shanks).
inline int64_t SqrtModPrime(int64_t a, int64_t p) {
  a = ResMod(a, p);
  if (a == 0) {
    return 0;
  }
  int64_t q = p - 1;
  int s = 0;
  while (q % 2 == 0) {
    q /= 2;
    s++;
  }
  int64_t z = LeastNonResidue(p);
  int64_t m = s;
  int64_t c = ResPow(z, q, p);
  int64_t t = ResPow(a, q, p);
  int64_t r = ResPow(a, (q + 1) / 2, p);
  while (t != 1) {
    int64_t i = 0, t2 = t;
    while (t2 != 1) {
      t2 = ResMul(t2, t2, p);
      i++;
    }
    int64_t b = c;
    for (int64_t j = 0; j < m - i - 1; j++) {
      b = ResMul(b, b, p);
    }
    m = i;
    c = ResMul(b, b, p);
    t = ResMul(t, c, p);
    r = ResMul(r, b, p);
  }
  return r;
}

// The integer residue of the p-integral rational x modulo M.
template <typename T, typename Tint> int64_t ResidueModInt(T const &x, int64_t M) {
  Tint M_T(static_cast<long>(M));
  Tint num = UniversalScalarConversion<Tint, T>(GetNumerator(x));
  Tint den = UniversalScalarConversion<Tint, T>(GetDenominator(x));
  int64_t num_r = UniversalScalarConversion<int64_t, Tint>(ResInt(num, M_T));
  int64_t den_r = UniversalScalarConversion<int64_t, Tint>(ResInt(den, M_T));
  return ResMul(num_r, ResInv(den_r, M), M);
}

inline int64_t EulerPhi(int64_t d) {
  int64_t res = 0;
  for (int64_t a = 1; a <= d; a++) {
    if (std::gcd(a, d) == 1) {
      res++;
    }
  }
  return res;
}

// Whether d (normalized) is a power of a prime: E/K is ramified there.
inline bool IsPrimePower(int64_t d) {
  if (d < 2) {
    return false;
  }
  int64_t p = 2;
  while (d % p != 0) {
    p++;
  }
  while (d % p == 0) {
    d /= p;
  }
  return d == 1;
}

// The order of q modulo d'.
inline int MultiplicativeOrder(int64_t q, int64_t m) {
  if (m == 1) {
    return 1;
  }
  int64_t x = q % m;
  int o = 1;
  while (x != 1) {
    x = (x * q) % m;
    o++;
  }
  return o;
}

template <typename Tint>
inline int LegendreModQ(Tint const &a, int64_t q) {
  return jacobi_symbol(a, UniversalScalarConversion<Tint, int64_t>(q));
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

inline bool IsFundamentalDiscriminant(int64_t D) {
  int64_t r = ((D % 4) + 4) % 4;
  if (r == 1) {
    return IsSquareFree(D);
  }
  if (r == 0) {
    int64_t k = D / 4;
    int64_t rk = ((k % 4) + 4) % 4;
    return (rk == 2 || rk == 3) && IsSquareFree(k);
  }
  return false;
}

// The residue in [0, mod) of the integer val, reduced exactly before the
// conversion to 64 bits (val itself may not fit: the images of a vector by
// a local isometry known modulo q^c have large coordinates)
template <typename T> int64_t ResidueOfInteger(T const &val, int64_t mod) {
  T mod_T = ScalarFromInt64<T>(mod);
  T fl = UniversalFloorScalarInteger<T, T>(val / mod_T);
  T r = val - mod_T * fl;
  return UniversalScalarConversion<int64_t, T>(r);
}

// The residue in (-mod/2, mod/2] of val: an entry smaller than mod / 2 is
// kept as it is
template <typename T>
int64_t CenteredResidueOfInteger(T const &val, int64_t mod) {
  int64_t r = ResidueOfInteger(val, mod);
  return (r > mod / 2) ? r - mod : r;
}

inline int64_t CenteredModP(int64_t x, int64_t mod) {
  int64_t r = x % mod;
  if (r < 0) {
    r += mod;
  }
  return (r > mod / 2) ? r - mod : r;
}

/*
  The square roots in the totally real K = Q(theta): d = y^2 forces y in
  den^{-1} O_K, O_K in disc(f)^{-1} Z[theta], and y(theta_j) = +/- the
  square root of d(theta_j) at the real roots theta_j of f. The roots and
  the square roots are approximated by rationals (PolynomialRealRoots.h),
  the Vandermonde system is solved exactly on those approximations, and
  for each sign pattern the coordinates of y are rounded and the candidate
  is checked exactly; so a square root returned is exact, and none is
  missed as the precision bounds the rounding error.
 */
template <typename Tint> size_t BitLength(Tint x) {
  if (x < 0) {
    x = -x;
  }
  size_t n = 0;
  Tint two(2);
  while (x > 0) {
    x = QuoInt(x, two);
    n++;
  }
  return n;
}

// The floor of the square root of a nonnegative integer.
template <typename Tint> Tint IntegerSquareRoot(Tint const &a) {
  Tint lo(0), hi(1);
  while (hi * hi <= a) {
    hi *= 2;
  }
  // lo^2 <= a < hi^2
  while (hi - lo > 1) {
    Tint mid = (lo + hi) / 2;
    if (mid * mid <= a) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

// The Gaussian binomial coefficient [a b]_p.
template <typename T> T GaussianBinomial(int a, int b, int64_t p) {
  if (b < 0 || b > a) {
    return T(0);
  }
  T p_T = ScalarFromInt64<T>(p);
  T res(1);
  for (int i = 0; i < b; i++) {
    res *= MyPowSigned(p_T, a - i) - T(1);
    res /= MyPowSigned(p_T, i + 1) - T(1);
  }
  return res;
}

inline int FloorDiv2(int a) { return (a >= 0) ? a / 2 : -((-a + 1) / 2); }

inline int MinusOnePower(int e) { return (ResMod(e, 2) == 0) ? 1 : -1; }
/*
  The same key for an integral (even at p = 2) Gram matrix in int64, by a
  Jordan splitting modulo p^K, K = v_p(det) + 5: pivoting on an entry of
  minimal valuation e, every update subtracts a product of a quotient known
  modulo p^{K-e} by an entry of valuation >= e, so the entries stay exact
  modulo p^K, which determines the units modulo 8 (or p) of every Jordan
  block. Returns nothing when p^{2K} would overflow; the caller then uses
  LocalSymbolKey.
 */
inline int64_t ValuationResidue(int64_t x, int64_t p, int K) {
  if (x == 0) {
    return K;
  }
  int v = 0;
  while (x % p == 0) {
    x /= p;
    v++;
  }
  return v;
}

// clang-format off
#endif  // SRC_GENERA_GENERAARITHMETIC_H_
// clang-format on
