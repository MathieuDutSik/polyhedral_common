// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENERAARITHMETIC_H_
#define SRC_GENERA_GENERAARITHMETIC_H_

// clang-format off
#include "ArithmeticFunctions.h"
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

// The value of a character, or a small signed key standing for one (a
// Legendre symbol, a unit residue modulo 8).
using T_character = int64_t;

// The Legendre (Kronecker) symbol of a residue a modulo the prime p.
inline int KroneckerResidue(T_residue a, T_small_prime p) {
  return kronecker_symbol<int64_t>(static_cast<int64_t>(a),
                                   static_cast<int64_t>(p));
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

// The integer residue of the p-integral rational x modulo M.
template <typename T, typename Tint>
T_residue ResidueModInt(T const &x, T_residue M) {
  Tint num = UniversalScalarConversion<Tint, T>(GetNumerator(x));
  Tint den = UniversalScalarConversion<Tint, T>(GetDenominator(x));
  return ResMul(ResModOf(num, M), ResInv(ResModOf(den, M), M), M);
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

// Whether t is a power of x modulo m (x prime to m; always for m <= 2).
inline bool IsPowerModulo(int64_t t, int64_t x, int64_t m) {
  if (m <= 2) {
    return true;
  }
  int64_t target = ((t % m) + m) % m;
  int64_t y = 1;
  int ord = MultiplicativeOrder(x, m);
  for (int k = 0; k < ord; k++) {
    if (y == target) {
      return true;
    }
    y = (y * (x % m)) % m;
  }
  return false;
}

// The Ramanujan sum c_n(j) = Tr_{Q(zeta_n) / Q}(zeta_n^j) = sum_{e | gcd(n,
// j)} mu(n / e) e.
inline int64_t RamanujanSum(int64_t n, int64_t j) {
  int64_t g = std::gcd(n, ((j % n) + n) % n);
  if (g == 0) {
    g = n;
  }
  int64_t res = 0;
  for (int64_t e = 1; e <= g; e++) {
    if (g % e == 0) {
      res += Moebius(n / e) * e;
    }
  }
  return res;
}

template <typename Tint>
inline int LegendreModQ(Tint const &a, int64_t q) {
  return jacobi_symbol(a, UniversalScalarConversion<Tint, int64_t>(q));
}

// (-1)^r x.
inline int64_t SignedPower(int r, int64_t x) { return (r % 2 == 0) ? x : -x; }

// The discriminant of Q(sqrt(x)) for x non-zero; 1 when x is a square.
inline int64_t FundamentalDiscriminant(int64_t x) {
  int64_t d = SquareFreePart(x);
  if (d == 1) {
    return 1;
  }
  int64_t r = ((d % 4) + 4) % 4;
  return (r == 1) ? d : 4 * d;
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

inline int FloorDiv2(int a) { return (a >= 0) ? a / 2 : -((-a + 1) / 2); }

inline int MinusOnePower(int e) { return (e % 2 == 0) ? 1 : -1; }
/*
  The same key for an integral (even at p = 2) Gram matrix in int64, by a
  Jordan splitting modulo p^K, K = v_p(det) + 5: pivoting on an entry of
  minimal valuation e, every update subtracts a product of a quotient known
  modulo p^{K-e} by an entry of valuation >= e, so the entries stay exact
  modulo p^K, which determines the units modulo 8 (or p) of every Jordan
  block. Returns nothing when p^{2K} would overflow; the caller then uses
  LocalSymbolKey.
 */
inline int ValuationResidue(T_residue x, T_small_prime p, int K) {
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
