// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_SIEGELSERIES_H_
#define SRC_GENERA_SIEGELSERIES_H_

// clang-format off
#include "MAT_Matrix.h"
#include "SiegelSpecialValues.h"
#include <cstdint>
#include <map>
#include <tuple>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_SIEGEL_SERIES
#endif

/*
  The local Siegel series of a half-integral matrix B over Z_p, and from it
  the local density alpha_p(S, T) when T is unimodular at p.

  Katsurada (H. Katsurada, An explicit formula for Siegel series, Amer. J.
  Math. 121 (1999), 415-452) writes the Siegel series as
  b_p(B, s) = gamma_p(B, X) F_p(B, X) with X = p^{-s}, F_p a polynomial,
  and gives a recursion for F_p along a Jordan decomposition of B
  (Theorems 4.1 and 4.2). The code below follows the Sage package
  siegel_series of S. Takemori (GPL-3, github.com/stakemori/siegel_series),
  which implements that recursion, translated to exact C++ arithmetic: the
  recursion is carried out on polynomials and rational functions in X, never
  at a point, since its coefficients have poles that only cancel in the sum.

  For T even unimodular at p of rank n = 2k (at p = 2 one of the two even
  unimodular classes), with S = 2B of rank m, the density in the convention
  of LocalDensity.h is the Siegel series at X = xi_T p^{-k}:

      alpha_p(S, T) = 2^{(m-1)[p=2]} gamma_p(B, X) F_p(B, X),
      gamma_p(B, X) = (1 - X) prod_{i=1}^{[m/2]} (1 - p^{2i} X^2)
                      / (1 - p^{m/2} xi_B X)          (last factor: m even),

  xi_T = kronecker((-1)^k det T, p), xi_B = kronecker(D_B, p) with D_B the
  fundamental discriminant of (-1)^{m/2} det S. This identity was pinned in
  the Sage prototype on exact representation numbers of E8 and of the
  rank-16 even unimodular genus. At p = 2 the factor 2^{m-1} is a
  convention: it makes the global constant of LocalDensity.h equal to 2;
  the density of the definition (exhaustive counting modulo 2^k) is the
  value without it.

  All the local invariants of Katsurada's recursion depend only on p-adic
  square classes, carried as (valuation, unit residue), the residue being
  taken modulo 8 at p = 2 and modulo p otherwise.
 */

//
// Univariate polynomials and rational functions over T
//

template <typename T> struct UnivariatePolynomial {
  std::vector<T> coeff;
  void normalize() {
    while (!coeff.empty() && coeff.back() == 0) {
      coeff.pop_back();
    }
  }
  int degree() const { return static_cast<int>(coeff.size()) - 1; }
  bool is_zero() const { return coeff.empty(); }
};

template <typename T> UnivariatePolynomial<T> PolyConstant(T const &c) {
  UnivariatePolynomial<T> P{{c}};
  P.normalize();
  return P;
}

// c * X^d
template <typename T> UnivariatePolynomial<T> PolyMonomial(T const &c, int d) {
  UnivariatePolynomial<T> P;
  if (c == 0) {
    return P;
  }
  P.coeff.assign(d + 1, T(0));
  P.coeff[d] = c;
  return P;
}

template <typename T>
UnivariatePolynomial<T> operator+(UnivariatePolynomial<T> const &A,
                                  UnivariatePolynomial<T> const &B) {
  UnivariatePolynomial<T> C;
  size_t len = std::max(A.coeff.size(), B.coeff.size());
  C.coeff.assign(len, T(0));
  for (size_t i = 0; i < A.coeff.size(); i++) {
    C.coeff[i] += A.coeff[i];
  }
  for (size_t i = 0; i < B.coeff.size(); i++) {
    C.coeff[i] += B.coeff[i];
  }
  C.normalize();
  return C;
}

template <typename T>
UnivariatePolynomial<T> operator-(UnivariatePolynomial<T> const &A,
                                  UnivariatePolynomial<T> const &B) {
  UnivariatePolynomial<T> C;
  size_t len = std::max(A.coeff.size(), B.coeff.size());
  C.coeff.assign(len, T(0));
  for (size_t i = 0; i < A.coeff.size(); i++) {
    C.coeff[i] += A.coeff[i];
  }
  for (size_t i = 0; i < B.coeff.size(); i++) {
    C.coeff[i] -= B.coeff[i];
  }
  C.normalize();
  return C;
}

template <typename T>
UnivariatePolynomial<T> operator*(UnivariatePolynomial<T> const &A,
                                  UnivariatePolynomial<T> const &B) {
  UnivariatePolynomial<T> C;
  if (A.is_zero() || B.is_zero()) {
    return C;
  }
  C.coeff.assign(A.coeff.size() + B.coeff.size() - 1, T(0));
  for (size_t i = 0; i < A.coeff.size(); i++) {
    for (size_t j = 0; j < B.coeff.size(); j++) {
      C.coeff[i + j] += A.coeff[i] * B.coeff[j];
    }
  }
  C.normalize();
  return C;
}

// P(a X)
template <typename T>
UnivariatePolynomial<T> PolyScaleVariable(UnivariatePolynomial<T> const &P,
                                          T const &a) {
  UnivariatePolynomial<T> Q = P;
  T pw(1);
  for (auto &c : Q.coeff) {
    c *= pw;
    pw *= a;
  }
  Q.normalize();
  return Q;
}

template <typename T>
T PolyEvaluate(UnivariatePolynomial<T> const &P, T const &x) {
  T val(0);
  for (size_t i = P.coeff.size(); i > 0; i--) {
    val = val * x + P.coeff[i - 1];
  }
  return val;
}

// A / B for a B that divides A exactly; the remainder is checked to vanish.
template <typename T>
UnivariatePolynomial<T> PolyExactDivision(UnivariatePolynomial<T> const &A,
                                          UnivariatePolynomial<T> const &B) {
  if (B.is_zero()) {
    std::cerr << "SSER: polynomial division by zero\n";
    throw TerminalException{1};
  }
  UnivariatePolynomial<T> R = A;
  UnivariatePolynomial<T> Q;
  int dB = B.degree();
  if (R.degree() >= dB) {
    Q.coeff.assign(R.degree() - dB + 1, T(0));
  }
  T lead = B.coeff[dB];
  while (!R.is_zero() && R.degree() >= dB) {
    int shift = R.degree() - dB;
    T c = R.coeff[R.degree()] / lead;
    Q.coeff[shift] = c;
    for (int i = 0; i <= dB; i++) {
      R.coeff[shift + i] -= c * B.coeff[i];
    }
    R.normalize();
  }
  if (!R.is_zero()) {
    std::cerr << "SSER: a division that should be exact left a remainder\n";
    throw TerminalException{1};
  }
  Q.normalize();
  return Q;
}

template <typename T> struct RationalFunction {
  UnivariatePolynomial<T> num;
  UnivariatePolynomial<T> den;
};

template <typename T>
RationalFunction<T> RatFromPoly(UnivariatePolynomial<T> const &P) {
  return {P, PolyConstant<T>(T(1))};
}

template <typename T>
RationalFunction<T> operator*(RationalFunction<T> const &A,
                              RationalFunction<T> const &B) {
  return {A.num * B.num, A.den * B.den};
}

template <typename T>
RationalFunction<T> operator+(RationalFunction<T> const &A,
                              RationalFunction<T> const &B) {
  return {A.num * B.den + B.num * A.den, A.den * B.den};
}

template <typename T>
RationalFunction<T> RatScaleVariable(RationalFunction<T> const &A,
                                     T const &a) {
  return {PolyScaleVariable(A.num, a), PolyScaleVariable(A.den, a)};
}

template <typename T>
UnivariatePolynomial<T> RatToPolynomial(RationalFunction<T> const &A) {
  return PolyExactDivision(A.num, A.den);
}

//
// p-adic square classes
//

inline int64_t ModPositive(int64_t a, int64_t M) { return ((a % M) + M) % M; }

inline int64_t InverseModInt64(int64_t a, int64_t M) {
  int64_t r0 = M, r1 = ModPositive(a, M), s0 = 0, s1 = 1;
  while (r1 != 0) {
    int64_t q = r0 / r1;
    int64_t r2 = r0 - q * r1;
    int64_t s2 = s0 - q * s1;
    r0 = r1;
    r1 = r2;
    s0 = s1;
    s1 = s2;
  }
  if (r0 != 1) {
    std::cerr << "SSER: " << a << " is not invertible modulo " << M << "\n";
    throw TerminalException{1};
  }
  return ModPositive(s0, M);
}

inline int64_t MulMod(int64_t a, int64_t b, int64_t M) {
  __int128 prod = static_cast<__int128>(ModPositive(a, M)) *
                  static_cast<__int128>(ModPositive(b, M));
  return static_cast<int64_t>(prod % M);
}

// The modulus of the unit residues at p.
inline int64_t UnitModulus(int64_t p) { return (p == 2) ? 8 : p; }

// A non-zero element of Q_p up to the information Hilbert symbols and
// characters need: p^val * u, u a unit known modulo UnitModulus(p).
struct PadicClass {
  int val;
  int64_t unit;
};

inline PadicClass PadicMul(PadicClass const &x, PadicClass const &y,
                           int64_t p) {
  return {x.val + y.val, MulMod(x.unit, y.unit, UnitModulus(p))};
}

inline PadicClass PadicFromInt64(int64_t x, int64_t p) {
  if (x == 0) {
    std::cerr << "SSER: PadicFromInt64 called on zero\n";
    throw TerminalException{1};
  }
  int val = 0;
  while (x % p == 0) {
    x /= p;
    val++;
  }
  return {val, ModPositive(x, UnitModulus(p))};
}

// The class of a non-zero rational, from its numerator and denominator.
template <typename T, typename Tint>
PadicClass PadicFromRational(T const &x, int64_t p) {
  if (x == 0) {
    std::cerr << "SSER: PadicFromRational called on zero\n";
    throw TerminalException{1};
  }
  Tint num = UniversalScalarConversion<Tint, T>(GetNumerator(x));
  Tint den = UniversalScalarConversion<Tint, T>(GetDenominator(x));
  Tint p_T(static_cast<long>(p));
  int val = 0;
  while (ResInt(num, p_T) == 0) {
    num /= p_T;
    val++;
  }
  while (ResInt(den, p_T) == 0) {
    den /= p_T;
    val--;
  }
  int64_t M = UnitModulus(p);
  Tint M_T(static_cast<long>(M));
  int64_t num_r =
      UniversalScalarConversion<int64_t, Tint>(ResInt(num, M_T));
  int64_t den_r =
      UniversalScalarConversion<int64_t, Tint>(ResInt(den, M_T));
  return {val, MulMod(num_r, InverseModInt64(den_r, M), M)};
}

// The Legendre symbol of a unit residue modulo the odd prime p.
inline int LegendreUnit(int64_t u, int64_t p) { return KroneckerSymbol(u, p); }

// The Hilbert symbol (x, y)_p.
inline int HilbertSymbolClass(PadicClass const &x, PadicClass const &y,
                              int64_t p) {
  int64_t a = x.val, b = y.val;
  if (p != 2) {
    int sign = 1;
    if ((a * b) % 2 != 0 && ((p - 1) / 2) % 2 != 0) {
      sign = -1;
    }
    if (b % 2 != 0) {
      sign *= LegendreUnit(x.unit, p);
    }
    if (a % 2 != 0) {
      sign *= LegendreUnit(y.unit, p);
    }
    return sign;
  }
  int64_t u = x.unit, v = y.unit;
  auto eps = [](int64_t w) -> int64_t { return ((w - 1) / 2) % 2; };
  auto omega = [](int64_t w) -> int64_t { return ((w * w - 1) / 8) % 2; };
  int64_t e = eps(u) * eps(v) + a * omega(v) + b * omega(u);
  return (ModPositive(e, 2) == 0) ? 1 : -1;
}

// The quadratic character chi_p of Katsurada: 0 when Q_p(sqrt x) / Q_p is
// ramified, +1 when x is a square, -1 when it is unramified non-trivial.
inline int ChiPadic(PadicClass const &x, int64_t p) {
  if (x.val % 2 != 0) {
    return 0;
  }
  if (p != 2) {
    return LegendreUnit(x.unit, p);
  }
  if (x.unit == 1) {
    return 1;
  }
  if (x.unit == 5) {
    return -1;
  }
  return 0;
}

//
// Jordan blocks
//

enum class SiegelBlockKind { Unit = 0, H = 1, Y = 2 };

// p^e * (unit) for a Unit block, 2^e times the binary form h = [[0, 1/2],
// [1/2, 0]] or y = [[1, 1/2], [1/2, 1]] for an H or Y block (p = 2 only).
struct SiegelJordanBlock {
  int e;
  SiegelBlockKind kind;
  int64_t unit;
};

struct SiegelJordanBlocks {
  int64_t p;
  std::vector<SiegelJordanBlock> blocks;
  int dim() const {
    int d = 0;
    for (auto const &b : blocks) {
      d += (b.kind == SiegelBlockKind::Unit) ? 1 : 2;
    }
    return d;
  }
};

inline std::tuple<int64_t, std::vector<std::tuple<int, int, int64_t>>>
SiegelBlocksKey(SiegelJordanBlocks const &q) {
  std::vector<std::tuple<int, int, int64_t>> v;
  for (auto const &b : q.blocks) {
    v.emplace_back(b.e, static_cast<int>(b.kind), b.unit);
  }
  return {q.p, v};
}

// The rational diagonal entries of the form the blocks describe (for the
// Hasse invariant), as p-adic classes.
inline std::vector<PadicClass> RationalDiagonal(SiegelJordanBlocks const &q) {
  std::vector<PadicClass> diag;
  for (auto const &b : q.blocks) {
    if (b.kind == SiegelBlockKind::Unit) {
      diag.push_back({b.e, b.unit});
    } else if (b.kind == SiegelBlockKind::H) {
      diag.push_back({b.e, 1});
      diag.push_back({b.e, 7});
    } else {
      diag.push_back({b.e, 1});
      diag.push_back({b.e, 3});
    }
  }
  return diag;
}

inline PadicClass GramDetClass(SiegelJordanBlocks const &q) {
  int64_t p = q.p;
  PadicClass det = {0, 1};
  for (auto const &b : q.blocks) {
    if (b.kind == SiegelBlockKind::Unit) {
      det = PadicMul(det, {b.e, b.unit}, p);
    } else if (b.kind == SiegelBlockKind::H) {
      // 2^{2e} * (-1/4)
      det = PadicMul(det, {2 * b.e - 2, 7}, p);
    } else {
      // 2^{2e} * (3/4)
      det = PadicMul(det, {2 * b.e - 2, 3}, p);
    }
  }
  return det;
}

inline int ContentOrder(SiegelJordanBlocks const &q) {
  int c = q.blocks[0].e;
  for (auto const &b : q.blocks) {
    c = std::min(c, b.e);
  }
  return c;
}

inline int HasseInvariantOMeara(SiegelJordanBlocks const &q) {
  int64_t p = q.p;
  std::vector<PadicClass> diag = RationalDiagonal(q);
  int res = 1;
  for (size_t i = 0; i < diag.size(); i++) {
    for (size_t j = i + 1; j < diag.size(); j++) {
      res *= HilbertSymbolClass(diag[i], diag[j], p);
    }
  }
  PadicClass minus_one = {0, ModPositive(-1, UnitModulus(p))};
  return res * HilbertSymbolClass(GramDetClass(q), minus_one, p);
}

inline PadicClass PadicSignedPower(int r, PadicClass const &x, int64_t p) {
  if (r % 2 == 0) {
    return x;
  }
  return {x.val, MulMod(x.unit, -1, UnitModulus(p))};
}

// The invariants of Katsurada, Section 3.
inline int XiInvariant(SiegelJordanBlocks const &q) {
  return ChiPadic(PadicSignedPower(q.dim() / 2, GramDetClass(q), q.p), q.p);
}

inline int EtaInvariant(SiegelJordanBlocks const &q) {
  int64_t p = q.p;
  int n = q.dim();
  PadicClass det = GramDetClass(q);
  int h1 = HilbertSymbolClass(det, PadicSignedPower((n - 1) / 2, det, p), p);
  PadicClass minus_one = {0, ModPositive(-1, UnitModulus(p))};
  int h2 = HilbertSymbolClass(minus_one, minus_one, p);
  int exp_h2 = ((n * n - 1) / 8) % 2;
  return HasseInvariantOMeara(q) * h1 * ((exp_h2 == 0) ? 1 : h2);
}

inline int SmallD(SiegelJordanBlocks const &q) {
  int v = GramDetClass(q).val;
  if (q.p == 2) {
    v += 2 * (q.dim() / 2);
  }
  return v;
}

inline int FloorDiv2(int a) { return (a >= 0) ? a / 2 : -((-a + 1) / 2); }

inline int DeltaInvariant(SiegelJordanBlocks const &q) {
  int db = SmallD(q);
  if (q.dim() % 2 == 0) {
    int dl = (q.p == 2) ? 1 : 0;
    return 2 * FloorDiv2(db + 1 - dl);
  }
  return db;
}

inline int XiDash(int xi) { return 1 + xi - xi * xi; }

inline int MinusOnePower(int e) { return (ModPositive(e, 2) == 0) ? 1 : -1; }

//
// Jordan decomposition over Z_p of a half-integral matrix
//

template <typename T>
MyMatrix<T> BracketAction(MyMatrix<T> const &A, MyMatrix<T> const &B) {
  return B.transpose() * A * B;
}

/*
  The entry of minimal order, the order of an off-diagonal entry being that
  of twice it. Among the entries of minimal order, a diagonal one is
  preferred for odd p and an off-diagonal one for p = 2; otherwise the first
  in the order (0,0), (0,1), ..., (1,1), ...
 */
template <typename T, typename Tint>
std::pair<int, int> FindMinOrderEntry(MyMatrix<T> const &S, int64_t p) {
  int n = S.rows();
  std::vector<std::tuple<int, int, int>> list;
  int min_val = 0;
  bool first = true;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      if (S(i, j) == 0) {
        continue;
      }
      T val = (i == j) ? S(i, j) : T(2) * S(i, j);
      int v = PadicFromRational<T, Tint>(val, p).val;
      list.emplace_back(v, i, j);
      if (first || v < min_val) {
        min_val = v;
        first = false;
      }
    }
  }
  if (first) {
    std::cerr << "SSER: Jordan decomposition of a zero matrix\n";
    throw TerminalException{1};
  }
  bool want_diag = (p != 2);
  for (auto const &[v, i, j] : list) {
    if (v == min_val && ((i == j) == want_diag)) {
      return {i, j};
    }
  }
  for (auto const &[v, i, j] : list) {
    if (v == min_val) {
      return {i, j};
    }
  }
  return {0, 0};
}

// Moves the (i0, i0) entry to (0, 0) and clears the first row and column.
template <typename T>
MyMatrix<T> JordanDiagonalStep(int i0, MyMatrix<T> const &Sin) {
  int n = Sin.rows();
  MyMatrix<T> P = IdentityMat<T>(n);
  P(0, 0) = 0;
  P(i0, i0) = 0;
  P(0, i0) = 1;
  P(i0, 0) = 1;
  if (i0 == 0) {
    P = IdentityMat<T>(n);
  }
  MyMatrix<T> S = BracketAction(Sin, P);
  MyMatrix<T> U = IdentityMat<T>(n);
  T a = S(0, 0);
  for (int j = 1; j < n; j++) {
    U(0, j) = -S(0, j) / a;
  }
  return BracketAction(S, U);
}

template <typename T>
MyMatrix<T> SubMatrixFrom(MyMatrix<T> const &S, int start) {
  int n = S.rows() - start;
  MyMatrix<T> R(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      R(i, j) = S(start + i, start + j);
    }
  }
  return R;
}

// p odd: diagonal entries p^{n_1} u_1, ... with n_1 <= n_2 <= ...
template <typename T, typename Tint>
std::vector<T> JordanDiagonalOdd(MyMatrix<T> S, int64_t p) {
  std::vector<T> acc;
  int n = S.rows();
  while (true) {
    if (n == 1) {
      acc.push_back(S(0, 0));
      return acc;
    }
    auto [i0, j0] = FindMinOrderEntry<T, Tint>(S, p);
    if (i0 != j0) {
      MyMatrix<T> U = IdentityMat<T>(n);
      U(j0, i0) = 1;
      S = BracketAction(S, U);
    }
    S = JordanDiagonalStep(i0, S);
    acc.push_back(S(0, 0));
    S = SubMatrixFrom(S, 1);
    n--;
  }
}

template <typename T, typename Tint>
SiegelJordanBlocks JordanBlocksOdd(MyMatrix<T> const &B, int64_t p) {
  std::vector<T> diag = JordanDiagonalOdd<T, Tint>(B, p);
  // the least quadratic non-residue modulo p
  int64_t nonres = 2;
  while (LegendreUnit(nonres, p) != -1) {
    nonres++;
  }
  std::vector<SiegelJordanBlock> res;
  size_t i = 0;
  while (i < diag.size()) {
    PadicClass c0 = PadicFromRational<T, Tint>(diag[i], p);
    int e = c0.val;
    int64_t prod = 1;
    size_t j = i;
    std::vector<SiegelJordanBlock> part;
    while (j < diag.size()) {
      PadicClass c = PadicFromRational<T, Tint>(diag[j], p);
      if (c.val != e) {
        break;
      }
      prod = MulMod(prod, c.unit, p);
      part.push_back({e, SiegelBlockKind::Unit, 1});
      j++;
    }
    if (LegendreUnit(prod, p) == -1) {
      part.back().unit = nonres;
    }
    res.insert(res.end(), part.begin(), part.end());
    i = j;
  }
  SiegelJordanBlocks q;
  q.p = p;
  q.blocks.assign(res.rbegin(), res.rend());
  return q;
}

// Whether a x^2 + b y^2 + c z^2 represents zero over Z_2, a, b, c odd.
inline bool DoesTwoAdicallyRepresentZero(int64_t a, int64_t b, int64_t c) {
  for (int64_t x : {0, 1, 4}) {
    for (int64_t y : {0, 1, 4}) {
      for (int64_t z : {0, 1, 4}) {
        bool some_odd = (x % 2 == 1) || (y % 2 == 1) || (z % 2 == 1);
        if (some_odd && ModPositive(a * x + b * y + c * z, 8) == 0) {
          return true;
        }
      }
    }
  }
  return false;
}

// An element of a list of 2-adic unit diagonal entries of one scale, or a
// binary form h / y produced from three of them.
struct TwoAdicItem {
  bool is_binary;
  int64_t unit;
  SiegelBlockKind kind;
};

// Rewrites diag(us) so that at most two unit diagonal entries remain, the
// others being traded for h / y forms of the next scale.
inline std::vector<TwoAdicItem>
TransJordanDec2(std::vector<TwoAdicItem> const &us) {
  if (us.size() <= 2) {
    return us;
  }
  int64_t u1 = us[0].unit, u2 = us[1].unit, u3 = us[2].unit;
  int64_t u = ModPositive(u1 * u2 * u3, 8);
  std::vector<TwoAdicItem> rest(us.begin() + 3, us.end());
  std::vector<TwoAdicItem> l = TransJordanDec2(rest);
  if (DoesTwoAdicallyRepresentZero(u1, u2, u3)) {
    l.push_back({false, ModPositive(-u, 8), SiegelBlockKind::Unit});
    l.push_back({true, 0, SiegelBlockKind::H});
  } else {
    l.push_back({false, ModPositive(3 * u, 8), SiegelBlockKind::Unit});
    l.push_back({true, 0, SiegelBlockKind::Y});
  }
  std::vector<TwoAdicItem> l_diag, l_non_diag;
  for (auto const &it : l) {
    if (it.is_binary) {
      l_non_diag.push_back(it);
    } else {
      l_diag.push_back(it);
    }
  }
  if (l_diag.size() <= 2) {
    return l;
  }
  std::vector<TwoAdicItem> l1 = TransJordanDec2(l_diag);
  l1.insert(l1.end(), l_non_diag.begin(), l_non_diag.end());
  return l1;
}

template <typename T, typename Tint>
SiegelJordanBlocks JordanBlocks2(MyMatrix<T> S) {
  int64_t p = 2;
  int n = S.rows();
  std::vector<MyMatrix<T>> mat_lst;
  while (true) {
    auto [i0, j0] = FindMinOrderEntry<T, Tint>(S, p);
    if (n == 1 || (n == 2 && i0 != j0)) {
      mat_lst.push_back(S);
      break;
    }
    if (i0 == j0) {
      S = JordanDiagonalStep(i0, S);
      MyMatrix<T> M(1, 1);
      M(0, 0) = S(0, 0);
      mat_lst.push_back(M);
      S = SubMatrixFrom(S, 1);
      n--;
    } else {
      // (i0, j0) -> (0, 1)
      MyMatrix<T> P = ZeroMatrix<T>(n, n);
      std::vector<int> order{i0, j0};
      for (int a = 0; a < n; a++) {
        if (a != i0 && a != j0) {
          order.push_back(a);
        }
      }
      for (int col = 0; col < n; col++) {
        P(order[col], col) = 1;
      }
      S = BracketAction(S, P);
      MyMatrix<T> U = IdentityMat<T>(n);
      T s00 = S(0, 0), s01 = S(0, 1), s11 = S(1, 1);
      T d = s00 * s11 - s01 * s01;
      for (int j = 2; j < n; j++) {
        T s0j = S(0, j), s1j = S(1, j);
        U(0, j) = (-s11 * s0j + s01 * s1j) / d;
        U(1, j) = (s01 * s0j - s00 * s1j) / d;
      }
      S = BracketAction(S, U);
      MyMatrix<T> M(2, 2);
      for (int a = 0; a < 2; a++) {
        for (int b = 0; b < 2; b++) {
          M(a, b) = S(a, b);
        }
      }
      mat_lst.push_back(M);
      S = SubMatrixFrom(S, 2);
      n -= 2;
    }
  }
  // the blocks before the trading of triples of units
  std::vector<SiegelJordanBlock> raw;
  for (auto const &M : mat_lst) {
    if (M.rows() == 1) {
      PadicClass c = PadicFromRational<T, Tint>(M(0, 0), p);
      raw.push_back({c.val, SiegelBlockKind::Unit, c.unit});
    } else {
      int e = PadicFromRational<T, Tint>(M(0, 1), p).val + 1;
      T scale(1);
      for (int i = 0; i < e - 1; i++) {
        scale *= T(2);
      }
      for (int i = e - 1; i < 0; i++) {
        scale /= T(2);
      }
      T det = (M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0)) / (scale * scale);
      PadicClass c = PadicFromRational<T, Tint>(det, p);
      bool is_y = (c.val == 0 && c.unit == 3);
      raw.push_back(
          {e, is_y ? SiegelBlockKind::Y : SiegelBlockKind::H, 0});
    }
  }
  // the binary blocks, then each scale of unit diagonal entries reduced to at
  // most two units, the traded triples becoming binary forms one scale up
  std::vector<SiegelJordanBlock> res;
  for (auto const &b : raw) {
    if (b.kind != SiegelBlockKind::Unit) {
      res.push_back(b);
    }
  }
  std::map<int, std::vector<TwoAdicItem>> by_scale;
  for (auto const &b : raw) {
    if (b.kind == SiegelBlockKind::Unit) {
      by_scale[b.e].push_back({false, b.unit, SiegelBlockKind::Unit});
    }
  }
  for (auto const &[a, us] : by_scale) {
    for (auto const &it : TransJordanDec2(us)) {
      if (it.is_binary) {
        res.push_back({a + 1, it.kind, 0});
      } else {
        res.push_back({a, SiegelBlockKind::Unit, it.unit});
      }
    }
  }
  // decreasing scale; within a scale the unit entries first
  std::map<int, std::vector<SiegelJordanBlock>> groups;
  for (auto const &b : res) {
    groups[b.e].push_back(b);
  }
  SiegelJordanBlocks q;
  q.p = 2;
  for (auto it = groups.rbegin(); it != groups.rend(); ++it) {
    for (auto const &b : it->second) {
      if (b.kind == SiegelBlockKind::Unit) {
        q.blocks.push_back(b);
      }
    }
    for (auto const &b : it->second) {
      if (b.kind != SiegelBlockKind::Unit) {
        q.blocks.push_back(b);
      }
    }
  }
  return q;
}

//
// Katsurada's recursion
//

template <typename T> T PowerInt(T const &a, int e) {
  T res(1);
  if (e >= 0) {
    for (int i = 0; i < e; i++) {
      res *= a;
    }
  } else {
    for (int i = 0; i < -e; i++) {
      res /= a;
    }
  }
  return res;
}

// 1 - c X^d
template <typename T> UnivariatePolynomial<T> OneMinus(T const &c, int d) {
  return PolyConstant<T>(T(1)) - PolyMonomial<T>(c, d);
}

inline int CheckedHalf(int a) {
  if (a % 2 != 0) {
    std::cerr << "SSER: an exponent of Katsurada's recursion is not "
                 "integral\n";
    throw TerminalException{1};
  }
  return a / 2;
}

inline void CheckNonNegative(int a) {
  if (a < 0) {
    std::cerr << "SSER: a negative power of X in Katsurada's recursion\n";
    throw TerminalException{1};
  }
}

template <typename T>
UnivariatePolynomial<T> SiegelSeriesDim1(SiegelJordanBlocks const &q) {
  T p_T = ScalarFromInt64<T>(q.p);
  UnivariatePolynomial<T> P;
  for (int a = 0; a <= ContentOrder(q); a++) {
    P = P + PolyMonomial<T>(PowerInt(p_T, a), a);
  }
  return P;
}

template <typename T>
RationalFunction<T> SiegelSeriesDim2Kernel(int64_t p, int a, int b) {
  if (b == 0) {
    return RatFromPoly(UnivariatePolynomial<T>());
  }
  a = std::min(a, b - 1);
  T p_T = ScalarFromInt64<T>(p);
  // r1 = (1 - (p^2 X)^{a+1}) / (1 - p^2 X) = sum_{i=0}^{a} (p^2 X)^i
  UnivariatePolynomial<T> r1;
  for (int i = 0; i <= a; i++) {
    r1 = r1 + PolyMonomial<T>(PowerInt(p_T, 2 * i), i);
  }
  UnivariatePolynomial<T> rn2 =
      PolyMonomial<T>(PowerInt(p_T, 3 * b + 1), 2 * b + 1) -
      PolyMonomial<T>(PowerInt(p_T, b + 2 * b - a), 2 * b - a);
  UnivariatePolynomial<T> rd2 =
      PolyMonomial<T>(p_T, 1) - PolyConstant<T>(T(1));
  // (r1 - rn2 / rd2) / (1 - p^3 X^2)
  RationalFunction<T> res = {r1 * rd2 - rn2,
                             rd2 * OneMinus<T>(PowerInt(p_T, 3), 2)};
  return res;
}

template <typename T>
UnivariatePolynomial<T> SiegelSeriesDim2(SiegelJordanBlocks const &q) {
  int64_t p = q.p;
  T p_T = ScalarFromInt64<T>(p);
  int c = ContentOrder(q);
  // det_4 = 4 * Gram_det, and the local class of the fundamental
  // discriminant of -det_4
  PadicClass det4 = PadicMul(GramDetClass(q), PadicFromInt64(4, p), p);
  PadicClass x = {det4.val, MulMod(det4.unit, -1, UnitModulus(p))};
  int val_fd;
  int kron_fd;
  if (p != 2) {
    if (x.val % 2 != 0) {
      val_fd = 1;
      kron_fd = 0;
    } else {
      val_fd = 0;
      kron_fd = LegendreUnit(x.unit, p);
    }
  } else {
    if (x.val % 2 != 0) {
      val_fd = 3;
      kron_fd = 0;
    } else if (x.unit % 4 == 1) {
      val_fd = 0;
      kron_fd = (x.unit == 1) ? 1 : -1;
    } else {
      val_fd = 2;
      kron_fd = 0;
    }
  }
  int f = CheckedHalf(det4.val - val_fd);
  RationalFunction<T> A = SiegelSeriesDim2Kernel<T>(p, c, f + 1);
  RationalFunction<T> B = SiegelSeriesDim2Kernel<T>(p, c, f);
  RationalFunction<T> coefB =
      RatFromPoly(PolyMonomial<T>(T(-kron_fd) * p_T, 1));
  return RatToPolynomial(A + coefB * B);
}

template <typename T>
RationalFunction<T> KatsuradaC1(SiegelJordanBlocks const &q,
                                SiegelJordanBlocks const &q2) {
  T p_T = ScalarFromInt64<T>(q.p);
  int n = q.dim();
  if (n % 2 == 0) {
    int xi = XiInvariant(q);
    return {OneMinus<T>(PowerInt(p_T, n / 2) * T(xi), 1),
            OneMinus<T>(PowerInt(p_T, n + 1), 2)};
  }
  int xi_tilde = (n > 1) ? XiInvariant(q2) : 1;
  return {PolyConstant<T>(T(1)),
          OneMinus<T>(PowerInt(p_T, (n + 1) / 2) * T(xi_tilde), 1)};
}

template <typename T>
RationalFunction<T> KatsuradaC0(SiegelJordanBlocks const &q,
                                SiegelJordanBlocks const &q2) {
  T p_T = ScalarFromInt64<T>(q.p);
  int n = q.dim();
  if (n % 2 == 0) {
    int xi = XiInvariant(q);
    int xi_dash = XiDash(xi);
    int eta_tilde = EtaInvariant(q2);
    int delta = DeltaInvariant(q);
    int delta_tilde = DeltaInvariant(q2);
    int expt = delta - delta_tilde + xi * xi;
    CheckNonNegative(expt);
    T coef = T(MinusOnePower(xi + 1) * xi_dash * eta_tilde) *
             PowerInt(p_T, (n / 2) * expt + CheckedHalf(delta));
    UnivariatePolynomial<T> num =
        OneMinus<T>(PowerInt(p_T, n / 2 + 1) * T(xi), 1) *
        PolyMonomial<T>(coef, expt);
    return {num, OneMinus<T>(PowerInt(p_T, n + 1), 2)};
  }
  int eta = EtaInvariant(q);
  int delta = DeltaInvariant(q);
  int xi_tilde = 1, xi_tilde_dash = 1, delta_tilde = 0;
  if (n > 1) {
    xi_tilde = XiInvariant(q2);
    xi_tilde_dash = XiDash(xi_tilde);
    delta_tilde = DeltaInvariant(q2);
  }
  int expt = delta - delta_tilde + 2 - xi_tilde * xi_tilde;
  CheckNonNegative(expt);
  T coef = T(MinusOnePower(xi_tilde) * xi_tilde_dash * eta) *
           PowerInt(p_T, ((n - 1) / 2) * expt +
                             CheckedHalf(2 * delta - delta_tilde + 2));
  return {PolyMonomial<T>(coef, expt),
          OneMinus<T>(PowerInt(p_T, (n + 1) / 2) * T(xi_tilde), 1)};
}

inline SiegelJordanBlocks BlocksTail(SiegelJordanBlocks const &q,
                                     size_t start) {
  SiegelJordanBlocks q2;
  q2.p = q.p;
  q2.blocks.assign(q.blocks.begin() + start, q.blocks.end());
  return q2;
}

inline SiegelJordanBlocks BlocksPrepend(std::vector<SiegelJordanBlock> const &b,
                                        SiegelJordanBlocks const &q) {
  SiegelJordanBlocks res;
  res.p = q.p;
  res.blocks = b;
  res.blocks.insert(res.blocks.end(), q.blocks.begin(), q.blocks.end());
  return res;
}

// The four coefficients of Katsurada's Theorem 4.2 (p = 2), for the binary
// block b1 (two units, or h / y) in front of q2.
template <typename T>
std::vector<RationalFunction<T>>
KatsuradaTheorem42Coefficients(std::vector<SiegelJordanBlock> const &b1,
                               SiegelJordanBlocks const &q2) {
  bool b1_unit = (b1[0].kind == SiegelBlockKind::Unit);
  int m = b1[0].e;
  int n = q2.dim() + 2;
  T two(2);
  SiegelJordanBlocks q = BlocksPrepend(b1, q2);
  SiegelJordanBlocks q3 =
      BlocksPrepend({{m, SiegelBlockKind::Unit, 1}}, q2);
  int delta = DeltaInvariant(q);
  int delta_tilde = DeltaInvariant(q3);
  int delta_hat = DeltaInvariant(q2);
  int sigma2;  // twice sigma
  if ((n % 2 == 0 && b1_unit && SmallD(q) % 2 == 1) ||
      (n % 2 == 0 && !b1_unit && XiInvariant(q2) == 0)) {
    sigma2 = 2 * delta_tilde - delta - delta_hat + 2;
  } else if (n % 2 == 1 && !b1_unit && SmallD(q3) % 2 == 0) {
    sigma2 = 4;
  } else {
    sigma2 = 0;
  }
  int sigma = CheckedHalf(sigma2);
  RationalFunction<T> c11, c10, c21, c20;
  if (n % 2 == 0) {
    int xi = XiInvariant(q);
    int xi_dash = XiDash(xi);
    int xi_hat = XiInvariant(q2);
    int xi_hat_dash = XiDash(xi_hat);
    int eta_tilde;
    if (b1_unit && SmallD(q2) % 2 == 0) {
      SiegelJordanBlocks qq =
          BlocksPrepend({{m, SiegelBlockKind::Unit, b1[1].unit}}, q2);
      eta_tilde = EtaInvariant(qq);
    } else if (!b1_unit && xi_hat != 0) {
      int sgn = MinusOnePower(((n - 1) * (n - 1) - 1) / 8);
      PadicClass two_m = {m, 1};
      PadicClass det2 = PadicSignedPower(n / 2 - 1, GramDetClass(q2), 2);
      eta_tilde =
          sgn * HasseInvariantOMeara(q2) * HilbertSymbolClass(two_m, det2, 2);
    } else {
      eta_tilde = 1;
    }
    c11 = {OneMinus<T>(PowerInt(two, n / 2) * T(xi), 1),
           OneMinus<T>(PowerInt(two, n + 1), 2)};
    int expt = delta - delta_tilde + xi * xi + sigma;
    CheckNonNegative(expt);
    T coef10 = T(MinusOnePower(xi + 1) * xi_dash * eta_tilde) *
               PowerInt(two, CheckedHalf(delta) + (n / 2) * expt);
    c10 = {OneMinus<T>(PowerInt(two, n / 2 + 1) * T(xi), 1) *
               PolyMonomial<T>(coef10, expt),
           OneMinus<T>(PowerInt(two, n + 1), 2)};
    c21 = {PolyConstant<T>(T(1)),
           OneMinus<T>(PowerInt(two, n / 2) * T(xi_hat), 1)};
    int expt2 = delta_tilde - delta_hat + 2 - xi_hat * xi_hat - sigma;
    CheckNonNegative(expt2);
    T coef20 = T(MinusOnePower(xi_hat) * xi_hat_dash * eta_tilde) *
               PowerInt(two, ((n - 2) / 2) * expt2 +
                                 CheckedHalf(2 * delta_tilde - delta_hat + 2 -
                                             2 * sigma));
    c20 = {PolyMonomial<T>(coef20, expt2),
           OneMinus<T>(PowerInt(two, n / 2) * T(xi_hat), 1)};
  } else {
    int eta = EtaInvariant(q);
    int eta_hat = EtaInvariant(q2);
    int xi_tilde = (!b1_unit && SmallD(q3) % 2 == 0) ? 1 : 0;
    c11 = {PolyConstant<T>(T(1)),
           OneMinus<T>(PowerInt(two, (n + 1) / 2) * T(xi_tilde), 1)};
    int expt = delta - delta_tilde + 2 - xi_tilde * xi_tilde + sigma;
    CheckNonNegative(expt);
    T coef10 = T(MinusOnePower(xi_tilde) * eta) *
               PowerInt(two, ((n - 1) / 2) * expt +
                                 CheckedHalf(2 * delta - delta_tilde + 2 +
                                             sigma));
    c10 = {PolyMonomial<T>(coef10, expt),
           OneMinus<T>(PowerInt(two, (n + 1) / 2) * T(xi_tilde), 1)};
    c21 = {OneMinus<T>(PowerInt(two, (n - 1) / 2) * T(xi_tilde), 1),
           OneMinus<T>(PowerInt(two, n), 2)};
    int expt2 = delta_tilde - delta_hat + xi_tilde * xi_tilde - sigma;
    CheckNonNegative(expt2);
    T coef20 = T(MinusOnePower(xi_tilde + 1) * eta_hat) *
               PowerInt(two, ((n - 1) / 2) * expt2 +
                                 CheckedHalf(delta_tilde - sigma));
    c20 = {OneMinus<T>(PowerInt(two, (n + 1) / 2) * T(xi_tilde), 1) *
               PolyMonomial<T>(coef20, expt2),
           OneMinus<T>(PowerInt(two, n), 2)};
  }
  return {c11, c10, c21, c20};
}

/*
  F_p(B, X) for the Jordan blocks q of B, memoised on the blocks: the
  recursion of a matrix runs through the tails of its own decomposition,
  and a sweep over many matrices meets the same tails again and again.
 */
template <typename T> struct SiegelSeriesCache {
  std::map<std::tuple<int64_t, std::vector<std::tuple<int, int, int64_t>>>,
           UnivariatePolynomial<T>>
      map;
};

template <typename T>
UnivariatePolynomial<T> SiegelSeriesPolynomialBlocks(SiegelJordanBlocks const &q,
                                                     SiegelSeriesCache<T> &cache) {
  auto key = SiegelBlocksKey(q);
  auto iter = cache.map.find(key);
  if (iter != cache.map.end()) {
    return iter->second;
  }
  UnivariatePolynomial<T> res;
  int64_t p = q.p;
  T p_T = ScalarFromInt64<T>(p);
  int dim = q.dim();
  if (dim == 1) {
    res = SiegelSeriesDim1<T>(q);
  } else if (dim == 2) {
    res = SiegelSeriesDim2<T>(q);
  } else {
    // Theorem 4.1 applies for odd p always, and for p = 2 when the first
    // scale has exactly one unit diagonal entry.
    bool use41 = (p != 2);
    if (p == 2) {
      int max_e = q.blocks[0].e;
      int n_unit = 0;
      for (auto const &b : q.blocks) {
        if (b.e == max_e && b.kind == SiegelBlockKind::Unit) {
          n_unit++;
        }
      }
      use41 = (n_unit == 1);
    }
    if (use41) {
      SiegelJordanBlocks q2 = BlocksTail(q, 1);
      UnivariatePolynomial<T> F2 = SiegelSeriesPolynomialBlocks(q2, cache);
      RationalFunction<T> sum =
          KatsuradaC1<T>(q, q2) * RatFromPoly(PolyScaleVariable(F2, p_T)) +
          KatsuradaC0<T>(q, q2) * RatFromPoly(F2);
      res = RatToPolynomial(sum);
    } else {
      std::vector<SiegelJordanBlock> b1;
      size_t start;
      if (q.blocks[0].kind != SiegelBlockKind::Unit) {
        b1 = {q.blocks[0]};
        start = 1;
      } else {
        b1 = {q.blocks[0], q.blocks[1]};
        start = 2;
      }
      SiegelJordanBlocks q2 = BlocksTail(q, start);
      std::vector<RationalFunction<T>> c =
          KatsuradaTheorem42Coefficients<T>(b1, q2);
      RationalFunction<T> const &c11 = c[0], &c10 = c[1], &c21 = c[2],
                                &c20 = c[3];
      T two(2), four(4);
      UnivariatePolynomial<T> F2 = SiegelSeriesPolynomialBlocks(q2, cache);
      RationalFunction<T> coef0 = c11 * RatScaleVariable(c21, two);
      RationalFunction<T> coef1 = c11 * RatScaleVariable(c20, two) + c10 * c21;
      RationalFunction<T> coef2 = c10 * c20;
      RationalFunction<T> sum =
          coef0 * RatFromPoly(PolyScaleVariable(F2, four)) +
          coef1 * RatFromPoly(PolyScaleVariable(F2, two)) +
          coef2 * RatFromPoly(F2);
      res = RatToPolynomial(sum);
    }
  }
  cache.map[key] = res;
  return res;
}

// The Jordan blocks of a half-integral matrix B over Z_p.
template <typename T, typename Tint>
SiegelJordanBlocks SiegelJordanBlocksOf(MyMatrix<T> const &B, int64_t p) {
  if (p == 2) {
    return JordanBlocks2<T, Tint>(B);
  }
  return JordanBlocksOdd<T, Tint>(B, p);
}

// F_p(B, X) for B = S/2, S an even Gram matrix (positive definite).
template <typename T, typename Tint>
UnivariatePolynomial<T> SiegelSeriesPolynomial(MyMatrix<T> const &S, int64_t p,
                                               SiegelSeriesCache<T> &cache) {
  MyMatrix<T> B = S / T(2);
  return SiegelSeriesPolynomialBlocks(SiegelJordanBlocksOf<T, Tint>(B, p),
                                      cache);
}

/*
  alpha_p(S, T) for T unimodular at p of rank 2k whose character at p is
  xi_T = kronecker((-1)^k det T, p) (at p = 2: T even unimodular, +1 for the
  class H^k and -1 for H^{k-1} + Y).
 */
template <typename T, typename Tint>
T AlphaUnimodularTarget(MyMatrix<T> const &S, int64_t p, int k, int xi_T,
                        SiegelSeriesCache<T> &cache) {
  int m = S.rows();
  T p_T = ScalarFromInt64<T>(p);
  T X = T(xi_T) / PowerInt(p_T, k);
  UnivariatePolynomial<T> F = SiegelSeriesPolynomial<T, Tint>(S, p, cache);
  T val = PolyEvaluate(F, X) * (T(1) - X);
  for (int i = 1; i <= m / 2; i++) {
    val *= T(1) - PowerInt(p_T, 2 * i) * X * X;
  }
  if (m % 2 == 0) {
    Tint detS_int = UniversalScalarConversion<Tint, T>(DeterminantMat(S));
    int64_t detS = UniversalScalarConversion<int64_t, Tint>(detS_int);
    int64_t DB = FundamentalDiscriminant(SignedPower(m / 2, detS));
    int xi_B = QuadraticCharacter(DB, p);
    val /= T(1) - PowerInt(p_T, m / 2) * T(xi_B) * X;
  }
  if (p == 2) {
    val *= PowerInt(T(2), m - 1);
  }
  return val;
}

// clang-format off
#endif  // SRC_GENERA_SIEGELSERIES_H_
// clang-format on
