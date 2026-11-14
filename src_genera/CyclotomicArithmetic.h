// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_CYCLOTOMICARITHMETIC_H_
#define SRC_GENERA_CYCLOTOMICARITHMETIC_H_

// clang-format off
#include "GeneraArithmetic.h"
#include "SiegelSpecialValues.h"
#include "MAT_Matrix.h"
#include <cstdint>
#include <vector>
// clang-format on

/*
  Cyclotomic fields and Dirichlet characters: the cyclotomic polynomials,
  the group of the characters modulo d and their conductors, elements of
  Q(zeta_e) and their norms, the generalized Bernoulli numbers and the
  products of Dirichlet L-values at negative integers.
 */

// The polynomials have integer coefficients, by increasing degree.
inline std::vector<int64_t>
ExactMonicPolynomialDivision(std::vector<int64_t> num,
                             std::vector<int64_t> const &den) {
  int dn = static_cast<int>(num.size()) - 1;
  int dd = static_cast<int>(den.size()) - 1;
  std::vector<int64_t> quo(dn - dd + 1, 0);
  for (int k = dn - dd; k >= 0; k--) {
    int64_t c = num[k + dd];
    quo[k] = c;
    for (int j = 0; j <= dd; j++) {
      num[k + j] -= c * den[j];
    }
  }
  for (int j = 0; j < dd; j++) {
    if (num[j] != 0) {
      std::cerr << "CYCLO: the division of polynomials is not exact\n";
      throw TerminalException{1};
    }
  }
  return quo;
}

inline std::vector<int64_t> CyclotomicPolynomial(int d) {
  std::vector<int64_t> num(d + 1, 0);
  num[0] = -1;
  num[d] = 1;
  for (int e = 1; e < d; e++) {
    if (d % e == 0) {
      num = ExactMonicPolynomialDivision(num, CyclotomicPolynomial(e));
    }
  }
  return num;
}

// Q(zeta_d) = Q(zeta_{d/2}) for d = 2 mod 4.
inline int64_t NormalizeCyclotomicIndex(int64_t d) {
  return (d % 4 == 2) ? d / 2 : d;
}

/*
  The Dirichlet characters modulo d: the units, a decomposition of (Z/d)^*
  as a product of cyclic groups of orders ord[j] with the discrete
  logarithms of every unit, and the exponent e. A character is a vector t
  with 0 <= t[j] < ord[j]; chi(a) = zeta_e^{sum_j t[j] log_j(a) e / ord[j]}.
 */
struct DirichletGroup {
  int64_t d;
  int64_t e;
  std::vector<int64_t> units;
  std::vector<int64_t> ord;
  // logs[a]: the discrete logarithms of the unit a (empty for a non-unit)
  std::vector<std::vector<int64_t>> logs;
};

inline DirichletGroup ComputeDirichletGroup(int64_t d) {
  DirichletGroup G;
  G.d = d;
  G.logs.assign(d, {});
  // the cyclic factors, one or two per prime power q^a dividing d, each by
  // a generator modulo q^a lifted to Z/d by the CRT (1 at the other primes)
  std::vector<std::pair<int64_t, int64_t>> factors;
  {
    int64_t x = d;
    for (int64_t q = 2; q <= x; q++) {
      if (x % q == 0) {
        int64_t qa = 1;
        while (x % q == 0) {
          x /= q;
          qa *= q;
        }
        factors.push_back({q, qa});
      }
    }
  }
  // generators modulo q^a with their orders
  struct Cyc {
    int64_t qa;
    int64_t gen;
    int64_t order;
  };
  std::vector<Cyc> cycs;
  for (auto &[q, qa] : factors) {
    if (q != 2) {
      int64_t phi = qa / q * (q - 1);
      for (int64_t g = 2; g < qa; g++) {
        if (std::gcd(g, qa) != 1) {
          continue;
        }
        int64_t x = 1, o = 0;
        do {
          x = (x * g) % qa;
          o++;
        } while (x != 1);
        if (o == phi) {
          cycs.push_back({qa, g, phi});
          break;
        }
      }
    } else if (qa == 4) {
      cycs.push_back({4, 3, 2});
    } else if (qa >= 8) {
      cycs.push_back({qa, qa - 1, 2});
      cycs.push_back({qa, 5, qa / 4});
    }
  }
  int64_t e = 1;
  for (auto &c : cycs) {
    G.ord.push_back(c.order);
    e = std::lcm(e, c.order);
  }
  G.e = e;
  for (int64_t a = 1; a < d || (d == 1 && a == 1); a++) {
    if (std::gcd(a, d) != 1) {
      continue;
    }
    G.units.push_back(a);
  }
  // discrete logarithms, factor by factor
  for (int64_t a : G.units) {
    std::vector<int64_t> lg;
    size_t j = 0;
    while (j < cycs.size()) {
      int64_t qa = cycs[j].qa;
      int64_t r = a % qa;
      if (j + 1 < cycs.size() && cycs[j + 1].qa == qa) {
        // 2^a with a >= 3: r = (-1)^s 5^t
        int64_t s = 0;
        int64_t rr = r;
        if (rr % 4 == 3) {
          s = 1;
          rr = qa - rr;
        }
        int64_t x = 1, t = 0;
        while (x != rr) {
          x = (x * 5) % qa;
          t++;
        }
        lg.push_back(s);
        lg.push_back(t);
        j += 2;
      } else {
        int64_t x = 1, t = 0;
        while (x != r) {
          x = (x * cycs[j].gen) % qa;
          t++;
        }
        lg.push_back(t);
        j += 1;
      }
    }
    G.logs[a] = lg;
  }
  return G;
}

// The exponent of zeta_e of chi(a), a unit modulo d.
inline int64_t CharacterExponent(DirichletGroup const &G,
                                 std::vector<int64_t> const &t, int64_t a) {
  std::vector<int64_t> const &lg = G.logs[a];
  int64_t s = 0;
  for (size_t j = 0; j < t.size(); j++) {
    s += t[j] * lg[j] * (G.e / G.ord[j]);
  }
  return s % G.e;
}

inline std::vector<std::vector<int64_t>>
ListCharacters(DirichletGroup const &G) {
  std::vector<std::vector<int64_t>> res{{}};
  for (int64_t o : G.ord) {
    std::vector<std::vector<int64_t>> next;
    for (auto &t : res) {
      for (int64_t x = 0; x < o; x++) {
        std::vector<int64_t> t2 = t;
        t2.push_back(x);
        next.push_back(t2);
      }
    }
    res = std::move(next);
  }
  return res;
}

/*
  The conductor f of chi and the exponents of the primitive character on
  the residues modulo f (-1 for the non-units).
 */
inline std::pair<int64_t, std::vector<int64_t>>
PrimitiveCharacter(DirichletGroup const &G, std::vector<int64_t> const &t) {
  int64_t d = G.d;
  for (int64_t f = 1; f <= d; f++) {
    if (d % f != 0) {
      continue;
    }
    bool trivial = true;
    for (int64_t a : G.units) {
      if ((a - 1) % f == 0 && CharacterExponent(G, t, a) != 0) {
        trivial = false;
        break;
      }
    }
    if (!trivial) {
      continue;
    }
    std::vector<int64_t> val(f, -1);
    for (int64_t a = 0; a < f; a++) {
      if (std::gcd(a, f) != 1 && f > 1) {
        continue;
      }
      // a lift of a modulo f to a unit modulo d
      for (int64_t b = a; b < d + f; b += f) {
        int64_t bb = b % d;
        if (std::gcd(bb, d) == 1 || d == 1) {
          val[a] = (d == 1) ? 0 : CharacterExponent(G, t, bb);
          break;
        }
      }
    }
    if (f == 1) {
      val[0] = 0;
    }
    return {f, val};
  }
  std::cerr << "CYCLO: no conductor found\n";
  throw TerminalException{1};
}

/*
  An element of Q(zeta_e) as its coefficients on 1, zeta, ..., zeta^{e-1}
  (not reduced), and the reduction modulo Phi_e to the power basis of
  degree phi(e).
 */
template <typename T>
std::vector<T> ReduceCyclotomic(std::vector<T> const &x, int64_t e) {
  std::vector<int64_t> Phi = CyclotomicPolynomial(static_cast<int>(e));
  int deg = static_cast<int>(Phi.size()) - 1;
  std::vector<T> r = x;
  for (int k = static_cast<int>(r.size()) - 1; k >= deg; k--) {
    T c = r[k];
    if (c != 0) {
      for (int j = 0; j <= deg; j++) {
        r[k - deg + j] -= c * ScalarFromInt64<T>(Phi[j]);
      }
    }
  }
  r.resize(deg);
  return r;
}

// N_{Q(zeta_e)/Q}(x), x in the power basis, as the determinant of the
// multiplication by x.
template <typename T> T NormCyclotomic(std::vector<T> const &x, int64_t e) {
  int deg = static_cast<int>(x.size());
  MyMatrix<T> M(deg, deg);
  std::vector<T> cur = x;
  for (int j = 0; j < deg; j++) {
    for (int i = 0; i < deg; i++) {
      M(j, i) = cur[i];
    }
    // multiply by zeta
    std::vector<T> next(deg + 1, T(0));
    for (int i = 0; i < deg; i++) {
      next[i + 1] = cur[i];
    }
    cur = ReduceCyclotomic(next, e);
  }
  return DeterminantMat(M);
}

/*
  B_{k, chi} for the primitive character of conductor f with the exponents
  val, as an element of Q(zeta_e) (not reduced):
      B_{k,chi} = sum_{j=0}^{k} C(k,j) B_j f^{j-1} S_{k-j},
      S_t = sum_{a=1}^{f} chi(a) a^t,
  with the convention B_{1, 1} = 1/2 of the trivial character.
 */
template <typename T>
std::vector<T> GeneralizedBernoulliCyclotomic(int k, int64_t f,
                                              std::vector<int64_t> const &val,
                                              int64_t e) {
  std::vector<T> res(e, T(0));
  std::vector<T> B = BernoulliNumbers<T>(k);
  if (f == 1) {
    res[0] = (k == 1) ? T(1) / T(2) : B[k];
    return res;
  }
  std::vector<std::vector<T>> S(k + 1, std::vector<T>(e, T(0)));
  for (int64_t a = 1; a <= f; a++) {
    int64_t ex = val[a % f];
    if (ex < 0) {
      continue;
    }
    T a_T = ScalarFromInt64<T>(a);
    T pw(1);
    for (int t = 0; t <= k; t++) {
      S[t][ex] += pw;
      pw *= a_T;
    }
  }
  T f_T = ScalarFromInt64<T>(f);
  T fpow = T(1) / f_T;
  for (int j = 0; j <= k; j++) {
    T c = BinomialCoefficient<T>(k, j) * B[j] * fpow;
    for (int64_t u = 0; u < e; u++) {
      res[u] += c * S[k - j][u];
    }
    fpow *= f_T;
  }
  return res;
}

// The exact r-th root of a positive rational which is an r-th power.
template <typename T>
inline T RationalRoot(T const &x, int r) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  Tint num = GetNumerator_z(x);
  Tint den = GetDenominator_z(x);
  std::optional<Tint> rn = exact_integer_root(num, static_cast<size_t>(r));
  std::optional<Tint> rd = exact_integer_root(den, static_cast<size_t>(r));
  if (!rn || !rd) {
    std::cerr << "CYCLO: the product of the norms is not a power\n";
    throw TerminalException{1};
  }
  T res = UniversalScalarConversion<T, Tint>(*rn);
  res /= UniversalScalarConversion<T, Tint>(*rd);
  return res;
}

/*
  prod L(1 - k, chi) over the characters chi modulo d with chi(-1) =
  (-1)^parity, taken primitive.
 */
template <typename T>
inline T ProductDirichletLValues(int64_t d, int k, int parity,
                                         std::vector<int64_t> const &H = {},
                                         std::vector<int64_t> const &Hsub = {}) {
  DirichletGroup G = ComputeDirichletGroup(d);
  int64_t e = G.e;
  int64_t phi_e = EulerPhi(e);
  T real_prod(1);
  T norm_prod(1);
  size_t n_char = 0;
  size_t n_complex = 0;
  for (auto &t : ListCharacters(G)) {
    // the characters of the subfield fixed by H: trivial on H
    bool on_field = true;
    for (int64_t h : H) {
      if (CharacterExponent(G, t, h) != 0) {
        on_field = false;
      }
    }
    if (!on_field) {
      continue;
    }
    // with Hsub (a supergroup of H): only the characters not trivial on
    // it, those of the field fixed by H that are not of the field fixed by
    // Hsub
    if (!Hsub.empty()) {
      bool on_sub = true;
      for (int64_t h : Hsub) {
        if (CharacterExponent(G, t, h) != 0) {
          on_sub = false;
        }
      }
      if (on_sub) {
        continue;
      }
    }
    int64_t ex_m1 = (d <= 2) ? 0 : CharacterExponent(G, t, d - 1);
    int par = (ex_m1 == 0) ? 0 : 1;
    if (par != parity) {
      continue;
    }
    n_char++;
    std::pair<int64_t, std::vector<int64_t>> prim = PrimitiveCharacter(G, t);
    std::vector<T> Bk =
        GeneralizedBernoulliCyclotomic<T>(k, prim.first, prim.second, e);
    // real when chi^2 is trivial, the values being then +-1
    bool is_real = true;
    for (size_t j = 0; j < t.size(); j++) {
      if ((2 * t[j]) % G.ord[j] != 0) {
        is_real = false;
      }
    }
    if (is_real) {
      T val = Bk[0];
      if (e % 2 == 0) {
        val -= Bk[e / 2];
      }
      real_prod *= -val / T(k);
    } else {
      n_complex++;
      std::vector<T> Bred = ReduceCyclotomic(Bk, e);
      T kk(k);
      T kpow = MyPow<T>(kk, phi_e);
      // N(-B / k) = N(B) / k^{phi(e)}, phi(e) even
      norm_prod *= NormCyclotomic(Bred, e) / kpow;
    }
  }
#ifdef DEBUG_CYCLOTOMIC_HERMITIAN_MASS
  std::cerr << "CYCLO: d=" << d << " k=" << k << " parity=" << parity
            << " characters=" << n_char << " complex=" << n_complex << "\n";
#endif
  (void)n_char;
  if (n_complex == 0) {
    return real_prod;
  }
  // the non-real values come in conjugate pairs: their product is positive
  return real_prod * RationalRoot(norm_prod, static_cast<int>(phi_e));
}

// clang-format off
#endif  // SRC_GENERA_CYCLOTOMICARITHMETIC_H_
// clang-format on
