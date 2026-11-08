// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_CYCLOTOMICHERMITIANMASS_H_
#define SRC_GENERA_CYCLOTOMICHERMITIANMASS_H_

// clang-format off
#include "SiegelSpecialValues.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include <optional>
#include <cstdint>
#include <numeric>
#include <map>
#include <set>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_CYCLOTOMIC_HERMITIAN_MASS
#endif

/*
  Masses of Hermitian lattices over E = Q(zeta_d) relative to its maximal
  real subfield K = Q(zeta_d)^+, n = [K : Q] = phi(d) / 2, exactly.

  The standard mass of a lattice of rank m is (Hecke.jl,
  QuadForm/Herm/Mass.jl, after Gan-Hanke-Yu)

      2^{1 - n m} prod_{i even, i <= m} zeta_K(1 - i)
                  prod_{i odd, i <= m} L(1 - i, chi_{E/K}),

  and the mass is its absolute value times the local factors at the bad
  primes. Hecke evaluates the L-values numerically and rounds against a
  bound on the denominator. Here E is abelian, so they are products of
  Dirichlet L-values: zeta_K(s) = prod_{chi even} L(s, chi), L(s, chi_{E/K})
  = prod_{chi odd} L(s, chi), over the characters modulo d taken primitive,
  and L(1 - k, chi) = -B_{k, chi} / k. The values of the characters are
  powers of zeta_e, e the exponent of (Z/d)^*. For a real character B_{k,
  chi} is rational. The non-real characters come in conjugate pairs, so
  their product is positive, and it is the phi(e)-th root of the product
  of the norms N_{Q(zeta_e)/Q}(B_{k, chi}), each a determinant of a
  multiplication matrix.

  Rank one. A lattice of rank one is a fractional ideal a with h(x, y) =
  c x y-bar; its unitary group is the group of the roots of unity of E and
  its genus is determined by E alone. The local factor is 1 at the primes
  unramified in E/K (a single Jordan component), and E/K is ramified at a
  finite prime only when d (made odd or divisible by 4) is a prime power,
  at the prime over it, where the factor is 1/2. So the mass is
  2^{1-n} |L(0, chi_{E/K})| times 1/2 for a prime power, which is h^- / w
  by the class number formula (h^- the relative class number, w the number
  of roots of unity).
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

inline int64_t CycloGcd(int64_t a, int64_t b) {
  return std::gcd(a, b);
}

inline int64_t EulerPhi(int64_t d) {
  int64_t res = 0;
  for (int64_t a = 1; a <= d; a++) {
    if (CycloGcd(a, d) == 1) {
      res++;
    }
  }
  return res;
}

// Q(zeta_d) = Q(zeta_{d/2}) for d = 2 mod 4.
inline int64_t NormalizeCyclotomicIndex(int64_t d) {
  return (d % 4 == 2) ? d / 2 : d;
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
        if (CycloGcd(g, qa) != 1) {
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
    if (CycloGcd(a, d) != 1) {
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
      if (CycloGcd(a, f) != 1 && f > 1) {
        continue;
      }
      // a lift of a modulo f to a unit modulo d
      for (int64_t b = a; b < d + f; b += f) {
        int64_t bb = b % d;
        if (CycloGcd(bb, d) == 1 || d == 1) {
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
    if (c == 0) {
      continue;
    }
    for (int j = 0; j <= deg; j++) {
      r[k - deg + j] -= c * ScalarFromInt64<T>(Phi[j]);
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
inline mpq_class RationalRoot(mpq_class const &x, int r) {
  mpz_class num = x.get_num();
  mpz_class den = x.get_den();
  mpz_class rn, rd;
  int ex1 = mpz_root(rn.get_mpz_t(), num.get_mpz_t(), r);
  int ex2 = mpz_root(rd.get_mpz_t(), den.get_mpz_t(), r);
  if (ex1 == 0 || ex2 == 0) {
    std::cerr << "CYCLO: the product of the norms is not a power\n";
    throw TerminalException{1};
  }
  return mpq_class(rn, rd);
}

/*
  prod L(1 - k, chi) over the characters chi modulo d with chi(-1) =
  (-1)^parity, taken primitive.
 */
inline mpq_class ProductDirichletLValues(int64_t d, int k, int parity,
                                         std::vector<int64_t> const &H = {},
                                         std::vector<int64_t> const &Hsub = {}) {
  using T = mpq_class;
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
      T kpow(1);
      for (int64_t i = 0; i < phi_e; i++) {
        kpow *= kk;
      }
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

// The standard mass of a Hermitian lattice of rank m over Q(zeta_d), in
// absolute value.
inline mpq_class StandardHermitianMassCyclotomic(int64_t d_in, int m) {
  int64_t d = NormalizeCyclotomicIndex(d_in);
  int64_t n = EulerPhi(d) / 2;
  mpq_class res(1);
  int64_t ex = 1 - n * m;
  mpq_class two(2);
  for (int64_t i = 0; i < (ex >= 0 ? ex : -ex); i++) {
    if (ex >= 0) {
      res *= two;
    } else {
      res /= two;
    }
  }
  for (int i = 1; i <= m; i++) {
    res *= ProductDirichletLValues(d, i, i % 2);
  }
  return (res < 0) ? mpq_class(-res) : res;
}

// The mass of the genus of a Hermitian lattice of rank one over
// Q(zeta_d), which depends only on d.
inline mpq_class RankOneHermitianMassCyclotomic(int64_t d_in) {
  int64_t d = NormalizeCyclotomicIndex(d_in);
  mpq_class res = StandardHermitianMassCyclotomic(d, 1);
  if (IsPrimePower(d)) {
    res /= 2;
  }
  return res;
}

/*
  Abelian CM fields. E is the subfield of Q(zeta_n) fixed by a subgroup H
  of (Z/n)^* not containing -1, so that E is CM with maximal real subfield
  K the field fixed by <H, -1>. Its Dirichlet characters are those modulo n
  trivial on H, the even ones being those of K, so the standard mass is
  the formula above with the characters restricted: zeta_K(s) and
  L(s, chi_{E/K}) are the products of the L(s, chi) over the even and the
  odd characters of E. The cyclotomic case is H = {1}.
 */
struct AbelianCMField {
  int64_t n;
  std::vector<int64_t> H;
};

inline AbelianCMField MakeAbelianCMField(int64_t n,
                                         std::vector<int64_t> const &gens) {
  // the subgroup generated by gens
  std::set<int64_t> H{1 % std::max<int64_t>(n, 1)};
  std::vector<int64_t> queue{1 % std::max<int64_t>(n, 1)};
  while (!queue.empty()) {
    int64_t a = queue.back();
    queue.pop_back();
    for (int64_t g : gens) {
      int64_t b = ((a * (g % n)) % n + n) % n;
      if (CycloGcd(b, n) != 1) {
        std::cerr << "CYCLO: " << g << " is not a unit modulo " << n << "\n";
        throw TerminalException{1};
      }
      if (H.insert(b).second) {
        queue.push_back(b);
      }
    }
  }
  if (n > 2 && H.count(n - 1) == 1) {
    std::cerr << "CYCLO: -1 is in H, the field is not CM\n";
    throw TerminalException{1};
  }
  return {n, std::vector<int64_t>(H.begin(), H.end())};
}

inline int64_t AbelianFieldDegree(AbelianCMField const &F) {
  return EulerPhi(F.n) / static_cast<int64_t>(F.H.size());
}

inline mpq_class StandardHermitianMassAbelian(AbelianCMField const &F, int m) {
  int64_t n = AbelianFieldDegree(F) / 2;
  mpq_class res(1);
  int64_t ex = 1 - n * m;
  mpq_class two(2);
  for (int64_t i = 0; i < (ex >= 0 ? ex : -ex); i++) {
    if (ex >= 0) {
      res *= two;
    } else {
      res /= two;
    }
  }
  for (int i = 1; i <= m; i++) {
    res *= ProductDirichletLValues(F.n, i, i % 2, F.H);
  }
  return (res < 0) ? mpq_class(-res) : res;
}

/*
  The decomposition of a rational prime q in the field of a group X of
  Dirichlet characters: e = [X : X_unr] with X_unr the characters of
  conductor prime to q, f the order of q on X_unr (the lcm of the orders of
  the chi(q)), and g = |X| / (e f).
 */
struct PrimeDecompositionData {
  int64_t e;
  int64_t f;
  int64_t g;
};

inline PrimeDecompositionData
PrimeDecompositionCharacters(DirichletGroup const &G,
                             std::vector<std::vector<int64_t>> const &X,
                             int64_t q) {
  int64_t n_unr = 0;
  int64_t f = 1;
  for (auto const &t : X) {
    std::pair<int64_t, std::vector<int64_t>> prim = PrimitiveCharacter(G, t);
    int64_t cond = prim.first;
    if (cond % q == 0) {
      continue;
    }
    n_unr++;
    int64_t ex = prim.second[q % cond];
    int64_t ord = G.e / std::gcd(G.e, ex);
    f = std::lcm(f, ord);
  }
  int64_t sz = static_cast<int64_t>(X.size());
  int64_t e = sz / n_unr;
  return {e, f, sz / (e * f)};
}

inline std::vector<std::vector<int64_t>>
FieldCharacters(DirichletGroup const &G, std::vector<int64_t> const &H,
                bool even_only) {
  std::vector<std::vector<int64_t>> res;
  for (auto &t : ListCharacters(G)) {
    bool ok = true;
    for (int64_t h : H) {
      if (CharacterExponent(G, t, h) != 0) {
        ok = false;
      }
    }
    if (ok && even_only && G.d > 2 && CharacterExponent(G, t, G.d - 1) != 0) {
      ok = false;
    }
    if (ok) {
      res.push_back(t);
    }
  }
  return res;
}

/*
  The finite primes of K ramified in E: over a rational prime q they all
  ramify when e_q(E) = 2 e_q(K), and there are g_q(K) of them. Returned as
  the pairs (q, number of primes of K over q ramified in E).
 */
inline std::vector<std::pair<int64_t, int64_t>>
RamifiedPrimesCM(AbelianCMField const &F) {
  DirichletGroup G = ComputeDirichletGroup(F.n);
  std::vector<std::vector<int64_t>> XE = FieldCharacters(G, F.H, false);
  std::vector<std::vector<int64_t>> XK = FieldCharacters(G, F.H, true);
  std::vector<std::pair<int64_t, int64_t>> res;
  int64_t x = F.n;
  for (int64_t q = 2; q <= x; q++) {
    if (x % q != 0) {
      continue;
    }
    while (x % q == 0) {
      x /= q;
    }
    PrimeDecompositionData dE = PrimeDecompositionCharacters(G, XE, q);
    PrimeDecompositionData dK = PrimeDecompositionCharacters(G, XK, q);
    if (dE.e == 2 * dK.e) {
      res.push_back({q, dK.g});
    }
  }
  return res;
}

/*
  The mass of the genus of a Hermitian lattice of rank one over E: the
  standard mass times 1/2 for each finite prime of K ramified in E (the
  local factor of a single Jordan block at a ramified prime), checked
  against Hecke.jl on Q(sqrt(-D)), D = 1, 2, 6, 7, 15, 23.
 */
inline mpq_class RankOneHermitianMassAbelian(AbelianCMField const &F) {
  mpq_class res = StandardHermitianMassAbelian(F, 1);
  for (auto const &pr : RamifiedPrimesCM(F)) {
    for (int64_t i = 0; i < pr.second; i++) {
      res /= 2;
    }
  }
  return res;
}

/*
  Higher rank: the local factors.

  The Hermitian lattice is given as a Z-lattice with Gram matrix B of the
  trace form b = Tr_{E/Q} h and the action of zeta_d (an isometry S on the
  row vectors, Phi_d(S) = 0). At a rational prime q, Z_q[x] / Phi_d is the
  product of the completions O_P of O_E at the primes P over q, cut out by
  idempotents: Phi_d = Phi_{d'}^{phi(q^a)} modulo q for d = q^a d', the
  factors g_i of Phi_{d'} modulo q (all of degree f, the order of q modulo
  d') give the primes P_i, and the idempotent of P_i, congruent to 1 modulo
  g_i^{phi(q^a)} and to 0 modulo the others, is lifted to Z / q^N by e <-
  3 e^2 - 2 e^3. The complex conjugation maps P_i to the prime of the
  reciprocal polynomial of g_i: the prime of K below is split in E when
  the two differ, and otherwise inert, or ramified when d is a power of q.

  A uniformizer of O_P is q when q does not divide d, and 1 - zeta^{d'}
  otherwise. The Jordan decomposition of h at P has the scales s_i with
  ranks m_i; the trace dual and the Hermitian dual differ by the different,
  L^{#b} = D^{-1} L^{#h}, so e_P L^{#b} / e_P L is the O_P-module
  (+)_i (O_P / P^{s_i + delta})^{m_i}, delta = v_P(D_{E/Q}). Its invariants
  are read off the lengths of the pi^k e_P L^{#b} + L over L, computed as
  indices of Z-lattices. For a split prime the module of P alone gives the
  elementary divisors of the pairing between the lattices at P and P-bar.

  The local factor at a prime of K unramified in E is that of Hecke.jl
  (QuadForm/Herm/Mass.jl, local_factor): 1 for a single Jordan block,
  otherwise q_K^N |GL or U_m(q_K)| / prod_i |GL or U_{m_i}(q_K)|, with N =
  sum_i s_i m_i (sum_{j<i} m_j - sum_{j>i} m_j) for the scales in
  increasing order, and the group orders normalized as prod_j (1 -
  (+-q_K)^{-j}). It does not change when all the scales are shifted, so
  delta does not matter there. The ramified prime, when d is a prime
  power, also needs the norm classes of the determinants of the blocks of
  even scale, and is not treated: the mass is then not computed.
 */

using PolyQ = std::vector<int64_t>;

// a b modulo q without overflow, for 0 <= a, b < q < 2^62
inline int64_t MulModQ(int64_t a, int64_t b, int64_t q) {
  __int128 p = static_cast<__int128>(a) * static_cast<__int128>(b);
  int64_t r = static_cast<int64_t>(p % q);
  return (r < 0) ? r + q : r;
}

inline void PolyTrim(PolyQ &a) {
  while (!a.empty() && a.back() == 0) {
    a.pop_back();
  }
}

inline PolyQ PolyReduce(PolyQ a, int64_t q) {
  for (auto &x : a) {
    x = ((x % q) + q) % q;
  }
  PolyTrim(a);
  return a;
}

inline PolyQ PolyMulQ(PolyQ const &a, PolyQ const &b, int64_t q) {
  if (a.empty() || b.empty()) {
    return {};
  }
  PolyQ c(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      c[i + j] = (c[i + j] + MulModQ(a[i], b[j], q)) % q;
    }
  }
  PolyTrim(c);
  return c;
}

inline PolyQ PolySubQ(PolyQ a, PolyQ const &b, int64_t q) {
  if (a.size() < b.size()) {
    a.resize(b.size(), 0);
  }
  for (size_t i = 0; i < b.size(); i++) {
    a[i] = ((a[i] - b[i]) % q + q) % q;
  }
  PolyTrim(a);
  return a;
}

inline int64_t InvModQ(int64_t a, int64_t q) {
  int64_t r0 = q, r1 = ((a % q) + q) % q, s0 = 0, s1 = 1;
  while (r1 != 0) {
    int64_t t = r0 / r1;
    int64_t r2 = r0 - t * r1;
    r0 = r1;
    r1 = r2;
    int64_t s2 = static_cast<int64_t>((static_cast<__int128>(s0) - static_cast<__int128>(t) * s1) % q);
    s0 = s1;
    s1 = s2;
  }
  return ((s0 % q) + q) % q;
}

// quotient and remainder modulo q, b nonzero
inline std::pair<PolyQ, PolyQ> PolyDivModQ(PolyQ a, PolyQ const &b, int64_t q) {
  PolyTrim(a);
  if (a.size() < b.size()) {
    return {{}, a};
  }
  int64_t inv = InvModQ(b.back(), q);
  PolyQ quo(a.size() - b.size() + 1, 0);
  for (size_t k = a.size() - b.size() + 1; k-- > 0;) {
    int64_t c = MulModQ(a[k + b.size() - 1], inv, q);
    quo[k] = c;
    for (size_t j = 0; j < b.size(); j++) {
      a[k + j] = ((a[k + j] - MulModQ(c, b[j], q)) % q + q) % q;
    }
  }
  PolyTrim(a);
  PolyTrim(quo);
  return {quo, a};
}

inline PolyQ PolyMonicQ(PolyQ a, int64_t q) {
  if (a.empty()) {
    return a;
  }
  int64_t inv = InvModQ(a.back(), q);
  for (auto &x : a) {
    x = MulModQ(x, inv, q);
  }
  return a;
}

inline PolyQ PolyGcdQ(PolyQ a, PolyQ b, int64_t q) {
  PolyTrim(a);
  PolyTrim(b);
  while (!b.empty()) {
    PolyQ r = PolyDivModQ(a, b, q).second;
    a = b;
    b = r;
  }
  return PolyMonicQ(a, q);
}

// u with u a = gcd modulo b (a and b coprime): the inverse of a modulo b
inline PolyQ PolyInverseModQ(PolyQ const &a, PolyQ const &b, int64_t q) {
  PolyQ r0 = b, r1 = PolyDivModQ(a, b, q).second;
  PolyQ s0 = {}, s1 = {1};
  while (!r1.empty()) {
    std::pair<PolyQ, PolyQ> qr = PolyDivModQ(r0, r1, q);
    PolyQ s2 = PolySubQ(s0, PolyMulQ(qr.first, s1, q), q);
    r0 = r1;
    r1 = qr.second;
    s0 = s1;
    s1 = s2;
  }
  // r0 is a nonzero constant
  if (r0.size() != 1) {
    std::cerr << "CYCLO: the polynomials are not coprime\n";
    throw TerminalException{1};
  }
  int64_t inv = InvModQ(r0[0], q);
  for (auto &x : s0) {
    x = MulModQ(x, inv, q);
  }
  return PolyDivModQ(s0, b, q).second;
}

inline PolyQ PolyPowModQ(PolyQ base, uint64_t e, PolyQ const &m, int64_t q) {
  PolyQ res = {1};
  base = PolyDivModQ(base, m, q).second;
  while (e > 0) {
    if (e & 1) {
      res = PolyDivModQ(PolyMulQ(res, base, q), m, q).second;
    }
    base = PolyDivModQ(PolyMulQ(base, base, q), m, q).second;
    e >>= 1;
  }
  return res;
}

/*
  The irreducible factors (monic) of the squarefree g modulo q, all of
  degree f, by the equal degree splitting of Cantor and Zassenhaus (q odd)
  or by the trace map (q = 2). The random polynomials are replaced by a
  deterministic sequence, which splits as well.
 */
inline std::vector<PolyQ> EqualDegreeFactorization(PolyQ const &g, int f,
                                                   int64_t q) {
  int deg = static_cast<int>(g.size()) - 1;
  if (deg == f) {
    return {PolyMonicQ(g, q)};
  }
  for (int64_t c = 1;; c++) {
    // r = x^{1 + c mod deg} + c, a varied family of polynomials
    PolyQ r(2 + static_cast<size_t>(c % deg), 0);
    r.back() = 1;
    r[0] = c % q;
    r = PolyReduce(r, q);
    PolyQ h;
    if (q == 2) {
      PolyQ t = PolyDivModQ(r, g, q).second;
      h = t;
      for (int i = 1; i < f; i++) {
        t = PolyDivModQ(PolyMulQ(t, t, q), g, q).second;
        h = PolySubQ(h, PolySubQ({}, t, q), q);
      }
    } else {
      // r^{(q^f - 1)/2} = prod_{i < f} (r^{(q-1)/2})^{q^i}, which keeps the
      // exponents within 64 bits
      PolyQ t = PolyPowModQ(r, static_cast<uint64_t>((q - 1) / 2), g, q);
      h = {1};
      for (int i = 0; i < f; i++) {
        h = PolyDivModQ(PolyMulQ(h, t, q), g, q).second;
        t = PolyPowModQ(t, static_cast<uint64_t>(q), g, q);
      }
      h = PolySubQ(h, {1}, q);
    }
    PolyQ dgc = PolyGcdQ(g, h, q);
    int dd = static_cast<int>(dgc.size()) - 1;
    if (dd > 0 && dd < deg) {
      std::vector<PolyQ> res = EqualDegreeFactorization(dgc, f, q);
      std::vector<PolyQ> res2 =
          EqualDegreeFactorization(PolyDivModQ(g, dgc, q).first, f, q);
      res.insert(res.end(), res2.begin(), res2.end());
      return res;
    }
    if (c > 100000) {
      std::cerr << "CYCLO: the equal degree splitting failed\n";
      throw TerminalException{1};
    }
  }
}

inline PolyQ ReciprocalPolyQ(PolyQ const &g, int64_t q) {
  PolyQ r(g.rbegin(), g.rend());
  return PolyMonicQ(PolyReduce(r, q), q);
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

/*
  Polynomials over Z (mpz), reduced modulo Phi_d (monic) and with the
  coefficients modulo Q = q^N.
 */
using PolyZ = std::vector<mpz_class>;

inline PolyZ PolyZReduce(PolyZ a, std::vector<int64_t> const &Phi,
                         mpz_class const &Q) {
  int deg = static_cast<int>(Phi.size()) - 1;
  for (int k = static_cast<int>(a.size()) - 1; k >= deg; k--) {
    mpz_class c = a[k];
    if (c == 0) {
      continue;
    }
    for (int j = 0; j <= deg; j++) {
      a[k - deg + j] -= c * mpz_class(static_cast<long>(Phi[j]));
    }
  }
  a.resize(deg, mpz_class(0));
  for (auto &x : a) {
    x %= Q;
    if (x < 0) {
      x += Q;
    }
  }
  return a;
}

inline PolyZ PolyZMul(PolyZ const &a, PolyZ const &b) {
  PolyZ c(a.size() + b.size() - 1, mpz_class(0));
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      c[i + j] += a[i] * b[j];
    }
  }
  return c;
}

/*
  The Jordan data of h at a prime of K over q: the splitting type (0 split,
  1 inert, 2 ramified), the residue degree of the prime of K and the list
  of (scale, rank) in increasing scales.
 */
struct HermitianLocalJordan {
  int type;
  int f_K;
  std::vector<std::pair<int, int>> blocks;
};

template <typename T>
int ValuationQRational(T const &x, int64_t q) {
  // v_q of a nonzero rational
  mpq_class y = UniversalScalarConversion<mpq_class, T>(x);
  mpz_class num = y.get_num(), den = y.get_den();
  if (num < 0) {
    num = -num;
  }
  int v = 0;
  mpz_class qz(static_cast<long>(q));
  while (num % qz == 0) {
    num /= qz;
    v++;
  }
  while (den % qz == 0) {
    den /= qz;
    v--;
  }
  return v;
}

inline int RankModQ(std::vector<std::vector<int64_t>> M, int64_t q) {
  int rows = static_cast<int>(M.size());
  int cols = rows ? static_cast<int>(M[0].size()) : 0;
  int rank = 0;
  for (int c = 0; c < cols && rank < rows; c++) {
    int piv = -1;
    for (int i = rank; i < rows; i++) {
      if (M[i][c] % q != 0) {
        piv = i;
        break;
      }
    }
    if (piv < 0) {
      continue;
    }
    std::swap(M[piv], M[rank]);
    int64_t inv = InvModQ(M[rank][c], q);
    for (int i = 0; i < rows; i++) {
      if (i != rank && M[i][c] % q != 0) {
        int64_t t = MulModQ(M[i][c], inv, q);
        for (int j = 0; j < cols; j++) {
          M[i][j] = ((M[i][j] - MulModQ(t, M[rank][j], q)) % q + q) % q;
        }
      }
    }
    rank++;
  }
  return rank;
}

// v_q of the index [L + span(rows of M) : L], M rational, L = Z^r.
inline int ValuationIndexOverlattice(MyMatrix<mpq_class> const &M, int64_t q) {
  int r = M.cols();
  mpz_class D(1);
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < r; j++) {
      D = lcm(D, M(i, j).get_den());
    }
  }
  MyMatrix<mpz_class> G(M.rows() + r, r);
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < r; j++) {
      mpq_class v = M(i, j) * D;
      G(i, j) = v.get_num();
    }
  }
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < r; j++) {
      G(M.rows() + i, j) = (i == j) ? D : mpz_class(0);
    }
  }
  MyMatrix<mpz_class> Bz = GetZbasis(G);
  mpq_class idx = mpq_class(mpz_class(1)) / mpq_class(DeterminantMat(Bz));
  for (int i = 0; i < r; i++) {
    idx *= mpq_class(D);
  }
  return ValuationQRational(idx, q);
}

inline std::vector<HermitianLocalJordan>
HermitianJordanAtPrime(MyMatrix<mpq_class> const &B, MyMatrix<mpq_class> const &S,
                       int64_t d, int64_t q) {
  int r = B.rows();
  int64_t qa = 1;
  int64_t dp = d;
  while (dp % q == 0) {
    dp /= q;
    qa *= q;
  }
  int a = 0;
  for (int64_t x = qa; x > 1; x /= q) {
    a++;
  }
  int64_t e_ram = (qa == 1) ? 1 : qa / q * (q - 1);
  int64_t delta = (qa == 1) ? 0 : a * e_ram - qa / q;
  int f = MultiplicativeOrder(q, dp);
  std::vector<int64_t> PhiD = CyclotomicPolynomial(static_cast<int>(d));
  std::vector<int64_t> PhiDp = CyclotomicPolynomial(static_cast<int>(dp));
  PolyQ phid_q = PolyReduce(PolyQ(PhiD.begin(), PhiD.end()), q);
  std::vector<PolyQ> gs =
      EqualDegreeFactorization(PolyReduce(PolyQ(PhiDp.begin(), PhiDp.end()), q),
                               f, q);
  // the precision: beyond the exponent of the discriminant group at q
  int vdet = ValuationQRational(DeterminantMat(B), q);
  int N = vdet + 4;
  mpz_class Q(1);
  for (int i = 0; i < N; i++) {
    Q *= static_cast<long>(q);
  }
  // the dual basis, the uniformizer and the idempotent of each factor
  MyMatrix<mpq_class> Binv = Inverse(B);
  MyMatrix<mpq_class> Id = IdentityMat<mpq_class>(r);
  MyMatrix<mpq_class> Pi;
  if (qa == 1) {
    Pi = Id * mpq_class(static_cast<long>(q));
  } else {
    MyMatrix<mpq_class> Sd = Id;
    for (int64_t i = 0; i < dp; i++) {
      Sd = Sd * S;
    }
    Pi = Id - Sd;
  }
  auto idempotent_matrix = [&](PolyQ const &g) -> MyMatrix<mpq_class> {
    PolyQ P = {1};
    for (int64_t i = 0; i < e_ram; i++) {
      P = PolyMulQ(P, g, q);
    }
    PolyQ Qp = PolyDivModQ(phid_q, P, q).first;
    // e = Qp * (Qp^{-1} mod P): 1 modulo P, 0 modulo Qp
    PolyQ e0 = PolyMulQ(Qp, PolyInverseModQ(Qp, P, q), q);
    PolyZ ez;
    for (int64_t c : e0) {
      ez.push_back(mpz_class(static_cast<long>(c)));
    }
    ez = PolyZReduce(ez, PhiD, Q);
    for (int it = 0; it < 64; it++) {
      PolyZ e2 = PolyZReduce(PolyZMul(ez, ez), PhiD, Q);
      PolyZ e3 = PolyZReduce(PolyZMul(e2, ez), PhiD, Q);
      PolyZ en(e2.size());
      for (size_t i = 0; i < e2.size(); i++) {
        en[i] = 3 * e2[i] - 2 * e3[i];
      }
      en = PolyZReduce(en, PhiD, Q);
      if (en == ez) {
        break;
      }
      ez = en;
    }
    MyMatrix<mpq_class> E = ZeroMatrix<mpq_class>(r, r);
    MyMatrix<mpq_class> pw = Id;
    for (size_t i = 0; i < ez.size(); i++) {
      E += pw * mpq_class(ez[i]);
      pw = pw * S;
    }
    return E;
  };
  std::vector<HermitianLocalJordan> res;
  std::vector<bool> done(gs.size(), false);
  for (size_t i = 0; i < gs.size(); i++) {
    if (done[i]) {
      continue;
    }
    PolyQ gr = ReciprocalPolyQ(gs[i], q);
    size_t j = i;
    for (size_t u = 0; u < gs.size(); u++) {
      if (gs[u] == gr) {
        j = u;
      }
    }
    done[i] = true;
    done[j] = true;
    HermitianLocalJordan hj;
    if (j != i) {
      hj.type = 0;
      hj.f_K = f;
    } else if (qa > 1 && dp == 1) {
      hj.type = 2;
      hj.f_K = f;
    } else {
      hj.type = 1;
      hj.f_K = f / 2;
    }
    MyMatrix<mpq_class> E = idempotent_matrix(gs[i]);
    // the O_P-rank: the rank modulo q of E divided by [O_P : Z_q]
    std::vector<std::vector<int64_t>> Emod(r, std::vector<int64_t>(r));
    mpz_class qz(static_cast<long>(q));
    for (int u = 0; u < r; u++) {
      for (int v = 0; v < r; v++) {
        mpz_class x = E(u, v).get_num() % qz;
        if (x < 0) {
          x += qz;
        }
        Emod[u][v] = x.get_si();
      }
    }
    int rk = RankModQ(Emod, q);
    int m = rk / static_cast<int>(e_ram * f);
    // len[k] = log_q |pi^k X|, X = e_P L^# / e_P L, until it vanishes
    std::vector<int> len;
    MyMatrix<mpq_class> Mk = Binv * E;
    while (true) {
      len.push_back(ValuationIndexOverlattice(Mk, q));
      if (len.back() == 0 || len.size() > static_cast<size_t>(4 * vdet + 8)) {
        break;
      }
      Mk = Mk * Pi;
    }
    // c[k] = #{invariants t > k} = dim_{F_q}(pi^k X / pi^{k+1} X) / f
    std::vector<int> c;
    for (size_t k = 0; k + 1 < len.size(); k++) {
      c.push_back((len[k] - len[k + 1]) / f);
    }
    std::map<int, int> counts;
    int n_pos = 0;
    for (size_t j = 1; j <= c.size(); j++) {
      int nb = c[j - 1] - ((j < c.size()) ? c[j] : 0);
      if (nb > 0) {
        counts[static_cast<int>(j)] += nb;
        n_pos += nb;
      }
    }
    if (m - n_pos > 0) {
      counts[0] += m - n_pos;
    }
    for (auto &kv : counts) {
      hj.blocks.push_back({kv.first - static_cast<int>(delta), kv.second});
    }
    res.push_back(hj);
  }
  return res;
}

inline mpq_class HermitianGroupOrderNormalized(int type, int m, int64_t qK) {
  // GL (split) or U (inert): prod_{j=1}^{m} (1 - (+-qK)^{-j})
  mpq_class res(1);
  mpq_class base = (type == 0) ? mpq_class(static_cast<long>(qK))
                               : mpq_class(-static_cast<long>(qK));
  mpq_class pw(1);
  for (int j = 1; j <= m; j++) {
    pw /= base;
    res *= (1 - pw);
  }
  return res;
}

inline mpq_class HermitianLocalFactorUnramified(HermitianLocalJordan const &hj,
                                                int64_t q) {
  if (hj.blocks.size() <= 1) {
    return mpq_class(1);
  }
  int64_t qK = 1;
  for (int i = 0; i < hj.f_K; i++) {
    qK *= q;
  }
  int m = 0;
  for (auto &b : hj.blocks) {
    m += b.second;
  }
  mpq_class f = HermitianGroupOrderNormalized(hj.type, m, qK);
  int64_t N = 0;
  int before = 0;
  for (size_t i = 0; i < hj.blocks.size(); i++) {
    int s = hj.blocks[i].first;
    int mi = hj.blocks[i].second;
    int after = m - before - mi;
    N += static_cast<int64_t>(s) * mi * (before - after);
    f /= HermitianGroupOrderNormalized(hj.type, mi, qK);
    before += mi;
  }
  mpq_class qpow(1);
  mpq_class qK_T(static_cast<long>(qK));
  for (int64_t i = 0; i < (N >= 0 ? N : -N); i++) {
    if (N >= 0) {
      qpow *= qK_T;
    } else {
      qpow /= qK_T;
    }
  }
  return qpow * f;
}

/*
  The ramified prime, d = q^a with q odd: E_P / Q_q is totally ramified, P
  is the only prime over q, its residue field is F_q and v_P(x) =
  v_q(N_{E/Q}(x)). The Hermitian form is computed in E exactly: an
  O_P-basis x_1, ..., x_m of L (x) Z_q is taken among the vectors
  independent in L / pi L (pi = 1 - zeta), and h(x, y) = sum_t b(zeta^t x,
  y) zeta_t^v with (zeta_t^v) the basis dual to (zeta^t) for Tr_{E/Q}. The
  Jordan splitting splits off a block of rank one when an entry of minimal
  valuation is on the diagonal and of rank two otherwise; the elimination
  coefficients are integral at P. A block of even scale s has a
  discriminant (-1)^{m_i (m_i - 1) / 2} det = N(pi)^{s m_i / 2} u, and it
  is a norm from E_P iff the residue of the unit u, its value at zeta = 1
  modulo q, is a square.
  The local factor is then that of Hecke.jl at a ramified non-dyadic
  prime, with the Sp and O^+- orders of the residue field F_q.
 */
using ElemE = std::vector<mpq_class>;

struct CyclotomicField {
  // d > 0: Q(zeta_d) on the powers of zeta; d = 0: a field Q(theta) on the
  // powers of theta, with Phi the minimal polynomial of theta and conj the
  // complex conjugate of theta on the same basis
  int64_t d;
  int deg;
  std::vector<int64_t> Phi;
  std::vector<mpq_class> conj;
  // v_P = v_q(N(a)) / val_div at a prime P alone over q, of residue
  // degree val_div
  int val_div = 1;
};

inline ElemE EReduce(std::vector<mpq_class> a, CyclotomicField const &K) {
  int deg = K.deg;
  for (int k = static_cast<int>(a.size()) - 1; k >= deg; k--) {
    mpq_class c = a[k];
    if (c == 0) {
      continue;
    }
    for (int j = 0; j <= deg; j++) {
      a[k - deg + j] -= c * mpq_class(static_cast<long>(K.Phi[j]));
    }
  }
  a.resize(deg, mpq_class(0));
  return a;
}

inline ElemE EMul(ElemE const &a, ElemE const &b, CyclotomicField const &K) {
  std::vector<mpq_class> c(a.size() + b.size(), mpq_class(0));
  for (size_t i = 0; i < a.size(); i++) {
    if (a[i] == 0) {
      continue;
    }
    for (size_t j = 0; j < b.size(); j++) {
      c[i + j] += a[i] * b[j];
    }
  }
  return EReduce(c, K);
}

inline ElemE ESub(ElemE a, ElemE const &b) {
  for (size_t i = 0; i < a.size(); i++) {
    a[i] -= b[i];
  }
  return a;
}

inline ElemE EZetaPower(int64_t t, CyclotomicField const &K) {
  if (K.d == 0) {
    // theta^t, t >= 0
    if (t < 0) {
      std::cerr << "CYCLO: a negative power of theta\n";
      throw TerminalException{1};
    }
    std::vector<mpq_class> c(std::max<int64_t>(t + 1, K.deg), mpq_class(0));
    c[t] = 1;
    return EReduce(c, K);
  }
  std::vector<mpq_class> c(K.d, mpq_class(0));
  c[((t % K.d) + K.d) % K.d] = 1;
  return EReduce(c, K);
}

inline ElemE EConj(ElemE const &a, CyclotomicField const &K) {
  if (K.d == 0) {
    // a(conj(theta)) by Horner
    ElemE res(K.deg, mpq_class(0));
    for (int i = K.deg - 1; i >= 0; i--) {
      res = EMul(res, K.conj, K);
      res[0] += a[i];
    }
    return res;
  }
  std::vector<mpq_class> c(K.d, mpq_class(0));
  for (int t = 0; t < K.deg; t++) {
    c[(K.d - t) % K.d] += a[t];
  }
  return EReduce(c, K);
}

inline MyMatrix<mpq_class> EMulMatrix(ElemE const &a, CyclotomicField const &K) {
  MyMatrix<mpq_class> M(K.deg, K.deg);
  for (int j = 0; j < K.deg; j++) {
    ElemE c = EMul(a, EZetaPower(j, K), K);
    for (int i = 0; i < K.deg; i++) {
      M(j, i) = c[i];
    }
  }
  return M;
}

inline bool EIsZero(ElemE const &a) {
  for (auto &x : a) {
    if (x != 0) {
      return false;
    }
  }
  return true;
}

inline ElemE EInv(ElemE const &a, CyclotomicField const &K) {
  // x a = 1: x M_a = e_0 with the rows of M_a the a zeta^j
  MyMatrix<mpq_class> M = EMulMatrix(a, K);
  MyMatrix<mpq_class> Minv = Inverse(M);
  ElemE x(K.deg);
  for (int i = 0; i < K.deg; i++) {
    x[i] = Minv(0, i);
  }
  return x;
}

// v_q(N(a)), large for zero
inline int EValuation(ElemE const &a, CyclotomicField const &K, int64_t q) {
  if (EIsZero(a)) {
    return 1 << 28;
  }
  return ValuationQRational(DeterminantMat(EMulMatrix(a, K)), q) / K.val_div;
}

inline int LegendreModQ(mpz_class a, int64_t q) {
  mpz_class qz(static_cast<long>(q));
  a %= qz;
  if (a < 0) {
    a += qz;
  }
  if (a == 0) {
    return 0;
  }
  mpz_class r;
  mpz_class e = (qz - 1) / 2;
  mpz_powm(r.get_mpz_t(), a.get_mpz_t(), e.get_mpz_t(), qz.get_mpz_t());
  return (r == 1) ? 1 : -1;
}

// the residue at zeta = 1 modulo q of a P-integral element
inline mpz_class EResidue(ElemE const &a, int64_t q) {
  mpz_class qz(static_cast<long>(q));
  mpz_class res(0);
  for (auto &x : a) {
    mpz_class num = x.get_num() % qz;
    mpz_class den = x.get_den() % qz;
    mpz_class inv;
    mpz_invert(inv.get_mpz_t(), den.get_mpz_t(), qz.get_mpz_t());
    res = (res + num * inv) % qz;
  }
  if (res < 0) {
    res += qz;
  }
  return res;
}

inline mpq_class HeckeGroupOrderSp(int m, int64_t q) {
  mpq_class res(1);
  for (int j = 2; j <= m; j += 2) {
    mpq_class pw(1);
    for (int i = 0; i < j; i++) {
      pw /= mpq_class(static_cast<long>(q));
    }
    res *= (1 - pw);
  }
  return res;
}

inline mpq_class HeckeGroupOrderO(int m, int eps, int64_t q) {
  auto qinv_pow = [&](int j) -> mpq_class {
    mpq_class pw(1);
    for (int i = 0; i < j; i++) {
      pw /= mpq_class(static_cast<long>(q));
    }
    return pw;
  };
  mpq_class res(2);
  if (m % 2 == 0) {
    int k = m / 2;
    res *= (1 - mpq_class(eps) * qinv_pow(k));
    for (int j = 2; j <= m - 2; j += 2) {
      res *= (1 - qinv_pow(j));
    }
  } else {
    for (int j = 2; j <= m; j += 2) {
      res *= (1 - qinv_pow(j));
    }
  }
  return res;
}

/*
  The Hermitian Gram matrix, in E, of an O_P-basis of L (x) Z_q at the
  ramified prime (d = q^a).
 */
struct HermitianLocalGram {
  CyclotomicField K;
  int m;
  std::vector<std::vector<ElemE>> H;
};

inline HermitianLocalGram
ComputeHermitianLocalGram(MyMatrix<mpq_class> const &B,
                          MyMatrix<mpq_class> const &S, int64_t d, int64_t q) {
  int r = B.rows();
  CyclotomicField K;
  K.d = d;
  K.Phi = CyclotomicPolynomial(static_cast<int>(d));
  K.deg = static_cast<int>(K.Phi.size()) - 1;
  int m = r / K.deg;
  // the dual basis of (zeta^t) for the trace: T_{st} = Tr(zeta^{s+t})
  auto trace = [&](ElemE const &x) -> mpq_class {
    // Tr(zeta^t) for t < deg, then linear
    mpq_class res(0);
    for (int t = 0; t < K.deg; t++) {
      if (x[t] == 0) {
        continue;
      }
      // Tr(zeta^t) is the Ramanujan sum c_d(t)
      mpq_class tr(0);
      int64_t g = CycloGcd(t, d);
      int64_t dd = d / g;
      // c_d(t) = mu(dd) phi(d) / phi(dd)
      int64_t mu = 1;
      {
        int64_t x2 = dd;
        for (int64_t p = 2; p <= x2; p++) {
          if (x2 % p == 0) {
            x2 /= p;
            if (x2 % p == 0) {
              mu = 0;
              break;
            }
            mu = -mu;
          }
        }
      }
      tr = mpq_class(static_cast<long>(mu * EulerPhi(d) / EulerPhi(dd)));
      res += x[t] * tr;
    }
    return res;
  };
  MyMatrix<mpq_class> T(K.deg, K.deg);
  for (int s2 = 0; s2 < K.deg; s2++) {
    for (int t = 0; t < K.deg; t++) {
      T(s2, t) = trace(EZetaPower(s2 + t, K));
    }
  }
  MyMatrix<mpq_class> Tinv = Inverse(T);
  // an O_P-basis among the standard vectors: independent modulo pi L
  MyMatrix<mpq_class> Id = IdentityMat<mpq_class>(r);
  MyMatrix<mpq_class> Pi = Id - S;
  std::vector<std::vector<int64_t>> rows;
  mpz_class qz(static_cast<long>(q));
  auto red = [&](mpq_class const &x) -> int64_t {
    mpz_class v = x.get_num() % qz;
    if (v < 0) {
      v += qz;
    }
    return v.get_si();
  };
  for (int i = 0; i < r; i++) {
    std::vector<int64_t> row(r);
    for (int j = 0; j < r; j++) {
      row[j] = red(Pi(i, j));
    }
    rows.push_back(row);
  }
  int base_rank = RankModQ(rows, q);
  std::vector<int> chosen;
  for (int k = 0; k < r && static_cast<int>(chosen.size()) < m; k++) {
    std::vector<std::vector<int64_t>> test = rows;
    for (int c : chosen) {
      std::vector<int64_t> e(r, 0);
      e[c] = 1;
      test.push_back(e);
    }
    std::vector<int64_t> e(r, 0);
    e[k] = 1;
    test.push_back(e);
    if (RankModQ(test, q) == base_rank + static_cast<int>(chosen.size()) + 1) {
      chosen.push_back(k);
    }
  }
  if (static_cast<int>(chosen.size()) != m) {
    std::cerr << "CYCLO: no O_P-basis found\n";
    throw TerminalException{1};
  }
  // h(x_i, x_j) = sum_t b(zeta^t x_i, x_j) zeta_t^v, zeta_t^v = sum_s
  // Tinv(t, s) zeta^s
  std::vector<MyMatrix<mpq_class>> Spow{Id};
  for (int t = 1; t < K.deg; t++) {
    Spow.push_back(Spow.back() * S);
  }
  std::vector<std::vector<ElemE>> H(m, std::vector<ElemE>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      ElemE val(K.deg, mpq_class(0));
      for (int t = 0; t < K.deg; t++) {
        // b(zeta^t x_i, x_j) with x_i the chosen standard vector
        mpq_class bt(0);
        for (int u = 0; u < r; u++) {
          bt += Spow[t](chosen[i], u) * B(u, chosen[j]);
        }
        for (int s2 = 0; s2 < K.deg; s2++) {
          val[s2] += bt * Tinv(t, s2);
        }
      }
      H[i][j] = val;
    }
  }
  return {K, m, H};
}

/*
  The Jordan splitting of a Hermitian form over E at the ramified prime P
  over q: a block of rank one when an entry of minimal valuation is on the
  diagonal, of rank two otherwise, the elimination coefficients being
  integral at P. Returns the blocks by scale (their v_P).
 */
using SmallGram = std::vector<std::vector<ElemE>>;

inline ElemE HermitianSmallDet(SmallGram const &G, CyclotomicField const &K) {
  if (G.size() == 1) {
    return G[0][0];
  }
  return ESub(EMul(G[0][0], G[1][1], K), EMul(G[0][1], G[1][0], K));
}

inline std::map<int, std::vector<SmallGram>>
HermitianJordanSplit(std::vector<std::vector<ElemE>> H,
                     CyclotomicField const &K, int64_t q) {
  std::map<int, std::vector<SmallGram>> blocks;
  while (!H.empty()) {
    int n = static_cast<int>(H.size());
    int vmin = 1 << 28, bi = -1, bj = -1;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        int v = EValuation(H[i][j], K, q);
        if (v < vmin || (v == vmin && i == j && bi != bj)) {
          vmin = v;
          bi = i;
          bj = j;
        }
      }
    }
    if (bi == bj) {
      ElemE hinv = EInv(H[bi][bi], K);
      for (int k = 0; k < n; k++) {
        if (k == bi) {
          continue;
        }
        ElemE c = EMul(H[k][bi], hinv, K);
        ElemE cc = EConj(c, K);
        for (int l = 0; l < n; l++) {
          H[k][l] = ESub(H[k][l], EMul(c, H[bi][l], K));
        }
        for (int l = 0; l < n; l++) {
          H[l][k] = ESub(H[l][k], EMul(H[l][bi], cc, K));
        }
      }
      blocks[vmin].push_back({{H[bi][bi]}});
      H.erase(H.begin() + bi);
      for (auto &row : H) {
        row.erase(row.begin() + bi);
      }
    } else {
      int i0 = bi, j0 = bj;
      ElemE a11 = H[i0][i0], a12 = H[i0][j0], a21 = H[j0][i0], a22 = H[j0][j0];
      ElemE det2 = ESub(EMul(a11, a22, K), EMul(a12, a21, K));
      ElemE dinv = EInv(det2, K);
      // the inverse of the 2x2 block: [[a22, -a12], [-a21, a11]] / det2
      for (int k = 0; k < n; k++) {
        if (k == i0 || k == j0) {
          continue;
        }
        // (c1, c2) = (H[k][i0], H[k][j0]) * block^{-1}
        ElemE x1 = H[k][i0], x2 = H[k][j0];
        ElemE c1 = EMul(ESub(EMul(x1, a22, K), EMul(x2, a21, K)), dinv, K);
        ElemE c2 = EMul(ESub(EMul(x2, a11, K), EMul(x1, a12, K)), dinv, K);
        ElemE cc1 = EConj(c1, K), cc2 = EConj(c2, K);
        for (int l = 0; l < n; l++) {
          H[k][l] = ESub(ESub(H[k][l], EMul(c1, H[i0][l], K)),
                         EMul(c2, H[j0][l], K));
        }
        for (int l = 0; l < n; l++) {
          H[l][k] = ESub(ESub(H[l][k], EMul(H[l][i0], cc1, K)),
                         EMul(H[l][j0], cc2, K));
        }
      }
      blocks[vmin].push_back({{a11, a12}, {a21, a22}});
      int hi = std::max(i0, j0), lo = std::min(i0, j0);
      H.erase(H.begin() + hi);
      H.erase(H.begin() + lo);
      for (auto &row : H) {
        row.erase(row.begin() + hi);
        row.erase(row.begin() + lo);
      }
    }
  }
  return blocks;
}

/*
  The dyadic ramified prime, d = 2^a with a >= 2, following the local
  factor of Hecke.jl there (QuadForm/Herm/Mass.jl): the factor of a maximal
  lattice when L is maximal integral at P, and otherwise the generic
  method. (Its formula for the scales 0 and 1 is never reached: the guard
  "g[2] != 0 || g[2] != 1" is always true.)

  The norm ideal n(L) of K is generated by the h(x_i, x_i) and the
  Tr_{E/K}(h(x_i, x_j) O_E), and v_p(Tr(c O_E)) = floor((v_P(c) + 2) / 2)
  since the different of E/K has exponent 2 at P; the valuations are kept
  in P-units. A lattice strictly containing L with the same norm level
  contains one of the L + R pi^{-1} y, y running over the lines of L / pi
  L, so L is maximal iff none of them keeps the norm level. For a maximal
  L of rank m the factor is 1/2 for m odd, and for m even 1 or (q^m - 1) /
  (2 q + 2) (q = 2) according as the discriminant (-1)^{m (m-1) / 2} det is
  a norm from E_P. The unit norms contain 1 + p^2 (the conductor of E_P /
  K_p is 2), so a unit u is a norm iff u = y y-bar modulo P^4 for a unit y
  of R / P^6.

  The generic method: a chain L = L_0 < L_1 < ... < L_k of lattices of
  index 2 with the norm level of L, L_k maximal; the factor is that of a
  maximal lattice of the same space times the product over i of the number
  of sublattices of index 2 of L_{i+1} isometric to L_i divided by the
  number of superlattices of index 2 of L_i isometric to L_{i+1}. The
  local isometry is decided by the genus symbol of Kirschmer (Hecke.jl,
  QuadForm/Herm/Legacy.jl): for the Jordan components in increasing scale
  S_i, the rank, S_i, whether the norm valuation of the component is S_i,
  the norm valuation v_i of the sum of the components scaled by
  (pi pi-bar)^{max(0, S_i - S_j)}, and the determinant of the first i
  components. Two lattices are isometric iff the first four agree, the
  last ratio of determinants is a norm and, for i < t, twice the normic
  defect of the i-th ratio (infinite for a norm, v_p(x) + 1 otherwise) is
  at least v_i + v_{i+1} - 2 S_i.
 */
struct DyadicContext {
  CyclotomicField K;
  ElemE one;
  ElemE pi;
  ElemE pi_inv;
  ElemE pipibar;
  ElemE pipibar_inv;
  std::vector<ElemE> pipow;
};

inline DyadicContext MakeDyadicContext(CyclotomicField const &K) {
  DyadicContext C;
  C.K = K;
  C.one = EZetaPower(0, K);
  C.pi = ESub(C.one, EZetaPower(1, K));
  C.pi_inv = EInv(C.pi, K);
  C.pipibar = EMul(C.pi, EConj(C.pi, K), K);
  C.pipibar_inv = EInv(C.pipibar, K);
  C.pipow = {C.one};
  for (int i = 1; i < 6; i++) {
    C.pipow.push_back(EMul(C.pipow.back(), C.pi, K));
  }
  return C;
}

inline int DyadicVal(ElemE const &x, DyadicContext const &C) {
  return EValuation(x, C.K, 2);
}

// in P-units: v_P of the ideal Tr(c O_E) O_E
inline int DyadicTraceVal(ElemE const &c, DyadicContext const &C) {
  int v = DyadicVal(c, C);
  if (v >= (1 << 27)) {
    return v;
  }
  int num = v + 2;
  int vp = (num >= 0) ? num / 2 : -((-num + 1) / 2);
  return 2 * vp;
}

inline int DyadicNormVal(std::vector<std::vector<ElemE>> const &G,
                         DyadicContext const &C) {
  int n = static_cast<int>(G.size());
  int res = 1 << 28;
  for (int i = 0; i < n; i++) {
    res = std::min(res, DyadicVal(G[i][i], C));
    for (int j = i + 1; j < n; j++) {
      res = std::min(res, DyadicTraceVal(G[i][j], C));
    }
  }
  return res;
}

// whether x (in K, nonzero) is a norm from E_P
inline bool DyadicIsLocalNorm(ElemE x, DyadicContext const &C) {
  CyclotomicField const &K = C.K;
  int w = DyadicVal(x, C);
  for (int i = 0; i < w / 2; i++) {
    x = EMul(x, C.pipibar_inv, K);
  }
  for (int i = 0; i < -w / 2; i++) {
    x = EMul(x, C.pipibar, K);
  }
  for (int64_t c = 0; c < 32; c++) {
    ElemE y = C.one;
    for (int i = 1; i < 6; i++) {
      if (c & (int64_t(1) << (i - 1))) {
        for (int t = 0; t < K.deg; t++) {
          y[t] += C.pipow[i][t];
        }
      }
    }
    if (DyadicVal(ESub(x, EMul(y, EConj(y, K), K)), C) >= 4) {
      return true;
    }
  }
  return false;
}

// U H U^*, the rows of U being the new basis in the old one
inline std::vector<std::vector<ElemE>>
HermitianCongruence(std::vector<std::vector<ElemE>> const &H,
                    std::vector<std::vector<ElemE>> const &U,
                    CyclotomicField const &K) {
  int m = static_cast<int>(H.size());
  std::vector<std::vector<ElemE>> UH(m, std::vector<ElemE>(m));
  for (int a = 0; a < m; a++) {
    for (int k = 0; k < m; k++) {
      ElemE sum(K.deg, mpq_class(0));
      for (int i = 0; i < m; i++) {
        if (!EIsZero(U[a][i])) {
          ElemE t = EMul(U[a][i], H[i][k], K);
          for (int u = 0; u < K.deg; u++) {
            sum[u] += t[u];
          }
        }
      }
      UH[a][k] = sum;
    }
  }
  std::vector<std::vector<ElemE>> R(m, std::vector<ElemE>(m));
  for (int a = 0; a < m; a++) {
    for (int b = 0; b < m; b++) {
      ElemE sum(K.deg, mpq_class(0));
      for (int k = 0; k < m; k++) {
        if (!EIsZero(U[b][k])) {
          ElemE t = EMul(UH[a][k], EConj(U[b][k], K), K);
          for (int u = 0; u < K.deg; u++) {
            sum[u] += t[u];
          }
        }
      }
      R[a][b] = sum;
    }
  }
  return R;
}

// L + R pi^{-1} y, y = sum_{i in mask} x_i
inline std::vector<std::vector<ElemE>>
DyadicOverlattice(std::vector<std::vector<ElemE>> const &H, int64_t mask,
                  DyadicContext const &C) {
  int m = static_cast<int>(H.size());
  ElemE zero(C.K.deg, mpq_class(0));
  std::vector<std::vector<ElemE>> U(m, std::vector<ElemE>(m, zero));
  int j = -1;
  for (int i = 0; i < m; i++) {
    U[i][i] = C.one;
    if (j < 0 && (mask & (int64_t(1) << i))) {
      j = i;
    }
  }
  for (int i = 0; i < m; i++) {
    U[j][i] = (mask & (int64_t(1) << i)) ? C.pi_inv : zero;
  }
  return HermitianCongruence(H, U, C.K);
}

// the kernel of the functional phi (mask) on L / pi L
inline std::vector<std::vector<ElemE>>
DyadicSublattice(std::vector<std::vector<ElemE>> const &H, int64_t mask,
                 DyadicContext const &C) {
  int m = static_cast<int>(H.size());
  ElemE zero(C.K.deg, mpq_class(0));
  std::vector<std::vector<ElemE>> U(m, std::vector<ElemE>(m, zero));
  int j = -1;
  for (int i = 0; i < m; i++) {
    if (mask & (int64_t(1) << i)) {
      j = i;
      break;
    }
  }
  for (int i = 0; i < m; i++) {
    if (i == j) {
      U[i][j] = C.pi;
    } else {
      U[i][i] = C.one;
      if (mask & (int64_t(1) << i)) {
        U[i][j] = ESub(zero, C.one);
      }
    }
  }
  return HermitianCongruence(H, U, C.K);
}

struct KirschmerSymbolEntry {
  int rank;
  int scale;
  bool normal;
  int v;
  ElemE det;
};

inline std::vector<KirschmerSymbolEntry>
KirschmerSymbol(std::vector<std::vector<ElemE>> const &H, DyadicContext const &C) {
  std::map<int, std::vector<SmallGram>> blocks = HermitianJordanSplit(H, C.K, 2);
  std::vector<int> scales;
  std::vector<int> nvals;
  std::vector<int> ranks;
  std::vector<ElemE> dets;
  for (auto &kv : blocks) {
    int nv = 1 << 28;
    int rk = 0;
    ElemE det = C.one;
    for (auto &G : kv.second) {
      nv = std::min(nv, DyadicNormVal(G, C));
      rk += static_cast<int>(G.size());
      det = EMul(det, HermitianSmallDet(G, C.K), C.K);
    }
    scales.push_back(kv.first);
    nvals.push_back(nv);
    ranks.push_back(rk);
    dets.push_back(det);
  }
  std::vector<KirschmerSymbolEntry> sym;
  ElemE det_partial = C.one;
  size_t t = scales.size();
  for (size_t i = 0; i < t; i++) {
    int v = 1 << 28;
    for (size_t j = 0; j < t; j++) {
      int shift = std::max(0, scales[i] - scales[j]);
      v = std::min(v, nvals[j] + 2 * shift);
    }
    det_partial = EMul(det_partial, dets[i], C.K);
    sym.push_back({ranks[i], scales[i], nvals[i] == scales[i], v, det_partial});
  }
  return sym;
}

inline bool DyadicIsLocallyIsometric(std::vector<KirschmerSymbolEntry> const &S1,
                                     std::vector<KirschmerSymbolEntry> const &S2,
                                     DyadicContext const &C) {
  size_t t = S1.size();
  if (t != S2.size()) {
    return false;
  }
  for (size_t i = 0; i < t; i++) {
    if (S1[i].rank != S2[i].rank || S1[i].scale != S2[i].scale ||
        S1[i].normal != S2[i].normal || S1[i].v != S2[i].v) {
      return false;
    }
  }
  auto ratio = [&](size_t i) -> ElemE {
    return EMul(S1[i].det, EInv(S2[i].det, C.K), C.K);
  };
  if (!DyadicIsLocalNorm(ratio(t - 1), C)) {
    return false;
  }
  for (size_t i = 0; i + 1 < t; i++) {
    ElemE x = ratio(i);
    if (DyadicIsLocalNorm(x, C)) {
      continue;
    }
    // normic defect v_p(x) + v_p(disc(O_E / O_K)) - 1, the discriminant
    // having valuation 2 at p
    int defect = DyadicVal(x, C) / 2 + 2 - 1;
    if (2 * defect < S1[i].v + S1[i + 1].v - 2 * S1[i].scale) {
      return false;
    }
  }
  return true;
}

inline std::optional<mpq_class>
HermitianLocalFactorDyadic(MyMatrix<mpq_class> const &B,
                           MyMatrix<mpq_class> const &S, int64_t d) {
  int64_t q = 2;
  HermitianLocalGram hlg = ComputeHermitianLocalGram(B, S, d, q);
  DyadicContext C = MakeDyadicContext(hlg.K);
  CyclotomicField const &K = C.K;
  int m = hlg.m;
  int64_t n_mask = int64_t(1) << m;
  // the factor of a maximal lattice of the same space
  auto maximal_factor = [&](std::vector<std::vector<ElemE>> const &H) -> mpq_class {
    if (m % 2 == 1) {
      return mpq_class(1, 2);
    }
    std::vector<std::vector<ElemE>> G = H;
    ElemE det = C.one;
    for (int c = 0; c < m; c++) {
      int piv = -1;
      for (int i = c; i < m; i++) {
        if (!EIsZero(G[i][c])) {
          piv = i;
          break;
        }
      }
      if (piv < 0) {
        std::cerr << "CYCLO: degenerate Hermitian form\n";
        throw TerminalException{1};
      }
      if (piv != c) {
        std::swap(G[piv], G[c]);
        for (auto &x : det) {
          x = -x;
        }
      }
      det = EMul(det, G[c][c], K);
      ElemE inv = EInv(G[c][c], K);
      for (int i = c + 1; i < m; i++) {
        ElemE t = EMul(G[i][c], inv, K);
        for (int j = c; j < m; j++) {
          G[i][j] = ESub(G[i][j], EMul(t, G[c][j], K));
        }
      }
    }
    if ((m * (m - 1) / 2) % 2 == 1) {
      for (auto &x : det) {
        x = -x;
      }
    }
    if (DyadicIsLocalNorm(det, C)) {
      return mpq_class(1);
    }
    mpq_class qm(1);
    for (int i = 0; i < m; i++) {
      qm *= 2;
    }
    return (qm - 1) / mpq_class(static_cast<long>(2 * q + 2));
  };
  int nL = DyadicNormVal(hlg.H, C);
  // an overlattice of index 2 keeping the norm level, if any
  auto find_over = [&](std::vector<std::vector<ElemE>> const &H)
      -> std::optional<std::vector<std::vector<ElemE>>> {
    for (int64_t mask = 1; mask < n_mask; mask++) {
      std::vector<std::vector<ElemE>> Hn = DyadicOverlattice(H, mask, C);
      if (DyadicNormVal(Hn, C) >= nL) {
        return Hn;
      }
    }
    return {};
  };
  std::vector<std::vector<std::vector<ElemE>>> chain{hlg.H};
  while (true) {
    std::optional<std::vector<std::vector<ElemE>>> opt = find_over(chain.back());
    if (!opt) {
      break;
    }
    chain.push_back(*opt);
    if (chain.size() > 64) {
      std::cerr << "CYCLO: the chain of overlattices does not end\n";
      throw TerminalException{1};
    }
  }
  mpq_class f = maximal_factor(hlg.H);
  for (size_t i = 0; i + 1 < chain.size(); i++) {
    std::vector<KirschmerSymbolEntry> sym_i = KirschmerSymbol(chain[i], C);
    std::vector<KirschmerSymbolEntry> sym_n = KirschmerSymbol(chain[i + 1], C);
    int64_t A = 0, Bc = 0;
    for (int64_t mask = 1; mask < n_mask; mask++) {
      if (DyadicIsLocallyIsometric(
              KirschmerSymbol(DyadicSublattice(chain[i + 1], mask, C), C), sym_i,
              C)) {
        A++;
      }
      if (DyadicIsLocallyIsometric(
              KirschmerSymbol(DyadicOverlattice(chain[i], mask, C), C), sym_n,
              C)) {
        Bc++;
      }
    }
    if (A == 0 || Bc == 0) {
      std::cerr << "CYCLO: a step of the chain has no isometric neighbour\n";
      throw TerminalException{1};
    }
    f *= mpq_class(static_cast<long>(A), static_cast<long>(Bc));
  }
  return f;
}

/*
  The local factor at the ramified prime.
 */
inline std::optional<mpq_class>
HermitianLocalFactorRamified(MyMatrix<mpq_class> const &B,
                             MyMatrix<mpq_class> const &S, int64_t d,
                             int64_t q) {
  if (q == 2) {
    return HermitianLocalFactorDyadic(B, S, d);
  }
  HermitianLocalGram hlg = ComputeHermitianLocalGram(B, S, d, q);
  CyclotomicField const &K = hlg.K;
  int m = hlg.m;
  std::vector<std::vector<ElemE>> H = hlg.H;
  // the Jordan splitting, the blocks of a scale gathered with their
  // determinant
  std::map<int, std::pair<int, ElemE>> blocks;  // scale -> (rank, det)
  for (auto &kv : HermitianJordanSplit(H, K, q)) {
    int rk = 0;
    ElemE det = EZetaPower(0, K);
    for (auto &G : kv.second) {
      rk += static_cast<int>(G.size());
      det = EMul(det, HermitianSmallDet(G, K), K);
    }
    blocks[kv.first] = {rk, det};
  }
  // the local factor of Hecke.jl at a ramified non-dyadic prime
  mpq_class f = HeckeGroupOrderSp(2 * (m / 2), q);
  int64_t N = 0;
  int after_total = m;
  ElemE pipibar = ESub(ESub(EZetaPower(0, K), EZetaPower(1, K)),
                       ESub(EZetaPower(-1, K), EZetaPower(0, K)));
  for (auto &kv : blocks) {
    int s = kv.first;
    int mi = kv.second.first;
    after_total -= mi;
    int ri = after_total;
    int64_t fl = (s >= 0) ? s / 2 : -((-s + 1) / 2);
    N -= fl * mi * mi;
    if (s % 2 != 0) {
      N -= static_cast<int64_t>(mi + 1) * (mi / 2);
      f /= HeckeGroupOrderSp(mi, q);
    } else {
      // det = (pi pi-bar)^{s mi / 2} u
      ElemE u = kv.second.second;
      int w = EValuation(u, K, q);
      ElemE pinv = EInv(pipibar, K);
      for (int i = 0; i < w / 2; i++) {
        u = EMul(u, pinv, K);
      }
      // the discriminant (-1)^{m_i (m_i - 1) / 2} det, as in Hecke.jl
      mpz_class res_u = EResidue(u, q);
      if ((mi * (mi - 1) / 2) % 2 == 1) {
        res_u = -res_u;
      }
      int leg = LegendreModQ(res_u, q);
      f /= HeckeGroupOrderO(mi, leg, q);
    }
    N -= static_cast<int64_t>(s) * mi * ri;
    N += static_cast<int64_t>(s) * mi * m / 2;
  }
  if (m % 2 == 0) {
    N += m / 2;
  }
  mpq_class qpow(1);
  for (int64_t i = 0; i < (N >= 0 ? N : -N); i++) {
    if (N >= 0) {
      qpow *= mpq_class(static_cast<long>(q));
    } else {
      qpow /= mpq_class(static_cast<long>(q));
    }
  }
  return qpow * f;
}

/*
  The mass of the genus of the Hermitian lattice over Q(zeta_d) given by
  the Gram matrix B of its trace form and the action S of zeta_d, or
  nothing at a ramified prime of E/K with rank >= 2.
 */
inline std::optional<mpq_class>
HermitianMassCyclotomicLattice(MyMatrix<mpq_class> const &B,
                               MyMatrix<mpq_class> const &S, int64_t d_in) {
  int64_t d = NormalizeCyclotomicIndex(d_in);
  int r = B.rows();
  int64_t phi = EulerPhi(d);
  int m = static_cast<int>(r / phi);
  if (m == 1) {
    return RankOneHermitianMassCyclotomic(d);
  }
  // For d_in = 2 d with d odd, Phi_{d_in}(S) = Phi_d(-S): -S is the action
  // of a primitive d-th root of unity.
  MyMatrix<mpq_class> Sd = S;
  if (d != d_in) {
    Sd = -S;
  }
  mpq_class res = StandardHermitianMassCyclotomic(d, m);
  mpq_class det = DeterminantMat(B);
  std::set<int64_t> primes;
  {
    mpz_class num = det.get_num();
    if (num < 0) {
      num = -num;
    }
    if (num > 1) {
      for (mpz_class const &p : FactorsInt(num)) {
        primes.insert(p.get_si());
      }
    }
    // the primes dividing d, where E/Q ramifies
    int64_t x = d;
    for (int64_t p = 2; p <= x; p++) {
      if (x % p == 0) {
        primes.insert(p);
        while (x % p == 0) {
          x /= p;
        }
      }
    }
  }
  for (int64_t q : primes) {
    for (auto &hj : HermitianJordanAtPrime(B, Sd, d, q)) {
      if (hj.type == 2) {
        std::optional<mpq_class> opt = HermitianLocalFactorRamified(B, Sd, d, q);
        if (!opt) {
          return {};
        }
        res *= *opt;
      } else {
        res *= HermitianLocalFactorUnramified(hj, q);
      }
    }
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_CYCLOTOMICHERMITIANMASS_H_
// clang-format on
