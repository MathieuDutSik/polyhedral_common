// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_HERMITIANGENUSSEED_H_
#define SRC_GENERA_HERMITIANGENUSSEED_H_

// clang-format off
#include "HermitianGenusSymbol.h"
#include <deque>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_HERMITIAN_GENUS_SEED
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_HERMITIAN_GENUS_SEED
#endif

/*
  A lattice in a genus of positive definite Hermitian lattices over the
  order O_E of a CM field, from its symbol (HermitianGenusSymbol.h), as
  representative(G) of Hecke.jl, at every prime (with the norm valuations
  at the dyadic ramified ones).

  1. The space. A Hermitian space of dimension m is determined, at a finite
     prime p of K_0, by the norm class of its determinant, and is definite
     at every real place when it is diag(1, ..., 1, a) with a totally
     positive. So a is searched among the small totally positive elements
     of O_{K_0} (the short vectors of Tr_{K_0}(x^2) on O_{K_0}, the fixed
     points of the conjugation) for a norm class -1 exactly at the primes
     where the symbol has its determinant -1 (the split primes have none);
     the product formula makes the number of those primes even, which the
     symbol must satisfy.
  2. A maximal lattice. From O_E^m with that form, the integral overlattices
     M + P^{-1} x of index N(P) are taken while there are some, at the
     primes above the divisors of 2 det, integral being for the norm: h(x,
     x) in O_{K_0} for all x, that is h(z_k, z_k) and Tr_{E/K_0} h(z_k, z_l)
     in O_{K_0} on a Z-basis. The norm-maximal lattices of a space are
     locally isometric (Kirschmer, Habilitation, Lemma 3.5.2), so every
     integral lattice of the space is, locally, isometric to a sublattice of
     M.
  3. The descent. At a prime p of K_0 where the local symbol of M is not
     that of the genus, a breadth first search over the sublattices of
     index N(P) (the kernels of the nonzero O_E-linear maps M -> O_E / P, at
     P and at conj(P)), one lattice being kept per local symbol: two
     lattices with the same symbol at p are locally isometric, so they have
     the same sublattices up to isometry and the search on the symbols is
     complete. The volume v_P(det) grows by one at each step, which bounds
     the search by that of the genus. The descent at p does not change the
     lattice at the other primes.
  The result is checked against the symbol at every prime.
 */

// The primes of O_E above each ell, kept so that their indices stay fixed.
template <typename T, typename Tint> struct HermPrimeCache {
  CMOrder<T, Tint> const *ord;
  std::map<int64_t, std::vector<HermIdealPrime<Tint>>> primes;
  std::vector<HermIdealPrime<Tint>> const &at(int64_t ell) {
    auto iter = primes.find(ell);
    if (iter == primes.end()) {
      iter = primes.emplace(ell, HermIdealPrimesAbove(*ord, ell)).first;
    }
    return iter->second;
  }
};

// A prime of K_0: the rational prime and the index of a prime of E above it.
struct HermK0Ref {
  int64_t ell;
  size_t i_P;
  bool operator<(HermK0Ref const &o) const {
    return std::make_pair(ell, i_P) < std::make_pair(o.ell, o.i_P);
  }
  bool operator==(HermK0Ref const &o) const {
    return ell == o.ell && i_P == o.i_P;
  }
};

// The representative i_P of the prime of K_0 under the prime j of E.
template <typename Tint>
size_t HermK0Representative(std::vector<HermIdealPrime<Tint>> const &l_prime,
                            size_t j) {
  return std::min(j, l_prime[j].conj_index);
}

/*
  The parsed symbol: for each prime, ell, the element named (theta
  coordinates, empty for ell alone) and the local symbol.
 */
template <typename T> struct HermSymbolInput {
  int m;
  std::vector<int64_t> l_ell;
  std::vector<std::vector<T>> l_elt;
  std::vector<HermLocalSymbol> l_sym;
};

inline std::string HermTrim(std::string const &s) {
  size_t a = s.find_first_not_of(" \t\n");
  if (a == std::string::npos) {
    return "";
  }
  size_t b = s.find_last_not_of(" \t\n");
  return s.substr(a, b - a + 1);
}

template <typename T>
HermSymbolInput<T> ParseHermSymbol(std::string const &str) {
  HermSymbolInput<T> res;
  std::vector<std::string> parts = STRING_Split(str, ";");
  if (parts.empty()) {
    std::cerr << "HERMSEED: empty symbol\n";
    throw TerminalException{1};
  }
  std::string head = HermTrim(parts[0]);
  if (head.substr(0, 2) != "m=") {
    std::cerr << "HERMSEED: the symbol must start with m=<rank>\n";
    throw TerminalException{1};
  }
  res.m = ParseScalar<int>(head.substr(2));
  for (size_t i = 1; i < parts.size(); i++) {
    std::string part = HermTrim(parts[i]);
    if (part.empty()) {
      continue;
    }
    // the prime name ends at the ':' following its closing parenthesis
    size_t pos_colon = part.find(':', part[0] == '(' ? part.find(')') : 0);
    std::string name = HermTrim(part.substr(0, pos_colon));
    std::string blocks = HermTrim(part.substr(pos_colon + 1));
    std::vector<T> elt;
    int64_t ell;
    if (name[0] == '(') {
      size_t pos_comma = name.find(',');
      ell = ParseScalar<int64_t>(name.substr(1, pos_comma - 1));
      size_t a = name.find('['), b = name.find(']');
      for (auto &c : STRING_Split(name.substr(a + 1, b - a - 1), ",")) {
        elt.push_back(ParseScalar<T>(HermTrim(c)));
      }
    } else {
      ell = ParseScalar<int64_t>(name);
    }
    HermLocalSymbol sym;
    size_t pos = 0;
    while ((pos = blocks.find('(', pos)) != std::string::npos) {
      size_t end = blocks.find(')', pos);
      std::vector<std::string> f =
          STRING_Split(blocks.substr(pos + 1, end - pos - 1), ",");
      if (f.size() != 3 && f.size() != 4) {
        std::cerr << "HERMSEED: a block is (s,r,+) or (s,r,-), with the "
                     "norm valuation (s,r,+,n) at a dyadic ramified prime\n";
        throw TerminalException{1};
      }
      std::string sign = HermTrim(f[2]);
      sym.blocks.push_back({ParseScalar<int>(HermTrim(f[0])),
                            ParseScalar<int>(HermTrim(f[1])),
                            sign == "+" ? 1 : -1});
      if (f.size() == 4) {
        sym.norms.push_back(ParseScalar<int>(HermTrim(f[3])));
      }
      pos = end + 1;
    }
    res.l_ell.push_back(ell);
    res.l_elt.push_back(elt);
    res.l_sym.push_back(sym);
  }
  return res;
}

// The prime of K_0 named by ell and the element elt (theta coordinates).
template <typename T, typename Tint>
HermK0Ref HermResolvePrime(HermPrimeCache<T, Tint> &cache, int64_t ell,
                           std::vector<T> const &elt) {
  CMOrder<T, Tint> const &ord = *cache.ord;
  std::vector<HermIdealPrime<Tint>> const &l_prime = cache.at(ell);
  std::vector<size_t> l_rep;
  for (size_t j = 0; j < l_prime.size(); j++) {
    if (HermK0Representative(l_prime, j) == j) {
      l_rep.push_back(j);
    }
  }
  if (elt.empty()) {
    if (l_rep.size() != 1) {
      std::cerr << "HERMSEED: several primes of K_0 above " << ell
                << ", one must be named by an element\n";
      throw TerminalException{1};
    }
    return {ell, l_rep[0]};
  }
  MyVector<T> th(ord.deg);
  for (int k = 0; k < ord.deg; k++) {
    th(k) = elt[k];
  }
  MyVector<T> a_e = Inverse(ord.basis_theta).transpose() * th;
  for (size_t j : l_rep) {
    if (HermResidueRat(l_prime[j], a_e).is_zero()) {
      return {ell, j};
    }
  }
  std::cerr << "HERMSEED: the element names no prime above " << ell << "\n";
  throw TerminalException{1};
}

/*
  The equality of local symbols, as == of Hecke.jl: the same scales and
  ranks and the same determinant; at a ramified prime that is not dyadic,
  also the same determinants of the blocks of even scale (those of odd
  scale are hyperbolic); at a dyadic ramified prime (Kirschmer, Theorem
  3.3.6), the same normal blocks (scale twice the norm), the same
  invariants n_i = min_j 2 max(0, s_i - s_j) + 2 n_j, and where the partial
  determinants differ 2 (e - 1) >= n_i + n_{i+1} - 2 s_i, e the valuation
  of the different (diff_val).
 */
inline bool HermLocalSymbolEquivalent(HermLocalSymbol const &a,
                                      HermLocalSymbol const &b, bool ramified,
                                      bool dyadic = false, int diff_val = 0) {
  size_t t = a.blocks.size();
  if (t != b.blocks.size()) {
    return false;
  }
  for (size_t i = 0; i < t; i++) {
    if (a.blocks[i][0] != b.blocks[i][0] || a.blocks[i][1] != b.blocks[i][1]) {
      return false;
    }
  }
  if (a.det() != b.det()) {
    return false;
  }
  if (!ramified) {
    return true;
  }
  if (!dyadic) {
    for (size_t i = 0; i < t; i++) {
      if (a.blocks[i][0] % 2 == 0 && a.blocks[i][2] != b.blocks[i][2]) {
        return false;
      }
    }
    return true;
  }
  if (a.norms.size() != t || b.norms.size() != t) {
    std::cerr << "HERMSEED: a dyadic ramified symbol needs the norm "
                 "valuations\n";
    throw TerminalException{1};
  }
  // the dyadic conditions of HermitianGenus.h
  return HermDyadicEquivalent(a, b, diff_val);
}

inline int HermLocalSymbolVolume(HermLocalSymbol const &a) {
  int v = 0;
  for (auto const &bl : a.blocks) {
    v += bl[0] * bl[1];
  }
  return v;
}

// Whether the lattice of basis A has its norm in O_{K_0}.
template <typename T, typename Tint>
bool HermIsNormIntegral(HermLatticeData<T, Tint> const &L) {
  int n = L.n();
  MyMatrix<T> C = UniversalMatrixConversion<T, Tint>(L.ord->conj).transpose();
  for (int k = 0; k < n; k++) {
    MyVector<Tint> ek = ZeroVector<Tint>(n);
    ek(k) = 1;
    if (!IsIntegralVector(L.herm(ek, ek))) {
      return false;
    }
    for (int l = k + 1; l < n; l++) {
      MyVector<Tint> el = ZeroVector<Tint>(n);
      el(l) = 1;
      MyVector<T> h = L.herm(ek, el);
      if (!IsIntegralVector(MyVector<T>(h + C * h))) {
        return false;
      }
    }
  }
  return true;
}

// The lattice M + I x for x in M and I a fractional ideal (rows on e), in
// the coordinates of M.
template <typename T, typename Tint>
MyMatrix<T> HermAddIdealTimes(HermLatticeData<T, Tint> const &L,
                              MyMatrix<T> const &I, MyVector<Tint> const &x) {
  int n = L.n();
  int d = L.ord->deg;
  MyMatrix<T> Gens(n + I.rows(), n);
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
      Gens(i, k) = (i == k) ? T(1) : T(0);
    }
  }
  MyVector<T> x_T = UniversalVectorConversion<T, Tint>(x);
  for (int r = 0; r < I.rows(); r++) {
    MyMatrix<T> Mc = ZeroMatrix<T>(n, n);
    for (int a = 0; a < d; a++) {
      if (I(r, a) != 0) {
        Mc += I(r, a) * UniversalMatrixConversion<T, Tint>(L.Omega[a]);
      }
    }
    MyVector<T> cx = Mc.transpose() * x_T;
    for (int k = 0; k < n; k++) {
      Gens(n + r, k) = cx(k);
    }
  }
  return GetZbasis(Gens);
}

/*
  The kernel of the O_E-linear map M -> O_E / P whose coordinates on the
  F_q-basis of the line space are c (codes), in the coordinates of M.
 */
template <typename T, typename Tint>
MyMatrix<T> HermHyperplaneSublattice(HermLatticeData<T, Tint> const &L,
                                     HermLineSpace const &ls,
                                     std::vector<int> const &c) {
  int n = L.n();
  int64_t ell = ls.ell();
  int f = ls.F.f;
  std::vector<std::vector<int64_t>> func;
  for (int k = 0; k < n; k++) {
    std::vector<int64_t> z(n, 0);
    z[k] = 1;
    std::vector<uint64_t> coord = ls.coordinates(z);
    uint64_t val = 0;
    for (int j = 0; j < ls.m; j++) {
      val = FqCodeAdd(ls.F, val,
                      FqCodeMul(ls.F, coord[j], static_cast<uint64_t>(c[j])));
    }
    std::vector<int64_t> row(f, 0);
    for (int t = 0; t < f; t++) {
      row[t] = static_cast<int64_t>(val % static_cast<uint64_t>(ell));
      val /= static_cast<uint64_t>(ell);
    }
    func.push_back(row);
  }
  Tint ell_T = UniversalScalarConversion<Tint, int64_t>(ell);
  MyMatrix<Tint> Ker = NullspaceMatMod(HermRowsToMatrix<Tint>(func, f), ell_T);
  MyMatrix<T> Gens(Ker.rows() + n, n);
  for (int i = 0; i < Ker.rows(); i++) {
    for (int k = 0; k < n; k++) {
      Gens(i, k) = UniversalScalarConversion<T, Tint>(Ker(i, k));
    }
  }
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
      Gens(Ker.rows() + i, k) =
          (i == k) ? UniversalScalarConversion<T, int64_t>(ell) : T(0);
    }
  }
  return GetZbasis(Gens);
}

/*
  The seed of the genus of the symbol: the realified ring of the space and
  the basis of the lattice in it.
 */
template <typename T, typename Tint> struct HermSeed {
  RealifiedRing<T, Tint> ring;
  MyMatrix<T> basis;
};

template <typename T, typename Tint>
HermSeed<T, Tint> HermSeedFromSymbol(CMOrder<T, Tint> const &ord,
                                     std::string const &str, std::ostream &os) {
  int d = ord.deg;
  HermSymbolInput<T> input = ParseHermSymbol<T>(str);
  int m = input.m;
  HermPrimeCache<T, Tint> cache;
  cache.ord = &ord;
  // the symbol, by prime of K_0
  std::map<HermK0Ref, HermLocalSymbol> target;
  std::set<HermK0Ref> non_norm;
  for (size_t i = 0; i < input.l_ell.size(); i++) {
    HermK0Ref ref = HermResolvePrime(cache, input.l_ell[i], input.l_elt[i]);
    HermIdealPrime<Tint> const &P = cache.at(ref.ell)[ref.i_P];
    int rk = 0;
    for (auto const &b : input.l_sym[i].blocks) {
      rk += b[1];
    }
    if (rk != m) {
      std::cerr << "HERMSEED: a local symbol is not of rank " << m << "\n";
      throw TerminalException{1};
    }
    target[ref] = input.l_sym[i];
    if (input.l_sym[i].det() == -1 && P.type != HermIdealType::Split) {
      non_norm.insert(ref);
    }
  }
  if (non_norm.size() % 2 != 0) {
    std::cerr << "HERMSEED: the determinant is not a norm at an odd number of "
                 "primes, there is no such positive definite lattice\n";
    throw TerminalException{1};
  }
  // 1. the space diag(1, ..., 1, a)
  MyMatrix<Tint> Cm = ord.conj - IdentityMat<Tint>(d);
  MyMatrix<Tint> B0 = NullspaceIntMat(Cm);
  MyMatrix<Tint> Q = ord.trace_norm_form();
  MyMatrix<T> Q0 = UniversalMatrixConversion<T, Tint>(MyMatrix<Tint>(B0 * Q * B0.transpose()));
  Tint disc = DeterminantMat(Q);
  auto prime_divisors = [&](Tint x) -> std::vector<int64_t> {
    std::vector<int64_t> res;
    x = T_abs(x);
    for (int64_t p = 2;; p++) {
      Tint p_T = UniversalScalarConversion<Tint, int64_t>(p);
      if (p_T * p_T > x) {
        break;
      }
      if (x % p_T == 0) {
        res.push_back(p);
        while (x % p_T == 0) {
          x /= p_T;
        }
      }
    }
    if (x > 1) {
      res.push_back(UniversalScalarConversion<int64_t, Tint>(x));
    }
    return res;
  };
  std::optional<MyVector<Tint>> opt_a;
  CVPSolver<T, Tint> solver(Q0, os);
  for (int bound = 4 * d; !opt_a && bound < 1000000; bound *= 2) {
    for (auto const &v : solver.at_most_norm_vectors(T(bound))) {
      for (int sgn : {1, -1}) {
        if (opt_a) {
          break;
        }
        MyVector<Tint> a = Tint(sgn) * (B0.transpose() * v);
        // totally positive: Tr(a y^2) > 0 on K_0
        MyMatrix<Tint> Ta = B0 * ord.mult_matrix(a) * Q * B0.transpose();
        // Q(b_i, b_j) = Tr(b_i conj(b_j)) = Tr(b_i b_j) on K_0, so Ta is the
        // form Tr(a y z)
        if (!IsPositiveDefinite(UniversalMatrixConversion<T, Tint>(Ta), os)) {
          continue;
        }
        Tint nrm = DeterminantMat(ord.mult_matrix(a));
        std::set<int64_t> l_ell;
        for (int64_t p : prime_divisors(nrm)) {
          l_ell.insert(p);
        }
        for (int64_t p : prime_divisors(disc)) {
          l_ell.insert(p);
        }
        for (auto const &ref : non_norm) {
          l_ell.insert(ref.ell);
        }
        bool ok = true;
        std::set<HermK0Ref> found;
        MyVector<T> a_T = UniversalVectorConversion<T, Tint>(a);
        for (int64_t ell : l_ell) {
          std::vector<HermIdealPrime<Tint>> const &l_prime = cache.at(ell);
          for (size_t j = 0; j < l_prime.size() && ok; j++) {
            if (HermK0Representative(l_prime, j) != j) {
              continue;
            }
            if (HermLocalNormClass(ord, l_prime, j, a_T) == -1) {
              found.insert({ell, j});
            }
          }
        }
        if (ok && found == non_norm) {
          opt_a = a;
        }
      }
      if (opt_a) {
        break;
      }
    }
  }
  if (!opt_a) {
    std::cerr << "HERMSEED: no element a for the space found\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> h(d, ZeroMatrix<T>(m, m));
  for (int i = 0; i < m - 1; i++) {
    h[0](i, i) = 1;
  }
  MyVector<T> a_th =
      ord.basis_theta.transpose() * UniversalVectorConversion<T, Tint>(*opt_a);
  for (int k = 0; k < d; k++) {
    h[k](m - 1, m - 1) = a_th(k);
  }
  HermSeed<T, Tint> seed;
  seed.ring = CMOrderRealifiedRing<T, Tint>(ord, m, h);
  int n = seed.ring.n();
  seed.basis = IdentityMat<T>(n);
#ifdef DEBUG_HERMITIAN_GENUS_SEED
  os << "HERMSEED: the space diag(1, ..., 1, a), a = " << StringVector(a_th)
     << "\n";
#endif
  auto lattice_of = [&](MyMatrix<T> const &A) -> HermLatticeData<T, Tint> {
    return HermLatticeOfBasis(ord, seed.ring, A);
  };
  // the primes to treat: those of the symbol and those above 2 det
  std::set<int64_t> l_ell_all;
  for (int64_t p : prime_divisors(T_abs(DeterminantMat(seed.ring.B)))) {
    l_ell_all.insert(p);
  }
  l_ell_all.insert(2);
  for (auto const &kv : target) {
    l_ell_all.insert(kv.first.ell);
  }
  // 2. a maximal lattice, at the primes of E above them
  for (int64_t ell : l_ell_all) {
    std::vector<HermIdealPrime<Tint>> const &l_prime = cache.at(ell);
    for (size_t j = 0; j < l_prime.size(); j++) {
      MyMatrix<T> P_inv = HermIdealInverse(ord, l_prime[j]);
      bool progress = true;
      while (progress) {
        progress = false;
        HermLatticeData<T, Tint> L = lattice_of(seed.basis);
        HermLineSpace ls = HermGetLineSpace(L, l_prime[j]);
        for (size_t idx = 0; idx < ls.n_line() && !progress; idx++) {
          std::vector<int64_t> z = ls.line_vector(idx);
          MyVector<Tint> x(n);
          for (int k = 0; k < n; k++) {
            x(k) = UniversalScalarConversion<Tint, int64_t>(z[k]);
          }
          MyMatrix<T> Bn = HermAddIdealTimes(L, P_inv, x);
          MyMatrix<T> A = Bn * seed.basis;
          // the norm integral lattices have an integral trace form
          MyMatrix<T> G_A = A * UniversalMatrixConversion<T, Tint>(seed.ring.B) * A.transpose();
          if (!IsIntegralMatrix(G_A)) {
            continue;
          }
          if (HermIsNormIntegral(lattice_of(A))) {
            seed.basis = A;
            progress = true;
          }
        }
      }
    }
  }
  // 3. the descent, prime of K_0 by prime of K_0
  HermLocalSymbol unimod;
  unimod.blocks.push_back({0, m, 1});
  for (int64_t ell : l_ell_all) {
    std::vector<HermIdealPrime<Tint>> const &l_prime = cache.at(ell);
    for (size_t j = 0; j < l_prime.size(); j++) {
      if (HermK0Representative(l_prime, j) != j) {
        continue;
      }
      HermK0Ref ref{ell, j};
      bool ramified = l_prime[j].type == HermIdealType::Ramified;
      bool dyadic = ramified && ell == 2;
      int diff_val = dyadic ? HermDifferentVal(ord, l_prime, j) : 0;
      HermLocalSymbol tgt = target.count(ref) ? target[ref] : unimod;
      if (dyadic && tgt.norms.empty()) {
        std::cerr << "HERMSEED: the symbol at the dyadic ramified prime "
                  << ell << " needs the norm valuations\n";
        throw TerminalException{1};
      }
      int vol_tgt = HermLocalSymbolVolume(tgt);
      std::vector<size_t> l_side{j};
      if (l_prime[j].type == HermIdealType::Split) {
        l_side.push_back(l_prime[j].conj_index);
      }
      HermLocalSymbol sym0 = HermLocalSymbolOfLattice(lattice_of(seed.basis), l_prime, j);
      if (HermLocalSymbolEquivalent(sym0, tgt, ramified, dyadic, diff_val)) {
        continue;
      }
      std::deque<MyMatrix<T>> queue{seed.basis};
      std::set<std::string> seen{HermLocalSymbolString(sym0)};
      std::optional<MyMatrix<T>> found;
      while (!queue.empty() && !found) {
        MyMatrix<T> A = queue.front();
        queue.pop_front();
        HermLatticeData<T, Tint> L = lattice_of(A);
        for (size_t side : l_side) {
          if (found) {
            break;
          }
          HermLineSpace ls = HermGetLineSpace(L, l_prime[side]);
          size_t n_line = ls.n_line();
          size_t n_try = std::min<size_t>(n_line, 3000);
          for (size_t t = 0; t < n_try && !found; t++) {
            size_t idx = (n_try == n_line)
                             ? t
                             : static_cast<size_t>(random_u64() % n_line);
            std::vector<int> c = ProjectiveLineOfIndex(
                ls.m, static_cast<int>(ls.F.q), idx);
            MyMatrix<T> A_sub = HermHyperplaneSublattice(L, ls, c) * A;
            HermLocalSymbol sym =
                HermLocalSymbolOfLattice(lattice_of(A_sub), l_prime, j);
            if (HermLocalSymbolVolume(sym) > vol_tgt) {
              continue;
            }
            if (HermLocalSymbolEquivalent(sym, tgt, ramified, dyadic,
                                          diff_val)) {
              found = A_sub;
              break;
            }
            std::string key = HermLocalSymbolString(sym);
            if (seen.insert(key).second) {
              queue.push_back(A_sub);
            }
          }
        }
      }
      if (!found) {
        std::cerr << "HERMSEED: no sublattice with the local symbol "
                  << HermLocalSymbolString(tgt) << " at " << ell
                  << ": the symbol is not that of an integral lattice of the "
                     "space\n";
        throw TerminalException{1};
      }
      seed.basis = *found;
    }
  }
  // 4. the check
  HermLatticeData<T, Tint> L = lattice_of(seed.basis);
  for (int64_t ell : l_ell_all) {
    std::vector<HermIdealPrime<Tint>> const &l_prime = cache.at(ell);
    for (size_t j = 0; j < l_prime.size(); j++) {
      if (HermK0Representative(l_prime, j) != j) {
        continue;
      }
      HermK0Ref ref{ell, j};
      HermLocalSymbol tgt = target.count(ref) ? target[ref] : unimod;
      HermLocalSymbol sym = HermLocalSymbolOfLattice(L, l_prime, j);
      bool ramified = l_prime[j].type == HermIdealType::Ramified;
      bool dyadic = ramified && ell == 2;
      int diff_val = dyadic ? HermDifferentVal(ord, l_prime, j) : 0;
      if (!HermLocalSymbolEquivalent(sym, tgt, ramified, dyadic, diff_val)) {
        std::cerr << "HERMSEED: the seed has the symbol "
                  << HermLocalSymbolString(sym) << " at " << ell
                  << " instead of " << HermLocalSymbolString(tgt) << "\n";
        throw TerminalException{1};
      }
    }
  }
  return seed;
}

// clang-format off
#endif  // SRC_GENERA_HERMITIANGENUSSEED_H_
// clang-format on
