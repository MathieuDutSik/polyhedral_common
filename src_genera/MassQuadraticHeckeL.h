// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSQUADRATICHECKEL_H_
#define SRC_GENERA_MASSQUADRATICHECKEL_H_

// clang-format off
#include "MassAbelianQuadraticLocal.h"
#include "RealFunctions.h"
#include <cmath>
#include <chrono>
#include <functional>
#include <map>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_QUADRATIC_HECKE_L
#endif

/*
  The quadratic Hecke character chi = chi_{E/K} of E = K(sqrt d) over the
  totally real abelian field K = Q(theta), and its L-function

      L(s, chi) = zeta_E(s) / zeta_K(s) = prod_p (1 - chi(p) N(p)^{-s})^{-1}.

  For each rational prime p the primes P of K over p are listed with their
  norm p^f and chi(P) = 1 (split in E), -1 (inert) or 0 (ramified). At a
  prime p not dividing 2 disc(f) N(d) they are read off the factors g of f
  modulo p: chi(P) is the quadratic character of d modulo (p, g(theta)).
  At the other primes the local arithmetic of LocalPrimeK is used: at an
  odd P, chi(P) = 0 when v_P(d) is odd and otherwise the character of the
  unit part; at a dyadic P of ramification index e, with d = pi^v u, the
  extension K_P(sqrt d) is ramified when v is odd (conductor exponent
  2e + 1) or when the quadratic defect delta of u is odd (exponent 2e + 1 -
  delta), unramified inert when delta = 2e and split when u is a square.
  The conductor f of chi is the product of the P^{c_P}, and N(f) = |D_E| /
  D_K^2.
 */
struct QuadHeckePrime {
  int f;
  int chi;
  // the exponent of P in the conductor
  int cond_exp;
};

template <typename T> struct QuadHeckeCharacter {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int n;
  // the parity: 1 when E is CM (chi odd at every real place), 0 when E is
  // totally real
  int delta;
  Tint D_K;
  Tint norm_conductor;
  // 2 and the primes of disc(f) and of N(d)
  std::set<int64_t> bad;
  std::map<int64_t, std::vector<QuadHeckePrime>> primes;
};

/*
  Arithmetic in F_p[x]/(f) for an odd p < 2^61 and f monic of degree n, on
  fixed arrays, in Montgomery form (x R mod p, R = 2^64): a product is one
  128-bit multiplication and a reduction without division, the sums stay
  below 2p < 2^62. The inner loop
  of the characters at the good primes.
 */
struct QuadHeckeModRing {
  int n;
  uint64_t p;
  // -p^{-1} mod 2^64, and R^2 mod p
  uint64_t pinv;
  uint64_t r2;
  // the coefficients f_0, ..., f_{n-1} of f, in Montgomery form
  std::vector<uint64_t> f;
  std::vector<uint64_t> tmp;
  void init(uint64_t p_in) {
    p = p_in;
    // Newton: x <- x (2 - p x) doubles the correct low bits of p^{-1}
    uint64_t x = p;
    for (int i = 0; i < 6; i++) {
      x *= 2 - p * x;
    }
    pinv = ~x + 1;
    unsigned __int128 R = (static_cast<unsigned __int128>(1) << 64) % p;
    r2 = static_cast<uint64_t>(R * R % p);
  }
  uint64_t redc(unsigned __int128 T) const {
    uint64_t m = static_cast<uint64_t>(T) * pinv;
    unsigned __int128 t = (T + static_cast<unsigned __int128>(m) * p) >> 64;
    uint64_t r = static_cast<uint64_t>(t);
    return (r >= p) ? r - p : r;
  }
  uint64_t to_mont(uint64_t a) const {
    return redc(static_cast<unsigned __int128>(a) * r2);
  }
  uint64_t from_mont(uint64_t a) const { return redc(a); }
  uint64_t mont_mul(uint64_t a, uint64_t b) const {
    return redc(static_cast<unsigned __int128>(a) * b);
  }
  // tmp (degree <= 2n - 2) reduced modulo f into out
  void reduce(std::vector<uint64_t> &out) {
    for (int k = 2 * n - 2; k >= n; k--) {
      uint64_t c = tmp[k];
      if (c == 0) {
        continue;
      }
      // x^k = x^{k-n} x^n = -x^{k-n} sum_i f_i x^i
      for (int i = 0; i < n; i++) {
        uint64_t prod = mont_mul(c, f[i]);
        uint64_t v = tmp[k - n + i];
        tmp[k - n + i] = (v >= prod) ? v - prod : v + p - prod;
      }
    }
    for (int i = 0; i < n; i++) {
      out[i] = tmp[i];
    }
  }
  uint64_t add(uint64_t a, uint64_t b) const {
    uint64_t v = a + b;
    return (v >= p) ? v - p : v;
  }
  // out = a b in F_p[x]/(f), all in Montgomery form (out may alias a or b)
  void mul(std::vector<uint64_t> const &a, std::vector<uint64_t> const &b,
           std::vector<uint64_t> &out) {
    std::fill(tmp.begin(), tmp.end(), 0);
    for (int i = 0; i < n; i++) {
      if (a[i] != 0) {
        for (int j = 0; j < n; j++) {
          tmp[i + j] = add(tmp[i + j], mont_mul(a[i], b[j]));
        }
      }
    }
    reduce(out);
  }
  // out = a^2, with the products a_i a_j, i < j, counted twice
  void sqr(std::vector<uint64_t> const &a, std::vector<uint64_t> &out) {
    std::fill(tmp.begin(), tmp.end(), 0);
    for (int i = 0; i < n; i++) {
      if (a[i] == 0) {
        continue;
      }
      tmp[2 * i] = add(tmp[2 * i], mont_mul(a[i], a[i]));
      uint64_t a2 = add(a[i], a[i]);
      for (int j = i + 1; j < n; j++) {
        tmp[i + j] = add(tmp[i + j], mont_mul(a2, a[j]));
      }
    }
    reduce(out);
  }
  // out = base^e, with a window of 4 bits
  void pow(std::vector<uint64_t> const &base, uint64_t e,
           std::vector<uint64_t> &out) {
    std::vector<std::vector<uint64_t>> tab(16, std::vector<uint64_t>(n, 0));
    tab[0][0] = to_mont(1);
    tab[1] = base;
    for (int i = 2; i < 16; i++) {
      mul(tab[i - 1], base, tab[i]);
    }
    int nbits = 0;
    for (uint64_t x = e; x > 0; x >>= 1) {
      nbits++;
    }
    int top = ((nbits + 3) / 4) * 4;
    std::vector<uint64_t> res = tab[0];
    for (int sh = top - 4; sh >= 0; sh -= 4) {
      for (int t = 0; t < 4; t++) {
        sqr(res, res);
      }
      int w = static_cast<int>((e >> sh) & 15);
      if (w != 0) {
        mul(res, tab[w], res);
      }
    }
    out = res;
  }
};

/*
  The characters at a good prime p (p prime to 2 disc(f) N(d), p > 2n,
  p < 2^61) without factoring f modulo p: K being Galois, the primes over p
  have the same residue degree f_p, the order of p in (Z/N)^* / H. In
  O_K / p = (+)_P F_P the element w = prod_{i < f_p} d^{p^i} has the norm
  N_{F_P/F_p}(d) as component at P, and u = w^{(p-1)/2} is chi(P) = +-1
  there, so that Tr(u) = f_p sum_P chi(P) modulo p, which determines the
  sum as |f_p sum_P chi(P)| <= n < p / 2.
 */
inline std::vector<QuadHeckePrime>
QuadHeckeGoodPrime(QuadHeckeModRing &R, std::vector<uint64_t> const &dq,
                   std::vector<uint64_t> const &traces, int f_p) {
  int n = R.n;
  uint64_t p = R.p;
  std::vector<uint64_t> dm(n);
  for (int i = 0; i < n; i++) {
    dm[i] = R.to_mont(dq[i]);
  }
  std::vector<uint64_t> w = dm, z = dm;
  for (int i = 1; i < f_p; i++) {
    R.pow(z, p, z);
    R.mul(w, z, w);
  }
  std::vector<uint64_t> u;
  R.pow(w, (p - 1) / 2, u);
  uint64_t tr = 0;
  for (int i = 0; i < n; i++) {
    uint64_t prod = static_cast<uint64_t>(
        static_cast<unsigned __int128>(R.from_mont(u[i])) * traces[i] % p);
    tr = (tr + prod) % p;
  }
  int64_t S = (tr > p / 2) ? static_cast<int64_t>(tr) - static_cast<int64_t>(p)
                           : static_cast<int64_t>(tr);
  int g = n / f_p;
#ifdef SANITY_CHECK_QUADRATIC_HECKE_L
  if (S % f_p != 0 || std::abs(S / f_p) > g) {
    std::cerr << "QHECKE: inconsistent trace " << S << " at p=" << p << "\n";
    throw TerminalException{1};
  }
#endif
  S /= f_p;
  std::vector<QuadHeckePrime> lp;
  for (int i = 0; i < (g + S) / 2; i++) {
    lp.push_back({f_p, 1, 0});
  }
  for (int i = 0; i < (g - S) / 2; i++) {
    lp.push_back({f_p, -1, 0});
  }
  return lp;
}

/*
  The characters at the good primes p > 2n of the field K = Q(theta) of
  conductor N with the subgroup H, for d on the powers of theta (made
  integral): the residue degree of p is its order in (Z/N)^* / H, and the
  sum of the chi(P) over the primes of degree 1, the coefficient c_1(p) of
  p^{-s} in L(s, chi), is 0 unless p splits completely.
 */
template <typename Tint> struct QuadHeckeGoodPrimes {
  int n;
  int64_t N;
  std::vector<char> inH;
  std::vector<Tint> phi;
  std::vector<Tint> dint;
  std::vector<Tint> trZ;
  QuadHeckeModRing R;
  std::vector<uint64_t> dq;
  std::vector<uint64_t> tr;
  int residue_degree(int64_t p) const {
    int f_p = 1;
    int64_t x = p % N;
    while (!inH[x]) {
      x = static_cast<int64_t>(static_cast<__int128>(x) * (p % N) % N);
      f_p++;
    }
    return f_p;
  }
  // phi, dint and trZ as machine integers, when they all fit
  bool small = false;
  std::vector<int64_t> phi_s, dint_s, trZ_s;
  void make_small() {
    Tint lim = UniversalScalarConversion<Tint, int64_t>(int64_t(1) << 62);
    small = true;
    for (auto const *v : {&phi, &dint, &trZ}) {
      for (auto const &x : *v) {
        if (T_abs(x) >= lim) {
          small = false;
        }
      }
    }
    if (small) {
      for (int i = 0; i < n; i++) {
        phi_s.push_back(UniversalScalarConversion<int64_t, Tint>(phi[i]));
        dint_s.push_back(UniversalScalarConversion<int64_t, Tint>(dint[i]));
        trZ_s.push_back(UniversalScalarConversion<int64_t, Tint>(trZ[i]));
      }
    }
  }
  void setup(int64_t p) {
    R.init(static_cast<uint64_t>(p));
    if (small) {
      auto residue = [&](int64_t c) -> uint64_t {
        int64_t r = c % p;
        return static_cast<uint64_t>((r < 0) ? r + p : r);
      };
      for (int i = 0; i < n; i++) {
        R.f[i] = R.to_mont(residue(phi_s[i]));
        dq[i] = residue(dint_s[i]);
        tr[i] = residue(trZ_s[i]);
      }
      return;
    }
    Tint pz = UniversalScalarConversion<Tint, int64_t>(p);
    auto residue = [&](Tint const &c) -> uint64_t {
      Tint r = ResInt(c, pz);
      return static_cast<uint64_t>(UniversalScalarConversion<int64_t, Tint>(r));
    };
    for (int i = 0; i < n; i++) {
      R.f[i] = R.to_mont(residue(phi[i]));
      dq[i] = residue(dint[i]);
      tr[i] = residue(trZ[i]);
    }
  }
  std::vector<QuadHeckePrime> local(int64_t p) {
    setup(p);
    return QuadHeckeGoodPrime(R, dq, tr, residue_degree(p));
  }
  int c1(int64_t p) {
    if (residue_degree(p) != 1) {
      return 0;
    }
    setup(p);
    int res = 0;
    for (auto const &P : QuadHeckeGoodPrime(R, dq, tr, 1)) {
      res += P.chi;
    }
    return res;
  }
};

template <typename T>
QuadHeckeGoodPrimes<typename underlying_z_ring<T>::ring_type>
MakeQuadHeckeGoodPrimes(AbelianField const &FK,
                        std::vector<typename underlying_z_ring<T>::ring_type> const &dint,
                        NumberField<T> const &K) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  QuadHeckeGoodPrimes<Tint> G;
  int n = K.deg;
  G.n = n;
  G.N = FK.n;
  G.inH.assign(FK.n, 0);
  for (int64_t h : FK.H) {
    G.inH[h] = 1;
  }
  for (int i = 0; i < n; i++) {
    G.phi.push_back(UniversalScalarConversion<Tint, int64_t>(K.Phi[i]));
  }
  G.dint = dint;
  std::vector<T> trT = PowerTraces<T>(K.Phi, n);
  for (int i = 0; i < n; i++) {
    G.trZ.push_back(GetNumerator_z(trT[i]));
  }
  G.R.n = n;
  G.R.f.resize(n);
  G.R.tmp.resize(2 * n);
  G.dq.resize(n);
  G.tr.resize(n);
  G.make_small();
  return G;
}

/*
  The characters at the primes p <= X and at the bad primes for the lattice (B,
  Omega) of rank m over K, with K and d from QuadDiscriminantK. At a bad prime p
  the discriminant is recomputed on a local generator at p (the same square
  class). Nothing when a bad prime has no local generator (a common index
  divisor). Each prime is passed to the callback, and stored in the
  result when store is set.
 */
template <typename T>
std::optional<QuadHeckeCharacter<T>> ComputeQuadHeckeCharacter(
    MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Omega,
    AbelianField const &FK, ElemE<T> const &d_in, NumberField<T> const &K,
    int m, int64_t X, bool store,
    std::function<void(int64_t, std::vector<QuadHeckePrime> const &)> const
        &f_prime) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  QuadHeckeCharacter<T> res;
  int n = K.deg;
  res.n = n;
  res.delta = (m % 4 == 0) ? 0 : 1;
  // d integral on the powers of theta: d den^2
  auto [dz, den] = ClearDenominators(d_in);
  std::vector<Tint> dint(n);
  for (int i = 0; i < n; i++) {
    dint[i] = dz[i] * den;
  }
  // the bad primes: 2, those of disc(f) and of N(d)
  std::set<int64_t> bad{2};
  {
    MyMatrix<T> Tm = PowerTraceMatrix<T>(K.Phi);
    ElemE<T> dT(n);
    for (int i = 0; i < n; i++) {
      dT[i] = UniversalScalarConversion<T, Tint>(dint[i]);
    }
    for (T const &x : {DeterminantMat(Tm), DeterminantMat(EMulMatrix(dT, K))}) {
      Tint num = T_abs(GetNumerator_z(x));
      if (num > 1) {
        for (Tint const &p : FactorsInt(num)) {
          bad.insert(UniversalScalarConversion<int64_t, Tint>(p));
        }
      }
    }
  }
  res.bad = bad;
  // D_K from the characters
  DirichletGroup G = ComputeDirichletGroup(FK.n);
  std::vector<std::vector<int64_t>> XK = FieldCharacters(G, FK.H, true);
  res.D_K = Tint(1);
  for (auto const &t : XK) {
    res.D_K *= UniversalScalarConversion<Tint, int64_t>(
        PrimitiveCharacter(G, t).first);
  }
  res.norm_conductor = Tint(1);
  // the primes up to X (a bit sieve), then the bad primes beyond X for the
  // conductor
  std::vector<bool> sieve(X + 1, true);
  std::vector<int64_t> plist_big;
  for (int64_t p : bad) {
    if (p > X) {
      plist_big.push_back(p);
    }
  }
  QuadHeckeGoodPrimes<Tint> GP = MakeQuadHeckeGoodPrimes(FK, dint, K);
  size_t i_big = 0;
  int64_t p_next = 2;
  auto next_prime = [&]() -> int64_t {
    while (p_next <= X) {
      int64_t p = p_next++;
      if (sieve[p]) {
        for (int64_t j = p * p; j <= X; j += p) {
          sieve[j] = false;
        }
        return p;
      }
    }
    if (i_big < plist_big.size()) {
      return plist_big[i_big++];
    }
    return 0;
  };
  for (int64_t p = next_prime(); p != 0; p = next_prime()) {
    std::vector<QuadHeckePrime> lp;
    if (bad.count(p) == 0 && p > 2 * n) {
      lp = GP.local(p);
    } else if (bad.count(p) == 0) {
      T_small_prime p_u = static_cast<T_small_prime>(p);
      PolyQ fq = PolyReduceSigned(K.Phi, p_u);
      std::vector<int64_t> dc;
      for (auto const &c : dint) {
        Tint r = ResInt(c, UniversalScalarConversion<Tint, int64_t>(p));
        dc.push_back(UniversalScalarConversion<int64_t, Tint>(r));
      }
      PolyQ dq = PolyReduceSigned(dc, p_u);
      for (auto const &[g, e] : FactorModQWithMultiplicity(fq)) {
        FiniteField F = MakeFiniteField(p_u, g);
        int chi = FqQuadraticCharacter(F, FqFromPoly(F, dq));
        lp.push_back({static_cast<int>(g.size()) - 1, chi, 0});
        (void)e;
      }
    } else {
      // the primes over p, each with a generator, its factor and the
      // discriminant on a K-basis of the lattice good at it
      struct BadPrime {
        LocalGenerator<T> lg;
        PolyQ g;
        int e;
        ElemE<T> dp;
      };
      std::vector<BadPrime> bps;
      auto disc_of = [&](LocalGenerator<T> const &lg,
                         std::vector<std::vector<ElemE<T>>> const &H)
          -> ElemE<T> {
        NumberField<T> Kp;
        Kp.d = 0;
        Kp.deg = n;
        Kp.Phi = lg.f;
        Kp.conj = lg.conj;
        ElemE<T> dp = DetK(H, Kp);
        if ((m * (m - 1) / 2) % 2 == 1) {
          dp = -dp;
        }
        return dp;
      };
      std::optional<LocalGenerator<T>> lg = FindLocalGenerator(B, Omega, p);
      if (lg) {
        ElemE<T> dp = disc_of(*lg, LocalOKGram(B, *lg, p));
        for (auto const &[g, e] :
             FactorModQWithMultiplicity(PolyReduceSigned(lg->f, p))) {
          bps.push_back({*lg, g, e, dp});
        }
      } else {
        // a common index divisor: a generator for each prime
        std::optional<std::vector<PrimeLocalGenerator<T>>> pgs =
            PrimeLocalGenerators(B, Omega, p,
                                 PrimeDecompositionCharacters(G, XK, p), false);
        if (!pgs) {
          return {};
        }
        for (auto const &pg : *pgs) {
          bps.push_back({pg.lg, pg.g, pg.e,
                         disc_of(pg.lg, LocalOKGramFromVectors(B, pg.lg,
                                                               pg.basis))});
        }
      }
      for (auto const &bp : bps) {
        PolyQ const &g = bp.g;
        int e = bp.e;
        ElemE<T> const &dp = bp.dp;
        LocalPrimeK<T> P = MakeLocalPrimeK(bp.lg, g, e, p);
        int fr = P.f_res;
        int v = LocalVal(dp, P);
        // the unit part, when v is even
        ElemE<T> u = (v % 2 == 0) ? LocalShift(dp, v, P) : dp;
        QuadHeckePrime hp{fr, 0, 0};
        if (p != 2) {
          if (v % 2 != 0) {
            hp.chi = 0;
            hp.cond_exp = 1;
          } else {
            hp.chi = FqQuadraticCharacter(P.F, LocalResidue(u, P));
          }
        } else {
          FiniteField const &F = P.F;
          if (v % 2 != 0) {
            hp.chi = 0;
            hp.cond_exp = 2 * e + 1;
          } else {
            int dl = QuadDefectUnit(u, P, F).first;
            if (dl == QD_INF) {
              hp.chi = 1;
            } else if (dl == 2 * e) {
              hp.chi = -1;
            } else {
              hp.chi = 0;
              hp.cond_exp = 2 * e + 1 - dl;
            }
          }
        }
        res.norm_conductor *= MyPow<Tint>(
            UniversalScalarConversion<Tint, int64_t>(p), hp.cond_exp * fr);
        lp.push_back(hp);
      }
    }
    f_prime(p, lp);
    if (store) {
      res.primes[p] = lp;
    }
  }
  return res;
}

/*
  The Dirichlet coefficients a_k, k <= X, of L(s, chi): a_k is the sum of
  chi(I) over the ideals I of K of norm k, multiplicative, with the local
  series prod_{P | p} 1 / (1 - chi(P) u^{f_P}) in u = p^{-s}.
 */
template <typename T>
std::vector<int64_t> QuadHeckeCoefficients(QuadHeckeCharacter<T> const &C,
                                           int64_t X) {
  std::vector<int64_t> a(X + 1, 0);
  a[1] = 1;
  // the smallest prime factor
  std::vector<int64_t> spf(X + 1, 0);
  for (int64_t i = 2; i <= X; i++) {
    if (spf[i] == 0) {
      for (int64_t j = i; j <= X; j += i) {
        if (spf[j] == 0) {
          spf[j] = i;
        }
      }
    }
  }
  std::map<int64_t, std::vector<int64_t>> local;
  for (auto const &[p, lp] : C.primes) {
    int kmax = 0;
    for (int64_t x = 1; x <= X / p; x *= p) {
      kmax++;
    }
    local[p] = QuadHeckeLocalSeries(lp, kmax);
  }
  for (int64_t k = 2; k <= X; k++) {
    int64_t p = spf[k];
    int64_t y = k;
    int e = 0;
    while (y % p == 0) {
      y /= p;
      e++;
    }
    auto it = local.find(p);
    int64_t cp = (it != local.end() && e < static_cast<int>(it->second.size()))
                     ? it->second[e]
                     : 0;
    a[k] = cp * a[y];
  }
  return a;
}

/*
  The value L(1 - r, chi) by the functional equation (Dokchitser, Computing
  special values of motivic L-functions, 2004). The completed function

      Lambda(s) = A^s gamma(s) L(s, chi),  gamma(s) = Gamma_R(s + delta)^n,
      Gamma_R(s) = pi^{-s/2} Gamma(s/2),  A = sqrt(|D_K| N(f)),

  with delta = 1 when E is CM and 0 when E is totally real, is entire and
  Lambda(s) = Lambda(1 - s), the root number of a quadratic character being
  1. With phi the inverse Mellin transform of gamma and G_s(t) = t^{-s}
  int_t^infty phi(x) x^s dx / x,

      Lambda(s) = sum_{k >= 1} a_k (G_s(k / A) + G_{1-s}(k / A)).

  Moving the line of integration of phi to the left,

      G_s(t) = gamma(s) t^{-s} - sum_{j >= 0} Res_{w = z_j} gamma(w) t^{-w} / (s
  - w)

  over the poles z_j = -delta - 2j of gamma, of order n. At w = z_j + eps,
  Gamma(-j + u) = ((-1)^j / j!) u^{-1} exp(S_j(u)) with S_j(u) = -gamma_E u
  + sum_{k >= 2} (-1)^k zeta(k) u^k / k + sum_{k >= 1} H_j^{(k)} u^k / k,
  H_j^{(k)} = sum_{i <= j} i^{-k}, so that the residue is

      t^{delta + 2j} C_j [eps^{n-1}] exp(n S_j(eps / 2) - eps (n/2) log pi)
                                     exp(-eps log t) sum_k eps^k / (s -
  z_j)^{k+1}

  with C_j = pi^{nj} 2^n ((-1)^j / j!)^n: a polynomial of degree n - 1 in
  log t. The series converges for every t; its terms reach about
  exp(n pi t^{2/n}) while G_s(t) is about exp(-n pi t^{2/n}), hence a
  working precision of about twice the target. The sum over k stops at
  t = k / A where n pi t^{2/n} exceeds the target.
 */
template <typename Tf> Tf QuadHeckeGamma(Tf const &s, int n, int delta) {
  Tf u = (s + real_from_int<Tf>(delta)) / real_from_int<Tf>(2);
  Tf g = real_exp(-u * real_log(real_pi<Tf>())) * real_gamma(u);
  return real_pow_int(g, n);
}

/*
  The kernel in a type Tf (double, DoubleDouble, mpreal or ArbReal, see
  RealFunctions.h), with the coefficients q[k][j] of t^{delta + 2j} (log
  t)^k of the series cut after J terms, and the bounds of what is cut:

  - By Cauchy on the circle |w - z_j| = 1 the j-th residue of gamma(w)
    t^{-w} / (s - w) is at most
      pi^{n(j+1/2)} (2 pi / Gamma(j + 1/2))^n |t^{-w}| / (|s - z_j| - 1)
    (res_const times |t^{-w}|), as (w + delta)/2 = -j + v with |v| = 1/2 and
    |Gamma(-j + v)| = |Gamma(1 + v)| / (|v| prod_{i <= j} |v - i|) <=
    sqrt(pi) / ((1/2) Gamma(j + 1/2) / sqrt(pi)). With t = e^u, u = x + iy,
    |t^{-w}| <= e^{(delta + 2j) x + |x| + |y|}. Consecutive bounds have the
    ratio pi^n e^{2x} / (j + 1/2)^n, so the tail is summed until the ratio
    is 1/2 and closed by a geometric series.
  - The majorant of |F(e^u)| on a region of the complex u, from the series
    with the absolute values of the coefficients and the tail above.
  For a ball type the tail is added to the radius of every value of F; for
  the other types it is not used in the values.
 */
template <typename Tf> struct QuadHeckeKernel {
  int n;
  int delta;
  int r;
  Tf gs;
  Tf g1s;
  // q[k][j]: the coefficient of t^{delta + 2j} (log t)^k
  std::vector<std::vector<Tf>> q;
  int J() const { return static_cast<int>(q[0].size()) - 1; }
  // the bound of the j-th residue without the factor |t^{-w}|, j > J
  Tf res_const(int j) const {
    Tf pi = real_pi<Tf>();
    Tf jh = real_from_int<Tf>(2 * j + 1) / real_from_int<Tf>(2);
    Tf base = real_from_int<Tf>(2) * pi / real_gamma(jh);
    Tf dist = real_from_int<Tf>(1 - r + delta + 2 * j) - real_from_int<Tf>(1);
    return real_exp(real_from_int<Tf>(n) * jh * real_log(pi)) *
           real_pow_int(base, n) * real_from_int<Tf>(2) / dist;
  }
  // an upper bound of the residues j > J for Re u in [., xe], |Re u| <= xa,
  // |Im u| <= y (+infinity when the bounds are not geometric soon enough)
  Tf tail(double xe, double xa, double y) const {
    Tf pi = real_pi<Tf>();
    Tf total = real_from_int<Tf>(0);
    Tf half = real_from_int<Tf>(1) / real_from_int<Tf>(2);
    Tf ex2 = real_exp(real_from_int<Tf>(2) * real_from_double<Tf>(xe));
    for (int j = J() + 1; j <= J() + 2000; j++) {
      Tf term = res_const(j) *
                real_exp(real_from_int<Tf>(delta + 2 * j) *
                             real_from_double<Tf>(xe) +
                         real_from_double<Tf>(xa) + real_from_double<Tf>(y));
      // the ratio to the next bound
      Tf ratio = real_pow_int(pi, n) * ex2 /
                 real_pow_int(real_from_int<Tf>(j) + half, n);
      if (real_certainly_le(ratio, half)) {
        total += real_from_int<Tf>(2) * term;
        return real_upper(total);
      }
      total += term;
    }
    // too far to the right (a large ellipse): no useful bound
    return real_infinity<Tf>();
  }
};

template <typename Tf>
QuadHeckeKernel<Tf> MakeQuadHeckeKernel(int n, int delta, int r, int J) {
  QuadHeckeKernel<Tf> ker;
  ker.n = n;
  ker.delta = delta;
  ker.r = r;
  Tf s = real_from_int<Tf>(r);
  Tf s1 = real_from_int<Tf>(1 - r);
  ker.gs = QuadHeckeGamma(s, n, delta);
  ker.g1s = QuadHeckeGamma(s1, n, delta);
  Tf pi = real_pi<Tf>();
  Tf logpi = real_log(pi);
  Tf euler = real_euler<Tf>();
  std::vector<Tf> zeta(n + 1, real_from_int<Tf>(0));
  for (int k = 2; k <= n; k++) {
    zeta[k] = real_zeta<Tf>(k);
  }
  std::vector<Tf> H(n + 1, real_from_int<Tf>(0));
  Tf C = real_pow_int(real_from_int<Tf>(2), n);
  Tf pin = real_pow_int(pi, n);
  std::vector<Tf> fact(n + 1, real_from_int<Tf>(1));
  for (int k = 1; k <= n; k++) {
    fact[k] = fact[k - 1] * real_from_int<Tf>(k);
  }
  ker.q.assign(n, std::vector<Tf>(J + 1, real_from_int<Tf>(0)));
  Tf half = real_from_int<Tf>(1) / real_from_int<Tf>(2);
  for (int j = 0; j <= J; j++) {
    if (j > 0) {
      Tf ij = real_from_int<Tf>(1) / real_from_int<Tf>(j);
      Tf pw = ij;
      for (int k = 1; k <= n; k++) {
        H[k] += pw;
        pw *= ij;
      }
      // C_j = C_{j-1} pi^n (-1/j)^n
      C *= pin * real_pow_int(-ij, n);
    }
    // A(eps) = n S_j(eps/2) - eps (n/2) log pi, up to eps^{n-1}
    std::vector<Tf> a(n, real_from_int<Tf>(0));
    Tf hk = half;
    for (int k = 1; k < n; k++) {
      Tf sk = (k == 1) ? Tf(H[1] - euler)
                       : Tf((((k % 2 == 0) ? zeta[k] : Tf(-zeta[k])) + H[k]) /
                            real_from_int<Tf>(k));
      a[k] = real_from_int<Tf>(n) * sk * hk;
      hk *= half;
    }
    if (n > 1) {
      a[1] -= logpi * real_from_int<Tf>(n) / real_from_int<Tf>(2);
    }
    // exp(A)
    std::vector<Tf> e(n, real_from_int<Tf>(0));
    e[0] = real_from_int<Tf>(1);
    for (int k = 1; k < n; k++) {
      Tf sum = real_from_int<Tf>(0);
      for (int i = 1; i <= k; i++) {
        sum += real_from_int<Tf>(i) * a[i] * e[k - i];
      }
      e[k] = sum / real_from_int<Tf>(k);
    }
    // sum_k eps^k (1/(s - z)^{k+1} + 1/(1 - s - z)^{k+1})
    Tf z = real_from_int<Tf>(-delta - 2 * j);
    Tf i1 = real_from_int<Tf>(1) / (s - z);
    Tf i2 = real_from_int<Tf>(1) / (s1 - z);
    std::vector<Tf> h(n, real_from_int<Tf>(0));
    Tf p1 = i1, p2 = i2;
    for (int k = 0; k < n; k++) {
      h[k] = p1 + p2;
      p1 *= i1;
      p2 *= i2;
    }
    for (int k = 0; k < n; k++) {
      // g_{n-1-k}
      Tf g = real_from_int<Tf>(0);
      for (int i = 0; i <= n - 1 - k; i++) {
        g += e[i] * h[n - 1 - k - i];
      }
      Tf sign = real_from_int<Tf>((k % 2 == 0) ? 1 : -1);
      ker.q[k][j] = C * g * sign / fact[k];
    }
  }
  return ker;
}

// F(e^u) = G_s(e^u) + G_{1-s}(e^u) for a real u (ud: u as a double), the
// tail of the series as an error for a ball type
template <typename Tf>
Tf QuadHeckeKernelValue(QuadHeckeKernel<Tf> const &ker, Tf const &u,
                        double ud) {
  Tf t2 = real_exp(real_from_int<Tf>(2) * u);
  Tf res = ker.gs * real_exp(-real_from_int<Tf>(ker.r) * u) +
           ker.g1s * real_exp(real_from_int<Tf>(ker.r - 1) * u);
  Tf corr = real_from_int<Tf>(0);
  Tf Lk = real_from_int<Tf>(1);
  int J = ker.J();
  for (int k = 0; k < ker.n; k++) {
    Tf poly = ker.q[k][J];
    for (int j = J - 1; j >= 0; j--) {
      poly = poly * t2 + ker.q[k][j];
    }
    corr += Lk * poly;
    Lk *= u;
  }
  if (ker.delta == 1) {
    corr *= real_exp(u);
  }
  res -= corr;
  if constexpr (RealIsBall<Tf>) {
    real_add_error(res, ker.tail(ud + 1e-9, std::abs(ud) + 1e-9, 0.0));
  }
  return res;
}

// An upper bound of |F(e^u)| for Re u in [xmin, xmax], |u| <= U, |Im u| <= y
template <typename Tf>
Tf QuadHeckeMajorant(QuadHeckeKernel<Tf> const &ker, double xmin, double xmax,
                     double U, double y) {
  Tf xm = real_from_double<Tf>(xmin), xM = real_from_double<Tf>(xmax);
  Tf Ub = real_from_double<Tf>(U);
  Tf M = real_abs(ker.gs) * real_exp(-real_from_int<Tf>(ker.r) * xm) +
         real_abs(ker.g1s) * real_exp(real_from_int<Tf>(ker.r - 1) * xM);
  Tf t2 = real_exp(real_from_int<Tf>(2) * xM);
  Tf Uk = real_from_int<Tf>(1);
  int J = ker.J();
  for (int k = 0; k < ker.n; k++) {
    Tf poly = real_upper(ker.q[k][J]);
    for (int j = J - 1; j >= 0; j--) {
      poly = poly * t2 + real_upper(ker.q[k][j]);
    }
    M += Uk * poly * real_exp(real_from_int<Tf>(ker.delta) * xM);
    Uk *= Ub;
  }
  M += ker.tail(xmax, std::max(std::abs(xmin), std::abs(xmax)), y);
  return real_upper(M);
}

/*
  The kernel on the terms k = 1, ..., X is evaluated through a Chebyshev
  interpolation p of degree D of F(e^u) on segments of the variable u =
  log t, the series being computed in Tf at the nodes only. A segment is
  halved until the error |F - p| is at most the tolerance:
  - for a ball type, the error is bounded: for F analytic in the
    Bernstein ellipse of parameter rho with |F| <= M there, |F - p_D| <= 4
    M rho^{-D} / (rho - 1) (Trefethen, Approximation theory and
    approximation practice, Thm 8.2), M from QuadHeckeMajorant, the best
    rho of a list being taken;
  - otherwise it is estimated by the last two coefficients |c_D| +
    |c_{D-1}|.
  Each segment keeps log2 max |F| at the nodes, from which the terms are
  evaluated in the cheapest type that has the bits they need.
 */
template <typename Tf> struct QuadHeckeInterp {
  std::vector<double> lo;
  std::vector<double> hi;
  std::vector<double> mid;
  std::vector<double> half;
  std::vector<std::vector<Tf>> coef;
  // the bound (a ball type) or the estimate of |F - p| on the segment
  std::vector<Tf> err;
  std::vector<double> log2_fmax;
};

template <typename Tf>
std::optional<QuadHeckeInterp<Tf>>
MakeQuadHeckeInterp(QuadHeckeKernel<Tf> const &ker, double u0, double u1,
                    double step, Tf const &tol, int D) {
  QuadHeckeInterp<Tf> res;
  Tf pi = real_pi<Tf>();
  // cos(pi k (2i + 1) / (2 (D + 1)))
  std::vector<std::vector<Tf>> cs(D + 1, std::vector<Tf>(D + 1));
  for (int k = 0; k <= D; k++) {
    for (int i = 0; i <= D; i++) {
      cs[k][i] = real_cos(pi * real_from_int<Tf>(k * (2 * i + 1)) /
                          real_from_int<Tf>(2 * (D + 1)));
    }
  }
  std::vector<double> rhos{1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 128};
  bool ok = true;
  std::function<void(double, double, int)> build = [&](double a, double b,
                                                       int depth) -> void {
    if (!ok) {
      return;
    }
    auto subdivide = [&]() -> void {
      if (depth > 40) {
        std::cerr << "QHECKE: no interpolation on [" << a << ", " << b
                  << "]\n";
        ok = false;
        return;
      }
      double mm = (a + b) / 2;
      build(a, mm, depth + 1);
      build(mm, b, depth + 1);
    };
    // slightly enlarged, for the rounding of u_k when assigning k
    double eps = 1e-9 * (1 + std::abs(a) + std::abs(b));
    double ae = a - eps, be = b + eps;
    double m = (ae + be) / 2, h = (be - ae) / 2;
    Tf err;
    if constexpr (RealIsBall<Tf>) {
      // the best bound over rho, before computing anything at the nodes
      bool have = false;
      for (double rho : rhos) {
        double ca = (rho + 1 / rho) / 2, cb = (rho - 1 / rho) / 2;
        Tf M = QuadHeckeMajorant(ker, m - h * ca, m + h * ca,
                                 std::abs(m) + h * ca, h * cb);
        Tf rb = real_from_double<Tf>(rho);
        Tf den = real_pow_int(rb, D) * (rb - real_from_int<Tf>(1));
        Tf e = real_upper(Tf(real_from_int<Tf>(4) * M / den));
        if (!have || real_certainly_le(e, err)) {
          err = e;
          have = true;
        }
      }
      if (!real_certainly_le(err, tol)) {
        subdivide();
        return;
      }
    }
    Tf mb = real_from_double<Tf>(m), hb = real_from_double<Tf>(h);
    std::vector<Tf> g(D + 1);
    double log2_fmax = -1e300;
    for (int i = 0; i <= D; i++) {
      // the node m + h cos(pi (2i+1) / (2(D+1)))
      Tf u = mb + hb * cs[1][i];
      g[i] = QuadHeckeKernelValue(ker, u, real_to_double(u));
      double gd = real_to_double(g[i]);
      if constexpr (!RealIsBall<Tf>) {
        if (!(gd > 0)) {
          std::cerr << "QHECKE: the kernel is not positive at u="
                    << real_to_double(u) << ", F=" << gd << "\n";
          throw TerminalException{1};
        }
      }
      log2_fmax = std::max(log2_fmax, std::log2(std::abs(gd) + 1e-300));
    }
    std::vector<Tf> c(D + 1);
    for (int k = 0; k <= D; k++) {
      Tf sum = real_from_int<Tf>(0);
      for (int i = 0; i <= D; i++) {
        sum += g[i] * cs[k][i];
      }
      c[k] = sum * real_from_int<Tf>(2) / real_from_int<Tf>(D + 1);
    }
    c[0] /= real_from_int<Tf>(2);
    if constexpr (!RealIsBall<Tf>) {
      err = real_abs(c[D]) + real_abs(c[D - 1]);
      if (!real_certainly_le(err, tol)) {
        subdivide();
        return;
      }
    }
    res.lo.push_back(a);
    res.hi.push_back(b);
    res.mid.push_back(m);
    res.half.push_back(h);
    res.coef.push_back(c);
    res.err.push_back(err);
    res.log2_fmax.push_back(log2_fmax);
  };
  for (double a = u0; a < u1; a += step) {
    build(a, std::min(a + step, u1), 0);
  }
  if (!ok) {
    return {};
  }
  return res;
}

// the segment of u: the last one starting at or below it
template <typename Tf>
size_t QuadHeckeSegment(QuadHeckeInterp<Tf> const &I, double ud) {
  size_t i = static_cast<size_t>(
      std::upper_bound(I.lo.begin(), I.lo.end(), ud) - I.lo.begin());
  return (i == 0) ? 0 : i - 1;
}

// p(u) on a segment of center m and half-width h, by Clenshaw, in a type
// Te holding the coefficients c
template <typename Te>
Te QuadHeckeChebyshevValue(std::vector<Te> const &c, double m, double h,
                           Te const &u) {
  Te x = (u - real_from_double<Te>(m)) / real_from_double<Te>(h);
  Te x2 = real_from_int<Te>(2) * x;
  Te b1 = real_from_int<Te>(0), b2 = real_from_int<Te>(0);
  for (int k = static_cast<int>(c.size()) - 1; k >= 1; k--) {
    Te b0 = x2 * b1 - b2 + c[k];
    b2 = b1;
    b1 = b0;
  }
  return x * b1 - b2 + c[0];
}

// The coefficients of the interpolation in the type Te, at its working
// precision
template <typename Te, typename Tf>
std::vector<std::vector<Te>>
QuadHeckeInterpCoefficients(QuadHeckeInterp<Tf> const &I) {
  std::vector<std::vector<Te>> res;
  for (auto const &c : I.coef) {
    std::vector<Te> ce;
    for (auto const &x : c) {
      ce.push_back(real_convert<Te>(x));
    }
    res.push_back(ce);
  }
  return res;
}

// The local series prod_{P | p} 1 / (1 - chi(P) u^{f_P}) up to u^kmax
inline std::vector<int64_t>
QuadHeckeLocalSeries(std::vector<QuadHeckePrime> const &lp, int kmax) {
  std::vector<int64_t> c(kmax + 1, 0);
  c[0] = 1;
  for (auto const &P : lp) {
    std::vector<int64_t> nc(kmax + 1, 0);
    for (int k = 0; k <= kmax; k++) {
      int64_t s = 0;
      int64_t pw = 1;
      for (int i = 0; i * P.f <= k; i++) {
        s += pw * c[k - i * P.f];
        pw *= P.chi;
        if (pw == 0) {
          break;
        }
      }
      nc[k] = s;
    }
    c = nc;
  }
  return c;
}

/*
  The data of L(s, chi) for the lattice: the character, and for the
  Dirichlet coefficients a_k, k <= X, the local series at the primes p <= S,
  S = max(floor(sqrt X), 2n), the a_j for j <= S, and the coefficients
  c_1(q) of q^{-s} at the bad primes q > S; at the other primes q > S,
  c_1(q) is computed when needed (QuadHeckeGoodPrimes). Everything is of
  size O(sqrt X).
 */
template <typename T> struct QuadHeckeLData {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> B;
  std::vector<MyMatrix<T>> Omega;
  AbelianField FK;
  ElemE<T> d;
  NumberField<T> K;
  int m;
  int n;
  int delta;
  Tint A2;
  std::set<int64_t> bad;
  QuadHeckeGoodPrimes<Tint> GP;
  int64_t X;
  int64_t S;
  std::vector<std::pair<int64_t, std::vector<int64_t>>> small;
  std::vector<int64_t> a_small;
  std::map<int64_t, int> c1_bad;
  // log p for the small primes, log j for j <= S
  std::vector<DoubleDouble> log_small;
  std::vector<DoubleDouble> log_j;
};

template <typename T>
std::optional<QuadHeckeLData<T>>
MakeQuadHeckeLData(MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Omega,
                   AbelianField const &FK, ElemE<T> const &d,
                   NumberField<T> const &K, int m) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  auto f_none = [](int64_t, std::vector<QuadHeckePrime> const &) -> void {};
  std::optional<QuadHeckeCharacter<T>> C =
      ComputeQuadHeckeCharacter(B, Omega, FK, d, K, m, 2, false, f_none);
  if (!C) {
    return {};
  }
  auto [dz, den] = ClearDenominators(d);
  std::vector<Tint> dint;
  for (auto const &c : dz) {
    dint.push_back(c * den);
  }
  QuadHeckeLData<T> res{B,
                        Omega,
                        FK,
                        d,
                        K,
                        m,
                        C->n,
                        C->delta,
                        C->D_K * C->norm_conductor,
                        C->bad,
                        MakeQuadHeckeGoodPrimes(FK, dint, K),
                        0,
                        0,
                        {},
                        {},
                        {},
                        {},
                        {}};
  return res;
}

template <typename T> void QuadHeckeLExtend(QuadHeckeLData<T> &Ld, int64_t X) {
  if (X <= Ld.X) {
    return;
  }
  int64_t S = 1;
  while ((S + 1) * (S + 1) <= X) {
    S++;
  }
  S = std::max<int64_t>(S, 2 * Ld.n);
  Ld.small.clear();
  Ld.c1_bad.clear();
  auto f_prime = [&](int64_t p, std::vector<QuadHeckePrime> const &lp) -> void {
    if (p <= S) {
      int kmax = 0;
      for (int64_t x = 1; x <= X / p; x *= p) {
        kmax++;
      }
      Ld.small.push_back({p, QuadHeckeLocalSeries(lp, kmax)});
    } else if (p <= X) {
      int c = 0;
      for (auto const &P : lp) {
        if (P.f == 1) {
          c += P.chi;
        }
      }
      Ld.c1_bad[p] = c;
    }
  };
  ComputeQuadHeckeCharacter(Ld.B, Ld.Omega, Ld.FK, Ld.d, Ld.K, Ld.m, S, false,
                            f_prime);
  // a_j for j <= S
  Ld.a_small.assign(S + 1, 1);
  Ld.a_small[0] = 0;
  std::vector<int64_t> rem(S + 1);
  for (int64_t j = 0; j <= S; j++) {
    rem[j] = j;
  }
  for (auto const &[p, c] : Ld.small) {
    for (int64_t j = p; j <= S; j += p) {
      int e = 0;
      while (rem[j] % p == 0) {
        rem[j] /= p;
        e++;
      }
      Ld.a_small[j] *= (e < static_cast<int>(c.size())) ? c[e] : 0;
    }
  }
  {
    mpfr_prec_t old = mpfr::mpreal::get_default_prec();
    mpfr::mpreal::set_default_prec(160);
    Ld.log_small.clear();
    for (auto const &pc : Ld.small) {
      Ld.log_small.push_back(DoubleDoubleFromMpreal(
          mpfr::log(mpfr::mpreal(static_cast<unsigned long>(pc.first)))));
    }
    Ld.log_j.assign(S + 1, DoubleDouble(0.0));
    for (int64_t j = 2; j <= S; j++) {
      Ld.log_j[j] = DoubleDoubleFromMpreal(
          mpfr::log(mpfr::mpreal(static_cast<unsigned long>(j))));
    }
    mpfr::mpreal::set_default_prec(old);
  }
  Ld.X = X;
  Ld.S = S;
}

/*
  f(k, a_k, log k) for the k <= X with a_k != 0, X <= Ld.X, in no particular
  order, log k in DoubleDouble: the S-smooth k by a depth-first walk over the
  local series of the primes p <= S (skipping the zero coefficients and what
  lies under them), and the k = j q with a prime q > S (at most one, as S^2
  >= X) and j <= X / q <= S from a segmented sieve over (S, X], the primes
  with c_1(q) = 0 being skipped.
 */
template <typename T, typename F>
void QuadHeckeForEachTerm(QuadHeckeLData<T> &Ld, int64_t X, F &&f) {
  int64_t S = Ld.S;
  std::function<void(size_t, int64_t, int64_t, DoubleDouble)> dfs =
      [&](size_t i0, int64_t k, int64_t ak, DoubleDouble lk) -> void {
    for (size_t i = i0; i < Ld.small.size(); i++) {
      int64_t p = Ld.small[i].first;
      if (k > X / p) {
        break;
      }
      std::vector<int64_t> const &c = Ld.small[i].second;
      int64_t kk = k;
      DoubleDouble lkk = lk;
      for (size_t e = 1; kk <= X / p; e++) {
        kk *= p;
        lkk += Ld.log_small[i];
        if (e < c.size() && c[e] != 0) {
          f(kk, ak * c[e], lkk);
          dfs(i + 1, kk, ak * c[e], lkk);
        }
      }
    }
  };
  f(int64_t(1), int64_t(1), DoubleDouble(0.0));
  dfs(0, 1, 1, DoubleDouble(0.0));
  // the primes q in (S, X]
  const int64_t block = int64_t(1) << 22;
  std::vector<char> is_p;
  for (int64_t lo = S + 1; lo <= X; lo += block) {
    int64_t hi = std::min(lo + block, X + 1);
    is_p.assign(hi - lo, 1);
    for (auto const &pc : Ld.small) {
      int64_t p = pc.first;
      if (p * p >= hi) {
        break;
      }
      int64_t first = std::max(p * p, ((lo + p - 1) / p) * p);
      for (int64_t x = first; x < hi; x += p) {
        is_p[x - lo] = 0;
      }
    }
    for (int64_t q = lo; q < hi; q++) {
      if (!is_p[q - lo]) {
        continue;
      }
      int c = 0;
      if (Ld.bad.count(q) == 1) {
        auto it = Ld.c1_bad.find(q);
        c = (it == Ld.c1_bad.end()) ? 0 : it->second;
      } else {
        c = Ld.GP.c1(q);
      }
      if (c == 0) {
        continue;
      }
      int64_t jmax = X / q;
      DoubleDouble lq = DoubleDoubleLog(DoubleDouble(q));
      for (int64_t j = 1; j <= jmax; j++) {
        int64_t aj = Ld.a_small[j];
        if (aj != 0) {
          f(j * q, aj * c, lq + Ld.log_j[j]);
        }
      }
    }
  }
}

/*
  L(1 - r, chi) in the type Tf: Lambda(1 - r) from the sum, divided by
  A^{1-r} gamma(1-r) (no pole there: (1 - r + delta)/2 is a half-integer
  as r = delta + 1 mod 2). Nothing when the number of terms exceeds the
  limit, or when the interpolation cannot be built.

  - For a ball type (ArbReal) the result is certified: the number X of
    terms comes from the bound of the terms k > X (as phi > 0, G_s(t) =
    t^{-s} int_t^infty phi(x) x^s dx / x <= gamma(sigma) t^{-sigma} for
    every sigma >= s, and |a_k| is at most the k-th coefficient of
    zeta(s)^n, so that for sigma > 2
      sum_{k > X} |a_k| F(k/A) <= 2 gamma(sigma) A^sigma X^{2-sigma}
      zeta(2)^n,
    sigma being optimized), and that bound, the interpolation errors and
    the rounding of every term are in the radius. The ball has an absolute
    radius of about 2^{-prec_bits} relative to the size of Lambda.
  - Otherwise X is the point where the kernel n pi t^{2/n} falls below the
    target, the terms where F < 2^{-abs_bits} being dropped, for about
    prec_bits relative bits. Each term is evaluated in double when it needs
    at most 45 bits, in DoubleDouble up to 85 and in Tf beyond; with
    all_double every term is in double precision, whatever it needs.

  The kernel and the interpolation are built in Tf at the working
  precision wprec: the series reaches about exp(n pi t^{2/n}) at t =
  t_max while F is about its inverse, hence about twice the target bits.
  For double or DoubleDouble this is only meaningful for the smallest
  conductors.
 */
template <typename Tf, typename T>
std::optional<Tf> QuadHeckeLNegative(QuadHeckeLData<T> &Ld, int prec_bits,
                                     bool all_double) {
  constexpr bool is_ball = RealIsBall<Tf>;
  int r = Ld.m / 2;
  int n = Ld.n;
  int delta = Ld.delta;
  double logA = 0.5 * std::log(real_from_integer<double>(Ld.A2));
  // |Lambda(1 - r)| = A^r gamma(r) |L(r)|, with L(r) of size about 1
  double log2_scale = r * logA / std::log(2.0) +
                      std::log2(std::abs(QuadHeckeGamma(
                          static_cast<double>(r), n, delta)));
  double target = std::max(20.0, prec_bits - log2_scale + (is_ball ? 20 : 10));
  double X_double = 0;
  double X_limit = 0;
  if constexpr (is_ball) {
    // the bound of the terms k > X, optimized over sigma
    double lz2 = n * std::log(M_PI * M_PI / 6);
    auto log_tail = [&](double logX) -> double {
      double best = 1e300;
      for (double sigma = std::max(3.0, r + 0.0); sigma < 2000;
           sigma *= 1.02) {
        double lg = n * (-(sigma + delta) / 2 * std::log(M_PI) +
                         std::lgamma((sigma + delta) / 2));
        double v =
            std::log(2.0) + lg + sigma * logA + (2 - sigma) * logX + lz2;
        best = std::min(best, v);
      }
      return best;
    };
    double goal = -(target + 8) * std::log(2.0);
    double logX = logA;
    while (log_tail(logX) > goal) {
      logX += 0.05;
    }
    X_double = std::ceil(std::exp(logX));
    X_limit = 4.0e8;
  } else {
    // the cutoff: n pi t^{2/n} - (r + n) log t > target log 2 + log(A t)
    double t_cut = 1.0;
    auto tail_ok = [&](double t) -> bool {
      return n * M_PI * std::pow(t, 2.0 / n) - (r + n) * std::log(t) >=
             (target + 40) * std::log(2.0) + logA + std::log(t) + 10;
    };
    while (!tail_ok(t_cut)) {
      t_cut *= 1.05;
    }
    X_double = std::ceil(t_cut * std::exp(logA)) + 1;
    X_limit = 1.0e12;
  }
  // each term to 2^{-abs_bits}, the X terms then add up to 2^{-target}
  double abs_bits = target + std::log2(X_double) + (is_ball ? 24 : 6);
  double u0 = -logA - 0.01;
  double u1 = std::log(X_double) - logA + 0.01;
  // J: the residue bounds at the largest Re u of the ellipses for a ball
  // type, at u1 otherwise; the series reaches its peak there
  double xJ = is_ball ? u1 + 0.5 : u1;
  int J = 1;
  while (true) {
    double jh = J + 0.5;
    double lt = n * jh * std::log(M_PI) + n * std::log(2 * M_PI) -
                n * std::lgamma(jh) + (delta + 2 * J) * xJ + std::abs(xJ);
    double ratio = n * std::log(M_PI) + 2 * xJ - n * std::log(jh);
    if (J > 8 && lt / std::log(2.0) < -abs_bits - 20 && ratio < -1) {
      break;
    }
    J++;
  }
  // the peak of the series, at t = e^{xJ}
  double log2_peak = n * M_PI * std::exp(2.0 * xJ / n) / std::log(2.0);
  double log2_Fmax = std::max(0.0, log2_scale) + 8;
  long wprec =
      static_cast<long>(abs_bits + std::max(log2_peak, log2_Fmax) + 64);
  long eprec = static_cast<long>(abs_bits + log2_Fmax + 32);
  int D = is_ball ? std::clamp(
                        static_cast<int>((log2_peak + abs_bits) / 3) + 10, 30,
                        240)
                  : 30;
  double step = is_ball ? 0.25 : 1.0;
  RealPrecisionGuard<Tf> g_w(wprec);
  QuadHeckeKernel<Tf> ker = MakeQuadHeckeKernel<Tf>(n, delta, r, J);
  Tf tol = real_pow_int(real_from_int<Tf>(2), -static_cast<long>(abs_bits));
  if constexpr (!is_ball) {
    // F is decreasing for large t: the terms where F < 2^{-abs_bits} are
    // dropped (within the margin of abs_bits), and the series is not
    // accurate enough there to keep F positive
    auto small = [&](double u) -> bool {
      return QuadHeckeKernelValue(ker, real_from_double<Tf>(u), u) < tol;
    };
    if (small(u1)) {
      double lo = 0.0, hi = u1;
      if (!small(lo)) {
        for (int it = 0; it < 50; it++) {
          double mid = (lo + hi) / 2;
          if (small(mid)) {
            hi = mid;
          } else {
            lo = mid;
          }
        }
        u1 = hi;
        X_double = std::min(X_double, std::floor(std::exp(u1 + logA)));
      }
    }
  }
  if (X_double > X_limit) {
    std::cerr << "QHECKE: L(1-r) needs " << X_double
              << " terms (A^2 = " << Ld.A2 << "), beyond the limit " << X_limit
              << "\n";
    return {};
  }
  int64_t X = static_cast<int64_t>(X_double);
#ifdef DEBUG_QHECKE
  std::chrono::steady_clock::time_point tc0 = std::chrono::steady_clock::now();
#endif
  std::optional<QuadHeckeInterp<Tf>> I =
      MakeQuadHeckeInterp(ker, u0, u1, step, tol, D);
  if (!I) {
    return {};
  }
  QuadHeckeLExtend(Ld, X);
#ifdef DEBUG_QHECKE
  std::chrono::steady_clock::time_point tc1 = std::chrono::steady_clock::now();
  std::cerr << "QHECKE: data up to " << X << " (S=" << Ld.S << "), "
            << I->lo.size() << " segments, in "
            << std::chrono::duration<double>(tc1 - tc0).count() << " s\n";
#endif
  // the terms
  RealPrecisionGuard<Tf> g_e(eprec);
  std::vector<std::vector<Tf>> coef_e = QuadHeckeInterpCoefficients<Tf>(*I);
  std::vector<std::vector<double>> coef_d =
      QuadHeckeInterpCoefficients<double>(*I);
  std::vector<std::vector<DoubleDouble>> coef_dd =
      QuadHeckeInterpCoefficients<DoubleDouble>(*I);
  // 0: Tf, 1: DoubleDouble, 2: double, from the relative bits needed
  std::vector<int> band;
  for (double l2 : I->log2_fmax) {
    double needed = abs_bits + l2 + 8;
    band.push_back(is_ball ? 0
                           : (all_double || needed <= 45)
                                 ? 2
                                 : ((needed <= 85) ? 1 : 0));
  }
  Tf logAe = real_log(real_from_integer<Tf>(Ld.A2)) / real_from_int<Tf>(2);
  double logA_d = real_to_double(logAe);
  DoubleDouble logA_dd = real_convert<DoubleDouble>(logAe);
  Tf Lambda = real_from_int<Tf>(0);
  // the double and DoubleDouble terms, summed in DoubleDouble
  DoubleDouble sum_dd(0.0);
  std::vector<int64_t> abs_sum(I->lo.size(), 0);
  bool inside = true;
  int64_t n_tf = 0, n_dd = 0, n_dbl = 0;
  QuadHeckeForEachTerm(Ld, X,
                       [&](int64_t k, int64_t ak, DoubleDouble lk) -> void {
    double ud = lk.hi - logA_d;
    size_t seg = QuadHeckeSegment(*I, ud);
    if (band[seg] == 2) {
      double v = static_cast<double>(ak) *
                 QuadHeckeChebyshevValue(coef_d[seg], I->mid[seg],
                                         I->half[seg], ud);
      sum_dd += DoubleDouble(v);
      n_dbl++;
    } else if (band[seg] == 1) {
      DoubleDouble u = lk - logA_dd;
      sum_dd += DoubleDouble(static_cast<double>(ak)) *
                QuadHeckeChebyshevValue(coef_dd[seg], I->mid[seg],
                                        I->half[seg], u);
      n_dd++;
    } else {
      if constexpr (is_ball) {
        if (ud < I->lo[seg] - 1e-9 || ud > I->hi[seg] + 1e-9) {
          inside = false;
          return;
        }
        abs_sum[seg] += (ak < 0) ? -ak : ak;
      }
      Tf u = real_log(real_from_int<Tf>(k)) - logAe;
      Lambda += real_from_int<Tf>(ak) *
                QuadHeckeChebyshevValue(coef_e[seg], I->mid[seg],
                                        I->half[seg], u);
      n_tf++;
    }
  });
  if (!inside) {
    std::cerr << "QHECKE: a term outside the segments\n";
    return {};
  }
  Lambda += real_convert<Tf>(sum_dd);
  if constexpr (is_ball) {
    // the interpolation errors
    Tf err = real_from_int<Tf>(0);
    for (size_t i = 0; i < abs_sum.size(); i++) {
      err += real_from_int<Tf>(abs_sum[i]) * I->err[i];
    }
    // the terms k > X, with sigma optimized as above
    double best = 1e300, sbest = 3;
    for (double sigma = std::max(3.0, r + 0.0); sigma < 2000; sigma *= 1.02) {
      double lg = n * (-(sigma + delta) / 2 * std::log(M_PI) +
                       std::lgamma((sigma + delta) / 2));
      double v = lg + sigma * logA + (2 - sigma) * std::log(X_double);
      if (v < best) {
        best = v;
        sbest = sigma;
      }
    }
    Tf sg = real_from_double<Tf>(sbest);
    Tf tail = real_from_int<Tf>(2) * QuadHeckeGamma(sg, n, delta) *
              real_exp(sg * logAe) *
              real_exp((real_from_int<Tf>(2) - sg) *
                       real_log(real_from_int<Tf>(X))) *
              real_pow_int(real_zeta<Tf>(2), n);
    err += real_upper(tail);
    real_add_error(Lambda, real_upper(err));
  }
  Tf s1 = real_from_int<Tf>(1 - r);
  Tf res = Lambda / (real_exp(logAe * s1) * QuadHeckeGamma(s1, n, delta));
#ifdef DEBUG_QHECKE
  std::chrono::steady_clock::time_point tc2 = std::chrono::steady_clock::now();
  std::cerr << "QHECKE: " << n_tf << " terms in Tf, " << n_dd
            << " in DoubleDouble, " << n_dbl << " in double, in "
            << std::chrono::duration<double>(tc2 - tc1).count() << " s\n";
  std::cerr << "QHECKE: n=" << n << " delta=" << delta << " r=" << r
            << " A^2=" << Ld.A2 << " X=" << X << " J=" << J << " D=" << D
            << " wprec=" << wprec << " eprec=" << eprec
            << " L(1-r)=" << real_to_double(res) << "\n";
#endif
  (void)n_tf;
  (void)n_dd;
  (void)n_dbl;
  return res;
}

/*
  An integer mu divisible by the order of every finite subgroup of
  GL_m(K) (Hecke.jl, _minkowski_multiple, after Guralnick and Lorenz): for
  a prime P of K over p with residue field F_q, the reduction modulo P is
  injective on the elements of order prime to p, so the p'-part of the
  order divides |GL_m(F_q)| = q^{m(m-1)/2} prod_{i=1}^m (q^i - 1). The
  candidate starts at the bound from p = 2 for the odd part and from p = 3
  for the 2-part, and is reduced by the primes p in turn, its p-part being
  kept, until it has not changed for 20 primes. Every candidate is a valid
  multiple.
 */
template <typename Tint>
Tint MinkowskiMultiple(AbelianField const &FK, int m) {
  DirichletGroup G = ComputeDirichletGroup(FK.n);
  std::vector<std::vector<int64_t>> XK = FieldCharacters(G, FK.H, true);
  auto min_f = [&](int64_t p) -> int64_t {
    // the residue degree is the same at all the primes over p
    return PrimeDecompositionCharacters(G, XK, p).f;
  };
  auto glf = [&](int64_t p) -> Tint {
    Tint q(1);
    Tint pz = UniversalScalarConversion<Tint, int64_t>(p);
    q *= MyPow<Tint>(pz, min_f(p));
    Tint res(1);
    Tint qi(1);
    for (int i = 1; i <= m; i++) {
      qi *= q;
      res *= qi - 1;
    }
    return res;
  };
  auto ppart = [](Tint x, int64_t p) -> Tint {
    Tint pz = UniversalScalarConversion<Tint, int64_t>(p);
    Tint res(1);
    while (x != 0 && ResInt(x, pz) == 0) {
      x = QuoInt(x, pz);
      res *= pz;
    }
    return res;
  };
  Tint cand = glf(2) * ppart(glf(3), 2);
  int stab = 0;
  int64_t p = 3;
  auto is_prime = [](int64_t x) -> bool {
    for (int64_t i = 2; i * i <= x; i++) {
      if (x % i == 0) {
        return false;
      }
    }
    return x >= 2;
  };
  while (stab <= 20) {
    Tint old = cand;
    Tint bound = glf(p) * ppart(cand, p);
    cand = GcdPair(cand, bound);
    if (old == cand) {
      stab++;
    }
    p++;
    while (!is_prime(p)) {
      p++;
    }
  }
  return cand;
}

/*
  The standard mass when E is not abelian: 2^{-nr} prod_{i odd <= m-3}
  zeta_K(-i), exact, times L(1 - r, chi_{E/K}) computed numerically; the
  mass sm * local is made exact as in Hecke.jl: with mu from
  MinkowskiMultiple, the mass of the genus times mu is an integer (its
  terms are the 1/|O(L_i)|), so b = mu num(local) sm is an integer. The
  L-value is computed to the bits of b plus a margin and b rounded; it is
  accepted when a second computation with 16 more bits rounds to the same
  integer, both within 2^{-8} of it. With the method Arb, b is moreover
  certified: the Arb enclosure of L(1 - r) (QuadHeckeLNegative<ArbReal>)
  times mu num(local) sm_exact must contain b and no other integer,
  otherwise the mass is not returned. Nothing when the character cannot
  be computed (a common index divisor).
 */
template <typename T>
std::optional<T> QuadraticMassNonAbelianE(
    MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Omega,
    AbelianField const &FK, ElemE<T> const &d, NumberField<T> const &K,
    int m, T const &local, QuadLValueMethod method) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  using mpfr::mpreal;
  int r = m / 2;
  int64_t deg = AbelianFieldDegree(FK);
  T sm_exact(1);
  for (int64_t i = 0; i < deg * r; i++) {
    sm_exact /= 2;
  }
  for (int i = 1; i <= m - 3; i += 2) {
    sm_exact *= DedekindZetaAbelian<T>(FK, i + 1);
  }
  Tint mu = MinkowskiMultiple<Tint>(FK, m);
  T factor = sm_exact * UniversalScalarConversion<T, Tint>(mu) *
             UniversalScalarConversion<T, Tint>(GetNumerator_z(local));
  std::optional<QuadHeckeLData<T>> Ld =
      MakeQuadHeckeLData(B, Omega, FK, d, K, m);
  if (!Ld) {
    return {};
  }
  bool all_double = (method == QuadLValueMethod::Double);
  // the size of b, for which a few bits of L(1 - r) are enough
  std::optional<mpreal> L0 = QuadHeckeLNegative<mpreal>(*Ld, 24, all_double);
  if (!L0) {
    return {};
  }
  mpreal::set_default_prec(256);
  double log2_b =
      mpfr::log2(mpfr::abs(mpreal_from_rational(factor, 256) * *L0)).toDouble();
  int bits = std::max(48, static_cast<int>(std::ceil(log2_b)) + 24);
  if (all_double) {
    // about 50 bits in double precision: b of at most about 26 bits
    if (log2_b + 24 > 50) {
      std::cerr << "QHECKE: the integer to round has " << log2_b
                << " bits, too many for lvalue=double\n";
      return {};
    }
    bits = 42;
  }
  auto rounded = [&](int prec) -> std::optional<Tint> {
    std::optional<mpreal> L = QuadHeckeLNegative<mpreal>(*Ld, prec, all_double);
    if (!L) {
      return {};
    }
    mpfr_prec_t wp = prec + 64;
    mpreal bf = mpreal_from_rational(factor, wp) * *L;
    Tint b = mpreal_nearest_integer<Tint>(bf);
    if (mpfr::abs(bf - mpreal_from_integer(b, wp)) >
        mpfr::pow(mpreal(2, wp), -8)) {
#ifdef DEBUG_QHECKE
      std::cerr << "QHECKE: b = " << bf.toString(40)
                << " is not close to an integer\n";
#endif
      return {};
    }
    return b;
  };
  for (int iter = 0; iter < (all_double ? 1 : 4); iter++) {
    std::optional<Tint> b1 = rounded(bits);
    std::optional<Tint> b2 = rounded(bits + (all_double ? 8 : 16));
    if (b1 && b2 && *b1 == *b2) {
#ifdef ENABLE_FLINT_SUPPORT
      // the integer certified by the Arb enclosure of L(1 - r)
      bool certified = (method != QuadLValueMethod::Arb);
      for (int cb = bits; cb <= bits + 64 && !certified; cb += 32) {
        std::optional<ArbReal> La = QuadHeckeLNegative<ArbReal>(*Ld, cb, false);
        if (!La) {
          break;
        }
        ArbPrecisionGuard g(cb + 128);
        ArbReal bb = arb_from_rational(factor) * *La;
        std::optional<Tint> bc = arb_unique_integer<Tint>(bb);
# ifdef DEBUG_QHECKE
        std::cerr << "QHECKE_ARB: b in " << bb.str(20) << "\n";
# endif
        if (bc) {
          if (*bc != *b1) {
            std::cerr << "QHECKE_ARB: the certified b=" << *bc
                      << " differs from the numerical b=" << *b1 << "\n";
            throw TerminalException{1};
          }
          certified = true;
        }
      }
      if (!certified) {
        return {};
      }
#else
      if (method == QuadLValueMethod::Arb) {
        std::cerr << "QHECKE: the Arb certification needs flint\n";
        throw TerminalException{1};
      }
#endif
      T sm = UniversalScalarConversion<T, Tint>(*b1) / factor * sm_exact;
      // sm = b / (mu num(local)), the sign of the L-value included
      if (sm < 0) {
        sm = -sm;
      }
      return sm * local;
    }
    bits *= 2;
  }
  return {};
}

// clang-format off
#endif  // SRC_GENERA_MASSQUADRATICHECKEL_H_
// clang-format on
