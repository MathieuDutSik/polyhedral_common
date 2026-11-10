// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_QUADRATICHECKEL_H_
#define SRC_GENERA_QUADRATICHECKEL_H_

// clang-format off
#include "AbelianQuadraticMass.h"
#include "MprealConversion.h"
#include <cmath>
#include <chrono>
#include <functional>
#include <map>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

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
      if (a[i] == 0) {
        continue;
      }
      for (int j = 0; j < n; j++) {
        tmp[i + j] = add(tmp[i + j], mont_mul(a[i], b[j]));
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
  if (S % f_p != 0 || std::abs(S / f_p) > g) {
    std::cerr << "QHECKE: inconsistent trace " << S << " at p=" << p << "\n";
    throw TerminalException{1};
  }
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
MakeQuadHeckeGoodPrimes(AbelianCMField const &FK,
                        std::vector<typename underlying_z_ring<T>::ring_type> const &dint,
                        CyclotomicField<T> const &K) {
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
    AbelianCMField const &FK, ElemE<T> const &d_in, CyclotomicField<T> const &K,
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
    std::vector<T> pt = PowerTraces<T>(K.Phi, 2 * n);
    MyMatrix<T> Tm(n, n);
    for (int a = 0; a < n; a++) {
      for (int b = 0; b < n; b++) {
        Tm(a, b) = pt[a + b];
      }
    }
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
      PolyQ fq = PolyReduce(PolyQ(K.Phi.begin(), K.Phi.end()), p);
      PolyQ dq;
      for (auto const &c : dint) {
        Tint r = ResInt(c, UniversalScalarConversion<Tint, int64_t>(p));
        dq.push_back(UniversalScalarConversion<int64_t, Tint>(r));
      }
      for (auto const &[g, e] : FactorModQWithMultiplicity(fq, p)) {
        PolyQ r = PolyDivModQ(PolyReduce(dq, p), g, p).second;
        PolyTrim(r);
        int chi = ResidueQuadraticCharacter(r, g, p);
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
        CyclotomicField<T> Kp;
        Kp.d = 0;
        Kp.deg = n;
        Kp.Phi = lg.f;
        Kp.conj = lg.conj;
        ElemE<T> dp = DetK(H, Kp);
        if ((m * (m - 1) / 2) % 2 == 1) {
          dp = EScal(dp, T(-1));
        }
        return dp;
      };
      std::optional<LocalGenerator<T>> lg = FindLocalGenerator(B, Omega, p);
      if (lg) {
        ElemE<T> dp = disc_of(*lg, LocalOKGram(B, *lg, p));
        for (auto const &[g, e] :
             FactorModQWithMultiplicity(PolyQ(lg->f.begin(), lg->f.end()), p)) {
          bps.push_back({*lg, g, e, dp});
        }
      } else {
        // a common index divisor: a generator for each prime
        std::optional<std::vector<QuadPrimeGenerator<T>>> pgs =
            QuadPrimeGenerators(B, Omega, p,
                                PrimeDecompositionCharacters(G, XK, p));
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
            hp.chi = ResidueQuadraticCharacter(LocalResidue(u, P), g, p);
          }
        } else {
          ResidueFieldGF F{P.q, P.g, P.f_res};
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
        for (int i = 0; i < hp.cond_exp * fr; i++) {
          res.norm_conductor *= UniversalScalarConversion<Tint, int64_t>(p);
        }
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
struct QuadHeckeKernel {
  int n;
  int delta;
  mpfr_prec_t wprec;
  mpfr::mpreal s;
  mpfr::mpreal gs;
  mpfr::mpreal g1s;
  // q[k][j]: the coefficient of t^{delta + 2j} (log t)^k
  std::vector<std::vector<mpfr::mpreal>> q;
};

inline mpfr::mpreal QuadHeckeGamma(mpfr::mpreal const &s, int n, int delta) {
  mpfr::mpreal u = (s + delta) / 2;
  mpfr::mpreal g = mpfr::pow(mpfr::const_pi(s.get_prec()), -u) * mpfr::gamma(u);
  return mpfr::pow(g, n);
}

// log2 of C_j t^{2j}, in double precision
inline double QuadHeckeTermLog2(int n, int j, double t) {
  return (n * j * std::log(M_PI) + n * std::log(2.0) -
          n * std::lgamma(j + 1.0) + 2.0 * j * std::log(t)) /
         std::log(2.0);
}

inline QuadHeckeKernel MakeQuadHeckeKernel(int n, int delta, int s_int,
                                           double t_max, mpfr_prec_t wprec) {
  using mpfr::mpreal;
  mpreal::set_default_prec(wprec);
  QuadHeckeKernel ker;
  ker.n = n;
  ker.delta = delta;
  ker.wprec = wprec;
  ker.s = mpreal(s_int);
  mpreal s1 = 1 - ker.s;
  ker.gs = QuadHeckeGamma(ker.s, n, delta);
  ker.g1s = QuadHeckeGamma(s1, n, delta);
  // the number of terms: C_j t_max^{2j} below 2^{-wprec}, past the maximum
  int J = 1;
  while (J < 8 ||
         QuadHeckeTermLog2(n, J, t_max) > -static_cast<double>(wprec) - 20 ||
         J < 3.0 * M_PI * std::pow(t_max, 2.0 / n)) {
    J++;
  }
  mpreal pi = mpfr::const_pi(wprec);
  mpreal logpi = mpfr::log(pi);
  mpreal euler = mpfr::const_euler(wprec);
  std::vector<mpreal> zeta(n + 1);
  for (int k = 2; k <= n; k++) {
    zeta[k] = mpfr::zeta(mpreal(k));
  }
  std::vector<mpreal> H(n + 1, mpreal(0));
  mpreal C = mpfr::pow(mpreal(2), n);
  mpreal pin = mpfr::pow(pi, n);
  ker.q.assign(n, std::vector<mpreal>(J + 1, mpreal(0)));
  std::vector<mpreal> fact(n + 1);
  fact[0] = 1;
  for (int k = 1; k <= n; k++) {
    fact[k] = fact[k - 1] * k;
  }
  for (int j = 0; j <= J; j++) {
    if (j > 0) {
      mpreal ij = 1 / mpreal(j);
      mpreal pw = ij;
      for (int k = 1; k <= n; k++) {
        H[k] += pw;
        pw *= ij;
      }
      // C_j = C_{j-1} pi^n (-1/j)^n
      C *= pin * mpfr::pow(-ij, n);
    }
    // A(eps) = n S_j(eps/2) - eps (n/2) log pi, up to eps^{n-1}
    std::vector<mpreal> a(n, mpreal(0));
    mpreal half = mpreal(1) / 2;
    mpreal hk = half;
    for (int k = 1; k < n; k++) {
      mpreal sk = (k == 1)
                      ? mpreal(-euler + H[1])
                      : mpreal(((k % 2 == 0) ? zeta[k] : -zeta[k]) + H[k]) / k;
      a[k] = n * sk * hk;
      hk *= half;
    }
    if (n > 1) {
      a[1] -= logpi * n / 2;
    }
    // exp(A)
    std::vector<mpreal> e(n, mpreal(0));
    e[0] = 1;
    for (int k = 1; k < n; k++) {
      mpreal sum = 0;
      for (int i = 1; i <= k; i++) {
        sum += i * a[i] * e[k - i];
      }
      e[k] = sum / k;
    }
    // sum_k eps^k (1/(s - z)^{k+1} + 1/(1 - s - z)^{k+1})
    mpreal z = -delta - 2 * j;
    mpreal i1 = 1 / (ker.s - z);
    mpreal i2 = 1 / (s1 - z);
    std::vector<mpreal> h(n);
    mpreal p1 = i1, p2 = i2;
    for (int k = 0; k < n; k++) {
      h[k] = p1 + p2;
      p1 *= i1;
      p2 *= i2;
    }
    for (int k = 0; k < n; k++) {
      // g_{n-1-k}
      mpreal g = 0;
      for (int i = 0; i <= n - 1 - k; i++) {
        g += e[i] * h[n - 1 - k - i];
      }
      mpreal sign = (k % 2 == 0) ? mpreal(1) : mpreal(-1);
      ker.q[k][j] = C * g * sign / fact[k];
    }
  }
  return ker;
}

// G_s(t) + G_{1-s}(t)
inline mpfr::mpreal QuadHeckeKernelValue(QuadHeckeKernel const &ker,
                                         mpfr::mpreal const &t) {
  using mpfr::mpreal;
  mpreal L = mpfr::log(t);
  mpreal t2 = t * t;
  mpreal res =
      ker.gs * mpfr::pow(t, -ker.s) + ker.g1s * mpfr::pow(t, ker.s - 1);
  mpreal corr = 0;
  mpreal Lk = 1;
  int J = static_cast<int>(ker.q[0].size()) - 1;
  for (int k = 0; k < ker.n; k++) {
    mpreal poly = ker.q[k][J];
    for (int j = J - 1; j >= 0; j--) {
      poly = poly * t2 + ker.q[k][j];
    }
    corr += Lk * poly;
    Lk *= L;
  }
  if (ker.delta == 1) {
    corr *= t;
  }
  return res - corr;
}

/*
  Double-double numbers hi + lo, |lo| <= ulp(hi) / 2, about 106 bits, for
  the terms of the sum that need more than double precision but much less
  than the working precision (Dekker, Knuth; Hida, Li and Bailey).
 */
struct QuadDD {
  double hi;
  double lo;
};

inline QuadDD QuadDDTwoSum(double a, double b) {
  double s = a + b;
  double bb = s - a;
  double e = (a - (s - bb)) + (b - bb);
  return {s, e};
}

inline QuadDD QuadDDQuickTwoSum(double a, double b) {
  double s = a + b;
  return {s, b - (s - a)};
}

inline QuadDD QuadDDAdd(QuadDD a, QuadDD b) {
  QuadDD s = QuadDDTwoSum(a.hi, b.hi);
  QuadDD t = QuadDDTwoSum(a.lo, b.lo);
  s.lo += t.hi;
  s = QuadDDQuickTwoSum(s.hi, s.lo);
  s.lo += t.lo;
  return QuadDDQuickTwoSum(s.hi, s.lo);
}

inline QuadDD QuadDDMul(QuadDD a, QuadDD b) {
  double p = a.hi * b.hi;
  double e = std::fma(a.hi, b.hi, -p);
  e += a.hi * b.lo + a.lo * b.hi;
  return QuadDDQuickTwoSum(p, e);
}

inline QuadDD QuadDDFromMpreal(mpfr::mpreal const &x) {
  double hi = x.toDouble();
  double lo = (x - mpfr::mpreal(hi)).toDouble();
  return QuadDDQuickTwoSum(hi, lo);
}

inline mpfr::mpreal QuadDDToMpreal(QuadDD const &x) {
  return mpfr::mpreal(x.hi) + mpfr::mpreal(x.lo);
}

// exp(x) for |x| < 700: x = k log 2 + r, exp(r / 2^10) by Taylor, then 10
// squarings kept as expm1 (2 e + e^2), and the scaling by 2^k
inline QuadDD QuadDDExp(QuadDD x) {
  const QuadDD ln2{0.6931471805599453, 2.3190468138462996e-17};
  double k = std::nearbyint(x.hi / ln2.hi);
  QuadDD r = QuadDDAdd(x, QuadDDMul(ln2, QuadDD{-k, 0.0}));
  r.hi = std::ldexp(r.hi, -10);
  r.lo = std::ldexp(r.lo, -10);
  // expm1(r) = r + r^2/2 + ... + r^9/9!
  QuadDD e = r;
  QuadDD term = r;
  for (int i = 2; i <= 9; i++) {
    term = QuadDDMul(term, r);
    term.hi /= i;
    term.lo /= i;
    term = QuadDDQuickTwoSum(term.hi, term.lo);
    e = QuadDDAdd(e, term);
  }
  for (int i = 0; i < 10; i++) {
    // expm1(2y) = 2 expm1(y) + expm1(y)^2
    QuadDD two_e{2 * e.hi, 2 * e.lo};
    e = QuadDDAdd(two_e, QuadDDMul(e, e));
  }
  QuadDD res = QuadDDAdd(QuadDD{1.0, 0.0}, e);
  int ki = static_cast<int>(k);
  return {std::ldexp(res.hi, ki), std::ldexp(res.lo, ki)};
}

/*
  The kernel F(t) = G_s(t) + G_{1-s}(t) is positive (phi is), and is
  evaluated on the terms k = 1, ..., X through a Chebyshev interpolation
  of g(u) = log F(e^u) on segments of the variable u = log t: the series
  above is only computed at the nodes, at the working precision, and a
  segment is halved until the last two Chebyshev coefficients are below
  the tolerance, which bounds the error of F by 2^{-abs_bits}. A term only
  needs F to abs_bits + log2 F relative bits, so a segment is evaluated in
  double precision when that is at most 45 bits, in double-double when it
  is at most 85, and with MPFR otherwise: the terms are mostly in the
  first two bands.
 */
struct QuadHeckeInterp {
  mpfr_prec_t eprec;
  std::vector<double> lo;
  std::vector<mpfr::mpreal> mid;
  std::vector<mpfr::mpreal> half;
  std::vector<std::vector<mpfr::mpreal>> coef;
  // 0: MPFR, 1: double-double, 2: double
  std::vector<int> band;
  std::vector<double> mid_d;
  std::vector<double> half_d;
  std::vector<std::vector<double>> coef_d;
  std::vector<QuadDD> mid_dd;
  std::vector<QuadDD> half_dd;
  std::vector<std::vector<QuadDD>> coef_dd;
};

inline QuadHeckeInterp MakeQuadHeckeInterp(QuadHeckeKernel const &ker,
                                           double u0, double u1,
                                           double abs_bits, mpfr_prec_t eprec) {
  using mpfr::mpreal;
  const int D = 30;
  QuadHeckeInterp res;
  res.eprec = eprec;
  mpreal::set_default_prec(ker.wprec);
  mpreal pi = mpfr::const_pi(ker.wprec);
  std::function<void(double, double, int)> build = [&](double a, double b,
                                                       int depth) -> void {
    mpreal::set_default_prec(ker.wprec);
    mpreal m = (mpreal(a) + mpreal(b)) / 2;
    mpreal h = (mpreal(b) - mpreal(a)) / 2;
    std::vector<mpreal> g(D + 1);
    double gmax = -1e300;
    for (int i = 0; i <= D; i++) {
      mpreal x = mpfr::cos(pi * (2 * i + 1) / (2 * (D + 1)));
      mpreal F = QuadHeckeKernelValue(ker, mpfr::exp(m + h * x));
      if (F <= 0) {
        std::cerr << "QHECKE: the kernel is not positive at u=" << m + h * x
                  << ", F=" << F << "\n";
        throw TerminalException{1};
      }
      g[i] = mpfr::log(F);
      gmax = std::max(gmax, g[i].toDouble());
    }
    std::vector<mpreal> c(D + 1);
    for (int k = 0; k <= D; k++) {
      mpreal sum = 0;
      for (int i = 0; i <= D; i++) {
        sum += g[i] * mpfr::cos(pi * k * (2 * i + 1) / (2 * (D + 1)));
      }
      c[k] = sum * 2 / (D + 1);
    }
    c[0] /= 2;
    // the error on g allowed: 2^{-abs_bits} / max F
    double log2_tol = -abs_bits - gmax / std::log(2.0) - 4;
    double log2_tail =
        mpfr::log2(mpfr::abs(c[D]) + mpfr::abs(c[D - 1]) +
                   mpreal(2, ker.wprec) *
                       mpfr::pow(mpreal(2), -static_cast<long>(ker.wprec)))
            .toDouble();
    if (log2_tail > log2_tol) {
      if (depth > 40) {
        std::cerr << "QHECKE: no convergence of the interpolation on [" << a
                  << ", " << b << "]\n";
        throw TerminalException{1};
      }
      double mm = (a + b) / 2;
      build(a, mm, depth + 1);
      build(mm, b, depth + 1);
      return;
    }
    auto at_eprec = [&](mpreal x) -> mpreal {
      x.setPrecision(static_cast<int>(eprec), MPFR_RNDN);
      return x;
    };
    res.lo.push_back(a);
    res.mid.push_back(at_eprec(m));
    res.half.push_back(at_eprec(h));
    std::vector<mpreal> ce(D + 1);
    for (int k = 0; k <= D; k++) {
      ce[k] = at_eprec(c[k]);
    }
    res.coef.push_back(ce);
    std::vector<double> cd(D + 1);
    for (int k = 0; k <= D; k++) {
      cd[k] = c[k].toDouble();
    }
    res.coef_d.push_back(cd);
    res.mid_d.push_back((a + b) / 2);
    res.half_d.push_back((b - a) / 2);
    std::vector<QuadDD> cdd(D + 1);
    for (int k = 0; k <= D; k++) {
      cdd[k] = QuadDDFromMpreal(c[k]);
    }
    res.coef_dd.push_back(cdd);
    res.mid_dd.push_back(QuadDDFromMpreal(m));
    res.half_dd.push_back(QuadDDFromMpreal(h));
    // the relative bits needed, with a margin
    double needed = abs_bits + gmax / std::log(2.0) + 8;
    res.band.push_back((needed <= 45) ? 2 : ((needed <= 85) ? 1 : 0));
  };
  for (double a = u0; a < u1; a += 1.0) {
    build(a, std::min(a + 1.0, u1), 0);
  }
  return res;
}

// the segment of u: the last one starting at or below it
inline size_t QuadHeckeSegment(QuadHeckeInterp const &I, double ud) {
  size_t i = static_cast<size_t>(
      std::upper_bound(I.lo.begin(), I.lo.end(), ud) - I.lo.begin());
  return (i == 0) ? 0 : i - 1;
}

// F(e^u) on the segment seg
inline mpfr::mpreal QuadHeckeInterpValue(QuadHeckeInterp const &I,
                                         mpfr::mpreal const &u, size_t seg) {
  using mpfr::mpreal;
  std::vector<mpreal> const &c = I.coef[seg];
  mpreal x = (u - I.mid[seg]) / I.half[seg];
  mpreal x2 = 2 * x;
  mpreal b1 = 0, b2 = 0;
  for (int k = static_cast<int>(c.size()) - 1; k >= 1; k--) {
    mpreal b0 = x2 * b1 - b2 + c[k];
    b2 = b1;
    b1 = b0;
  }
  return mpfr::exp(x * b1 - b2 + c[0]);
}

// the same in double-double
inline QuadDD QuadHeckeInterpValueDD(QuadHeckeInterp const &I, QuadDD u,
                                     size_t seg) {
  std::vector<QuadDD> const &c = I.coef_dd[seg];
  // x = (u - m) / h, h a power of 2 times a double
  QuadDD diff = QuadDDAdd(u, QuadDD{-I.mid_dd[seg].hi, -I.mid_dd[seg].lo});
  double inv = 1.0 / I.half_dd[seg].hi;
  QuadDD x = QuadDDMul(diff, QuadDD{inv, 0.0});
  // one Newton step for the division: x += (diff - x h) / h
  QuadDD xh = QuadDDMul(x, I.half_dd[seg]);
  QuadDD r = QuadDDAdd(diff, QuadDD{-xh.hi, -xh.lo});
  x = QuadDDAdd(x, QuadDDMul(r, QuadDD{inv, 0.0}));
  QuadDD x2{2 * x.hi, 2 * x.lo};
  QuadDD b1{0, 0}, b2{0, 0};
  for (int k = static_cast<int>(c.size()) - 1; k >= 1; k--) {
    QuadDD b0 = QuadDDAdd(QuadDDAdd(QuadDDMul(x2, b1), QuadDD{-b2.hi, -b2.lo}),
                          c[k]);
    b2 = b1;
    b1 = b0;
  }
  QuadDD g = QuadDDAdd(QuadDDAdd(QuadDDMul(x, b1), QuadDD{-b2.hi, -b2.lo}),
                       c[0]);
  return QuadDDExp(g);
}

// the same in double precision
inline double QuadHeckeInterpValueDouble(QuadHeckeInterp const &I, double u,
                                         size_t seg) {
  std::vector<double> const &c = I.coef_d[seg];
  double x = (u - I.mid_d[seg]) / I.half_d[seg];
  double x2 = 2 * x;
  double b1 = 0, b2 = 0;
  for (int k = static_cast<int>(c.size()) - 1; k >= 1; k--) {
    double b0 = x2 * b1 - b2 + c[k];
    b2 = b1;
    b1 = b0;
  }
  return std::exp(x * b1 - b2 + c[0]);
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
  AbelianCMField FK;
  ElemE<T> d;
  CyclotomicField<T> K;
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
  std::vector<QuadDD> log_small;
  std::vector<QuadDD> log_j;
};

template <typename T>
std::optional<QuadHeckeLData<T>>
MakeQuadHeckeLData(MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Omega,
                   AbelianCMField const &FK, ElemE<T> const &d,
                   CyclotomicField<T> const &K, int m) {
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
      Ld.log_small.push_back(QuadDDFromMpreal(
          mpfr::log(mpfr::mpreal(static_cast<unsigned long>(pc.first)))));
    }
    Ld.log_j.assign(S + 1, QuadDD{0, 0});
    for (int64_t j = 2; j <= S; j++) {
      Ld.log_j[j] = QuadDDFromMpreal(
          mpfr::log(mpfr::mpreal(static_cast<unsigned long>(j))));
    }
    mpfr::mpreal::set_default_prec(old);
  }
  Ld.X = X;
  Ld.S = S;
}

// log q to double-double precision, q < 2^53: y = log q in double, then
// log q = y + log(q e^{-y}) with q e^{-y} = 1 + delta, |delta| < 2^{-50}
inline QuadDD QuadDDLogInteger(int64_t q) {
  double y = std::log(static_cast<double>(q));
  QuadDD e = QuadDDExp(QuadDD{-y, 0.0});
  QuadDD t = QuadDDMul(e, QuadDD{static_cast<double>(q), 0.0});
  QuadDD delta = QuadDDAdd(t, QuadDD{-1.0, 0.0});
  return QuadDDAdd(QuadDD{y, 0.0}, delta);
}

/*
  f(k, a_k, log k) for the k <= X with a_k != 0, X <= Ld.X, in no
  particular order, log k in double-double: the S-smooth k by a depth-first walk over the local series of the
  primes p <= S (skipping the zero coefficients and what lies under them),
  and the k = j q with a prime q > S (at most one, as S^2 >= X) and j <= X /
  q <= S from a segmented sieve over (S, X], the primes with c_1(q) = 0
  being skipped.
 */
template <typename T, typename F>
void QuadHeckeForEachTerm(QuadHeckeLData<T> &Ld, int64_t X, F &&f) {
  int64_t S = Ld.S;
  std::function<void(size_t, int64_t, int64_t, QuadDD)> dfs =
      [&](size_t i0, int64_t k, int64_t ak, QuadDD lk) -> void {
    for (size_t i = i0; i < Ld.small.size(); i++) {
      int64_t p = Ld.small[i].first;
      if (k > X / p) {
        break;
      }
      std::vector<int64_t> const &c = Ld.small[i].second;
      int64_t kk = k;
      QuadDD lkk = lk;
      for (size_t e = 1; kk <= X / p; e++) {
        kk *= p;
        lkk = QuadDDAdd(lkk, Ld.log_small[i]);
        if (e < c.size() && c[e] != 0) {
          f(kk, ak * c[e], lkk);
          dfs(i + 1, kk, ak * c[e], lkk);
        }
      }
    }
  };
  f(int64_t(1), int64_t(1), QuadDD{0.0, 0.0});
  dfs(0, 1, 1, QuadDD{0.0, 0.0});
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
      QuadDD lq = QuadDDLogInteger(q);
      for (int64_t j = 1; j <= jmax; j++) {
        int64_t aj = Ld.a_small[j];
        if (aj != 0) {
          f(j * q, aj * c, QuadDDAdd(lq, Ld.log_j[j]));
        }
      }
    }
  }
}

/*
  L(1 - r, chi) to about prec_bits relative bits: Lambda(1 - r) from the
  sum, divided by A^{1-r} gamma(1-r) (no pole there: (1 - r + delta)/2 is
  a half-integer as r = delta + 1 mod 2). Nothing when the number of terms
  exceeds the limit (10^12).
 */
template <typename T>
std::optional<mpfr::mpreal> QuadHeckeLNegative(QuadHeckeLData<T> &Ld,
                                               int prec_bits) {
  using mpfr::mpreal;
  int r = Ld.m / 2;
  int n = Ld.n;
  int delta = Ld.delta;
  mpfr_prec_t old_prec = mpreal::get_default_prec();
  mpfr_prec_t p0 = 64;
  mpreal::set_default_prec(p0);
  mpreal A2 = mpreal_from_integer(Ld.A2, p0);
  double logA = 0.5 * mpfr::log(A2).toDouble();
  // |Lambda(1 - r)| = A^r gamma(r) |L(r)|, with L(r) of size about 1
  double log2_scale =
      (r * logA) / std::log(2.0) +
      mpfr::log2(mpfr::abs(QuadHeckeGamma(mpreal(r, p0), n, delta))).toDouble();
  double target = std::max(20.0, prec_bits - log2_scale + 10);
  // the cutoff: n pi t^{2/n} - (r + n) log t > target log 2 + log(A t)
  double t_max = 1.0;
  auto tail_ok = [&](double t) -> bool {
    return n * M_PI * std::pow(t, 2.0 / n) - (r + n) * std::log(t) >=
           (target + 40) * std::log(2.0) + logA + std::log(t) + 10;
  };
  while (!tail_ok(t_max)) {
    t_max *= 1.05;
  }
  double X_double = std::ceil(t_max * std::exp(logA)) + 1;
  // each term to 2^{-abs_bits}, the X terms then add up to 2^{-target}
  double abs_bits = target + std::log2(X_double) + 6;
  double log2_peak = n * M_PI * std::pow(t_max, 2.0 / n) / std::log(2.0);
  double log2_Fmax = std::max(0.0, log2_scale) + 8;
  mpfr_prec_t wprec =
      static_cast<mpfr_prec_t>(abs_bits + std::max(log2_peak, log2_Fmax) + 64);
  mpfr_prec_t eprec = static_cast<mpfr_prec_t>(abs_bits + log2_Fmax + 32);
  QuadHeckeKernel ker = MakeQuadHeckeKernel(n, delta, r, t_max, wprec);
  double u0 = -logA - 0.01;
  double u1 = std::log(X_double) - logA + 0.01;
  {
    // F is decreasing for large t: the terms where F < 2^{-abs_bits} are
    // dropped (within the margin of abs_bits), and the series is not
    // accurate enough there to keep F positive
    mpreal::set_default_prec(wprec);
    mpreal thr = mpfr::pow(mpreal(2), -abs_bits);
    auto small = [&](double u) -> bool {
      return QuadHeckeKernelValue(ker, mpfr::exp(mpreal(u))) < thr;
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
  // beyond this the sum is out of reach: the mass is left not computed
  const double X_limit = 1.0e12;
  if (X_double > X_limit) {
    std::cerr << "QHECKE: L(1-r) needs " << X_double
              << " terms (A^2 = " << Ld.A2 << "), beyond the limit " << X_limit
              << "\n";
    mpreal::set_default_prec(old_prec);
    return {};
  }
  int64_t X = static_cast<int64_t>(X_double);
#ifdef DEBUG_QHECKE
  std::chrono::steady_clock::time_point tc0 = std::chrono::steady_clock::now();
#endif
  QuadHeckeLExtend(Ld, X);
  QuadHeckeInterp I = MakeQuadHeckeInterp(ker, u0, u1, abs_bits, eprec);
#ifdef DEBUG_QHECKE
  std::chrono::steady_clock::time_point tc1 = std::chrono::steady_clock::now();
  std::cerr << "QHECKE: data up to " << X << " (S=" << Ld.S << ") in "
            << std::chrono::duration<double>(tc1 - tc0).count() << " s\n";
#endif
  mpreal::set_default_prec(eprec);
  mpreal logAe = mpfr::log(mpreal_from_integer(Ld.A2, eprec)) / 2;
  double logA_d = logAe.toDouble();
  mpreal Lambda = 0;
  QuadDD logA_dd = QuadDDFromMpreal(logAe);
  QuadDD mlogA_dd{-logA_dd.hi, -logA_dd.lo};
  // the double and double-double terms, summed in double-double
  QuadDD sum_dd{0, 0};
  int64_t n_mpfr = 0, n_dd = 0, n_dbl = 0;
  QuadHeckeForEachTerm(Ld, X, [&](int64_t k, int64_t ak, QuadDD lk) -> void {
    double ud = lk.hi - logA_d;
    size_t seg = QuadHeckeSegment(I, ud);
    if (I.band[seg] == 2) {
      double v =
          static_cast<double>(ak) * QuadHeckeInterpValueDouble(I, ud, seg);
      sum_dd = QuadDDAdd(sum_dd, QuadDD{v, 0.0});
      n_dbl++;
    } else if (I.band[seg] == 1) {
      QuadDD u = QuadDDAdd(lk, mlogA_dd);
      QuadDD v = QuadDDMul(QuadHeckeInterpValueDD(I, u, seg),
                           QuadDD{static_cast<double>(ak), 0.0});
      sum_dd = QuadDDAdd(sum_dd, v);
      n_dd++;
    } else {
      mpreal u = mpfr::log(mpreal(static_cast<unsigned long>(k))) - logAe;
      Lambda += ak * QuadHeckeInterpValue(I, u, seg);
      n_mpfr++;
    }
  });
  Lambda += QuadDDToMpreal(sum_dd);
#ifdef DEBUG_QHECKE
  std::chrono::steady_clock::time_point tc2 = std::chrono::steady_clock::now();
  std::cerr << "QHECKE: " << n_mpfr << " terms with MPFR, " << n_dd
            << " in double-double, " << n_dbl << " in double, in "
            << std::chrono::duration<double>(tc2 - tc1).count() << " s\n";
#endif
  mpreal s1 = mpreal(1 - r);
  mpreal res = Lambda / (mpfr::exp(logAe * s1) * QuadHeckeGamma(s1, n, delta));
#ifdef DEBUG_QHECKE
  std::cerr << "QHECKE: n=" << n << " delta=" << delta << " r=" << r
            << " A^2=" << Ld.A2 << " t_max=" << t_max << " X=" << X
            << " J=" << ker.q[0].size() << " segments=" << I.lo.size()
            << " wprec=" << wprec << " eprec=" << eprec
            << " L(1-r)=" << res.toString(30) << "\n";
#endif
  mpreal::set_default_prec(old_prec);
  return res;
}

// The certified evaluation with Arb, when flint is there
#ifdef ENABLE_FLINT_SUPPORT
#include "QuadraticHeckeLArb.h"
#endif

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
Tint MinkowskiMultiple(AbelianCMField const &FK, int m) {
  DirichletGroup G = ComputeDirichletGroup(FK.n);
  std::vector<std::vector<int64_t>> XK = FieldCharacters(G, FK.H, true);
  auto min_f = [&](int64_t p) -> int64_t {
    // the residue degree is the same at all the primes over p
    return PrimeDecompositionCharacters(G, XK, p).f;
  };
  auto glf = [&](int64_t p) -> Tint {
    Tint q(1);
    Tint pz = UniversalScalarConversion<Tint, int64_t>(p);
    for (int64_t i = 0; i < min_f(p); i++) {
      q *= pz;
    }
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
  certified: the Arb enclosure of L(1 - r) (QuadraticHeckeLArb.h) times mu
  num(local) sm_exact must contain b and no other integer, otherwise the
  mass is not returned. Nothing when the character cannot be computed (a common
  index divisor).
 */
template <typename T>
std::optional<T> QuadraticMassNonAbelianE(
    MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Omega,
    AbelianCMField const &FK, ElemE<T> const &d, CyclotomicField<T> const &K,
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
  // the size of b, for which a few bits of L(1 - r) are enough
  std::optional<mpreal> L0 = QuadHeckeLNegative(*Ld, 24);
  if (!L0) {
    return {};
  }
  mpreal::set_default_prec(256);
  double log2_b =
      mpfr::log2(mpfr::abs(mpreal_from_rational(factor, 256) * *L0)).toDouble();
  int bits = std::max(48, static_cast<int>(std::ceil(log2_b)) + 24);
  auto rounded = [&](int prec) -> std::optional<Tint> {
    std::optional<mpreal> L = QuadHeckeLNegative(*Ld, prec);
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
  for (int iter = 0; iter < 4; iter++) {
    std::optional<Tint> b1 = rounded(bits);
    std::optional<Tint> b2 = rounded(bits + 16);
    if (b1 && b2 && *b1 == *b2) {
#ifdef ENABLE_FLINT_SUPPORT
      // the integer certified by the Arb enclosure of L(1 - r)
      bool certified = (method == QuadLValueMethod::Mpfr);
      for (int cb = bits; cb <= bits + 64 && !certified; cb += 32) {
        std::optional<ArbReal> La = QuadHeckeLNegativeArb(*Ld, cb);
        if (!La) {
          break;
        }
        ArbPrecisionGuard g(cb + 128);
        ArbReal bb = arb_from_rational(factor) * *La;
        std::optional<Tint> bc = arb_unique_integer<Tint>(bb);
#ifdef DEBUG_QHECKE
        std::cerr << "QHECKE_ARB: b in " << bb.str(20) << "\n";
#endif
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
#endif  // SRC_GENERA_QUADRATICHECKEL_H_
// clang-format on
