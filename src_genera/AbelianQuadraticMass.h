// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ABELIANQUADRATICMASS_H_
#define SRC_GENERA_ABELIANQUADRATICMASS_H_

// clang-format off
#include "AbelianHermitianMass.h"
#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <vector>
// clang-format on

/*
  Masses of quadratic lattices over a totally real abelian field K, the
  subfield of Q(zeta_n) fixed by a subgroup H containing -1 (Hecke.jl,
  QuadForm/Quad/MassQuad.jl, after Siegel and Kirschmer). The standard mass
  of a lattice of rank m, r = floor(m / 2), is

      2^{-[K:Q] r} prod_{i odd, i <= m - 2} zeta_K(-i)            (m odd)
      2^{-[K:Q] r} prod_{i odd, i <= m - 3} zeta_K(-i) L_r        (m even)

  with L_r = zeta_K(1 - r) when the discriminant (-1)^{m(m-1)/2} det is a
  square in K, and L(1 - r, chi_{E/K}) for E = K(sqrt disc) otherwise. As
  K is abelian, zeta_K(1 - k) is the product of the L(1 - k, chi) over the
  Dirichlet characters trivial on H; when E is abelian too, given as the
  subfield of Q(zeta_n) fixed by H_E (of index 2 in H), L(1 - r, chi_{E/K})
  is the product over the characters trivial on H_E but not on H.
 */

inline mpq_class DedekindZetaAbelian(AbelianCMField const &F, int k) {
  // zeta_K(1 - k): all the characters of K, even as K is real
  return ProductDirichletLValues(F.n, k, 0, F.H);
}

// The quadratic field data: a totally real field as (n, H) with -1 in H.
inline AbelianCMField MakeAbelianRealField(int64_t n,
                                           std::vector<int64_t> const &gens) {
  std::vector<int64_t> g = gens;
  if (n > 2) {
    g.push_back(n - 1);
  }
  // MakeAbelianCMField refuses -1: build the subgroup directly
  std::set<int64_t> Hs{1 % std::max<int64_t>(n, 1)};
  std::vector<int64_t> queue{1 % std::max<int64_t>(n, 1)};
  while (!queue.empty()) {
    int64_t a = queue.back();
    queue.pop_back();
    for (int64_t x : g) {
      int64_t b = ((a * (x % n)) % n + n) % n;
      if (Hs.insert(b).second) {
        queue.push_back(b);
      }
    }
  }
  return {n, std::vector<int64_t>(Hs.begin(), Hs.end())};
}

/*
  The standard mass; for m even with a non-square discriminant, FE gives
  E = K(sqrt disc) as the subfield of Q(zeta_n) fixed by H_E (same n).
  Nothing when E is needed and not given.
 */
inline std::optional<mpq_class>
StandardQuadraticMassAbelian(AbelianCMField const &FK, int m, bool disc_square,
                             std::optional<AbelianCMField> const &FE) {
  if (m == 0) {
    return mpq_class(1);
  }
  if (m == 1) {
    return mpq_class(1, 2);
  }
  int64_t deg = AbelianFieldDegree(FK);
  int r = m / 2;
  mpq_class res(1);
  for (int64_t i = 0; i < deg * r; i++) {
    res /= 2;
  }
  int top = (m % 2 == 1) ? m - 2 : m - 3;
  for (int i = 1; i <= top; i += 2) {
    res *= DedekindZetaAbelian(FK, i + 1);
  }
  if (m % 2 == 0) {
    if (disc_square) {
      res *= DedekindZetaAbelian(FK, r);
    } else {
      if (!FE || FE->n != FK.n) {
        return {};
      }
      // the characters of E not of K, whatever their parity
      mpq_class prod(1);
      for (int parity = 0; parity <= 1; parity++) {
        prod *= ProductDirichletLValues(FK.n, r, parity, FE->H, FK.H);
      }
      res *= prod;
    }
  }
  return (res < 0) ? mpq_class(-res) : res;
}

/*
  The local factors (Hecke.jl, local_factor of MassQuad.jl). The lattice
  is given as for the Hermitian masses: the Gram matrix B of the trace form
  b = Tr_{K/Q} phi on the Z-lattice X, and the actions Omega of a Z-basis
  of O_K. At a prime q, with theta_q a generator of K with q prime to the
  index, phi is written over K on m vectors of X spanning X (x) Z_q over
  O_K (x) Z_q (LocalOKGram), and each prime p = (q, g(theta_q)) of K is
  handled on that Gram matrix with the exact arithmetic of K, v_p and the
  residues being read through the idempotent of p (LocalPrimeK), whatever
  the ramification and the number of primes over q.

  At an odd prime the factor only depends on the Jordan components: their
  scales s_i = v_p, ranks m_i and the square classes of their determinants
  divided by pi^{s_i m_i} (QuadJordanBlock). At a dyadic prime unramified
  over 2 it is Cho's formula (LocalFactorCho), and at a dyadic prime
  ramified over 2 see LocalFactorDyadicRamified. The common index divisors
  are not computed.
 */
struct QuadJordanBlock {
  int s;
  int m;
  int chi;
};

// Hecke's local factor at an odd prime from the Jordan components, the
// residue field having qf elements.
inline mpq_class QuadLocalFactorOdd(std::vector<QuadJordanBlock> blocks,
                                    int64_t qf) {
  std::sort(blocks.begin(), blocks.end(),
            [](QuadJordanBlock const &a, QuadJordanBlock const &b) {
              return a.s < b.s;
            });
  if (blocks.size() <= 1) {
    return mpq_class(1);
  }
  int chi_m1 = (qf % 4 == 1) ? 1 : -1;
  auto sign_class = [&](int k) -> int { return (k % 2 == 0) ? 1 : chi_m1; };
  int m = 0;
  int64_t vd = 0;
  int chid = 1;
  for (auto const &bl : blocks) {
    m += bl.m;
    vd += static_cast<int64_t>(bl.s) * bl.m;
    chid *= bl.chi;
  }
  chid *= sign_class(m * (m - 1) / 2);
  mpq_class f;
  if (m % 2 == 1) {
    f = HeckeGroupOrderO(m, 1, qf);
  } else if (vd % 2 != 0) {
    f = HeckeGroupOrderO(m - 1, 1, qf);
  } else {
    f = HeckeGroupOrderO(m, chid, qf);
  }
  // twice the exponent N
  int64_t N2 = 0;
  int after = m;
  for (auto const &bl : blocks) {
    int64_t s = bl.s, mi = bl.m;
    after -= bl.m;
    int64_t ri = after;
    int sq = ((s * mi) % 2 != 0) ? -1 : sign_class(bl.m * (bl.m - 1) / 2) * bl.chi;
    f /= HeckeGroupOrderO(bl.m, sq, qf);
    N2 -= s * mi * (mi + 1) + 2 * s * mi * ri;
    N2 += s * mi * (m + 1);
  }
  if (m % 2 == 0 && vd % 2 != 0) {
    N2 += 1 - m;
  }
  if (N2 % 2 != 0) {
    std::cerr << "ABQUAD: the exponent N is not integral\n";
    throw TerminalException{1};
  }
  int64_t N = N2 / 2;
  mpq_class qpow(1);
  for (int64_t i = 0; i < (N >= 0 ? N : -N); i++) {
    qpow *= mpq_class(static_cast<long>(qf));
  }
  return (N >= 0) ? mpq_class(qpow * f) : mpq_class(f / qpow);
}

/*
  The arithmetic at a prime p = (q, g(theta)) of K, theta = theta_q with q
  prime to the index: the elements of K are kept exactly on the powers of
  theta, v_p(x) is v_q of the norm of eps x + 1 - eps over f (eps the
  idempotent of p in Z_q[theta], known modulo q^N for an N raised until
  the norm is nonzero modulo q^N), and the residue of x in F_q[x] / g is
  read on eps x. This holds whatever the other primes over q are.
 */
struct LocalPrimeK {
  CyclotomicField K;
  int64_t q;
  PolyQ g;
  int e;
  int f_res;
  int64_t qf;
  ElemE pi;
  std::map<int, PolyZ> eps;
};

inline LocalPrimeK MakeLocalPrimeK(LocalGenerator const &lg, PolyQ const &g,
                                   int e, int64_t q) {
  LocalPrimeK P;
  int deg = static_cast<int>(lg.f.size()) - 1;
  P.K.d = 0;
  P.K.deg = deg;
  P.K.Phi = lg.f;
  P.K.conj = lg.conj;
  P.q = q;
  P.g = g;
  P.e = e;
  P.f_res = static_cast<int>(g.size()) - 1;
  P.qf = 1;
  for (int i = 0; i < P.f_res; i++) {
    if (P.qf > (int64_t(1) << 62) / q) {
      std::cerr << "ABQUAD: the residue field of order " << q << "^"
                << P.f_res << " is too large\n";
      throw TerminalException{1};
    }
    P.qf *= q;
  }
  P.pi = ElemE(deg, mpq_class(0));
  if (e == 1) {
    P.pi[0] = mpq_class(static_cast<long>(q));
  } else {
    ElemE pw = EZetaPower(0, P.K);
    for (size_t i = 0; i < g.size(); i++) {
      for (int t = 0; t < deg; t++) {
        P.pi[t] += mpq_class(static_cast<long>(g[i])) * pw[t];
      }
      pw = EMul(pw, EZetaPower(1, P.K), P.K);
    }
  }
  return P;
}

inline PolyZ const &LocalEps(LocalPrimeK &P, int N) {
  auto it = P.eps.find(N);
  if (it == P.eps.end()) {
    it = P.eps.emplace(N, AbelianPrimeIdempotentPoly(P.K.Phi, P.g, P.e, P.q, N))
             .first;
  }
  return it->second;
}

// x = y / d with y in Z[theta] and d > 0
inline std::pair<PolyZ, mpz_class> ClearDenominators(ElemE const &x) {
  mpz_class d(1);
  for (auto const &c : x) {
    d = lcm(d, c.get_den());
  }
  PolyZ y;
  for (auto const &c : x) {
    mpq_class t = c * mpq_class(d);
    y.push_back(t.get_num());
  }
  return {y, d};
}

inline int ValuationMpz(mpz_class x, int64_t q) {
  int v = 0;
  mpz_class qz(static_cast<long>(q));
  while (x != 0 && x % qz == 0) {
    x /= qz;
    v++;
  }
  return v;
}

inline int LocalVal(ElemE const &x, LocalPrimeK &P) {
  if (EIsZero(x)) {
    return 1 << 28;
  }
  auto [y, d] = ClearDenominators(x);
  int k = ValuationMpz(d, P.q);
  int N = P.f_res * (P.e * k + 4) + 4;
  while (true) {
    mpz_class Q = PowInt(P.q, N);
    PolyZ const &eps = LocalEps(P, N);
    PolyZ w = PolyZReduce(PolyZMul(eps, y), P.K.Phi, Q);
    for (size_t i = 0; i < eps.size(); i++) {
      w[i] -= eps[i];
    }
    w[0] += 1;
    w = PolyZReduce(w, P.K.Phi, Q);
    ElemE we;
    for (auto const &c : w) {
      we.push_back(mpq_class(c));
    }
    mpq_class det = DeterminantMat(EMulMatrix(we, P.K));
    mpz_class dn = det.get_num();
    if (dn % Q != 0) {
      int v = ValuationMpz(dn, P.q);
      return v / P.f_res - P.e * k;
    }
    N *= 2;
  }
}

// The residue in F_q[x] / g of x with v_p(x) >= 0.
inline PolyQ LocalResidue(ElemE const &x, LocalPrimeK &P) {
  auto [y, d] = ClearDenominators(x);
  int k = ValuationMpz(d, P.q);
  mpz_class qk = PowInt(P.q, k);
  mpz_class Q = qk * mpz_class(static_cast<long>(P.q));
  PolyZ z = PolyZReduce(PolyZMul(LocalEps(P, k + 1), y), P.K.Phi, Q);
  mpz_class qz(static_cast<long>(P.q));
  mpz_class dd = (d / qk) % qz;
  int64_t inv = InvModQ(dd.get_si(), P.q);
  PolyQ res;
  for (auto const &c : z) {
    if (c % qk != 0) {
      std::cerr << "ABQUAD: the residue of a non integral element\n";
      throw TerminalException{1};
    }
    mpz_class t = (c / qk) % qz;
    res.push_back(MulModQ(t.get_si(), inv, P.q));
  }
  PolyTrim(res);
  if (res.empty()) {
    return {};
  }
  PolyQ r = PolyDivModQ(res, P.g, P.q).second;
  PolyTrim(r);
  return r;
}

/*
  The Gram matrix over K of phi on m vectors of X spanning X (x) Z_q over
  O_K (x) Z_q: standard vectors first, then small combinations.
 */
inline std::vector<std::vector<ElemE>>
LocalOKGram(MyMatrix<mpq_class> const &B, LocalGenerator const &lg, int64_t q) {
  int r = B.rows();
  int deg = static_cast<int>(lg.f.size()) - 1;
  int m = r / deg;
  mpz_class qz(static_cast<long>(q));
  auto red = [&](mpq_class const &x) -> int64_t {
    mpz_class num = x.get_num() % qz;
    mpz_class den = x.get_den() % qz;
    if (num < 0) {
      num += qz;
    }
    if (den < 0) {
      den += qz;
    }
    return MulModQ(num.get_si(), InvModQ(den.get_si(), q), q);
  };
  std::vector<MyMatrix<mpq_class>> Spow{IdentityMat<mpq_class>(r)};
  for (int t = 1; t < deg; t++) {
    Spow.push_back(Spow.back() * lg.S);
  }
  auto orbit_rows = [&](MyVector<mpq_class> const &v) {
    std::vector<std::vector<int64_t>> rows;
    for (int t = 0; t < deg; t++) {
      MyVector<mpq_class> w = Spow[t].transpose() * v;
      std::vector<int64_t> row(r);
      for (int j = 0; j < r; j++) {
        row[j] = red(w(j));
      }
      rows.push_back(row);
    }
    return rows;
  };
  std::vector<MyVector<mpq_class>> chosen;
  std::vector<std::vector<int64_t>> acc;
  auto try_add = [&](MyVector<mpq_class> const &v) -> void {
    std::vector<std::vector<int64_t>> test = acc;
    for (auto &row : orbit_rows(v)) {
      test.push_back(row);
    }
    if (RankModQ(test, q) == static_cast<int>(test.size())) {
      chosen.push_back(v);
      acc = test;
    }
  };
  for (int k = 0; k < r && static_cast<int>(chosen.size()) < m; k++) {
    MyVector<mpq_class> v = ZeroVector<mpq_class>(r);
    v(k) = 1;
    try_add(v);
  }
  uint64_t seed = 1;
  for (int it = 0; it < 1000 && static_cast<int>(chosen.size()) < m; it++) {
    MyVector<mpq_class> v = ZeroVector<mpq_class>(r);
    for (int j = 0; j < r; j++) {
      seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
      v(j) = static_cast<long>((seed >> 33) % 5) - 2;
    }
    try_add(v);
  }
  if (static_cast<int>(chosen.size()) != m) {
    std::cerr << "ABQUAD: no local O_K-basis found\n";
    throw TerminalException{1};
  }
  std::vector<mpq_class> pt = PowerTraces(lg.f, 2 * deg);
  MyMatrix<mpq_class> T(deg, deg);
  for (int a = 0; a < deg; a++) {
    for (int b = 0; b < deg; b++) {
      T(a, b) = pt[a + b];
    }
  }
  MyMatrix<mpq_class> Tinv = Inverse(T);
  std::vector<std::vector<ElemE>> H(m, std::vector<ElemE>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      ElemE val(deg, mpq_class(0));
      for (int t = 0; t < deg; t++) {
        MyVector<mpq_class> w = Spow[t].transpose() * chosen[i];
        mpq_class bt = w.dot(B * chosen[j]);
        for (int s2 = 0; s2 < deg; s2++) {
          val[s2] += bt * Tinv(t, s2);
        }
      }
      H[i][j] = val;
    }
  }
  return H;
}

/*
  The Jordan decomposition of Hecke (jordan_decomposition of Lattices.jl)
  on a Gram matrix over K at p: pieces of rank one, or of rank two at a
  dyadic prime when an off-diagonal entry has a smaller valuation than the
  diagonal ones (normalised to pi^v); consecutive pieces of the same scale
  make a component. The components keep their pieces, along which their
  Gram matrices are block diagonal.
 */
struct LocalJordanComponent {
  int s;
  std::vector<std::vector<std::vector<ElemE>>> pieces;
  int rank() const {
    int r = 0;
    for (auto const &G : pieces) {
      r += static_cast<int>(G.size());
    }
    return r;
  }
};

inline ElemE EAddE(ElemE a, ElemE const &b) {
  for (size_t i = 0; i < a.size(); i++) {
    a[i] += b[i];
  }
  return a;
}

inline ElemE EScal(ElemE a, mpq_class const &c) {
  for (auto &x : a) {
    x *= c;
  }
  return a;
}

inline std::vector<LocalJordanComponent>
LocalJordanDecomposition(std::vector<std::vector<ElemE>> G, LocalPrimeK &P) {
  CyclotomicField const &K = P.K;
  bool even = (P.q == 2);
  int n = static_cast<int>(G.size());
  auto swap_idx = [&](int a, int b) {
    if (a == b) {
      return;
    }
    std::swap(G[a], G[b]);
    for (auto &row : G) {
      std::swap(row[a], row[b]);
    }
  };
  // x_l <- x_l + c x_k
  auto add_mult = [&](int l, int k, ElemE const &c) {
    for (int u = 0; u < n; u++) {
      G[l][u] = EAddE(G[l][u], EMul(c, G[k][u], K));
    }
    for (int u = 0; u < n; u++) {
      G[u][l] = EAddE(G[u][l], EMul(G[u][k], c, K));
    }
  };
  auto scale_idx = [&](int k, ElemE const &c) {
    for (int u = 0; u < n; u++) {
      G[k][u] = EMul(c, G[k][u], K);
    }
    for (int u = 0; u < n; u++) {
      G[u][k] = EMul(G[u][k], c, K);
    }
  };
  std::vector<LocalJordanComponent> res;
  int oldval = 1 << 29;
  int k = 0;
  while (k < n) {
    int mval = 1 << 28, p1 = -1, p2 = -1;
    for (int i = k; i < n; i++) {
      int v = LocalVal(G[i][i], P);
      if (v < mval) {
        mval = v;
        p1 = i;
        p2 = i;
      }
    }
    for (int i = k; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        int v = LocalVal(G[i][j], P);
        if (v < mval) {
          mval = v;
          p1 = i;
          p2 = j;
        }
      }
    }
    if (p1 < 0) {
      std::cerr << "ABQUAD: degenerate Gram matrix\n";
      throw TerminalException{1};
    }
    if (mval != oldval) {
      res.push_back({mval, {}});
      oldval = mval;
    }
    if (even && p1 != p2) {
      swap_idx(p1, k);
      swap_idx(p2 == k ? p1 : p2, k + 1);
      // x_k <- pi^v / T12 x_k
      ElemE T12 = G[k][k + 1];
      int v12 = LocalVal(T12, P);
      ElemE c = EInv(T12, K);
      for (int i = 0; i < v12; i++) {
        c = EMul(c, P.pi, K);
      }
      ElemE pinv = EInv(P.pi, K);
      for (int i = 0; i < -v12; i++) {
        c = EMul(c, pinv, K);
      }
      scale_idx(k, c);
      ElemE T11 = G[k][k], T22 = G[k + 1][k + 1];
      T12 = G[k][k + 1];
      ElemE d = ESub(EMul(T11, T22, K), EMul(T12, T12, K));
      ElemE dinv = EInv(d, K);
      for (int l = k + 2; l < n; l++) {
        ElemE tl = ESub(EMul(T12, G[k + 1][l], K), EMul(T22, G[k][l], K));
        ElemE ul = ESub(EMul(T12, G[k][l], K), EMul(T11, G[k + 1][l], K));
        ElemE c1 = EMul(tl, dinv, K), c2 = EMul(ul, dinv, K);
        add_mult(l, k, c1);
        add_mult(l, k + 1, c2);
      }
      res.back().pieces.push_back(
          {{G[k][k], G[k][k + 1]}, {G[k + 1][k], G[k + 1][k + 1]}});
      k += 2;
    } else {
      if (p1 != p2) {
        add_mult(p1, p2, EZetaPower(0, K));
      }
      swap_idx(p1, k);
      ElemE nrm = G[k][k];
      ElemE ninv = EInv(nrm, K);
      for (int l = k + 1; l < n; l++) {
        ElemE c = EScal(EMul(G[k][l], ninv, K), mpq_class(-1));
        add_mult(l, k, c);
      }
      res.back().pieces.push_back({{G[k][k]}});
      k += 1;
    }
  }
  return res;
}

inline ElemE LocalComponentDet(LocalJordanComponent const &C,
                               CyclotomicField const &K) {
  ElemE det = EZetaPower(0, K);
  for (auto const &G : C.pieces) {
    if (G.size() == 1) {
      det = EMul(det, G[0][0], K);
    } else {
      det = EMul(det, ESub(EMul(G[0][0], G[1][1], K), EMul(G[0][1], G[1][0], K)),
                 K);
    }
  }
  return det;
}

// x pi^{-k}
inline ElemE LocalShift(ElemE x, int k, LocalPrimeK const &P) {
  if (k > 0) {
    ElemE pinv = EInv(P.pi, P.K);
    for (int i = 0; i < k; i++) {
      x = EMul(x, pinv, P.K);
    }
  }
  for (int i = 0; i < -k; i++) {
    x = EMul(x, P.pi, P.K);
  }
  return x;
}

// The square class (1 or -1) in F_q[x] / g of a residue, q odd.
inline int ResidueSquareClass(PolyQ const &u, LocalPrimeK const &P) {
  if (u.empty()) {
    std::cerr << "ABQUAD: the square class of zero\n";
    throw TerminalException{1};
  }
  PolyQ y = PolyPowModQ(u, static_cast<uint64_t>((P.qf - 1) / 2), P.g, P.q);
  PolyTrim(y);
  return (y.size() == 1 && y[0] == 1) ? 1 : -1;
}

// The local factor at an odd prime p from the Gram matrix over K.
inline mpq_class LocalFactorOddGram(std::vector<std::vector<ElemE>> const &H,
                                    LocalPrimeK &P) {
  std::vector<QuadJordanBlock> blocks;
  for (auto const &C : LocalJordanDecomposition(H, P)) {
    int mi = C.rank();
    int chi = 1;
    if ((C.s * mi) % 2 == 0) {
      ElemE u = LocalShift(LocalComponentDet(C, P.K), C.s * mi, P);
      chi = ResidueSquareClass(LocalResidue(u, P), P);
    }
    blocks.push_back({C.s, mi, chi});
  }
  return QuadLocalFactorOdd(blocks, P.qf);
}

/*
  The residue field F = F_q[x] / g: elements as reduced polynomials.
 */
struct ResidueFieldGF {
  int64_t q;
  PolyQ g;
  int f;
};

inline PolyQ GFReduce(PolyQ a, ResidueFieldGF const &F) {
  a = PolyReduce(a, F.q);
  if (a.empty()) {
    return a;
  }
  PolyQ r = PolyDivModQ(a, F.g, F.q).second;
  PolyTrim(r);
  return r;
}

inline PolyQ GFAdd(PolyQ a, PolyQ const &b, ResidueFieldGF const &F) {
  if (a.size() < b.size()) {
    a.resize(b.size(), 0);
  }
  for (size_t i = 0; i < b.size(); i++) {
    a[i] = (a[i] + b[i]) % F.q;
  }
  PolyTrim(a);
  return a;
}

inline PolyQ GFNeg(PolyQ a, ResidueFieldGF const &F) {
  for (auto &x : a) {
    x = (F.q - x) % F.q;
  }
  PolyTrim(a);
  return a;
}

inline PolyQ GFMul(PolyQ const &a, PolyQ const &b, ResidueFieldGF const &F) {
  return GFReduce(PolyMulQ(a, b, F.q), F);
}

inline PolyQ GFInv(PolyQ const &a, ResidueFieldGF const &F) {
  return GFReduce(PolyInverseModQ(a, F.g, F.q), F);
}

// a^(q^k)
inline PolyQ GFFrobenius(PolyQ a, int k, ResidueFieldGF const &F) {
  for (int i = 0; i < k; i++) {
    PolyQ b = {1};
    for (int64_t j = 0; j < F.q; j++) {
      b = GFMul(b, a, F);
    }
    a = b;
  }
  return a;
}

// the trace to F_q
inline int64_t GFTrace(PolyQ const &a, ResidueFieldGF const &F) {
  PolyQ t;
  PolyQ x = a;
  for (int i = 0; i < F.f; i++) {
    t = GFAdd(t, x, F);
    x = GFFrobenius(x, 1, F);
  }
  if (t.size() > 1) {
    std::cerr << "ABQUAD: the trace is not in the prime field\n";
    throw TerminalException{1};
  }
  return t.empty() ? 0 : t[0];
}

using GFMat = std::vector<std::vector<PolyQ>>;

// A basis of {v : v M = 0} for M with nr rows.
inline GFMat GFLeftKernel(GFMat const &M, size_t nr, size_t nc,
                          ResidueFieldGF const &F) {
  // row reduce [M | I]
  GFMat A(nr, std::vector<PolyQ>(nc + nr));
  for (size_t i = 0; i < nr; i++) {
    for (size_t j = 0; j < nc; j++) {
      A[i][j] = M[i][j];
    }
    A[i][nc + i] = {1};
  }
  size_t row = 0;
  for (size_t col = 0; col < nc && row < nr; col++) {
    size_t piv = row;
    while (piv < nr && A[piv][col].empty()) {
      piv++;
    }
    if (piv == nr) {
      continue;
    }
    std::swap(A[piv], A[row]);
    PolyQ inv = GFInv(A[row][col], F);
    for (auto &x : A[row]) {
      x = GFMul(x, inv, F);
    }
    for (size_t i = 0; i < nr; i++) {
      if (i != row && !A[i][col].empty()) {
        PolyQ c = A[i][col];
        for (size_t j = 0; j < nc + nr; j++) {
          A[i][j] = GFAdd(A[i][j], GFNeg(GFMul(c, A[row][j], F), F), F);
        }
      }
    }
    row++;
  }
  GFMat res;
  for (size_t i = row; i < nr; i++) {
    res.push_back(std::vector<PolyQ>(A[i].begin() + nc, A[i].end()));
  }
  return res;
}

inline PolyQ GFBil(std::vector<PolyQ> const &x, GFMat const &Q,
                   std::vector<PolyQ> const &y, ResidueFieldGF const &F) {
  PolyQ res;
  for (size_t i = 0; i < x.size(); i++) {
    if (x[i].empty()) {
      continue;
    }
    for (size_t j = 0; j < y.size(); j++) {
      if (!y[j].empty() && !Q[i][j].empty()) {
        res = GFAdd(res, GFMul(GFMul(x[i], Q[i][j], F), y[j], F), F);
      }
    }
  }
  return res;
}

/*
  A quadratic form over F of characteristic 2 given by an upper triangular
  matrix Q (Q(x) = x Q x^T): the dimension of the quotient by its radical
  (the vectors of the radical of the polar form where Q vanishes) and, for
  that dimension even, the type: +1 when the Arf invariant has trace 0.
 */
inline std::pair<int, int> QuadFormTypeChar2(GFMat const &Q, ResidueFieldGF const &F) {
  size_t n = Q.size();
  if (n == 0) {
    return {0, 1};
  }
  GFMat Bf(n, std::vector<PolyQ>(n));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      Bf[i][j] = GFAdd(Q[i][j], Q[j][i], F);
    }
  }
  GFMat N = GFLeftKernel(Bf, n, n, F);
  bool qzero = true;
  for (auto const &v : N) {
    if (!GFBil(v, Q, v, F).empty()) {
      qzero = false;
    }
  }
  int dim_rad = static_cast<int>(N.size()) - (qzero ? 0 : 1);
  int mi = static_cast<int>(n) - dim_rad;
  if (mi % 2 == 1 || mi == 0) {
    return {mi, 1};
  }
  // a complement of N: standard vectors independent modulo N
  GFMat basis = N;
  GFMat comp;
  for (size_t k = 0; k < n && basis.size() < n; k++) {
    std::vector<PolyQ> e(n);
    e[k] = {1};
    GFMat test = basis;
    test.push_back(e);
    // independent iff the left kernel of test (as rows) is trivial
    if (GFLeftKernel(test, test.size(), n, F).empty()) {
      basis = test;
      comp.push_back(e);
    }
  }
  // symplectic reduction on comp
  PolyQ arf;
  while (!comp.empty()) {
    std::vector<PolyQ> x = comp[0];
    size_t j = 1;
    while (j < comp.size() && GFBil(x, Bf, comp[j], F).empty()) {
      j++;
    }
    if (j == comp.size()) {
      std::cerr << "ABQUAD: degenerate polar form\n";
      throw TerminalException{1};
    }
    PolyQ c = GFInv(GFBil(x, Bf, comp[j], F), F);
    std::vector<PolyQ> y = comp[j];
    for (auto &t : y) {
      t = GFMul(t, c, F);
    }
    // b(x, y) = 1
    arf = GFAdd(arf, GFMul(GFBil(x, Q, x, F), GFBil(y, Q, y, F), F), F);
    GFMat rest;
    for (size_t k = 1; k < comp.size(); k++) {
      if (k == j) {
        continue;
      }
      std::vector<PolyQ> z = comp[k];
      // z <- z + b(z, y) x + b(z, x) y (char 2)
      PolyQ a1 = GFBil(z, Bf, y, F), a2 = GFBil(z, Bf, x, F);
      for (size_t u = 0; u < n; u++) {
        z[u] = GFAdd(z[u], GFAdd(GFMul(a1, x[u], F), GFMul(a2, y[u], F), F), F);
      }
      rest.push_back(z);
    }
    comp = rest;
  }
  return {mi, (GFTrace(arf, F) == 0) ? 1 : -1};
}

/*
  The square class data of an element d of K at a dyadic prime p
  unramified over 2 (Hecke: is_local_square, and the valuation at p of the
  relative discriminant of K(sqrt d) / K): returns 0 for a square, else
  that valuation (0 unramified, 2 or 3 ramified). With d = 2^v u, v odd
  gives 3; otherwise u = w^2 (1 + 2 t) with w lifting the square root of
  the residue of u; t not divisible by p gives 2, and t = 2 t' gives a
  square or the unramified extension as the trace of t' to F_2 is 0 or 1.
 */
inline int DyadicUnramifiedSquareClass(ElemE const &d, LocalPrimeK &P,
                                       ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  int v = LocalVal(d, P);
  if (v % 2 != 0) {
    return 3;
  }
  ElemE u = LocalShift(d, v, P);
  PolyQ ur = LocalResidue(u, P);
  PolyQ wr = GFFrobenius(ur, F.f - 1, F);
  ElemE w(K.deg, mpq_class(0));
  for (size_t i = 0; i < wr.size(); i++) {
    w[i] = mpq_class(static_cast<long>(wr[i]));
  }
  ElemE x = EMul(u, EInv(EMul(w, w, K), K), K);
  ElemE x1 = x;
  x1[0] -= 1;
  int vx = LocalVal(x1, P);
  if (vx == 1) {
    return 2;
  }
  if (vx >= 3) {
    return 0;
  }
  PolyQ t = LocalResidue(LocalShift(x1, 2, P), P);
  return (GFTrace(t, F) == 0) ? 0 : -1;
}

/*
  The local factor at a dyadic prime p unramified over 2, after Cho
  (Hecke: _local_factor_cho), from the Jordan decomposition. For each scale
  s, the residues of the diagonal of the Gram matrix scaled by
  2^{2 (s - S_j)} below s and divided by 2^s define the space of the
  vectors where they vanish, with the quadratic form x A x^T / 2^{s+1} on
  the square roots of their coordinates; its quotient by the radical gives
  the group O_{m_s}(q) of type +/-. The parity type of a component is
  whether its norm has the valuation of its scale.
 */
inline mpq_class LocalFactorCho(std::vector<std::vector<ElemE>> const &H,
                                LocalPrimeK &P) {
  CyclotomicField const &K = P.K;
  ResidueFieldGF F{P.q, P.g, P.f_res};
  int64_t q = P.qf;
  std::vector<LocalJordanComponent> comps = LocalJordanDecomposition(H, P);
  int m = static_cast<int>(H.size());
  size_t nc = comps.size();
  std::vector<int> S(nc), M(nc);
  std::vector<bool> PT(nc);
  for (size_t i = 0; i < nc; i++) {
    S[i] = comps[i].s;
    M[i] = comps[i].rank();
    bool pt = false;
    for (auto const &G : comps[i].pieces) {
      for (size_t a = 0; a < G.size(); a++) {
        if (LocalVal(G[a][a], P) == S[i]) {
          pt = true;
        }
      }
    }
    PT[i] = pt;
  }
  // the flattened Gram matrix and the component of each index
  std::vector<std::vector<ElemE>> AG0(m, std::vector<ElemE>(m, ElemE(K.deg, mpq_class(0))));
  std::vector<size_t> comp_of(m);
  {
    int pos = 0;
    for (size_t i = 0; i < nc; i++) {
      for (auto const &G : comps[i].pieces) {
        for (size_t a = 0; a < G.size(); a++) {
          for (size_t b = 0; b < G.size(); b++) {
            AG0[pos + a][pos + b] = G[a][b];
          }
          comp_of[pos + a] = i;
        }
        pos += static_cast<int>(G.size());
      }
    }
  }
  auto pow2 = [&](int k) -> mpq_class {
    mpq_class r(1);
    for (int i = 0; i < (k >= 0 ? k : -k); i++) {
      r *= 2;
    }
    return (k >= 0) ? r : mpq_class(1 / r);
  };
  // the forms V_s: dimension and type
  std::vector<std::pair<int, int>> V(nc);
  for (size_t is = 0; is < nc; is++) {
    int s = S[is];
    std::vector<std::vector<ElemE>> AG = AG0;
    for (int a = 0; a < m; a++) {
      for (int b = 0; b < m; b++) {
        int Sj = S[comp_of[a]];
        if (Sj < s) {
          AG[a][b] = EScal(AG[a][b], pow2(2 * (s - Sj)));
        }
      }
    }
    GFMat Dcol(m, std::vector<PolyQ>(1));
    for (int a = 0; a < m; a++) {
      Dcol[a][0] = LocalResidue(EScal(AG[a][a], pow2(-s)), P);
    }
    GFMat Bk = GFLeftKernel(Dcol, m, 1, F);
    size_t nb = Bk.size();
    std::vector<std::vector<ElemE>> BK(nb, std::vector<ElemE>(m, ElemE(K.deg, mpq_class(0))));
    for (size_t i = 0; i < nb; i++) {
      for (int a = 0; a < m; a++) {
        PolyQ r = GFFrobenius(Bk[i][a], F.f - 1, F);
        for (size_t t = 0; t < r.size(); t++) {
          BK[i][a][t] = mpq_class(static_cast<long>(r[t]));
        }
      }
    }
    GFMat Q(nb, std::vector<PolyQ>(nb));
    for (size_t i = 0; i < nb; i++) {
      for (size_t j = i; j < nb; j++) {
        ElemE val(K.deg, mpq_class(0));
        for (int a = 0; a < m; a++) {
          if (EIsZero(BK[i][a])) {
            continue;
          }
          for (int b = 0; b < m; b++) {
            if (EIsZero(BK[j][b]) || EIsZero(AG[a][b])) {
              continue;
            }
            val = EAddE(val, EMul(EMul(BK[i][a], AG[a][b], K), BK[j][b], K));
          }
        }
        val = EScal(val, pow2(-(s + 1)));
        if (i < j) {
          val = EScal(val, mpq_class(2));
        }
        Q[i][j] = LocalResidue(val, P);
      }
    }
    V[is] = QuadFormTypeChar2(Q, F);
  }
  auto in_S = [&](int x) -> int {
    for (size_t i = 0; i < nc; i++) {
      if (S[i] == x) {
        return static_cast<int>(i);
      }
    }
    return -1;
  };
  int n_alpha = 0;
  for (size_t i = 0; i < nc; i++) {
    int below = in_S(S[i] - 1), above = in_S(S[i] + 1);
    if (V[i].first != 0 && (below < 0 || !PT[below]) && (above < 0 || !PT[above])) {
      n_alpha++;
    }
  }
  int n_beta = 0;
  for (size_t i = 0; i < nc; i++) {
    if (!PT[i]) {
      continue;
    }
    int idx = in_S(S[i] + 2);
    if (idx < 0 || !PT[idx]) {
      n_beta++;
    }
  }
  mpq_class res = pow2(n_alpha + n_beta);
  int64_t N = 0;
  for (size_t i = 0; i < nc; i++) {
    if (PT[i]) {
      N++;
    } else {
      N += M[i];
    }
    if (i + 1 < nc && PT[i] && PT[i + 1] && S[i + 1] == S[i] + 1) {
      N--;
    }
  }
  for (size_t i = 0; i < nc; i++) {
    int64_t after = 0;
    for (size_t j = i + 1; j < nc; j++) {
      after += M[j];
    }
    N += static_cast<int64_t>(S[i]) * M[i] * (M[i] + 1) / 2;
    N += static_cast<int64_t>(S[i]) * M[i] * after;
  }
  for (auto const &v : V) {
    if (v.first == 0) {
      continue;
    }
    res *= HeckeGroupOrderO(v.first, v.second, q) / 2;
  }
  auto qpow = [&](int64_t k) -> mpq_class {
    mpq_class r(1);
    for (int64_t i = 0; i < (k >= 0 ? k : -k); i++) {
      r *= mpq_class(static_cast<long>(q));
    }
    return (k >= 0) ? r : mpq_class(1 / r);
  };
  mpq_class beta = qpow(N) * res / 2;
  // twice the exponent
  int64_t exp2 = 0;
  for (size_t i = 0; i < nc; i++) {
    exp2 += static_cast<int64_t>(m + 1) * S[i] * M[i];
  }
  mpq_class Hord;
  if (m % 2 == 1) {
    exp2 += m + 1;
    Hord = HeckeGroupOrderO(m, 1, q);
  } else {
    exp2 += 2 * m;
    ElemE d = EZetaPower(0, K);
    for (auto const &C : comps) {
      d = EMul(d, LocalComponentDet(C, K), K);
    }
    if ((m * (m - 1) / 2) % 2 == 1) {
      d = EScal(d, mpq_class(-1));
    }
    int sc = DyadicUnramifiedSquareClass(d, P, F);
    if (sc == 0) {
      Hord = HeckeGroupOrderO(m, 1, q);
    } else if (sc == -1) {
      Hord = HeckeGroupOrderO(m, -1, q);
    } else {
      Hord = HeckeGroupOrderO(m - 1, 1, q);
      exp2 += static_cast<int64_t>(sc) * (1 - m);
    }
  }
  if (exp2 % 2 != 0) {
    std::cerr << "ABQUAD: the exponent of Cho's factor is not integral\n";
    throw TerminalException{1};
  }
  return qpow(exp2 / 2) * Hord / 2 / beta;
}

/*
  The dyadic primes ramified over 2 (Hecke: local_factor of MassQuad.jl,
  with the quadratic defects and the Hilbert symbols of NumField/Hilbert.jl
  by Kirschmer, and his local isometry test of Quad/Legacy.jl). After a
  rescaling by a power of pi: the formula of a maximal lattice when L_p is
  maximal, that of a unimodular lattice when L_p is modular, and otherwise
  the chain L = L_0 < ... < L_k of integral overlattices of index p, L_k
  maximal, with the factor of the maximal lattice times the product over i
  of the number of sublattices of index p of L_{i+1} isometric to L_i
  divided by the number of superlattices of index p of L_i isometric to
  L_{i+1}. The lattices are kept by their Gram matrices over K on local
  bases; the factors do not depend on the rescaling.
 */
using KMat = std::vector<std::vector<ElemE>>;
constexpr int QD_INF = 1 << 28;

inline ElemE LiftResidue(PolyQ const &r, CyclotomicField const &K) {
  ElemE x(K.deg, mpq_class(0));
  for (size_t i = 0; i < r.size(); i++) {
    x[i] = mpq_class(static_cast<long>(r[i]));
  }
  return x;
}

inline bool EIsOne(ElemE const &a) {
  ElemE b = a;
  b[0] -= 1;
  return EIsZero(b);
}

inline ElemE EOne(CyclotomicField const &K) { return EZetaPower(0, K); }

// (the defect of the unit a, a in the same square class with a - 1 of that
// valuation), at a dyadic prime
inline std::pair<int, ElemE> QuadDefectUnit(ElemE a, LocalPrimeK &P,
                                            ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  auto sqrt_res = [&](ElemE const &x) -> ElemE {
    return LiftResidue(GFFrobenius(LocalResidue(x, P), F.f - 1, F), K);
  };
  ElemE s = sqrt_res(a);
  a = EMul(a, EInv(EMul(s, s, K), K), K);
  auto val1 = [&](ElemE const &x) -> int {
    ElemE y = x;
    y[0] -= 1;
    return EIsZero(y) ? QD_INF : LocalVal(y, P);
  };
  int w = val1(a);
  int ee = 2 * P.e;
  while (w < ee && w % 2 == 0) {
    ElemE am1 = a;
    am1[0] -= 1;
    ElemE t = sqrt_res(LocalShift(am1, w, P));
    ElemE c = EOne(K);
    ElemE pw = EOne(K);
    for (int i = 0; i < w / 2; i++) {
      pw = EMul(pw, P.pi, K);
    }
    c = EAddE(c, EMul(t, pw, K));
    a = EMul(a, EInv(EMul(c, c, K), K), K);
    w = val1(a);
  }
  if (w > ee) {
    return {QD_INF, EOne(K)};
  }
  if (w == ee) {
    ElemE am1 = a;
    am1[0] -= 1;
    PolyQ c = LocalResidue(EScal(am1, mpq_class(1, 4)), P);
    if (GFTrace(c, F) == 0) {
      return {QD_INF, EOne(K)};
    }
  }
  return {w, a};
}

inline int QuadDefect(ElemE const &a, LocalPrimeK &P, ResidueFieldGF const &F) {
  if (EIsZero(a)) {
    return QD_INF;
  }
  int v = LocalVal(a, P);
  if (v % 2 != 0) {
    return v;
  }
  int w = QuadDefectUnit(LocalShift(a, v, P), P, F).first;
  return (w == QD_INF) ? QD_INF : w + v;
}

// The Hilbert symbol (a, b)_p at the dyadic prime p (Kirschmer).
inline int HilbertSymbolDyadic(ElemE a, ElemE b, LocalPrimeK &P,
                               ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  int v = LocalVal(a, P);
  int w = LocalVal(b, P);
  if (v % 2 != 0) {
    if (w % 2 == 0) {
      std::swap(a, b);
      std::swap(v, w);
    } else {
      a = EScal(EMul(a, b, K), mpq_class(-1));
      v = v + w;
    }
  }
  a = LocalShift(a, v, P);
  auto [vd, a2] = QuadDefectUnit(a, P, F);
  a = a2;
  int wf = (w >= 0) ? w / 2 : -((-w + 1) / 2);
  b = LocalShift(b, 2 * wf, P);
  w = w - 2 * wf;
  if (vd != QD_INF && vd % 2 != 0 && w == 0) {
    ElemE oma = EScal(a, mpq_class(-1));
    oma[0] += 1;
    b = LocalShift(EMul(oma, b, K), vd - 1, P);
    w = 1;
  }
  while (vd != QD_INF && vd % 2 != 0) {
    ElemE t = b;
    for (int i = 0; i < vd - 1; i++) {
      t = EMul(t, P.pi, K);
    }
    ElemE am1 = a;
    am1[0] -= 1;
    ElemE sq = LiftResidue(
        GFFrobenius(LocalResidue(EMul(am1, EInv(t, K), K), P), F.f - 1, F), K);
    ElemE fac = EScal(EMul(EMul(sq, sq, K), t, K), mpq_class(-1));
    fac[0] += 1;
    a = EMul(a, fac, K);
    auto pr = QuadDefectUnit(a, P, F);
    vd = pr.first;
    a = pr.second;
  }
  return (vd != QD_INF && w % 2 != 0) ? -1 : 1;
}

// A diagonal form isometric to the nondegenerate G.
inline std::vector<ElemE> DiagonalizeK(KMat F, CyclotomicField const &K) {
  int n = static_cast<int>(F.size());
  std::vector<ElemE> res;
  for (int i = 0; i < n; i++) {
    if (EIsZero(F[i][i])) {
      int j = -1;
      for (int t = i + 1; t < n; t++) {
        if (!EIsZero(F[t][t])) {
          j = t;
          break;
        }
      }
      if (j >= 0) {
        std::swap(F[i], F[j]);
        for (auto &row : F) {
          std::swap(row[i], row[j]);
        }
      } else {
        for (int t = i + 1; t < n; t++) {
          if (!EIsZero(F[i][t])) {
            j = t;
            break;
          }
        }
        if (j < 0) {
          std::cerr << "ABQUAD: a degenerate form\n";
          throw TerminalException{1};
        }
        // x_i <- x_i + c x_j, c = 1 / (2 F_ij)
        ElemE c = EInv(EScal(F[i][j], mpq_class(2)), K);
        for (int u = 0; u < n; u++) {
          F[i][u] = EAddE(F[i][u], EMul(c, F[j][u], K));
        }
        for (int u = 0; u < n; u++) {
          F[u][i] = EAddE(F[u][i], EMul(F[u][j], c, K));
        }
      }
    }
    ElemE dinv = EInv(F[i][i], K);
    for (int t = i + 1; t < n; t++) {
      ElemE c = EScal(EMul(F[t][i], dinv, K), mpq_class(-1));
      for (int u = 0; u < n; u++) {
        F[t][u] = EAddE(F[t][u], EMul(c, F[i][u], K));
      }
      for (int u = 0; u < n; u++) {
        F[u][t] = EAddE(F[u][t], EMul(F[u][i], c, K));
      }
    }
    res.push_back(F[i][i]);
  }
  return res;
}

inline int HasseInvariantDyadic(std::vector<ElemE> const &D, LocalPrimeK &P,
                                ResidueFieldGF const &F) {
  int h = 1;
  for (size_t i = 0; i < D.size(); i++) {
    for (size_t j = i + 1; j < D.size(); j++) {
      h *= HilbertSymbolDyadic(D[i], D[j], P, F);
    }
  }
  return h;
}

inline ElemE DetK(KMat const &G, CyclotomicField const &K) {
  ElemE d = EOne(K);
  for (auto const &x : DiagonalizeK(G, K)) {
    d = EMul(d, x, K);
  }
  return d;
}

// The Witt invariant of the space of G (Hecke: _hasse_witt).
inline int WittInvariantDyadic(KMat const &G, LocalPrimeK &P,
                               ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  std::vector<ElemE> D = DiagonalizeK(G, K);
  int h = HasseInvariantDyadic(D, P, F);
  ElemE d = EOne(K);
  for (auto const &x : D) {
    d = EMul(d, x, K);
  }
  int n = static_cast<int>(D.size()) % 8;
  ElemE c = EOne(K);
  if (n == 3 || n == 4) {
    c = EScal(d, mpq_class(-1));
  } else if (n == 5 || n == 6) {
    c = EScal(EOne(K), mpq_class(-1));
  } else if (n == 7 || n == 0) {
    c = d;
  }
  return h * HilbertSymbolDyadic(EScal(EOne(K), mpq_class(-1)), c, P, F);
}

inline KMat ScaleKMat(KMat G, ElemE const &c, CyclotomicField const &K) {
  for (auto &row : G) {
    for (auto &x : row) {
      x = EMul(x, c, K);
    }
  }
  return G;
}

inline ElemE PiPower(int k, LocalPrimeK const &P) {
  return LocalShift(EOne(P.K), -k, P);
}

inline int NormValK(KMat const &G, LocalPrimeK &P) {
  int v = QD_INF;
  for (size_t i = 0; i < G.size(); i++) {
    v = std::min(v, LocalVal(G[i][i], P));
    for (size_t j = i + 1; j < G.size(); j++) {
      int w = LocalVal(G[i][j], P);
      if (w != QD_INF) {
        v = std::min(v, w + P.e);
      }
    }
  }
  return v;
}

inline int ScaleValK(KMat const &G, LocalPrimeK &P) {
  int v = QD_INF;
  for (auto const &row : G) {
    for (auto const &x : row) {
      v = std::min(v, LocalVal(x, P));
    }
  }
  return v;
}

inline mpq_class QPowRat(int64_t q, int64_t k) {
  mpq_class r(1);
  for (int64_t i = 0; i < (k >= 0 ? k : -k); i++) {
    r *= mpq_class(static_cast<long>(q));
  }
  return (k >= 0) ? r : mpq_class(1 / r);
}

// Hecke: _local_factor_maximal, on the space of G
inline mpq_class LocalFactorMaximalDyadic(KMat const &G, LocalPrimeK &P,
                                          ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  int m = static_cast<int>(G.size());
  int r = m / 2;
  if (m == 1) {
    return mpq_class(1);
  }
  int w = WittInvariantDyadic(G, P, F);
  ElemE d = DetK(G, K);
  if ((m * (m - 1) / 2) % 2 == 1) {
    d = EScal(d, mpq_class(-1));
  }
  int64_t q = P.qf;
  int vd = LocalVal(d, P);
  if (m % 2 == 1) {
    if (vd % 2 != 0) {
      return (QPowRat(q, r) + w) / 2;
    } else if (w == -1) {
      return (QPowRat(q, m - 1) - 1) / mpq_class(static_cast<long>(q + 1)) / 2;
    }
  } else {
    int qd = QuadDefect(d, P, F);
    if (vd % 2 != 0) {
      return mpq_class(1, 2);
    } else if (qd == QD_INF) {
      if (w == -1) {
        return (QPowRat(q, r - 1) - 1) * (QPowRat(q, r) - 1) /
               mpq_class(static_cast<long>(q + 1)) / 2;
      }
    } else if (qd % 2 != 0) {
      return mpq_class(1, 2);
    } else if (w == -1) {
      return (QPowRat(q, r - 1) + 1) * (QPowRat(q, r) + 1) /
             mpq_class(static_cast<long>(q + 1)) / 2;
    }
  }
  return mpq_class(1);
}

/*
  The genus symbol of Kirschmer at a dyadic prime: for the Jordan
  components G_i of scales s_i, the norm generators a_i, u_i = v(a_i), the
  weights w_i and the f_i.
 */
struct QuadKirschmerSymbol {
  std::vector<int> dims, sL, wL, uL, fL;
  std::vector<ElemE> aL;
  std::vector<std::vector<ElemE>> diag;
  std::vector<ElemE> dets;
};

inline QuadKirschmerSymbol KirschmerSymbolDyadic(KMat const &G, LocalPrimeK &P,
                                             ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  std::vector<LocalJordanComponent> comps = LocalJordanDecomposition(G, P);
  QuadKirschmerSymbol S;
  int t = static_cast<int>(comps.size());
  int e = P.e;
  // the diagonal entries of each component
  std::vector<std::vector<ElemE>> Dc(t);
  for (int i = 0; i < t; i++) {
    S.dims.push_back(comps[i].rank());
    S.sL.push_back(comps[i].s);
    for (auto const &Gp : comps[i].pieces) {
      for (size_t a = 0; a < Gp.size(); a++) {
        Dc[i].push_back(Gp[a][a]);
      }
    }
    // a diagonalisation of each piece, for the Hasse invariants
    std::vector<ElemE> dg;
    for (auto const &Gp : comps[i].pieces) {
      for (auto const &x : DiagonalizeK(Gp, K)) {
        dg.push_back(x);
      }
    }
    S.diag.push_back(dg);
    S.dets.push_back(LocalComponentDet(comps[i], K));
  }
  for (int i = 0; i < t; i++) {
    std::vector<ElemE> D;
    for (int j = 0; j < t; j++) {
      for (auto const &x : Dc[j]) {
        D.push_back(j < i ? EMul(x, PiPower(2 * (S.sL[i] - S.sL[j]), P), K)
                          : x);
      }
    }
    int mv = QD_INF, pos = -1;
    for (size_t k = 0; k < D.size(); k++) {
      int v = EIsZero(D[k]) ? QD_INF : LocalVal(D[k], P);
      if (v < mv) {
        mv = v;
        pos = static_cast<int>(k);
      }
    }
    ElemE a = (e + S.sL[i] <= mv) ? PiPower(e + S.sL[i], P) : D[pos];
    S.aL.push_back(a);
    int u = LocalVal(a, P);
    S.uL.push_back(u);
    int w = e + S.sL[i];
    ElemE ainv = EInv(a, K);
    for (auto const &x : D) {
      int qd = QuadDefect(EMul(x, ainv, K), P, F);
      if (qd != QD_INF) {
        w = std::min(w, u + qd);
      }
    }
    S.wL.push_back(w);
  }
  for (int k = 0; k + 1 < t; k++) {
    int ex = S.uL[k] + S.uL[k + 1];
    int val;
    if (ex % 2 != 0) {
      val = ex;
    } else {
      val = std::min({QuadDefect(EMul(S.aL[k], S.aL[k + 1], K), P, F),
                      S.uL[k] + S.wL[k + 1], S.uL[k + 1] + S.wL[k],
                      e + ex / 2 + S.sL[k]});
    }
    S.fL.push_back(val - 2 * S.sL[k]);
  }
  return S;
}

// Kirschmer's test (Hecke: is_locally_isometric_kirschmer) for two
// lattices of the same quadratic space at a dyadic prime.
inline bool IsLocallyIsometricDyadic(QuadKirschmerSymbol const &S1,
                                     QuadKirschmerSymbol const &S2, LocalPrimeK &P,
                                     ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  if (S1.dims != S2.dims || S1.sL != S2.sL || S1.wL != S2.wL ||
      S1.uL != S2.uL) {
    return false;
  }
  int e = P.e;
  for (size_t k = 0; k < S1.aL.size(); k++) {
    ElemE b = EMul(S1.aL[k], EInv(S2.aL[k], K), K);
    int qd = QuadDefect(b, P, F);
    if (qd < S1.wL[k] - S2.uL[k]) {
      return false;
    }
  }
  auto comp_hasse = [&](size_t i, ElemE const &a) -> bool {
    std::vector<ElemE> D1, D2;
    for (size_t j = 0; j <= i; j++) {
      for (auto const &x : S1.diag[j]) {
        D1.push_back(x);
      }
      for (auto const &x : S2.diag[j]) {
        D2.push_back(x);
      }
    }
    D2.push_back(a);
    ElemE prod = EOne(K);
    for (auto const &x : D1) {
      prod = EMul(prod, x, K);
    }
    for (auto const &x : D2) {
      prod = EMul(prod, x, K);
    }
    D1.push_back(prod);
    return HasseInvariantDyadic(D1, P, F) == HasseInvariantDyadic(D2, P, F);
  };
  ElemE d1 = EOne(K), d2 = EOne(K);
  for (size_t i = 0; i < S1.fL.size(); i++) {
    d1 = EMul(d1, S1.dets[i], K);
    d2 = EMul(d2, S2.dets[i], K);
    ElemE dq = EMul(d1, EInv(d2, K), K);
    if (LocalVal(dq, P) != 0) {
      return false;
    }
    if (QuadDefect(dq, P, F) < S1.fL[i]) {
      return false;
    }
    if (S1.fL[i] > 2 * e + S1.uL[i + 1] - S2.wL[i + 1] &&
        !comp_hasse(i, S2.aL[i + 1])) {
      return false;
    }
    if (S1.fL[i] > 2 * e + S1.uL[i] - S2.wL[i] && !comp_hasse(i, S2.aL[i])) {
      return false;
    }
  }
  return true;
}

// The vectors of k^m with first nonzero coordinate 1.
inline std::vector<std::vector<PolyQ>> ProjectivePointsGF(int m,
                                                          ResidueFieldGF const &F) {
  std::vector<PolyQ> elts;
  int64_t qf = 1;
  for (int i = 0; i < F.f; i++) {
    qf *= F.q;
  }
  for (int64_t u = 0; u < qf; u++) {
    PolyQ a;
    int64_t t = u;
    for (int i = 0; i < F.f; i++) {
      a.push_back(t % F.q);
      t /= F.q;
    }
    PolyTrim(a);
    elts.push_back(a);
  }
  std::vector<std::vector<PolyQ>> res;
  for (int lead = 0; lead < m; lead++) {
    int nfree = m - lead - 1;
    int64_t tot = 1;
    for (int i = 0; i < nfree; i++) {
      tot *= qf;
    }
    for (int64_t u = 0; u < tot; u++) {
      std::vector<PolyQ> v(m);
      v[lead] = {1};
      int64_t t = u;
      for (int i = lead + 1; i < m; i++) {
        v[i] = elts[t % qf];
        t /= qf;
      }
      res.push_back(v);
    }
  }
  return res;
}

inline KMat TransformKMat(KMat const &T, KMat const &G, CyclotomicField const &K) {
  size_t n = G.size();
  KMat TG(n, std::vector<ElemE>(n, ElemE(K.deg, mpq_class(0))));
  for (size_t i = 0; i < n; i++) {
    for (size_t k = 0; k < n; k++) {
      if (EIsZero(T[i][k])) {
        continue;
      }
      for (size_t j = 0; j < n; j++) {
        TG[i][j] = EAddE(TG[i][j], EMul(T[i][k], G[k][j], K));
      }
    }
  }
  KMat res(n, std::vector<ElemE>(n, ElemE(K.deg, mpq_class(0))));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < n; k++) {
        if (!EIsZero(T[j][k])) {
          res[i][j] = EAddE(res[i][j], EMul(TG[i][k], T[j][k], K));
        }
      }
    }
  }
  return res;
}

// L + pi^{-1} y with y = sum c_i b_i
inline KMat SuperlatticeGram(KMat const &G, std::vector<PolyQ> const &c,
                             LocalPrimeK const &P) {
  CyclotomicField const &K = P.K;
  size_t n = G.size();
  KMat T(n, std::vector<ElemE>(n, ElemE(K.deg, mpq_class(0))));
  size_t lead = 0;
  while (c[lead].empty()) {
    lead++;
  }
  ElemE pinv = EInv(P.pi, K);
  for (size_t i = 0; i < n; i++) {
    if (i == lead) {
      for (size_t j = 0; j < n; j++) {
        T[i][j] = EMul(LiftResidue(c[j], K), pinv, K);
      }
    } else {
      T[i][i] = EOne(K);
    }
  }
  return TransformKMat(T, G, K);
}

// {x in L : sum lambda_i x_i = 0 mod p}
inline KMat SublatticeGram(KMat const &G, std::vector<PolyQ> const &lambda,
                           LocalPrimeK const &P) {
  CyclotomicField const &K = P.K;
  size_t n = G.size();
  KMat T(n, std::vector<ElemE>(n, ElemE(K.deg, mpq_class(0))));
  size_t lead = 0;
  while (lambda[lead].empty()) {
    lead++;
  }
  for (size_t i = 0; i < n; i++) {
    if (i == lead) {
      T[i][i] = P.pi;
    } else {
      T[i][i] = EOne(K);
      T[i][lead] = EScal(LiftResidue(lambda[i], K), mpq_class(-1));
    }
  }
  return TransformKMat(T, G, K);
}

// An integral overlattice of index p of G (of integral norm), if any.
inline std::optional<KMat> IntegralOverlattice(KMat const &G, LocalPrimeK &P,
                                               ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  size_t n = G.size();
  for (auto const &c : ProjectivePointsGF(static_cast<int>(n), F)) {
    std::vector<ElemE> y(n, ElemE(K.deg, mpq_class(0)));
    // y G: the products B(y, b_j)
    std::vector<ElemE> cl;
    for (auto const &x : c) {
      cl.push_back(LiftResidue(x, K));
    }
    bool ok = true;
    ElemE Qy(K.deg, mpq_class(0));
    for (size_t j = 0; j < n && ok; j++) {
      ElemE bj(K.deg, mpq_class(0));
      for (size_t i = 0; i < n; i++) {
        if (!EIsZero(cl[i])) {
          bj = EAddE(bj, EMul(cl[i], G[i][j], K));
        }
      }
      if (!EIsZero(bj) && LocalVal(bj, P) < 1 - P.e) {
        ok = false;
      }
      if (!EIsZero(cl[j])) {
        Qy = EAddE(Qy, EMul(bj, cl[j], K));
      }
    }
    if (!ok) {
      continue;
    }
    if (!EIsZero(Qy) && LocalVal(Qy, P) < 2) {
      continue;
    }
    return SuperlatticeGram(G, c, P);
  }
  return {};
}

// Hecke: _local_factor_unimodular, G unimodular of integral norm
inline mpq_class LocalFactorUnimodularDyadic(KMat const &G, LocalPrimeK &P,
                                             ResidueFieldGF const &F) {
  CyclotomicField const &K = P.K;
  QuadKirschmerSymbol S = KirschmerSymbolDyadic(G, P, F);
  if (S.sL.size() != 1 || S.sL[0] != 0) {
    std::cerr << "ABQUAD: the lattice is not unimodular\n";
    throw TerminalException{1};
  }
  int d = S.dims[0];
  int b = S.wL[0];
  int a = S.uL[0];
  int64_t q = P.qf;
  int e = P.e;
  auto qp = [&](mpq_class const &x) -> mpq_class {
    if (x.get_den() != 1) {
      std::cerr << "ABQUAD: a non integral exponent\n";
      throw TerminalException{1};
    }
    return QPowRat(q, x.get_num().get_si());
  };
  auto eps = [&](int c) -> int { return (c == QD_INF) ? 1 : -1; };
  mpq_class lf;
  if (d % 2 == 0) {
    int r = d / 2;
    ElemE disc = S.dets[0];
    if (r % 2 == 1) {
      disc = EScal(disc, mpq_class(-1));
    }
    int c = QuadDefect(disc, P, F);
    mpq_class rh = mpq_class(r) - mpq_class(1, 2);
    if (d == 2) {
      if (a < b && b == 2 && c == 2 * e) {
        lf = qp((e - a - 1) / 2) * (mpq_class(static_cast<long>(q)) - eps(c));
      } else if (b == e && a + e + 1 <= c && c < 2 * e) {
        lf = 2 * qp((c - e - a) / 2);
      } else {
        lf = 1;
      }
    } else if ((a + b) % 2 == 0) {
      if (e == a) {
        lf = 1;
      } else if (c >= 2 * e) {
        lf = qp(mpq_class(e - a) * rh - r) * (QPowRat(q, r) - eps(c));
      } else if (a + e + 1 <= c) {
        lf = 2 * qp(mpq_class(c - e - a - 1) * rh);
      } else {
        lf = qp(mpq_class((c - e - a - 1) * (r - 1)));
      }
    } else {
      if (c == a + b) {
        lf = 1;
      } else if (c >= 2 * e) {
        int hs = (c == QD_INF || a % 2 == 0) ? 1 : -1;
        int cc;
        if (WittInvariantDyadic(G, P, F) == hs) {
          cc = c;
        } else {
          cc = (c == QD_INF) ? 2 * e : QD_INF;
        }
        if (e == b && cc == QD_INF) {
          int ee = eps(c);
          lf = qp(mpq_class(e - a - 1) * rh) * (QPowRat(q, r) - ee) *
               (QPowRat(q, r - 1) + ee) / 2;
        } else if (e == b) {
          lf = qp(mpq_class(e - a - 1) * rh) * mpq_class(static_cast<long>(q + 1));
        } else if (cc == QD_INF && c == cc) {
          int ee = eps(c);
          lf = qp(mpq_class(2 * e - b - a - 3) * rh + r) * (QPowRat(q, r) - ee) *
               (QPowRat(q, 2 * r - 2) - ee) / 2;
        } else if (cc == QD_INF) {
          // as in Hecke: (2 e - b - a - 3 (r - 1/2) + r)
          lf = qp(mpq_class(2 * e - b - a) - 3 * rh + r) *
               mpq_class(static_cast<long>(q + 1));
        } else if (c == cc) {
          int ee = eps(c);
          lf = qp(mpq_class(2 * e - b - a - 3) * rh + r) *
               (QPowRat(q, r - 1) + 1) * (QPowRat(q, r) - ee) *
               (QPowRat(q, r - 1) + ee) / 2;
        } else {
          lf = qp(mpq_class(2 * e - b - a - 3) * rh + r) *
               (QPowRat(q, r - 1) + 1) * mpq_class(static_cast<long>(q + 1));
        }
      } else if (b == e) {
        lf = 2 * qp(mpq_class(c - e - a) * rh);
      } else {
        lf = (QPowRat(q, 2 * r - 2) - 1) * qp(mpq_class(c - a - b - 2) * rh + 1);
      }
    }
  } else {
    if (a != 0) {
      std::cerr << "ABQUAD: odd unimodular with a != 0\n";
      throw TerminalException{1};
    }
    int r = (d - 3) / 2;
    int w = WittInvariantDyadic(G, P, F);
    if (e == b) {
      lf = 1;
    } else if (e % 2 == 1) {
      lf = (QPowRat(q, r + 1) - w) * QPowRat(q, (r + 1) * (e - b - 1));
    } else {
      lf = (w == 1 ? mpq_class((QPowRat(q, 2 * r + 2) - 1) / 2)
                   : mpq_class(static_cast<long>(q + 1))) *
           QPowRat(q, (r + 1) * (e - b - 1));
    }
  }
  ElemE two = EScal(EOne(K), mpq_class(2));
  return LocalFactorMaximalDyadic(ScaleKMat(G, two, K), P, F) * lf;
}

inline mpq_class LocalFactorDyadicRamified(KMat const &H, LocalPrimeK &P) {
  CyclotomicField const &K = P.K;
  ResidueFieldGF F{P.q, P.g, P.f_res};
  int m = static_cast<int>(H.size());
  if (m == 1) {
    return mpq_class(1);
  }
  KMat L = ScaleKMat(H, PiPower(-NormValK(H, P), P), K);
  std::optional<KMat> over = IntegralOverlattice(L, P, F);
  if (!over) {
    return LocalFactorMaximalDyadic(L, P, F);
  }
  int sv = ScaleValK(H, P);
  if (sv * m == LocalVal(DetK(H, K), P)) {
    return LocalFactorUnimodularDyadic(ScaleKMat(H, PiPower(-sv, P), K), P, F);
  }
  std::vector<KMat> chain{L};
  while (over) {
    chain.push_back(*over);
    over = IntegralOverlattice(chain.back(), P, F);
  }
  mpq_class f = LocalFactorMaximalDyadic(L, P, F);
  for (size_t i = 0; i + 1 < chain.size(); i++) {
    QuadKirschmerSymbol Si = KirschmerSymbolDyadic(chain[i], P, F);
    QuadKirschmerSymbol Sn = KirschmerSymbolDyadic(chain[i + 1], P, F);
    int64_t n_sub = 0, n_super = 0;
    for (auto const &lam : ProjectivePointsGF(m, F)) {
      KMat M = SublatticeGram(chain[i + 1], lam, P);
      if (IsLocallyIsometricDyadic(Si, KirschmerSymbolDyadic(M, P, F), P, F)) {
        n_sub++;
      }
    }
    for (auto const &c : ProjectivePointsGF(m, F)) {
      KMat M = SuperlatticeGram(chain[i], c, P);
      if (IsLocallyIsometricDyadic(Sn, KirschmerSymbolDyadic(M, P, F), P, F)) {
        n_super++;
      }
    }
    if (n_super == 0) {
      std::cerr << "ABQUAD: no superlattice isometric to the next in the chain\n";
      throw TerminalException{1};
    }
    f *= mpq_class(static_cast<long>(n_sub)) / mpq_class(static_cast<long>(n_super));
  }
  return f;
}

/*
  The product of the local factors at the primes of K over q, or nothing
  when not computed: the odd primes, and the dyadic primes unramified over
  2 (Cho) and those ramified over 2.
 */
inline std::optional<mpq_class>
QuadLocalFactorAtPrime(MyMatrix<mpq_class> const &B,
                       std::vector<MyMatrix<mpq_class>> const &Omega,
                       int64_t q) {
  std::optional<LocalGenerator> lg = FindLocalGenerator(B, Omega, q);
  if (!lg) {
    return {};
  }
  std::vector<std::pair<PolyQ, int>> fac =
      FactorModQWithMultiplicity(PolyQ(lg->f.begin(), lg->f.end()), q);
  int m = B.rows() / (static_cast<int>(lg->f.size()) - 1);
  if (m == 1) {
    return mpq_class(1);
  }
  std::vector<std::vector<ElemE>> H = LocalOKGram(B, *lg, q);
  mpq_class res(1);
  for (auto const &[g, e] : fac) {
    LocalPrimeK P = MakeLocalPrimeK(*lg, g, e, q);
    if (q != 2) {
      res *= LocalFactorOddGram(H, P);
    } else if (e == 1) {
      res *= LocalFactorCho(H, P);
    } else {
      res *= LocalFactorDyadicRamified(H, P);
    }
  }
  return res;
}

/*
  The local factors: q -> the factor, or nothing for q not computed. The
  odd primes where phi is unimodular up to scaling at every prime of K
  over q are left out (their factor is 1); 2 is always listed.
 */
inline std::map<int64_t, std::optional<mpq_class>>
QuadraticLocalFactors(MyMatrix<mpq_class> const &B,
                      std::vector<MyMatrix<mpq_class>> const &Omega,
                      AbelianCMField const &F) {
  int deg = static_cast<int>(Omega.size());
  int m = B.rows() / deg;
  DirichletGroup G = ComputeDirichletGroup(F.n);
  std::vector<std::vector<int64_t>> XK = FieldCharacters(G, F.H, true);
  // a rescaling with phi integral: b(X, X) in Z gives phi in the inverse
  // different
  mpz_class den(1);
  for (int i = 0; i < B.rows(); i++) {
    for (int j = 0; j < B.cols(); j++) {
      den = lcm(den, B(i, j).get_den());
    }
  }
  mpz_class dK(1);
  std::set<int64_t> primes{2};
  {
    int64_t x = F.n;
    for (int64_t p = 2; p <= x; p++) {
      if (x % p == 0) {
        int v = DiscriminantValuationCharacters(G, XK, p);
        for (int i = 0; i < v; i++) {
          dK *= static_cast<long>(p);
        }
        while (x % p == 0) {
          x /= p;
        }
      }
    }
  }
  mpq_class detB = DeterminantMat(B) * mpq_class(den * dK);
  for (int i = 1; i < B.rows(); i++) {
    detB *= mpq_class(den * dK);
  }
  mpz_class num = detB.get_num();
  if (num < 0) {
    num = -num;
  }
  if (num > 1) {
    for (mpz_class const &p : FactorsInt(num)) {
      if (p > mpz_class(1L << 31)) {
        std::cerr << "ABQUAD: the prime " << p << " is too large\n";
        throw TerminalException{1};
      }
      primes.insert(p.get_si());
    }
  }
  std::map<int64_t, std::optional<mpq_class>> res;
  for (int64_t q : primes) {
    if (q != 2) {
      // the sum over the primes p of K over q of f_p v_p(det phi)
      int c = ValuationQRational(detB, q) -
              m * DiscriminantValuationCharacters(G, XK, q);
      if (c == 0) {
        continue;
      }
    }
    res[q] = QuadLocalFactorAtPrime(B, Omega, q);
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_ABELIANQUADRATICMASS_H_
// clang-format on
