// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ABELIANHERMITIANMASS_H_
#define SRC_GENERA_ABELIANHERMITIANMASS_H_

// clang-format off
#include "NumberFieldElement.h"
#include "ClassicalGroupOrders.h"
#include "PolynomialModQ.h"
#include "MAT_MatrixModPrime.h"
#include "CyclotomicHermitianMass.h"
#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  Masses of Hermitian lattices over an abelian CM field E that need not be
  cyclotomic: E the subfield of Q(zeta_n) fixed by H (AbelianCMField), K
  its maximal real subfield.

  The lattice is given as a Z-lattice X with the Gram matrix B of the trace
  form b = Tr_{E/Q} h, and the action of a Z-basis omega_1, ..., omega_n of
  O_E (matrices on the row vectors). The standard mass is that of
  StandardHermitianMassAbelian. The local factors are computed prime by
  prime with the machinery of the cyclotomic case, on a generator theta_q
  of E in O_E chosen so that q does not divide [O_E : Z[theta_q]]: then
  Z_q[theta_q] = O_E (x) Z_q, the primes of E over q are the (q, g(theta_q))
  for the irreducible factors g of the minimal polynomial f modulo q
  (Dedekind), with the multiplicity of g as ramification index, and
  g(theta_q) is a uniformizer when that index exceeds 1. Such a theta_q is
  found among small combinations of the omega_i, v_q(disc f) = v_q(disc
  O_E) being the test; a prime dividing every index (a common index
  divisor) leaves the mass not computed.

  The primes of E over q are typed by the complex conjugation (the adjoint
  of theta_q for b, a polynomial in theta_q): split over K when conjugate
  to another one, otherwise ramified in E / K when their ramification index
  is twice that of K at q (from the characters), and inert otherwise.
  Unramified primes of K, including those over primes ramified in E / Q,
  get the local factor of HermitianLocalFactorUnramified, the different of
  E / Q being (f'(theta_q)) locally. A ramified prime P of K is treated
  on the Hermitian form over E of vectors spanning X (x) Z_q over O_E (x)
  Z_q, with v_P and the residues read through the idempotent of P
  (LocalPrimeK), so whatever the other primes over q: for q odd the
  Jordan splitting, the residue field F_q[x] / g of order q^f for the
  square classes, and the group orders of Hecke at q^{f_K}; for q = 2 the
  chain method with the relative different exponent.
 */

// Tr(theta^k), k = 0, ..., 2 deg, from the minimal polynomial by Newton.
template <typename T>
inline std::vector<T> PowerTraces(std::vector<int64_t> const &f,
                                          int kmax) {
  int deg = static_cast<int>(f.size()) - 1;
  // f = x^deg + c_{deg-1} x^{deg-1} + ... + c_0
  std::vector<T> p(kmax + 1, T(0));
  p[0] = deg;
  for (int k = 1; k <= kmax; k++) {
    T s(0);
    // Newton: p_k + c_{deg-1} p_{k-1} + ... = -k c_{deg-k} for k <= deg
    for (int i = 1; i <= std::min(k - 1, deg); i++) {
      s += UniversalScalarConversion<T, int64_t>(f[deg - i]) * p[k - i];
    }
    if (k <= deg) {
      s += UniversalScalarConversion<T, int64_t>(k) *
           UniversalScalarConversion<T, int64_t>(f[deg - k]);
    }
    p[k] = -s;
  }
  return p;
}

/*
  The minimal polynomial of a matrix S acting on an E-vector space (rows)
  through a generator of E: the first dependence of v, v S, v S^2, ... for
  a nonzero v, which is the minimal polynomial of the element. Also the
  coordinates of a vector w = v M in v, v S, ..., v S^{deg-1} (M in E).
 */
template <typename T> struct KrylovData {
  std::vector<T> minpoly;
  MyMatrix<T> basis;
};

template <typename T>
inline KrylovData<T> KrylovMinimalPolynomial(MyMatrix<T> const &S) {
  int r = S.rows();
  MyVector<T> v = ZeroVector<T>(r);
  v(0) = 1;
  std::vector<MyVector<T>> pw{v};
  for (int k = 1; k <= r; k++) {
    MyVector<T> w = S.transpose() * pw.back();
    std::vector<MyVector<T>> rows = pw;
    rows.push_back(w);
    MyMatrix<T> M = MatrixFromVectorFamily(rows);
    MyMatrix<T> ker = NullspaceMat(M);
    if (ker.rows() > 0) {
      // sum_i ker(0, i) v S^i = 0, normalized at i = k
      std::vector<T> mp(k + 1);
      for (int i = 0; i <= k; i++) {
        mp[i] = ker(0, i) / ker(0, k);
      }
      MyMatrix<T> basis(k, r);
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < r; j++) {
          basis(i, j) = pw[i](j);
        }
      }
      return {mp, basis};
    }
    pw.push_back(w);
  }
  std::cerr << "ABHERM: no minimal polynomial\n";
  throw TerminalException{1};
}

/*
  A generator theta of E in O_E with q prime to the index, its matrix S on
  X, its minimal polynomial f and the conjugate of theta as a polynomial in
  theta. Nothing when no small combination works.
 */
template <typename T> struct LocalGenerator {
  MyMatrix<T> S;
  std::vector<int64_t> f;
  std::vector<T> conj;
};

template <typename T>
inline std::optional<LocalGenerator<T>>
FindLocalGenerator(MyMatrix<T> const &B,
                   std::vector<MyMatrix<T>> const &Omega, int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int deg = static_cast<int>(Omega.size());
  int r = B.rows();
  int m = r / deg;
  // disc O_E = det(Tr(omega_i omega_j)), Tr_{E/Q} = Tr on X / m
  MyMatrix<T> Tr(deg, deg);
  for (int i = 0; i < deg; i++) {
    for (int j = 0; j < deg; j++) {
      MyMatrix<T> P = Omega[i] * Omega[j];
      T t(0);
      for (int u = 0; u < r; u++) {
        t += P(u, u);
      }
      Tr(i, j) = t / UniversalScalarConversion<T, int64_t>(m);
    }
  }
  int vdisc = RationalValuation(DeterminantMat(Tr), q);
  MyMatrix<T> Binv = Inverse(B);
  // combinations with coefficients in [-2, 2], by increasing size
  std::vector<std::vector<int>> cands;
  int n_tot = MyPow<int>(5, deg);
  for (int u = 0; u < n_tot; u++) {
    std::vector<int> c(deg);
    int t = u;
    for (int i = 0; i < deg; i++) {
      c[i] = t % 5 - 2;
      t /= 5;
    }
    cands.push_back(c);
  }
  std::stable_sort(cands.begin(), cands.end(),
                   [](std::vector<int> const &a, std::vector<int> const &b) {
                     int na = 0, nb = 0;
                     for (int x : a) {
                       na += std::abs(x);
                     }
                     for (int x : b) {
                       nb += std::abs(x);
                     }
                     return na < nb;
                   });
  for (auto const &c : cands) {
    MyMatrix<T> S = ZeroMatrix<T>(r, r);
    for (int i = 0; i < deg; i++) {
      if (c[i] != 0) {
        S += Omega[i] * T(c[i]);
      }
    }
    KrylovData<T> kd = KrylovMinimalPolynomial(S);
    if (static_cast<int>(kd.minpoly.size()) != deg + 1) {
      continue;
    }
    std::vector<int64_t> f;
    bool integral = true;
    for (auto const &x : kd.minpoly) {
      if (GetDenominator_z(x) != 1) {
        integral = false;
      }
      f.push_back(UniversalScalarConversion<int64_t, Tint>(GetNumerator_z(x)));
    }
    if (!integral) {
      continue;
    }
    // disc f = det(Tr(theta^{i+j}))
    std::vector<T> pt = PowerTraces<T>(f, 2 * deg);
    MyMatrix<T> D(deg, deg);
    for (int i = 0; i < deg; i++) {
      for (int j = 0; j < deg; j++) {
        D(i, j) = pt[i + j];
      }
    }
    if (RationalValuation(DeterminantMat(D), q) != vdisc) {
      continue;
    }
    // the conjugate: S* = (B^{-1} S B)^T on the rows, v S* in the Krylov
    // basis of v
    MyMatrix<T> Sstar = (Binv * S * B).transpose();
    MyVector<T> v = GetMatrixRow(kd.basis, 0);
    MyVector<T> w = Sstar.transpose() * v;
    // w = sum_i a_i v S^i
    MyMatrix<T> M = kd.basis;
    std::optional<MyVector<T>> sol = SolutionMat(M, w);
    if (!sol) {
      std::cerr << "ABHERM: the conjugate of theta is not in Q(theta)\n";
      throw TerminalException{1};
    }
    std::vector<T> conj(deg);
    for (int i = 0; i < deg; i++) {
      conj[i] = (*sol)(i);
    }
    return LocalGenerator<T>{S, f, conj};
  }
  return {};
}

// The matrix of g(S) for a polynomial g with integer coefficients, or
// residues read as integers.
template <typename T, typename Tcoef>
inline MyMatrix<T> PolyOfMatrixQ(std::vector<Tcoef> const &g,
                                 MyMatrix<T> const &S) {
  int r = S.rows();
  MyMatrix<T> res = ZeroMatrix<T>(r, r);
  for (size_t i = g.size(); i-- > 0;) {
    res = res * S;
    T g_i = UniversalScalarConversion<T, Tcoef>(g[i]);
    for (int j = 0; j < r; j++) {
      res(j, j) += g_i;
    }
  }
  return res;
}

/*
  A prime of E over q: the factor g, its multiplicity (the ramification
  index over q), the conjugate factor, the type over K (0 split, 1 inert,
  2 ramified), f_K = log_q of the residue field of the prime of K below,
  and the idempotent of the prime on X modulo q^N.
 */
struct AbelianPrime {
  PolyQ g;
  int e;
  size_t conj_index;
  int type;
  int f_K;
};

template <typename T>
inline std::vector<AbelianPrime>
AbelianPrimesOverQ(LocalGenerator<T> const &lg, int64_t q, int64_t eK) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  std::vector<std::pair<PolyQ, int>> fac =
      FactorModQWithMultiplicity(PolyReduceSigned(lg.f, q), q);
  // the conjugate of the prime of g: the factor dividing g(conj(x))
  auto compose_conj = [&](PolyQ const &g) -> PolyQ {
    // g(c(x)) modulo f and q, c the conjugate of theta (q-integral)
    PolyQ c;
    for (auto const &x : lg.conj) {
      c.push_back(ResidueModInt<T, Tint>(x, static_cast<T_residue>(q)));
    }
    PolyQ fq = PolyReduceSigned(lg.f, q);
    PolyQ res = {0};
    for (size_t i = g.size(); i-- > 0;) {
      res = PolyDivModQ(PolyMulQ(res, c, q), fq, q).second;
      if (res.empty()) {
        res = {0};
      }
      res[0] = (res[0] + g[i]) % q;
      PolyTrim(res);
      if (res.empty()) {
        res = {0};
      }
    }
    return res;
  };
  std::vector<AbelianPrime> res;
  for (size_t i = 0; i < fac.size(); i++) {
    PolyQ gc = compose_conj(fac[i].first);
    size_t j = i;
    for (size_t u = 0; u < fac.size(); u++) {
      PolyQ rr = PolyDivModQ(gc, fac[u].first, q).second;
      PolyTrim(rr);
      if (rr.empty()) {
        j = u;
        break;
      }
    }
    AbelianPrime ap;
    ap.g = fac[i].first;
    ap.e = fac[i].second;
    ap.conj_index = j;
    int deg_g = static_cast<int>(fac[i].first.size()) - 1;
    if (j != i) {
      ap.type = 0;
      ap.f_K = deg_g;
    } else if (ap.e == 2 * eK) {
      ap.type = 2;
      ap.f_K = deg_g;
    } else {
      ap.type = 1;
      ap.f_K = deg_g / 2;
    }
    res.push_back(ap);
  }
  return res;
}

/*
  The idempotent on X of the prime (q, g(theta)) of ramification index e,
  lifted modulo q^N: the CRT idempotent of g^e in F_q[x] / f, refined by e
  <- 3 e^2 - 2 e^3 modulo (f, q^N), evaluated at theta.
 */
template <typename Tint>
inline PolyZ<Tint> AbelianPrimeIdempotentPoly(std::vector<int64_t> const &f,
                                        PolyQ const &g, int e, int64_t q,
                                        int N) {
  Tint Q = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(q), N);
  PolyQ fq = PolyReduceSigned(f, q);
  PolyQ P = {1};
  for (int i = 0; i < e; i++) {
    P = PolyMulQ(P, g, q);
  }
  PolyQ Qp = PolyDivModQ(fq, P, q).first;
  PolyQ e0 = PolyMulQ(Qp, PolyInverseModQ(Qp, P, q), q);
  PolyZ<Tint> ez = PolyZOfPolyQ<Tint>(e0);
  ez = PolyZReduce(ez, f, Q);
  for (int it = 0; it < 64; it++) {
    PolyZ<Tint> e2 = PolyZReduce(PolyZMul(ez, ez), f, Q);
    PolyZ<Tint> e3 = PolyZReduce(PolyZMul(e2, ez), f, Q);
    PolyZ<Tint> en(e2.size());
    for (size_t k = 0; k < e2.size(); k++) {
      en[k] = 3 * e2[k] - 2 * e3[k];
    }
    en = PolyZReduce(en, f, Q);
    if (en == ez) {
      break;
    }
    ez = en;
  }
  return ez;
}

template <typename T>
inline MyMatrix<T> AbelianPrimeIdempotent(LocalGenerator<T> const &lg,
                                                  PolyQ const &g, int e,
                                                  int64_t q, int N) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int r = lg.S.rows();
  PolyZ<Tint> ez = AbelianPrimeIdempotentPoly<Tint>(lg.f, g, e, q, N);
  MyMatrix<T> E = ZeroMatrix<T>(r, r);
  MyMatrix<T> pw = IdentityMat<T>(r);
  for (size_t k = 0; k < ez.size(); k++) {
    E += pw * T(ez[k]);
    pw = pw * lg.S;
  }
  return E;
}

/*
  The Jordan data at a prime of E over q unramified over K (split: one of
  the pair), as HermitianJordanAtPrime does for Q(zeta_d).
 */
template <typename T>
inline HermitianLocalJordan
AbelianJordanAtPrime(MyMatrix<T> const &B, LocalGenerator<T> const &lg,
                     std::vector<AbelianPrime> const &primes, size_t i_prime,
                     int64_t q) {
  int r = B.rows();
  AbelianPrime const &ap = primes[i_prime];
  int f_res = static_cast<int>(ap.g.size()) - 1;
  int vdet = RationalValuation(DeterminantMat(B), q);
  std::vector<int64_t> f = lg.f;
  // the precision covers the discriminant group and the different
  int vdiscf = 0;
  {
    std::vector<T> pt = PowerTraces<T>(f, 2 * (static_cast<int>(f.size()) - 1));
    int deg = static_cast<int>(f.size()) - 1;
    MyMatrix<T> D(deg, deg);
    for (int a = 0; a < deg; a++) {
      for (int b = 0; b < deg; b++) {
        D(a, b) = pt[a + b];
      }
    }
    vdiscf = RationalValuation(DeterminantMat(D), q);
  }
  int N = vdet + vdiscf + 6;
  MyMatrix<T> Id = IdentityMat<T>(r);
  MyMatrix<T> E = AbelianPrimeIdempotent(lg, ap.g, ap.e, q, N);
  // the uniformizer: q, or g(theta) when the prime is ramified over q
  MyMatrix<T> Pi =
      (ap.e == 1) ? MyMatrix<T>(Id * UniversalScalarConversion<T, int64_t>(q))
                  : PolyOfMatrixQ(ap.g, lg.S);
  // delta = v_P(f'(theta)): v_q(det(f'(S) E + 1 - E)) on X is m f v_P
  int deg = static_cast<int>(f.size()) - 1;
  int m_total = r / deg;
  std::vector<int64_t> fprime;
  for (size_t k = 1; k < f.size(); k++) {
    fprime.push_back(static_cast<int64_t>(k) * f[k]);
  }
  MyMatrix<T> FP = PolyOfMatrixQ(fprime, lg.S);
  // the O_P-rank m
  T_small_prime q_u = static_cast<T_small_prime>(q);
  ModMatrix Emod = ModMatrixOfMatrix(E, q_u);
  int rk = static_cast<int>(RankModPrime(Emod, q_u));
  int m = rk / (ap.e * f_res);
  int delta = 0;
  {
    MyMatrix<T> M = FP * E + Id - E;
    int vv = RationalValuation(DeterminantMat(M), q);
    // on X the prime part has Z_q-rank m e f, and det = N(f')^m locally
    delta = vv / (m * f_res);
  }
  (void)m_total;
  HermitianLocalJordan hj;
  hj.type = ap.type;
  hj.f_K = ap.f_K;
  MyMatrix<T> Binv = Inverse(B);
  std::vector<int> len;
  MyMatrix<T> Mk = Binv * E;
  while (true) {
    len.push_back(ValuationIndexOverlattice(Mk, q));
    if (len.back() == 0 || len.size() > static_cast<size_t>(4 * (vdet + vdiscf) + 8)) {
      break;
    }
    Mk = Mk * Pi;
  }
  std::vector<int> c;
  for (size_t k = 0; k + 1 < len.size(); k++) {
    c.push_back((len[k] - len[k + 1]) / f_res);
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
    hj.blocks.push_back({kv.first - delta, kv.second});
  }
  return hj;
}

/*
  The arithmetic at a prime p = (q, g(theta)) of K, theta = theta_q with q
  prime to the index: the elements of K are kept exactly on the powers of
  theta, v_p(x) is v_q of the norm of eps x + 1 - eps over f (eps the
  idempotent of p in Z_q[theta], known modulo q^N for an N raised until
  the norm is nonzero modulo q^N), and the residue of x in F_q[x] / g is
  read on eps x. This holds whatever the other primes over q are. The
  same holds for a prime P of a CM field E, and LocalOKGram then gives the
  Hermitian form over E (h linear in its first argument).
 */
template <typename T> struct LocalPrimeK {
  using Tint = typename underlying_z_ring<T>::ring_type;
  CyclotomicField<T> K;
  int64_t q;
  PolyQ g;
  int e;
  int f_res;
  // the order of the residue field, an abstract size
  Tint qf;
  ElemE<T> pi;
  ElemE<T> pi_inv;
  std::map<int, PolyZ<Tint>> eps;
};

template <typename T>
inline LocalPrimeK<T> MakeLocalPrimeK(LocalGenerator<T> const &lg, PolyQ const &g,
                                   int e, int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  LocalPrimeK<T> P;
  int deg = static_cast<int>(lg.f.size()) - 1;
  P.K.d = 0;
  P.K.deg = deg;
  P.K.Phi = lg.f;
  P.K.conj = lg.conj;
  P.q = q;
  P.g = g;
  P.e = e;
  P.f_res = static_cast<int>(g.size()) - 1;
  P.qf = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(q), P.f_res);
  if (e == 1) {
    P.pi = ElemE<T>(deg, T(0));
    P.pi[0] = UniversalScalarConversion<T, int64_t>(q);
  } else {
    P.pi = EPolyOfTheta(g, P.K);
  }
  P.pi_inv = EInv(P.pi, P.K);
  return P;
}

template <typename T>
inline PolyZ<typename underlying_z_ring<T>::ring_type> const &LocalEps(LocalPrimeK<T> &P, int N) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  auto it = P.eps.find(N);
  if (it == P.eps.end()) {
    it = P.eps.emplace(N, AbelianPrimeIdempotentPoly<Tint>(P.K.Phi, P.g, P.e, P.q, N))
             .first;
  }
  return it->second;
}

// x = y / d with y in Z[theta] and d > 0
template <typename T>
inline std::pair<PolyZ<typename underlying_z_ring<T>::ring_type>, typename underlying_z_ring<T>::ring_type>
ClearDenominators(ElemE<T> const &x) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  Tint d(1);
  for (auto const &c : x) {
    d = LCMpair(d, GetDenominator_z(c));
  }
  PolyZ<Tint> y;
  for (auto const &c : x) {
    T t = c * T(d);
    y.push_back(GetNumerator_z(t));
  }
  return {y, d};
}

template <typename T>
inline int LocalVal(ElemE<T> const &x, LocalPrimeK<T> &P) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  if (EIsZero(x)) {
    return 1 << 28;
  }
  auto [y, d] = ClearDenominators(x);
  int k = IntegerValuation(d, P.q);
  int N = P.f_res * (P.e * k + 4) + 4;
  while (true) {
    Tint Q = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(P.q), N);
    PolyZ<Tint> const &eps = LocalEps(P, N);
    PolyZ<Tint> w = PolyZReduce(PolyZMul(eps, y), P.K.Phi, Q);
    for (size_t i = 0; i < eps.size(); i++) {
      w[i] -= eps[i];
    }
    w[0] += 1;
    w = PolyZReduce(w, P.K.Phi, Q);
    ElemE<T> we;
    for (auto const &c : w) {
      we.push_back(T(c));
    }
    T det = DeterminantMat(EMulMatrix(we, P.K));
    Tint dn = GetNumerator_z(det);
    if (ResInt(dn, Q) != 0) {
      int v = IntegerValuation(dn, P.q);
      return v / P.f_res - P.e * k;
    }
    N *= 2;
  }
}

// The residue in F_q[x] / g of x with v_p(x) >= 0.
template <typename T>
inline PolyQ LocalResidue(ElemE<T> const &x, LocalPrimeK<T> &P) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  auto [y, d] = ClearDenominators(x);
  int k = IntegerValuation(d, P.q);
  Tint qk = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(P.q), k);
  Tint Q = qk * UniversalScalarConversion<Tint, int64_t>(P.q);
  PolyZ<Tint> z = PolyZReduce(PolyZMul(LocalEps(P, k + 1), y), P.K.Phi, Q);
  T_small_prime q_u = static_cast<T_small_prime>(P.q);
  T_residue inv = ResInv(ResModOf(QuoInt(d, qk), q_u), q_u);
  PolyQ res;
  for (auto const &c : z) {
    if (ResInt(c, qk) != 0) {
      std::cerr << "ABQUAD: the residue of a non integral element\n";
      throw TerminalException{1};
    }
    res.push_back(ResMul(ResModOf(QuoInt(c, qk), q_u), inv, q_u));
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
  The Gram matrix over K of phi on the vectors v_i of X: phi(v_i, v_j) on
  the powers of theta, from the traces Tr(theta^t phi(v_i, v_j)) = b(theta^t
  v_i, v_j), theta being any generator of K.
 */
template <typename T>
inline std::vector<std::vector<ElemE<T>>>
LocalOKGramFromVectors(MyMatrix<T> const &B, LocalGenerator<T> const &lg,
                       std::vector<MyVector<T>> const &chosen) {
  int r = B.rows();
  int deg = static_cast<int>(lg.f.size()) - 1;
  int m = static_cast<int>(chosen.size());
  std::vector<MyMatrix<T>> Spow{IdentityMat<T>(r)};
  for (int t = 1; t < deg; t++) {
    Spow.push_back(Spow.back() * lg.S);
  }
  std::vector<T> pt = PowerTraces<T>(lg.f, 2 * deg);
  MyMatrix<T> Tm(deg, deg);
  for (int a = 0; a < deg; a++) {
    for (int b = 0; b < deg; b++) {
      Tm(a, b) = pt[a + b];
    }
  }
  MyMatrix<T> Tinv = Inverse(Tm);
  std::vector<std::vector<ElemE<T>>> H(m, std::vector<ElemE<T>>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      ElemE<T> val(deg, T(0));
      for (int t = 0; t < deg; t++) {
        MyVector<T> w = Spow[t].transpose() * chosen[i];
        T bt = w.dot(B * chosen[j]);
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
  The Gram matrix over K of phi on m vectors of X spanning X (x) Z_q over
  O_K (x) Z_q: standard vectors first, then small combinations.
 */
template <typename T> struct LocalOKBasisGram {
  std::vector<MyVector<T>> basis;
  std::vector<std::vector<ElemE<T>>> H;
};

template <typename T>
inline LocalOKBasisGram<T>
LocalOKGramBasis(MyMatrix<T> const &B, LocalGenerator<T> const &lg,
                 int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int r = B.rows();
  int deg = static_cast<int>(lg.f.size()) - 1;
  int m = r / deg;
  T_small_prime q_u = static_cast<T_small_prime>(q);
  auto red = [&](T const &x) -> T_residue {
    return ResidueModInt<T, Tint>(x, q_u);
  };
  std::vector<MyMatrix<T>> Spow{IdentityMat<T>(r)};
  for (int t = 1; t < deg; t++) {
    Spow.push_back(Spow.back() * lg.S);
  }
  auto orbit_rows = [&](MyVector<T> const &v) -> ModMatrix {
    ModMatrix rows(deg, r);
    for (int t = 0; t < deg; t++) {
      MyVector<T> w = Spow[t].transpose() * v;
      for (int j = 0; j < r; j++) {
        rows(t, j) = red(w(j));
      }
    }
    return rows;
  };
  std::vector<MyVector<T>> chosen;
  ModMatrix acc(0, r);
  auto try_add = [&](MyVector<T> const &v) -> void {
    ModMatrix rows = orbit_rows(v);
    ModMatrix test(acc.rows() + rows.rows(), r);
    test.topRows(acc.rows()) = acc;
    test.bottomRows(rows.rows()) = rows;
    if (RankModPrime(test, q_u) == static_cast<size_t>(test.rows())) {
      chosen.push_back(v);
      acc = test;
    }
  };
  for (int k = 0; k < r && static_cast<int>(chosen.size()) < m; k++) {
    MyVector<T> v = ZeroVector<T>(r);
    v(k) = 1;
    try_add(v);
  }
  uint64_t seed = 1;
  for (int it = 0; it < 1000 && static_cast<int>(chosen.size()) < m; it++) {
    MyVector<T> v = ZeroVector<T>(r);
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
  return {chosen, LocalOKGramFromVectors(B, lg, chosen)};
}

template <typename T>
inline std::vector<std::vector<ElemE<T>>>
LocalOKGram(MyMatrix<T> const &B, LocalGenerator<T> const &lg, int64_t q) {
  return LocalOKGramBasis(B, lg, q).H;
}

// x pi^{-k}
template <typename T>
inline ElemE<T> LocalShift(ElemE<T> x, int k, LocalPrimeK<T> const &P) {
  if (k > 0) {
    ElemE<T> pinv = P.pi_inv.empty() ? EInv(P.pi, P.K) : P.pi_inv;
    for (int i = 0; i < k; i++) {
      x = EMul(x, pinv, P.K);
    }
  }
  for (int i = 0; i < -k; i++) {
    x = EMul(x, P.pi, P.K);
  }
  return x;
}

/*
  The local factor at a ramified prime P of K over the odd prime q: the
  formula of HermitianLocalFactorRamified with the residue field of P, of
  order q^f. The Hermitian Gram matrix over E is that of m vectors of X
  spanning X (x) Z_q over O_E (x) Z_q, and v_P and the residues are those
  of LocalPrimeK, so P need not be alone over q.
 */
template <typename T> struct AbelianLocalGram {
  CyclotomicField<T> K;
  int m;
  std::vector<std::vector<ElemE<T>>> H;
  std::shared_ptr<LocalPrimeK<T>> P;
};

template <typename T>
inline AbelianLocalGram<T> AbelianLocalHermitianGram(MyMatrix<T> const &B,
                                                  LocalGenerator<T> const &lg,
                                                  AbelianPrime const &ap,
                                                  int64_t q) {
  std::shared_ptr<LocalPrimeK<T>> P =
      std::make_shared<LocalPrimeK<T>>(MakeLocalPrimeK(lg, ap.g, ap.e, q));
  CyclotomicField<T> K = P->K;
  K.val_div = static_cast<int>(ap.g.size()) - 1;
  K.valuation = [P](ElemE<T> const &x) -> int { return LocalVal(x, *P); };
  std::vector<std::vector<ElemE<T>>> H = LocalOKGram(B, lg, q);
  int m = static_cast<int>(H.size());
  return {K, m, H, P};
}

template <typename T>
inline T AbelianLocalFactorRamified(MyMatrix<T> const &B,
                                            LocalGenerator<T> const &lg,
                                            AbelianPrime const &ap, int64_t q) {
  AbelianLocalGram<T> alg = AbelianLocalHermitianGram(B, lg, ap, q);
  CyclotomicField<T> const &K = alg.K;
  int m = alg.m;
  std::vector<std::vector<ElemE<T>>> const &H = alg.H;
  std::map<int, std::pair<int, ElemE<T>>> blocks;
  for (auto &kv : HermitianJordanSplit(H, K, q)) {
    int rk = 0;
    ElemE<T> det = EZetaPower(0, K);
    for (auto &G : kv.second) {
      rk += static_cast<int>(G.size());
      det = EMul(det, HermitianSmallDet(G, K), K);
    }
    blocks[kv.first] = {rk, det};
  }
  T qK = MyPow<T>(UniversalScalarConversion<T, int64_t>(q), ap.f_K);
  T fac = HeckeGroupOrderSp<T>(2 * (m / 2), qK);
  int64_t Nexp = 0;
  int after_total = m;
  // pi = g(theta), and pi pi-bar
  ElemE<T> pi = EPolyOfTheta(ap.g, K);
  ElemE<T> pipibar = EMul(pi, EConj(pi, K), K);
  ElemE<T> pinv = EInv(pipibar, K);
  // the square class of a P-unit in the residue field F_q[x] / g
  auto is_square = [&](ElemE<T> const &u) -> int {
    PolyQ ur = LocalResidue(u, *alg.P);
    if (ur.empty()) {
      return 0;
    }
    return ResidueQuadraticCharacter(ur, ap.g, q);
  };
  for (auto &kv : blocks) {
    int s = kv.first;
    int mi = kv.second.first;
    after_total -= mi;
    int ri = after_total;
    int64_t fl = (s >= 0) ? s / 2 : -((-s + 1) / 2);
    Nexp -= fl * mi * mi;
    if (s % 2 != 0) {
      Nexp -= static_cast<int64_t>(mi + 1) * (mi / 2);
      fac /= HeckeGroupOrderSp<T>(mi, qK);
    } else {
      // divided by (pi pi-bar)^{w / 2}, w of either sign (a block of
      // negative scale when the form has denominators at P)
      ElemE<T> u = kv.second.second;
      int w = EValuation(u, K, q);
      for (int i = 0; i < w / 2; i++) {
        u = EMul(u, pinv, K);
      }
      for (int i = 0; i < -w / 2; i++) {
        u = EMul(u, pipibar, K);
      }
      if ((mi * (mi - 1) / 2) % 2 == 1) {
        for (auto &x : u) {
          x = -x;
        }
      }
      fac /= HeckeGroupOrderO<T>(mi, is_square(u), qK);
    }
    int64_t smi = static_cast<int64_t>(s) * mi;
    Nexp -= smi * ri;
    Nexp += smi * m / 2;
  }
  if (m % 2 == 0) {
    Nexp += m / 2;
  }
  T qpow(1);
  for (int64_t i = 0; i < (Nexp >= 0 ? Nexp : -Nexp); i++) {
    if (Nexp >= 0) {
      qpow *= qK;
    } else {
      qpow /= qK;
    }
  }
  return qpow * fac;
}

/*
  The dyadic ramified prime P of E over the prime p of K, following HermitianLocalFactorDyadic with what was particular to Q(zeta_{2^a})
  made general: the exponent d of the relative different of E / K at P
  (2 there), which is also the conductor exponent of E_P / K_p, so that
  v_p(Tr(c O_E)) = floor((v_P(c) + d) / 2), a unit u of K is a norm iff
  u = y y-bar modulo p^d for a unit y of O_P / P^{2d}, and the normic
  defect of a non-norm x is v_p(x) + d - 1; and the residue field
  F_P = F_2[x] / g of order 2^f, whose elements are the sum b_j theta^j,
  j < f, b_j in {0, 1}, over which the lines of L / pi L are taken.
 */
template <typename T> struct GenericDyadic {
  CyclotomicField<T> K;
  ElemE<T> one;
  ElemE<T> pi;
  ElemE<T> pi_inv;
  ElemE<T> pipibar;
  ElemE<T> pipibar_inv;
  int d;
  int f;
  // the representatives of F_P, the zero first
  std::vector<ElemE<T>> reps;
};

template <typename T>
inline int GDVal(ElemE<T> const &x, GenericDyadic<T> const &C) {
  return EValuation(x, C.K, 2);
}

template <typename T>
inline int GDTraceVal(ElemE<T> const &c, GenericDyadic<T> const &C) {
  int v = GDVal(c, C);
  if (v >= (1 << 27)) {
    return v;
  }
  int num = v + C.d;
  int vp = (num >= 0) ? num / 2 : -((-num + 1) / 2);
  return 2 * vp;
}

template <typename T>
inline int GDNormVal(std::vector<std::vector<ElemE<T>>> const &G,
                     GenericDyadic<T> const &C) {
  int n = static_cast<int>(G.size());
  int res = 1 << 28;
  for (int i = 0; i < n; i++) {
    res = std::min(res, GDVal(G[i][i], C));
    for (int j = i + 1; j < n; j++) {
      res = std::min(res, GDTraceVal(G[i][j], C));
    }
  }
  return res;
}

template <typename T>
inline bool GDIsLocalNorm(ElemE<T> x, GenericDyadic<T> const &C) {
  CyclotomicField<T> const &K = C.K;
  int w = GDVal(x, C);
  for (int i = 0; i < w / 2; i++) {
    x = EMul(x, C.pipibar_inv, K);
  }
  for (int i = 0; i < -w / 2; i++) {
    x = EMul(x, C.pipibar, K);
  }
  size_t nr = C.reps.size();
  int len = 2 * C.d;
  // y = sum_{i < len} a_i pi^i, a_0 a nonzero representative
  std::vector<ElemE<T>> pipow{C.one};
  for (int i = 1; i < len; i++) {
    pipow.push_back(EMul(pipow.back(), C.pi, K));
  }
  std::vector<size_t> idx(len, 0);
  idx[0] = 1;
  while (true) {
    ElemE<T> y(K.deg, T(0));
    for (int i = 0; i < len; i++) {
      if (idx[i] == 0) {
        continue;
      }
      ElemE<T> t = EMul(C.reps[idx[i]], pipow[i], K);
      for (int u = 0; u < K.deg; u++) {
        y[u] += t[u];
      }
    }
    if (GDVal(ESub(x, EMul(y, EConj(y, K), K)), C) >= 2 * C.d) {
      return true;
    }
    int i = 0;
    while (i < len) {
      idx[i]++;
      if (idx[i] < nr) {
        break;
      }
      idx[i] = (i == 0) ? 1 : 0;
      i++;
    }
    if (i == len) {
      return false;
    }
  }
}

// the lines of (O_P / P)^m: the first nonzero coordinate is 1
inline std::vector<std::vector<size_t>> GDLines(int m, size_t nr) {
  std::vector<std::vector<size_t>> res;
  for (int lead = 0; lead < m; lead++) {
    size_t n_tail = 1;
    for (int i = lead + 1; i < m; i++) {
      n_tail *= nr;
    }
    for (size_t u = 0; u < n_tail; u++) {
      std::vector<size_t> a(m, 0);
      a[lead] = 1;
      size_t t = u;
      for (int i = lead + 1; i < m; i++) {
        a[i] = t % nr;
        t /= nr;
      }
      res.push_back(a);
    }
  }
  return res;
}

// L + R pi^{-1} y, y = sum_i a_i x_i, a_lead = 1
template <typename T>
inline std::vector<std::vector<ElemE<T>>>
GDOverlattice(std::vector<std::vector<ElemE<T>>> const &H,
              std::vector<size_t> const &a, GenericDyadic<T> const &C) {
  int m = static_cast<int>(H.size());
  ElemE<T> zero(C.K.deg, T(0));
  std::vector<std::vector<ElemE<T>>> U(m, std::vector<ElemE<T>>(m, zero));
  int j = -1;
  for (int i = 0; i < m; i++) {
    U[i][i] = C.one;
    if (j < 0 && a[i] != 0) {
      j = i;
    }
  }
  for (int i = 0; i < m; i++) {
    U[j][i] = (a[i] != 0) ? EMul(C.reps[a[i]], C.pi_inv, C.K) : zero;
  }
  return HermitianCongruence(H, U, C.K);
}

// the kernel of x -> sum_i a_i c_i modulo pi, a_lead = 1
template <typename T>
inline std::vector<std::vector<ElemE<T>>>
GDSublattice(std::vector<std::vector<ElemE<T>>> const &H,
             std::vector<size_t> const &a, GenericDyadic<T> const &C) {
  int m = static_cast<int>(H.size());
  ElemE<T> zero(C.K.deg, T(0));
  std::vector<std::vector<ElemE<T>>> U(m, std::vector<ElemE<T>>(m, zero));
  int j = -1;
  for (int i = 0; i < m; i++) {
    if (a[i] != 0) {
      j = i;
      break;
    }
  }
  for (int i = 0; i < m; i++) {
    if (i == j) {
      U[i][j] = C.pi;
    } else {
      U[i][i] = C.one;
      if (a[i] != 0) {
        U[i][j] = ESub(zero, C.reps[a[i]]);
      }
    }
  }
  return HermitianCongruence(H, U, C.K);
}

template <typename T> struct GDSymbolEntry {
  int rank;
  int scale;
  bool normal;
  int v;
  ElemE<T> det;
};

template <typename T>
inline std::vector<GDSymbolEntry<T>>
GDKirschmerSymbol(std::vector<std::vector<ElemE<T>>> const &H,
                  GenericDyadic<T> const &C) {
  std::map<int, std::vector<SmallGram<T>>> blocks = HermitianJordanSplit(H, C.K, 2);
  std::vector<int> scales, nvals, ranks;
  std::vector<ElemE<T>> dets;
  for (auto &kv : blocks) {
    int nv = 1 << 28;
    int rk = 0;
    ElemE<T> det = C.one;
    for (auto &G : kv.second) {
      nv = std::min(nv, GDNormVal(G, C));
      rk += static_cast<int>(G.size());
      det = EMul(det, HermitianSmallDet(G, C.K), C.K);
    }
    scales.push_back(kv.first);
    nvals.push_back(nv);
    ranks.push_back(rk);
    dets.push_back(det);
  }
  std::vector<GDSymbolEntry<T>> sym;
  ElemE<T> det_partial = C.one;
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

template <typename T>
inline bool GDIsLocallyIsometric(std::vector<GDSymbolEntry<T>> const &S1,
                                 std::vector<GDSymbolEntry<T>> const &S2,
                                 GenericDyadic<T> const &C) {
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
  auto ratio = [&](size_t i) -> ElemE<T> {
    return EMul(S1[i].det, EInv(S2[i].det, C.K), C.K);
  };
  if (!GDIsLocalNorm(ratio(t - 1), C)) {
    return false;
  }
  for (size_t i = 0; i + 1 < t; i++) {
    ElemE<T> x = ratio(i);
    if (GDIsLocalNorm(x, C)) {
      continue;
    }
    int defect = GDVal(x, C) / 2 + C.d - 1;
    if (2 * defect < S1[i].v + S1[i + 1].v - 2 * S1[i].scale) {
      return false;
    }
  }
  return true;
}

/*
  The local factor at a dyadic ramified prime P, of relative different
  exponent d: the chain method of
  HermitianLocalFactorDyadic over the residue field of P.
 */
template <typename T>
inline T AbelianLocalFactorDyadic(MyMatrix<T> const &B,
                                          LocalGenerator<T> const &lg,
                                          AbelianPrime const &ap, int d) {
  AbelianLocalGram<T> alg = AbelianLocalHermitianGram(B, lg, ap, 2);
  GenericDyadic<T> C;
  C.K = alg.K;
  CyclotomicField<T> const &K = C.K;
  int deg = K.deg;
  int m = alg.m;
  C.one = EZetaPower(0, K);
  C.pi = EPolyOfTheta(ap.g, K);
  C.pi_inv = EInv(C.pi, K);
  C.pipibar = EMul(C.pi, EConj(C.pi, K), K);
  C.pipibar_inv = EInv(C.pipibar, K);
  C.d = d;
  C.f = static_cast<int>(ap.g.size()) - 1;
  for (int64_t u = 0; u < (int64_t(1) << C.f); u++) {
    ElemE<T> x(deg, T(0));
    for (int j = 0; j < C.f; j++) {
      if (u & (int64_t(1) << j)) {
        ElemE<T> t = EZetaPower(j, K);
        for (int k = 0; k < deg; k++) {
          x[k] += t[k];
        }
      }
    }
    C.reps.push_back(x);
  }
  int64_t qK = int64_t(1) << ap.f_K;
  std::vector<std::vector<size_t>> lines = GDLines(m, C.reps.size());
  auto maximal_factor = [&](std::vector<std::vector<ElemE<T>>> const &H) -> T {
    if (m % 2 == 1) {
      return RationalFraction<T>(1, 2);
    }
    std::vector<std::vector<ElemE<T>>> G = H;
    ElemE<T> det = C.one;
    for (int c = 0; c < m; c++) {
      int piv = -1;
      for (int i = c; i < m; i++) {
        if (!EIsZero(G[i][c])) {
          piv = i;
          break;
        }
      }
      if (piv < 0) {
        std::cerr << "ABHERM: degenerate Hermitian form\n";
        throw TerminalException{1};
      }
      if (piv != c) {
        std::swap(G[piv], G[c]);
        for (auto &x : det) {
          x = -x;
        }
      }
      det = EMul(det, G[c][c], K);
      ElemE<T> inv = EInv(G[c][c], K);
      for (int i = c + 1; i < m; i++) {
        ElemE<T> t = EMul(G[i][c], inv, K);
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
    if (GDIsLocalNorm(det, C)) {
      return T(1);
    }
    T qm = MyPow<T>(UniversalScalarConversion<T, int64_t>(qK), m);
    return (qm - 1) / UniversalScalarConversion<T, int64_t>(2 * qK + 2);
  };
  int nL = GDNormVal(alg.H, C);
  auto find_over = [&](std::vector<std::vector<ElemE<T>>> const &H)
      -> std::optional<std::vector<std::vector<ElemE<T>>>> {
    for (auto const &a : lines) {
      std::vector<std::vector<ElemE<T>>> Hn = GDOverlattice(H, a, C);
      if (GDNormVal(Hn, C) >= nL) {
        return Hn;
      }
    }
    return {};
  };
  std::vector<std::vector<std::vector<ElemE<T>>>> chain{alg.H};
  while (true) {
    std::optional<std::vector<std::vector<ElemE<T>>>> opt = find_over(chain.back());
    if (!opt) {
      break;
    }
    chain.push_back(*opt);
    if (chain.size() > 64) {
      std::cerr << "ABHERM: the chain of overlattices does not end\n";
      throw TerminalException{1};
    }
  }
  T fac = maximal_factor(alg.H);
  for (size_t i = 0; i + 1 < chain.size(); i++) {
    std::vector<GDSymbolEntry<T>> sym_i = GDKirschmerSymbol(chain[i], C);
    std::vector<GDSymbolEntry<T>> sym_n = GDKirschmerSymbol(chain[i + 1], C);
    int64_t A = 0, Bc = 0;
    for (auto const &a : lines) {
      if (GDIsLocallyIsometric(GDKirschmerSymbol(GDSublattice(chain[i + 1], a, C), C),
                               sym_i, C)) {
        A++;
      }
      if (GDIsLocallyIsometric(GDKirschmerSymbol(GDOverlattice(chain[i], a, C), C),
                               sym_n, C)) {
        Bc++;
      }
    }
    if (A == 0 || Bc == 0) {
      std::cerr << "ABHERM: a step of the chain has no isometric neighbor\n";
      throw TerminalException{1};
    }
    fac *= RationalFraction<T>(A, Bc);
  }
  return fac;
}

// sum over the characters of v_q(conductor): v_q(|disc|) of their field
inline int DiscriminantValuationCharacters(DirichletGroup const &G,
                                           std::vector<std::vector<int64_t>> const &X,
                                           int64_t q) {
  int res = 0;
  for (auto const &t : X) {
    int64_t c = PrimitiveCharacter(G, t).first;
    while (c % q == 0) {
      c /= q;
      res++;
    }
  }
  return res;
}

/*
  The Jordan data from an idempotent: the prime cut out by E on X, of
  residue degree f_res and ramification index e_ram over q, with the
  uniformizer Pi and delta = v_P(different of E / Q).
 */
template <typename T>
inline HermitianLocalJordan
JordanFromIdempotent(MyMatrix<T> const &B, MyMatrix<T> const &E,
                     MyMatrix<T> const &Pi, int f_res, int e_ram,
                     int delta, int type, int f_K, int64_t q, int bound) {
  T_small_prime q_u = static_cast<T_small_prime>(q);
  ModMatrix Emod = ModMatrixOfMatrix(E, q_u);
  int m = static_cast<int>(RankModPrime(Emod, q_u)) / (e_ram * f_res);
  HermitianLocalJordan hj;
  hj.type = type;
  hj.f_K = f_K;
  MyMatrix<T> Binv = Inverse(B);
  std::vector<int> len;
  MyMatrix<T> Mk = Binv * E;
  while (true) {
    len.push_back(ValuationIndexOverlattice(Mk, q));
    if (len.back() == 0 || len.size() > static_cast<size_t>(bound)) {
      break;
    }
    Mk = Mk * Pi;
  }
  std::vector<int> c;
  for (size_t k = 0; k + 1 < len.size(); k++) {
    c.push_back((len[k] - len[k + 1]) / f_res);
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
    hj.blocks.push_back({kv.first - delta, kv.second});
  }
  return hj;
}

/*
  The idempotents on X of the primes over q of a field given by the
  actions Omega of a Z-basis of its ring of integers, q unramified: those
  of several elements, lifted by E <- 3 E^2 - 2 E^3 modulo q^N, N =
  v_q(det B) + 6, and refined by their products until there are n_primes
  of them (fewer when the refinement does not get there).
 */
template <typename T>
inline std::vector<MyMatrix<T>>
PrimeIdempotentsByRefinement(MyMatrix<T> const &B,
                             std::vector<MyMatrix<T>> const &Omega,
                             int64_t q, int64_t n_primes) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int r = B.rows();
  int deg = static_cast<int>(Omega.size());
  int vdet = RationalValuation(DeterminantMat(B), q);
  int N = vdet + 6;
  Tint Q = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(q), N);
  Tint qz = UniversalScalarConversion<Tint, int64_t>(q);
  auto reduce_mod = [&](MyMatrix<T> M, Tint const &mod) {
    for (int i = 0; i < M.rows(); i++) {
      for (int j = 0; j < M.cols(); j++) {
        Tint num = GetNumerator_z(M(i, j)), den = GetDenominator_z(M(i, j));
        Tint inv = mod_inv(ResInt(den, mod), mod);
        Tint prod = num * inv;
        Tint v = ResInt(prod, mod);
        M(i, j) = UniversalScalarConversion<T, Tint>(v);
      }
    }
    return M;
  };
  auto lift_idempotent = [&](MyMatrix<T> E) {
    for (int it = 0; it < 64; it++) {
      MyMatrix<T> E2 = E * E;
      MyMatrix<T> En =
          reduce_mod(MyMatrix<T>(E2 * T(3) - E2 * E * T(2)), Q);
      if (En == E) {
        break;
      }
      E = En;
    }
    return E;
  };
  auto is_zero_mod_q = [&](MyMatrix<T> const &M) -> bool {
    MyMatrix<T> R = reduce_mod(M, qz);
    for (int i = 0; i < R.rows(); i++) {
      for (int j = 0; j < R.cols(); j++) {
        if (R(i, j) != 0) {
          return false;
        }
      }
    }
    return true;
  };
  std::vector<MyMatrix<T>> idem{IdentityMat<T>(r)};
  for (int i0 = 0; i0 < deg && static_cast<int64_t>(idem.size()) < n_primes; i0++) {
    for (int i1 = 0; i1 < deg && static_cast<int64_t>(idem.size()) < n_primes;
         i1++) {
      MyMatrix<T> S = Omega[i0] + Omega[i1] * T(i1 + 1);
      KrylovData<T> kd = KrylovMinimalPolynomial(S);
      std::vector<int64_t> f;
      for (auto const &x : kd.minpoly) {
        if (GetDenominator_z(x) != 1) {
          f.clear();
          break;
        }
        f.push_back(UniversalScalarConversion<int64_t, Tint>(GetNumerator_z(x)));
      }
      if (f.empty()) {
        continue;
      }
      PolyQ fq = PolyReduceSigned(f, q);
      std::vector<std::pair<PolyQ, int>> fac = FactorModQWithMultiplicity(fq, q);
      if (fac.size() <= 1) {
        continue;
      }
      std::vector<MyMatrix<T>> next;
      for (auto const &eps : idem) {
        for (auto const &[g, mult] : fac) {
          PolyQ P = {1};
          for (int k = 0; k < mult; k++) {
            P = PolyMulQ(P, g, q);
          }
          PolyQ Qp = PolyDivModQ(fq, P, q).first;
          PolyQ e0 = PolyMulQ(Qp, PolyInverseModQ(Qp, P, q), q);
          MyMatrix<T> Eg = lift_idempotent(
              reduce_mod(PolyOfMatrixQ(e0, S), Q));
          MyMatrix<T> prod = reduce_mod(MyMatrix<T>(eps * Eg), Q);
          if (!is_zero_mod_q(prod)) {
            next.push_back(lift_idempotent(prod));
          }
        }
      }
      idem = next;
    }
  }
  return idem;
}

/*
  The local factors at a prime q unramified in E / Q that divides every
  index [O_E : Z[theta]] (a common index divisor, as 2 in Q(sqrt 5,
  sqrt -7)): the idempotents of the g primes over q are obtained by
  refining those of several elements of O_E, each lifted on X by E <- 3 E^2
  - 2 E^3 modulo q^N; the uniformizer is q and the different is a unit.
 */
template <typename T>
inline std::optional<T>
UnramifiedFactorsByIdempotents(MyMatrix<T> const &B,
                               std::vector<MyMatrix<T>> const &Omega,
                               int64_t q, PrimeDecompositionData const &dE) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int r = B.rows();
  int vdet = RationalValuation(DeterminantMat(B), q);
  int N = vdet + 6;
  Tint Q = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(q), N);
  Tint qz = UniversalScalarConversion<Tint, int64_t>(q);
  auto reduce_mod = [&](MyMatrix<T> M, Tint const &mod) {
    for (int i = 0; i < M.rows(); i++) {
      for (int j = 0; j < M.cols(); j++) {
        Tint num = GetNumerator_z(M(i, j)), den = GetDenominator_z(M(i, j));
        Tint inv = mod_inv(ResInt(den, mod), mod);
        Tint prod = num * inv;
        Tint v = ResInt(prod, mod);
        M(i, j) = UniversalScalarConversion<T, Tint>(v);
      }
    }
    return M;
  };
  auto is_zero_mod_q = [&](MyMatrix<T> const &M) -> bool {
    MyMatrix<T> R = reduce_mod(M, qz);
    for (int i = 0; i < R.rows(); i++) {
      for (int j = 0; j < R.cols(); j++) {
        if (R(i, j) != 0) {
          return false;
        }
      }
    }
    return true;
  };
  std::vector<MyMatrix<T>> idem =
      PrimeIdempotentsByRefinement(B, Omega, q, dE.g);
  if (static_cast<int64_t>(idem.size()) != dE.g) {
    return {};
  }
  // the conjugation: E* = (B^{-1} E B)^T
  MyMatrix<T> Binv = Inverse(B);
  std::vector<int> conj_of(idem.size(), -1);
  for (size_t i = 0; i < idem.size(); i++) {
    MyMatrix<T> Es = reduce_mod((Binv * idem[i] * B).transpose(), Q);
    for (size_t j = 0; j < idem.size(); j++) {
      if (is_zero_mod_q(MyMatrix<T>(Es - idem[j]))) {
        conj_of[i] = static_cast<int>(j);
      }
    }
    if (conj_of[i] < 0) {
      return {};
    }
  }
  T res(1);
  std::set<size_t> done;
  MyMatrix<T> Pi =
      IdentityMat<T>(r) * UniversalScalarConversion<T, int64_t>(q);
  for (size_t i = 0; i < idem.size(); i++) {
    if (done.count(i)) {
      continue;
    }
    size_t j = static_cast<size_t>(conj_of[i]);
    done.insert(i);
    done.insert(j);
    int f_res = static_cast<int>(dE.f);
    int type = (j != i) ? 0 : 1;
    int f_K = (j != i) ? f_res : f_res / 2;
    res *= HermitianLocalFactorUnramified<T>(
        JordanFromIdempotent(B, idem[i], Pi, f_res, 1, 0, type, f_K, q,
                             4 * vdet + 8),
        q);
  }
  return res;
}

/*
  The mass of the genus of the Hermitian lattice over E given by the Gram
  matrix B of its trace form and the actions Omega of a Z-basis of O_E, or
  nothing in the cases not treated (see above).
 */
template <typename T>
inline std::optional<T>
HermitianMassAbelianLattice(MyMatrix<T> const &B,
                            std::vector<MyMatrix<T>> const &Omega,
                            AbelianCMField const &F) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int deg = static_cast<int>(AbelianFieldDegree(F));
  if (static_cast<int>(Omega.size()) != deg) {
    std::cerr << "ABHERM: " << Omega.size() << " basis elements for a field "
              << "of degree " << deg << "\n";
    throw TerminalException{1};
  }
  int r = B.rows();
  int m = r / deg;
  if (m == 1) {
    return RankOneHermitianMassAbelian<T>(F);
  }
  T res = StandardHermitianMassAbelian<T>(F, m);
  std::set<int64_t> primes;
  {
    Tint num = GetNumerator_z(DeterminantMat(B));
    if (num < 0) {
      num = -num;
    }
    if (num > 1) {
      for (Tint const &p : FactorsInt(num)) {
        primes.insert(UniversalScalarConversion<int64_t, Tint>(p));
      }
    }
    int64_t x = F.n;
    for (int64_t p = 2; p <= x; p++) {
      if (x % p == 0) {
        primes.insert(p);
        while (x % p == 0) {
          x /= p;
        }
      }
    }
  }
  DirichletGroup G = ComputeDirichletGroup(F.n);
  std::vector<std::vector<int64_t>> XK = FieldCharacters(G, F.H, true);
  std::vector<std::vector<int64_t>> XE = FieldCharacters(G, F.H, false);
  for (int64_t q : primes) {
    std::optional<LocalGenerator<T>> lg = FindLocalGenerator(B, Omega, q);
    if (!lg) {
      PrimeDecompositionData dE = PrimeDecompositionCharacters(G, XE, q);
      if (dE.e != 1) {
        return {};
      }
      std::optional<T> fac =
          UnramifiedFactorsByIdempotents(B, Omega, q, dE);
      if (!fac) {
        return {};
      }
      res *= *fac;
      continue;
    }
    int64_t eK = (F.n % q == 0) ? PrimeDecompositionCharacters(G, XK, q).e : 1;
    std::vector<AbelianPrime> aps = AbelianPrimesOverQ(*lg, q, eK);
    std::set<size_t> done;
    for (size_t i = 0; i < aps.size(); i++) {
      if (done.count(i)) {
        continue;
      }
      done.insert(i);
      done.insert(aps[i].conj_index);
      if (aps[i].type == 2) {
        if (q == 2) {
          // d = v_p(d(E / K)): the g primes over 2 are all ramified over
          // K, so v_2 |disc E| - 2 v_2 |disc K| = g f_K d
          int vE = DiscriminantValuationCharacters(G, XE, 2);
          int vK = DiscriminantValuationCharacters(G, XK, 2);
          int g_primes = static_cast<int>(aps.size());
          int dexp = (vE - 2 * vK) / (g_primes * aps[i].f_K);
          res *= AbelianLocalFactorDyadic(B, *lg, aps[i], dexp);
        } else {
          res *= AbelianLocalFactorRamified(B, *lg, aps[i], q);
        }
      } else {
        res *= HermitianLocalFactorUnramified<T>(
            AbelianJordanAtPrime(B, *lg, aps, i, q), q);
      }
    }
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_ABELIANHERMITIANMASS_H_
// clang-format on
