// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ABELIANHERMITIANMASS_H_
#define SRC_GENERA_ABELIANHERMITIANMASS_H_

// clang-format off
#include "CyclotomicHermitianMass.h"
#include <map>
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
  E / Q being (f'(theta_q)) locally. A ramified prime of K is treated when
  it is the only prime of E over q and q is odd: the Jordan splitting over
  E with v_P = v_q(N) / f, the residue field F_q[x] / g of order q^f for
  the square classes, and the group orders of Hecke at q^{f_K}. The other
  ramified cases (dyadic, or several primes of E over q) leave the mass
  not computed.
 */

// The distinct irreducible factors of f modulo q with their multiplicities.
inline std::vector<std::pair<PolyQ, int>> FactorModQWithMultiplicity(PolyQ f,
                                                                      int64_t q) {
  f = PolyMonicQ(PolyReduce(f, q), q);
  std::vector<std::pair<PolyQ, int>> res;
  PolyQ rem = f;
  int deg = static_cast<int>(f.size()) - 1;
  PolyQ xpol{0, 1};
  for (int k = 1; k <= deg && rem.size() > 1; k++) {
    // x^{q^k} - x modulo rem: the product of the irreducible factors of
    // degree dividing k, the smaller ones being already removed
    PolyQ y = xpol;
    for (int i = 0; i < k; i++) {
      y = PolyPowModQ(y, static_cast<uint64_t>(q), rem, q);
    }
    PolyQ g = PolyGcdQ(PolySubQ(y, xpol, q), rem, q);
    PolyTrim(g);
    if (g.size() <= 1) {
      continue;
    }
    for (auto const &fac : EqualDegreeFactorization(PolyMonicQ(g, q), k, q)) {
      int mult = 0;
      while (true) {
        std::pair<PolyQ, PolyQ> qr = PolyDivModQ(rem, fac, q);
        PolyQ r = qr.second;
        PolyTrim(r);
        if (!r.empty()) {
          break;
        }
        rem = qr.first;
        mult++;
      }
      res.push_back({fac, mult});
    }
  }
  return res;
}

// Tr(theta^k), k = 0, ..., 2 deg, from the minimal polynomial by Newton.
inline std::vector<mpq_class> PowerTraces(std::vector<int64_t> const &f,
                                          int kmax) {
  int deg = static_cast<int>(f.size()) - 1;
  // f = x^deg + c_{deg-1} x^{deg-1} + ... + c_0
  std::vector<mpq_class> p(kmax + 1, mpq_class(0));
  p[0] = deg;
  for (int k = 1; k <= kmax; k++) {
    mpq_class s(0);
    // Newton: p_k + c_{deg-1} p_{k-1} + ... = -k c_{deg-k} for k <= deg
    for (int i = 1; i <= std::min(k - 1, deg); i++) {
      s += mpq_class(static_cast<long>(f[deg - i])) * p[k - i];
    }
    if (k <= deg) {
      s += mpq_class(static_cast<long>(k)) *
           mpq_class(static_cast<long>(f[deg - k]));
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
struct KrylovData {
  std::vector<mpq_class> minpoly;
  MyMatrix<mpq_class> basis;
};

inline KrylovData KrylovMinimalPolynomial(MyMatrix<mpq_class> const &S) {
  int r = S.rows();
  MyVector<mpq_class> v = ZeroVector<mpq_class>(r);
  v(0) = 1;
  std::vector<MyVector<mpq_class>> pw{v};
  for (int k = 1; k <= r; k++) {
    MyVector<mpq_class> w = S.transpose() * pw.back();
    MyMatrix<mpq_class> M(k + 1, r);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < r; j++) {
        M(i, j) = pw[i](j);
      }
    }
    for (int j = 0; j < r; j++) {
      M(k, j) = w(j);
    }
    MyMatrix<mpq_class> ker = NullspaceMat(M);
    if (ker.rows() > 0) {
      // sum_i ker(0, i) v S^i = 0, normalized at i = k
      std::vector<mpq_class> mp(k + 1);
      for (int i = 0; i <= k; i++) {
        mp[i] = ker(0, i) / ker(0, k);
      }
      MyMatrix<mpq_class> basis(k, r);
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
struct LocalGenerator {
  MyMatrix<mpq_class> S;
  std::vector<int64_t> f;
  std::vector<mpq_class> conj;
};

inline std::optional<LocalGenerator>
FindLocalGenerator(MyMatrix<mpq_class> const &B,
                   std::vector<MyMatrix<mpq_class>> const &Omega, int64_t q) {
  int deg = static_cast<int>(Omega.size());
  int r = B.rows();
  int m = r / deg;
  // disc O_E = det(Tr(omega_i omega_j)), Tr_{E/Q} = Tr on X / m
  MyMatrix<mpq_class> Tr(deg, deg);
  for (int i = 0; i < deg; i++) {
    for (int j = 0; j < deg; j++) {
      MyMatrix<mpq_class> P = Omega[i] * Omega[j];
      mpq_class t(0);
      for (int u = 0; u < r; u++) {
        t += P(u, u);
      }
      Tr(i, j) = t / mpq_class(static_cast<long>(m));
    }
  }
  int vdisc = ValuationQRational(DeterminantMat(Tr), q);
  MyMatrix<mpq_class> Binv = Inverse(B);
  // combinations with coefficients in [-2, 2], by increasing size
  std::vector<std::vector<int>> cands;
  int n_tot = 1;
  for (int i = 0; i < deg; i++) {
    n_tot *= 5;
  }
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
    MyMatrix<mpq_class> S = ZeroMatrix<mpq_class>(r, r);
    for (int i = 0; i < deg; i++) {
      if (c[i] != 0) {
        S += Omega[i] * mpq_class(c[i]);
      }
    }
    KrylovData kd = KrylovMinimalPolynomial(S);
    if (static_cast<int>(kd.minpoly.size()) != deg + 1) {
      continue;
    }
    std::vector<int64_t> f;
    bool integral = true;
    for (auto const &x : kd.minpoly) {
      if (x.get_den() != 1) {
        integral = false;
      }
      f.push_back(x.get_num().get_si());
    }
    if (!integral) {
      continue;
    }
    // disc f = det(Tr(theta^{i+j}))
    std::vector<mpq_class> pt = PowerTraces(f, 2 * deg);
    MyMatrix<mpq_class> D(deg, deg);
    for (int i = 0; i < deg; i++) {
      for (int j = 0; j < deg; j++) {
        D(i, j) = pt[i + j];
      }
    }
    if (ValuationQRational(DeterminantMat(D), q) != vdisc) {
      continue;
    }
    // the conjugate: S* = (B^{-1} S B)^T on the rows, v S* in the Krylov
    // basis of v
    MyMatrix<mpq_class> Sstar = (Binv * S * B).transpose();
    MyVector<mpq_class> v = GetMatrixRow(kd.basis, 0);
    MyVector<mpq_class> w = Sstar.transpose() * v;
    // w = sum_i a_i v S^i
    MyMatrix<mpq_class> M = kd.basis;
    std::optional<MyVector<mpq_class>> sol = SolutionMat(M, w);
    if (!sol) {
      std::cerr << "ABHERM: the conjugate of theta is not in Q(theta)\n";
      throw TerminalException{1};
    }
    std::vector<mpq_class> conj(deg);
    for (int i = 0; i < deg; i++) {
      conj[i] = (*sol)(i);
    }
    return LocalGenerator{S, f, conj};
  }
  return {};
}

inline mpz_class PowInt(int64_t q, int64_t e) {
  mpz_class r(1);
  for (int64_t i = 0; i < e; i++) {
    r *= static_cast<long>(q);
  }
  return r;
}

// The matrix of g(S) for a polynomial g with integer coefficients.
inline MyMatrix<mpq_class> PolyOfMatrixQ(PolyQ const &g, MyMatrix<mpq_class> const &S) {
  int r = S.rows();
  MyMatrix<mpq_class> res = ZeroMatrix<mpq_class>(r, r);
  for (size_t i = g.size(); i-- > 0;) {
    res = res * S;
    for (int j = 0; j < r; j++) {
      res(j, j) += mpq_class(static_cast<long>(g[i]));
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

inline std::vector<AbelianPrime>
AbelianPrimesOverQ(LocalGenerator const &lg, int64_t q, int64_t eK) {
  std::vector<std::pair<PolyQ, int>> fac =
      FactorModQWithMultiplicity(PolyQ(lg.f.begin(), lg.f.end()), q);
  // the conjugate of the prime of g: the factor dividing g(conj(x))
  auto compose_conj = [&](PolyQ const &g) -> PolyQ {
    // g(c(x)) modulo f and q, c the conjugate of theta (q-integral)
    mpz_class qz(static_cast<long>(q));
    PolyQ c;
    for (auto const &x : lg.conj) {
      mpz_class num = x.get_num() % qz;
      mpz_class den = x.get_den() % qz;
      if (num < 0) {
        num += qz;
      }
      if (den < 0) {
        den += qz;
      }
      int64_t v = MulModQ(num.get_si(), InvModQ(den.get_si(), q), q);
      c.push_back(v);
    }
    PolyQ fq = PolyReduce(PolyQ(lg.f.begin(), lg.f.end()), q);
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
  The Jordan data at a prime of E over q unramified over K (split: one of
  the pair), as HermitianJordanAtPrime does for Q(zeta_d).
 */
inline HermitianLocalJordan
AbelianJordanAtPrime(MyMatrix<mpq_class> const &B, LocalGenerator const &lg,
                     std::vector<AbelianPrime> const &primes, size_t i_prime,
                     int64_t q) {
  int r = B.rows();
  AbelianPrime const &ap = primes[i_prime];
  int f_res = static_cast<int>(ap.g.size()) - 1;
  int vdet = ValuationQRational(DeterminantMat(B), q);
  std::vector<int64_t> f = lg.f;
  // the precision covers the discriminant group and the different
  int vdiscf = 0;
  {
    std::vector<mpq_class> pt = PowerTraces(f, 2 * (static_cast<int>(f.size()) - 1));
    int deg = static_cast<int>(f.size()) - 1;
    MyMatrix<mpq_class> D(deg, deg);
    for (int a = 0; a < deg; a++) {
      for (int b = 0; b < deg; b++) {
        D(a, b) = pt[a + b];
      }
    }
    vdiscf = ValuationQRational(DeterminantMat(D), q);
  }
  int N = vdet + vdiscf + 6;
  mpz_class Q = PowInt(q, N);
  PolyQ fq = PolyReduce(PolyQ(f.begin(), f.end()), q);
  // the idempotent of the prime, lifted modulo q^N
  PolyQ P = {1};
  for (int i = 0; i < ap.e; i++) {
    P = PolyMulQ(P, ap.g, q);
  }
  PolyQ Qp = PolyDivModQ(fq, P, q).first;
  PolyQ e0 = PolyMulQ(Qp, PolyInverseModQ(Qp, P, q), q);
  PolyZ ez;
  for (int64_t c : e0) {
    ez.push_back(mpz_class(static_cast<long>(c)));
  }
  ez = PolyZReduce(ez, f, Q);
  for (int it = 0; it < 64; it++) {
    PolyZ e2 = PolyZReduce(PolyZMul(ez, ez), f, Q);
    PolyZ e3 = PolyZReduce(PolyZMul(e2, ez), f, Q);
    PolyZ en(e2.size());
    for (size_t k = 0; k < e2.size(); k++) {
      en[k] = 3 * e2[k] - 2 * e3[k];
    }
    en = PolyZReduce(en, f, Q);
    if (en == ez) {
      break;
    }
    ez = en;
  }
  MyMatrix<mpq_class> Id = IdentityMat<mpq_class>(r);
  MyMatrix<mpq_class> E = ZeroMatrix<mpq_class>(r, r);
  {
    MyMatrix<mpq_class> pw = Id;
    for (size_t k = 0; k < ez.size(); k++) {
      E += pw * mpq_class(ez[k]);
      pw = pw * lg.S;
    }
  }
  // the uniformizer: q, or g(theta) when the prime is ramified over q
  MyMatrix<mpq_class> Pi =
      (ap.e == 1) ? MyMatrix<mpq_class>(Id * mpq_class(static_cast<long>(q)))
                  : PolyOfMatrixQ(ap.g, lg.S);
  // delta = v_P(f'(theta)): v_q(det(f'(S) E + 1 - E)) on X is m f v_P
  int deg = static_cast<int>(f.size()) - 1;
  int m_total = r / deg;
  PolyQ fprime;
  for (size_t k = 1; k < f.size(); k++) {
    fprime.push_back(static_cast<int64_t>(k) * f[k]);
  }
  MyMatrix<mpq_class> FP = PolyOfMatrixQ(fprime, lg.S);
  // the O_P-rank m
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
  int m = rk / (ap.e * f_res);
  int delta = 0;
  {
    MyMatrix<mpq_class> M = FP * E + Id - E;
    int vv = ValuationQRational(DeterminantMat(M), q);
    // on X the prime part has Z_q-rank m e f, and det = N(f')^m locally
    delta = vv / (m * f_res);
  }
  (void)m_total;
  HermitianLocalJordan hj;
  hj.type = ap.type;
  hj.f_K = ap.f_K;
  MyMatrix<mpq_class> Binv = Inverse(B);
  std::vector<int> len;
  MyMatrix<mpq_class> Mk = Binv * E;
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
  The local factor at a ramified prime of K, P alone over q, q odd: the
  formula of HermitianLocalFactorRamified with the residue field of order
  q^f.
 */
struct AbelianLocalGram {
  CyclotomicField K;
  int m;
  std::vector<std::vector<ElemE>> H;
};

// The Hermitian Gram matrix over E of an O_P-basis of X at the ramified
// prime P alone over q, chosen among the standard vectors.
inline AbelianLocalGram AbelianLocalHermitianGram(MyMatrix<mpq_class> const &B,
                                                  LocalGenerator const &lg,
                                                  AbelianPrime const &ap,
                                                  int64_t q) {
  int r = B.rows();
  int deg = static_cast<int>(lg.f.size()) - 1;
  int m = r / deg;
  int f_res = static_cast<int>(ap.g.size()) - 1;
  CyclotomicField K;
  K.d = 0;
  K.deg = deg;
  K.Phi = lg.f;
  K.conj = lg.conj;
  K.val_div = f_res;
  // the traces of the powers of theta, and the dual basis
  std::vector<mpq_class> pt = PowerTraces(lg.f, 2 * deg);
  MyMatrix<mpq_class> T(deg, deg);
  for (int s2 = 0; s2 < deg; s2++) {
    for (int t = 0; t < deg; t++) {
      T(s2, t) = pt[s2 + t];
    }
  }
  MyMatrix<mpq_class> Tinv = Inverse(T);
  // an O_P-basis among the standard vectors: their O_P / P-spans
  // independent modulo pi X, the residue field having degree f_res
  MyMatrix<mpq_class> Pi = PolyOfMatrixQ(ap.g, lg.S);
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
  std::vector<std::vector<int64_t>> rows;
  for (int i = 0; i < r; i++) {
    std::vector<int64_t> row(r);
    for (int j = 0; j < r; j++) {
      row[j] = red(Pi(i, j));
    }
    rows.push_back(row);
  }
  int base_rank = RankModQ(rows, q);
  std::vector<MyMatrix<mpq_class>> Spow{IdentityMat<mpq_class>(r)};
  for (int t = 1; t < deg; t++) {
    Spow.push_back(Spow.back() * lg.S);
  }
  std::vector<int> chosen;
  std::vector<std::vector<int64_t>> acc = rows;
  for (int k = 0; k < r && static_cast<int>(chosen.size()) < m; k++) {
    std::vector<std::vector<int64_t>> test = acc;
    for (int t = 0; t < f_res; t++) {
      std::vector<int64_t> row(r);
      for (int j = 0; j < r; j++) {
        row[j] = red(Spow[t](k, j));
      }
      test.push_back(row);
    }
    if (RankModQ(test, q) ==
        base_rank + (static_cast<int>(chosen.size()) + 1) * f_res) {
      chosen.push_back(k);
      acc = test;
    }
  }
  if (static_cast<int>(chosen.size()) != m) {
    std::cerr << "ABHERM: no O_P-basis found\n";
    throw TerminalException{1};
  }
  std::vector<std::vector<ElemE>> H(m, std::vector<ElemE>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      ElemE val(deg, mpq_class(0));
      for (int t = 0; t < deg; t++) {
        mpq_class bt(0);
        for (int u = 0; u < r; u++) {
          bt += Spow[t](chosen[i], u) * B(u, chosen[j]);
        }
        for (int s2 = 0; s2 < deg; s2++) {
          val[s2] += bt * Tinv(t, s2);
        }
      }
      H[i][j] = val;
    }
  }
  return {K, m, H};
}

inline mpq_class AbelianLocalFactorRamified(MyMatrix<mpq_class> const &B,
                                            LocalGenerator const &lg,
                                            AbelianPrime const &ap, int64_t q) {
  AbelianLocalGram alg = AbelianLocalHermitianGram(B, lg, ap, q);
  CyclotomicField const &K = alg.K;
  int m = alg.m;
  int deg = K.deg;
  std::vector<std::vector<ElemE>> const &H = alg.H;
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
  std::map<int, std::pair<int, ElemE>> blocks;
  for (auto &kv : HermitianJordanSplit(H, K, q)) {
    int rk = 0;
    ElemE det = EZetaPower(0, K);
    for (auto &G : kv.second) {
      rk += static_cast<int>(G.size());
      det = EMul(det, HermitianSmallDet(G, K), K);
    }
    blocks[kv.first] = {rk, det};
  }
  int64_t qK = 1;
  for (int i = 0; i < ap.f_K; i++) {
    qK *= q;
  }
  mpq_class fac = HeckeGroupOrderSp(2 * (m / 2), qK);
  int64_t Nexp = 0;
  int after_total = m;
  // pi = g(theta), and pi pi-bar
  ElemE pi(deg, mpq_class(0));
  {
    ElemE pw = EZetaPower(0, K);
    for (size_t i = 0; i < ap.g.size(); i++) {
      for (int t = 0; t < deg; t++) {
        pi[t] += mpq_class(static_cast<long>(ap.g[i])) * pw[t];
      }
      pw = EMul(pw, EZetaPower(1, K), K);
    }
  }
  ElemE pipibar = EMul(pi, EConj(pi, K), K);
  ElemE pinv = EInv(pipibar, K);
  // the square class of a P-unit in the residue field F_q[x] / g
  auto is_square = [&](ElemE const &u) -> int {
    PolyQ uq;
    for (auto const &x : u) {
      uq.push_back(red(x));
    }
    PolyQ ur = PolyDivModQ(uq, ap.g, q).second;
    PolyTrim(ur);
    if (ur.empty()) {
      return 0;
    }
    PolyQ y = PolyPowModQ(ur, static_cast<uint64_t>((qK - 1) / 2), ap.g, q);
    PolyTrim(y);
    return (y.size() == 1 && y[0] == 1) ? 1 : -1;
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
      fac /= HeckeGroupOrderSp(mi, qK);
    } else {
      ElemE u = kv.second.second;
      int w = EValuation(u, K, q);
      for (int i = 0; i < w / 2; i++) {
        u = EMul(u, pinv, K);
      }
      if ((mi * (mi - 1) / 2) % 2 == 1) {
        for (auto &x : u) {
          x = -x;
        }
      }
      fac /= HeckeGroupOrderO(mi, is_square(u), qK);
    }
    Nexp -= static_cast<int64_t>(s) * mi * ri;
    Nexp += static_cast<int64_t>(s) * mi * m / 2;
  }
  if (m % 2 == 0) {
    Nexp += m / 2;
  }
  mpq_class qpow(1);
  for (int64_t i = 0; i < (Nexp >= 0 ? Nexp : -Nexp); i++) {
    if (Nexp >= 0) {
      qpow *= mpq_class(static_cast<long>(qK));
    } else {
      qpow /= mpq_class(static_cast<long>(qK));
    }
  }
  return qpow * fac;
}

/*
  The dyadic ramified prime P of E over the prime p of K, P alone over 2,
  following HermitianLocalFactorDyadic with what was particular to Q(zeta_{2^a})
  made general: the exponent d of the relative different of E / K at P
  (2 there), which is also the conductor exponent of E_P / K_p, so that
  v_p(Tr(c O_E)) = floor((v_P(c) + d) / 2), a unit u of K is a norm iff
  u = y y-bar modulo p^d for a unit y of O_P / P^{2d}, and the normic
  defect of a non-norm x is v_p(x) + d - 1; and the residue field
  F_P = F_2[x] / g of order 2^f, whose elements are the sum b_j theta^j,
  j < f, b_j in {0, 1}, over which the lines of L / pi L are taken.
 */
struct GenericDyadic {
  CyclotomicField K;
  ElemE one;
  ElemE pi;
  ElemE pi_inv;
  ElemE pipibar;
  ElemE pipibar_inv;
  int d;
  int f;
  // the representatives of F_P, the zero first
  std::vector<ElemE> reps;
};

inline int GDVal(ElemE const &x, GenericDyadic const &C) {
  return EValuation(x, C.K, 2);
}

inline int GDTraceVal(ElemE const &c, GenericDyadic const &C) {
  int v = GDVal(c, C);
  if (v >= (1 << 27)) {
    return v;
  }
  int num = v + C.d;
  int vp = (num >= 0) ? num / 2 : -((-num + 1) / 2);
  return 2 * vp;
}

inline int GDNormVal(std::vector<std::vector<ElemE>> const &G,
                     GenericDyadic const &C) {
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

inline bool GDIsLocalNorm(ElemE x, GenericDyadic const &C) {
  CyclotomicField const &K = C.K;
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
  std::vector<ElemE> pipow{C.one};
  for (int i = 1; i < len; i++) {
    pipow.push_back(EMul(pipow.back(), C.pi, K));
  }
  std::vector<size_t> idx(len, 0);
  idx[0] = 1;
  while (true) {
    ElemE y(K.deg, mpq_class(0));
    for (int i = 0; i < len; i++) {
      if (idx[i] == 0) {
        continue;
      }
      ElemE t = EMul(C.reps[idx[i]], pipow[i], K);
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
inline std::vector<std::vector<ElemE>>
GDOverlattice(std::vector<std::vector<ElemE>> const &H,
              std::vector<size_t> const &a, GenericDyadic const &C) {
  int m = static_cast<int>(H.size());
  ElemE zero(C.K.deg, mpq_class(0));
  std::vector<std::vector<ElemE>> U(m, std::vector<ElemE>(m, zero));
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
inline std::vector<std::vector<ElemE>>
GDSublattice(std::vector<std::vector<ElemE>> const &H,
             std::vector<size_t> const &a, GenericDyadic const &C) {
  int m = static_cast<int>(H.size());
  ElemE zero(C.K.deg, mpq_class(0));
  std::vector<std::vector<ElemE>> U(m, std::vector<ElemE>(m, zero));
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

struct GDSymbolEntry {
  int rank;
  int scale;
  bool normal;
  int v;
  ElemE det;
};

inline std::vector<GDSymbolEntry>
GDKirschmerSymbol(std::vector<std::vector<ElemE>> const &H,
                  GenericDyadic const &C) {
  std::map<int, std::vector<SmallGram>> blocks = HermitianJordanSplit(H, C.K, 2);
  std::vector<int> scales, nvals, ranks;
  std::vector<ElemE> dets;
  for (auto &kv : blocks) {
    int nv = 1 << 28;
    int rk = 0;
    ElemE det = C.one;
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
  std::vector<GDSymbolEntry> sym;
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

inline bool GDIsLocallyIsometric(std::vector<GDSymbolEntry> const &S1,
                                 std::vector<GDSymbolEntry> const &S2,
                                 GenericDyadic const &C) {
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
  if (!GDIsLocalNorm(ratio(t - 1), C)) {
    return false;
  }
  for (size_t i = 0; i + 1 < t; i++) {
    ElemE x = ratio(i);
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
  The local factor at the dyadic ramified prime P alone over 2, of
  relative different exponent d: the chain method of
  HermitianLocalFactorDyadic over the residue field of P.
 */
inline mpq_class AbelianLocalFactorDyadic(MyMatrix<mpq_class> const &B,
                                          LocalGenerator const &lg,
                                          AbelianPrime const &ap, int d) {
  AbelianLocalGram alg = AbelianLocalHermitianGram(B, lg, ap, 2);
  GenericDyadic C;
  C.K = alg.K;
  CyclotomicField const &K = C.K;
  int deg = K.deg;
  int m = alg.m;
  C.one = EZetaPower(0, K);
  C.pi = ElemE(deg, mpq_class(0));
  {
    ElemE pw = C.one;
    for (size_t i = 0; i < ap.g.size(); i++) {
      for (int t = 0; t < deg; t++) {
        C.pi[t] += mpq_class(static_cast<long>(ap.g[i])) * pw[t];
      }
      pw = EMul(pw, EZetaPower(1, K), K);
    }
  }
  C.pi_inv = EInv(C.pi, K);
  C.pipibar = EMul(C.pi, EConj(C.pi, K), K);
  C.pipibar_inv = EInv(C.pipibar, K);
  C.d = d;
  C.f = static_cast<int>(ap.g.size()) - 1;
  for (int64_t u = 0; u < (int64_t(1) << C.f); u++) {
    ElemE x(deg, mpq_class(0));
    for (int j = 0; j < C.f; j++) {
      if (u & (int64_t(1) << j)) {
        ElemE t = EZetaPower(j, K);
        for (int k = 0; k < deg; k++) {
          x[k] += t[k];
        }
      }
    }
    C.reps.push_back(x);
  }
  int64_t qK = int64_t(1) << ap.f_K;
  std::vector<std::vector<size_t>> lines = GDLines(m, C.reps.size());
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
    if (GDIsLocalNorm(det, C)) {
      return mpq_class(1);
    }
    mpq_class qm(1);
    for (int i = 0; i < m; i++) {
      qm *= mpq_class(static_cast<long>(qK));
    }
    return (qm - 1) / mpq_class(static_cast<long>(2 * qK + 2));
  };
  int nL = GDNormVal(alg.H, C);
  auto find_over = [&](std::vector<std::vector<ElemE>> const &H)
      -> std::optional<std::vector<std::vector<ElemE>>> {
    for (auto const &a : lines) {
      std::vector<std::vector<ElemE>> Hn = GDOverlattice(H, a, C);
      if (GDNormVal(Hn, C) >= nL) {
        return Hn;
      }
    }
    return {};
  };
  std::vector<std::vector<std::vector<ElemE>>> chain{alg.H};
  while (true) {
    std::optional<std::vector<std::vector<ElemE>>> opt = find_over(chain.back());
    if (!opt) {
      break;
    }
    chain.push_back(*opt);
    if (chain.size() > 64) {
      std::cerr << "ABHERM: the chain of overlattices does not end\n";
      throw TerminalException{1};
    }
  }
  mpq_class fac = maximal_factor(alg.H);
  for (size_t i = 0; i + 1 < chain.size(); i++) {
    std::vector<GDSymbolEntry> sym_i = GDKirschmerSymbol(chain[i], C);
    std::vector<GDSymbolEntry> sym_n = GDKirschmerSymbol(chain[i + 1], C);
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
      std::cerr << "ABHERM: a step of the chain has no isometric neighbour\n";
      throw TerminalException{1};
    }
    fac *= mpq_class(static_cast<long>(A), static_cast<long>(Bc));
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
inline HermitianLocalJordan
JordanFromIdempotent(MyMatrix<mpq_class> const &B, MyMatrix<mpq_class> const &E,
                     MyMatrix<mpq_class> const &Pi, int f_res, int e_ram,
                     int delta, int type, int f_K, int64_t q, int bound) {
  int r = B.rows();
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
  int m = RankModQ(Emod, q) / (e_ram * f_res);
  HermitianLocalJordan hj;
  hj.type = type;
  hj.f_K = f_K;
  MyMatrix<mpq_class> Binv = Inverse(B);
  std::vector<int> len;
  MyMatrix<mpq_class> Mk = Binv * E;
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
  The local factors at a prime q unramified in E / Q that divides every
  index [O_E : Z[theta]] (a common index divisor, as 2 in Q(sqrt 5,
  sqrt -7)): the idempotents of the g primes over q are obtained by
  refining those of several elements of O_E, each lifted on X by E <- 3 E^2
  - 2 E^3 modulo q^N; the uniformizer is q and the different is a unit.
 */
inline std::optional<mpq_class>
UnramifiedFactorsByIdempotents(MyMatrix<mpq_class> const &B,
                               std::vector<MyMatrix<mpq_class>> const &Omega,
                               int64_t q, PrimeDecompositionData const &dE) {
  int r = B.rows();
  int deg = static_cast<int>(Omega.size());
  int vdet = ValuationQRational(DeterminantMat(B), q);
  int N = vdet + 6;
  mpz_class Q = PowInt(q, N);
  mpz_class qz(static_cast<long>(q));
  auto reduce_mod = [&](MyMatrix<mpq_class> M, mpz_class const &mod) {
    for (int i = 0; i < M.rows(); i++) {
      for (int j = 0; j < M.cols(); j++) {
        mpz_class num = M(i, j).get_num(), den = M(i, j).get_den();
        mpz_class inv;
        mpz_invert(inv.get_mpz_t(), den.get_mpz_t(), mod.get_mpz_t());
        mpz_class v = (num * inv) % mod;
        if (v < 0) {
          v += mod;
        }
        M(i, j) = mpq_class(v);
      }
    }
    return M;
  };
  auto lift_idempotent = [&](MyMatrix<mpq_class> E) {
    for (int it = 0; it < 64; it++) {
      MyMatrix<mpq_class> E2 = E * E;
      MyMatrix<mpq_class> En =
          reduce_mod(MyMatrix<mpq_class>(E2 * mpq_class(3) - E2 * E * mpq_class(2)), Q);
      if (En == E) {
        break;
      }
      E = En;
    }
    return E;
  };
  auto is_zero_mod_q = [&](MyMatrix<mpq_class> const &M) -> bool {
    MyMatrix<mpq_class> R = reduce_mod(M, qz);
    for (int i = 0; i < R.rows(); i++) {
      for (int j = 0; j < R.cols(); j++) {
        if (R(i, j) != 0) {
          return false;
        }
      }
    }
    return true;
  };
  std::vector<MyMatrix<mpq_class>> idem{IdentityMat<mpq_class>(r)};
  for (int i0 = 0; i0 < deg && static_cast<int64_t>(idem.size()) < dE.g; i0++) {
    for (int i1 = 0; i1 < deg && static_cast<int64_t>(idem.size()) < dE.g;
         i1++) {
      MyMatrix<mpq_class> S = Omega[i0] + Omega[i1] * mpq_class(i1 + 1);
      KrylovData kd = KrylovMinimalPolynomial(S);
      std::vector<int64_t> f;
      for (auto const &x : kd.minpoly) {
        if (x.get_den() != 1) {
          f.clear();
          break;
        }
        f.push_back(x.get_num().get_si());
      }
      if (f.empty()) {
        continue;
      }
      PolyQ fq = PolyReduce(PolyQ(f.begin(), f.end()), q);
      std::vector<std::pair<PolyQ, int>> fac = FactorModQWithMultiplicity(fq, q);
      if (fac.size() <= 1) {
        continue;
      }
      std::vector<MyMatrix<mpq_class>> next;
      for (auto const &eps : idem) {
        for (auto const &[g, mult] : fac) {
          PolyQ P = {1};
          for (int k = 0; k < mult; k++) {
            P = PolyMulQ(P, g, q);
          }
          PolyQ Qp = PolyDivModQ(fq, P, q).first;
          PolyQ e0 = PolyMulQ(Qp, PolyInverseModQ(Qp, P, q), q);
          MyMatrix<mpq_class> Eg = lift_idempotent(
              reduce_mod(PolyOfMatrixQ(e0, S), Q));
          MyMatrix<mpq_class> prod = reduce_mod(MyMatrix<mpq_class>(eps * Eg), Q);
          if (!is_zero_mod_q(prod)) {
            next.push_back(lift_idempotent(prod));
          }
        }
      }
      idem = next;
    }
  }
  if (static_cast<int64_t>(idem.size()) != dE.g) {
    return {};
  }
  // the conjugation: E* = (B^{-1} E B)^T
  MyMatrix<mpq_class> Binv = Inverse(B);
  std::vector<int> conj_of(idem.size(), -1);
  for (size_t i = 0; i < idem.size(); i++) {
    MyMatrix<mpq_class> Es = reduce_mod((Binv * idem[i] * B).transpose(), Q);
    for (size_t j = 0; j < idem.size(); j++) {
      if (is_zero_mod_q(MyMatrix<mpq_class>(Es - idem[j]))) {
        conj_of[i] = static_cast<int>(j);
      }
    }
    if (conj_of[i] < 0) {
      return {};
    }
  }
  mpq_class res(1);
  std::set<size_t> done;
  MyMatrix<mpq_class> Pi =
      IdentityMat<mpq_class>(r) * mpq_class(static_cast<long>(q));
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
    res *= HermitianLocalFactorUnramified(
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
inline std::optional<mpq_class>
HermitianMassAbelianLattice(MyMatrix<mpq_class> const &B,
                            std::vector<MyMatrix<mpq_class>> const &Omega,
                            AbelianCMField const &F) {
  int deg = static_cast<int>(AbelianFieldDegree(F));
  if (static_cast<int>(Omega.size()) != deg) {
    std::cerr << "ABHERM: " << Omega.size() << " basis elements for a field "
              << "of degree " << deg << "\n";
    throw TerminalException{1};
  }
  int r = B.rows();
  int m = r / deg;
  if (m == 1) {
    return RankOneHermitianMassAbelian(F);
  }
  mpq_class res = StandardHermitianMassAbelian(F, m);
  std::set<int64_t> primes;
  {
    mpz_class num = DeterminantMat(B).get_num();
    if (num < 0) {
      num = -num;
    }
    if (num > 1) {
      for (mpz_class const &p : FactorsInt(num)) {
        primes.insert(p.get_si());
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
    std::optional<LocalGenerator> lg = FindLocalGenerator(B, Omega, q);
    if (!lg) {
      PrimeDecompositionData dE = PrimeDecompositionCharacters(G, XE, q);
      if (dE.e != 1) {
        return {};
      }
      std::optional<mpq_class> fac =
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
        if (aps.size() > 1) {
          return {};
        }
        if (q == 2) {
          // d = v_p(d(E / K)) = (v_2 |disc E| - 2 v_2 |disc K|) / f_K, p
          // alone over 2
          int vE = DiscriminantValuationCharacters(G, XE, 2);
          int vK = DiscriminantValuationCharacters(G, XK, 2);
          int dexp = (vE - 2 * vK) / aps[i].f_K;
          res *= AbelianLocalFactorDyadic(B, *lg, aps[i], dexp);
        } else {
          res *= AbelianLocalFactorRamified(B, *lg, aps[i], q);
        }
      } else {
        res *= HermitianLocalFactorUnramified(
            AbelianJordanAtPrime(B, *lg, aps, i, q), q);
      }
    }
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_ABELIANHERMITIANMASS_H_
// clang-format on
