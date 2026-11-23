// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_PRIMEDECOMPOSITION_H_
#define SRC_GENERA_PRIMEDECOMPOSITION_H_

// clang-format off
#include "NumberField.h"
#include "PolynomialModQ.h"
#include "MAT_MatrixModPrime.h"
#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  The local arithmetic of a number field E at a rational prime q, E acting
  on a Z-lattice X through the matrices of a Z-basis of O_E:
  - a generator theta_q of E in O_E with q prime to the index
    [O_E : Z[theta_q]] (LocalGenerator), so that Z_q[theta_q] = O_E (x) Z_q;
  - the primes of E over q by Dedekind, the factors g of the minimal
    polynomial modulo q, with their ramification index, the conjugate
    factor and the type over the real subfield (AbelianPrime);
  - at a prime P = (q, g(theta_q)), whatever the other primes over q, the
    valuation v_P, the residues in F_q[x] / g and a uniformizer, through
    the idempotent of P (LocalPrimeK);
  - PrimeValuation, the valuation at a prime given either by the norm (P
    alone over q) or by LocalPrimeK;
  - at a common index divisor q (no theta_q), a generator theta_P for each
    prime P over q, regular at P (PrimeLocalGenerators).
 */

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
  The complex conjugate of theta (acting by S on X, of Krylov data kd) as
  a polynomial in theta: the adjoint S* = (B^{-1} S B)^T for the trace
  form B, on the rows, and v S* in the Krylov basis of v.
 */
template <typename T>
inline std::vector<T> ConjugateOfGenerator(MyMatrix<T> const &B,
                                           MyMatrix<T> const &S,
                                           KrylovData<T> const &kd) {
  MyMatrix<T> Sstar = (Inverse(B) * S * B).transpose();
  MyVector<T> v = GetMatrixRow(kd.basis, 0);
  MyVector<T> w = Sstar.transpose() * v;
  // w = sum_i a_i v S^i
  std::optional<MyVector<T>> sol = SolutionMat(kd.basis, w);
  if (!sol) {
    std::cerr << "PRIME_DECOMP: the conjugate of theta is not in Q(theta)\n";
    throw TerminalException{1};
  }
  int deg = kd.basis.rows();
  std::vector<T> conj(deg);
  for (int i = 0; i < deg; i++) {
    conj[i] = (*sol)(i);
  }
  return conj;
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
    MyMatrix<T> D = PowerTraceMatrix<T>(f);
    if (RationalValuation(DeterminantMat(D), q) != vdisc) {
      continue;
    }
    return LocalGenerator<T>{S, f, ConjugateOfGenerator(B, S, kd)};
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

template <typename T>
inline MyMatrix<T> PolyOfMatrixQ(PolyQ const &g, MyMatrix<T> const &S) {
  return PolyOfMatrixQ(g.c, S);
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
  NumberField<T> K;
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
  The valuation v_P of E at a prime P over q: v_q(N(a)) / div when P is
  the only prime over q (div its residue degree, or 1 when the
  ramification is counted in the norm, as at 1 - zeta in Q(zeta_d)), or
  read through the idempotent of P (LocalPrimeK) otherwise.
 */
template <typename T> struct PrimeValuation {
  int64_t q;
  int div = 1;
  // not owned, kept alive by the caller
  LocalPrimeK<T> *local = nullptr;
};

template <typename T>
inline int PrimeVal(ElemE<T> const &a, NumberField<T> const &K,
                    PrimeValuation<T> const &v) {
  if (EIsZero(a)) {
    return 1 << 28;
  }
  if (v.local) {
    return LocalVal(a, *v.local);
  }
  return ENormValuation(a, K, v.q) / v.div;
}

/*
  Dedekind's criterion: with f = prod g_i^{e_i} modulo q (f monic, the g_i
  monic) and F = (f - prod G_i^{e_i}) / q for lifts G_i of the g_i,
  Z_q[theta] is maximal at the prime (q, g(theta)) iff e = 1 or g does not
  divide F modulo q.
 */
template <typename Tint>
inline bool DedekindRegularAt(std::vector<int64_t> const &f,
                              std::vector<std::pair<PolyQ, int>> const &fac,
                              PolyQ const &g, int e, int64_t q) {
  if (e == 1) {
    return true;
  }
  PolyZ<Tint> prod{Tint(1)};
  for (auto const &[gi, ei] : fac) {
    PolyZ<Tint> gz = PolyZOfPolyQ<Tint>(gi);
    for (int k = 0; k < ei; k++) {
      prod = PolyZMul(prod, gz);
    }
  }
  Tint qz = UniversalScalarConversion<Tint, int64_t>(q);
  PolyQ Fq;
  for (size_t k = 0; k < f.size(); k++) {
    Tint c = UniversalScalarConversion<Tint, int64_t>(f[k]);
    if (k < prod.size()) {
      c -= prod[k];
    }
    // exact: f = prod modulo q
    Tint d = QuoInt(c, qz);
    Tint r = ResInt(d, qz);
    if (r < 0) {
      r += qz;
    }
    Fq.push_back(UniversalScalarConversion<int64_t, Tint>(r));
  }
  PolyTrim(Fq);
  if (Fq.empty()) {
    return false;
  }
  PolyQ rem = PolyDivModQ(Fq, g, q).second;
  PolyTrim(rem);
  return !rem.empty();
}

/*
  A common index divisor q: no generator of O_K has q prime to its index,
  but at a single prime P over q it suffices that the factor of f_theta
  modulo q belonging to P be g^e with Z_q[theta] maximal at P (Dedekind's
  criterion, automatic for e = 1), g coprime to the other factors: then the
  idempotent of P lies in Z_q[theta] and O_{K_P} = Z_q[theta] eps_P. So for
  each prime P a generator theta_P is searched among the small elements of
  O_K, P being recognized by its idempotent modulo q; the local arithmetic
  of P is that of LocalPrimeK for (theta_P, g, e), and the Gram matrix over
  K is taken on m vectors v_i with O_{K_P}-span eps_P X, that is with the
  theta^t eps_P v_i, t < e f, of F_q-rank m e f modulo q. Nothing when the
  g primes are not all found. For a totally real field the conjugate of
  theta is theta; for a CM field (cm) it is computed (ConjugateOfGenerator).
 */
template <typename T> struct PrimeLocalGenerator {
  LocalGenerator<T> lg;
  PolyQ g;
  int e;
  std::vector<MyVector<T>> basis;
};

template <typename T>
inline std::optional<std::vector<PrimeLocalGenerator<T>>>
PrimeLocalGenerators(MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Omega,
                     int64_t q, PrimeDecompositionData const &dK, bool cm) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int e_ram = static_cast<int>(dK.e);
  int deg = static_cast<int>(Omega.size());
  int r = B.rows();
  int m = r / deg;
  int f_res = static_cast<int>(dK.f);
  T_small_prime q_u = static_cast<T_small_prime>(q);
  auto red = [&](T const &x) -> T_residue {
    return ResidueModInt<T, Tint>(x, q_u);
  };
  // the small elements of O_K, by increasing size
  std::vector<std::vector<int>> cands;
  int n_tot = MyPow<int>(7, deg);
  for (int u = 0; u < n_tot; u++) {
    std::vector<int> c(deg);
    int t = u;
    for (int i = 0; i < deg; i++) {
      c[i] = t % 7 - 3;
      t /= 7;
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
  std::set<std::vector<int64_t>> found;
  std::vector<PrimeLocalGenerator<T>> res;
  for (auto const &c : cands) {
    if (static_cast<int64_t>(found.size()) == dK.g) {
      break;
    }
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
    // the conjugate of theta: theta itself for a totally real field
    std::vector<T> conj(deg, T(0));
    if (cm) {
      conj = ConjugateOfGenerator(B, S, kd);
    } else if (deg > 1) {
      conj[1] = T(1);
    }
    LocalGenerator<T> lg{S, f, conj};
    std::vector<std::pair<PolyQ, int>> fac =
        FactorModQWithMultiplicity(PolyReduceSigned(f, q), q);
    for (auto const &[g, mult] : fac) {
      if (mult != e_ram || static_cast<int>(g.size()) - 1 != f_res ||
          !DedekindRegularAt<Tint>(f, fac, g, mult, q)) {
        continue;
      }
      MyMatrix<T> E = AbelianPrimeIdempotent(lg, g, e_ram, q, 1);
      std::vector<int64_t> key;
      for (int i = 0; i < r; i++) {
        for (int j = 0; j < r; j++) {
          key.push_back(red(E(i, j)));
        }
      }
      if (found.count(key) == 1) {
        continue;
      }
      found.insert(key);
      // m vectors with O_{K_P}-span eps_P X
      int ef = e_ram * f_res;
      std::vector<MyMatrix<T>> Spow{IdentityMat<T>(r)};
      for (int t = 1; t < ef; t++) {
        Spow.push_back(Spow.back() * S);
      }
      std::vector<MyVector<T>> chosen;
      ModMatrix acc(0, r);
      auto try_add = [&](MyVector<T> const &v) -> void {
        ModMatrix test(acc.rows() + ef, r);
        test.topRows(acc.rows()) = acc;
        for (int t = 0; t < ef; t++) {
          MyVector<T> w = E.transpose() * (Spow[t].transpose() * v);
          for (int j = 0; j < r; j++) {
            test(acc.rows() + t, j) = red(w(j));
          }
        }
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
      if (static_cast<int>(chosen.size()) != m) {
        std::cerr << "PRIME_DECOMP: no local basis at a prime over the common index "
                     "divisor "
                  << q << "\n";
        throw TerminalException{1};
      }
      res.push_back({lg, g, e_ram, chosen});
    }
  }
  if (static_cast<int64_t>(found.size()) != dK.g) {
    return {};
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_PRIMEDECOMPOSITION_H_
// clang-format on
