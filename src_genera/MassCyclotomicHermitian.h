// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSCYCLOTOMICHERMITIAN_H_
#define SRC_GENERA_MASSCYCLOTOMICHERMITIAN_H_

// clang-format off
#include "CyclotomicArithmetic.h"
#include "ScalarNumberField.h"
#include "ScalarAbelianField.h"
#include "PrimeDecomposition.h"
#include "LocalJordan.h"
#include "ClassicalGroupOrders.h"
#include "ScalarPolynomialModQ.h"
#include "GeneraArithmetic.h"
#include "SiegelSpecialValues.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include "NumberTheoryModular.h"
#include "PolynomialRealRoots.h"
#include "MAT_MatrixModPrime.h"
#include <functional>
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


// The standard mass of a Hermitian lattice of rank m over Q(zeta_d), in
// absolute value.
template <typename T>
inline T StandardHermitianMassCyclotomic(int64_t d_in, int m) {
  int64_t d = NormalizeCyclotomicIndex(d_in);
  int64_t n = EulerPhi(d) / 2;
  T res(1);
  int64_t ex = 1 - n * m;
  T two(2);
  for (int64_t i = 0; i < (ex >= 0 ? ex : -ex); i++) {
    if (ex >= 0) {
      res *= two;
    } else {
      res /= two;
    }
  }
  for (int i = 1; i <= m; i++) {
    res *= ProductDirichletLValues<T>(d, i, i % 2);
  }
  return (res < 0) ? T(-res) : res;
}

// The mass of the genus of a Hermitian lattice of rank one over
// Q(zeta_d), which depends only on d.
template <typename T>
inline T RankOneHermitianMassCyclotomic(int64_t d_in) {
  int64_t d = NormalizeCyclotomicIndex(d_in);
  T res = StandardHermitianMassCyclotomic<T>(d, 1);
  if (PrimePowerDecomposition(d)) {
    res /= 2;
  }
  return res;
}

/*
  The standard mass over an abelian CM field (ScalarAbelianField.h): its
  Dirichlet characters are those modulo n trivial on H, the even ones being
  those of K, so the formula above holds with the characters restricted:
  zeta_K(s) and L(s, chi_{E/K}) are the products of the L(s, chi) over the
  even and the odd characters of E. The cyclotomic case is H = {1}.
 */
template <typename T>
inline T StandardHermitianMassAbelian(AbelianField const &F, int m) {
  int64_t n = AbelianFieldDegree(F) / 2;
  T res(1);
  int64_t ex = 1 - n * m;
  T two(2);
  for (int64_t i = 0; i < (ex >= 0 ? ex : -ex); i++) {
    if (ex >= 0) {
      res *= two;
    } else {
      res /= two;
    }
  }
  for (int i = 1; i <= m; i++) {
    res *= ProductDirichletLValues<T>(F.n, i, i % 2, F.H);
  }
  return (res < 0) ? T(-res) : res;
}

/*
  The finite primes of K ramified in E: over a rational prime q they all
  ramify when e_q(E) = 2 e_q(K), and there are g_q(K) of them. Returned as
  the pairs (q, number of primes of K over q ramified in E).
 */
inline std::vector<std::pair<int64_t, int64_t>>
RamifiedPrimesCM(AbelianField const &F) {
  std::vector<std::pair<int64_t, int64_t>> res;
  int64_t x = F.n;
  for (int64_t q = 2; q <= x; q++) {
    if (x % q != 0) {
      continue;
    }
    while (x % q == 0) {
      x /= q;
    }
    CMPrimeDecomposition dec = DecomposePrimeCM(F, q);
    if (dec.type == 2) {
      res.push_back({q, dec.K.g});
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
template <typename T>
inline T RankOneHermitianMassAbelian(AbelianField const &F) {
  T res = StandardHermitianMassAbelian<T>(F, 1);
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
  even scale (HermitianLocalFactorRamified).
 */


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

// v_q of the index [L + span(rows of M) : L], M rational, L = Z^r.
template <typename T>
inline int ValuationIndexOverlattice(MyMatrix<T> const &M, int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int r = M.cols();
  Tint D(1);
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < r; j++) {
      D = LCMpair(D, GetDenominator_z(M(i, j)));
    }
  }
  MyMatrix<Tint> G(M.rows() + r, r);
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < r; j++) {
      T v = M(i, j) * D;
      G(i, j) = GetNumerator_z(v);
    }
  }
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < r; j++) {
      G(M.rows() + i, j) = (i == j) ? D : Tint(0);
    }
  }
  MyMatrix<Tint> Bz = GetZbasis(G);
  T idx = T(Tint(1)) / T(DeterminantMat(Bz));
  idx *= MyPow<T>(T(D), r);
  return RationalValuation(idx, q);
}

template <typename T>
inline std::vector<HermitianLocalJordan>
HermitianJordanAtPrime(MyMatrix<T> const &B, MyMatrix<T> const &S,
                       int64_t d, int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
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
  PolyQ phid_q = PolyReduceSigned(PhiD, q);
  std::vector<PolyQ> gs =
      EqualDegreeFactorization(PolyReduceSigned(PhiDp, q), f);
  // the precision: beyond the exponent of the discriminant group at q
  int vdet = RationalValuation(DeterminantMat(B), q);
  int N = vdet + 4;
  Tint Q = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(q), N);
  // the dual basis, the uniformizer and the idempotent of each factor
  MyMatrix<T> Binv = Inverse(B);
  MyMatrix<T> Id = IdentityMat<T>(r);
  MyMatrix<T> Pi;
  if (qa == 1) {
    Pi = Id * UniversalScalarConversion<T, int64_t>(q);
  } else {
    MyMatrix<T> Sd = Id;
    for (int64_t i = 0; i < dp; i++) {
      Sd *= S;
    }
    Pi = Id - Sd;
  }
  auto idempotent_matrix = [&](PolyQ const &g) -> MyMatrix<T> {
    PolyQ P(g.q, {1});
    for (int64_t i = 0; i < e_ram; i++) {
      P = P * g;
    }
    PolyQ Qp = phid_q / P;
    // e = Qp * (Qp^{-1} mod P): 1 modulo P, 0 modulo Qp
    PolyQ e0 = Qp * PolyInverseModQ(Qp, P);
    PolyZ<Tint> ez = PolyZOfPolyQ<Tint>(e0);
    ez = PolyZReduce(ez, PhiD, Q);
    for (int it = 0; it < 64; it++) {
      PolyZ<Tint> e2 = PolyZReduce(PolyZMul(ez, ez), PhiD, Q);
      PolyZ<Tint> e3 = PolyZReduce(PolyZMul(e2, ez), PhiD, Q);
      PolyZ<Tint> en(e2.size());
      for (size_t i = 0; i < e2.size(); i++) {
        en[i] = 3 * e2[i] - 2 * e3[i];
      }
      en = PolyZReduce(en, PhiD, Q);
      if (en == ez) {
        break;
      }
      ez = en;
    }
    MyMatrix<T> E = ZeroMatrix<T>(r, r);
    MyMatrix<T> pw = Id;
    for (size_t i = 0; i < ez.size(); i++) {
      E += pw * T(ez[i]);
      pw *= S;
    }
    return E;
  };
  std::vector<HermitianLocalJordan> res;
  std::vector<bool> done(gs.size(), false);
  for (size_t i = 0; i < gs.size(); i++) {
    if (done[i]) {
      continue;
    }
    PolyQ gr = ReciprocalPolyQ(gs[i]);
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
    MyMatrix<T> E = idempotent_matrix(gs[i]);
    // the O_P-rank: the rank modulo q of E divided by [O_P : Z_q]
    int rk = static_cast<int>(RankModPrime(ModMatrixOfMatrix(E, q), q));
    int m = rk / static_cast<int>(e_ram * f);
    // len[k] = log_q |pi^k X|, X = e_P L^# / e_P L, until it vanishes
    std::vector<int> len;
    MyMatrix<T> Mk = Binv * E;
    while (true) {
      len.push_back(ValuationIndexOverlattice(Mk, q));
      if (len.back() == 0 || len.size() > static_cast<size_t>(4 * vdet + 8)) {
        break;
      }
      Mk *= Pi;
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

template <typename T>
inline T HermitianLocalFactorUnramified(HermitianLocalJordan const &hj,
                                                int64_t q) {
  if (hj.blocks.size() <= 1) {
    return T(1);
  }
  T qK = MyPow<T>(UniversalScalarConversion<T, int64_t>(q), hj.f_K);
  int m = 0;
  for (auto &b : hj.blocks) {
    m += b.second;
  }
  T f = HermitianGroupOrderNormalized<T>(hj.type, m, qK);
  int64_t N = 0;
  int before = 0;
  for (size_t i = 0; i < hj.blocks.size(); i++) {
    int s = hj.blocks[i].first;
    int mi = hj.blocks[i].second;
    int after = m - before - mi;
    N += static_cast<int64_t>(s) * mi * (before - after);
    f /= HermitianGroupOrderNormalized<T>(hj.type, mi, qK);
    before += mi;
  }
  T qpow(1);
  T qK_T = qK;
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

/*
  The Hermitian Gram matrix, in E, of an O_P-basis of L (x) Z_q at the
  ramified prime (d = q^a).
 */
template <typename T> struct HermitianLocalGram {
  NumberField<T> K;
  int m;
  std::vector<std::vector<ElemE<T>>> H;
};

template <typename T>
inline HermitianLocalGram<T>
ComputeHermitianLocalGram(MyMatrix<T> const &B,
                          MyMatrix<T> const &S, int64_t d, int64_t q) {
  int r = B.rows();
  NumberField<T> K;
  K.d = d;
  K.Phi = CyclotomicPolynomial(static_cast<int>(d));
  K.deg = static_cast<int>(K.Phi.size()) - 1;
  int m = r / K.deg;
  // the dual basis of (zeta^t) for the trace: T_{st} = Tr(zeta^{s+t})
  auto trace = [&](ElemE<T> const &x) -> T {
    // Tr(zeta^t) for t < deg, then linear
    T res(0);
    for (int t = 0; t < K.deg; t++) {
      if (x[t] == 0) {
        continue;
      }
      // Tr(zeta^t) is the Ramanujan sum c_d(t)
      T tr(0);
      int64_t g = std::gcd(t, d);
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
      tr = UniversalScalarConversion<T, int64_t>(mu * EulerPhi(d) / EulerPhi(dd));
      res += x[t] * tr;
    }
    return res;
  };
  MyMatrix<T> Tm(K.deg, K.deg);
  for (int s2 = 0; s2 < K.deg; s2++) {
    for (int t = 0; t < K.deg; t++) {
      Tm(s2, t) = trace(EZetaPower(s2 + t, K));
    }
  }
  MyMatrix<T> Tinv = Inverse(Tm);
  // an O_P-basis among the standard vectors: independent modulo pi L
  MyMatrix<T> Id = IdentityMat<T>(r);
  MyMatrix<T> Pi = Id - S;
  ModMatrix rows = ModMatrixOfMatrix(Pi, q);
  int base_rank = static_cast<int>(RankModPrime(rows, q));
  std::vector<int> chosen;
  for (int k = 0; k < r && static_cast<int>(chosen.size()) < m; k++) {
    // the rows of Pi and the unit vectors of chosen and of k
    int n_c = static_cast<int>(chosen.size());
    ModMatrix test = ZeroMatrix<T_residue>(r + n_c + 1, r);
    test.topRows(r) = rows;
    for (int i_c = 0; i_c < n_c; i_c++) {
      test(r + i_c, chosen[i_c]) = 1;
    }
    test(r + n_c, k) = 1;
    if (static_cast<int>(RankModPrime(test, q)) == base_rank + n_c + 1) {
      chosen.push_back(k);
    }
  }
  if (static_cast<int>(chosen.size()) != m) {
    std::cerr << "CYCLO: no O_P-basis found\n";
    throw TerminalException{1};
  }
  // h(x_i, x_j) = sum_t b(zeta^t x_i, x_j) zeta_t^v, zeta_t^v = sum_s
  // Tinv(t, s) zeta^s
  std::vector<MyMatrix<T>> Spow{Id};
  for (int t = 1; t < K.deg; t++) {
    Spow.push_back(Spow.back() * S);
  }
  std::vector<std::vector<ElemE<T>>> H(m, std::vector<ElemE<T>>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      ElemE<T> val(K.deg, T(0));
      for (int t = 0; t < K.deg; t++) {
        // b(zeta^t x_i, x_j) with x_i the chosen standard vector
        T bt(0);
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
template <typename T> struct DyadicContext {
  NumberField<T> K;
  // the valuation at the prime over 2, alone (1 - zeta)
  PrimeValuation<T> v{2, 1, nullptr};
  ElemE<T> one;
  ElemE<T> pi;
  ElemE<T> pi_inv;
  ElemE<T> pipibar;
  ElemE<T> pipibar_inv;
  std::vector<ElemE<T>> pipow;
};

template <typename T>
inline DyadicContext<T> MakeDyadicContext(NumberField<T> const &K) {
  DyadicContext<T> C;
  C.K = K;
  C.one = EZetaPower(0, K);
  C.pi = C.one - EZetaPower(1, K);
  C.pi_inv = EInv(C.pi, K);
  C.pipibar = EMul(C.pi, EConj(C.pi, K), K);
  C.pipibar_inv = EInv(C.pipibar, K);
  C.pipow = {C.one};
  for (int i = 1; i < 6; i++) {
    C.pipow.push_back(EMul(C.pipow.back(), C.pi, K));
  }
  return C;
}

template <typename T>
inline int DyadicVal(ElemE<T> const &x, DyadicContext<T> const &C) {
  return PrimeVal(x, C.K, C.v);
}

// in P-units: v_P of the ideal Tr(c O_E) O_E
template <typename T>
inline int DyadicTraceVal(ElemE<T> const &c, DyadicContext<T> const &C) {
  int v = DyadicVal(c, C);
  if (v >= (1 << 27)) {
    return v;
  }
  int num = v + 2;
  int vp = (num >= 0) ? num / 2 : -((-num + 1) / 2);
  return 2 * vp;
}

template <typename T>
inline int DyadicNormVal(std::vector<std::vector<ElemE<T>>> const &G,
                         DyadicContext<T> const &C) {
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
template <typename T>
inline bool DyadicIsLocalNorm(ElemE<T> x, DyadicContext<T> const &C) {
  NumberField<T> const &K = C.K;
  int w = DyadicVal(x, C);
  for (int i = 0; i < w / 2; i++) {
    x = EMul(x, C.pipibar_inv, K);
  }
  for (int i = 0; i < -w / 2; i++) {
    x = EMul(x, C.pipibar, K);
  }
  for (int64_t c = 0; c < 32; c++) {
    ElemE<T> y = C.one;
    for (int i = 1; i < 6; i++) {
      if (c & (int64_t(1) << (i - 1))) {
        for (int t = 0; t < K.deg; t++) {
          y[t] += C.pipow[i][t];
        }
      }
    }
    if (DyadicVal(x - EMul(y, EConj(y, K), K), C) >= 4) {
      return true;
    }
  }
  return false;
}

// U H U^*, the rows of U being the new basis in the old one
template <typename T>
inline std::vector<std::vector<ElemE<T>>>
HermitianCongruence(std::vector<std::vector<ElemE<T>>> const &H,
                    std::vector<std::vector<ElemE<T>>> const &U,
                    NumberField<T> const &K) {
  int m = static_cast<int>(H.size());
  std::vector<std::vector<ElemE<T>>> UH(m, std::vector<ElemE<T>>(m));
  for (int a = 0; a < m; a++) {
    for (int k = 0; k < m; k++) {
      ElemE<T> sum(K.deg, T(0));
      for (int i = 0; i < m; i++) {
        if (!U[a][i].is_zero()) {
          ElemE<T> t = EMul(U[a][i], H[i][k], K);
          for (int u = 0; u < K.deg; u++) {
            sum[u] += t[u];
          }
        }
      }
      UH[a][k] = sum;
    }
  }
  std::vector<std::vector<ElemE<T>>> R(m, std::vector<ElemE<T>>(m));
  for (int a = 0; a < m; a++) {
    for (int b = 0; b < m; b++) {
      ElemE<T> sum(K.deg, T(0));
      for (int k = 0; k < m; k++) {
        if (!U[b][k].is_zero()) {
          ElemE<T> t = EMul(UH[a][k], EConj(U[b][k], K), K);
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
template <typename T>
inline std::vector<std::vector<ElemE<T>>>
DyadicOverlattice(std::vector<std::vector<ElemE<T>>> const &H, int64_t mask,
                  DyadicContext<T> const &C) {
  int m = static_cast<int>(H.size());
  ElemE<T> zero(C.K.deg, T(0));
  std::vector<std::vector<ElemE<T>>> U(m, std::vector<ElemE<T>>(m, zero));
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
template <typename T>
inline std::vector<std::vector<ElemE<T>>>
DyadicSublattice(std::vector<std::vector<ElemE<T>>> const &H, int64_t mask,
                 DyadicContext<T> const &C) {
  int m = static_cast<int>(H.size());
  ElemE<T> zero(C.K.deg, T(0));
  std::vector<std::vector<ElemE<T>>> U(m, std::vector<ElemE<T>>(m, zero));
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
        U[i][j] = zero - C.one;
      }
    }
  }
  return HermitianCongruence(H, U, C.K);
}

template <typename T> struct KirschmerSymbolEntry {
  int rank;
  int scale;
  bool normal;
  int v;
  ElemE<T> det;
};

template <typename T>
inline std::vector<KirschmerSymbolEntry<T>>
KirschmerSymbol(std::vector<std::vector<ElemE<T>>> const &H, DyadicContext<T> const &C) {
  std::vector<LocalJordanComponent<T>> comps =
      LocalJordanSplit(H, C.K, C.v, JordanForm::Hermitian, nullptr);
  std::vector<int> scales;
  std::vector<int> nvals;
  std::vector<int> ranks;
  std::vector<ElemE<T>> dets;
  for (auto const &comp : comps) {
    int nv = 1 << 28;
    int rk = 0;
    ElemE<T> det = C.one;
    for (auto &G : comp.pieces) {
      nv = std::min(nv, DyadicNormVal(G, C));
      rk += static_cast<int>(G.size());
      det = EMul(det, JordanPieceDet(G, C.K), C.K);
    }
    scales.push_back(comp.s);
    nvals.push_back(nv);
    ranks.push_back(rk);
    dets.push_back(det);
  }
  std::vector<KirschmerSymbolEntry<T>> sym;
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
inline bool DyadicIsLocallyIsometric(std::vector<KirschmerSymbolEntry<T>> const &S1,
                                     std::vector<KirschmerSymbolEntry<T>> const &S2,
                                     DyadicContext<T> const &C) {
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
  if (!DyadicIsLocalNorm(ratio(t - 1), C)) {
    return false;
  }
  for (size_t i = 0; i + 1 < t; i++) {
    ElemE<T> x = ratio(i);
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

template <typename T>
inline std::optional<T>
HermitianLocalFactorDyadic(MyMatrix<T> const &B,
                           MyMatrix<T> const &S, int64_t d) {
  int64_t q = 2;
  HermitianLocalGram<T> hlg = ComputeHermitianLocalGram(B, S, d, q);
  DyadicContext<T> C = MakeDyadicContext(hlg.K);
  NumberField<T> const &K = C.K;
  int m = hlg.m;
  int64_t n_mask = int64_t(1) << m;
  // the factor of a maximal lattice of the same space
  auto maximal_factor = [&](std::vector<std::vector<ElemE<T>>> const &H) -> T {
    if (m % 2 == 1) {
      return RationalFraction<T>(1, 2);
    }
    std::vector<std::vector<ElemE<T>>> G = H;
    ElemE<T> det = C.one;
    for (int c = 0; c < m; c++) {
      int piv = -1;
      for (int i = c; i < m; i++) {
        if (!G[i][c].is_zero()) {
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
      ElemE<T> inv = EInv(G[c][c], K);
      for (int i = c + 1; i < m; i++) {
        ElemE<T> t = EMul(G[i][c], inv, K);
        for (int j = c; j < m; j++) {
          G[i][j] -= EMul(t, G[c][j], K);
        }
      }
    }
    if ((m * (m - 1) / 2) % 2 == 1) {
      for (auto &x : det) {
        x = -x;
      }
    }
    if (DyadicIsLocalNorm(det, C)) {
      return T(1);
    }
    T qm = MyPow<T>(2, m);
    return (qm - 1) / UniversalScalarConversion<T, int64_t>(2 * q + 2);
  };
  int nL = DyadicNormVal(hlg.H, C);
  // an overlattice of index 2 keeping the norm level, if any
  auto find_over = [&](std::vector<std::vector<ElemE<T>>> const &H)
      -> std::optional<std::vector<std::vector<ElemE<T>>>> {
    for (int64_t mask = 1; mask < n_mask; mask++) {
      std::vector<std::vector<ElemE<T>>> Hn = DyadicOverlattice(H, mask, C);
      if (DyadicNormVal(Hn, C) >= nL) {
        return Hn;
      }
    }
    return {};
  };
  std::vector<std::vector<std::vector<ElemE<T>>>> chain{hlg.H};
  while (true) {
    std::optional<std::vector<std::vector<ElemE<T>>>> opt = find_over(chain.back());
    if (!opt) {
      break;
    }
    chain.push_back(*opt);
    if (chain.size() > 64) {
      std::cerr << "CYCLO: the chain of overlattices does not end\n";
      throw TerminalException{1};
    }
  }
  T f = maximal_factor(hlg.H);
  for (size_t i = 0; i + 1 < chain.size(); i++) {
    std::vector<KirschmerSymbolEntry<T>> sym_i = KirschmerSymbol(chain[i], C);
    std::vector<KirschmerSymbolEntry<T>> sym_n = KirschmerSymbol(chain[i + 1], C);
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
      std::cerr << "CYCLO: a step of the chain has no isometric neighbor\n";
      throw TerminalException{1};
    }
    f *= RationalFraction<T>(A, Bc);
  }
  return f;
}

/*
  The local factor at the ramified prime.
 */
template <typename T>
inline std::optional<T>
HermitianLocalFactorRamified(MyMatrix<T> const &B,
                             MyMatrix<T> const &S, int64_t d,
                             int64_t q) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  if (q == 2) {
    return HermitianLocalFactorDyadic(B, S, d);
  }
  HermitianLocalGram<T> hlg = ComputeHermitianLocalGram(B, S, d, q);
  NumberField<T> const &K = hlg.K;
  int m = hlg.m;
  std::vector<std::vector<ElemE<T>>> H = hlg.H;
  // the Jordan splitting, the blocks of a scale gathered with their
  // determinant
  std::map<int, std::pair<int, ElemE<T>>> blocks;  // scale -> (rank, det)
  // the prime of E over q is alone (1 - zeta for q | d), its valuation
  // counted in the norm
  PrimeValuation<T> vP{q, 1, nullptr};
  for (auto const &comp :
       LocalJordanSplit(H, K, vP, JordanForm::Hermitian, nullptr)) {
    blocks[comp.s] = {comp.rank(), LocalComponentDet(comp, K)};
  }
  // the local factor of Hecke.jl at a ramified non-dyadic prime
  T qT = UniversalScalarConversion<T, int64_t>(q);
  T f = HeckeGroupOrderSp<T>(2 * (m / 2), qT);
  int64_t N = 0;
  int after_total = m;
  // (1 - zeta)(1 - zeta^{-1}) = 2 - zeta - zeta^{-1}
  ElemE<T> one = EZetaPower(0, K);
  ElemE<T> pipibar = one + one - EZetaPower(1, K) - EZetaPower(-1, K);
  for (auto &kv : blocks) {
    int s = kv.first;
    int mi = kv.second.first;
    after_total -= mi;
    int ri = after_total;
    int64_t fl = (s >= 0) ? s / 2 : -((-s + 1) / 2);
    N -= fl * mi * mi;
    if (s % 2 != 0) {
      N -= static_cast<int64_t>(mi + 1) * (mi / 2);
      f /= HeckeGroupOrderSp<T>(mi, qT);
    } else {
      // det = (pi pi-bar)^{s mi / 2} u
      ElemE<T> u = kv.second.second;
      int w = PrimeVal(u, K, vP);
      // divided by (pi pi-bar)^{w / 2}, w of either sign (a block of
      // negative scale when the form has denominators at P)
      ElemE<T> pinv = EInv(pipibar, K);
      for (int i = 0; i < w / 2; i++) {
        u = EMul(u, pinv, K);
      }
      for (int i = 0; i < -w / 2; i++) {
        u = EMul(u, pipibar, K);
      }
      // the discriminant (-1)^{m_i (m_i - 1) / 2} det, as in Hecke.jl
      Tint res_u = EResidue(u, q);
      if ((mi * (mi - 1) / 2) % 2 == 1) {
        res_u = -res_u;
      }
      int leg = LegendreModQ<Tint>(res_u, q);
      f /= HeckeGroupOrderO<T>(mi, leg, qT);
    }
    int64_t smi = static_cast<int64_t>(s) * mi;
    N -= smi * ri;
    N += smi * m / 2;
  }
  if (m % 2 == 0) {
    N += m / 2;
  }
  T qpow(1);
  for (int64_t i = 0; i < (N >= 0 ? N : -N); i++) {
    if (N >= 0) {
      qpow *= qT;
    } else {
      qpow /= qT;
    }
  }
  return qpow * f;
}

/*
  The mass of the genus of the Hermitian lattice over Q(zeta_d) given by
  the Gram matrix B of its trace form and the action S of zeta_d, or
  nothing at a ramified prime of E/K with rank >= 2.
 */
template <typename T>
inline std::optional<T>
HermitianMassCyclotomicLattice(MyMatrix<T> const &B,
                               MyMatrix<T> const &S, int64_t d_in) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int64_t d = NormalizeCyclotomicIndex(d_in);
  int r = B.rows();
  int64_t phi = EulerPhi(d);
  int m = static_cast<int>(r / phi);
  if (m == 1) {
    return RankOneHermitianMassCyclotomic<T>(d);
  }
  // For d_in = 2 d with d odd, Phi_{d_in}(S) = Phi_d(-S): -S is the action
  // of a primitive d-th root of unity.
  MyMatrix<T> Sd = S;
  if (d != d_in) {
    Sd = -S;
  }
  T res = StandardHermitianMassCyclotomic<T>(d, m);
  T det = DeterminantMat(B);
  std::set<int64_t> primes;
  {
    Tint num = GetNumerator_z(det);
    if (num < 0) {
      num = -num;
    }
    if (num > 1) {
      for (Tint const &p : FactorsInt(num)) {
        primes.insert(UniversalScalarConversion<int64_t, Tint>(p));
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
        std::optional<T> opt = HermitianLocalFactorRamified(B, Sd, d, q);
        if (!opt) {
          return {};
        }
        res *= *opt;
      } else {
        res *= HermitianLocalFactorUnramified<T>(hj, q);
      }
    }
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_MASSCYCLOTOMICHERMITIAN_H_
// clang-format on
