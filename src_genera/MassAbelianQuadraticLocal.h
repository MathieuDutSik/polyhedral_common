// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSABELIANQUADRATICLOCAL_H_
#define SRC_GENERA_MASSABELIANQUADRATICLOCAL_H_

// clang-format off
#include "MAT_MatrixModPrime.h"
#include "ScalarNumberField.h"
#include "ScalarAbelianField.h"
#include "LocalGramOverOK.h"
#include "LocalJordan.h"
#include "ScalarFiniteField.h"
#include "ClassicalGroupOrders.h"
#include "ArithmeticFunctions.h"
#include "ResidueArithmetic.h"
#include "MassAbelianHermitian.h"
#include "ClassicLLL.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_ABELIAN_QUADRATIC_MASS
#endif

/*
  The quadratic masses over a totally real abelian field K, their local
  part (MassAbelianQuadratic.h for the global mass of a lattice): the
  standard mass from the Dedekind zeta values of K, the local factors at
  the odd and the dyadic primes, the square classes, the Hilbert symbols
  and the Witt invariants at the dyadic primes, the fields K and E =
  K(sqrt d), and the choice of the L-value method.

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

template <typename T>
inline T DedekindZetaAbelian(AbelianField const &F, int k) {
  // zeta_K(1 - k): all the characters of K, even as K is real
  return ProductDirichletLValues<T>(F.n, k, 0, F.H);
}

/*
  The standard mass; for m even with a non-square discriminant, FE gives
  E = K(sqrt disc) as the subfield of Q(zeta_n) fixed by H_E (same n).
  Nothing when E is needed and not given.
 */
template <typename T>
inline std::optional<T>
StandardQuadraticMassAbelian(AbelianField const &FK, int m, bool disc_square,
                             std::optional<AbelianField> const &FE) {
  if (m == 0) {
    return T(1);
  }
  if (m == 1) {
    return T(1) / T(2);
  }
  int64_t deg = AbelianFieldDegree(FK);
  int r = m / 2;
  T res(1);
  for (int64_t i = 0; i < deg * r; i++) {
    res /= 2;
  }
  int top = (m % 2 == 1) ? m - 2 : m - 3;
  for (int i = 1; i <= top; i += 2) {
    res *= DedekindZetaAbelian<T>(FK, i + 1);
  }
  if (m % 2 == 0) {
    if (disc_square) {
      res *= DedekindZetaAbelian<T>(FK, r);
    } else {
      if (!FE || FE->n != FK.n) {
        return {};
      }
      // the characters of E not of K, whatever their parity
      T prod(1);
      for (int parity = 0; parity <= 1; parity++) {
        prod *= ProductDirichletLValues<T>(FK.n, r, parity, FE->H, FK.H);
      }
      res *= prod;
    }
  }
  return (res < 0) ? T(-res) : res;
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
template <typename T>
inline T QuadLocalFactorOdd(std::vector<QuadJordanBlock> blocks,
                            T const &qf) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  std::sort(blocks.begin(), blocks.end(),
            [](QuadJordanBlock const &a, QuadJordanBlock const &b) {
              return a.s < b.s;
            });
  if (blocks.size() <= 1) {
    return T(1);
  }
  int chi_m1 = (ResInt(GetNumerator_z(qf), Tint(4)) == 1) ? 1 : -1;
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
  T f;
  if (m % 2 == 1) {
    f = HeckeGroupOrderO<T>(m, 1, qf);
  } else if (vd % 2 != 0) {
    f = HeckeGroupOrderO<T>(m - 1, 1, qf);
  } else {
    f = HeckeGroupOrderO<T>(m, chid, qf);
  }
  // twice the exponent N
  int64_t N2 = 0;
  int after = m;
  for (auto const &bl : blocks) {
    int64_t s = bl.s, mi = bl.m;
    after -= bl.m;
    int64_t ri = after;
    int sq = ((s * mi) % 2 != 0) ? -1 : sign_class(bl.m * (bl.m - 1) / 2) * bl.chi;
    f /= HeckeGroupOrderO<T>(bl.m, sq, qf);
    N2 -= s * mi * (mi + 1) + 2 * s * mi * ri;
    N2 += s * mi * (m + 1);
  }
  if (m % 2 == 0 && vd % 2 != 0) {
    N2 += 1 - m;
  }
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
  if (N2 % 2 != 0) {
    std::cerr << "ABQUAD: the exponent N is not integral\n";
    throw TerminalException{1};
  }
#endif
  int64_t N = N2 / 2;
  T qpow = MyPow<T>(qf, (N >= 0 ? N : -N));
  return (N >= 0) ? T(qpow * f) : T(f / qpow);
}

// The square class (1 or -1) in F_q[x] / g of a residue, q odd.
template <typename T>
inline int ResidueSquareClass(FqElt const &u, LocalPrimeK<T> const &P) {
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
  if (u.is_zero()) {
    std::cerr << "ABQUAD: the square class of zero\n";
    throw TerminalException{1};
  }
#endif
  return FqQuadraticCharacter(P.F, u);
}

// The local factor at an odd prime p from the Gram matrix over K.
template <typename T>
inline T LocalFactorOddGram(std::vector<std::vector<ElemE<T>>> const &H,
                                    LocalPrimeK<T> &P) {
  std::vector<QuadJordanBlock> blocks;
  for (auto const &C : LocalJordanSplit(H, P.K, PrimeValuation<T>{P.q, 1, &P},
                                     JordanForm::Quadratic, &P.pi)) {
    int mi = C.rank();
    int chi = 1;
    if ((C.s * mi) % 2 == 0) {
      ElemE<T> u = LocalShift(LocalComponentDet(C, P.K), C.s * mi, P);
      chi = ResidueSquareClass(LocalResidue(u, P), P);
    }
    blocks.push_back({C.s, mi, chi});
  }
  return QuadLocalFactorOdd<T>(
      blocks, UniversalScalarConversion<T, typename underlying_z_ring<T>::ring_type>(P.qf));
}

// The matrices over a residue field F (ScalarFiniteField.h).
using FqMat = std::vector<std::vector<FqElt>>;

// A basis of {v : v M = 0} for M with nr rows.
inline FqMat FqLeftKernel(FiniteField const &F, FqMat const &M, size_t nr, size_t nc) {
  // row reduce [M | I]
  FqMat A(nr, std::vector<FqElt>(nc + nr));
  for (size_t i = 0; i < nr; i++) {
    for (size_t j = 0; j < nc; j++) {
      A[i][j] = M[i][j];
    }
    A[i][nc + i] = FqOne(F);
  }
  size_t row = 0;
  for (size_t col = 0; col < nc && row < nr; col++) {
    size_t piv = row;
    while (piv < nr && A[piv][col].is_zero()) {
      piv++;
    }
    if (piv == nr) {
      continue;
    }
    std::swap(A[piv], A[row]);
    FqElt inv = FqInv(F, A[row][col]);
    for (auto &x : A[row]) {
      x = FqMul(F, x, inv);
    }
    for (size_t i = 0; i < nr; i++) {
      if (i != row && !A[i][col].is_zero()) {
        FqElt c = A[i][col];
        for (size_t j = 0; j < nc + nr; j++) {
          A[i][j] = FqAdd(F, A[i][j], FqNeg(F, FqMul(F, c, A[row][j])));
        }
      }
    }
    row++;
  }
  FqMat res;
  for (size_t i = row; i < nr; i++) {
    res.push_back(std::vector<FqElt>(A[i].begin() + nc, A[i].end()));
  }
  return res;
}

inline FqElt FqBil(FiniteField const &F, std::vector<FqElt> const &x, FqMat const &Q,
                   std::vector<FqElt> const &y) {
  FqElt res;
  for (size_t i = 0; i < x.size(); i++) {
    if (x[i].is_zero()) {
      continue;
    }
    for (size_t j = 0; j < y.size(); j++) {
      if (!y[j].is_zero() && !Q[i][j].is_zero()) {
        res = FqAdd(F, res, FqMul(F, FqMul(F, x[i], Q[i][j]), y[j]));
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
inline std::pair<int, int> QuadFormTypeChar2(FqMat const &Q, FiniteField const &F) {
  size_t n = Q.size();
  if (n == 0) {
    return {0, 1};
  }
  FqMat Bf(n, std::vector<FqElt>(n));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      Bf[i][j] = FqAdd(F, Q[i][j], Q[j][i]);
    }
  }
  FqMat N = FqLeftKernel(F, Bf, n, n);
  bool qzero = true;
  for (auto const &v : N) {
    if (!FqBil(F, v, Q, v).is_zero()) {
      qzero = false;
    }
  }
  int dim_rad = static_cast<int>(N.size()) - (qzero ? 0 : 1);
  int mi = static_cast<int>(n) - dim_rad;
  if (mi % 2 == 1 || mi == 0) {
    return {mi, 1};
  }
  // a complement of N: standard vectors independent modulo N
  FqMat basis = N;
  FqMat comp;
  for (size_t k = 0; k < n && basis.size() < n; k++) {
    std::vector<FqElt> e(n);
    e[k] = FqOne(F);
    FqMat test = basis;
    test.push_back(e);
    // independent iff the left kernel of test (as rows) is trivial
    if (FqLeftKernel(F, test, test.size(), n).empty()) {
      basis = test;
      comp.push_back(e);
    }
  }
  // symplectic reduction on comp
  FqElt arf;
  while (!comp.empty()) {
    std::vector<FqElt> x = comp[0];
    size_t j = 1;
    while (j < comp.size() && FqBil(F, x, Bf, comp[j]).is_zero()) {
      j++;
    }
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
    if (j == comp.size()) {
      std::cerr << "ABQUAD: degenerate polar form\n";
      throw TerminalException{1};
    }
#endif
    FqElt c = FqInv(F, FqBil(F, x, Bf, comp[j]));
    std::vector<FqElt> y = comp[j];
    for (auto &t : y) {
      t = FqMul(F, t, c);
    }
    // b(x, y) = 1
    arf = FqAdd(F, arf, FqMul(F, FqBil(F, x, Q, x), FqBil(F, y, Q, y)));
    FqMat rest;
    for (size_t k = 1; k < comp.size(); k++) {
      if (k == j) {
        continue;
      }
      std::vector<FqElt> z = comp[k];
      // z <- z + b(z, y) x + b(z, x) y (char 2)
      FqElt a1 = FqBil(F, z, Bf, y), a2 = FqBil(F, z, Bf, x);
      for (size_t u = 0; u < n; u++) {
        z[u] = FqAdd(F, z[u], FqAdd(F, FqMul(F, a1, x[u]), FqMul(F, a2, y[u])));
      }
      rest.push_back(z);
    }
    comp = rest;
  }
  return {mi, (FqTrace(F, arf) == 0) ? 1 : -1};
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
template <typename T>
inline int DyadicUnramifiedSquareClass(ElemE<T> const &d, LocalPrimeK<T> &P,
                                       FiniteField const &F) {
  NumberField<T> const &K = P.K;
  int v = LocalVal(d, P);
  if (v % 2 != 0) {
    return 3;
  }
  ElemE<T> u = LocalShift(d, v, P);
  FqElt ur = LocalResidue(u, P);
  FqElt wr = FqFrobenius(F, ur, F.f - 1);
  ElemE<T> w(K.deg, T(0));
  for (size_t i = 0; i < wr.size(); i++) {
    w[i] = UniversalScalarConversion<T, int64_t>(wr[i]);
  }
  ElemE<T> x = EMul(u, EInv(EMul(w, w, K), K), K);
  ElemE<T> x1 = x;
  x1[0] -= 1;
  int vx = LocalVal(x1, P);
  if (vx == 1) {
    return 2;
  }
  if (vx >= 3) {
    return 0;
  }
  FqElt t = LocalResidue(LocalShift(x1, 2, P), P);
  return (FqTrace(F, t) == 0) ? 0 : -1;
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
template <typename T>
inline T LocalFactorCho(std::vector<std::vector<ElemE<T>>> const &H,
                                LocalPrimeK<T> &P) {
  NumberField<T> const &K = P.K;
  FiniteField F = MakeFiniteField(static_cast<T_small_prime>(P.q), P.g);
  T q = UniversalScalarConversion<T, typename underlying_z_ring<T>::ring_type>(P.qf);
  std::vector<LocalJordanComponent<T>> comps = LocalJordanSplit(H, P.K, PrimeValuation<T>{P.q, 1, &P},
                                     JordanForm::Quadratic, &P.pi);
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
  std::vector<std::vector<ElemE<T>>> AG0(m, std::vector<ElemE<T>>(m, ElemE<T>(K.deg, T(0))));
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
  auto pow2 = [&](int k) -> T {
    T r = MyPow<T>(2, (k >= 0 ? k : -k));
    return (k >= 0) ? r : T(1 / r);
  };
  // the forms V_s: dimension and type
  std::vector<std::pair<int, int>> V(nc);
  for (size_t is = 0; is < nc; is++) {
    int s = S[is];
    std::vector<std::vector<ElemE<T>>> AG = AG0;
    for (int a = 0; a < m; a++) {
      for (int b = 0; b < m; b++) {
        int Sj = S[comp_of[a]];
        if (Sj < s) {
          AG[a][b] *= pow2(2 * (s - Sj));
        }
      }
    }
    FqMat Dcol(m, std::vector<FqElt>(1));
    for (int a = 0; a < m; a++) {
      Dcol[a][0] = LocalResidue(AG[a][a] * pow2(-s), P);
    }
    FqMat Bk = FqLeftKernel(F, Dcol, m, 1);
    size_t nb = Bk.size();
    std::vector<std::vector<ElemE<T>>> BK(nb, std::vector<ElemE<T>>(m, ElemE<T>(K.deg, T(0))));
    for (size_t i = 0; i < nb; i++) {
      for (int a = 0; a < m; a++) {
        FqElt r = FqFrobenius(F, Bk[i][a], F.f - 1);
        for (size_t t = 0; t < r.size(); t++) {
          BK[i][a][t] = UniversalScalarConversion<T, int64_t>(r[t]);
        }
      }
    }
    FqMat Q(nb, std::vector<FqElt>(nb));
    for (size_t i = 0; i < nb; i++) {
      for (size_t j = i; j < nb; j++) {
        ElemE<T> val(K.deg, T(0));
        for (int a = 0; a < m; a++) {
          if (BK[i][a].is_zero()) {
            continue;
          }
          for (int b = 0; b < m; b++) {
            if (BK[j][b].is_zero() || AG[a][b].is_zero()) {
              continue;
            }
            val += EMul(EMul(BK[i][a], AG[a][b], K), BK[j][b], K);
          }
        }
        val *= pow2(-(s + 1));
        if (i < j) {
          val *= T(2);
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
  T res = pow2(n_alpha + n_beta);
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
  // after: the total rank of the components after the i-th
  int64_t after = 0;
  for (size_t i = 0; i < nc; i++) {
    after += M[i];
  }
  for (size_t i = 0; i < nc; i++) {
    after -= M[i];
    int64_t sm = static_cast<int64_t>(S[i]) * M[i];
    N += sm * (M[i] + 1) / 2 + sm * after;
  }
  for (auto const &v : V) {
    if (v.first == 0) {
      continue;
    }
    res *= HeckeGroupOrderO<T>(v.first, v.second, q) / 2;
  }
  auto qpow = [&](int64_t k) -> T {
    T r = MyPow<T>(q, (k >= 0 ? k : -k));
    return (k >= 0) ? r : T(1 / r);
  };
  T beta = qpow(N) * res / 2;
  // twice the exponent
  int64_t exp2 = 0;
  for (size_t i = 0; i < nc; i++) {
    exp2 += static_cast<int64_t>(m + 1) * S[i] * M[i];
  }
  T Hord;
  if (m % 2 == 1) {
    exp2 += m + 1;
    Hord = HeckeGroupOrderO<T>(m, 1, q);
  } else {
    exp2 += 2 * m;
    ElemE<T> d = EZetaPower(0, K);
    for (auto const &C : comps) {
      d = EMul(d, LocalComponentDet(C, K), K);
    }
    if ((m * (m - 1) / 2) % 2 == 1) {
      d = -d;
    }
    int sc = DyadicUnramifiedSquareClass(d, P, F);
    if (sc == 0) {
      Hord = HeckeGroupOrderO<T>(m, 1, q);
    } else if (sc == -1) {
      Hord = HeckeGroupOrderO<T>(m, -1, q);
    } else {
      Hord = HeckeGroupOrderO<T>(m - 1, 1, q);
      exp2 += static_cast<int64_t>(sc) * (1 - m);
    }
  }
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
  if (exp2 % 2 != 0) {
    std::cerr << "ABQUAD: the exponent of Cho's factor is not integral\n";
    throw TerminalException{1};
  }
#endif
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
template <typename T> using KMat = std::vector<std::vector<ElemE<T>>>;
constexpr int QD_INF = 1 << 28;

template <typename T>
inline ElemE<T> LiftResidue(FqElt const &r, NumberField<T> const &K) {
  ElemE<T> x(K.deg, T(0));
  for (size_t i = 0; i < r.size(); i++) {
    x[i] = UniversalScalarConversion<T, int64_t>(r[i]);
  }
  return x;
}

// (the defect of the unit a, a in the same square class with a - 1 of that
// valuation), at a dyadic prime
template <typename T>
inline std::pair<int, ElemE<T>> QuadDefectUnit(ElemE<T> a, LocalPrimeK<T> &P,
                                            FiniteField const &F) {
  NumberField<T> const &K = P.K;
  auto sqrt_res = [&](ElemE<T> const &x) -> ElemE<T> {
    return LiftResidue(FqFrobenius(F, LocalResidue(x, P), F.f - 1), K);
  };
  ElemE<T> s = sqrt_res(a);
  a = EMul(a, EInv(EMul(s, s, K), K), K);
  auto val1 = [&](ElemE<T> const &x) -> int {
    ElemE<T> y = x;
    y[0] -= 1;
    return y.is_zero() ? QD_INF : LocalVal(y, P);
  };
  int w = val1(a);
  int ee = 2 * P.e;
  while (w < ee && w % 2 == 0) {
    ElemE<T> am1 = a;
    am1[0] -= 1;
    ElemE<T> t = sqrt_res(LocalShift(am1, w, P));
    ElemE<T> c = EOne(K);
    ElemE<T> pw = EOne(K);
    for (int i = 0; i < w / 2; i++) {
      pw = EMul(pw, P.pi, K);
    }
    c += EMul(t, pw, K);
    a = EMul(a, EInv(EMul(c, c, K), K), K);
    w = val1(a);
  }
  if (w > ee) {
    return {QD_INF, EOne(K)};
  }
  if (w == ee) {
    ElemE<T> am1 = a;
    am1[0] -= 1;
    FqElt c = LocalResidue(am1 * (T(1) / T(4)), P);
    if (FqTrace(F, c) == 0) {
      return {QD_INF, EOne(K)};
    }
  }
  return {w, a};
}

template <typename T>
inline int QuadDefect(ElemE<T> const &a, LocalPrimeK<T> &P, FiniteField const &F) {
  if (a.is_zero()) {
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
template <typename T>
inline int HilbertSymbolDyadic(ElemE<T> a, ElemE<T> b, LocalPrimeK<T> &P,
                               FiniteField const &F) {
  NumberField<T> const &K = P.K;
  int v = LocalVal(a, P);
  int w = LocalVal(b, P);
  if (v % 2 != 0) {
    if (w % 2 == 0) {
      std::swap(a, b);
      std::swap(v, w);
    } else {
      a = -EMul(a, b, K);
      v += w;
    }
  }
  a = LocalShift(a, v, P);
  auto [vd, a2] = QuadDefectUnit(a, P, F);
  a = a2;
  int wf = (w >= 0) ? w / 2 : -((-w + 1) / 2);
  b = LocalShift(b, 2 * wf, P);
  w -= 2 * wf;
  if (vd != QD_INF && vd % 2 != 0 && w == 0) {
    ElemE<T> oma = -a;
    oma[0] += 1;
    b = LocalShift(EMul(oma, b, K), vd - 1, P);
    w = 1;
  }
  while (vd != QD_INF && vd % 2 != 0) {
    ElemE<T> t = b;
    for (int i = 0; i < vd - 1; i++) {
      t = EMul(t, P.pi, K);
    }
    ElemE<T> am1 = a;
    am1[0] -= 1;
    ElemE<T> sq = LiftResidue(
        FqFrobenius(F, LocalResidue(EMul(am1, EInv(t, K), K), P), F.f - 1), K);
    ElemE<T> fac = -EMul(EMul(sq, sq, K), t, K);
    fac[0] += 1;
    a = EMul(a, fac, K);
    auto pr = QuadDefectUnit(a, P, F);
    vd = pr.first;
    a = pr.second;
  }
  return (vd != QD_INF && w % 2 != 0) ? -1 : 1;
}

// A diagonal form isometric to the nondegenerate G.
template <typename T>
inline std::vector<ElemE<T>> DiagonalizeK(KMat<T> F, NumberField<T> const &K) {
  int n = static_cast<int>(F.size());
  std::vector<ElemE<T>> res;
  for (int i = 0; i < n; i++) {
    if (F[i][i].is_zero()) {
      int j = -1;
      for (int t = i + 1; t < n; t++) {
        if (!F[t][t].is_zero()) {
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
          if (!F[i][t].is_zero()) {
            j = t;
            break;
          }
        }
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
        if (j < 0) {
          std::cerr << "ABQUAD: a degenerate form\n";
          throw TerminalException{1};
        }
#endif
        // x_i <- x_i + c x_j, c = 1 / (2 F_ij)
        ElemE<T> c = EInv(F[i][j] * T(2), K);
        for (int u = 0; u < n; u++) {
          F[i][u] += EMul(c, F[j][u], K);
        }
        for (int u = 0; u < n; u++) {
          F[u][i] += EMul(F[u][j], c, K);
        }
      }
    }
    ElemE<T> dinv = EInv(F[i][i], K);
    for (int t = i + 1; t < n; t++) {
      ElemE<T> c = -EMul(F[t][i], dinv, K);
      for (int u = 0; u < n; u++) {
        F[t][u] += EMul(c, F[i][u], K);
      }
      for (int u = 0; u < n; u++) {
        F[u][t] += EMul(F[u][i], c, K);
      }
    }
    res.push_back(F[i][i]);
  }
  return res;
}

template <typename T>
inline int HasseInvariantDyadic(std::vector<ElemE<T>> const &D, LocalPrimeK<T> &P,
                                FiniteField const &F) {
  int h = 1;
  for (size_t i = 0; i < D.size(); i++) {
    for (size_t j = i + 1; j < D.size(); j++) {
      h *= HilbertSymbolDyadic(D[i], D[j], P, F);
    }
  }
  return h;
}

template <typename T>
inline ElemE<T> DetK(KMat<T> const &G, NumberField<T> const &K) {
  ElemE<T> d = EOne(K);
  for (auto const &x : DiagonalizeK(G, K)) {
    d = EMul(d, x, K);
  }
  return d;
}

// The Witt invariant of the space of G (Hecke: _hasse_witt).
template <typename T>
inline int WittInvariantDyadic(KMat<T> const &G, LocalPrimeK<T> &P,
                               FiniteField const &F) {
  NumberField<T> const &K = P.K;
  std::vector<ElemE<T>> D = DiagonalizeK(G, K);
  int h = HasseInvariantDyadic(D, P, F);
  ElemE<T> d = EOne(K);
  for (auto const &x : D) {
    d = EMul(d, x, K);
  }
  int n = static_cast<int>(D.size()) % 8;
  ElemE<T> c = EOne(K);
  if (n == 3 || n == 4) {
    c = -d;
  } else if (n == 5 || n == 6) {
    c = -EOne(K);
  } else if (n == 7 || n == 0) {
    c = d;
  }
  return h * HilbertSymbolDyadic(-EOne(K), c, P, F);
}

template <typename T>
inline KMat<T> ScaleKMat(KMat<T> G, ElemE<T> const &c, NumberField<T> const &K) {
  for (auto &row : G) {
    for (auto &x : row) {
      x = EMul(x, c, K);
    }
  }
  return G;
}

template <typename T>
inline ElemE<T> PiPower(int k, LocalPrimeK<T> const &P) {
  return LocalShift(EOne(P.K), -k, P);
}

template <typename T>
inline int NormValK(KMat<T> const &G, LocalPrimeK<T> &P) {
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

template <typename T>
inline int ScaleValK(KMat<T> const &G, LocalPrimeK<T> &P) {
  int v = QD_INF;
  for (auto const &row : G) {
    for (auto const &x : row) {
      v = std::min(v, LocalVal(x, P));
    }
  }
  return v;
}

template <typename T>
inline T QPowRat(T const &q, int64_t k) {
  T r = MyPow<T>(q, (k >= 0 ? k : -k));
  return (k >= 0) ? r : T(1 / r);
}

// Hecke: _local_factor_maximal, on the space of G
template <typename T>
inline T LocalFactorMaximalDyadic(KMat<T> const &G, LocalPrimeK<T> &P,
                                          FiniteField const &F) {
  NumberField<T> const &K = P.K;
  int m = static_cast<int>(G.size());
  int r = m / 2;
  if (m == 1) {
    return T(1);
  }
  int w = WittInvariantDyadic(G, P, F);
  ElemE<T> d = DetK(G, K);
  if ((m * (m - 1) / 2) % 2 == 1) {
    d = -d;
  }
  T q = UniversalScalarConversion<T, typename underlying_z_ring<T>::ring_type>(P.qf);
  int vd = LocalVal(d, P);
  if (m % 2 == 1) {
    if (vd % 2 != 0) {
      return (QPowRat(q, r) + w) / 2;
    } else if (w == -1) {
      return (QPowRat(q, m - 1) - 1) / T(q + 1) / 2;
    }
  } else {
    int qd = QuadDefect(d, P, F);
    if (vd % 2 != 0) {
      return T(1) / T(2);
    } else if (qd == QD_INF) {
      if (w == -1) {
        return (QPowRat(q, r - 1) - 1) * (QPowRat(q, r) - 1) /
               T(q + 1) / 2;
      }
    } else if (qd % 2 != 0) {
      return T(1) / T(2);
    } else if (w == -1) {
      return (QPowRat(q, r - 1) + 1) * (QPowRat(q, r) + 1) /
             T(q + 1) / 2;
    }
  }
  return T(1);
}

/*
  The genus symbol of Kirschmer at a dyadic prime: for the Jordan
  components G_i of scales s_i, the norm generators a_i, u_i = v(a_i), the
  weights w_i and the f_i.
 */
template <typename T> struct QuadKirschmerSymbol {
  std::vector<int> dims, sL, wL, uL, fL;
  std::vector<ElemE<T>> aL;
  std::vector<std::vector<ElemE<T>>> diag;
  std::vector<ElemE<T>> dets;
};

template <typename T>
inline QuadKirschmerSymbol<T> KirschmerSymbolDyadic(KMat<T> const &G, LocalPrimeK<T> &P,
                                             FiniteField const &F) {
  NumberField<T> const &K = P.K;
  std::vector<LocalJordanComponent<T>> comps = LocalJordanSplit(G, P.K, PrimeValuation<T>{P.q, 1, &P},
                                     JordanForm::Quadratic, &P.pi);
  QuadKirschmerSymbol<T> S;
  int t = static_cast<int>(comps.size());
  int e = P.e;
  // the diagonal entries of each component
  std::vector<std::vector<ElemE<T>>> Dc(t);
  for (int i = 0; i < t; i++) {
    S.dims.push_back(comps[i].rank());
    S.sL.push_back(comps[i].s);
    for (auto const &Gp : comps[i].pieces) {
      for (size_t a = 0; a < Gp.size(); a++) {
        Dc[i].push_back(Gp[a][a]);
      }
    }
    // a diagonalization of each piece, for the Hasse invariants
    std::vector<ElemE<T>> dg;
    for (auto const &Gp : comps[i].pieces) {
      for (auto const &x : DiagonalizeK(Gp, K)) {
        dg.push_back(x);
      }
    }
    S.diag.push_back(dg);
    S.dets.push_back(LocalComponentDet(comps[i], K));
  }
  for (int i = 0; i < t; i++) {
    std::vector<ElemE<T>> D;
    for (int j = 0; j < t; j++) {
      for (auto const &x : Dc[j]) {
        D.push_back(j < i ? EMul(x, PiPower(2 * (S.sL[i] - S.sL[j]), P), K)
                          : x);
      }
    }
    int mv = QD_INF, pos = -1;
    for (size_t k = 0; k < D.size(); k++) {
      int v = D[k].is_zero() ? QD_INF : LocalVal(D[k], P);
      if (v < mv) {
        mv = v;
        pos = static_cast<int>(k);
      }
    }
    ElemE<T> a = (e + S.sL[i] <= mv) ? PiPower(e + S.sL[i], P) : D[pos];
    S.aL.push_back(a);
    int u = LocalVal(a, P);
    S.uL.push_back(u);
    int w = e + S.sL[i];
    ElemE<T> ainv = EInv(a, K);
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
template <typename T>
inline bool IsLocallyIsometricDyadic(QuadKirschmerSymbol<T> const &S1,
                                     QuadKirschmerSymbol<T> const &S2, LocalPrimeK<T> &P,
                                     FiniteField const &F) {
  NumberField<T> const &K = P.K;
  if (S1.dims != S2.dims || S1.sL != S2.sL || S1.wL != S2.wL ||
      S1.uL != S2.uL) {
    return false;
  }
  int e = P.e;
  for (size_t k = 0; k < S1.aL.size(); k++) {
    ElemE<T> b = EMul(S1.aL[k], EInv(S2.aL[k], K), K);
    int qd = QuadDefect(b, P, F);
    if (qd < S1.wL[k] - S2.uL[k]) {
      return false;
    }
  }
  auto comp_hasse = [&](size_t i, ElemE<T> const &a) -> bool {
    std::vector<ElemE<T>> D1, D2;
    for (size_t j = 0; j <= i; j++) {
      for (auto const &x : S1.diag[j]) {
        D1.push_back(x);
      }
      for (auto const &x : S2.diag[j]) {
        D2.push_back(x);
      }
    }
    D2.push_back(a);
    ElemE<T> prod = EOne(K);
    for (auto const &x : D1) {
      prod = EMul(prod, x, K);
    }
    for (auto const &x : D2) {
      prod = EMul(prod, x, K);
    }
    D1.push_back(prod);
    return HasseInvariantDyadic(D1, P, F) == HasseInvariantDyadic(D2, P, F);
  };
  ElemE<T> d1 = EOne(K), d2 = EOne(K);
  for (size_t i = 0; i < S1.fL.size(); i++) {
    d1 = EMul(d1, S1.dets[i], K);
    d2 = EMul(d2, S2.dets[i], K);
    ElemE<T> dq = EMul(d1, EInv(d2, K), K);
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
inline std::vector<std::vector<FqElt>> FqAllProjectivePoints(FiniteField const &F, int m) {
  std::vector<FqElt> elts;
  FqCheckCodes(F);
  int64_t qf = static_cast<int64_t>(F.q);
  for (int64_t u = 0; u < qf; u++) {
    elts.push_back(FqDecode(F, static_cast<uint64_t>(u)));
  }
  std::vector<std::vector<FqElt>> res;
  for (int lead = 0; lead < m; lead++) {
    int nfree = m - lead - 1;
    int64_t tot = MyPow<int64_t>(qf, nfree);
    for (int64_t u = 0; u < tot; u++) {
      std::vector<FqElt> v(m);
      v[lead] = FqOne(F);
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

template <typename T>
inline KMat<T> TransformKMat(KMat<T> const &Tm, KMat<T> const &G, NumberField<T> const &K) {
  size_t n = G.size();
  KMat<T> TG(n, std::vector<ElemE<T>>(n, ElemE<T>(K.deg, T(0))));
  for (size_t i = 0; i < n; i++) {
    for (size_t k = 0; k < n; k++) {
      if (Tm[i][k].is_zero()) {
        continue;
      }
      for (size_t j = 0; j < n; j++) {
        TG[i][j] += EMul(Tm[i][k], G[k][j], K);
      }
    }
  }
  KMat<T> res(n, std::vector<ElemE<T>>(n, ElemE<T>(K.deg, T(0))));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < n; k++) {
        if (!Tm[j][k].is_zero()) {
          res[i][j] += EMul(TG[i][k], Tm[j][k], K);
        }
      }
    }
  }
  return res;
}

// L + pi^{-1} y with y = sum c_i b_i
template <typename T>
inline KMat<T> SuperlatticeGram(KMat<T> const &G, std::vector<FqElt> const &c,
                             LocalPrimeK<T> const &P) {
  NumberField<T> const &K = P.K;
  size_t n = G.size();
  KMat<T> Tm(n, std::vector<ElemE<T>>(n, ElemE<T>(K.deg, T(0))));
  size_t lead = 0;
  while (c[lead].is_zero()) {
    lead++;
  }
  ElemE<T> pinv = EInv(P.pi, K);
  for (size_t i = 0; i < n; i++) {
    if (i == lead) {
      for (size_t j = 0; j < n; j++) {
        Tm[i][j] = EMul(LiftResidue(c[j], K), pinv, K);
      }
    } else {
      Tm[i][i] = EOne(K);
    }
  }
  return TransformKMat(Tm, G, K);
}

// {x in L : sum lambda_i x_i = 0 mod p}
template <typename T>
inline KMat<T> SublatticeGram(KMat<T> const &G, std::vector<FqElt> const &lambda,
                           LocalPrimeK<T> const &P) {
  NumberField<T> const &K = P.K;
  size_t n = G.size();
  KMat<T> Tm(n, std::vector<ElemE<T>>(n, ElemE<T>(K.deg, T(0))));
  size_t lead = 0;
  while (lambda[lead].is_zero()) {
    lead++;
  }
  for (size_t i = 0; i < n; i++) {
    if (i == lead) {
      Tm[i][i] = P.pi;
    } else {
      Tm[i][i] = EOne(K);
      Tm[i][lead] = -LiftResidue(lambda[i], K);
    }
  }
  return TransformKMat(Tm, G, K);
}

// An integral overlattice of index p of G (of integral norm), if any.
template <typename T>
inline std::optional<KMat<T>> IntegralOverlattice(KMat<T> const &G, LocalPrimeK<T> &P,
                                               FiniteField const &F) {
  NumberField<T> const &K = P.K;
  size_t n = G.size();
  for (auto const &c : FqAllProjectivePoints(F, static_cast<int>(n))) {
    std::vector<ElemE<T>> y(n, ElemE<T>(K.deg, T(0)));
    // y G: the products B(y, b_j)
    std::vector<ElemE<T>> cl;
    for (auto const &x : c) {
      cl.push_back(LiftResidue(x, K));
    }
    bool ok = true;
    ElemE<T> Qy(K.deg, T(0));
    for (size_t j = 0; j < n && ok; j++) {
      ElemE<T> bj(K.deg, T(0));
      for (size_t i = 0; i < n; i++) {
        if (!cl[i].is_zero()) {
          bj += EMul(cl[i], G[i][j], K);
        }
      }
      if (!bj.is_zero() && LocalVal(bj, P) < 1 - P.e) {
        ok = false;
      }
      if (!cl[j].is_zero()) {
        Qy += EMul(bj, cl[j], K);
      }
    }
    if (!ok) {
      continue;
    }
    if (!Qy.is_zero() && LocalVal(Qy, P) < 2) {
      continue;
    }
    return SuperlatticeGram(G, c, P);
  }
  return {};
}

// Hecke: _local_factor_unimodular, G unimodular of integral norm
template <typename T>
inline T LocalFactorUnimodularDyadic(KMat<T> const &G, LocalPrimeK<T> &P,
                                             FiniteField const &F) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  NumberField<T> const &K = P.K;
  QuadKirschmerSymbol<T> S = KirschmerSymbolDyadic(G, P, F);
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
  if (S.sL.size() != 1 || S.sL[0] != 0) {
    std::cerr << "ABQUAD: the lattice is not unimodular\n";
    throw TerminalException{1};
  }
#endif
  int d = S.dims[0];
  int b = S.wL[0];
  int a = S.uL[0];
  T q = UniversalScalarConversion<T, typename underlying_z_ring<T>::ring_type>(P.qf);
  int e = P.e;
  auto qp = [&](T const &x) -> T {
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
    if (GetDenominator_z(x) != 1) {
      std::cerr << "ABQUAD: a non integral exponent\n";
      throw TerminalException{1};
    }
#endif
    return QPowRat(q, UniversalScalarConversion<int64_t, Tint>(GetNumerator_z(x)));
  };
  auto eps = [&](int c) -> int { return (c == QD_INF) ? 1 : -1; };
  T lf;
  if (d % 2 == 0) {
    int r = d / 2;
    ElemE<T> disc = S.dets[0];
    if (r % 2 == 1) {
      disc = -disc;
    }
    int c = QuadDefect(disc, P, F);
    T rh = T(r) - T(1) / T(2);
    if (d == 2) {
      if (a < b && b == 2 && c == 2 * e) {
        lf = qp((e - a - 1) / 2) * (q - eps(c));
      } else if (b == e && a + e + 1 <= c && c < 2 * e) {
        lf = 2 * qp((c - e - a) / 2);
      } else {
        lf = 1;
      }
    } else if ((a + b) % 2 == 0) {
      if (e == a) {
        lf = 1;
      } else if (c >= 2 * e) {
        lf = qp(T(e - a) * rh - r) * (QPowRat(q, r) - eps(c));
      } else if (a + e + 1 <= c) {
        lf = 2 * qp(T(c - e - a - 1) * rh);
      } else {
        lf = qp(T((c - e - a - 1) * (r - 1)));
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
          lf = qp(T(e - a - 1) * rh) * (QPowRat(q, r) - ee) *
               (QPowRat(q, r - 1) + ee) / 2;
        } else if (e == b) {
          lf = qp(T(e - a - 1) * rh) * T(q + 1);
        } else if (cc == QD_INF && c == cc) {
          int ee = eps(c);
          lf = qp(T(2 * e - b - a - 3) * rh + r) * (QPowRat(q, r) - ee) *
               (QPowRat(q, 2 * r - 2) - ee) / 2;
        } else if (cc == QD_INF) {
          // as in Hecke: (2 e - b - a - 3 (r - 1/2) + r)
          lf = qp(T(2 * e - b - a) - 3 * rh + r) *
               T(q + 1);
        } else if (c == cc) {
          int ee = eps(c);
          lf = qp(T(2 * e - b - a - 3) * rh + r) *
               (QPowRat(q, r - 1) + 1) * (QPowRat(q, r) - ee) *
               (QPowRat(q, r - 1) + ee) / 2;
        } else {
          lf = qp(T(2 * e - b - a - 3) * rh + r) *
               (QPowRat(q, r - 1) + 1) * T(q + 1);
        }
      } else if (b == e) {
        lf = 2 * qp(T(c - e - a) * rh);
      } else {
        lf = (QPowRat(q, 2 * r - 2) - 1) * qp(T(c - a - b - 2) * rh + 1);
      }
    }
  } else {
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
    if (a != 0) {
      std::cerr << "ABQUAD: odd unimodular with a != 0\n";
      throw TerminalException{1};
    }
#endif
    int r = (d - 3) / 2;
    int w = WittInvariantDyadic(G, P, F);
    if (e == b) {
      lf = 1;
    } else if (e % 2 == 1) {
      lf = (QPowRat(q, r + 1) - w) * QPowRat(q, (r + 1) * (e - b - 1));
    } else {
      lf = (w == 1 ? T((QPowRat(q, 2 * r + 2) - 1) / 2)
                   : T(q + 1)) *
           QPowRat(q, (r + 1) * (e - b - 1));
    }
  }
  ElemE<T> two = EOne(K) * T(2);
  return LocalFactorMaximalDyadic(ScaleKMat(G, two, K), P, F) * lf;
}

template <typename T>
inline T LocalFactorDyadicRamified(KMat<T> const &H, LocalPrimeK<T> &P) {
  NumberField<T> const &K = P.K;
  FiniteField F = MakeFiniteField(static_cast<T_small_prime>(P.q), P.g);
  int m = static_cast<int>(H.size());
  if (m == 1) {
    return T(1);
  }
  KMat<T> L = ScaleKMat(H, PiPower(-NormValK(H, P), P), K);
  std::optional<KMat<T>> over = IntegralOverlattice(L, P, F);
  if (!over) {
    return LocalFactorMaximalDyadic(L, P, F);
  }
  int sv = ScaleValK(H, P);
  if (sv * m == LocalVal(DetK(H, K), P)) {
    return LocalFactorUnimodularDyadic(ScaleKMat(H, PiPower(-sv, P), K), P, F);
  }
  std::vector<KMat<T>> chain{L};
  while (over) {
    chain.push_back(*over);
    over = IntegralOverlattice(chain.back(), P, F);
  }
  T f = LocalFactorMaximalDyadic(L, P, F);
  for (size_t i = 0; i + 1 < chain.size(); i++) {
    QuadKirschmerSymbol<T> Si = KirschmerSymbolDyadic(chain[i], P, F);
    QuadKirschmerSymbol<T> Sn = KirschmerSymbolDyadic(chain[i + 1], P, F);
    int64_t n_sub = 0, n_super = 0;
    for (auto const &lam : FqAllProjectivePoints(F, m)) {
      KMat<T> M = SublatticeGram(chain[i + 1], lam, P);
      if (IsLocallyIsometricDyadic(Si, KirschmerSymbolDyadic(M, P, F), P, F)) {
        n_sub++;
      }
    }
    for (auto const &c : FqAllProjectivePoints(F, m)) {
      KMat<T> M = SuperlatticeGram(chain[i], c, P);
      if (IsLocallyIsometricDyadic(Sn, KirschmerSymbolDyadic(M, P, F), P, F)) {
        n_super++;
      }
    }
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
    if (n_super == 0) {
      std::cerr << "ABQUAD: no superlattice isometric to the next in the chain\n";
      throw TerminalException{1};
    }
#endif
    f *= UniversalScalarConversion<T, int64_t>(n_sub) / UniversalScalarConversion<T, int64_t>(n_super);
  }
  return f;
}

/*
  A common index divisor q (no generator of O_K with q prime to its index)
  unramified in K, with primes of residue degree 1: the idempotents of the
  primes over q are refined from several elements of O_K
  (PrimeIdempotentsByRefinement), and K_p = Q_q, so the local lattice at p
  is the Z_q-lattice eps_p X with the trace form, handled over Q. Other
  common index divisors are not computed.
 */
template <typename T>
inline std::optional<T>
QuadLocalFactorsByIdempotents(MyMatrix<T> const &B,
                              std::vector<MyMatrix<T>> const &Omega,
                              int64_t q, PrimeDecompositionData const &dK) {
  if (dK.e != 1 || dK.f != 1) {
    return {};
  }
  std::vector<MyMatrix<T>> idem =
      PrimeIdempotentsByRefinement(B, Omega, q, dK.g);
  if (static_cast<int64_t>(idem.size()) != dK.g) {
    return {};
  }
  // Q as Q(theta), theta = 0, with the prime (q, x)
  LocalGenerator<T> lgQ;
  lgQ.S = ZeroMatrix<T>(1, 1);
  lgQ.f = {0, 1};
  lgQ.conj = {T(0)};
  T_small_prime q_u = static_cast<T_small_prime>(q);
  T res(1);
  for (auto const &E : idem) {
    // the rows of E independent modulo q: a Z_q-basis of eps_p X
    ModMatrix Emod = ModMatrixOfMatrix(E, q_u);
    ModMatrix acc(0, E.cols());
    std::vector<int> chosen;
    for (int i = 0; i < E.rows(); i++) {
      ModMatrix test(acc.rows() + 1, E.cols());
      test.topRows(acc.rows()) = acc;
      test.row(acc.rows()) = Emod.row(i);
      if (RankModPrime(test, q_u) == static_cast<size_t>(test.rows())) {
        chosen.push_back(i);
        acc = test;
      }
    }
    int m = static_cast<int>(chosen.size());
    if (m == 1) {
      continue;
    }
    MyMatrix<T> R(m, E.cols());
    for (int k = 0; k < m; k++) {
      for (int j = 0; j < E.cols(); j++) {
        R(k, j) = E(chosen[k], j);
      }
    }
    MyMatrix<T> Gq = R * B * R.transpose();
    std::vector<std::vector<ElemE<T>>> H(m, std::vector<ElemE<T>>(m));
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        H[i][j] = ElemE<T>{Gq(i, j)};
      }
    }
    LocalPrimeK<T> P = MakeLocalPrimeK(lgQ, PolyQ::X(static_cast<T_small_prime>(q)), 1, q);
    res *= (q == 2) ? LocalFactorCho(H, P) : LocalFactorOddGram(H, P);
  }
  return res;
}

// The local factors at a common index divisor q, at each prime over q:
// LocalFactorOddGram at an odd q, Cho's at q = 2 unramified,
// LocalFactorDyadicRamified at q = 2 ramified.
template <typename T>
inline std::optional<T>
QuadLocalFactorsByPrimeGenerators(MyMatrix<T> const &B,
                                  std::vector<MyMatrix<T>> const &Omega,
                                  int64_t q, PrimeDecompositionData const &dK) {
  int m = B.rows() / static_cast<int>(Omega.size());
  std::optional<std::vector<PrimeLocalGenerator<T>>> pgs =
      PrimeLocalGenerators(B, Omega, q, dK, false);
  if (!pgs) {
    return {};
  }
  T res(1);
  if (m == 1) {
    return res;
  }
  for (auto const &pg : *pgs) {
    std::vector<std::vector<ElemE<T>>> H =
        LocalOKGramFromVectors(B, pg.lg, pg.basis);
    LocalPrimeK<T> P = MakeLocalPrimeK(pg.lg, pg.g, pg.e, q);
    if (q != 2) {
      res *= LocalFactorOddGram(H, P);
    } else if (pg.e == 1) {
      res *= LocalFactorCho(H, P);
    } else {
      res *= LocalFactorDyadicRamified(H, P);
    }
  }
  return res;
}

/*
  The product of the local factors at the primes of K over q, or nothing
  when not computed: the odd primes, and the dyadic primes unramified over
  2 (Cho) and those ramified over 2.
 */
template <typename T>
inline std::optional<T>
QuadLocalFactorAtPrime(MyMatrix<T> const &B,
                       std::vector<MyMatrix<T>> const &Omega,
                       int64_t q, PrimeDecompositionData const &dK) {
  std::optional<LocalGenerator<T>> lg = FindLocalGenerator(B, Omega, q);
  if (!lg) {
    if (dK.e == 1 && dK.f == 1) {
      return QuadLocalFactorsByIdempotents(B, Omega, q, dK);
    }
    return QuadLocalFactorsByPrimeGenerators(B, Omega, q, dK);
  }
  std::vector<std::pair<PolyQ, int>> fac =
      FactorModQWithMultiplicity(PolyReduceSigned(lg->f, q));
  int m = B.rows() / (static_cast<int>(lg->f.size()) - 1);
  if (m == 1) {
    return T(1);
  }
  std::vector<std::vector<ElemE<T>>> H = LocalOKGram(B, *lg, q);
  T res(1);
  for (auto const &[g, e] : fac) {
    LocalPrimeK<T> P = MakeLocalPrimeK(*lg, g, e, q);
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
template <typename T>
inline std::map<int64_t, std::optional<T>>
QuadraticLocalFactors(MyMatrix<T> const &B,
                      std::vector<MyMatrix<T>> const &Omega,
                      AbelianField const &F) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int deg = static_cast<int>(Omega.size());
  int m = B.rows() / deg;
  DirichletGroup G = ComputeDirichletGroup(F.n);
  std::vector<std::vector<int64_t>> XK = FieldCharacters(G, F.H, true);
  // a rescaling with phi integral: b(X, X) in Z gives phi in the inverse
  // different
  Tint den(1);
  for (int i = 0; i < B.rows(); i++) {
    for (int j = 0; j < B.cols(); j++) {
      den = LCMpair(den, GetDenominator_z(B(i, j)));
    }
  }
  Tint dK(1);
  std::set<int64_t> primes{2};
  {
    int64_t x = F.n;
    for (int64_t p = 2; p <= x; p++) {
      if (x % p == 0) {
        int v = DiscriminantValuationCharacters(G, XK, p);
        dK *= MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), v);
        while (x % p == 0) {
          x /= p;
        }
      }
    }
  }
  T detB = DeterminantMat(B) * T(den * dK);
  for (int i = 1; i < B.rows(); i++) {
    detB *= T(den * dK);
  }
  Tint num = GetNumerator_z(detB);
  if (num < 0) {
    num = -num;
  }
  if (num > 1) {
    for (Tint const &p : FactorsInt(num)) {
      if (p > Tint(1L << 31)) {
        std::cerr << "ABQUAD: the prime " << p << " is too large\n";
        throw TerminalException{1};
      }
      primes.insert(UniversalScalarConversion<int64_t, Tint>(p));
    }
  }
  std::map<int64_t, std::optional<T>> res;
  for (int64_t q : primes) {
    if (q != 2) {
      // the sum over the primes p of K over q of f_p v_p(det phi)
      int c = RationalValuation(detB, q) -
              m * DiscriminantValuationCharacters(G, XK, q);
      if (c == 0) {
        continue;
      }
    }
    res[q] = QuadLocalFactorAtPrime(B, Omega, q,
                                    PrimeDecompositionCharacters(G, XK, q));
  }
  return res;
}

template <typename T>
inline std::optional<ElemE<T>> SqrtInK(ElemE<T> const &d, NumberField<T> const &K) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int n = K.deg;
  if (d.is_zero()) {
    return d;
  }
  auto [dz, den] = ClearDenominators(d);
  // dd = d den^2, integral
  ElemE<T> dd(n, T(0));
  for (int i = 0; i < n; i++) {
    Tint prod = dz[i] * den;
    dd[i] = UniversalScalarConversion<T, Tint>(prod);
  }
  MyMatrix<T> Tm = PowerTraceMatrix<T>(K.Phi);
  T discf = DeterminantMat(Tm);
  Tint D = T_abs(GetNumerator_z(discf));
  // the precision, in bits, from the sizes of the data
  size_t bits_c = 1;
  for (auto const &x : dd) {
    bits_c = std::max(bits_c, BitLength(GetNumerator_z(x)));
  }
  size_t bits_f = 1;
  for (auto const &x : K.Phi) {
    bits_f = std::max(bits_f, BitLength(UniversalScalarConversion<Tint, int64_t>(x)));
  }
  size_t need = BitLength(D) + bits_c + 2 * n * (bits_f + 1) + 4 * BitLength(Tint(n + 1)) + 64;
  size_t prec = std::max<size_t>(256, 2 * need);
  T eps(1);
  for (size_t i = 0; i < prec; i++) {
    eps /= 2;
  }
  std::vector<T> poly = UniversalStdVectorScalarConversion<T, int64_t>(K.Phi);
  std::vector<T> th = real_roots_approximation(poly, eps);
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
  if (static_cast<int>(th.size()) != n) {
    std::cerr << "ABQUAD: the field is not totally real\n";
    throw TerminalException{1};
  }
#endif
  std::vector<T> sq;
  for (int j = 0; j < n; j++) {
    T v(0);
    for (int i = n - 1; i >= 0; i--) {
      v *= th[j];
      v += dd[i];
    }
    if (v <= 0) {
      return {};
    }
    sq.push_back(square_root_approximation(v, eps));
  }
  // the inverse of the Vandermonde matrix V(j, i) = th_j^i
  MyMatrix<T> V(n, n);
  for (int j = 0; j < n; j++) {
    T pw(1);
    for (int i = 0; i < n; i++) {
      V(j, i) = pw;
      pw *= th[j];
    }
  }
  MyMatrix<T> Vinv = Inverse(V);
  T DT = UniversalScalarConversion<T, Tint>(D);
  for (uint64_t pat = 0; pat < (uint64_t(1) << (n - 1)); pat++) {
    ElemE<T> y(n, T(0));
    for (int i = 0; i < n; i++) {
      T c(0);
      for (int j = 0; j < n; j++) {
        bool neg = (j > 0) && ((pat >> (j - 1)) & 1);
        T term = Vinv(i, j) * sq[j];
        if (neg) {
          c -= term;
        } else {
          c += term;
        }
      }
      T cD = c * DT;
      Tint zi = UniversalNearestScalarInteger<Tint, T>(cD);
      T yi = UniversalScalarConversion<T, Tint>(zi);
      yi /= DT;
      y[i] = yi;
    }
    if (EMul(y, y, K) == dd) {
      T inv_den = T(1) / UniversalScalarConversion<T, Tint>(den);
      return y * inv_den;
    }
  }
  return {};
}

// The field (n, H) seen in Q(zeta_N), n | N.
inline AbelianField LiftAbelianField(AbelianField const &F, int64_t N) {
  std::set<int64_t> Hs(F.H.begin(), F.H.end());
  std::vector<int64_t> H;
  for (int64_t a = 1; a < N; a++) {
    if (std::gcd(a, N) == 1 && Hs.count(a % F.n) == 1) {
      H.push_back(a);
    }
  }
  return {N, H};
}

/*
  E = K(sqrt d) abelian over Q but not of the form K(sqrt D), D rational.

  - The automorphisms of K = Q(theta) (totally real, Galois) are the roots
    g(theta) of f in K. Writing g(theta_1) = theta_j at the real roots,
    the integer relation between 1, theta_1, ..., theta_1^{n-1} and
    theta_j is found by LLL on the Gram matrix I + W^2 v v^T at a working
    precision, and each candidate is checked exactly by f(g) = 0 in K.
  - E is Galois over Q iff sigma(d) / d is a square s_sigma^2 in K for
    every sigma, the lifts sending sqrt d to s_sigma sqrt d; it is then
    abelian iff sigma(s_tau) s_sigma = tau(s_sigma) s_tau for all sigma,
    tau. Both are exact (SqrtInK).
  - An abelian E has a conductor dividing N = lcm(n, l^{max(v_l(n), 1)} for
    the odd primes l dividing 2 n N(d), 2^{max(v_2(n) + 1, 3)}): a
    character chi of E has chi^2 a character of K, and at an odd l the
    square of a character of conductor l^k, k >= 2, keeps that conductor.
    So E is the field fixed by one of the subgroups of index 2 of the lift
    H_N of H to (Z/N)^*, containing -1 iff E is real. A prime l splitting
    completely in K (l mod N in H_N) is in the subgroup of E iff d is a
    square modulo a prime of K over l; such primes remove the other
    candidates one by one, all the computations being exact, and since E
    is among them the one left is E.
 */
template <typename T>
ElemE<T> EApplyAut(ElemE<T> const &x, ElemE<T> const &g,
                   NumberField<T> const &K) {
  ElemE<T> res(K.deg, T(0));
  for (int i = K.deg - 1; i >= 0; i--) {
    res = EMul(res, g, K);
    res[0] += x[i];
  }
  return res;
}

template <typename T>
std::vector<ElemE<T>> AutomorphismsTotallyReal(NumberField<T> const &K) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  int n = K.deg;
  auto is_root = [&](ElemE<T> const &g) -> bool {
    ElemE<T> v(n, T(0));
    for (int i = static_cast<int>(K.Phi.size()) - 1; i >= 0; i--) {
      v = EMul(v, g, K);
      v[0] += UniversalScalarConversion<T, int64_t>(K.Phi[i]);
    }
    return v.is_zero();
  };
  std::vector<T> poly = UniversalStdVectorScalarConversion<T, int64_t>(K.Phi);
  for (size_t prec : {256, 512, 1024, 2048, 4096}) {
    T eps(1);
    for (size_t i = 0; i < prec; i++) {
      eps /= 2;
    }
    std::vector<T> th = real_roots_approximation(poly, eps);
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
    if (static_cast<int>(th.size()) != n) {
      std::cerr << "ABQUAD: the field is not totally real\n";
      throw TerminalException{1};
    }
#endif
    T W2 = MyPow<T>(2, prec);
    std::vector<ElemE<T>> res;
    for (int j = 0; j < n; j++) {
      std::vector<T> v(n + 1);
      T pw(1);
      for (int k = 0; k < n; k++) {
        v[k] = pw;
        pw *= th[0];
      }
      v[n] = th[j];
      MyMatrix<T> G(n + 1, n + 1);
      for (int a = 0; a <= n; a++) {
        for (int b = 0; b <= n; b++) {
          T val = W2 * v[a];
          val *= v[b];
          if (a == b) {
            val += 1;
          }
          G(a, b) = val;
        }
      }
      LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(G, std::cerr);
      for (int r = 0; r <= n; r++) {
        Tint w = red.Pmat(r, n);
        if (w == 0) {
          continue;
        }
        ElemE<T> g(n, T(0));
        T wT = UniversalScalarConversion<T, Tint>(w);
        for (int k = 0; k < n; k++) {
          T c = UniversalScalarConversion<T, Tint>(red.Pmat(r, k));
          c /= wT;
          g[k] = -c;
        }
        if (is_root(g)) {
          if (std::find(res.begin(), res.end(), g) == res.end()) {
            res.push_back(g);
          }
          break;
        }
      }
    }
    if (static_cast<int>(res.size()) == n) {
      return res;
    }
  }
  std::cerr << "ABQUAD: the automorphisms of K are not found\n";
  throw TerminalException{1};
}

template <typename T>
bool IsAbelianQuadraticExtension(ElemE<T> const &d, NumberField<T> const &K,
                                 std::vector<ElemE<T>> const &auts) {
  ElemE<T> dinv = EInv(d, K);
  std::vector<ElemE<T>> sv;
  for (auto const &g : auts) {
    std::optional<ElemE<T>> opt =
        SqrtInK(EMul(EApplyAut(d, g, K), dinv, K), K);
    if (!opt) {
      return false;
    }
    sv.push_back(*opt);
  }
  for (size_t a = 0; a < auts.size(); a++) {
    for (size_t b = a + 1; b < auts.size(); b++) {
      ElemE<T> lhs = EMul(EApplyAut(sv[b], auts[a], K), sv[a], K);
      ElemE<T> rhs = EMul(EApplyAut(sv[a], auts[b], K), sv[b], K);
      if (lhs != rhs) {
        return false;
      }
    }
  }
  return true;
}

// The subgroups of index 2 of the finite subgroup Hl of (Z/N)^*.
inline std::vector<std::vector<int64_t>>
IndexTwoSubgroups(std::vector<int64_t> const &Hl, int64_t N) {
  auto closure = [&](std::vector<int64_t> const &gens) -> std::set<int64_t> {
    std::set<int64_t> S{1 % N};
    std::vector<int64_t> queue{1 % N};
    while (!queue.empty()) {
      int64_t a = queue.back();
      queue.pop_back();
      for (int64_t g : gens) {
        int64_t b = (a * g) % N;
        if (S.insert(b).second) {
          queue.push_back(b);
        }
      }
    }
    return S;
  };
  std::vector<int64_t> squares;
  for (int64_t h : Hl) {
    squares.push_back((h * h) % N);
  }
  // a basis of Hl modulo the squares
  std::vector<int64_t> basis;
  std::vector<int64_t> gens = squares;
  std::set<int64_t> cur = closure(gens);
  for (int64_t h : Hl) {
    if (cur.count(h) == 0) {
      basis.push_back(h);
      gens.push_back(h);
      cur = closure(gens);
    }
  }
  size_t s = basis.size();
  // the coordinates of each element over F_2
  std::set<int64_t> Ssq = closure(squares);
  std::map<int64_t, uint64_t> coord;
  for (uint64_t c = 0; c < (uint64_t(1) << s); c++) {
    int64_t rep = 1 % N;
    for (size_t i = 0; i < s; i++) {
      if ((c >> i) & 1) {
        rep = (rep * basis[i]) % N;
      }
    }
    for (int64_t x : Ssq) {
      coord[(rep * x) % N] = c;
    }
  }
  std::vector<std::vector<int64_t>> res;
  for (uint64_t lam = 1; lam < (uint64_t(1) << s); lam++) {
    std::vector<int64_t> Hp;
    for (int64_t h : Hl) {
      if (__builtin_popcountll(coord[h] & lam) % 2 == 0) {
        Hp.push_back(h);
      }
    }
    res.push_back(Hp);
  }
  return res;
}

/*
  E = K(sqrt d) for the discriminant d of the quadratic lattice (B,
  Omega), when E is abelian over Q: (N, H_E). Nothing when E is not
  abelian. The conductor is bounded from the primes ps, which must
  contain those ramified in E (and those of 2 n): by default the primes
  of N(d), which may carry primes unramified in E when d was computed in
  a basis local at another prime.
 */
template <typename T>
std::optional<AbelianField> AbelianQuadraticExtensionFieldPrimes(
    MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Omega,
    AbelianField const &FK, ElemE<T> const &d, NumberField<T> const &K,
    int m, std::set<int64_t> ps) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  std::vector<ElemE<T>> auts = AutomorphismsTotallyReal(K);
  if (!IsAbelianQuadraticExtension(d, K, auts)) {
    return {};
  }
  for (Tint const &p : FactorsInt(UniversalScalarConversion<Tint, int64_t>(2 * FK.n))) {
    ps.insert(UniversalScalarConversion<int64_t, Tint>(p));
  }
  auto vp = [](int64_t x, int64_t p) -> int {
    int v = 0;
    while (x % p == 0) {
      x /= p;
      v++;
    }
    return v;
  };
  int64_t N = FK.n;
  for (int64_t p : ps) {
    int64_t pk = 1;
    int k = (p == 2) ? std::max(vp(FK.n, 2) + 1, 3) : std::max(vp(FK.n, p), 1);
    pk *= MyPow<int64_t>(p, k);
    N = std::lcm(N, pk);
  }
  AbelianField FK2 = LiftAbelianField(FK, N);
  bool real = (m % 4 == 0);
  std::vector<std::vector<int64_t>> cands;
  for (auto const &Hp : IndexTwoSubgroups(FK2.H, N)) {
    bool has_m1 = std::find(Hp.begin(), Hp.end(), N - 1) != Hp.end();
    if (has_m1 == real) {
      cands.push_back(Hp);
    }
  }
  std::set<int64_t> HN(FK2.H.begin(), FK2.H.end());
  auto is_prime = [](int64_t x) -> bool {
    for (int64_t p = 2; p * p <= x; p++) {
      if (x % p == 0) {
        return false;
      }
    }
    return x >= 2;
  };
  for (int64_t l = 3; cands.size() > 1 && l < 10000000; l += 2) {
    if (!is_prime(l) || N % l == 0 || ps.count(l) == 1 ||
        HN.count(l % N) == 0) {
      continue;
    }
    std::optional<LocalGenerator<T>> lg = FindLocalGenerator(B, Omega, l);
    if (!lg) {
      continue;
    }
    NumberField<T> Kl;
    Kl.d = 0;
    Kl.deg = K.deg;
    Kl.Phi = lg->f;
    Kl.conj = lg->conj;
    ElemE<T> dl = DetK(LocalOKGram(B, *lg, l), Kl);
    if ((m * (m - 1) / 2) % 2 == 1) {
      dl = -dl;
    }
    std::vector<std::pair<PolyQ, int>> fac =
        FactorModQWithMultiplicity(PolyReduceSigned(lg->f, l));
    LocalPrimeK<T> P = MakeLocalPrimeK(*lg, fac[0].first, 1, l);
    if (LocalVal(dl, P) != 0) {
      continue;
    }
    FqElt res = LocalResidue(dl, P);
    bool split = (FqQuadraticCharacter(P.F, res) == 1);
    std::vector<std::vector<int64_t>> next;
    for (auto const &Hp : cands) {
      bool in = std::binary_search(Hp.begin(), Hp.end(), l % N);
      if (in == split) {
        next.push_back(Hp);
      }
    }
    cands = next;
  }
#ifdef SANITY_CHECK_ABELIAN_QUADRATIC_MASS
  if (cands.size() != 1) {
    std::cerr << "ABQUAD: " << cands.size()
              << " candidates for the abelian field E\n";
    throw TerminalException{1};
  }
#endif
  return AbelianField{N, cands[0]};
}

/*
  The evaluation of L(1 - r, chi_{E/K}) when E is not abelian
  (MassQuadraticHeckeL.h): Mpfr rounds an MPFR value checked at two
  precisions, Arb moreover certifies the integer with Arb balls (proven
  error bounds, about twice the time), and Double evaluates every term of
  the sum in double precision (the kernel itself is still built with
  MPFR), which is enough when the integer to round has at most about 26
  bits and gives up otherwise. Arb needs a build with flint.
 */
enum class QuadLValueMethod { Double, Mpfr, Arb };

inline QuadLValueMethod ParseQuadLValueMethod(std::string const &s) {
  if (s == "double") {
    return QuadLValueMethod::Double;
  }
  if (s == "mpfr") {
    return QuadLValueMethod::Mpfr;
  }
  if (s == "arb") {
#ifdef ENABLE_FLINT_SUPPORT
    return QuadLValueMethod::Arb;
#else
    std::cerr << "ABQUAD: lvalue=arb needs flint, this program is built "
                 "without flint (make ENABLE_FLINT_SUPPORT=1)\n";
    throw TerminalException{1};
#endif
  }
  std::cerr << "ABQUAD: the L-value method " << s
            << " is not double, mpfr or arb\n";
  throw TerminalException{1};
}

/*
  The option lvalue=double / lvalue=mpfr / lvalue=arb of the command line
  programs, found
  anywhere among the arguments and removed from them; mpfr when absent.
 */
inline QuadLValueMethod ExtractQuadLValueOption(int &argc, char *argv[]) {
  QuadLValueMethod res = QuadLValueMethod::Mpfr;
  int j = 0;
  for (int i = 0; i < argc; i++) {
    std::string a = argv[i];
    if (a.rfind("lvalue=", 0) == 0) {
      res = ParseQuadLValueMethod(a.substr(7));
    } else {
      argv[j++] = argv[i];
    }
  }
  argc = j;
  return res;
}

// clang-format off
#endif  // SRC_GENERA_MASSABELIANQUADRATICLOCAL_H_
// clang-format on
