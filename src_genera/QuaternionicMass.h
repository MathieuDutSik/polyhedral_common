// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_QUATERNIONICMASS_H_
#define SRC_GENERA_QUATERNIONICMASS_H_

// clang-format off
#include "ResidueArithmeticLegacy.h"
#include "MultiplicityLocalFactor.h"
#include "AbelianQuadraticMass.h"
#include "SiegelSpecialValues.h"
#include "ClassicLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include "MAT_MatrixModPrime.h"
#include <algorithm>
#include <chrono>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_QUATERNIONIC_MASS
#endif

/*
  The mass of a component V_rho of a finite group G whose End_G(V_rho) is
  M_r(D), D a totally definite quaternion algebra over a totally real field
  K of degree c (Albert type III: the adjoint involution is of symplectic
  type, and the centralizer C_rho = U is an inner form of Sp_{2r} over K,
  simply connected).

  The method compares volumes. For two lattices X and S in D-Hermitian
  spaces of the same rank (all positive definite ones are isometric up to
  a K-similitude, which does not change U),

      mass(X) / mass(S) = prod_p vol(U(S_p)) / vol(U(X_p)),

  for one Haar measure carried over by the similitude. For a lattice Z
  with Lie lattice u(Z) = {Y in u : Y Z subset Z}, u the skew-adjoint
  D-linear endomorphisms, the Cayley transform Y -> (1 + Y)(1 - Y)^{-1}
  is a bijection of p^c u(Z) onto the congruence subgroup Gamma_Z(p^c) of
  U(Z_p) (c >= 1, p odd) and of p^{c-1} u(Z) onto Gamma_Z(p^c) (c >= 2,
  p = 2), with Jacobian of valuation 0 and dim u respectively, so that in
  both cases

      vol(U(Z_p)) = |U(Z_p) mod p^c| p^{-c dim u} |disc u(Z)|_p^{1/2},

  the volume of u(Z) being measured by the trace form tr(Y Y') on the
  natural representation, which the similitude preserves. |U(Z_p) mod
  p^c| is the order of the image found by LocalOrthogonalImageGenerators
  with the commutant of D, certified generators and exact order; it is
  also computed modulo p^{c+1} under SANITY_CHECK_QUATERNIONIC_MASS, and
  the ratio p^{dim u} (the Cayley transform) is checked.

  S is the principal lattice O^r, O a maximal order of D, with the form
  Tr_{K/Q}(lambda trd(sum x_i y_i^*)), lambda totally positive in K chosen
  to make S unimodular at as many primes as possible. Its mass is Prasad's
  formula (Hashimoto-Ibukiyama over Q):

      mass(S) = prod_{i=1}^r 2^{-c} |zeta_K(1 - 2i)|
                prod_{P | d(D)} prod_{i=1}^r (N P^i + (-1)^i).

  The image searches also impose, modulo p, the integral forms B(w x, y)
  for w in D (a Z-basis of their integral points: w runs over the trace
  dual of the order), which fix the Hermitian form h modulo p O where the
  form b = Tr trd h alone fixes it modulo a larger ideal (at 2 for the
  Hurwitz order, h is then only known modulo the prime above 2).

  The product of the volume ratios runs over the p dividing det X, det S
  and [O : O_X], O_X the multiplier ring of X; elsewhere X_p and S_p are
  unimodular over the maximal order, hence similar.

  X is L_rho when G spans D (V_rho free of rank r over D); when G spans
  M_s(D), s > 1, X = L cap e V_rho for a primitive self-adjoint idempotent
  e of the span of G, D acting by the e g e (FiniteGroupComponents), and
  the multiplicity factors [U(X_p) : Stab(L_p)] enter as in the split
  case. Nothing is returned beyond the bounds of the searches.
 */

struct QuaternionicMassBounds {
  // the F_p-space searched for an overorder at p
  size_t max_overorder = size_t(1) << 24;
  MultiplicityFactorBounds image;
};

/*
  A Z-basis (rows) of the integer x with x Eq = 0, for the commutant
  systems: many unknowns (k^2), many more conditions than their rank, and
  a kernel of small dimension. The integral elimination of the whole
  system suffers from the growth of the coefficients, so the kernel is
  found modulo p and lifted:
  - the reduced echelon form modulo p of the conditions independent
    modulo p (ColumnRankProfileModPrime), giving the kernel basis K with
    the identity on the free unknowns;
  - its entries as fractions by rational reconstruction, and K Eq = 0
    checked exactly. The rank over Q is at least the rank modulo p, so the
    rational kernel has dimension at most that of K, and K, of full rank,
    spans it;
  - the saturation: the integer points of span_Q(K) are the y K with y
    integral (the coordinates on the free unknowns) and y K integral, that
    is y (D K) = 0 modulo the common denominator D of K, one column at a
    time on a basis of the y found so far.
  When the reconstruction fails or does not check, the integral
  elimination of the full system, which is also the comparison under
  SANITY_CHECK_QUATERNIONIC_MASS for the small systems.
 */
template <typename Tint>
MyMatrix<Tint> IntegralKernelOfConditions(MyMatrix<Tint> const &Eq) {
  int N = Eq.rows();
  int n_cond = Eq.cols();
  T_small_prime const p = 2147483629;
  auto fallback = [&]() -> MyMatrix<Tint> {
    return NullspaceIntTrMat(MyMatrix<Tint>(Eq.transpose()));
  };
  std::vector<int> cols = ColumnRankProfileModPrime(Eq);
  int r = static_cast<int>(cols.size());
  // the reduced echelon form modulo p of the selected conditions (rows)
  ModMatrix A(r, N);
  for (int j = 0; j < r; j++) {
    for (int i = 0; i < N; i++) {
      A(j, i) = ResModOf(Eq(i, cols[j]), p);
    }
  }
  RowSpaceModP rs = ComputeRowSpaceModP(std::move(A), p);
  ModMatrix const &E = rs.basis;
  std::vector<int> pivcol(rs.pivots.begin(), rs.pivots.end());
  std::vector<bool> is_pivot(N, false);
  for (int col : pivcol) {
    is_pivot[col] = true;
  }
  std::vector<int> freecol;
  for (int u = 0; u < N; u++) {
    if (!is_pivot[u]) {
      freecol.push_back(u);
    }
  }
  int d = static_cast<int>(freecol.size());
  if (d == 0) {
    return MyMatrix<Tint>(0, N);
  }
  // the fractions n/m with |n|, m below sqrt(p/2) from a residue
  int64_t const bound = 32767;
  auto reconstruct = [&](T_residue a)
      -> std::optional<std::pair<int64_t, int64_t>> {
    int64_t r0 = ResToSigned(p), r1 = ResToSigned(a), t0 = 0, t1 = 1;
    while (r1 > bound) {
      int64_t q = r0 / r1;
      int64_t r2 = r0 - q * r1;
      int64_t t2 = t0 - q * t1;
      r0 = r1;
      r1 = r2;
      t0 = t1;
      t1 = t2;
    }
    if (t1 == 0 || t1 > bound || t1 < -bound) {
      return {};
    }
    if (t1 < 0) {
      return std::pair<int64_t, int64_t>{-r1, -t1};
    }
    return std::pair<int64_t, int64_t>{r1, t1};
  };
  // K(f, pivcol[i]) = -E[i][freecol[f]], K(f, freecol[f]) = 1
  std::vector<std::vector<std::pair<int64_t, int64_t>>> Kfrac(
      d, std::vector<std::pair<int64_t, int64_t>>(pivcol.size()));
  Tint D(1);
  for (int f = 0; f < d; f++) {
    for (size_t i = 0; i < pivcol.size(); i++) {
      T_residue a = ResNeg(E(static_cast<int>(i), freecol[f]), p);
      std::optional<std::pair<int64_t, int64_t>> opt = reconstruct(a);
      if (!opt) {
        return fallback();
      }
      Kfrac[f][i] = *opt;
      Tint den = UniversalScalarConversion<Tint, int64_t>(opt->second);
      D = LCMpair(D, den);
    }
  }
  MyMatrix<Tint> KD = ZeroMatrix<Tint>(d, N);
  for (int f = 0; f < d; f++) {
    KD(f, freecol[f]) = D;
    for (size_t i = 0; i < pivcol.size(); i++) {
      auto const &[num, den] = Kfrac[f][i];
      if (num != 0) {
        Tint num_T = UniversalScalarConversion<Tint, int64_t>(num);
        Tint den_T = UniversalScalarConversion<Tint, int64_t>(den);
        KD(f, pivcol[i]) = num_T * (D / den_T);
      }
    }
  }
  // KD Eq = 0, by the nonzero entries of the columns of Eq
  for (int c = 0; c < n_cond; c++) {
    std::vector<int> nz;
    for (int i = 0; i < N; i++) {
      if (Eq(i, c) != 0) {
        nz.push_back(i);
      }
    }
    for (int f = 0; f < d; f++) {
      Tint s(0);
      for (int i : nz) {
        s += KD(f, i) * Eq(i, c);
      }
      if (s != 0) {
        return fallback();
      }
    }
  }
  if (D == 1) {
    return KD;
  }
  // a basis B (rows) of the y in Z^d with y KD = 0 modulo D, one column
  // at a time: the values v = B c modulo D brought to (g, 0, ..., 0) by
  // unimodular row operations (Euclid), then the row of g multiplied by
  // D / g
  MyMatrix<Tint> B = IdentityMat<Tint>(d);
  for (int u = 0; u < N; u++) {
    if (!is_pivot[u]) {
      // the column of a free unknown is D times a unit vector
      continue;
    }
    std::vector<Tint> v(d);
    bool all_zero = true;
    for (int i = 0; i < d; i++) {
      Tint s(0);
      for (int f = 0; f < d; f++) {
        if (B(i, f) != 0 && KD(f, u) != 0) {
          s += B(i, f) * KD(f, u);
        }
      }
      v[i] = ResInt(s, D);
      if (v[i] != 0) {
        all_zero = false;
      }
    }
    if (all_zero) {
      continue;
    }
    int i0 = -1;
    for (int i = 0; i < d; i++) {
      if (v[i] == 0) {
        continue;
      }
      if (i0 == -1) {
        i0 = i;
        continue;
      }
      // Euclid on (v[i0], v[i]), the same operations on the rows
      while (v[i] != 0) {
        Tint q = QuoInt(v[i0], v[i]);
        for (int f = 0; f < d; f++) {
          B(i0, f) -= q * B(i, f);
        }
        v[i0] -= q * v[i];
        std::swap(v[i0], v[i]);
        B.row(i0).swap(B.row(i));
      }
    }
    Tint g = GcdPair(v[i0], D);
    Tint mult = D / g;
    for (int f = 0; f < d; f++) {
      B(i0, f) *= mult;
    }
  }
  MyMatrix<Tint> res = B * KD;
  for (int r2 = 0; r2 < res.rows(); r2++) {
    for (int u = 0; u < N; u++) {
      res(r2, u) /= D;
    }
  }
#ifdef SANITY_CHECK_QUATERNIONIC_MASS
  if (N <= 400) {
    MyMatrix<Tint> ref = fallback();
    if (ComputeRowHermiteNormalForm(ref).second !=
        ComputeRowHermiteNormalForm(res).second) {
      std::cerr << "QUAT_MASS: IntegralKernelOfConditions differs from the "
                   "integral elimination\n";
      throw TerminalException{1};
    }
  }
#endif
  return res;
}

template <typename T>
std::vector<MyMatrix<T>> QuatZBasis(std::vector<MyMatrix<T>> const &L, int k) {
  MyMatrix<T> F(static_cast<int>(L.size()), k * k);
  for (size_t i = 0; i < L.size(); i++) {
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        F(static_cast<int>(i), u * k + v) = L[i](u, v);
      }
    }
  }
  MyMatrix<T> Bz = GetZbasis(F);
  std::vector<MyMatrix<T>> res;
  for (int i = 0; i < Bz.rows(); i++) {
    MyMatrix<T> Y(k, k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Y(u, v) = Bz(i, u * k + v);
      }
    }
    res.push_back(Y);
  }
  return res;
}

// A Z-basis of span_Q(L) cap M_k(Z)
template <typename T, typename Tint>
std::vector<MyMatrix<T>> QuatIntegralPoints(std::vector<MyMatrix<T>> const &L,
                                            int k) {
  int N = k * k;
  MyMatrix<T> R(static_cast<int>(L.size()), N);
  for (size_t i = 0; i < L.size(); i++) {
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        R(static_cast<int>(i), u * k + v) = L[i](u, v);
      }
    }
  }
  // the y with R y = 0, then the integer x orthogonal to them
  MyMatrix<T> perp = NullspaceTrMat(R);
  std::vector<MyMatrix<T>> res;
  if (perp.rows() == 0) {
    for (int w = 0; w < N; w++) {
      MyMatrix<T> Y = ZeroMatrix<T>(k, k);
      Y(w / k, w % k) = T(1);
      res.push_back(Y);
    }
    return res;
  }
  MyMatrix<Tint> perp_i =
      UniversalMatrixConversion<Tint, T>(RemoveFractionMatrix(perp));
  MyMatrix<Tint> ker = NullspaceIntTrMat(perp_i);
  for (int r = 0; r < ker.rows(); r++) {
    MyMatrix<T> Y(k, k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Y(u, v) = UniversalScalarConversion<T, Tint>(ker(r, u * k + v));
      }
    }
    res.push_back(Y);
  }
  return res;
}

/*
  A Z-basis of the integral Y (columns) with Y C = C Y for the C of gens
  and, when skew, Y^T B + B Y = 0.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<T>> QuatIntegralCommutant(MyMatrix<T> const &B,
                                               std::vector<MyMatrix<T>> const &gens,
                                               bool skew) {
  int k = B.rows();
  int kk = k * k;
  int n_cond = kk * static_cast<int>(gens.size() + (skew ? 1 : 0));
  MyMatrix<Tint> Eq = ZeroMatrix<Tint>(kk, n_cond);
  int col = 0;
  for (auto const &Cr : gens) {
    MyMatrix<Tint> C = UniversalMatrixConversion<Tint, T>(Cr);
    // (Y C - C Y)(u, v) = sum_w Y(u, w) C(w, v) - C(u, w) Y(w, v)
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        for (int w = 0; w < k; w++) {
          Eq(u * k + w, col) += C(w, v);
          Eq(w * k + v, col) -= C(u, w);
        }
        col++;
      }
    }
  }
  if (skew) {
    MyMatrix<Tint> Bi = UniversalMatrixConversion<Tint, T>(B);
    // (Y^T B + B Y)(u, v) = sum_w Y(w, u) B(w, v) + B(u, w) Y(w, v)
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        for (int w = 0; w < k; w++) {
          Eq(w * k + u, col) += Bi(w, v);
          Eq(w * k + v, col) += Bi(u, w);
        }
        col++;
      }
    }
  }
  MyMatrix<Tint> ker = IntegralKernelOfConditions(Eq);
  std::vector<MyMatrix<T>> res;
  for (int r = 0; r < ker.rows(); r++) {
    MyMatrix<T> Y(k, k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Y(u, v) = UniversalScalarConversion<T, Tint>(ker(r, u * k + v));
      }
    }
    res.push_back(Y);
  }
  return res;
}

// positive definite, by the leading principal minors
template <typename T> bool QuatIsPositiveDefinite(MyMatrix<T> const &G) {
  int n = G.rows();
  for (int s = 1; s <= n; s++) {
    MyMatrix<T> M(s, s);
    for (int i = 0; i < s; i++) {
      for (int j = 0; j < s; j++) {
        M(i, j) = G(i, j);
      }
    }
    if (DeterminantMat(M) <= T(0)) {
      return false;
    }
  }
  return true;
}

template <typename T> T QuatTrace(MyMatrix<T> const &M) {
  T t(0);
  for (int i = 0; i < M.rows(); i++) {
    t += M(i, i);
  }
  return t;
}

// the discriminant of a lattice of matrices for the trace form tr(Y Y')
template <typename T> T QuatTraceDiscriminant(std::vector<MyMatrix<T>> const &L) {
  int d = static_cast<int>(L.size());
  MyMatrix<T> G(d, d);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      G(i, j) = QuatTrace(MyMatrix<T>(L[i] * L[j]));
    }
  }
  return DeterminantMat(G);
}

/*
  The decomposition (e, f, g) of p in the subfield of Q(zeta_n) fixed by
  the subgroup H of (Z/n)^*: with n = p^a n', the inertia group is the
  image of the u = 1 mod n', the decomposition group adds a u = p mod n',
  u = 1 mod p^a.
 */
inline std::tuple<int, int, int> AbelianPrimeDecomposition(int64_t n,
                                                           std::vector<int64_t> const &H,
                                                           int64_t p) {
  if (n <= 2) {
    return {1, 1, 1};
  }
  int64_t pa = 1, np = n;
  while (np % p == 0) {
    np /= p;
    pa *= p;
  }
  auto closure = [&](std::vector<int64_t> const &gens) -> size_t {
    std::set<int64_t> S{1};
    std::vector<int64_t> queue{1};
    while (!queue.empty()) {
      int64_t a = queue.back();
      queue.pop_back();
      for (int64_t g : gens) {
        int64_t b = (a * g) % n;
        if (S.insert(b).second) {
          queue.push_back(b);
        }
      }
    }
    return S.size();
  };
  std::vector<int64_t> inertia(H.begin(), H.end());
  for (int64_t u = 1; u < n; u++) {
    if (std::gcd(u, n) == 1 && (u - 1) % np == 0) {
      inertia.push_back(u);
    }
  }
  std::vector<int64_t> decomp = inertia;
  if (np > 1) {
    for (int64_t u = 1; u < n; u++) {
      if (std::gcd(u, n) == 1 && ((u - p) % np + np) % np == 0 &&
          (u - 1) % pa == 0) {
        decomp.push_back(u);
        break;
      }
    }
  }
  size_t h = H.size();
  size_t phi = static_cast<size_t>(EulerPhi(n));
  int ee = static_cast<int>(closure(inertia) / h);
  int dd = static_cast<int>(closure(decomp) / h);
  int ff = dd / ee;
  int gg = static_cast<int>(phi / h) / dd;
  return {ee, ff, gg};
}

// Tr_{K/Q}(trd(x y)) for x, y in D acting on F of dimension 4 c r
template <typename T>
T QuatTrTrd(MyMatrix<T> const &X, MyMatrix<T> const &Y, int r) {
  return QuatTrace(MyMatrix<T>(X * Y)) / ScalarFromInt64<T>(2 * r);
}

template <typename T>
T QuatOrderDisc(std::vector<MyMatrix<T>> const &O, int r) {
  int d = static_cast<int>(O.size());
  MyMatrix<T> G(d, d);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      G(i, j) = QuatTrTrd(O[i], O[j], r);
    }
  }
  return T_abs(DeterminantMat(G));
}

/*
  A maximal order containing the order O (Z-basis, matrices on F). At
  each p dividing |disc O| / d_K^4: if O_p is not maximal there is an
  order O' containing O with x in O' \ O, p x in O, and such an x lies in
  O^# (the dual for Tr trd); the x of (O^# cap p^{-1} O) / O are tried, O'
  being the ring they generate with O when it is an order (its
  discriminant stays a multiple of d_K^4).
 */
template <typename T, typename Tint>
std::optional<std::vector<MyMatrix<T>>>
QuatMaximalOrder(std::vector<MyMatrix<T>> O, int r, T const &dK,
                 QuaternionicMassBounds const &bounds,
                 [[maybe_unused]] std::ostream &os) {
  int k = O[0].rows();
  int d = static_cast<int>(O.size());
  T dK4 = dK * dK * dK * dK;
  auto t_start = std::chrono::steady_clock::now();
  (void)t_start;
  auto ring_closure = [&](std::vector<MyMatrix<T>> M)
      -> std::optional<std::vector<MyMatrix<T>>> {
    for (int iter = 0; iter < 64; iter++) {
      std::vector<MyMatrix<T>> prods = M;
      for (auto const &X1 : M) {
        for (auto const &X2 : M) {
          prods.push_back(X1 * X2);
        }
      }
      std::vector<MyMatrix<T>> next = QuatZBasis(prods, k);
      if (static_cast<int>(next.size()) != d) {
        return {};
      }
      T disc = QuatOrderDisc(next, r);
      T q = disc / dK4;
      if (!IsInteger(q) || q < T(1)) {
        return {};
      }
      if (disc == QuatOrderDisc(M, r)) {
        return next;
      }
      M = next;
    }
    return {};
  };
  while (true) {
    T q = QuatOrderDisc(O, r) / dK4;
    if (!IsInteger(q)) {
      std::cerr << "QUAT_MASS: the discriminant of the order is not a "
                   "multiple of d_K^4\n";
      throw TerminalException{1};
    }
    Tint q_i = UniversalScalarConversion<Tint, T>(q);
    std::vector<Tint> primes;
    if (q_i > 1) {
      primes = FactorsInt(q_i);
      std::sort(primes.begin(), primes.end());
      primes.erase(std::unique(primes.begin(), primes.end()), primes.end());
    }
    bool enlarged = false;
    // the coordinates of an element of D in the basis of O
    MyMatrix<T> Ob(d, k * k);
    for (int a = 0; a < d; a++) {
      for (int u = 0; u < k; u++) {
        for (int v = 0; v < k; v++) {
          Ob(a, u * k + v) = O[a](u, v);
        }
      }
    }
    auto coordsO = [&](MyMatrix<T> const &Y) -> MyVector<T> {
      MyVector<T> y(k * k);
      for (int u = 0; u < k; u++) {
        for (int v = 0; v < k; v++) {
          y(u * k + v) = Y(u, v);
        }
      }
      std::optional<MyVector<T>> sol = SolutionMat(Ob, y);
      if (!sol) {
        std::cerr << "QUAT_MASS: an element is not in the span of the order\n";
        throw TerminalException{1};
      }
      return *sol;
    };
    /*
      First the left and right orders of the two-sided ideal I = O cap p
      O^# (in coordinates the v with G v = 0 mod p, plus p Z^d, G the Gram
      matrix of Tr trd): orders containing O, found by linear algebra (the
      c with M_b c integral for the elements b of a basis of I is R^{-1}
      Z^d, R a Z-basis of the rows of the M_b). When neither is larger
      than O at p, the elements of (O^# cap p^{-1} O) / O are tried.
     */
#ifdef DEBUG_QUATERNIONIC_MASS
    os << "QUAT_MASS: time " << "loop start" << " " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count() << "\n";
#endif
    for (Tint const &p_i : primes) {
      int64_t p = UniversalScalarConversion<int64_t, Tint>(p_i);
      MyMatrix<T> Gt(d, d);
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          Gt(i, j) = QuatTrTrd(O[i], O[j], r);
        }
      }
      T_small_prime p_u = static_cast<T_small_prime>(p);
      RowSpaceModP rsI = ComputeRowSpaceModP(ModMatrixOfMatrix(Gt, p_u), p_u);
      // the kernel vectors, one per free column in increasing order
      ModMatrix kerI = FreeKernelBasisModP(rsI, d, p_u);
      std::vector<bool> pivI(d, false);
      for (size_t pv : rsI.pivots) {
        pivI[pv] = true;
      }
      std::vector<MyVector<T>> Igen;
      int i_kerI = 0;
      for (int f = 0; f < d; f++) {
        if (!pivI[f]) {
          MyVector<T> v(d);
          for (int a = 0; a < d; a++) {
            v(a) = UniversalScalarConversion<T, T_residue>(kerI(i_kerI, a));
          }
          i_kerI++;
          Igen.push_back(v);
        }
        MyVector<T> w = ZeroVector<T>(d);
        w(f) = ScalarFromInt64<T>(p);
        Igen.push_back(w);
      }
      MyMatrix<T> Iall(static_cast<int>(Igen.size()), d);
      for (size_t i2 = 0; i2 < Igen.size(); i2++) {
        for (int a = 0; a < d; a++) {
          Iall(static_cast<int>(i2), a) = Igen[i2](a);
        }
      }
      MyMatrix<T> IB = GetZbasis(Iall);
      if (IB.rows() == d && T_abs(DeterminantMat(IB)) > T(1)) {
        MyMatrix<T> IBinv = Inverse(IB);
        std::vector<MyMatrix<T>> Ielts;
        for (int i2 = 0; i2 < d; i2++) {
          MyMatrix<T> Y = ZeroMatrix<T>(k, k);
          for (int a = 0; a < d; a++) {
            Y += IB(i2, a) * O[a];
          }
          Ielts.push_back(Y);
        }
        for (int side = 0; side < 2 && !enlarged; side++) {
          // the rows of the M_b: the I-coordinates of w_a b (left) or b w_a
          MyMatrix<T> Rows(d * d, d);
          for (int ib = 0; ib < d; ib++) {
            for (int a = 0; a < d; a++) {
              MyMatrix<T> Y = (side == 0) ? MyMatrix<T>(O[a] * Ielts[ib])
                                          : MyMatrix<T>(Ielts[ib] * O[a]);
              MyVector<T> cI = IBinv.transpose() * coordsO(Y);
              for (int t = 0; t < d; t++) {
                Rows(ib * d + t, a) = cI(t);
              }
            }
          }
          MyMatrix<T> R = GetZbasis(Rows);
          if (R.rows() != d) {
            continue;
          }
          MyMatrix<T> Rinv = Inverse(R);
          if (T_abs(DeterminantMat(Rinv)) <= T(1)) {
            continue;
          }
          std::vector<MyMatrix<T>> On;
          for (int col = 0; col < d; col++) {
            MyMatrix<T> Y = ZeroMatrix<T>(k, k);
            for (int a = 0; a < d; a++) {
              Y += Rinv(a, col) * O[a];
            }
            On.push_back(Y);
          }
#ifdef DEBUG_QUATERNIONIC_MASS
          os << "QUAT_MASS: the " << (side == 0 ? "left" : "right")
             << " order of O cap p O^# at p=" << p << " has index "
             << T_abs(DeterminantMat(Rinv)) << " over O\n";
#endif
          std::optional<std::vector<MyMatrix<T>>> opt =
              ring_closure(QuatZBasis(On, k));
#ifdef DEBUG_QUATERNIONIC_MASS
          os << "QUAT_MASS: its ring closure " << (opt ? "is an order" : "fails")
             << "\n";
#endif
          if (opt && QuatOrderDisc(*opt, r) < QuatOrderDisc(O, r)) {
#ifdef DEBUG_QUATERNIONIC_MASS
            os << "QUAT_MASS: the " << (side == 0 ? "left" : "right")
               << " order of O cap p O^# at p=" << p << " enlarges O\n";
#endif
            O = *opt;
            enlarged = true;
          }
        }
      }
      if (enlarged) {
        break;
      }
    }
    if (enlarged) {
      continue;
    }
#ifdef DEBUG_QUATERNIONIC_MASS
    os << "QUAT_MASS: time " << "idealizer done" << " " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count() << "\n";
#endif
    for (Tint const &p_i : primes) {
      int64_t p = UniversalScalarConversion<int64_t, Tint>(p_i);
      // whether O contains O_K: its integral center has discriminant d_K
      // (needed by the filter y^2 = 0 mod p O, x^2 = trd(x) x - nrd(x), and
      // by the restriction to the y with a y in p O below)
      bool has_OK = false;
      std::vector<MyMatrix<T>> Zs;
      {
        MyMatrix<T> Eq(d, k * k * d);
        for (int i = 0; i < d; i++) {
          for (int j = 0; j < d; j++) {
            MyMatrix<T> C = O[i] * O[j] - O[j] * O[i];
            for (int u = 0; u < k; u++) {
              for (int v = 0; v < k; v++) {
                Eq(i, j * k * k + u * k + v) = C(u, v);
              }
            }
          }
        }
        MyMatrix<Tint> Eq_i =
            UniversalMatrixConversion<Tint, T>(RemoveFractionMatrix(Eq));
        MyMatrix<Tint> Zc = NullspaceIntMat(Eq_i);
        int cc = Zc.rows();
        MyMatrix<T> Gz(cc, cc);
        for (int t = 0; t < cc; t++) {
          MyMatrix<T> Z = ZeroMatrix<T>(k, k);
          for (int i = 0; i < d; i++) {
            Z += UniversalScalarConversion<T, Tint>(Zc(t, i)) * O[i];
          }
          Zs.push_back(Z);
        }
        for (int t1 = 0; t1 < cc; t1++) {
          for (int t2 = 0; t2 < cc; t2++) {
            Gz(t1, t2) = QuatTrace(MyMatrix<T>(Zs[t1] * Zs[t2])) /
                         ScalarFromInt64<T>(k / cc);
          }
        }
        has_OK = (4 * r * cc == k) && (T_abs(DeterminantMat(Gz)) == dK);
      }
      // the x = sum a_i w_i / p in O^#: sum_i a_i Tr trd(w_i w_j) = 0 mod p
      T_small_prime p_u = static_cast<T_small_prime>(p);
      MyMatrix<T> GmT(d, d);
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          GmT(j, i) = QuatTrTrd(O[i], O[j], r);
        }
      }
      ModMatrix Gm = ModMatrixOfMatrix(GmT, p_u);
      /*
        With O_K in O: an overorder O' at p has an x outside O with a x in
        O, a the product of the primes of K over p (O'/O is an O_K-module
        killed by a power of p), so y = p x satisfies a y in p O. a / p O_K
        is the nilradical of O_K / p O_K, the kernel of the F_p-linear map
        z -> z^{p^m}, p^m >= [K:Q]. The space of the y then has dimension
        4 sum_P f_P over F_p instead of 4 [K:Q].
       */
      if (has_OK) {
        int cc = static_cast<int>(Zs.size());
        MyMatrix<T> Zb(cc, k * k);
        for (int t = 0; t < cc; t++) {
          for (int u = 0; u < k; u++) {
            for (int v = 0; v < k; v++) {
              Zb(t, u * k + v) = Zs[t](u, v);
            }
          }
        }
        auto coordsZ = [&](MyMatrix<T> const &Y) -> MyVector<T> {
          MyVector<T> y(k * k);
          for (int u = 0; u < k; u++) {
            for (int v = 0; v < k; v++) {
              y(u * k + v) = Y(u, v);
            }
          }
          std::optional<MyVector<T>> sol = SolutionMat(Zb, y);
          if (!sol) {
            std::cerr << "QUAT_MASS: a power of a central element is not "
                         "central\n";
            throw TerminalException{1};
          }
          return *sol;
        };
        int64_t pm = p;
        while (pm < cc) {
          pm *= p;
        }
        // the Frobenius z -> z^{p^m} on O_K / p O_K, columns
        MyMatrix<T> FrT(cc, cc);
        for (int t = 0; t < cc; t++) {
          MyMatrix<T> Zp = IdentityMat<T>(k);
          for (int64_t e2 = 0; e2 < pm; e2++) {
            Zp = Zp * Zs[t];
          }
          MyVector<T> cz = coordsZ(Zp);
          for (int u = 0; u < cc; u++) {
            FrT(u, t) = cz(u);
          }
        }
        RowSpaceModP rf =
            ComputeRowSpaceModP(ModMatrixOfMatrix(FrT, p_u), p_u);
        ModMatrix kerF = FreeKernelBasisModP(rf, cc, p_u);
        for (int i_f = 0; i_f < kerF.rows(); i_f++) {
          MyMatrix<T> nz = ZeroMatrix<T>(k, k);
          for (int t = 0; t < cc; t++) {
            if (kerF(i_f, t) != 0) {
              nz += UniversalScalarConversion<T, T_residue>(kerF(i_f, t)) *
                    Zs[t];
            }
          }
          // the coordinates of nz y in O, 0 modulo p
          std::vector<MyVector<T>> cols_n;
          for (int i = 0; i < d; i++) {
            cols_n.push_back(coordsO(MyMatrix<T>(nz * O[i])));
          }
          for (int t = 0; t < d; t++) {
            ModVector row(d);
            for (int i = 0; i < d; i++) {
              row[i] = ResModOf(GetNumerator_z(cols_n[i](t)), p_u);
            }
            AppendModRow(Gm, row);
          }
        }
      }
      ModMatrix ker = FreeKernelBasisModP(ComputeRowSpaceModP(Gm, p_u), d, p_u);
      size_t n_ker = static_cast<size_t>(ker.rows());
#ifdef DEBUG_QUATERNIONIC_MASS
    os << "QUAT_MASS: time " << "kernel computed" << " " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count() << "\n";
#endif
      size_t n_el = 1;
      for (size_t i = 0; i < n_ker; i++) {
        n_el *= static_cast<size_t>(p);
        if (n_el > bounds.max_overorder) {
#ifdef DEBUG_QUATERNIONIC_MASS
          os << "QUAT_MASS: the overorder search at p=" << p
             << " is beyond the bounds\n";
#endif
          return {};
        }
      }
      // x = y / p, y = sum a_i w_i: two cheap necessary conditions in
      // integers, Tr trd(y^2) = 0 mod p^2 (x^2 is integral) and y^2 = 0
      // mod p O (p^2 x^2 is in p O'), from Tr trd(w_i w_j) and the
      // structure constants of O modulo p
      std::vector<std::vector<int64_t>> Tij(d, std::vector<int64_t>(d));
      std::vector<std::vector<std::vector<int64_t>>> Sij(
          d, std::vector<std::vector<int64_t>>(d, std::vector<int64_t>(d)));
      std::vector<std::vector<std::vector<int64_t>>> C2(
          d, std::vector<std::vector<int64_t>>(d, std::vector<int64_t>(d)));
      int64_t p2 = p * p;
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          Tij[i][j] = ResidueOfInteger(
                               QuatTrTrd(O[i], O[j], r), p2);
          MyVector<T> cij = coordsO(MyMatrix<T>(O[i] * O[j]));
          for (int t = 0; t < d; t++) {
            Sij[i][j][t] = ResidueOfInteger(cij(t), p);
            C2[i][j][t] = ResidueOfInteger(cij(t), p2);
          }
        }
      }
      // and x^2 in O' subset O^#: Tr trd(y^2 w_k) = 0 mod p^2 for every k,
      // from the Tr trd(w_i w_j w_k)
      std::vector<std::vector<std::vector<int64_t>>> T3(
          d, std::vector<std::vector<int64_t>>(d, std::vector<int64_t>(d)));
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          MyMatrix<T> Wij = O[i] * O[j];
          for (int kk = 0; kk < d; kk++) {
            T3[i][j][kk] = ResidueOfInteger(
                                    QuatTrTrd(Wij, O[kk], r), p2);
          }
        }
      }
#ifdef DEBUG_QUATERNIONIC_MASS
    os << "QUAT_MASS: time " << "filters precomputed" << " " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count() << "\n";
#endif
      std::vector<int64_t> coef(n_ker, 0), a(d);
      size_t n_q2 = 0, n_t3 = 0, n_pair = 0, n_sq = 0, n_closure = 0;
      for (size_t s = 1; s < n_el && !enlarged; s++) {
        size_t tt = s;
        for (size_t i = 0; i < n_ker; i++) {
          coef[i] = static_cast<int64_t>(tt % static_cast<size_t>(p));
          tt /= static_cast<size_t>(p);
        }
        // x and lambda x (lambda prime to p) generate the same order at p:
        // one coefficient vector per line, the last nonzero one being 1
        {
          size_t i_last = n_ker;
          for (size_t i = n_ker; i-- > 0;) {
            if (coef[i] != 0) {
              i_last = i;
              break;
            }
          }
          if (i_last == n_ker || coef[i_last] != 1) {
            continue;
          }
        }
        for (int i = 0; i < d; i++) {
          int64_t ai = 0;
          for (size_t j = 0; j < n_ker; j++) {
            ai = (ai + coef[j] * ResToSigned(ker(static_cast<int>(j), i))) % p;
          }
          a[i] = ai;
        }
        int64_t q2 = 0;
        for (int i = 0; i < d; i++) {
          if (a[i] != 0) {
            for (int j = 0; j < d; j++) {
              q2 = (q2 + a[i] * a[j] % p2 * Tij[i][j]) % p2;
            }
          }
        }
        if (q2 != 0) {
          continue;
        }
        n_q2++;
        bool t3_zero = true;
        for (int kk = 0; kk < d && t3_zero; kk++) {
          int64_t v = 0;
          for (int i = 0; i < d; i++) {
            if (a[i] != 0) {
              for (int j = 0; j < d; j++) {
                v = (v + a[i] * a[j] % p2 * T3[i][j][kk]) % p2;
              }
            }
          }
          if (v != 0) {
            t3_zero = false;
          }
        }
        if (!t3_zero) {
          continue;
        }
        n_t3++;
        // x w_k in O', so Tr trd((x w_k)(x w_l)) is an integer: with v_k the
        // coordinates of y w_k, v_k^T T v_l = 0 mod p^2 for all k, l
        {
          std::vector<std::vector<int64_t>> vk(d, std::vector<int64_t>(d, 0));
          for (int kk = 0; kk < d; kk++) {
            for (int i = 0; i < d; i++) {
              if (a[i] != 0) {
                for (int t2 = 0; t2 < d; t2++) {
                  vk[kk][t2] = (vk[kk][t2] + a[i] * C2[i][kk][t2]) % p2;
                }
              }
            }
          }
          bool pair_zero = true;
          for (int k1 = 0; k1 < d && pair_zero; k1++) {
            std::vector<int64_t> Tv(d, 0);
            for (int u = 0; u < d; u++) {
              for (int v = 0; v < d; v++) {
                Tv[u] = (Tv[u] + Tij[u][v] * vk[k1][v]) % p2;
              }
            }
            for (int k2 = 0; k2 < d && pair_zero; k2++) {
              int64_t w = 0;
              for (int u = 0; u < d; u++) {
                w = (w + vk[k2][u] * Tv[u]) % p2;
              }
              if (w != 0) {
                pair_zero = false;
              }
            }
          }
          if (!pair_zero) {
            continue;
          }
        }
        n_pair++;
        bool sq_zero = true;
        for (int t = 0; t < d && sq_zero && has_OK; t++) {
          int64_t v = 0;
          for (int i = 0; i < d; i++) {
            if (a[i] != 0) {
              for (int j = 0; j < d; j++) {
                v = (v + a[i] * a[j] % p * Sij[i][j][t]) % p;
              }
            }
          }
          if (v != 0) {
            sq_zero = false;
          }
        }
        if (!sq_zero) {
          continue;
        }
        n_sq++;
        MyMatrix<T> x = ZeroMatrix<T>(k, k);
        for (int i = 0; i < d; i++) {
          if (a[i] != 0) {
            x += ScalarFromInt64<T>(a[i]) * O[i];
          }
        }
        x /= ScalarFromInt64<T>(p);
        n_closure++;
        std::vector<MyMatrix<T>> M = O;
        M.push_back(x);
        std::optional<std::vector<MyMatrix<T>>> opt =
            ring_closure(QuatZBasis(M, k));
        if (opt) {
          O = *opt;
          enlarged = true;
        }
      }
#ifdef DEBUG_QUATERNIONIC_MASS
      os << "QUAT_MASS: overorder candidates at p=" << p << ": " << n_el
         << " (space of dimension " << n_ker << "), passing Tr y^2 "
         << n_q2 << ", Tr y^2 w " << n_t3 << ", (y w)(y w') " << n_pair
         << ", y^2 " << n_sq << ", closures " << n_closure << "\n";
#endif
      (void)n_q2;
      (void)n_t3;
      (void)n_pair;
      (void)n_sq;
      (void)n_closure;
      if (enlarged) {
        break;
      }
    }
    if (!enlarged) {
      return O;
    }
  }
}

/*
  |U(Z_p) mod p^c| and mod p^{c+1}, for the lattice of Gram matrix B with
  the action of D by the integral matrices gens (columns), and the
  integral commutant Lam of gens.
 */
template <typename T, typename Tint>
std::optional<std::pair<Tint, Tint>>
QuatLocalImageOrders(MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &gens,
                     std::vector<MyMatrix<T>> const &Lam,
                     std::vector<MyMatrix<T>> const &Forms, int64_t p, int c,
                     MultiplicityFactorBounds const &bounds,
                     [[maybe_unused]] std::ostream &os) {
  int m = B.rows();
  int e = 0;
  MyMatrix<T> Binv = Inverse(B);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      if (Binv(i, j) != T(0)) {
        e = std::max(e, RationalValuation(T(GetDenominator(Binv(i, j))), p));
      }
    }
  }
  int Wn = MultiplicityNewtonPrecision(p);
  Tint pW = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), Wn);
  T pe(1);
  pe *= MyPow<T>(UniversalScalarConversion<T, int64_t>(p), e);
  ModRowsI64 A_mod(m, std::vector<int64_t>(m)), Be(m, std::vector<int64_t>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      Tint av = ResInt(UniversalScalarConversion<Tint, T>(B(i, j)), pW);
      if (av < 0) {
        av += pW;
      }
      A_mod[i][j] = UniversalScalarConversion<int64_t, Tint>(av);
      T x = Binv(i, j) * pe;
      Tint den = GetDenominator_z(x);
      Tint inv = mod_inv(ResInt(den, pW), pW);
      Tint v = ResInt(Tint(GetNumerator_z(x) * inv), pW);
      if (v < 0) {
        v += pW;
      }
      Be[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
    }
  }
  // the residues modulo p^W in [0, p^W): the searches only use the
  // matrices modulo p^W, and their entries may not fit in 64 bits
  auto to_mod = [&](MyMatrix<T> const &Y) -> ModRowsI64 {
    ModRowsI64 Ym(m, std::vector<int64_t>(m));
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        Tint v = ResInt(UniversalScalarConversion<Tint, T>(Y(i, j)), pW);
        if (v < 0) {
          v += pW;
        }
        Ym[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
      }
    }
    return Ym;
  };
  std::vector<ModRowsI64> Om, LamM, FormsM;
  for (auto const &N : Forms) {
    FormsM.push_back(to_mod(N));
  }
  for (auto const &g : gens) {
    Om.push_back(to_mod(g));
  }
  for (auto const &Y : Lam) {
    LamM.push_back(to_mod(Y));
  }
  auto exact = [&](LocalImageOrder const &o) -> Tint {
    Tint v(1);
    for (size_t x : o.orbits) {
      v *= UniversalScalarConversion<Tint, int64_t>(static_cast<int64_t>(x));
    }
    v *= MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), o.layer_exp);
    return v;
  };
  std::vector<Tint> res;
#ifdef SANITY_CHECK_QUATERNIONIC_MASS
  int c_top = c + 1;
#else
  int c_top = c;
#endif
  std::vector<MyMatrix<Tint>> LamInt =
      UniversalStdVectorMatrixConversion<Tint, T>(Lam);
  std::optional<std::vector<ModRowsI64>> Lie = ResidueMatrices(
      IntegralLieLatticeColumns<T, Tint>(B, LamInt), ResidueModulusOf(pW));
  for (int cc = c; cc <= c_top; cc++) {
    LocalImageOrder ord;
    std::optional<std::vector<ModRowsI64>> opt = LocalOrthogonalImageGenerators(
        A_mod, Be, Om, LamM, FormsM, Lie, p, cc, e, bounds, &ord);
#ifdef DEBUG_QUATERNIONIC_MASS
    os << "QUAT_MASS: p=" << p << " c=" << cc << " e=" << e << " image "
       << (opt ? "found" : "not found") << ", order "
       << (opt ? exact(ord) : Tint(0)) << "\n";
#endif
    if (!opt) {
      return {};
    }
    res.push_back(exact(ord));
  }
  // without the sanity check the order modulo p^{c+1} is the theorem's
  return std::make_pair(res[0], res.size() > 1 ? res[1] : Tint(0));
}

/*
  The type of the component, given in the columns convention: the Gram
  matrix B, End_G (a Q-basis) and its center. r >= 1 when End_G = M_r(D)
  with a symplectic-type adjoint involution over a real center (Albert
  type III), 0 otherwise.
 */
template <typename T>
int QuaternionicRank(MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &End,
                     std::vector<MyMatrix<T>> const &Center) {
  int k = B.rows();
  MyMatrix<T> Binv = Inverse(B);
  auto adj = [&](MyMatrix<T> const &Y) -> MyMatrix<T> {
    return Binv * Y.transpose() * B;
  };
  for (auto const &Z : Center) {
    if (adj(Z) != Z) {
      return 0;
    }
  }
  int c = static_cast<int>(Center.size());
  int e = static_cast<int>(End.size());
  if (c == 0 || e % (4 * c) != 0) {
    return 0;
  }
  int r = static_cast<int>(std::llround(std::sqrt(static_cast<double>(e / (4 * c)))));
  if (4 * c * r * r != e) {
    return 0;
  }
  // the dimension of the skew-adjoint part of End_G: c r (2 r + 1) for
  // the symplectic type, c r (2 r - 1) for the orthogonal one
  MyMatrix<T> M(e, k * k);
  for (int i = 0; i < e; i++) {
    MyMatrix<T> S = End[i] + adj(End[i]);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        M(i, u * k + v) = S(u, v);
      }
    }
  }
  int skew = e - RankMat(M);
  if (skew != c * r * (2 * r + 1)) {
    return 0;
  }
  return r;
}

/*
  The mass of the quaternionic component, columns convention: B the Gram
  matrix of X, Dspan matrices spanning D acting on X over Q (the group
  elements when G spans D, X = L_rho; the e g e restricted to X = L cap e
  V_rho otherwise), FK the center K (when c > 1) and d_K. The multiplier
  ring O_X of X in D, the center and the D-linear maps all come from
  Dspan.
 */
template <typename T, typename Tint>
std::optional<T> QuaternionicComponentMass(
    MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Dspan,
    std::optional<AbelianCMField> const &FK, T const &dK, int r,
    QuaternionicMassBounds const &bounds, std::ostream &os) {
  int k = B.rows();
  auto c_start = std::chrono::steady_clock::now();
  (void)c_start;
  // the multiplier ring of X in D
  std::vector<MyMatrix<T>> O0 = QuatIntegralPoints<T, Tint>(Dspan, k);
  std::vector<MyMatrix<T>> const &gens = O0;
  // the center of D: the elements of the span of O0 commuting with O0
  std::vector<MyMatrix<T>> Center;
  {
    int d0 = static_cast<int>(O0.size());
    int kk = k * k;
    MyMatrix<T> Eq(d0, kk * d0);
    for (int i = 0; i < d0; i++) {
      for (int j = 0; j < d0; j++) {
        MyMatrix<T> C = O0[i] * O0[j] - O0[j] * O0[i];
        for (int u = 0; u < k; u++) {
          for (int v = 0; v < k; v++) {
            Eq(i, j * kk + u * k + v) = C(u, v);
          }
        }
      }
    }
    MyMatrix<T> ker = NullspaceMat(Eq);
    for (int t = 0; t < ker.rows(); t++) {
      MyMatrix<T> Z = ZeroMatrix<T>(k, k);
      for (int i = 0; i < d0; i++) {
        Z += ker(t, i) * O0[i];
      }
      Center.push_back(Z);
    }
  }
  int c = static_cast<int>(Center.size());
  int dimD = 4 * c;
  int dim_u = c * r * (2 * r + 1);
  MyMatrix<T> Binv = Inverse(B);
  auto adj = [&](MyMatrix<T> const &Y) -> MyMatrix<T> {
    return Binv * Y.transpose() * B;
  };
  if (static_cast<int>(O0.size()) != dimD || k != dimD * r) {
    std::cerr << "QUAT_MASS: the span has dimension " << O0.size()
              << " for a center of degree " << c << " and rank " << k
              << ", not that of D\n";
    throw TerminalException{1};
  }
  // the adjoint is the canonical involution: x + x^* is central
  for (auto const &w : O0) {
    MyMatrix<T> t = w + adj(w);
    for (auto const &g : gens) {
      if (t * g != g * t) {
        std::cerr << "QUAT_MASS: the adjoint is not the canonical involution\n";
        throw TerminalException{1};
      }
    }
  }
  T dK4 = dK * dK * dK * dK;
#ifdef DEBUG_QUATERNIONIC_MASS
  os << "QUAT_MASS: ctime before maximal order " << std::chrono::duration<double>(std::chrono::steady_clock::now() - c_start).count() << "\n";
#endif
  std::optional<std::vector<MyMatrix<T>>> opt_O =
      QuatMaximalOrder<T, Tint>(O0, r, dK, bounds, os);
  if (!opt_O) {
    return {};
  }
  // a Z-basis of O starting with 1, so that the first vector of each
  // copy of O in S = O^r generates a free module modulo p (with a basis
  // from the Hermite form it may be a zero divisor modulo p, which leaves
  // few columns forced in the image search)
  std::vector<MyMatrix<T>> O;
  {
    int dO = static_cast<int>(opt_O->size());
    MyMatrix<T> Ob(dO, k * k);
    for (int a = 0; a < dO; a++) {
      for (int u = 0; u < k; u++) {
        for (int v = 0; v < k; v++) {
          Ob(a, u * k + v) = (*opt_O)[a](u, v);
        }
      }
    }
    MyMatrix<T> Id = IdentityMat<T>(k);
    MyVector<T> idv(k * k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        idv(u * k + v) = Id(u, v);
      }
    }
    std::optional<MyVector<T>> c1 = SolutionMat(Ob, idv);
    if (!c1) {
      std::cerr << "QUAT_MASS: 1 is not in the span of the order\n";
      throw TerminalException{1};
    }
    MyVector<Tint> c1i(dO);
    for (int a = 0; a < dO; a++) {
      if (!IsInteger((*c1)(a))) {
        std::cerr << "QUAT_MASS: 1 is not in the order\n";
        throw TerminalException{1};
      }
      c1i(a) = UniversalScalarConversion<Tint, T>((*c1)(a));
    }
    MyMatrix<Tint> comp = ComplementToBasis(c1i);
    O.push_back(Id);
    for (int r2 = 0; r2 < comp.rows(); r2++) {
      MyMatrix<T> w = ZeroMatrix<T>(k, k);
      for (int a = 0; a < dO; a++) {
        w += UniversalScalarConversion<T, Tint>(comp(r2, a)) * (*opt_O)[a];
      }
      O.push_back(w);
    }
    if (QuatOrderDisc(O, r) != QuatOrderDisc(*opt_O, r)) {
      std::cerr << "QUAT_MASS: the new basis of O is not a basis\n";
      throw TerminalException{1};
    }
  }
#ifdef DEBUG_QUATERNIONIC_MASS
  os << "QUAT_MASS: ctime after maximal order " << std::chrono::duration<double>(std::chrono::steady_clock::now() - c_start).count() << "\n";
#endif
  // the ramified primes: |disc O| / d_K^4 = N(d(D))^2
  T qmax = QuatOrderDisc(O, r) / dK4;
  Tint qmax_i = UniversalScalarConversion<Tint, T>(qmax);
  std::vector<Tint> ram;
  if (qmax_i > 1) {
    ram = FactorsInt(qmax_i);
    std::sort(ram.begin(), ram.end());
    ram.erase(std::unique(ram.begin(), ram.end()), ram.end());
  }
  // the standard mass
  T mass_S(1);
  {
    std::vector<T> Bern = BernoulliNumbers<T>(2 * r);
    for (int i = 1; i <= r; i++) {
      T z = FK ? DedekindZetaAbelian<T>(*FK, 2 * i)
               : T(-Bern[2 * i] / ScalarFromInt64<T>(2 * i));
      mass_S *= T_abs(z);
      for (int j = 0; j < c; j++) {
        mass_S /= ScalarFromInt64<T>(2);
      }
    }
  }
  int n_ram_places = c;
  for (Tint const &p_i : ram) {
    int64_t p = UniversalScalarConversion<int64_t, Tint>(p_i);
    int v = RationalValuation(qmax, p);
    int f = 1, g = 1;
    if (FK) {
      auto [ee, ff, gg] = AbelianPrimeDecomposition(FK->n, FK->H, p);
      (void)ee;
      f = ff;
      g = gg;
    }
    if (v % (2 * f) != 0 || v / (2 * f) > g) {
      std::cerr << "QUAT_MASS: v_p(N(d(D))^2)=" << v << " at p=" << p
                << " does not fit f=" << f << " g=" << g << "\n";
      throw TerminalException{1};
    }
    int g_ram = v / (2 * f);
    n_ram_places += g_ram;
    T q = MyPow<T>(UniversalScalarConversion<T, int64_t>(p), f);
    for (int s = 0; s < g_ram; s++) {
      T qi(1);
      for (int i = 1; i <= r; i++) {
        qi *= q;
        mass_S *= qi + ScalarFromInt64<T>((i % 2 == 0) ? 1 : -1);
      }
    }
  }
  if (n_ram_places % 2 != 0) {
    std::cerr << "QUAT_MASS: an odd number of ramified places\n";
    throw TerminalException{1};
  }
  // S = O^r: the left multiplications by the generators and the Gram
  // matrices Tr trd(lambda w_a w_b^*), lambda in K
  int d = dimD;
  MyMatrix<T> Obasis(d, k * k);
  for (int a = 0; a < d; a++) {
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        Obasis(a, u * k + v) = O[a](u, v);
      }
    }
  }
  auto coords = [&](MyMatrix<T> const &Y) -> MyVector<T> {
    MyVector<T> y(k * k);
    for (int u = 0; u < k; u++) {
      for (int v = 0; v < k; v++) {
        y(u * k + v) = Y(u, v);
      }
    }
    std::optional<MyVector<T>> sol = SolutionMat(Obasis, y);
    if (!sol) {
      std::cerr << "QUAT_MASS: an element is not in the span of the order\n";
      throw TerminalException{1};
    }
    return *sol;
  };
  int kS = d * r;
  std::vector<MyMatrix<T>> gensS;
  for (auto const &g : gens) {
    MyMatrix<T> Lg = ZeroMatrix<T>(kS, kS);
    for (int a = 0; a < d; a++) {
      MyVector<T> s = coords(MyMatrix<T>(g * O[a]));
      for (int b = 0; b < d; b++) {
        if (!IsInteger(s(b))) {
          std::cerr << "QUAT_MASS: a generator is not in the maximal order\n";
          throw TerminalException{1};
        }
        for (int i = 0; i < r; i++) {
          Lg(i * d + b, i * d + a) = s(b);
        }
      }
    }
    gensS.push_back(Lg);
  }
#ifdef DEBUG_QUATERNIONIC_MASS
  os << "QUAT_MASS: ctime gensS built " << std::chrono::duration<double>(std::chrono::steady_clock::now() - c_start).count() << "\n";
#endif
  auto gram_S = [&](MyMatrix<T> const &lambda) -> MyMatrix<T> {
    MyMatrix<T> G = ZeroMatrix<T>(kS, kS);
    for (int a = 0; a < d; a++) {
      for (int b = 0; b < d; b++) {
        T v = QuatTrTrd(MyMatrix<T>(lambda * O[a]), adj(O[b]), r);
        for (int i = 0; i < r; i++) {
          G(i * d + a, i * d + b) = v;
        }
      }
    }
    return G;
  };
  // the lambda with an integral Gram matrix: a lattice, by the integral
  // points of the span of the Gram matrices of a basis of K
  std::vector<MyMatrix<T>> GramCenter;
  for (auto const &Z : Center) {
    GramCenter.push_back(gram_S(Z));
  }
  std::vector<MyMatrix<T>> GramInt = QuatIntegralPoints<T, Tint>(GramCenter, kS);
  MyMatrix<T> Gflat(c, kS * kS);
  for (int t = 0; t < c; t++) {
    for (int u = 0; u < kS; u++) {
      for (int v = 0; v < kS; v++) {
        Gflat(t, u * kS + v) = GramCenter[t](u, v);
      }
    }
  }
  std::vector<MyMatrix<T>> LambdaBasis;
  for (auto const &Gi : GramInt) {
    MyVector<T> y(kS * kS);
    for (int u = 0; u < kS; u++) {
      for (int v = 0; v < kS; v++) {
        y(u * kS + v) = Gi(u, v);
      }
    }
    std::optional<MyVector<T>> mu = SolutionMat(Gflat, y);
    if (!mu) {
      std::cerr << "QUAT_MASS: no lambda for an integral Gram matrix\n";
      throw TerminalException{1};
    }
    MyMatrix<T> lam = ZeroMatrix<T>(k, k);
    for (int t = 0; t < c; t++) {
      lam += (*mu)(t) * Center[t];
    }
    LambdaBasis.push_back(lam);
  }
  // LLL-reduced for the trace form tr(lambda lambda'), positive definite
  // on the totally real K, so that the small combinations include
  // totally positive elements
  {
    int nb = static_cast<int>(LambdaBasis.size());
    MyMatrix<T> Gl(nb, nb);
    for (int a = 0; a < nb; a++) {
      for (int b = 0; b < nb; b++) {
        Gl(a, b) = QuatTrace(MyMatrix<T>(LambdaBasis[a] * LambdaBasis[b]));
      }
    }
    LLLreduction<T, Tint> red = LLLreducedBasis<T, Tint>(Gl, os);
    std::vector<MyMatrix<T>> reduced;
    for (int a = 0; a < nb; a++) {
      MyMatrix<T> lam = ZeroMatrix<T>(k, k);
      for (int b = 0; b < nb; b++) {
        lam += UniversalScalarConversion<T, Tint>(red.Pmat(a, b)) * LambdaBasis[b];
      }
      reduced.push_back(lam);
    }
    LambdaBasis = reduced;
  }
  // the primes of the numerator and of the denominator of x
  auto prime_set = [&](T const &x) -> std::set<int64_t> {
    std::set<int64_t> S;
    for (Tint xi : {Tint(T_abs(GetNumerator_z(x))), Tint(GetDenominator_z(x))}) {
      if (xi > 1) {
        for (Tint const &q : FactorsInt(xi)) {
          S.insert(UniversalScalarConversion<int64_t, Tint>(q));
        }
      }
    }
    return S;
  };
  T detX = DeterminantMat(B);
  std::set<int64_t> primesX = prime_set(detX);
  std::optional<MyMatrix<T>> best_G;
  size_t best_extra = 0;
  T best_det(0);
  {
    // lambda = 1 (always admissible: O is an order) and the combinations
    // with coefficients in [-3, 3] of the reduced basis
    int nb = static_cast<int>(LambdaBasis.size());
    int bound = 3;
    std::vector<MyMatrix<T>> cands{IdentityMat<T>(k)};
    std::vector<int> cf(nb, -bound);
    while (true) {
      MyMatrix<T> lam = ZeroMatrix<T>(k, k);
      for (int t = 0; t < nb; t++) {
        lam += ScalarFromInt64<T>(cf[t]) * LambdaBasis[t];
      }
      cands.push_back(lam);
      int t = 0;
      while (t < nb && cf[t] == bound) {
        cf[t] = -bound;
        t++;
      }
      if (t == nb) {
        break;
      }
      cf[t]++;
    }
    // gram_S is linear in lambda: the Gram matrices of the reduced basis,
    // then their combinations
    std::vector<MyMatrix<T>> GramLB;
    for (auto const &lam : LambdaBasis) {
      GramLB.push_back(gram_S(lam));
    }
    MyMatrix<T> Gram1 = gram_S(IdentityMat<T>(k));
    std::vector<MyMatrix<T>> cand_gram{Gram1};
    {
      std::vector<int> cg(nb, -bound);
      while (true) {
        MyMatrix<T> G = ZeroMatrix<T>(kS, kS);
        for (int t = 0; t < nb; t++) {
          G += ScalarFromInt64<T>(cg[t]) * GramLB[t];
        }
        cand_gram.push_back(G);
        int t = 0;
        while (t < nb && cg[t] == bound) {
          cg[t] = -bound;
          t++;
        }
        if (t == nb) {
          break;
        }
        cg[t]++;
      }
    }
    for (size_t ic = 0; ic < cands.size(); ic++) {
      MyMatrix<T> const &lam = cands[ic];
      MyMatrix<T> G = cand_gram[ic];
      if (!IsZeroMatrix(lam) && QuatIsPositiveDefinite(G)) {
        T det = DeterminantMat(G);
        size_t extra = 0;
        for (int64_t q : prime_set(det)) {
          if (primesX.count(q) == 0) {
            extra++;
          }
        }
        if (!best_G || extra < best_extra ||
            (extra == best_extra && det < best_det)) {
          best_G = G;
          best_extra = extra;
          best_det = det;
        }
      }
    }
  }
  if (!best_G) {
    std::cerr << "QUAT_MASS: no totally positive lambda found\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_QUATERNIONIC_MASS
  os << "QUAT_MASS: ctime lambda chosen " << std::chrono::duration<double>(std::chrono::steady_clock::now() - c_start).count() << "\n";
#endif
  MyMatrix<T> const &BS = *best_G;
  // the adjoint of the left multiplication by w on S is that by w^*
  for (size_t ig = 0; ig < gens.size(); ig++) {
    MyMatrix<T> wbar = adj(gens[ig]);
    MyMatrix<T> Lbar = ZeroMatrix<T>(kS, kS);
    for (int a = 0; a < d; a++) {
      MyVector<T> sv = coords(MyMatrix<T>(wbar * O[a]));
      for (int b = 0; b < d; b++) {
        for (int i = 0; i < r; i++) {
          Lbar(i * d + b, i * d + a) = sv(b);
        }
      }
    }
    if (MyMatrix<T>(gensS[ig].transpose() * BS) != MyMatrix<T>(BS * Lbar)) {
      std::cerr << "QUAT_MASS: the adjoint on S is not the canonical "
                   "involution\n";
      throw TerminalException{1};
    }
  }
  // the Lie lattices and their discriminants
  std::vector<MyMatrix<T>> uX = QuatIntegralCommutant<T, Tint>(B, gens, true);
  std::vector<MyMatrix<T>> uS = QuatIntegralCommutant<T, Tint>(BS, gensS, true);
  if (static_cast<int>(uX.size()) != dim_u ||
      static_cast<int>(uS.size()) != dim_u) {
    std::cerr << "QUAT_MASS: the Lie lattices have ranks " << uX.size()
              << " and " << uS.size() << ", not " << dim_u << "\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_QUATERNIONIC_MASS
  os << "QUAT_MASS: ctime Lie lattices " << std::chrono::duration<double>(std::chrono::steady_clock::now() - c_start).count() << "\n";
#endif
  T discX = QuatTraceDiscriminant(uX);
  T discS = QuatTraceDiscriminant(uS);
  std::vector<MyMatrix<T>> LamX = QuatIntegralCommutant<T, Tint>(B, gens, false);
  std::vector<MyMatrix<T>> LamS = QuatIntegralCommutant<T, Tint>(BS, gensS, false);
  // the image searches commute with a Z-basis of the multiplier ring of X
  // and of O on S: modulo p the generators of G may span less (Q8 spans
  // the Lipschitz order, of index 2 in the Hurwitz order)
  std::vector<MyMatrix<T>> OmS;
  for (auto const &w : O) {
    MyMatrix<T> Lw = ZeroMatrix<T>(kS, kS);
    for (int a = 0; a < d; a++) {
      MyVector<T> sv = coords(MyMatrix<T>(w * O[a]));
      for (int b = 0; b < d; b++) {
        for (int i = 0; i < r; i++) {
          Lw(i * d + b, i * d + a) = sv(b);
        }
      }
    }
    OmS.push_back(Lw);
  }
  // the integral forms B(w x, y) = x^T C_w^T B y, w in D
  auto forms = [&](MyMatrix<T> const &Bz, std::vector<MyMatrix<T>> const &W)
      -> std::vector<MyMatrix<T>> {
    std::vector<MyMatrix<T>> L;
    for (auto const &w : W) {
      L.push_back(MyMatrix<T>(w.transpose() * Bz));
    }
    return QuatIntegralPoints<T, Tint>(L, Bz.rows());
  };
  std::vector<MyMatrix<T>> FormsX = forms(B, O0);
  std::vector<MyMatrix<T>> FormsS = forms(BS, OmS);
  // the primes where X and S may differ
  std::set<int64_t> P = primesX;
  for (int64_t q : prime_set(DeterminantMat(BS))) {
    P.insert(q);
  }
  for (int64_t q : prime_set(QuatOrderDisc(O0, r) / QuatOrderDisc(O, r))) {
    P.insert(q);
  }
  for (int64_t q : prime_set(discX / discS)) {
    if (P.count(q) == 0) {
      std::cerr << "QUAT_MASS: the Lie lattices differ at p=" << q
                << " outside the bad primes\n";
      throw TerminalException{1};
    }
  }
#ifdef DEBUG_QUATERNIONIC_MASS
  os << "QUAT_MASS: c=" << c << " r=" << r << " |disc O0|="
     << QuatOrderDisc(O0, r) << " |disc O|=" << QuatOrderDisc(O, r)
     << " mass(S)=" << mass_S << " det X=" << detX << " det S=" << DeterminantMat(BS)
     << " disc u(X)=" << discX << " disc u(S)=" << discS << "\n";
#endif
  T mass = mass_S;
  // [O : O_X], whose primes are those where O_X is not maximal
  T index_O = QuatOrderDisc(O0, r) / QuatOrderDisc(O, r);
  for (int64_t p : P) {
    /*
      r = 1 and O_X maximal at p: X_p = O_p x (left ideals of a local
      maximal order are principal), h(a x, b x) = a h(x, x) b^* with h(x,
      x) central, and the right multiplication by x^{-1} is a D-linear
      similitude of X_p onto S_p, so U(X_p) and U(S_p) are conjugate and
      the ratio is 1, without the image search (whose first column would
      run over the orbit, of size about p^{dim u}). The Lie lattices then
      have the same discriminant at p, which is checked.
     */
    if (r == 1 && RationalValuation(index_O, p) == 0) {
      if (RationalValuation(discX, p) != RationalValuation(discS, p)) {
        std::cerr << "QUAT_MASS: rank one over the maximal order at p=" << p
                  << " but the Lie lattices differ there\n";
        throw TerminalException{1};
      }
#ifdef DEBUG_QUATERNIONIC_MASS
      os << "QUAT_MASS: p=" << p << " rank one over the maximal order, ratio 1\n";
#endif
      continue;
    }
    int cc = (p == 2) ? 2 : 1;
    std::optional<std::pair<Tint, Tint>> nX =
        QuatLocalImageOrders<T, Tint>(B, O0, LamX, FormsX, p, cc, bounds.image,
                                      os);
    if (!nX) {
      return {};
    }
    std::optional<std::pair<Tint, Tint>> nS =
        QuatLocalImageOrders<T, Tint>(BS, OmS, LamS, FormsS, p, cc,
                                      bounds.image, os);
    if (!nS) {
      return {};
    }
    Tint pdim = MyPow<Tint>(UniversalScalarConversion<Tint, int64_t>(p), dim_u);
#ifdef SANITY_CHECK_QUATERNIONIC_MASS
    if (nX->second != nX->first * pdim || nS->second != nS->first * pdim) {
      std::cerr << "QUAT_MASS: the images at p=" << p << " do not grow by p^"
                << dim_u << " from p^" << cc << " to p^" << cc + 1 << "\n";
      throw TerminalException{1};
    }
#endif
    int vX = RationalValuation(discX, p);
    int vS = RationalValuation(discS, p);
    if ((vX - vS) % 2 != 0) {
      std::cerr << "QUAT_MASS: the Lie discriminants differ by an odd power of "
                << p << "\n";
      throw TerminalException{1};
    }
    T ratio = UniversalScalarConversion<T, Tint>(nS->first) /
              UniversalScalarConversion<T, Tint>(nX->first);
    int ex = (vX - vS) / 2;
    for (int i = 0; i < (ex >= 0 ? ex : -ex); i++) {
      if (ex >= 0) {
        ratio *= UniversalScalarConversion<T, int64_t>(p);
      } else {
        ratio /= UniversalScalarConversion<T, int64_t>(p);
      }
    }
#ifdef DEBUG_QUATERNIONIC_MASS
    os << "QUAT_MASS: p=" << p << " |U(X_p) mod p^c|=" << nX->first
       << " |U(S_p) mod p^c|=" << nS->first << " v_p(disc u)=" << vX << ","
       << vS << " ratio=" << ratio << "\n";
#endif
    mass *= ratio;
  }
  return mass;
}

// clang-format off
#endif  // SRC_GENERA_QUATERNIONICMASS_H_
// clang-format on
