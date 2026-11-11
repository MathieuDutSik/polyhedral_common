// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_QUATERNIONICMASS_H_
#define SRC_GENERA_QUATERNIONICMASS_H_

// clang-format off
#include "MultiplicityLocalFactor.h"
#include "AbelianQuadraticMass.h"
#include "SiegelSpecialValues.h"
#include "ClassicLLL.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "factorizations.h"
#include <algorithm>
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
  type, and the centraliser C_rho = U is an inner form of Sp_{2r} over K,
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
  size_t max_overorder = size_t(1) << 20;
  MultiplicityFactorBounds image;
};

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
  MyMatrix<Tint> ker = NullspaceIntTrMat(MyMatrix<Tint>(Eq.transpose()));
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
QuatMaximalOrder(std::vector<MyMatrix<T>> O, int r, T const &dK4,
                 QuaternionicMassBounds const &bounds,
                 [[maybe_unused]] std::ostream &os) {
  int k = O[0].rows();
  int d = static_cast<int>(O.size());
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
    for (Tint const &p_i : primes) {
      int64_t p = UniversalScalarConversion<int64_t, Tint>(p_i);
      // the x = sum a_i w_i / p in O^#: sum_i a_i Tr trd(w_i w_j) = 0 mod p
      ModMat Gm(d, std::vector<int64_t>(d));
      for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
          T v = QuatTrTrd(O[i], O[j], r);
          Gm[j][i] = ModP(UniversalScalarConversion<int64_t, T>(v), p);
        }
      }
      ModRowSpace rs = ModRowReduce(Gm, p);
      std::vector<bool> is_piv(d, false);
      for (size_t pv : rs.pivots) {
        is_piv[pv] = true;
      }
      ModMat ker;
      for (int f = 0; f < d; f++) {
        if (is_piv[f]) {
          continue;
        }
        std::vector<int64_t> v(d, 0);
        v[f] = 1;
        for (size_t i = 0; i < rs.pivots.size(); i++) {
          v[rs.pivots[i]] = ModP(-rs.basis[i][f], p);
        }
        ker.push_back(v);
      }
      size_t n_el = 1;
      for (size_t i = 0; i < ker.size(); i++) {
        n_el *= static_cast<size_t>(p);
        if (n_el > bounds.max_overorder) {
#ifdef DEBUG_QUATERNIONIC_MASS
          os << "QUAT_MASS: the overorder search at p=" << p
             << " is beyond the bounds\n";
#endif
          return {};
        }
      }
      std::vector<int64_t> coef(ker.size(), 0);
      for (size_t s = 1; s < n_el && !enlarged; s++) {
        size_t tt = s;
        for (size_t i = 0; i < ker.size(); i++) {
          coef[i] = static_cast<int64_t>(tt % static_cast<size_t>(p));
          tt /= static_cast<size_t>(p);
        }
        MyMatrix<T> x = ZeroMatrix<T>(k, k);
        for (int i = 0; i < d; i++) {
          int64_t a = 0;
          for (size_t j = 0; j < ker.size(); j++) {
            a = (a + coef[j] * ker[j][i]) % p;
          }
          if (a != 0) {
            x += ScalarFromInt64<T>(a) * O[i];
          }
        }
        x /= ScalarFromInt64<T>(p);
        // Tr_{K/Q} trd(x^2) is an integer for x integral
        if (!IsInteger(QuatTrTrd(x, x, r))) {
          continue;
        }
        std::vector<MyMatrix<T>> M = O;
        M.push_back(x);
        std::optional<std::vector<MyMatrix<T>>> opt =
            ring_closure(QuatZBasis(M, k));
        if (opt) {
          O = *opt;
          enlarged = true;
        }
      }
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
        e = std::max(e, ValuationOfRational(T(GetDenominator(Binv(i, j))), p));
      }
    }
  }
  int Wn = MultiplicityNewtonPrecision(p);
  Tint pW(1);
  for (int i = 0; i < Wn; i++) {
    pW *= UniversalScalarConversion<Tint, int64_t>(p);
  }
  T pe(1);
  for (int i = 0; i < e; i++) {
    pe *= UniversalScalarConversion<T, int64_t>(p);
  }
  ModMat A_mod(m, std::vector<int64_t>(m)), Be(m, std::vector<int64_t>(m));
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
  auto to_mod = [&](MyMatrix<T> const &Y) -> ModMat {
    ModMat Ym(m, std::vector<int64_t>(m));
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
  std::vector<ModMat> Om, LamM, FormsM;
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
    for (size_t i = 0; i < o.layer_exp; i++) {
      v *= UniversalScalarConversion<Tint, int64_t>(p);
    }
    return v;
  };
  std::vector<Tint> res;
#ifdef SANITY_CHECK_QUATERNIONIC_MASS
  int c_top = c + 1;
#else
  int c_top = c;
#endif
  for (int cc = c; cc <= c_top; cc++) {
    LocalImageOrder ord;
    std::optional<std::vector<ModMat>> opt = LocalOrthogonalImageGenerators(
        A_mod, Be, Om, LamM, FormsM, p, cc, e, bounds, &ord);
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
  matrix B, End_G (a Q-basis) and its centre. r >= 1 when End_G = M_r(D)
  with a symplectic-type adjoint involution over a real centre (Albert
  type III), 0 otherwise.
 */
template <typename T>
int QuaternionicRank(MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &End,
                     std::vector<MyMatrix<T>> const &Centre) {
  int k = B.rows();
  MyMatrix<T> Binv = Inverse(B);
  auto adj = [&](MyMatrix<T> const &Y) -> MyMatrix<T> {
    return Binv * Y.transpose() * B;
  };
  for (auto const &Z : Centre) {
    if (adj(Z) != Z) {
      return 0;
    }
  }
  int c = static_cast<int>(Centre.size());
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
  V_rho otherwise), FK the centre K (when c > 1) and d_K. The multiplier
  ring O_X of X in D, the centre and the D-linear maps all come from
  Dspan.
 */
template <typename T, typename Tint>
std::optional<T> QuaternionicComponentMass(
    MyMatrix<T> const &B, std::vector<MyMatrix<T>> const &Dspan,
    std::optional<AbelianCMField> const &FK, T const &dK, int r,
    QuaternionicMassBounds const &bounds, std::ostream &os) {
  int k = B.rows();
  // the multiplier ring of X in D
  std::vector<MyMatrix<T>> O0 = QuatIntegralPoints<T, Tint>(Dspan, k);
  std::vector<MyMatrix<T>> const &gens = O0;
  // the centre of D: the elements of the span of O0 commuting with O0
  std::vector<MyMatrix<T>> Centre;
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
      Centre.push_back(Z);
    }
  }
  int c = static_cast<int>(Centre.size());
  int dimD = 4 * c;
  int dim_u = c * r * (2 * r + 1);
  MyMatrix<T> Binv = Inverse(B);
  auto adj = [&](MyMatrix<T> const &Y) -> MyMatrix<T> {
    return Binv * Y.transpose() * B;
  };
  if (static_cast<int>(O0.size()) != dimD || k != dimD * r) {
    std::cerr << "QUAT_MASS: the span has dimension " << O0.size()
              << " for a centre of degree " << c << " and rank " << k
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
  std::optional<std::vector<MyMatrix<T>>> opt_O =
      QuatMaximalOrder<T, Tint>(O0, r, dK4, bounds, os);
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
    int v = ValuationOfRational(qmax, p);
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
    T q(1);
    for (int i = 0; i < f; i++) {
      q *= UniversalScalarConversion<T, int64_t>(p);
    }
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
  std::vector<MyMatrix<T>> GramCentre;
  for (auto const &Z : Centre) {
    GramCentre.push_back(gram_S(Z));
  }
  std::vector<MyMatrix<T>> GramInt = QuatIntegralPoints<T, Tint>(GramCentre, kS);
  MyMatrix<T> Gflat(c, kS * kS);
  for (int t = 0; t < c; t++) {
    for (int u = 0; u < kS; u++) {
      for (int v = 0; v < kS; v++) {
        Gflat(t, u * kS + v) = GramCentre[t](u, v);
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
      lam += (*mu)(t) * Centre[t];
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
    for (auto const &lam : cands) {
      MyMatrix<T> G = gram_S(lam);
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
  for (int64_t p : P) {
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
    Tint pdim(1);
    for (int i = 0; i < dim_u; i++) {
      pdim *= UniversalScalarConversion<Tint, int64_t>(p);
    }
#ifdef SANITY_CHECK_QUATERNIONIC_MASS
    if (nX->second != nX->first * pdim || nS->second != nS->first * pdim) {
      std::cerr << "QUAT_MASS: the images at p=" << p << " do not grow by p^"
                << dim_u << " from p^" << cc << " to p^" << cc + 1 << "\n";
      throw TerminalException{1};
    }
#endif
    int vX = ValuationOfRational(discX, p);
    int vS = ValuationOfRational(discS, p);
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
