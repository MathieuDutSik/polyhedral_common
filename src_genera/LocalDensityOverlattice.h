// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_LOCALDENSITYOVERLATTICE_H_
#define SRC_GENERA_LOCALDENSITYOVERLATTICE_H_

// clang-format off
#include "LocalDensityBadPrimes.h"
#include "MAT_MatrixInt.h"
#include <cstdint>
#include <map>
#include <optional>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE
#endif

/*
  alpha_p(S, L) for a lattice L that is not unimodular at p, from Siegel
  series only.

  Split L_p = U + P with U the unimodular Jordan component (even unimodular
  at p = 2) and P the rest, and let R = P(-1), the form P negated, whose
  discriminant form is minus that of L. Gluing L + R along the graph of that
  isometry gives a lattice N unimodular at p (even unimodular at p = 2) of
  rank n + rank P, with det N = (-1)^{rank P} det U up to squares. Every
  primitive embedding of R in N has a complement of rank n, of discriminant
  form that of L and of determinant that of L, hence isometric to L; and the
  vectors of N orthogonal to a primitive copy of R are, modulo p^k, exactly
  the complement modulo p^k. Counting the representations of S + R by N
  whose R-part is primitive therefore gives
      alpha_p(S, L) = alpha_p^{R-prim}(S + R, N) / alpha_p^*(R, N),
  the normalisations matching since (n + rank R) - rank R = n.

  Both are primitive densities. Sorting the representations of S + L0 by
  the saturation Lambda of their L0-part (an integral, even at p = 2,
  overlattice of L0 of p-power index):
      alpha_p(S + L0, N) = sum_{Lambda >= L0} p^{(m+1-n) v[Lambda:L0]}
                           alpha_p^{Lambda-prim}(S + Lambda, N),
  m = rank S (the Jacobian of X_R = Y G is |det G|^{n_N} on the vectors and
  |det G|^{rank R + 1 + m} on the Gram matrices). The full densities are
  Siegel series (AlphaUnimodularTarget), and the relation is inverted over
  the finite poset of the overlattices of R.
 */

// The unimodular / non-unimodular split of a Gram matrix at p.
template <typename T> struct UnimodularSplit {
  // the rank of the unimodular component, and its determinant class
  int r_U;
  int64_t detU_unit;
  // the Gram matrix of the rest, integral, in the normalisation of the input
  MyMatrix<T> P;
};

template <typename T, typename Tint>
UnimodularSplit<T> ComputeUnimodularSplit(MyMatrix<T> const &Tgram,
                                          int64_t p) {
  UnimodularSplit<T> res;
  res.r_U = 0;
  int64_t M = UnitModulus(p);
  res.detU_unit = 1;
  std::vector<MyMatrix<T>> blocksP;
  if (p == 2) {
    SiegelJordanBlocks q = JordanBlocks2<T, Tint>(Tgram / T(2));
    for (auto const &b : q.blocks) {
      if (b.kind == SiegelBlockKind::Unit) {
        // 2^e u in T / 2, so 2^{e+1} u in T: never unimodular for even T
        MyMatrix<T> G(1, 1);
        G(0, 0) = PowerInt(T(2), b.e + 1) * ScalarFromInt64<T>(b.unit);
        blocksP.push_back(G);
        continue;
      }
      MyMatrix<T> G(2, 2);
      T scale = PowerInt(T(2), b.e);
      if (b.kind == SiegelBlockKind::H) {
        G(0, 0) = 0;
        G(1, 1) = 0;
      } else {
        G(0, 0) = 2 * scale;
        G(1, 1) = 2 * scale;
      }
      G(0, 1) = scale;
      G(1, 0) = scale;
      if (b.e == 0) {
        res.r_U += 2;
        res.detU_unit = MulMod(res.detU_unit,
                               (b.kind == SiegelBlockKind::H) ? 7 : 3, M);
      } else {
        blocksP.push_back(G);
      }
    }
  } else {
    std::vector<PadicClass> D = OddDiagonalClasses<T, Tint>(Tgram, p);
    for (auto const &c : D) {
      // the diagonal entry 2 c of T
      if (c.val == 0) {
        res.r_U++;
        res.detU_unit = MulMod(res.detU_unit, MulMod(2, c.unit, M), M);
      } else {
        MyMatrix<T> G(1, 1);
        G(0, 0) = T(2) * PowerInt(ScalarFromInt64<T>(p), c.val) *
                  ScalarFromInt64<T>(c.unit);
        blocksP.push_back(G);
      }
    }
  }
  int mP = 0;
  for (auto const &G : blocksP) {
    mP += G.rows();
  }
  res.P = ZeroMatrix<T>(mP, mP);
  int pos = 0;
  for (auto const &G : blocksP) {
    for (int i = 0; i < G.rows(); i++) {
      for (int j = 0; j < G.rows(); j++) {
        res.P(pos + i, pos + j) = G(i, j);
      }
    }
    pos += G.rows();
  }
  return res;
}

// The integral (even at p = 2) overlattices of p-power index of the lattice
// with Gram matrix R, as bases (rows) in the coordinates of R.
template <typename T> struct OverlatticeSet {
  std::vector<MyMatrix<T>> basis;
  std::vector<MyMatrix<T>> gram;
  std::vector<int> val;
};

template <typename T> bool LatticeContained(MyMatrix<T> const &B1,
                                            MyMatrix<T> const &B2) {
  return IsIntegralMatrix(MyMatrix<T>(B1 * Inverse(B2)));
}

template <typename T, typename Tint>
OverlatticeSet<T> EnumerateOverlattices(MyMatrix<T> const &R, int64_t p) {
  int m = R.rows();
  bool even = (p == 2);
  T p_T = ScalarFromInt64<T>(p);
  OverlatticeSet<T> os;
  os.basis.push_back(IdentityMat<T>(m));
  os.gram.push_back(R);
  os.val.push_back(0);
  // the non-zero vectors of (Z/p)^m with first non-zero coordinate 1
  std::vector<std::vector<int64_t>> dirs;
  {
    std::vector<int64_t> c(m, 0);
    int64_t total = 1;
    for (int i = 0; i < m; i++) {
      total *= p;
    }
    for (int64_t idx = 1; idx < total; idx++) {
      int64_t w = idx;
      for (int i = 0; i < m; i++) {
        c[i] = w % p;
        w /= p;
      }
      int first = 0;
      while (c[first] == 0) {
        first++;
      }
      if (c[first] == 1) {
        dirs.push_back(c);
      }
    }
  }
  for (size_t pos = 0; pos < os.basis.size(); pos++) {
    MyMatrix<T> B = os.basis[pos];
    for (auto const &c : dirs) {
      MyVector<T> v = ZeroVector<T>(m);
      for (int i = 0; i < m; i++) {
        if (c[i] != 0) {
          for (int j = 0; j < m; j++) {
            v(j) += ScalarFromInt64<T>(c[i]) * B(i, j);
          }
        }
      }
      v /= p_T;
      MyVector<T> Rv = R * v;
      bool ok = true;
      for (int i = 0; i < m && ok; i++) {
        T prod(0);
        for (int j = 0; j < m; j++) {
          prod += B(i, j) * Rv(j);
        }
        if (!IsInteger(prod)) {
          ok = false;
        }
      }
      if (!ok) {
        continue;
      }
      T nv = v.dot(Rv);
      if (!IsInteger(nv) || (even && !IsInteger(nv / T(2)))) {
        continue;
      }
      MyMatrix<T> Ext(m + 1, m);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          Ext(i, j) = B(i, j);
        }
      }
      for (int j = 0; j < m; j++) {
        Ext(m, j) = v(j);
      }
      MyMatrix<T> Bnew = GetZbasis(Ext);
      int val_new = os.val[pos] + 1;
      bool known = false;
      for (size_t i = 0; i < os.basis.size(); i++) {
        if (os.val[i] == val_new && LatticeContained(Bnew, os.basis[i])) {
          known = true;
          break;
        }
      }
      if (!known) {
        os.basis.push_back(Bnew);
        os.gram.push_back(Bnew * R * Bnew.transpose());
        os.val.push_back(val_new);
      }
    }
  }
  return os;
}

template <typename T>
MyMatrix<T> BlockDiagonal2(MyMatrix<T> const &A, MyMatrix<T> const &B) {
  int a = A.rows(), b = B.rows();
  MyMatrix<T> M = ZeroMatrix<T>(a + b, a + b);
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < a; j++) {
      M(i, j) = A(i, j);
    }
  }
  for (int i = 0; i < b; i++) {
    for (int j = 0; j < b; j++) {
      M(a + i, a + j) = B(i, j);
    }
  }
  return M;
}

/*
  The density of the definition at a unimodular target. At p = 2,
  AlphaUnimodularTarget is 2^{m-1} times it, the convention under which the
  global constant of LocalDensity.h is c = 2; the gluing below needs the
  densities themselves, since its two sides have different ranks.
 */
template <typename T, typename Tint>
T TrueAlphaUnimodularTarget(MyMatrix<T> const &S, int64_t p, int k, int xi_T,
                            SiegelSeriesCache<T> &cache) {
  T val = AlphaUnimodularTarget<T, Tint>(S, p, k, xi_T, cache);
  if (p == 2) {
    val /= PowerInt(T(2), S.rows() - 1);
  }
  return val;
}

/*
  The primitive densities alpha^{Lambda-prim}(S + Lambda, N) for all the
  overlattices Lambda of R, by inversion of
      f(L0) = sum_{Lambda >= L0} p^{c v[Lambda:L0]} g(Lambda),  c = m+1-n.
  Returns g(R).
 */
template <typename T, typename Tint>
T PrimitiveDensityOverlattices(MyMatrix<T> const &S, OverlatticeSet<T> const &os,
                               int64_t p, int k_N, int xi_N, int n,
                               SiegelSeriesCache<T> &cache) {
  int m = S.rows();
  int c = m + 1 - n;
  T p_T = ScalarFromInt64<T>(p);
  size_t K = os.basis.size();
  std::vector<T> f(K), g(K);
  for (size_t i = 0; i < K; i++) {
    MyMatrix<T> G = (m == 0) ? os.gram[i] : BlockDiagonal2(S, os.gram[i]);
    f[i] = TrueAlphaUnimodularTarget<T, Tint>(G, p, k_N, xi_N, cache);
  }
  // from the largest overlattices down
  std::vector<size_t> order(K);
  for (size_t i = 0; i < K; i++) {
    order[i] = i;
  }
  std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b) { return os.val[a] > os.val[b]; });
  for (size_t idx : order) {
    T val = f[idx];
    for (size_t j = 0; j < K; j++) {
      if (os.val[j] > os.val[idx] &&
          LatticeContained(os.basis[idx], os.basis[j])) {
        val -= PowerInt(p_T, c * (os.val[j] - os.val[idx])) * g[j];
      }
    }
    g[idx] = val;
  }
  return g[0];
}

/*
  alpha_p(S, L) by the gluing into a unimodular lattice, in the convention of
  AlphaUnimodularTarget (at p = 2, 2^{m-1} times the density of the
  definition). Nothing when the rank of N would be odd (only possible at odd
  p; the Siegel series bridge is pinned for even rank).
 */
template <typename T, typename Tint>
std::optional<T> OverlatticeDensity(MyMatrix<T> const &S, MyMatrix<T> const &Tgram,
                                    int64_t p, SiegelSeriesCache<T> &cache) {
  int n = Tgram.rows();
  UnimodularSplit<T> us = ComputeUnimodularSplit<T, Tint>(Tgram, p);
  int mP = us.P.rows();
  int n_N = n + mP;
  if (n_N % 2 != 0) {
    return {};
  }
  // det N = (-1)^{mP} det U up to squares
  int64_t detN = SignedPower(mP, us.detU_unit);
  int xi_N = KroneckerSymbol(SignedPower(n_N / 2, detN), p);
  if (mP == 0) {
    return AlphaUnimodularTarget<T, Tint>(S, p, n / 2, xi_N, cache);
  }
  // the result is returned in the convention of AlphaUnimodularTarget
  T conv = (p == 2) ? PowerInt(T(2), S.rows() - 1) : T(1);
  MyMatrix<T> R = -us.P;
  OverlatticeSet<T> os = EnumerateOverlattices<T, Tint>(R, p);
  MyMatrix<T> empty(0, 0);
  T num = PrimitiveDensityOverlattices<T, Tint>(S, os, p, n_N / 2, xi_N, n,
                                                cache);
  T den = PrimitiveDensityOverlattices<T, Tint>(empty, os, p, n_N / 2, xi_N, n,
                                                cache);
  if (den == 0) {
    std::cerr << "LDOV: R has no primitive embedding in N\n";
    throw TerminalException{1};
  }
  return conv * num / den;
}


/*
  What the gluing needs of a genus at p, computed once: the overlattices of
  R, the class of N, and the denominator alpha^*(R, N).
 */
template <typename T> struct OverlatticeGenusData {
  int64_t p;
  int n;
  int n_N;
  int xi_N;
  bool unimodular;
  OverlatticeSet<T> os;
  T den;
};

template <typename T, typename Tint>
std::optional<OverlatticeGenusData<T>>
ComputeOverlatticeGenusData(MyMatrix<T> const &Tgram, int64_t p,
                            SiegelSeriesCache<T> &cache) {
  OverlatticeGenusData<T> data;
  data.p = p;
  data.n = Tgram.rows();
  UnimodularSplit<T> us = ComputeUnimodularSplit<T, Tint>(Tgram, p);
  int mP = us.P.rows();
  data.n_N = data.n + mP;
  if (data.n_N % 2 != 0) {
    return {};
  }
  int64_t detN = SignedPower(mP, us.detU_unit);
  data.xi_N = KroneckerSymbol(SignedPower(data.n_N / 2, detN), p);
  data.unimodular = (mP == 0);
  if (!data.unimodular) {
    data.os = EnumerateOverlattices<T, Tint>(-us.P, p);
    MyMatrix<T> empty(0, 0);
    data.den = PrimitiveDensityOverlattices<T, Tint>(
        empty, data.os, p, data.n_N / 2, data.xi_N, data.n, cache);
    if (data.den == 0) {
      std::cerr << "LDOV: R has no primitive embedding in N\n";
      throw TerminalException{1};
    }
  }
  return data;
}

// alpha_p(S, L) from the genus data, in the convention of
// AlphaUnimodularTarget.
template <typename T, typename Tint>
T OverlatticeDensityFromData(MyMatrix<T> const &S,
                             OverlatticeGenusData<T> const &data,
                             SiegelSeriesCache<T> &cache) {
  int64_t p = data.p;
  if (data.unimodular) {
    return AlphaUnimodularTarget<T, Tint>(S, p, data.n / 2, data.xi_N, cache);
  }
  T conv = (p == 2) ? PowerInt(T(2), S.rows() - 1) : T(1);
  T num = PrimitiveDensityOverlattices<T, Tint>(S, data.os, p, data.n_N / 2,
                                                data.xi_N, data.n, cache);
  return conv * num / data.den;
}

/*
  alpha_p(S, T) by exhaustive counting modulo p^k, for tiny cases only (the
  loop is over p^{k n m} matrices): the density of the definition, as an
  oracle (at p = 2 the other functions differ from it by 2^{m-1}).
  Returns the value at level k and fails unless it equals the value at k-1.
 */
template <typename T>
std::optional<T> BruteForceDensity(MyMatrix<T> const &S, MyMatrix<T> const &Tgram,
                                   int64_t p, int k) {
  int n = Tgram.rows();
  int m = S.rows();
  auto count_at = [&](int kk) -> T {
    int64_t q = PowerInt64(p, kk);
    int64_t q2 = 2 * q;
    int64_t total = 1;
    for (int i = 0; i < n * m; i++) {
      total *= q;
    }
    std::vector<int64_t> Ti(n * n), Si(m * m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        Ti[i * n + j] = UniversalScalarConversion<int64_t, T>(Tgram(i, j));
      }
    }
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        Si[i * m + j] = UniversalScalarConversion<int64_t, T>(S(i, j));
      }
    }
    std::vector<int64_t> X(n * m);
    int64_t cnt = 0;
    for (int64_t idx = 0; idx < total; idx++) {
      int64_t w = idx;
      for (int i = 0; i < n * m; i++) {
        X[i] = w % q;
        w /= q;
      }
      bool ok = true;
      for (int a = 0; a < m && ok; a++) {
        for (int b = a; b < m && ok; b++) {
          int64_t val = 0;
          for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
              val += X[i * m + a] * Ti[i * n + j] * X[j * m + b];
            }
          }
          int64_t mod = (a == b) ? q2 : q;
          if (ModPositive(val - Si[a * m + b], mod) != 0) {
            ok = false;
          }
        }
      }
      if (ok) {
        cnt++;
      }
    }
    T p_T = ScalarFromInt64<T>(p);
    return ScalarFromInt64<T>(cnt) / PowerInt(p_T, kk * (n * m - m * (m + 1) / 2));
  };
  T v1 = count_at(k - 1);
  T v2 = count_at(k);
  if (v1 != v2) {
    return {};
  }
  return v2;
}


/*
  M_S / mass(G) for the genus G of the even lattice Tgram of even rank n and
  a positive definite even S of rank m with n - m >= 3, from the local
  densities alone:
      2 * (generic part) * prod_{p | 2 det S det T} alpha_p(S, T),
  alpha_p from the Siegel series where T is unimodular at p, from the unit
  chain and the reduction formula at the odd primes dividing det T, and
  from the gluing into an even unimodular lattice at p = 2 when det T is
  even.
 */
template <typename T> struct GenusDensityDataFull {
  GenusDensityData odd;
  std::optional<OverlatticeGenusData<T>> two_adic;
};

template <typename T, typename Tint>
GenusDensityDataFull<T>
ComputeGenusDensityDataFull(MyMatrix<T> const &Tgram,
                            SiegelSeriesCache<T> &cache) {
  GenusDensityDataFull<T> data;
  data.odd = ComputeGenusDensityData<T, Tint>(Tgram);
  if (data.odd.detT % 2 == 0) {
    data.two_adic = ComputeOverlatticeGenusData<T, Tint>(Tgram, 2, cache);
    if (!data.two_adic) {
      std::cerr << "LDOV: the gluing at p = 2 gave a lattice of odd rank\n";
      throw TerminalException{1};
    }
  }
  return data;
}

template <typename T, typename Tint>
T MassRatioFromDensitiesFull(MyMatrix<T> const &S,
                             GenusDensityDataFull<T> const &data,
                             LocalDensityCache<T> &cache) {
  int n = data.odd.n;
  int m = S.rows();
  if (n % 2 != 0) {
    std::cerr << "LDOV: the mass ratio needs an even rank\n";
    throw TerminalException{1};
  }
  if (n - m < 3) {
    std::cerr << "LDOV: the mass ratio needs n - m >= 3\n";
    throw TerminalException{1};
  }
  int64_t detS = DeterminantInt64<T, Tint>(S);
  int64_t detT = data.odd.detT;
  int k = n / 2;
  T val = T(2) * GenericDensityPart<T>(n, m, detS, detT);
  for (int64_t p : BadPrimes(detS, detT)) {
    if (detT % p != 0) {
      int xi_T = KroneckerSymbol(SignedPower(k, detT), p);
      val *= AlphaUnimodularTarget<T, Tint>(S, p, k, xi_T, cache.siegel);
    } else if (p == 2) {
      val *= OverlatticeDensityFromData<T, Tint>(S, *data.two_adic, cache.siegel);
    } else {
      val *= DensityOddDiagonal<T>(OddDiagonalClasses<T, Tint>(S, p),
                                   data.odd.odd_diagonal.at(p), p, cache.bad);
    }
  }
  return val;
}

// clang-format off
#endif  // SRC_GENERA_LOCALDENSITYOVERLATTICE_H_
// clang-format on
