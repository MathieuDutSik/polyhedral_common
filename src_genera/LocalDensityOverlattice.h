// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_LOCALDENSITYOVERLATTICE_H_
#define SRC_GENERA_LOCALDENSITYOVERLATTICE_H_

// clang-format off
#include "LocalDensityBadPrimes.h"
#include "MAT_MatrixInt.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <unordered_set>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_LOCAL_DENSITY_OVERLATTICE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE
#endif

//
// Local isometry classes
//

// The Jordan invariants at an odd prime: for every level, its rank and the
// Legendre symbol of the product of its units (of T / 2; the factor 2 is a
// unit and changes the classes of both sides alike).
inline std::vector<std::array<int, 3>>
OddLocalSymbolFromClasses(std::vector<PadicClass> const &classes, int64_t p) {
  std::map<int, std::pair<int, int64_t>> levels;
  for (auto const &c : classes) {
    auto &lv = levels[c.val];
    if (lv.first == 0) {
      lv.second = 1;
    }
    lv.first++;
    lv.second = MulMod(lv.second, c.unit, p);
  }
  std::vector<std::array<int, 3>> res;
  for (auto const &kv : levels) {
    res.push_back({kv.first, kv.second.first, LegendreUnit(kv.second.second, p)});
  }
  return res;
}

template <typename T, typename Tint>
std::vector<std::array<int, 3>> OddLocalSymbol(MyMatrix<T> const &Gram,
                                               int64_t p) {
  return OddLocalSymbolFromClasses(OddDiagonalClasses<T, Tint>(Gram, p), p);
}

/*
  The canonical 2-adic symbol of an even lattice (Conway-Sloane, Sphere
  Packings, Lattices and Groups, ch. 15, Sec. 7.5-7.6; the reduction follows
  Sage's canonical_2_adic_reduction): one quintuple (scale exponent, rank,
  determinant, type, oddity) per Jordan component, the determinants reduced
  to +-1, the oddities fused over the compartments and the signs walked to
  the front of every train. Two 2-adic lattices are isometric iff their
  canonical symbols are equal.
 */
inline std::vector<std::array<int, 5>>
CanonicalTwoAdicReduction(std::vector<std::array<int, 5>> sym) {
  for (auto &c : sym) {
    c[2] = (c[2] == 1 || c[2] == 7) ? 1 : -1;
  }
  // compartments: maximal runs of odd components of consecutive scales
  std::vector<std::vector<int>> compartments;
  {
    int i = 0, r = static_cast<int>(sym.size());
    while (i < r) {
      if (sym[i][3] == 1) {
        int v = sym[i][0];
        std::vector<int> c;
        while (i < r && sym[i][3] == 1 && sym[i][0] == v) {
          c.push_back(i);
          i++;
          v++;
        }
        compartments.push_back(c);
      } else {
        i++;
      }
    }
  }
  for (auto const &compart : compartments) {
    int oddity = 0;
    for (int i : compart) {
      oddity += sym[i][4];
      sym[i][4] = 0;
    }
    sym[compart[0]][4] = ModPositive(oddity, 8);
  }
  // trains
  std::vector<std::vector<int>> trains;
  {
    std::vector<std::array<int, 5>> ext = sym;
    ext.push_back({sym.back()[0] + 1, 0, 1, 0, 0});
    std::vector<int> new_train{0};
    for (size_t i = 1; i + 1 < ext.size(); i++) {
      auto const &prev = ext[i - 1];
      auto const &cur = ext[i];
      if (cur[0] - prev[0] > 2) {
        trains.push_back(new_train);
        new_train = {static_cast<int>(i)};
      } else if (cur[0] - prev[0] == 2 && cur[3] * prev[3] == 0) {
        trains.push_back(new_train);
        new_train = {static_cast<int>(i)};
      } else if (prev[3] == 0 && cur[3] == 0) {
        trains.push_back(new_train);
        new_train = {static_cast<int>(i)};
      } else {
        new_train.push_back(static_cast<int>(i));
      }
    }
    trains.push_back(new_train);
  }
  for (auto const &train : trains) {
    int t = static_cast<int>(train.size());
    for (int i = 0; i < t - 1; i++) {
      int t1 = train[t - i - 1];
      if (sym[t1][2] == -1) {
        sym[t1][2] = 1;
        sym[t1 - 1][2] *= -1;
        for (auto const &compart : compartments) {
          bool in = std::find(compart.begin(), compart.end(), t1 - 1) != compart.end() ||
                    std::find(compart.begin(), compart.end(), t1) != compart.end();
          if (in) {
            sym[compart[0]][4] = ModPositive(sym[compart[0]][4] + 4, 8);
          }
        }
      }
    }
  }
  return sym;
}

// The 2-adic symbol of an even Gram matrix, from one Jordan splitting:
// a unit u at level e of T / 2 is an odd component 2^{e+1} u of T; the
// binary blocks h (det -1 = 7) and y (det 3) of level e are components of
// scale 2^e. Neither adds to the oddity: <u> + y = <3u, -u, -u> and
// <u> + h = <u, 1, -1> have the oddity u of <u>. The determinants are kept
// mod 8 (Sage's symbol_tuple_list); they depend on the decomposition, the
// canonical symbol below does not.
inline std::vector<std::array<int, 5>>
TwoAdicSymbolFromSplitting(std::vector<SiegelJordanBlock> const &blocks) {
  std::map<int, std::array<int, 5>> comp;
  for (auto const &b : blocks) {
    int s = (b.kind == SiegelBlockKind::Unit) ? b.e + 1 : b.e;
    auto iter = comp.find(s);
    if (iter == comp.end()) {
      iter = comp.emplace(s, std::array<int, 5>{s, 0, 1, 0, 0}).first;
    }
    auto &c = iter->second;
    if (b.kind == SiegelBlockKind::Unit) {
      c[1] += 1;
      c[2] = static_cast<int>(MulMod(c[2], b.unit, 8));
      c[3] = 1;
      c[4] = ModPositive(c[4] + b.unit, 8);
    } else if (b.kind == SiegelBlockKind::H) {
      c[1] += 2;
      c[2] = static_cast<int>(MulMod(c[2], 7, 8));
    } else {
      c[1] += 2;
      c[2] = static_cast<int>(MulMod(c[2], 3, 8));
    }
  }
  std::vector<std::array<int, 5>> sym;
  for (auto &kv : comp) {
    if (kv.second[3] == 0) {
      kv.second[4] = 0;
    }
    sym.push_back(kv.second);
  }
  return sym;
}

template <typename T, typename Tint>
std::vector<std::array<int, 5>> TwoAdicSymbolRaw(MyMatrix<T> const &Gram) {
  return TwoAdicSymbolFromSplitting(JordanSplitting2<T, Tint>(Gram / T(2)));
}

template <typename T, typename Tint>
std::vector<std::array<int, 5>> TwoAdicSymbol(MyMatrix<T> const &Gram) {
  return CanonicalTwoAdicReduction(TwoAdicSymbolRaw<T, Tint>(Gram));
}

// The local symbol at p flattened to integers: equal iff locally isometric
template <typename T, typename Tint>
std::vector<int> LocalSymbolKey(MyMatrix<T> const &Gram, int64_t p) {
  std::vector<int> key;
  if (p == 2) {
    for (auto const &c : TwoAdicSymbol<T, Tint>(Gram)) {
      key.insert(key.end(), c.begin(), c.end());
    }
  } else {
    for (auto const &c : OddLocalSymbol<T, Tint>(Gram, p)) {
      key.insert(key.end(), c.begin(), c.end());
    }
  }
  return key;
}

/*
  The same key for an integral (even at p = 2) Gram matrix in int64, by a
  Jordan splitting modulo p^K, K = v_p(det) + 5: pivoting on an entry of
  minimal valuation e, every update subtracts a product of a quotient known
  modulo p^{K-e} by an entry of valuation >= e, so the entries stay exact
  modulo p^K, which determines the units modulo 8 (or p) of every Jordan
  block. Returns nothing when p^{2K} would overflow; the caller then uses
  LocalSymbolKey.
 */
inline int64_t ValuationResidue(int64_t x, int64_t p, int K) {
  if (x == 0) {
    return K;
  }
  int v = 0;
  while (x % p == 0) {
    x /= p;
    v++;
  }
  return v;
}

inline int64_t InverseModPower(int64_t u, int64_t mod) {
  // u odd (p = 2) or prime to p: extended Euclid
  int64_t a = ModPositive(u, mod), b = mod, x0 = 1, x1 = 0;
  while (b != 0) {
    int64_t q = a / b;
    int64_t t = a - q * b;
    a = b;
    b = t;
    t = x0 - q * x1;
    x0 = x1;
    x1 = t;
  }
  return ModPositive(x0, mod);
}

inline std::optional<std::vector<int>>
LocalSymbolKeyModular(MyMatrix<int64_t> const &Gram, int64_t p, int v_det) {
  int n = Gram.rows();
  int K = v_det + 5;
  int64_t M = 1;
  for (int i = 0; i < K; i++) {
    if (M > (int64_t(1) << 30) / p) {
      return {};
    }
    M *= p;
  }
  std::vector<std::vector<int64_t>> A(n, std::vector<int64_t>(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      A[i][j] = ModPositive(Gram(i, j), M);
    }
  }
  std::vector<int> active(n);
  for (int i = 0; i < n; i++) {
    active[i] = i;
  }
  auto remove = [&](int i) -> void {
    active.erase(std::find(active.begin(), active.end(), i));
  };
  int64_t pe_div;
  std::vector<SiegelJordanBlock> blocks2;
  std::vector<PadicClass> classes;
  int64_t inv2 = (p == 2) ? 0 : InverseModPower(2, p);
  while (!active.empty()) {
    int best_v = K + 1, bi = -1, bj = -1;
    for (int i : active) {
      for (int j : active) {
        if (j < i) {
          continue;
        }
        int v = static_cast<int>(ValuationResidue(A[i][j], p, K));
        // a diagonal pivot when it has the minimal valuation
        if (v < best_v || (v == best_v && i == j && bi != bj)) {
          best_v = v;
          bi = i;
          bj = j;
        }
      }
    }
    if (best_v >= K - 3) {
      // a degenerate matrix, or precision exhausted
      return {};
    }
    int e = best_v;
    pe_div = 1;
    for (int i = 0; i < e; i++) {
      pe_div *= p;
    }
    if (bi != bj && p != 2) {
      // x_i <- x_i + x_j makes the diagonal of valuation e
      int i = bi, j = bj;
      for (int k : active) {
        A[i][k] = ModPositive(A[i][k] + A[j][k], M);
      }
      for (int k : active) {
        A[k][i] = (k == i) ? A[k][i] : A[i][k];
      }
      A[i][i] = ModPositive(A[i][i] + A[j][i], M);
      bj = bi;
      if (ValuationResidue(A[i][i], p, K) != e) {
        std::cerr << "LDOV: the diagonal pivot has the wrong valuation\n";
        throw TerminalException{1};
      }
    }
    if (bi == bj) {
      int i = bi;
      int64_t u = A[i][i] / pe_div;
      int64_t uinv = InverseModPower(u, M / pe_div);
      remove(i);
      for (int j : active) {
        int64_t fj = ModPositive((A[j][i] / pe_div) * uinv, M / pe_div);
        for (int k : active) {
          A[j][k] = ModPositive(A[j][k] - fj * A[i][k], M);
        }
      }
      if (p == 2) {
        blocks2.push_back({e - 1, SiegelBlockKind::Unit, ModPositive(u, 8)});
      } else {
        classes.push_back({e, MulMod(ModPositive(u, p), inv2, p)});
      }
    } else {
      int i = bi, j = bj;
      int64_t a = A[i][i], b = A[i][j], c = A[j][j];
      int64_t pe2 = pe_div * pe_div;
      int64_t d = a * c - b * b;
      int64_t d_red = ModPositive(d / pe2, M / pe_div);
      int64_t dinv = InverseModPower(d_red, M / pe_div);
      remove(i);
      remove(j);
      std::vector<std::pair<int64_t, int64_t>> coef;
      for (int k : active) {
        int64_t num_a = c * A[i][k] - b * A[j][k];
        int64_t num_b = a * A[j][k] - b * A[i][k];
        int64_t al = ModPositive((num_a / pe2) % (M / pe_div) * dinv, M / pe_div);
        int64_t be = ModPositive((num_b / pe2) % (M / pe_div) * dinv, M / pe_div);
        coef.push_back({al, be});
      }
      for (size_t x = 0; x < active.size(); x++) {
        int kx = active[x];
        for (int l : active) {
          A[kx][l] = ModPositive(A[kx][l] - coef[x].first * A[i][l] -
                                     coef[x].second * A[j][l],
                                 M);
        }
      }
      bool is_y = (ModPositive(d_red, 8) == 3);
      blocks2.push_back({e, is_y ? SiegelBlockKind::Y : SiegelBlockKind::H, 0});
    }
  }
  std::vector<int> key;
  if (p == 2) {
    for (auto const &c : CanonicalTwoAdicReduction(TwoAdicSymbolFromSplitting(blocks2))) {
      key.insert(key.end(), c.begin(), c.end());
    }
  } else {
    for (auto const &c : OddLocalSymbolFromClasses(classes, p)) {
      key.insert(key.end(), c.begin(), c.end());
    }
  }
  return key;
}

template <typename T, typename Tint>
bool LocallyIsometric(MyMatrix<T> const &G1, MyMatrix<T> const &G2, int64_t p) {
  return LocalSymbolKey<T, Tint>(G1, p) == LocalSymbolKey<T, Tint>(G2, p);
}

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


// A basis of the left kernel {c : c G = 0 mod p} of an integral G, by
// Gaussian elimination over F_p.
template <typename T, typename Tint>
std::vector<std::vector<int64_t>> KernelBasisModP(MyMatrix<T> const &G,
                                                  int64_t p) {
  int m = G.rows();
  Tint p_T(static_cast<long>(p));
  // the columns of G mod p, as the rows of the system c G = 0
  std::vector<std::vector<int64_t>> A(m, std::vector<int64_t>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      Tint g = UniversalScalarConversion<Tint, T>(G(j, i));
      A[i][j] = UniversalScalarConversion<int64_t, Tint>(ResInt(g, p_T));
    }
  }
  // row reduction of A (m equations in the m unknowns c)
  std::vector<int> pivot_col;
  int row = 0;
  for (int col = 0; col < m && row < m; col++) {
    int piv = -1;
    for (int r = row; r < m; r++) {
      if (A[r][col] % p != 0) {
        piv = r;
        break;
      }
    }
    if (piv == -1) {
      continue;
    }
    std::swap(A[row], A[piv]);
    int64_t inv = InverseModInt64(A[row][col], p);
    for (int c = 0; c < m; c++) {
      A[row][c] = MulMod(A[row][c], inv, p);
    }
    for (int r = 0; r < m; r++) {
      if (r != row && A[r][col] != 0) {
        int64_t f = A[r][col];
        for (int c = 0; c < m; c++) {
          A[r][c] = ModPositive(A[r][c] - MulMod(f, A[row][c], p), p);
        }
      }
    }
    pivot_col.push_back(col);
    row++;
  }
  std::vector<int> free_col;
  for (int col = 0; col < m; col++) {
    if (std::find(pivot_col.begin(), pivot_col.end(), col) == pivot_col.end()) {
      free_col.push_back(col);
    }
  }
  int t = static_cast<int>(free_col.size());
  // a basis of the kernel, one vector per free unknown
  std::vector<std::vector<int64_t>> basis;
  for (int f = 0; f < t; f++) {
    std::vector<int64_t> v(m, 0);
    v[free_col[f]] = 1;
    for (size_t r = 0; r < pivot_col.size(); r++) {
      v[pivot_col[r]] = ModPositive(-A[r][free_col[f]], p);
    }
    basis.push_back(v);
  }
  return basis;
}

// The vectors c of (Z/p)^m with c G = 0 mod p, one per line (first non-zero
// coordinate 1).
template <typename T, typename Tint>
std::vector<std::vector<int64_t>> ProjectiveKernelModP(MyMatrix<T> const &G,
                                                       int64_t p) {
  int m = G.rows();
  std::vector<std::vector<int64_t>> basis = KernelBasisModP<T, Tint>(G, p);
  int t = static_cast<int>(basis.size());
  std::vector<std::vector<int64_t>> dirs;
  int64_t total = 1;
  for (int i = 0; i < t; i++) {
    total *= p;
  }
  for (int64_t idx = 1; idx < total; idx++) {
    std::vector<int64_t> coef(t);
    int64_t w = idx;
    for (int i = 0; i < t; i++) {
      coef[i] = w % p;
      w /= p;
    }
    int first = 0;
    while (coef[first] == 0) {
      first++;
    }
    if (coef[first] != 1) {
      continue;
    }
    std::vector<int64_t> c(m, 0);
    for (int i = 0; i < t; i++) {
      for (int j = 0; j < m; j++) {
        c[j] = ModPositive(c[j] + coef[i] * basis[i][j], p);
      }
    }
    dirs.push_back(c);
  }
  return dirs;
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
  std::vector<std::unordered_set<MyMatrix<Tint>>> seen;
  for (size_t pos = 0; pos < os.basis.size(); pos++) {
    MyMatrix<T> B = os.basis[pos];
    // v = c B / p lies in the dual of Lambda iff c G_Lambda = 0 mod p, so the
    // extensions of index p come from the left kernel of G_Lambda mod p
    std::vector<std::vector<int64_t>> dirs =
        ProjectiveKernelModP<T, Tint>(os.gram[pos], p);
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
      // an overlattice of index p^val lies in p^{-val} Z^m: the Hermite
      // normal form of p^val Bnew identifies it
      MyMatrix<T> Bscaled = Bnew;
      for (int i = 0; i < val_new; i++) {
        Bscaled *= p_T;
      }
      MyMatrix<Tint> key = ComputeRowHermiteNormalForm(
                               UniversalMatrixConversion<Tint, T>(Bscaled))
                               .second;
      if (static_cast<int>(seen.size()) <= val_new) {
        seen.resize(val_new + 1);
      }
      if (seen[val_new].insert(key).second) {
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
  The same g(R) from the elementary overlattices alone. The interval
  [R, Lambda] of the poset of overlattices is the lattice of subgroups of
  Lambda / R, whose Moebius function vanishes unless Lambda / R = (Z/p)^k,
  where it is (-1)^k p^{k(k-1)/2} (P. Hall), so
      g(R) = sum_k (-1)^k p^{k(k-1)/2 + c k}
                   sum_{Lambda / R = (Z/p)^k} f(Lambda),
  and f(Lambda), a Siegel series of S + Lambda, depends only on the local
  class of Lambda: the elementary overlattices are tabulated once, by
  (k, class), with a count and one representative.

  They are R + W, W a subspace of V = D(R)[p] (the vectors c / p, c in the
  left kernel of the Gram matrix mod p): R + W is integral iff
  c G d^t = 0 mod p^2 on W, and even iff moreover c G c^t = 0 mod 2 p^2
  (automatic at odd p). Both depend on c modulo p only. The subspaces are
  enumerated in reduced row echelon form, each from its canonical parent,
  the span of its rows but the last.
 */
template <typename T> struct ElementaryOverlatticeClass {
  int k;
  int64_t count;
  MyMatrix<T> gram;
};

template <typename T, typename Tint>
std::vector<ElementaryOverlatticeClass<T>>
ElementaryOverlatticeClasses(MyMatrix<T> const &R, int64_t p) {
  int m = R.rows();
  std::vector<std::vector<int64_t>> kb = KernelBasisModP<T, Tint>(R, p);
  int t = static_cast<int>(kb.size());
  int64_t p2 = p * p;
  int64_t mod_q = (p == 2) ? 2 * p2 : p2;
  std::vector<std::vector<int64_t>> Ri(m, std::vector<int64_t>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      Ri[i][j] = UniversalScalarConversion<int64_t, T>(R(i, j));
    }
  }
  auto form = [&](std::vector<int64_t> const &x,
                  std::vector<int64_t> const &y) -> int64_t {
    int64_t s = 0;
    for (int i = 0; i < m; i++) {
      if (x[i] != 0) {
        for (int j = 0; j < m; j++) {
          s += x[i] * Ri[i][j] * y[j];
        }
      }
    }
    return s;
  };
  // q on the basis mod 2p^2 (p^2 at odd p), b mod p^2
  std::vector<int64_t> Qd(t);
  std::vector<std::vector<int64_t>> Bm(t, std::vector<int64_t>(t));
  for (int i = 0; i < t; i++) {
    for (int j = 0; j < t; j++) {
      Bm[i][j] = ModPositive(form(kb[i], kb[j]), p2);
    }
    Qd[i] = ModPositive(form(kb[i], kb[i]), mod_q);
  }
  auto q_of = [&](std::vector<int64_t> const &a) -> int64_t {
    int64_t s = 0;
    for (int i = 0; i < t; i++) {
      if (a[i] == 0) {
        continue;
      }
      s += a[i] * a[i] * Qd[i];
      for (int j = i + 1; j < t; j++) {
        s += 2 * a[i] * a[j] * Bm[i][j];
      }
    }
    return ModPositive(s, mod_q);
  };
  auto b_of = [&](std::vector<int64_t> const &a,
                  std::vector<int64_t> const &a2) -> int64_t {
    int64_t s = 0;
    for (int i = 0; i < t; i++) {
      if (a[i] != 0) {
        for (int j = 0; j < t; j++) {
          s += a[i] * a2[j] * Bm[i][j];
        }
      }
    }
    return ModPositive(s, p2);
  };
  // the free column of each kernel vector: kb[i] has a 1 there and 0 at
  // the free columns of the others
  std::vector<int> free_col(t);
  for (int i = 0; i < t; i++) {
    free_col[i] = -1;
    for (int c = 0; c < m && free_col[i] < 0; c++) {
      bool ok = (kb[i][c] == 1);
      for (int i2 = 0; i2 < t && ok; i2++) {
        if (i2 != i && kb[i2][c] != 0) {
          ok = false;
        }
      }
      if (ok) {
        free_col[i] = c;
      }
    }
    if (free_col[i] < 0) {
      std::cerr << "LDOV: kernel basis without a free column\n";
      throw TerminalException{1};
    }
  }
  std::map<std::pair<int, std::vector<int>>, size_t> index;
  std::vector<ElementaryOverlatticeClass<T>> classes;
  int v_det_R = 0;
  {
    int64_t detR = DeterminantInt64<T, Tint>(R);
    while (detR % p == 0) {
      detR /= p;
      v_det_R++;
    }
  }
  // A basis of R + W: the rows of W lifted to c_r / p replace the basis
  // vectors e_i at the free columns of their pivots (a row of W in reduced
  // echelon form has a 1 at its pivot and 0 at the other pivots, so the
  // lift has a 1 at that free column and 0 at the others). The Gram matrix
  // is integral and computed in int64.
  auto record = [&](std::vector<std::vector<int64_t>> const &rows,
                    std::vector<int> const &pivots) -> void {
    int k = static_cast<int>(rows.size());
    std::vector<std::vector<int64_t>> lifts(k, std::vector<int64_t>(m, 0));
    std::vector<bool> replaced(m, false);
    for (int r = 0; r < k; r++) {
      for (int i = 0; i < t; i++) {
        if (rows[r][i] != 0) {
          for (int j = 0; j < m; j++) {
            lifts[r][j] += rows[r][i] * kb[i][j];
          }
        }
      }
      replaced[free_col[pivots[r]]] = true;
    }
    // each basis vector: (is it a lift, index)
    std::vector<std::pair<bool, int>> vecs;
    for (int i = 0; i < m; i++) {
      if (!replaced[i]) {
        vecs.push_back({false, i});
      }
    }
    for (int r = 0; r < k; r++) {
      vecs.push_back({true, r});
    }
    std::vector<int64_t> Rx(m);
    MyMatrix<int64_t> gram(m, m);
    for (int a = 0; a < m; a++) {
      // R times the vector a (times p for a lift)
      for (int j = 0; j < m; j++) {
        Rx[j] = vecs[a].first ? 0 : Ri[vecs[a].second][j];
      }
      if (vecs[a].first) {
        for (int j = 0; j < m; j++) {
          for (int i = 0; i < m; i++) {
            Rx[j] += lifts[vecs[a].second][i] * Ri[i][j];
          }
        }
      }
      for (int b = 0; b < m; b++) {
        int64_t val = vecs[b].first ? 0 : Rx[vecs[b].second];
        if (vecs[b].first) {
          for (int j = 0; j < m; j++) {
            val += lifts[vecs[b].second][j] * Rx[j];
          }
        }
        int64_t div = (vecs[a].first ? p : 1) * (vecs[b].first ? p : 1);
        gram(a, b) = val / div;
      }
    }
    std::optional<std::vector<int>> opt =
        LocalSymbolKeyModular(gram, p, v_det_R - 2 * k);
    std::vector<int> sym =
        opt ? *opt
            : LocalSymbolKey<T, Tint>(UniversalMatrixConversion<T, int64_t>(gram),
                                      p);
#ifdef SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE_EXTENSIVE
    if (opt && *opt != LocalSymbolKey<T, Tint>(
                           UniversalMatrixConversion<T, int64_t>(gram), p)) {
      std::cerr << "LDOV: the modular local symbol disagrees\n";
      throw TerminalException{1};
    }
#endif
    std::pair<int, std::vector<int>> key{k, sym};
    auto iter = index.find(key);
    if (iter == index.end()) {
      index[key] = classes.size();
      classes.push_back({k, 1, UniversalMatrixConversion<T, int64_t>(gram)});
    } else {
      classes[iter->second].count++;
    }
  };
  std::vector<std::vector<int64_t>> rows;
  std::vector<int> pivots;
#ifdef TIMINGS_LOCAL_DENSITY_OVERLATTICE
  MicrosecondTime time;
  int64_t n_subspace = 0, n_candidate = 0;
  int64_t time_record = 0;
#endif
  std::function<void()> rec = [&]() -> void {
#ifdef TIMINGS_LOCAL_DENSITY_OVERLATTICE
    n_subspace++;
    MicrosecondTime time_rec;
#endif
    record(rows, pivots);
#ifdef TIMINGS_LOCAL_DENSITY_OVERLATTICE
    time_record += time_rec.eval_int64();
#endif
    int start = pivots.empty() ? 0 : pivots.back() + 1;
    for (int j = start; j < t; j++) {
      bool blocked = false;
      for (auto const &r : rows) {
        if (r[j] != 0) {
          blocked = true;
        }
      }
      if (blocked) {
        continue;
      }
      // the free coordinates of the new row: after j, off the pivots
      std::vector<int> free_pos;
      for (int f = j + 1; f < t; f++) {
        if (std::find(pivots.begin(), pivots.end(), f) == pivots.end()) {
          free_pos.push_back(f);
        }
      }
      int64_t n_assign = 1;
      for (size_t i = 0; i < free_pos.size(); i++) {
        n_assign *= p;
      }
      std::vector<int64_t> a(t, 0);
      a[j] = 1;
      for (int64_t idx = 0; idx < n_assign; idx++) {
        int64_t w = idx;
        for (int f : free_pos) {
          a[f] = w % p;
          w /= p;
        }
#ifdef TIMINGS_LOCAL_DENSITY_OVERLATTICE
        n_candidate++;
#endif
        if (q_of(a) != 0) {
          continue;
        }
        bool ortho = true;
        for (auto const &r : rows) {
          if (b_of(a, r) != 0) {
            ortho = false;
            break;
          }
        }
        if (!ortho) {
          continue;
        }
        rows.push_back(a);
        pivots.push_back(j);
        rec();
        rows.pop_back();
        pivots.pop_back();
      }
    }
  };
  rec();
#ifdef TIMINGS_LOCAL_DENSITY_OVERLATTICE
  std::cerr << "LDOV: ElementaryOverlatticeClasses t=" << t << " subspaces="
            << n_subspace << " candidates=" << n_candidate
            << " classes=" << classes.size() << " record=" << time_record
            << "us total=" << time << "\n";
#endif
  return classes;
}

template <typename T, typename Tint>
T PrimitiveDensityElementary(
    MyMatrix<T> const &S,
    std::vector<ElementaryOverlatticeClass<T>> const &classes, int64_t p,
    int k_N, int xi_N, int n, SiegelSeriesCache<T> &cache) {
  int m = S.rows();
  int c = m + 1 - n;
  T p_T = ScalarFromInt64<T>(p);
  T g(0);
  for (auto const &cls : classes) {
    MyMatrix<T> G = (m == 0) ? cls.gram : BlockDiagonal2(S, cls.gram);
    T f = TrueAlphaUnimodularTarget<T, Tint>(G, p, k_N, xi_N, cache);
    T coef = PowerInt(p_T, cls.k * (cls.k - 1) / 2 + c * cls.k) *
             ScalarFromInt64<T>(cls.count);
    if (cls.k % 2 == 1) {
      coef = -coef;
    }
    g += coef * f;
  }
  return g;
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
  std::vector<ElementaryOverlatticeClass<T>> classes =
      ElementaryOverlatticeClasses<T, Tint>(R, p);
  MyMatrix<T> empty(0, 0);
  T num = PrimitiveDensityElementary<T, Tint>(S, classes, p, n_N / 2, xi_N, n,
                                              cache);
  T den = PrimitiveDensityElementary<T, Tint>(empty, classes, p, n_N / 2,
                                              xi_N, n, cache);
#ifdef SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE_EXTENSIVE
  std::cerr << "LDOV: SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE_EXTENSIVE, the "
               "Moebius inversion over all the overlattices\n";
  OverlatticeSet<T> os = EnumerateOverlattices<T, Tint>(R, p);
  T num_all = PrimitiveDensityOverlattices<T, Tint>(S, os, p, n_N / 2, xi_N, n,
                                                    cache);
  T den_all = PrimitiveDensityOverlattices<T, Tint>(empty, os, p, n_N / 2,
                                                    xi_N, n, cache);
  if (num_all != num || den_all != den) {
    std::cerr << "LDOV: elementary overlattices give " << num << " / " << den
              << ", all overlattices " << num_all << " / " << den_all << "\n";
    throw TerminalException{1};
  }
  std::cerr << "LDOV: SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE_EXTENSIVE, done\n";
#endif
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
  std::vector<ElementaryOverlatticeClass<T>> classes;
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
    data.classes = ElementaryOverlatticeClasses<T, Tint>(-us.P, p);
    MyMatrix<T> empty(0, 0);
    data.den = PrimitiveDensityElementary<T, Tint>(
        empty, data.classes, p, data.n_N / 2, data.xi_N, data.n, cache);
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
  T num = PrimitiveDensityElementary<T, Tint>(S, data.classes, p,
                                              data.n_N / 2, data.xi_N, data.n,
                                              cache);
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


//
// M_R for R of full rank
//

// For every S and p, the number of even overlattices of S of each index p^v
// and local symbol at p: the enumeration does not depend on the genus.
template <typename T> struct FullRankOverlatticeCache {
  using SymbolCount = std::map<std::pair<int, std::vector<int>>, int64_t>;
  std::map<int64_t, std::unordered_map<MyMatrix<T>, SymbolCount>> counts;
  // the counts at one index p^h, by (p, h), then S
  std::map<std::pair<int64_t, int>,
           std::unordered_map<MyMatrix<T>, std::map<std::vector<int>, int64_t>>>
      counts_at_index;
};

template <typename T, typename Tint>
typename FullRankOverlatticeCache<T>::SymbolCount const &
OverlatticeSymbolCount(MyMatrix<T> const &S, int64_t p,
                       FullRankOverlatticeCache<T> &cache) {
  auto &by_S = cache.counts[p];
  auto iter = by_S.find(S);
  if (iter != by_S.end()) {
    return iter->second;
  }
  typename FullRankOverlatticeCache<T>::SymbolCount sc;
  OverlatticeSet<T> os = EnumerateOverlattices<T, Tint>(S, p);
  for (size_t i = 0; i < os.basis.size(); i++) {
    sc[{os.val[i], LocalSymbolKey<T, Tint>(os.gram[i], p)}]++;
  }
  return by_S.emplace(S, std::move(sc)).first->second;
}

// The canonical row Hermite normal form (positive pivots, the entries above
// a pivot reduced to [0, pivot)) of an integral matrix of full column rank
// m, as its m x m upper triangle; small entries, int64.
inline std::vector<int64_t>
HermiteNormalFormInt64(std::vector<std::vector<int64_t>> M, int m) {
  int nr = static_cast<int>(M.size());
  int row = 0;
  for (int c = 0; c < m && row < nr; c++) {
    int piv = -1;
    for (int i = row; i < nr; i++) {
      if (M[i][c] != 0) {
        piv = i;
        break;
      }
    }
    if (piv < 0) {
      continue;
    }
    std::swap(M[row], M[piv]);
    for (int i = row + 1; i < nr; i++) {
      while (M[i][c] != 0) {
        int64_t q = M[row][c] / M[i][c];
        for (int j = c; j < m; j++) {
          M[row][j] -= q * M[i][j];
        }
        std::swap(M[row], M[i]);
      }
    }
    if (M[row][c] < 0) {
      for (int j = c; j < m; j++) {
        M[row][j] = -M[row][j];
      }
    }
    for (int i = 0; i < row; i++) {
      int64_t q = M[i][c] / M[row][c];
      if (M[i][c] - q * M[row][c] < 0) {
        q--;
      }
      if (q != 0) {
        for (int j = c; j < m; j++) {
          M[i][j] -= q * M[row][j];
        }
      }
    }
    row++;
  }
  if (row != m) {
    std::cerr << "LDOV: HermiteNormalFormInt64 expects a full rank matrix\n";
    throw TerminalException{1};
  }
  std::vector<int64_t> key;
  for (int i = 0; i < m; i++) {
    for (int j = i; j < m; j++) {
      key.push_back(M[i][j]);
    }
  }
  return key;
}

/*
  The even (integral at odd p) overlattices of index p^h of the lattice
  with Gram matrix R, counted by local symbol. Each is carried as the
  Hermite normal form of P B, P = p^h and B a basis in the coordinates of
  R (integral, since p^h kills Lambda / R), level by level, each level
  deduplicated; the extensions of Lambda come from the left kernel of its
  Gram matrix mod p, as in EnumerateOverlattices, but only up to the index
  needed, in int64, and the class from the modular Jordan splitting.
 */
template <typename T, typename Tint>
std::map<std::vector<int>, int64_t>
OverlatticeClassCountsAtIndex(MyMatrix<T> const &R, int64_t p, int h) {
  int m = R.rows();
  std::vector<std::vector<int64_t>> Ri(m, std::vector<int64_t>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      Ri[i][j] = UniversalScalarConversion<int64_t, T>(R(i, j));
    }
  }
  int64_t P = 1;
  for (int i = 0; i < h; i++) {
    P *= p;
  }
  int64_t P2 = P * P;
  int64_t mod_even = (p == 2) ? 2 : 1;
  auto unpack = [&](std::vector<int64_t> const &key) {
    std::vector<std::vector<int64_t>> M(m, std::vector<int64_t>(m, 0));
    size_t pos = 0;
    for (int i = 0; i < m; i++) {
      for (int j = i; j < m; j++) {
        M[i][j] = key[pos++];
      }
    }
    return M;
  };
  // the products x R y^t of the rows of M, divided by P^2
  auto gram_of = [&](std::vector<std::vector<int64_t>> const &M) {
    MyMatrix<int64_t> G(m, m);
    std::vector<std::vector<int64_t>> MR(m, std::vector<int64_t>(m, 0));
    for (int i = 0; i < m; i++) {
      for (int k = 0; k < m; k++) {
        if (M[i][k] != 0) {
          for (int j = 0; j < m; j++) {
            MR[i][j] += M[i][k] * Ri[k][j];
          }
        }
      }
    }
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) {
        int64_t v = 0;
        for (int k = 0; k < m; k++) {
          v += MR[i][k] * M[j][k];
        }
        G(i, j) = v / P2;
      }
    }
    return std::make_pair(G, MR);
  };
  std::vector<std::vector<int64_t>> frontier;
  {
    std::vector<std::vector<int64_t>> M0(m, std::vector<int64_t>(m, 0));
    for (int i = 0; i < m; i++) {
      M0[i][i] = P;
    }
    frontier.push_back(HermiteNormalFormInt64(M0, m));
  }
  for (int level = 0; level < h; level++) {
    std::unordered_set<std::vector<int64_t>, std::hash<std::vector<int64_t>>>
        seen;
    std::vector<std::vector<int64_t>> next;
    for (auto const &key : frontier) {
      std::vector<std::vector<int64_t>> M = unpack(key);
      auto [G, MR] = gram_of(M);
      for (auto const &c : ProjectiveKernelModP<int64_t, int64_t>(G, p)) {
        std::vector<int64_t> v(m, 0);
        for (int i = 0; i < m; i++) {
          if (c[i] != 0) {
            for (int j = 0; j < m; j++) {
              v[j] += c[i] * M[i][j];
            }
          }
        }
        bool ok = true;
        for (int j = 0; j < m && ok; j++) {
          if (v[j] % p != 0) {
            ok = false;
          }
          v[j] /= p;
        }
        if (!ok) {
          continue;
        }
        // integrality with Lambda and evenness of the new vector
        std::vector<int64_t> vR(m, 0);
        for (int k = 0; k < m; k++) {
          if (v[k] != 0) {
            for (int j = 0; j < m; j++) {
              vR[j] += v[k] * Ri[k][j];
            }
          }
        }
        for (int i = 0; i < m && ok; i++) {
          int64_t s = 0;
          for (int j = 0; j < m; j++) {
            s += vR[j] * M[i][j];
          }
          if (s % P2 != 0) {
            ok = false;
          }
        }
        int64_t nv = 0;
        for (int j = 0; j < m; j++) {
          nv += vR[j] * v[j];
        }
        if (!ok || nv % (P2 * mod_even) != 0) {
          continue;
        }
        std::vector<std::vector<int64_t>> Mext = M;
        Mext.push_back(v);
        std::vector<int64_t> knew = HermiteNormalFormInt64(Mext, m);
        if (seen.insert(knew).second) {
          next.push_back(std::move(knew));
        }
      }
    }
    frontier = std::move(next);
  }
  int v_det = 0;
  {
    int64_t d = DeterminantInt64<T, Tint>(R);
    while (d % p == 0) {
      d /= p;
      v_det++;
    }
  }
  std::map<std::vector<int>, int64_t> counts;
  for (auto const &key : frontier) {
    MyMatrix<int64_t> G = gram_of(unpack(key)).first;
    std::optional<std::vector<int>> opt =
        LocalSymbolKeyModular(G, p, v_det - 2 * h);
    std::vector<int> sym =
        opt ? *opt
            : LocalSymbolKey<T, Tint>(UniversalMatrixConversion<T, int64_t>(G), p);
    counts[sym]++;
  }
  return counts;
}

/*
  M_S = sum_{L in G} r(S, L) / |Aut L| for S of full rank n: an embedding
  S -> L identifies L with an overlattice of S, and counting the pairs (L,
  embedding) with weight 1 / |Aut L| gives exactly the number of even
  overlattices of S, in S tensor Q, that lie in the genus G. An overlattice
  is given prime by prime, so the count is the product over the primes p
  dividing det S of the number of overlattices of p-power index (the index
  sqrt(det S / det T)_p) locally isometric to T at p.
 */
template <typename T, typename Tint>
T FullRankRepresentationMass(MyMatrix<T> const &S, MyMatrix<T> const &Tgram,
                             FullRankOverlatticeCache<T> &cache) {
  int64_t detS = DeterminantInt64<T, Tint>(S);
  int64_t detT = DeterminantInt64<T, Tint>(Tgram);
  if (detT == 0 || detS % detT != 0) {
    return T(0);
  }
  int64_t ratio = detS / detT;
  T count(1);
  for (int64_t p : BadPrimes(1, detS)) {
    int v = 0;
    int64_t w = ratio;
    while (w % p == 0) {
      w /= p;
      v++;
    }
    if (v % 2 != 0) {
      return T(0);
    }
    if (detS % p != 0) {
      // p = 2 with det S odd: both even unimodular at 2, of the same class
      continue;
    }
    auto &by_S = cache.counts_at_index[{p, v / 2}];
    auto it_S = by_S.find(S);
    if (it_S == by_S.end()) {
      it_S = by_S.emplace(S, OverlatticeClassCountsAtIndex<T, Tint>(S, p, v / 2))
                 .first;
    }
    std::vector<int> keyT = LocalSymbolKey<T, Tint>(Tgram, p);
    auto iter = it_S->second.find(keyT);
    int64_t n_p = (iter == it_S->second.end()) ? 0 : iter->second;
#ifdef SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE_EXTENSIVE
    std::cerr << "LDOV: SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE_EXTENSIVE, the "
                 "full rank count against all the overlattices\n";
    auto const &sc = OverlatticeSymbolCount<T, Tint>(S, p, cache);
    auto iter_all = sc.find({v / 2, keyT});
    int64_t n_all = (iter_all == sc.end()) ? 0 : iter_all->second;
    if (n_all != n_p) {
      std::cerr << "LDOV: " << n_p << " overlattices at the index, " << n_all
                << " among all\n";
      throw TerminalException{1};
    }
    std::cerr << "LDOV: SANITY_CHECK_LOCAL_DENSITY_OVERLATTICE_EXTENSIVE, done\n";
#endif
    if (n_p == 0) {
      return T(0);
    }
    count *= ScalarFromInt64<T>(n_p);
  }
  return count;
}

/*
  M_S / mass(G) for the genus G of the even lattice Tgram of even rank n and
  a positive definite even S of rank m with n - m >= 1, from the local
  densities alone:
      c * (generic part) * prod_{p | 2 det S det T} alpha_p(S, T),
  alpha_p from the Siegel series where T is unimodular at p, from the unit
  chain and the reduction formula at the odd primes dividing det T, and
  from the gluing into an even unimodular lattice at p = 2 when det T is
  even, with c = 2, or c = 1 when n - m = 1. For n = m see
  FullRankRepresentationMass.
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

// alpha_p(S, T) for T the genus of the data, in the conventions of
// AlphaUnimodularTarget, at a prime p | 2 det S det T.
template <typename T, typename Tint>
T LocalDensityAtPrime(MyMatrix<T> const &S, GenusDensityDataFull<T> const &data,
                      int64_t p, LocalDensityCache<T> &cache) {
  int64_t detT = data.odd.detT;
  if (detT % p != 0) {
    int k = data.odd.n / 2;
    int xi_T = KroneckerSymbol(SignedPower(k, detT), p);
    return AlphaUnimodularTarget<T, Tint>(S, p, k, xi_T, cache.siegel);
  }
  if (p == 2) {
    return OverlatticeDensityFromData<T, Tint>(S, *data.two_adic, cache.siegel);
  }
  return DensityOddDiagonal<T>(OddDiagonalClasses<T, Tint>(S, p),
                               data.odd.odd_diagonal.at(p), p, cache.bad);
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
  if (n - m < 1) {
    std::cerr << "LDOV: the mass ratio needs n - m >= 1\n";
    throw TerminalException{1};
  }
  int64_t detS = DeterminantInt64<T, Tint>(S);
  int64_t detT = data.odd.detT;
  // Siegel's constant: 2 in the conventions of LocalDensity.h, halved when
  // n = m + 1 (checked on the certified representation masses of eleven odd
  // and sixteen even rank-14 genera, ranks 12 and 13)
  T eps = (n - m == 1) ? T(1) : T(2);
  T val = eps * GenericDensityPart<T>(n, m, detS, detT);
  for (int64_t p : BadPrimes(detS, detT)) {
    val *= LocalDensityAtPrime<T, Tint>(S, data, p, cache);
  }
  return val;
}


// clang-format off
#endif  // SRC_GENERA_LOCALDENSITYOVERLATTICE_H_
// clang-format on
