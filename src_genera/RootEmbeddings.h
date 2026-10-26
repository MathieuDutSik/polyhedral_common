// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ROOTEMBEDDINGS_H_
#define SRC_GENERA_ROOTEMBEDDINGS_H_

// clang-format off
#include "SiegelSpecialValues.h"
#include <algorithm>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

/*
  The number r(R, T) of isometric embeddings of a root lattice R into a root
  lattice T (ADE types, simply laced).

  An embedding sends the roots of R to roots of T, so r(R, L) = r(R, rho(L))
  for any lattice L with root system rho(L). With R = R_1 + R', R_1
  irreducible, an embedding is one of R_1 with image S, a connected root
  subsystem of T isometric to R_1, together with one of R' into the roots of
  T orthogonal to S. The |Aut R_1| embeddings with a given image leave the
  same orthogonal root system, so

      r(R_1 + R', T) = |Aut R_1| sum_S r(R', Phi_S^perp),

  and only the multiset of the types of Phi_S^perp matters: the split table
  of T, which factors over the irreducible components of T since a connected
  S lies in one of them. The split table of an irreducible type is computed
  by a walk over its connected root subsystems.

  The image of an embedding is reflection-closed: the norm 2 vectors of the
  root lattice of a root system are exactly its roots, so the roots of T in
  the integral span of S are those of S. A subsystem is therefore closed
  under the reflections s_a(b) = b - (a.b) a, which is all the walk needs.

  The type of a connected subsystem with N positive roots is read off its
  Coxeter number h, every root meeting 2h - 4 others up to sign: the rank is
  2N / h, and (rank, N) separates A_r, D_r and E_6, E_7, E_8.
 */

// An ADE type: sorted (letter, rank) pairs, the empty one being the zero
// root system. D_2 and D_3 are written A_1 + A_1 and A_3.
using AdeType = std::vector<std::pair<char, int>>;

inline AdeType AdeCanonical(AdeType t) {
  std::sort(t.begin(), t.end());
  return t;
}

inline int AdeRank(AdeType const &t) {
  int r = 0;
  for (auto const &c : t) {
    r += c.second;
  }
  return r;
}

inline int64_t AdeNumberRootsIrreducible(char letter, int r) {
  if (letter == 'A') {
    return static_cast<int64_t>(r) * (r + 1);
  }
  if (letter == 'D') {
    return 2 * static_cast<int64_t>(r) * (r - 1);
  }
  if (r == 6) {
    return 72;
  }
  if (r == 7) {
    return 126;
  }
  return 240;
}

inline int64_t AdeNumberRoots(AdeType const &t) {
  int64_t n = 0;
  for (auto const &c : t) {
    n += AdeNumberRootsIrreducible(c.first, c.second);
  }
  return n;
}

inline std::string AdeTypeName(AdeType const &t) {
  if (t.empty()) {
    return "0";
  }
  std::string s;
  for (size_t i = 0; i < t.size(); i++) {
    if (i > 0) {
      s += "+";
    }
    s += t[i].first + std::to_string(t[i].second);
  }
  return s;
}

inline AdeType ParseAdeType(std::string const &name) {
  AdeType t;
  if (name == "0") {
    return t;
  }
  size_t pos = 0;
  while (pos < name.size()) {
    size_t next = name.find('+', pos);
    if (next == std::string::npos) {
      next = name.size();
    }
    std::string comp = name.substr(pos, next - pos);
    char letter = comp[0];
    int r = std::stoi(comp.substr(1));
    bool valid = (letter == 'A' && r >= 1) || (letter == 'D' && r >= 4) ||
                 (letter == 'E' && r >= 6 && r <= 8);
    if (!valid) {
      std::cerr << "ROOT: invalid ADE component " << comp << "\n";
      throw TerminalException{1};
    }
    t.push_back({letter, r});
    pos = next + 1;
  }
  return AdeCanonical(t);
}

// Every non-zero ADE type of rank at most max_rank.
inline std::vector<AdeType> AdeTypesUpToRank(int max_rank) {
  std::vector<std::pair<char, int>> comps;
  for (int r = 1; r <= max_rank; r++) {
    comps.push_back({'A', r});
  }
  for (int r = 4; r <= max_rank; r++) {
    comps.push_back({'D', r});
  }
  for (int r = 6; r <= 8 && r <= max_rank; r++) {
    comps.push_back({'E', r});
  }
  std::vector<AdeType> out;
  AdeType cur;
  std::function<void(size_t, int)> rec = [&](size_t start, int left) -> void {
    if (!cur.empty()) {
      out.push_back(AdeCanonical(cur));
    }
    for (size_t i = start; i < comps.size(); i++) {
      if (comps[i].second <= left) {
        cur.push_back(comps[i]);
        rec(i, left - comps[i].second);
        cur.pop_back();
      }
    }
  };
  rec(0, max_rank);
  return out;
}

// The Cartan matrix, E_n with the Bourbaki labelling (the chain
// 1-3-4-5-...-n and 2 on 4).
template <typename T> MyMatrix<T> CartanMatrixIrreducible(char letter, int r) {
  MyMatrix<T> M = ZeroMatrix<T>(r, r);
  for (int i = 0; i < r; i++) {
    M(i, i) = 2;
  }
  auto edge = [&](int a, int b) -> void {
    M(a, b) = -1;
    M(b, a) = -1;
  };
  if (letter == 'A') {
    for (int i = 0; i + 1 < r; i++) {
      edge(i, i + 1);
    }
  } else if (letter == 'D') {
    for (int i = 0; i + 2 < r; i++) {
      edge(i, i + 1);
    }
    edge(r - 3, r - 1);
  } else {
    edge(0, 2);
    edge(1, 3);
    for (int i = 2; i + 1 < r; i++) {
      edge(i, i + 1);
    }
  }
  return M;
}

template <typename T> MyMatrix<T> AdeGramMatrix(AdeType const &t) {
  int n = AdeRank(t);
  MyMatrix<T> M = ZeroMatrix<T>(n, n);
  int pos = 0;
  for (auto const &c : t) {
    MyMatrix<T> C = CartanMatrixIrreducible<T>(c.first, c.second);
    for (int i = 0; i < c.second; i++) {
      for (int j = 0; j < c.second; j++) {
        M(pos + i, pos + j) = C(i, j);
      }
    }
    pos += c.second;
  }
  return M;
}

// |Aut| of an irreducible root lattice: the Weyl group times the diagram
// automorphisms.
template <typename Tint> Tint AutOrderIrreducible(char letter, int r) {
  Tint fact(1);
  for (int i = 2; i <= r; i++) {
    fact *= ScalarFromInt64<Tint>(i);
  }
  if (letter == 'A') {
    return (r == 1) ? Tint(2) : Tint(2) * fact * ScalarFromInt64<Tint>(r + 1);
  }
  if (letter == 'D') {
    if (r == 4) {
      return ScalarFromInt64<Tint>(1152);
    }
    Tint pow2(1);
    for (int i = 0; i < r; i++) {
      pow2 *= Tint(2);
    }
    return pow2 * fact;
  }
  if (r == 6) {
    return ScalarFromInt64<Tint>(103680);
  }
  if (r == 7) {
    return ScalarFromInt64<Tint>(2903040);
  }
  return ScalarFromInt64<Tint>(696729600);
}

template <typename Tint> Tint AutOrderAde(AdeType const &t) {
  Tint res(1);
  size_t i = 0;
  while (i < t.size()) {
    size_t j = i;
    while (j < t.size() && t[j] == t[i]) {
      j++;
    }
    Tint a = AutOrderIrreducible<Tint>(t[i].first, t[i].second);
    for (size_t k = i; k < j; k++) {
      res *= a * ScalarFromInt64<Tint>(static_cast<int64_t>(k - i + 1));
    }
    i = j;
  }
  return res;
}

//
// The positive roots of an irreducible root system and the walk over its
// connected subsystems
//

using RootMask = std::vector<uint64_t>;

struct RootMaskHash {
  size_t operator()(RootMask const &m) const {
    size_t h = 0;
    for (uint64_t w : m) {
      h ^= std::hash<uint64_t>()(w) + 0x9e3779b97f4a7c15ULL + (h << 6) +
           (h >> 2);
    }
    return h;
  }
};

struct IrreducibleRootSystem {
  int rank;
  int n_pos;
  int n_words;
  // the positive roots in the basis of the simple roots
  std::vector<std::vector<int64_t>> coords;
  // inner products of the positive roots
  std::vector<std::vector<int>> ip;
  // reflect[a][b]: the positive root +-s_b(a), for ip[a][b] = +-1
  std::vector<std::vector<int>> reflect;
  // adj[a]: the roots b with ip[a][b] != 0, a included
  std::vector<RootMask> adj;
};

inline bool MaskHas(RootMask const &m, int i) {
  return (m[i >> 6] >> (i & 63)) & 1;
}

inline void MaskSet(RootMask &m, int i) { m[i >> 6] |= uint64_t(1) << (i & 63); }

inline std::vector<int> MaskElements(RootMask const &m, int n) {
  std::vector<int> res;
  for (int i = 0; i < n; i++) {
    if (MaskHas(m, i)) {
      res.push_back(i);
    }
  }
  return res;
}

inline IrreducibleRootSystem BuildIrreducibleRootSystem(char letter, int r) {
  MyMatrix<int> C = CartanMatrixIrreducible<int>(letter, r);
  // the positive roots in the basis of the simple roots, by increasing
  // height: b + a_i is a root when b.a_i = -1 (simply laced)
  std::vector<std::vector<int>> roots;
  std::map<std::vector<int>, int> index;
  for (int i = 0; i < r; i++) {
    std::vector<int> v(r, 0);
    v[i] = 1;
    index[v] = static_cast<int>(roots.size());
    roots.push_back(v);
  }
  auto dot = [&](std::vector<int> const &x, std::vector<int> const &y) -> int {
    int s = 0;
    for (int i = 0; i < r; i++) {
      if (x[i] != 0) {
        for (int j = 0; j < r; j++) {
          s += x[i] * C(i, j) * y[j];
        }
      }
    }
    return s;
  };
  for (size_t pos = 0; pos < roots.size(); pos++) {
    for (int i = 0; i < r; i++) {
      std::vector<int> e(r, 0);
      e[i] = 1;
      if (dot(roots[pos], e) == -1) {
        std::vector<int> w = roots[pos];
        w[i]++;
        if (index.count(w) == 0) {
          index[w] = static_cast<int>(roots.size());
          roots.push_back(w);
        }
      }
    }
  }
  IrreducibleRootSystem rs;
  rs.rank = r;
  rs.n_pos = static_cast<int>(roots.size());
  for (auto const &v : roots) {
    rs.coords.push_back(std::vector<int64_t>(v.begin(), v.end()));
  }
  rs.n_words = (rs.n_pos + 63) / 64;
  rs.ip.assign(rs.n_pos, std::vector<int>(rs.n_pos));
  rs.reflect.assign(rs.n_pos, std::vector<int>(rs.n_pos, -1));
  rs.adj.assign(rs.n_pos, RootMask(rs.n_words, 0));
  for (int a = 0; a < rs.n_pos; a++) {
    for (int b = 0; b < rs.n_pos; b++) {
      rs.ip[a][b] = dot(roots[a], roots[b]);
    }
  }
  for (int a = 0; a < rs.n_pos; a++) {
    for (int b = 0; b < rs.n_pos; b++) {
      int c = rs.ip[a][b];
      if (c != 0) {
        MaskSet(rs.adj[a], b);
      }
      if (c == 1 || c == -1) {
        std::vector<int> w(r);
        bool neg = false;
        for (int i = 0; i < r; i++) {
          w[i] = roots[a][i] - c * roots[b][i];
          if (w[i] < 0) {
            neg = true;
          }
        }
        if (neg) {
          for (int i = 0; i < r; i++) {
            w[i] = -w[i];
          }
        }
        rs.reflect[a][b] = index.at(w);
      }
    }
  }
  return rs;
}

// The reflection closure of a closed subsystem and one more root q.
inline RootMask RootClosure(IrreducibleRootSystem const &rs, RootMask mask,
                            int q) {
  std::vector<int> elems = MaskElements(mask, rs.n_pos);
  std::vector<int> fresh{q};
  MaskSet(mask, q);
  elems.push_back(q);
  while (!fresh.empty()) {
    int a = fresh.back();
    fresh.pop_back();
    size_t n_elem = elems.size();
    for (size_t k = 0; k < n_elem; k++) {
      int b = elems[k];
      int c = rs.reflect[a][b];
      if (c >= 0 && !MaskHas(mask, c)) {
        MaskSet(mask, c);
        elems.push_back(c);
        fresh.push_back(c);
      }
    }
  }
  return mask;
}

// The ADE type of the subsystem given by a set of positive roots.
inline AdeType SubsystemType(IrreducibleRootSystem const &rs,
                             std::vector<int> const &elems) {
  int m = static_cast<int>(elems.size());
  std::vector<int> comp(m, -1);
  AdeType t;
  for (int i = 0; i < m; i++) {
    if (comp[i] >= 0) {
      continue;
    }
    std::vector<int> stack{i};
    std::vector<int> members{i};
    comp[i] = i;
    while (!stack.empty()) {
      int a = stack.back();
      stack.pop_back();
      for (int b = 0; b < m; b++) {
        if (comp[b] < 0 && rs.ip[elems[a]][elems[b]] != 0) {
          comp[b] = i;
          stack.push_back(b);
          members.push_back(b);
        }
      }
    }
    int64_t N = static_cast<int64_t>(members.size());
    int64_t d = 0;
    for (int b : members) {
      if (b != i && rs.ip[elems[i]][elems[b]] != 0) {
        d++;
      }
    }
    // 2 h - 4 = d, rank = 2 N / h
    int64_t h2 = d + 4;
    int64_t rank = (4 * N) / h2;
    char letter = 0;
    if (2 * N == rank * (rank + 1)) {
      letter = 'A';
    } else if (N == rank * (rank - 1)) {
      letter = 'D';
    } else if ((rank == 6 && N == 36) || (rank == 7 && N == 63) ||
               (rank == 8 && N == 120)) {
      letter = 'E';
    }
    if (letter == 0 || 4 * N != rank * h2) {
      std::cerr << "ROOT: no ADE type with " << N << " positive roots and "
                << "Coxeter number " << h2 / 2 << "\n";
      throw TerminalException{1};
    }
    t.push_back({letter, static_cast<int>(rank)});
  }
  return AdeCanonical(t);
}

// A row Hermite normal form of the coordinates of a few roots (at most the
// rank, entries bounded by the highest root), in int64.
struct SmallHermiteForm {
  std::vector<std::vector<int64_t>> rows;
  std::vector<int> pivots;
};

inline SmallHermiteForm ComputeSmallHermiteForm(IrreducibleRootSystem const &rs,
                                                std::vector<int> const &basis) {
  std::vector<std::vector<int64_t>> M;
  for (int b : basis) {
    M.push_back(rs.coords[b]);
  }
  int nr = static_cast<int>(M.size());
  int nc = rs.rank;
  SmallHermiteForm H;
  int row = 0;
  for (int c = 0; c < nc && row < nr; c++) {
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
        for (int j = c; j < nc; j++) {
          M[row][j] -= q * M[i][j];
        }
        std::swap(M[row], M[i]);
      }
    }
    if (M[row][c] < 0) {
      for (int j = c; j < nc; j++) {
        M[row][j] = -M[row][j];
      }
    }
    H.rows.push_back(M[row]);
    H.pivots.push_back(c);
    row++;
  }
  return H;
}

// The representative of v modulo the span, pivot entries in [0, h), then
// the sign making the first non-zero entry positive.
inline std::vector<int64_t> ReduceModSpanUpToSign(SmallHermiteForm const &H,
                                                  std::vector<int64_t> v) {
  for (size_t i = 0; i < H.rows.size(); i++) {
    int c = H.pivots[i];
    int64_t h = H.rows[i][c];
    int64_t q = v[c] / h;
    if (v[c] - q * h < 0) {
      q--;
    }
    if (q != 0) {
      for (size_t j = c; j < v.size(); j++) {
        v[j] -= q * H.rows[i][j];
      }
    }
  }
  std::vector<int64_t> w = v;
  for (int64_t &x : w) {
    x = -x;
  }
  for (size_t i = 0; i < H.rows.size(); i++) {
    int c = H.pivots[i];
    int64_t h = H.rows[i][c];
    int64_t q = w[c] / h;
    if (w[c] - q * h < 0) {
      q--;
    }
    if (q != 0) {
      for (size_t j = c; j < w.size(); j++) {
        w[j] -= q * H.rows[i][j];
      }
    }
  }
  return std::min(v, w);
}

// For an irreducible type: for every irreducible R_1, the number of
// subsystems isometric to R_1 with each type of orthogonal root system.
using SplitTable = std::map<AdeType, std::map<AdeType, int64_t>>;

inline SplitTable ComputeSplitTable(char letter, int r) {
  IrreducibleRootSystem rs = BuildIrreducibleRootSystem(letter, r);
  int n = rs.n_pos;
  SplitTable table;
  auto record = [&](RootMask const &mask) -> void {
    std::vector<int> elems = MaskElements(mask, n);
    AdeType ty = SubsystemType(rs, elems);
    RootMask meet(rs.n_words, 0);
    for (int a : elems) {
      for (int w = 0; w < rs.n_words; w++) {
        meet[w] |= rs.adj[a][w];
      }
    }
    std::vector<int> perp;
    for (int q = 0; q < n; q++) {
      if (!MaskHas(meet, q)) {
        perp.push_back(q);
      }
    }
    table[ty][SubsystemType(rs, perp)]++;
  };
  // a subsystem is carried with a basis of its integral span, at most rank
  // roots
  std::unordered_set<RootMask, RootMaskHash> seen;
  std::vector<std::pair<RootMask, std::vector<int>>> frontier;
  for (int q = 0; q < n; q++) {
    RootMask m(rs.n_words, 0);
    MaskSet(m, q);
    seen.insert(m);
    frontier.push_back({m, {q}});
    record(m);
  }
  while (!frontier.empty()) {
    std::vector<std::pair<RootMask, std::vector<int>>> next;
    for (auto const &[mask, basis] : frontier) {
      RootMask meet(rs.n_words, 0);
      for (int a : MaskElements(mask, n)) {
        for (int w = 0; w < rs.n_words; w++) {
          meet[w] |= rs.adj[a][w];
        }
      }
      SmallHermiteForm H = ComputeSmallHermiteForm(rs, basis);
      // roots q, q' with q' = +-q modulo the span of the subsystem generate
      // the same child; one closure per class
      std::map<std::vector<int64_t>, int> classes;
      for (int q = 0; q < n; q++) {
        if (MaskHas(mask, q) || !MaskHas(meet, q)) {
          continue;
        }
        classes.emplace(ReduceModSpanUpToSign(H, rs.coords[q]), q);
      }
      for (auto const &kv : classes) {
        int q = kv.second;
        RootMask child = RootClosure(rs, mask, q);
        if (seen.insert(child).second) {
          record(child);
          std::vector<int> child_basis = basis;
          child_basis.push_back(q);
          next.push_back({std::move(child), std::move(child_basis)});
        }
      }
    }
    frontier = std::move(next);
  }
  return table;
}

template <typename Tint> struct RootEmbeddingCache {
  std::map<std::pair<char, int>, SplitTable> split;
  std::map<std::pair<AdeType, AdeType>, Tint> embeddings;
};

template <typename Tint>
SplitTable const &GetSplitTable(char letter, int r,
                                RootEmbeddingCache<Tint> &cache) {
  auto iter = cache.split.find({letter, r});
  if (iter != cache.split.end()) {
    return iter->second;
  }
  return cache.split.emplace(std::make_pair(letter, r),
                             ComputeSplitTable(letter, r))
      .first->second;
}

// r(R, T), the number of isometric embeddings of the root lattice R into
// the root lattice T.
template <typename Tint>
Tint RootEmbeddingNumber(AdeType const &R, AdeType const &T,
                         RootEmbeddingCache<Tint> &cache) {
  if (R.empty()) {
    return Tint(1);
  }
  if (AdeRank(R) > AdeRank(T) || AdeNumberRoots(R) > AdeNumberRoots(T)) {
    return Tint(0);
  }
  auto key = std::make_pair(R, T);
  auto iter = cache.embeddings.find(key);
  if (iter != cache.embeddings.end()) {
    return iter->second;
  }
  AdeType R1{R[0]};
  AdeType Rrest(R.begin() + 1, R.end());
  Tint total(0);
  for (size_t i = 0; i < T.size(); i++) {
    if (i > 0 && T[i] == T[i - 1]) {
      // an isomorphic component gives the same sum
      continue;
    }
    size_t mult = 1;
    while (i + mult < T.size() && T[i + mult] == T[i]) {
      mult++;
    }
    if (T[i].second < R1[0].second) {
      continue;
    }
    SplitTable const &sp = GetSplitTable(T[i].first, T[i].second, cache);
    auto it_r = sp.find(R1);
    if (it_r == sp.end()) {
      continue;
    }
    AdeType others;
    for (size_t j = 0; j < T.size(); j++) {
      if (j != i) {
        others.push_back(T[j]);
      }
    }
    for (auto const &[C, count] : it_r->second) {
      AdeType Tnew = others;
      Tnew.insert(Tnew.end(), C.begin(), C.end());
      Tint val = RootEmbeddingNumber(Rrest, AdeCanonical(Tnew), cache);
      total += ScalarFromInt64<Tint>(count) *
               ScalarFromInt64<Tint>(static_cast<int64_t>(mult)) * val;
    }
  }
  total *= AutOrderIrreducible<Tint>(R1[0].first, R1[0].second);
  cache.embeddings[key] = total;
  return total;
}

// clang-format off
#endif  // SRC_GENERA_ROOTEMBEDDINGS_H_
// clang-format on
