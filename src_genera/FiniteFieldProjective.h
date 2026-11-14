// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_FINITEFIELDPROJECTIVE_H_
#define SRC_GENERA_FINITEFIELDPROJECTIVE_H_

// clang-format off
#include "PolynomialModQ.h"
#include "GeneraArithmetic.h"
#include <cstdint>
#include <iostream>
#include <vector>
// clang-format on

/*
  Projective spaces over F_p and over F_q = F_p[theta]/mu, indexed without
  being enumerated, and the arithmetic of F_q: on coefficient vectors, and
  on the integer codes of the elements through tables for small q.
 */

/*
  The points of the projective space P^{n-1}(F_p), indexed WITHOUT being
  enumerated. There are (p^n - 1) / (p - 1) of them: 16383 for n = 14 and
  p = 2, but 1.5 x 10^9 at p = 5, where materializing the list (as an
  earlier version did) is itself a hundred gigabytes. The representatives
  are the vectors whose first non-zero coordinate is 1; the index runs
  through the blocks ipos = 0, ..., n-1 (position of that coordinate), each
  block of size p^{n-ipos-1}, the free coordinates read little-endian
  (v[ipos+1] least significant). ProjectiveLineOfIndex and
  IndexOfProjectiveLine are mutually inverse bijections for this order.
 */
inline size_t NumberProjectiveLines(int n, int p) {
  size_t const max_val = std::numeric_limits<size_t>::max();
  size_t p_s = static_cast<size_t>(p);
  size_t total = 0;
  size_t power = 1;
  for (int k = 0; k < n; k++) {
    if (total > max_val - power) {
      std::cerr << "GENUS: (p^n - 1)/(p - 1) overflows size_t for p=" << p
                << " n=" << n << "\n";
      throw TerminalException{1};
    }
    total += power;
    if (k + 1 < n) {
      if (power > max_val / p_s) {
        std::cerr << "GENUS: p^" << (k + 1) << " overflows size_t for p=" << p
                  << "\n";
        throw TerminalException{1};
      }
      power *= p_s;
    }
  }
  return total;
}

inline std::vector<int> ProjectiveLineOfIndex(int n, int p, size_t index) {
  size_t p_s = static_cast<size_t>(p);
  // block size p^{n-1} of the block ipos = 0
  size_t block = MyPow<size_t>(p_s, n - 1);
  size_t work = index;
  int ipos = 0;
  while (work >= block) {
    work -= block;
    block /= p_s;
    ipos++;
#ifdef SANITY_CHECK_GENUS_ENUMERATION
    if (ipos == n) {
      std::cerr << "GENUS: ProjectiveLineOfIndex, index out of range\n";
      throw TerminalException{1};
    }
#endif
  }
  std::vector<int> v(n, 0);
  v[ipos] = 1;
  for (int i = ipos + 1; i < n; i++) {
    v[i] = static_cast<int>(work % p_s);
    work /= p_s;
  }
  return v;
}

inline size_t IndexOfProjectiveLine(std::vector<int> const &v, int p) {
  size_t p_s = static_cast<size_t>(p);
  int n = static_cast<int>(v.size());
  int ipos = 0;
  while (ipos < n && ((v[ipos] % p) + p) % p == 0) {
    ipos++;
  }
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (ipos == n) {
    std::cerr << "GENUS: IndexOfProjectiveLine called on the zero vector\n";
    throw TerminalException{1};
  }
#endif
  int scal = ResInv(v[ipos], p);
  // start of the block: sum_{k < ipos} p^{n-k-1}
  size_t block = MyPow<size_t>(p_s, n - 1);
  size_t index = 0;
  for (int k = 0; k < ipos; k++) {
    index += block;
    block /= p_s;
  }
  // offset inside the block, little-endian in the normalized coordinates
  size_t mult = 1;
  for (int i = ipos + 1; i < n; i++) {
    int val = ((((v[i] % p) + p) % p) * scal) % p;
    index += mult * static_cast<size_t>(val);
    mult *= p_s;
  }
  return index;
}

// Kept only for the sanity check of the indexing, as an INDEPENDENT
// enumeration of the same order (small cases only); the enumeration itself
// never materializes the projective space.
inline std::vector<std::vector<int>> EnumerateProjectiveLines(int n, int p) {
  std::vector<std::vector<int>> ListLines;
  for (int ipos = 0; ipos < n; ipos++) {
    int n_free = n - ipos - 1;
    int64_t n_case = MyPow<int64_t>(p, n_free);
    for (int64_t idx = 0; idx < n_case; idx++) {
      std::vector<int> v(n, 0);
      v[ipos] = 1;
      int64_t work = idx;
      for (int i = ipos + 1; i < n; i++) {
        v[i] = static_cast<int>(work % p);
        work /= p;
      }
      ListLines.push_back(v);
    }
  }
  return ListLines;
}

#ifdef SANITY_CHECK_GENUS_ENUMERATION
// The indexing against the independent enumeration: exhaustive when the
// space is small, and for a large space the roundtrip index -> line ->
// index on the block boundaries and on a stride through the range.
inline void CheckProjectiveLineIndexing(int n, int p) {
  size_t n_line = NumberProjectiveLines(n, p);
  if (n_line <= 100000) {
    std::vector<std::vector<int>> ListLines = EnumerateProjectiveLines(n, p);
    if (ListLines.size() != n_line) {
      std::cerr << "GENUS: NumberProjectiveLines disagrees with the "
                   "enumeration\n";
      throw TerminalException{1};
    }
    for (size_t i = 0; i < n_line; i++) {
      if (ProjectiveLineOfIndex(n, p, i) != ListLines[i] ||
          IndexOfProjectiveLine(ListLines[i], p) != i) {
        std::cerr << "GENUS: projective line indexing broken at index " << i
                  << "\n";
        throw TerminalException{1};
      }
    }
    return;
  }
  auto roundtrip = [&](size_t i) -> void {
    std::vector<int> v = ProjectiveLineOfIndex(n, p, i);
    if (IndexOfProjectiveLine(v, p) != i) {
      std::cerr << "GENUS: projective line roundtrip broken at index " << i
                << "\n";
      throw TerminalException{1};
    }
  };
  size_t block_start = 0;
  size_t block = n_line - (n_line - 1) / static_cast<size_t>(p);
  // block = p^{n-1}; the boundaries of every block and their neighbors
  for (int ipos = 0; ipos < n; ipos++) {
    roundtrip(block_start);
    if (block_start + block - 1 < n_line) {
      roundtrip(block_start + block - 1);
    }
    block_start += block;
    block /= static_cast<size_t>(p);
  }
  size_t stride = 1 + n_line / 10000;
  for (size_t i = 0; i < n_line; i += stride) {
    roundtrip(i);
  }
}
#endif

/*
  F_q = F_p[theta] / mu: elements as coefficient vectors of length f, and
  an encoding as integers in [0, q).
 */
struct FiniteFieldModP {
  int64_t p;
  int f;
  PolyQ mu;
  int64_t q;
};

inline std::vector<int64_t> FqMul(FiniteFieldModP const &F,
                                  std::vector<int64_t> const &a,
                                  std::vector<int64_t> const &b) {
  PolyQ prod = PolyMulQ(a, b, F.p);
  std::pair<PolyQ, PolyQ> qr = PolyDivModQ(prod, F.mu, F.p);
  std::vector<int64_t> r = qr.second;
  r.resize(F.f, 0);
  return r;
}

inline std::vector<int64_t> FqAdd(FiniteFieldModP const &F,
                                  std::vector<int64_t> const &a,
                                  std::vector<int64_t> const &b) {
  std::vector<int64_t> r(F.f);
  for (int i = 0; i < F.f; i++) {
    r[i] = (a[i] + b[i]) % F.p;
  }
  return r;
}

inline bool FqIsZero(std::vector<int64_t> const &a) {
  return std::all_of(a.begin(), a.end(), [](int64_t x) { return x == 0; });
}

inline std::vector<int64_t> FqInverse(FiniteFieldModP const &F,
                                      std::vector<int64_t> const &a) {
  PolyQ a2 = a;
  PolyTrim(a2);
  PolyQ inv = PolyInverseModQ(a2, F.mu, F.p);
  inv.resize(F.f, 0);
  return inv;
}

inline int64_t FqEncode(FiniteFieldModP const &F,
                        std::vector<int64_t> const &a) {
  int64_t c = 0;
  for (int i = F.f; i-- > 0;) {
    c = c * F.p + a[i];
  }
  return c;
}

inline std::vector<int64_t> FqDecode(FiniteFieldModP const &F, int64_t c) {
  std::vector<int64_t> a(F.f);
  for (int i = 0; i < F.f; i++) {
    a[i] = c % F.p;
    c /= F.p;
  }
  return a;
}

/*
  The points of P^{m-1}(F_q), normalized with the first nonzero coordinate
  equal to 1, indexed by the position i of that coordinate (q^{m-1-i}
  points each) and the remaining coordinates in base q.
 */
inline size_t NumberProjectivePointsFq(int m, int64_t q) {
  // the points are indexed by size_t: a count beyond it is refused, as in
  // NumberProjectiveLines, rather than wrapped around
  size_t const max_val = std::numeric_limits<size_t>::max();
  size_t qs = static_cast<size_t>(q);
  size_t tot = 0, pw = 1;
  for (int i = 0; i < m; i++) {
    if (tot > max_val - pw) {
      std::cerr << "GENUS: |P^" << (m - 1) << "(F_" << q
                << ")| overflows size_t\n";
      throw TerminalException{1};
    }
    tot += pw;
    if (i + 1 < m) {
      if (pw > max_val / qs) {
        std::cerr << "GENUS: " << q << "^" << (i + 1) << " overflows size_t\n";
        throw TerminalException{1};
      }
      pw *= qs;
    }
  }
  return tot;
}

inline std::vector<std::vector<int64_t>>
ProjectivePointFqOfIndex(FiniteFieldModP const &F, int m, size_t idx) {
  // the blocks in the order i = m-1 (one point) up to i = 0
  size_t q = static_cast<size_t>(F.q);
  std::vector<size_t> block(m);
  size_t pw = 1;
  for (int i = m - 1; i >= 0; i--) {
    block[i] = pw;
    pw *= q;
  }
  int lead = m - 1;
  for (int i = m - 1; i >= 0; i--) {
    if (idx < block[i]) {
      lead = i;
      break;
    }
    idx -= block[i];
  }
  std::vector<std::vector<int64_t>> x(m, std::vector<int64_t>(F.f, 0));
  x[lead][0] = 1;
  for (int j = m - 1; j > lead; j--) {
    x[j] = FqDecode(F, static_cast<int64_t>(idx % q));
    idx /= q;
  }
  return x;
}

inline size_t IndexOfProjectivePointFq(FiniteFieldModP const &F, int m,
                                       std::vector<std::vector<int64_t>> x) {
  size_t q = static_cast<size_t>(F.q);
  int lead = 0;
  while (lead < m && FqIsZero(x[lead])) {
    lead++;
  }
  std::vector<int64_t> inv = FqInverse(F, x[lead]);
  for (int j = lead; j < m; j++) {
    x[j] = FqMul(F, x[j], inv);
  }
  size_t offset = 0, pw = 1;
  for (int i = m - 1; i > lead; i--) {
    offset += pw;
    pw *= q;
  }
  size_t idx = 0;
  for (int j = lead + 1; j < m; j++) {
    idx = idx * q + static_cast<size_t>(FqEncode(F, x[j]));
  }
  return offset + idx;
}

/*
  F_q for a small q as tables on the codes of FqEncode, and the projective
  points and the action on Hom of the orbit enumeration in those codes: the
  same indices as ProjectivePointFqOfIndex and IndexOfProjectivePointFq,
  without the polynomial arithmetic and the allocations of every product
  (the orbits of P^11(F_4), 5.6 10^6 points, under the centralizer of
  (A2^12, w) took about ten minutes per class).
 */
struct FqTables {
  int64_t q;
  std::vector<int> add;
  std::vector<int> mul;
  std::vector<int> inv;
};

inline size_t const fq_tables_max_q = 4096;

inline FqTables BuildFqTables(FiniteFieldModP const &F) {
  int64_t q = F.q;
  FqTables tb{q, std::vector<int>(q * q), std::vector<int>(q * q),
              std::vector<int>(q, 0)};
  for (int64_t a = 0; a < q; a++) {
    std::vector<int64_t> fa = FqDecode(F, a);
    for (int64_t b = 0; b < q; b++) {
      std::vector<int64_t> fb = FqDecode(F, b);
      tb.add[a * q + b] = static_cast<int>(FqEncode(F, FqAdd(F, fa, fb)));
      int c = static_cast<int>(FqEncode(F, FqMul(F, fa, fb)));
      tb.mul[a * q + b] = c;
      if (c == 1) {
        tb.inv[a] = static_cast<int>(b);
      }
    }
  }
  return tb;
}

// The point of index idx as codes, x of size m.
inline void ProjectivePointCodesOfIndex(FqTables const &tb, int m, size_t idx,
                                        std::vector<int> &x) {
  size_t q = static_cast<size_t>(tb.q);
  size_t pw = 1;
  // the block of the lead i holds q^{m-1-i} points, in the order i = m-1
  // (one point) up to i = 0
  int lead = m - 1;
  for (int i = m - 1; i >= 0; i--) {
    if (idx < pw) {
      lead = i;
      break;
    }
    idx -= pw;
    pw *= q;
  }
  for (int j = 0; j < lead; j++) {
    x[j] = 0;
  }
  x[lead] = 1;
  for (int j = m - 1; j > lead; j--) {
    x[j] = static_cast<int>(idx % q);
    idx /= q;
  }
}

// The index of the point of codes x, which is normalized in place.
inline size_t IndexOfProjectivePointCodes(FqTables const &tb, int m,
                                          std::vector<int> &x) {
  size_t q = static_cast<size_t>(tb.q);
  int lead = 0;
  while (lead < m && x[lead] == 0) {
    lead++;
  }
  int inv = tb.inv[x[lead]];
  for (int j = lead; j < m; j++) {
    x[j] = tb.mul[x[j] * tb.q + inv];
  }
  size_t offset = 0, pw = 1;
  for (int i = m - 1; i > lead; i--) {
    offset += pw;
    pw *= q;
  }
  size_t idx = 0;
  for (int j = lead + 1; j < m; j++) {
    idx = idx * q + static_cast<size_t>(x[j]);
  }
  return offset + idx;
}

// The action of ApplyActionOnHom in codes: cc[i * m + j] the code of c[i][j].
inline std::vector<int>
ActionOnHomCodes(FiniteFieldModP const &F,
                 std::vector<std::vector<std::vector<int64_t>>> const &c) {
  size_t m = c.size();
  std::vector<int> cc(m * m);
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < m; j++) {
      cc[i * m + j] = static_cast<int>(FqEncode(F, c[i][j]));
    }
  }
  return cc;
}

inline void ApplyActionOnHomCodes(FqTables const &tb, int m,
                                  std::vector<int> const &cc,
                                  std::vector<int> const &alpha,
                                  std::vector<int> &beta) {
  for (int j = 0; j < m; j++) {
    beta[j] = 0;
  }
  for (int i = 0; i < m; i++) {
    if (alpha[i] == 0) {
      continue;
    }
    int const *row = cc.data() + i * m;
    int const *mul_a = tb.mul.data() + alpha[i] * tb.q;
    for (int j = 0; j < m; j++) {
      beta[j] = tb.add[beta[j] * tb.q + mul_a[row[j]]];
    }
  }
}

// clang-format off
#endif  // SRC_GENERA_FINITEFIELDPROJECTIVE_H_
// clang-format on
