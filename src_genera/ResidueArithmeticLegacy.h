// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_RESIDUEARITHMETICLEGACY_H_
#define SRC_GENERA_RESIDUEARITHMETICLEGACY_H_

// clang-format off
#include "MAT_MatrixModPrime.h"
#include "ResidueArithmetic.h"
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  TRANSITIONAL. The signed residue helpers and the matrices of residues as
  rows of int64_t, under the names *I64, for the headers of src_genera not
  yet converted to T_residue, T_small_prime and ModMatrix = MyMatrix of
  ResidueArithmetic.h and MAT_MatrixModPrime.h. Every use is a place left
  to convert; the header goes when none remains.
 */

// a mod m in [0, m), for any a.
inline int64_t ResModI64(int64_t a, int64_t m) {
  int64_t r = a % m;
  return (r < 0) ? r + m : r;
}

// x mod m in [0, m), for an unbounded integer x.
template <typename Tint> int64_t ResModOfI64(Tint const &x, int64_t m) {
  Tint m_T = UniversalScalarConversion<Tint, int64_t>(m);
  Tint r = ResInt(x, m_T);
  if (r < Tint(0)) {
    r += m_T;
  }
  return UniversalScalarConversion<int64_t, Tint>(r);
}

// The operations on residues a, b in [0, m).
inline int64_t ResAddI64(int64_t a, int64_t b, int64_t m) {
  int64_t s = a + b;
  return (s >= m) ? s - m : s;
}

inline int64_t ResSubI64(int64_t a, int64_t b, int64_t m) {
  int64_t s = a - b;
  return (s < 0) ? s + m : s;
}

inline int64_t ResNegI64(int64_t a, int64_t m) { return (a == 0) ? 0 : m - a; }

// a b mod m in [0, m), for any signs of a and b.
inline int64_t ResMulI64(int64_t a, int64_t b, int64_t m) {
  __int128 prod = static_cast<__int128>(a) * static_cast<__int128>(b);
  int64_t r = static_cast<int64_t>(prod % static_cast<__int128>(m));
  return (r < 0) ? r + m : r;
}

// a^e mod m in [0, m), e >= 0, by binary exponentiation.
inline int64_t ResPowI64(int64_t a, int64_t e, int64_t m) {
  int64_t res = ResModI64(1, m);
  a = ResModI64(a, m);
  while (e > 0) {
    if (e & 1) {
      res = ResMulI64(res, a, m);
    }
    a = ResMulI64(a, a, m);
    e >>= 1;
  }
  return res;
}

// The inverse of a modulo m, an error when a is not invertible. The
// coefficients of the extended Euclid algorithm are bounded by m.
inline int64_t ResInvI64(int64_t a, int64_t m) {
  int64_t r0 = m, r1 = ResModI64(a, m);
  int64_t s0 = 0, s1 = 1;
  while (r1 != 0) {
    int64_t q = r0 / r1;
    int64_t r2 = r0 - q * r1;
    r0 = r1;
    r1 = r2;
    int64_t s2 = s0 - q * s1;
    s0 = s1;
    s1 = s2;
  }
  if (r0 != 1) {
    std::cerr << "RESIDUE: " << a << " is not invertible modulo " << m << "\n";
    throw TerminalException{1};
  }
  return ResModI64(s0, m);
}

using ModRowsI64 = std::vector<std::vector<int64_t>>;

// The conversions between the two representations.
inline ModRowsI64 ModRowsOfModMatrixI64(ModMatrix const &X) {
  ModRowsI64 R(X.rows(), std::vector<int64_t>(X.cols()));
  for (int i = 0; i < X.rows(); i++) {
    for (int j = 0; j < X.cols(); j++) {
      R[i][j] = ResToSigned(X(i, j));
    }
  }
  return R;
}

inline ModMatrix ModMatrixOfModRowsI64(ModRowsI64 const &R, int n_col) {
  ModMatrix X(static_cast<int>(R.size()), n_col);
  for (size_t i = 0; i < R.size(); i++) {
    for (int j = 0; j < n_col; j++) {
      X(static_cast<int>(i), j) = static_cast<T_residue>(R[i][j]);
    }
  }
  return X;
}

// X Y modulo mod.
inline ModRowsI64 ModRowsProductI64(ModRowsI64 const &X, ModRowsI64 const &Y,
                                  int64_t mod) {
  size_t n = X.size();
  size_t k = Y.size();
  size_t l = Y.empty() ? 0 : Y[0].size();
  ModRowsI64 Z(n, std::vector<int64_t>(l, 0));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < k; j++) {
      int64_t x = ResModI64(X[i][j], mod);
      if (x == 0) {
        continue;
      }
      for (size_t u = 0; u < l; u++) {
        Z[i][u] = ResAddI64(Z[i][u], ResMulI64(x, ResModI64(Y[j][u], mod), mod), mod);
      }
    }
  }
  return Z;
}

inline ModRowsI64 ModRowsIdentityI64(size_t n) {
  ModRowsI64 Id(n, std::vector<int64_t>(n, 0));
  for (size_t i = 0; i < n; i++) {
    Id[i][i] = 1;
  }
  return Id;
}

/*
  The reduced row echelon form of the row space modulo the prime p: the
  basis rows and their pivot columns. The coordinates of a vector of the
  row space in that basis are its entries at the pivots.
 */
struct RowSpaceModPI64 {
  ModRowsI64 basis;
  std::vector<size_t> pivots;
};

inline RowSpaceModPI64 ComputeRowSpaceModPI64(ModRowsI64 rows, int64_t p) {
  RowSpaceModPI64 res;
  size_t n_col = rows.empty() ? 0 : rows[0].size();
  size_t r = 0;
  for (auto &row : rows) {
    for (auto &x : row) {
      x = ResModI64(x, p);
    }
  }
  for (size_t col = 0; col < n_col && r < rows.size(); col++) {
    size_t piv = r;
    while (piv < rows.size() && rows[piv][col] == 0) {
      piv++;
    }
    if (piv == rows.size()) {
      continue;
    }
    std::swap(rows[r], rows[piv]);
    int64_t inv = ResInvI64(rows[r][col], p);
    for (size_t j = col; j < n_col; j++) {
      rows[r][j] = ResMulI64(rows[r][j], inv, p);
    }
    for (size_t i = 0; i < rows.size(); i++) {
      int64_t c = rows[i][col];
      if (i != r && c != 0) {
        for (size_t j = col; j < n_col; j++) {
          if (rows[r][j] != 0) {
            rows[i][j] = ResSubI64(rows[i][j], ResMulI64(c, rows[r][j], p), p);
          }
        }
      }
    }
    res.pivots.push_back(col);
    r++;
  }
  rows.resize(r);
  res.basis = std::move(rows);
  return res;
}

// The rank modulo the prime p.
inline size_t RankModPrimeI64(ModRowsI64 const &rows, int64_t p) {
  return ComputeRowSpaceModPI64(rows, p).pivots.size();
}

// The coordinates of a vector of the row space in its echelon basis.
inline std::vector<int64_t>
CoordinatesInRowSpaceModPI64(RowSpaceModPI64 const &rs,
                          std::vector<int64_t> const &v) {
  std::vector<int64_t> c(rs.pivots.size());
  for (size_t i = 0; i < rs.pivots.size(); i++) {
    c[i] = v[rs.pivots[i]];
  }
  return c;
}

// A basis of the x with x M = 0 modulo p (M with r rows).
inline ModRowsI64 LeftNullspaceModPI64(ModRowsI64 const &M, int64_t p) {
  size_t r = M.size();
  size_t c = M.empty() ? 0 : M[0].size();
  ModRowsI64 aug(r, std::vector<int64_t>(c + r, 0));
  for (size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < c; j++) {
      aug[i][j] = ResModI64(M[i][j], p);
    }
    aug[i][c + i] = 1;
  }
  RowSpaceModPI64 rs = ComputeRowSpaceModPI64(std::move(aug), p);
  ModRowsI64 res;
  for (size_t i = 0; i < rs.basis.size(); i++) {
    if (rs.pivots[i] >= c) {
      res.push_back(std::vector<int64_t>(rs.basis[i].begin() + c,
                                         rs.basis[i].end()));
    }
  }
  return res;
}

// The solution x of x B = v modulo p (B with independent rows), if any.
inline std::optional<std::vector<int64_t>>
SolveRowsModPI64(ModRowsI64 const &B, std::vector<int64_t> const &v, int64_t p) {
  size_t r = B.size();
  ModRowsI64 M = B;
  M.push_back(v);
  ModRowsI64 ker = LeftNullspaceModPI64(M, p);
  for (auto const &x : ker) {
    if (x[r] != 0) {
      int64_t inv = ResInvI64(ResNegI64(x[r], p), p);
      std::vector<int64_t> sol(r);
      for (size_t i = 0; i < r; i++) {
        sol[i] = ResMulI64(x[i], inv, p);
      }
      return sol;
    }
  }
  return {};
}


// clang-format off
#endif  // SRC_GENERA_RESIDUEARITHMETICLEGACY_H_
// clang-format on
