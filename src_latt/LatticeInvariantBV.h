// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEINVARIANTBV_H_
#define SRC_LATT_LATTICEINVARIANTBV_H_

// clang-format off
#include "Shvec_exact.h"
#include "Timings.h"
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_LATTICE_INVARIANT_BV
#endif

#ifdef TIMINGS
#define TIMINGS_LATTICE_INVARIANT_BV
#endif

/*
  An isometry invariant of a positive definite integral lattice in the
  manner of the BV invariant of Chenevier and Taibi (Unimodular lattices of
  rank 29 and related even genera of small determinant, arXiv:2601.19780;
  fast_marked_HBV in their eqfminim.c), at a cost cubic in the number N of
  short vectors.

  1. S = the vectors of norm at most M, up to sign, M being the least norm
     bound for which |S| >= min_vectors (an isometry invariant choice).
  2. A_ij = f(b(v_i, v_j)) for an even function f, so that A does not
     depend on the signs of the representatives:
     * parity: b mod 2, the choice of Chenevier-Taibi, natural for odd
       lattices; A is a bit matrix and (A^2)_ij = popcount(row_i & row_j),
       N^3 / 64 word operations;
     * absolute: |b|, as their absBVdual;
     * modular: the class of +-b mod q, coded by min(b mod q, -b mod q).
     The last two use a dense N x N product in int64.
  3. For every v_i, the sorted list of the pairs (A_ij, (A^2)_ij) over j
     is hashed; the sorted list of these N hashes is hashed again.
  The result, with M and N, is an isometry invariant: an isometry permutes S
  and preserves A. It is not complete. Two lattices with different
  invariants are not isometric, while equal invariants prove nothing (and
  Chenevier and Taibi tune the variant per genus to separate their classes).
 */

enum class BVFunction { Parity, Absolute, Modular };

struct LatticeBVOptions {
  BVFunction function = BVFunction::Parity;
  // the q of BVFunction::Modular
  int64_t modulus = 3;
  // the least number of vectors up to sign; 0 for twice the rank
  size_t min_vectors = 0;
  // the largest: beyond it the computation stops with an error rather than
  // take a cubic time on a huge set
  size_t max_vectors = 20000;
};

inline BVFunction ParseBVFunction(std::string const &str) {
  if (str == "parity") {
    return BVFunction::Parity;
  }
  if (str == "absolute") {
    return BVFunction::Absolute;
  }
  if (str == "modular") {
    return BVFunction::Modular;
  }
  std::cerr << "BV: unknown function " << str
            << ", allowed: parity, absolute, modular\n";
  throw TerminalException{1};
}

template <typename T> struct LatticeBVInvariant {
  T max_norm;
  size_t n_vect;
  uint64_t hash;
};

template <typename T>
bool operator<(LatticeBVInvariant<T> const &a, LatticeBVInvariant<T> const &b) {
  if (a.max_norm != b.max_norm) {
    return a.max_norm < b.max_norm;
  }
  if (a.n_vect != b.n_vect) {
    return a.n_vect < b.n_vect;
  }
  return a.hash < b.hash;
}

template <typename T>
bool operator==(LatticeBVInvariant<T> const &a,
                LatticeBVInvariant<T> const &b) {
  return a.max_norm == b.max_norm && a.n_vect == b.n_vect && a.hash == b.hash;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, LatticeBVInvariant<T> const &x) {
  os << "BV(M=" << x.max_norm << ", N=" << x.n_vect << ", h=" << x.hash
     << ")";
  return os;
}

// The mixing step of the hashes of Chenevier-Taibi.
inline uint64_t BVHashMix(uint64_t h, uint64_t x) {
  return (h ^ x) * UINT64_C(1111111111111111111);
}

inline uint64_t BVHashSeed() { return UINT64_C(13282407956253574712); }

// The hash of the sorted rows of the pairs (A_ij, (A^2)_ij), then of the
// sorted row hashes.
inline uint64_t BVHashRows(std::vector<std::vector<uint64_t>> &rows) {
  std::vector<uint64_t> hashes;
  hashes.reserve(rows.size());
  for (auto &row : rows) {
    std::sort(row.begin(), row.end());
    uint64_t h = BVHashSeed();
    for (uint64_t x : row) {
      h = BVHashMix(h, x);
    }
    hashes.push_back(h);
  }
  std::sort(hashes.begin(), hashes.end());
  uint64_t h = BVHashSeed();
  for (uint64_t x : hashes) {
    h = BVHashMix(h, x);
  }
  return h;
}

// A pair (a, s) as one word, a in the high bits: a < 2^24, s < 2^40.
inline uint64_t BVPair(uint64_t a, uint64_t s) { return (a << 40) | s; }

/*
  The vectors of S, up to sign, as rows, and M: the norm bounds are raised
  through the integers from the minimum until there are min_vectors of them.
 */
template <typename T, typename Tint>
std::pair<T, std::vector<MyVector<Tint>>>
BVShortVectors(MyMatrix<T> const &GramMat, LatticeBVOptions const &opt,
               std::ostream &os) {
  size_t min_vectors =
      (opt.min_vectors == 0) ? size_t(2 * GramMat.rows()) : opt.min_vectors;
  CVPSolver<T, Tint> solver(GramMat, os);
  T M = solver.shortest_vectors().min;
  while (true) {
    std::vector<MyVector<Tint>> SHV = solver.at_most_norm_vectors(M);
    if (SHV.size() > opt.max_vectors) {
      std::cerr << "BV: " << SHV.size() << " vectors of norm at most " << M
                << " above max_vectors=" << opt.max_vectors << "\n";
      throw TerminalException{1};
    }
    if (SHV.size() >= min_vectors) {
      return {M, std::move(SHV)};
    }
    M += T(1);
  }
}

template <typename T, typename Tint>
LatticeBVInvariant<T> LatticeBV(MyMatrix<T> const &GramMat,
                                LatticeBVOptions const &opt,
                                std::ostream &os) {
#ifdef TIMINGS_LATTICE_INVARIANT_BV
  MicrosecondTime time;
#endif
  if (!IsIntegralMatrix(GramMat)) {
    std::cerr << "BV: the Gram matrix must be integral\n";
    throw TerminalException{1};
  }
  int n = GramMat.rows();
  auto [M, SHV] = BVShortVectors<T, Tint>(GramMat, opt, os);
  size_t N = SHV.size();
#ifdef TIMINGS_LATTICE_INVARIANT_BV
  os << "|BV: short vectors|=" << time << " N=" << N << "\n";
#endif
  MyMatrix<int64_t> G = UniversalMatrixConversion<int64_t, T>(GramMat);
  // w_i = G v_i in int64: the inner products are bounded by M
  std::vector<std::vector<int64_t>> V(N, std::vector<int64_t>(n));
  std::vector<std::vector<int64_t>> W(N, std::vector<int64_t>(n, 0));
  for (size_t i = 0; i < N; i++) {
    for (int a = 0; a < n; a++) {
      V[i][a] = UniversalScalarConversion<int64_t, Tint>(SHV[i](a));
    }
    for (int a = 0; a < n; a++) {
      int64_t s = 0;
      for (int b = 0; b < n; b++) {
        s += G(a, b) * V[i][b];
      }
      W[i][a] = s;
    }
  }
  auto inner = [&](size_t i, size_t j) -> int64_t {
    int64_t s = 0;
    for (int a = 0; a < n; a++) {
      s += W[i][a] * V[j][a];
    }
    return s;
  };
  std::vector<std::vector<uint64_t>> rows(N, std::vector<uint64_t>(N));
  if (opt.function == BVFunction::Parity) {
    // A as bit rows of length N
    size_t n_word = (N + 63) / 64;
    std::vector<uint64_t> A(N * n_word, 0);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        if ((inner(i, j) & 1) != 0) {
          A[i * n_word + j / 64] |= uint64_t(1) << (j % 64);
        }
      }
    }
    for (size_t i = 0; i < N; i++) {
      uint64_t const *ri = A.data() + i * n_word;
      for (size_t j = 0; j < N; j++) {
        uint64_t const *rj = A.data() + j * n_word;
        uint64_t s = 0;
        for (size_t k = 0; k < n_word; k++) {
          s += static_cast<uint64_t>(__builtin_popcountll(ri[k] & rj[k]));
        }
        uint64_t a = (ri[j / 64] >> (j % 64)) & 1;
        rows[i][j] = BVPair(a, s);
      }
    }
  } else {
    int64_t q = opt.modulus;
    if (opt.function == BVFunction::Modular && q < 2) {
      std::cerr << "BV: the modulus must be at least 2\n";
      throw TerminalException{1};
    }
    auto f = [&](int64_t x) -> int64_t {
      if (opt.function == BVFunction::Absolute) {
        return (x < 0) ? -x : x;
      }
      int64_t r = ((x % q) + q) % q;
      return std::min(r, q - r);
    };
    std::vector<int64_t> A(N * N);
    for (size_t i = 0; i < N; i++) {
      for (size_t j = i; j < N; j++) {
        int64_t val = f(inner(i, j));
        A[i * N + j] = val;
        A[j * N + i] = val;
      }
    }
    for (size_t i = 0; i < N; i++) {
      int64_t const *ri = A.data() + i * N;
      for (size_t j = 0; j < N; j++) {
        // A is symmetric: (A^2)_ij is the product of the rows i and j
        int64_t const *rj = A.data() + j * N;
        int64_t s = 0;
        for (size_t k = 0; k < N; k++) {
          s += ri[k] * rj[k];
        }
        rows[i][j] = BVPair(static_cast<uint64_t>(ri[j]),
                            static_cast<uint64_t>(s));
      }
    }
  }
  LatticeBVInvariant<T> res{M, N, BVHashRows(rows)};
#ifdef TIMINGS_LATTICE_INVARIANT_BV
  os << "|BV: LatticeBV|=" << time << "\n";
#endif
#ifdef DEBUG_LATTICE_INVARIANT_BV
  os << "BV: " << res << "\n";
#endif
  return res;
}

// clang-format off
#endif  // SRC_LATT_LATTICEINVARIANTBV_H_
// clang-format on
