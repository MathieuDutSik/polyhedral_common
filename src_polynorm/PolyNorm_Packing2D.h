// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_PACKING2D_H_
#define SRC_POLYNORM_POLYNORM_PACKING2D_H_

// clang-format off
#include "PolyNorm_LatticeOpt.h"
#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_PACKING2D
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYNORM_PACKING2D
#endif

/*
  The densest lattice packing of a convex polygon P, exactly.

  The lattice packings of P are the admissible lattices of D = P - P, so
  the density is vol(P) / Delta(D) with Delta(D) the critical determinant.
  For a centrally symmetric convex disk D, Minkowski showed that a critical
  lattice has a basis b1, b2 with b1, b2 and b1 + b2 on the boundary of D,
  and conversely that a lattice with a basis having b1, b2, b1 + b2 on the
  boundary is admissible (the hexagon +-b1, +-b2, +-(b1 + b2) inscribed in
  D has no other lattice point; Gruber and Lekkerkerker, Geometry of
  Numbers, section 2.16; Mount and Silverman, J. Algorithms 1990).

  So Delta(D) = min |det(b1, b2)| over b1, b2, b1 + b2 on the boundary.
  For a polygon, put b1 on the edge j1 (b1 = A1 + s E1, s in [0, 1]), b2
  on the edge j2 (b2 = A2 + t E2) and b1 + b2 on the edge j3. The last
  condition is one linear equation in (s, t) and two inequalities (the
  neighbors of j3), so the feasible set is a segment, along which det is a
  quadratic polynomial. Its minimum is at an end of the segment or at the
  vertex of the parabola, which are rational for a rational polygon. The
  minimum over all triples of edges is Delta(D), and the computation is
  exact and global.
 */

template <typename T> struct PolyNormPacking2D {
  // The critical determinant of D and the basis (rows) of a critical
  // lattice.
  T Delta;
  MyMatrix<T> W;
  // The number of triples of edges whose segment is not empty.
  size_t n_feasible;
};

template <typename T>
T PolyNorm_Det2(MyVector<T> const &u, MyVector<T> const &v) {
  return u(0) * v(1) - u(1) * v(0);
}

template <typename T, typename Tint>
PolyNormPacking2D<T>
PolyNorm_ExactPacking2D(MyMatrix<T> const &LmatDiff,
                        [[maybe_unused]] std::ostream &os) {
  int m = LmatDiff.rows();
  if (LmatDiff.cols() != 2) {
    std::cerr << "POLYNORM_PACKING2D: the polygon should be 2-dimensional\n";
    throw TerminalException{1};
  }
  // The edges in angular order of their normals.
  std::vector<int> ord(m);
  for (int j = 0; j < m; j++) {
    ord[j] = j;
  }
  std::vector<double> angle(m);
  for (int j = 0; j < m; j++) {
    angle[j] = std::atan2(PolyNorm_ToDouble(LmatDiff(j, 1)),
                          PolyNorm_ToDouble(LmatDiff(j, 0)));
  }
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) -> bool { return angle[a] < angle[b]; });
  std::vector<MyVector<T>> ListNormal(m);
  for (int u = 0; u < m; u++) {
    ListNormal[u] = GetMatrixRow(LmatDiff, ord[u]);
  }
  // The vertex between the edges u and u+1.
  auto vertex = [&](int u, int v) -> MyVector<T> {
    MyMatrix<T> M(2, 2);
    for (int k = 0; k < 2; k++) {
      M(0, k) = ListNormal[u](k);
      M(1, k) = ListNormal[v](k);
    }
    MyVector<T> one(2);
    one(0) = 1;
    one(1) = 1;
    return Inverse(M) * one;
  };
  std::vector<MyVector<T>> ListA(m), ListE(m);
  for (int u = 0; u < m; u++) {
    int u_prev = (u + m - 1) % m;
    int u_next = (u + 1) % m;
    MyVector<T> A = vertex(u_prev, u);
    MyVector<T> B = vertex(u, u_next);
    ListA[u] = A;
    ListE[u] = B - A;
  }
  // compat(u, v): some point x of the edge v has ell_u(x) <= 0, required
  // for b2 when b1 is on u, since ell_u(b1 + b2) <= 1 = ell_u(b1).
  std::vector<std::vector<uint8_t>> compat(m, std::vector<uint8_t>(m, 0));
  for (int u = 0; u < m; u++) {
    for (int v = 0; v < m; v++) {
      T val1 = ListNormal[u].dot(ListA[v]);
      T val2 = ListNormal[u].dot(ListA[v] + ListE[v]);
      compat[u][v] = (val1 <= 0 || val2 <= 0) ? 1 : 0;
    }
  }
  std::optional<T> best;
  MyMatrix<T> W_best(2, 2);
  size_t n_feasible = 0;
  for (int u1 = 0; u1 < m; u1++) {
    for (int u2 = 0; u2 < m; u2++) {
      if (!compat[u1][u2] || !compat[u2][u1]) {
        continue;
      }
      for (int u3 = 0; u3 < m; u3++) {
        MyVector<T> const &d3 = ListNormal[u3];
        // b1 + b2 = A1 + A2 + s E1 + t E2 on the line ell_3 = 1:
        // alpha s + beta t = gamma.
        T alpha = d3.dot(ListE[u1]);
        T beta = d3.dot(ListE[u2]);
        T gamma = 1 - d3.dot(ListA[u1] + ListA[u2]);
        // Parametrization (s, t) = P + x Q.
        MyVector<T> P(2), Q(2);
        if (beta != 0) {
          P(0) = 0;
          P(1) = gamma / beta;
          Q(0) = 1;
          Q(1) = -alpha / beta;
        } else if (alpha != 0) {
          P(0) = gamma / alpha;
          P(1) = 0;
          Q(0) = 0;
          Q(1) = 1;
        } else {
          // Three parallel edges: ell_3(b1 + b2) is 2 or 0, never 1.
#ifdef SANITY_CHECK_POLYNORM_PACKING2D
          if (gamma == 0) {
            std::cerr << "POLYNORM_PACKING2D: degenerate triple of edges\n";
            throw TerminalException{1};
          }
#endif
          continue;
        }
        // The interval of x: the inequalities c0 + c1 x >= 0.
        std::vector<std::pair<T, T>> ListIneq;
        ListIneq.push_back({P(0), Q(0)});
        ListIneq.push_back({1 - P(0), -Q(0)});
        ListIneq.push_back({P(1), Q(1)});
        ListIneq.push_back({1 - P(1), -Q(1)});
        for (int w : {(u3 + m - 1) % m, (u3 + 1) % m}) {
          MyVector<T> const &dw = ListNormal[w];
          // 1 - ell_w(A1 + A2 + s E1 + t E2) >= 0.
          T c0 = 1 - dw.dot(ListA[u1] + ListA[u2]) - P(0) * dw.dot(ListE[u1]) -
                 P(1) * dw.dot(ListE[u2]);
          T c1 = -Q(0) * dw.dot(ListE[u1]) - Q(1) * dw.dot(ListE[u2]);
          ListIneq.push_back({c0, c1});
        }
        std::optional<T> lo, hi;
        bool feasible = true;
        for (auto &eIneq : ListIneq) {
          T const &c0 = eIneq.first;
          T const &c1 = eIneq.second;
          if (c1 == 0) {
            if (c0 < 0) {
              feasible = false;
            }
          } else if (c1 > 0) {
            T x = -c0 / c1;
            if (!lo || x > *lo) {
              lo = x;
            }
          } else {
            T x = -c0 / c1;
            if (!hi || x < *hi) {
              hi = x;
            }
          }
        }
        if (!feasible || !lo || !hi || *lo > *hi) {
          continue;
        }
        n_feasible++;
        // b1 = R1 + x S1, b2 = R2 + x S2.
        MyVector<T> R1 = ListA[u1] + P(0) * ListE[u1];
        MyVector<T> S1 = Q(0) * ListE[u1];
        MyVector<T> R2 = ListA[u2] + P(1) * ListE[u2];
        MyVector<T> S2 = Q(1) * ListE[u2];
        T qa = PolyNorm_Det2(S1, S2);
        T qb = PolyNorm_Det2(R1, S2) + PolyNorm_Det2(S1, R2);
        T qc = PolyNorm_Det2(R1, R2);
        auto qval = [&](T const &x) -> T { return qa * x * x + qb * x + qc; };
        T mid = (*lo + *hi) / 2;
        T val_mid = qval(mid);
        if (val_mid == 0) {
          // b1 and b2 dependent, which the boundary conditions exclude.
#ifdef SANITY_CHECK_POLYNORM_PACKING2D
          std::cerr << "POLYNORM_PACKING2D: degenerate basis on a segment\n";
          throw TerminalException{1};
#endif
          continue;
        }
        T sign = val_mid > 0 ? T(1) : T(-1);
        std::vector<T> ListX{*lo, *hi};
        if (sign * qa > 0) {
          T x = -qb / (2 * qa);
          if (*lo <= x && x <= *hi) {
            ListX.push_back(x);
          }
        }
        for (auto &x : ListX) {
          T val = sign * qval(x);
          if (!best || val < *best) {
            best = val;
            MyVector<T> b1 = R1 + x * S1;
            MyVector<T> b2 = R2 + x * S2;
            for (int k = 0; k < 2; k++) {
              W_best(0, k) = b1(k);
              W_best(1, k) = b2(k);
            }
          }
        }
      }
    }
  }
  if (!best) {
    std::cerr << "POLYNORM_PACKING2D: no feasible configuration found\n";
    throw TerminalException{1};
  }
  if (DeterminantMat(W_best) < 0) {
    for (int k = 0; k < 2; k++) {
      W_best(0, k) = -W_best(0, k);
    }
  }
#ifdef SANITY_CHECK_POLYNORM_PACKING2D
  PolyNormLatticePackingEval<T, Tint> eval =
      PolyNorm_EvaluatePacking<T, Tint>(LmatDiff, W_best, os);
  if (eval.alpha != 1 || eval.det != *best) {
    std::cerr << "POLYNORM_PACKING2D: the critical lattice is not admissible"
              << " alpha=" << eval.alpha << " det=" << eval.det
              << " Delta=" << *best << "\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_POLYNORM_PACKING2D
  os << "POLYNORM_PACKING2D: m=" << m << " n_feasible=" << n_feasible
     << " Delta=" << *best << "\n";
#endif
  return {*best, W_best, n_feasible};
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_PACKING2D_H_
// clang-format on
