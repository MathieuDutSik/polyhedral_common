// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ERDAHL_ERDAHL_FLIP_H_
#define SRC_ERDAHL_ERDAHL_FLIP_H_

// clang-format off
#include "erdahl_polyhedron.h"
#include <cmath>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ERDAHL_FLIP
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ERDAHL_FLIP
#endif

#ifdef TIMINGS
#define TIMINGS_ERDAHL_FLIP
#endif

/*
  The flipping algorithm (Algorithm "Flipping" of the paper, the GAP
  function InfDel_LiftingDelaunay).

  Take Delaunay polyhedra D1 subset D2 subset D3 of perfection ranks r+2,
  r+1, r relative to the space W, with functions f1, f2, f3 of W. Then
  Space_W(D1) = Space_W(D3) + R f1 + R f2 and the functions of W vanishing
  on D1 and nonnegative on D3 are, modulo Space_W(D3), the pencil
      a f1 + b f2  with  a f1(v) + b f2(v) >= 0 for v in D3.
  Since f1, f2 >= 0 on Z^n, and f2 = 0 < f1 on D2 - D1, this is a two
  dimensional cone in the half plane a >= 0. One extreme ray is (0, 1),
  that is f2 and D2. The other is (1, -lambda) with
      lambda = min { f1(v) / f2(v) : v in D3, f2(v) > 0 }
  and gives the Delaunay polyhedron D2' != D2 with D1 subset D2' subset D3.

  The minimum over the infinite set D3 is obtained iteratively: lambda is
  computed over a finite set of candidates, and g = f1 - lambda f2 is
  tested for nonnegativity on D3 by closest vector problems. A negative
  point v has f2(v) > 0 and f1(v)/f2(v) < lambda; it is added and the
  process repeated. The candidates also contain directions v of L(D3),
  whose ratio Quad(f1)[v] / Quad(f2)[v] is the limit of the ratio along v:
  when the answer is the value lambda_psd where the quadratic part of g on
  L(D3) stops being positive semidefinite, the points only approach it and
  the directions are needed (see erdahl_pencil_bound). The process
  terminates since the cone of functions of W vanishing on D1 and
  nonnegative on D3 is polyhedral.

  Everything stays in W since f1, f2, f3 are in W.
 */

/*
  Some points of D3 with f2 > 0, to start the iteration: the
  representatives of D3 and their translates by the basis vectors of L(D3)
  and their opposites, growing the translations until one point is found.
 */
template <typename T, typename Tint>
std::vector<MyVector<Tint>>
erdahl_flip_initial_points(DelaunayPolyhedron<T, Tint> const &D2,
                           DelaunayPolyhedron<T, Tint> const &D3) {
  std::vector<MyVector<Tint>> l_pt;
  int n = erdahl_dimension(D3);
  auto insert = [&](MyVector<Tint> const &e) -> void {
    if (EvaluationQuadForm<T, Tint>(D2.F, e) > 0) {
      l_pt.push_back(e);
    }
  };
  for (int i = 0; i < D3.EXT.rows(); i++) {
    insert(GetMatrixRow(D3.EXT, i));
  }
  Tint mult(1);
  while (l_pt.empty()) {
    for (int i = 0; i < D3.EXT.rows(); i++) {
      MyVector<Tint> e = GetMatrixRow(D3.EXT, i);
      for (int j = 0; j < D3.L.rows(); j++) {
        for (int sign = -1; sign <= 1; sign += 2) {
          MyVector<Tint> f = e;
          for (int k = 0; k < n; k++) {
            f(1 + k) += sign * mult * D3.L(j, k);
          }
          insert(f);
        }
      }
    }
    if (D3.L.rows() == 0 && l_pt.empty()) {
      std::cerr << "ERDAHL: D2 = D3, which is not allowed in the flip\n";
      throw TerminalException{1};
    }
    mult *= 2;
  }
  return l_pt;
}

/*
  The rational numbers p/q approximating x given by the convergents of its
  continued fraction, with q <= max_denom.
 */
template <typename T>
std::vector<T> erdahl_continued_fraction_convergents(double x, long max_denom) {
  std::vector<T> l_conv;
  long p_prev = 1, q_prev = 0;
  long p_curr = static_cast<long>(std::floor(x)), q_curr = 1;
  double frac = x - std::floor(x);
  l_conv.push_back(T(p_curr));
  for (int iter = 0; iter < 60; iter++) {
    if (std::abs(frac) < 1e-15) {
      break;
    }
    double inv = 1 / frac;
    long a = static_cast<long>(std::floor(inv));
    frac = inv - std::floor(inv);
    long p_next = a * p_curr + p_prev;
    long q_next = a * q_curr + q_prev;
    if (q_next > max_denom || q_next <= 0) {
      break;
    }
    p_prev = p_curr;
    q_prev = q_curr;
    p_curr = p_next;
    q_curr = q_next;
    l_conv.push_back(T(p_curr) / T(q_curr));
  }
  return l_conv;
}

/*
  The pencil Q1 - lambda Q2 of two positive semidefinite forms on Z^d is
  positive semidefinite exactly for lambda <= lambda_psd. For g = f1 -
  lambda f2 to be nonnegative on D3, lambda <= lambda_psd is needed for the
  quadratic parts restricted to L(D3). When lambda_psd is rational it is
  attained by integral directions v with Q1[v] = lambda_psd Q2[v], which
  are the eigen-directions of the GAP code (InfDel_GetEigenRationalConditions).
 */
template <typename T, typename Tint> struct ErdahlPencilBound {
  // False if Q1 - lambda Q2 is positive semidefinite for all lambda.
  bool bounded;
  // Whether lambda_psd is rational (and then equal to lambda).
  bool is_rational;
  T lambda;
  // Directions of Z^d realizing lambda_psd when rational.
  std::vector<MyVector<Tint>> directions;
};

template <typename T, typename Tint>
ErdahlPencilBound<T, Tint> erdahl_pencil_bound(MyMatrix<T> const &Q1,
                                               MyMatrix<T> const &Q2,
                                               std::ostream &os) {
  int d = Q1.rows();
  auto integral_kernel = [&](MyMatrix<T> const &M) -> MyMatrix<Tint> {
    MyMatrix<T> Mred = RemoveFractionMatrix(M);
    MyMatrix<Tint> M_int = UniversalMatrixConversion<Tint, T>(Mred);
    return NullspaceIntMat(M_int);
  };
  // On the kernel K1 of Q1 the pencil is - lambda Q2.
  MyMatrix<Tint> K1 = integral_kernel(Q1);
  std::vector<MyVector<Tint>> l_dir0;
  for (int i = 0; i < K1.rows(); i++) {
    MyVector<Tint> k = GetMatrixRow(K1, i);
    if (EvaluationQuadForm<T, Tint>(Q2, k) > 0) {
      l_dir0.push_back(k);
    }
  }
  if (!l_dir0.empty()) {
    return {true, true, T(0), l_dir0};
  }
  // Q2 vanishes on K1, hence Q2 K1 = 0, and the pencil lives on a
  // complement of K1 where Q1 is positive definite.
  MyMatrix<Tint> C1 = SubspaceCompletionInt(K1, d);
  int dc = C1.rows();
  if (dc == 0) {
    return {false, false, T(0), {}};
  }
  MyMatrix<T> C1_T = UniversalMatrixConversion<T, Tint>(C1);
  MyMatrix<T> Q1c = C1_T * Q1 * C1_T.transpose();
  MyMatrix<T> Q2c = C1_T * Q2 * C1_T.transpose();
  if (IsZeroMatrix(Q2c)) {
    return {false, false, T(0), {}};
  }
  // The largest mu with Q2c - mu Q1c singular, lambda_psd = 1 / mu.
  MyMatrix<double> Q1d = UniversalMatrixConversion<double, T>(Q1c);
  MyMatrix<double> Q2d = UniversalMatrixConversion<double, T>(Q2c);
  Eigen::GeneralizedSelfAdjointEigenSolver<MyMatrix<double>> eig(Q2d, Q1d);
  double mu_max = eig.eigenvalues().maxCoeff();
  for (auto &mu : erdahl_continued_fraction_convergents<T>(mu_max, 1000000000)) {
    if (mu <= 0) {
      continue;
    }
    MyMatrix<T> N = Q2c - mu * Q1c;
    if (RankMat(N) < dc && IsPositiveSemiDefinite(MyMatrix<T>(-N), os)) {
      MyMatrix<Tint> K = integral_kernel(N);
      std::vector<MyVector<Tint>> l_dir;
      for (int i = 0; i < K.rows(); i++) {
        MyVector<Tint> k = GetMatrixRow(K, i);
        l_dir.push_back(C1.transpose() * k);
      }
      return {true, true, T(1) / mu, l_dir};
    }
  }
  return {true, false, T(0), {}};
}

// The quadratic form of F restricted to the lattice L: L Quad(F) L^T.
template <typename T, typename Tint>
MyMatrix<T> erdahl_quad_on_lattice(MyMatrix<T> const &F,
                                   MyMatrix<Tint> const &L) {
  MyMatrix<T> L_T = UniversalMatrixConversion<T, Tint>(L);
  return L_T * erdahl_get_quad(F) * L_T.transpose();
}

/*
  With compute_function = false the function of the result is left empty,
  to be computed by erdahl_ensure_function if the result is kept.
 */
template <typename T, typename Tint>
DelaunayPolyhedron<T, Tint>
erdahl_flip(ErdahlFunctionSpace<T> const &W,
            DelaunayPolyhedron<T, Tint> const &D1,
            DelaunayPolyhedron<T, Tint> const &D2,
            DelaunayPolyhedron<T, Tint> const &D3, bool compute_function,
            std::ostream &os) {
#ifdef SANITY_CHECK_ERDAHL_FLIP
  if (!erdahl_is_subset(D1, D2) || !erdahl_is_subset(D2, D3)) {
    std::cerr << "ERDAHL: the flip requires D1 subset D2 subset D3\n";
    throw TerminalException{1};
  }
  int r3 = erdahl_perfection_rank(W, D3);
  int r2 = erdahl_perfection_rank(W, D2);
  int r1 = erdahl_perfection_rank(W, D1);
  if (r2 != r3 + 1 || r1 != r3 + 2) {
    std::cerr << "ERDAHL: the flip requires ranks r+2, r+1, r but got r1="
              << r1 << " r2=" << r2 << " r3=" << r3 << "\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS_ERDAHL_FLIP
  MicrosecondTime time_total;
  size_t n_iter = 0;
#endif
  ErdahlLatticeSet<Tint> ls3 = erdahl_lattice_set(D3);
  /*
    The candidates for the minimum ratio: points of D3 (rows (1, x)) with
    f2 > 0 and directions of L(D3) (rows (0, v)) with Quad(f2)[v] > 0,
    whose ratio is the limit of the ratio along the direction.
   */
  std::vector<MyVector<Tint>> l_cand = erdahl_flip_initial_points(D2, D3);
  MyMatrix<T> Q1L = erdahl_quad_on_lattice(D1.F, D3.L);
  MyMatrix<T> Q2L = erdahl_quad_on_lattice(D2.F, D3.L);
  ErdahlPencilBound<T, Tint> pb = erdahl_pencil_bound<T, Tint>(Q1L, Q2L, os);
  for (auto &v : pb.directions) {
    MyVector<Tint> vn = D3.L.transpose() * v;
    l_cand.push_back(erdahl_direction_row(vn));
  }
  auto is_psd = [&](T const &lambda) -> bool {
    return IsPositiveSemiDefinite(MyMatrix<T>(Q1L - lambda * Q2L), os);
  };
  auto get_lambda = [&]() -> T {
    bool is_first = true;
    T lambda(0);
    for (auto &e : l_cand) {
      T val2 = EvaluationQuadForm<T, Tint>(D2.F, e);
      if (val2 <= 0) {
        continue;
      }
      T ratio = EvaluationQuadForm<T, Tint>(D1.F, e) / val2;
      if (is_first || ratio < lambda) {
        lambda = ratio;
        is_first = false;
      }
    }
    if (is_first) {
      std::cerr << "ERDAHL: no candidate with f2 > 0 in the flip\n";
      throw TerminalException{1};
    }
    return lambda;
  };
  MyMatrix<T> G;
  while (true) {
#ifdef TIMINGS_ERDAHL_FLIP
    n_iter++;
#endif
    T lambda = get_lambda();
    if (!is_psd(lambda)) {
      /*
        lambda > lambda_psd, which is then irrational: otherwise its
        directions are candidates. So the answer is < lambda_psd, and a
        point of ratio below lambda_psd is searched with the functions
        g = f1 - mu f2 for mu < lambda_psd, which are bounded below on the
        cosets of D3. The interval [mu_lo, mu_hi] brackets lambda_psd.
       */
      T mu_lo(0);
      T mu_hi = lambda;
      bool found = false;
      while (!found) {
        T mu = (mu_lo + mu_hi) / 2;
        if (!is_psd(mu)) {
          mu_hi = mu;
          continue;
        }
        MyMatrix<T> Gmu = D1.F - mu * D2.F;
        std::vector<MyVector<Tint>> l_neg =
            erdahl_negative_points_on_set(Gmu, ls3, os);
        if (l_neg.empty()) {
          mu_lo = mu;
        } else {
          for (auto &e : l_neg) {
            l_cand.push_back(e);
          }
          found = true;
        }
      }
#ifdef DEBUG_ERDAHL_FLIP
      os << "ERDAHL: flip, irrational bound handled, new lambda="
         << get_lambda() << "\n";
#endif
      continue;
    }
    G = D1.F - lambda * D2.F;
    std::vector<MyVector<Tint>> l_neg =
        erdahl_negative_points_on_set(G, ls3, os);
#ifdef DEBUG_ERDAHL_FLIP
    os << "ERDAHL: flip, lambda=" << lambda << " |l_cand|=" << l_cand.size()
       << " |l_neg|=" << l_neg.size() << "\n";
#endif
    if (l_neg.empty()) {
      break;
    }
    for (auto &e : l_neg) {
      l_cand.push_back(e);
    }
  }
#ifdef TIMINGS_ERDAHL_FLIP
  os << "|ERDAHL: flip, ratio iterations n_iter=" << n_iter
     << " |l_cand|=" << l_cand.size() << "|=" << time_total << "\n";
#endif
  DelaunayPolyhedron<T, Tint> D2p =
      erdahl_polyhedron_extension<T, Tint>(W, G, D3, compute_function, os);
#ifdef TIMINGS_ERDAHL_FLIP
  os << "|ERDAHL: flip, zero set and canonical function, d3=" << D3.L.rows()
     << " |EXT(D2')|=" << D2p.EXT.rows() << "|=" << time_total << "\n";
#endif
#ifdef SANITY_CHECK_ERDAHL_FLIP
  if (!erdahl_is_subset(D1, D2p) || !erdahl_is_subset(D2p, D3)) {
    std::cerr << "ERDAHL: the flipped polyhedron is not between D1 and D3\n";
    throw TerminalException{1};
  }
  if (erdahl_is_equal(D2, D2p)) {
    std::cerr << "ERDAHL: the flip returned D2\n";
    throw TerminalException{1};
  }
  if (erdahl_perfection_rank(W, D2p) != r2) {
    std::cerr << "ERDAHL: the flipped polyhedron has the wrong rank\n";
    throw TerminalException{1};
  }
#endif
  return D2p;
}

// clang-format off
#endif  // SRC_ERDAHL_ERDAHL_FLIP_H_
// clang-format on
