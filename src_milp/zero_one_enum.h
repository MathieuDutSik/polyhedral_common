// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MILP_ZERO_ONE_ENUM_H_
#define SRC_MILP_ZERO_ONE_ENUM_H_

// clang-format off
#include "MAT_Matrix.h"
#include "Boost_bitset_kernel.h"
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ZERO_ONE_ENUM
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ZERO_ONE_ENUM
#endif

#ifdef TIMINGS
#define TIMINGS_ZERO_ONE_ENUM
#endif

/*
  Enumeration of the lattice vectors all of whose ambient coordinates
  are +-Fq, after the enumeration of A. Wassermann's solvediophant.

  The input is an integral basis B of d vectors of Z^N, and the wanted
  vectors are the v = t B, t in Z^d, with |v_l| = Fq for every l. The
  norm bound ||v||^2 <= N Fq^2 that follows is far weaker than the
  coordinate conditions themselves, and enumerating the ball rather
  than the cube is what makes a plain Fincke-Pohst hopeless here. The
  three prunings below carry the coordinate conditions into the tree.

  Write b*_0, ..., b*_{d-1} for the Gram-Schmidt vectors, c_j for their
  squared norms and mu for the coefficients. The search fixes
  t_{d-1}, ..., t_0 in that order and maintains, at level i,

     y_i = t_i + sum_{k>i} mu_{k,i} t_k
     cs_i = sum_{j>=i} c_j y_j^2        the squared norm of P_i v
     w_i  = sum_{j>=i} y_j b*_j         that same P_i v, in coordinates

  where P_i is the orthogonal projection on the span of b*_i, ..., that
  is on the orthogonal complement of b_0, ..., b_{i-1}.

  --- Norm. cs_i <= ||v||^2 <= Fd, the usual Fincke-Pohst bound, and it
      is what bounds the range of t_i.
  --- Hoelder. P_i being an orthogonal projection,
         cs_i = ||P_i v||^2 = <P_i v, v> <= ||v||_infinity ||P_i v||_1
      so cs_i <= Fq ||w_i||_1 holds for every prefix of a solution. It
      is a much sharper test than the norm one, since it is the
      coordinate bound and not the norm that it uses.
  --- Determined coordinates. Let first_nonzero[l] be the least j with
      B[j][l] nonzero. For i = first_nonzero[l] every vector of the
      span of b_0, ..., b_{i-1} has its l-th coordinate zero, so
      (P_i v)[l] = v[l]: the coordinate l is settled at level i, and
      |w_i[l]| must be exactly Fq. Half of the coordinates of the
      instances at hand are settled by the middle of the tree, so this
      is the pruning that does the work.
  --- Dual bounds. With d_i the dual basis of the b_j inside their
      span, t_i = <v, d_i>, whence
         |t_i| <= min(Fq ||d_i||_1, sqrt(Fd) ||d_i||_2)
      by Hoelder and by Cauchy-Schwarz. Cheap, and it caps the range of
      every coefficient before the search starts.

  The Gram-Schmidt and the tests are in double precision, as in
  solvediophant, with a tolerance eps: the prunings are relaxed by eps
  so that a solution is never cut by a rounding error, and every vector
  reaching the bottom is verified exactly in integers before being
  reported. So a reported solution is certain, while the completeness
  of the enumeration holds up to that tolerance.
*/

struct ZeroOneEnumOptions {
  // Absolute value that every ambient coordinate of a solution has
  double Fq = 1.0;
  // Tolerance of the floating point tests
  double eps = 1e-6;
  size_t max_node = std::numeric_limits<size_t>::max();
  bool use_hoelder = true;
  bool use_only_zeros = true;
  bool use_dual_bounds = true;
};

struct ZeroOneEnumResult {
  bool resolved;
  size_t n_node;
  size_t n_solution;
  size_t n_prune_hoelder;
  size_t n_prune_only_zeros;
  size_t n_prune_dual;
};

template <typename Tint, typename F> struct ZeroOneEnumSearch {
private:
  int d;
  int n_amb;
  MyMatrix<Tint> const &B;
  ZeroOneEnumOptions opt;
  F f;
  std::ostream &os;
  // Gram-Schmidt in double
  std::vector<std::vector<double>> Bd;
  std::vector<std::vector<double>> bstar;
  std::vector<std::vector<double>> mu;
  std::vector<double> cnorm;
  // Coordinates settled when the search reaches a given level
  std::vector<std::vector<int>> settled_at;
  std::vector<double> dual_bound;
  // Per level state
  std::vector<std::vector<double>> w;
  std::vector<double> cs;
  std::vector<double> sigma;
  std::vector<long> t;
  double Fd;
  size_t n_node;
  size_t n_solution;
  size_t n_prune_hoelder;
  size_t n_prune_only_zeros;
  size_t n_prune_dual;
  bool budget_hit;
  bool has_dead_coordinate;

  void GramSchmidt() {
    bstar.assign(d, std::vector<double>(n_amb, 0.0));
    mu.assign(d, std::vector<double>(d, 0.0));
    cnorm.assign(d, 0.0);
    for (int i = 0; i < d; i++) {
      for (int l = 0; l < n_amb; l++)
        bstar[i][l] = Bd[i][l];
      for (int j = 0; j < i; j++) {
        double scal = 0.0;
        for (int l = 0; l < n_amb; l++)
          scal += Bd[i][l] * bstar[j][l];
        mu[i][j] = scal / cnorm[j];
        for (int l = 0; l < n_amb; l++)
          bstar[i][l] -= mu[i][j] * bstar[j][l];
      }
      double norm = 0.0;
      for (int l = 0; l < n_amb; l++)
        norm += bstar[i][l] * bstar[i][l];
      cnorm[i] = norm;
    }
  }

  // first_nonzero[l], and the coordinates settled at each level
  void ComputeSettled() {
    settled_at.assign(d, std::vector<int>());
    has_dead_coordinate = false;
    for (int l = 0; l < n_amb; l++) {
      int first = -1;
      for (int i = 0; i < d; i++) {
        if (B(i, l) != 0) {
          first = i;
          break;
        }
      }
      if (first == -1) {
        // The coordinate is zero for every lattice vector, so it can
        // never be equal to Fq and there is no solution at all
        has_dead_coordinate = true;
        return;
      }
      settled_at[first].push_back(l);
    }
  }

  // The dual basis inside the span, from the inverse of the unit lower
  // triangular matrix of the mu
  void ComputeDualBounds() {
    dual_bound.assign(d, std::numeric_limits<double>::max());
    if (!opt.use_dual_bounds)
      return;
    std::vector<std::vector<double>> minv(d, std::vector<double>(d, 0.0));
    for (int i = 0; i < d; i++) {
      minv[i][i] = 1.0;
      for (int j = i - 1; j >= 0; j--) {
        double scal = 0.0;
        for (int k = j + 1; k <= i; k++)
          scal += mu[k][j] * minv[i][k];
        minv[i][j] = -scal;
      }
    }
    double sqrt_Fd = std::sqrt(Fd);
    std::vector<double> dvec(n_amb);
    for (int i = 0; i < d; i++) {
      for (int l = 0; l < n_amb; l++)
        dvec[l] = 0.0;
      for (int j = i; j < d; j++) {
        double coef = minv[j][i] / cnorm[j];
        if (coef == 0.0)
          continue;
        for (int l = 0; l < n_amb; l++)
          dvec[l] += coef * bstar[j][l];
      }
      double l1 = 0.0, l2 = 0.0;
      for (int l = 0; l < n_amb; l++) {
        l1 += std::fabs(dvec[l]);
        l2 += dvec[l] * dvec[l];
      }
      l2 = std::sqrt(l2);
      dual_bound[i] = std::min(opt.Fq * l1, sqrt_Fd * l2) * (1.0 + opt.eps);
    }
  }

  // All the coordinates are fixed: check in exact integers
  void TreatLeaf() {
    MyVector<Tint> v = ZeroVector<Tint>(n_amb);
    for (int i = 0; i < d; i++) {
      if (t[i] == 0)
        continue;
      Tint ti = UniversalScalarConversion<Tint, long>(t[i]);
      for (int l = 0; l < n_amb; l++)
        v(l) += ti * B(i, l);
    }
    Tint Fq_i = UniversalScalarConversion<Tint, long>(
        static_cast<long>(std::llround(opt.Fq)));
    for (int l = 0; l < n_amb; l++) {
      if (v(l) != Fq_i && v(l) != -Fq_i)
        return;
    }
    n_solution++;
    f(v);
  }

  void Explore(int level) {
    if (budget_hit)
      return;
    if (level < 0) {
      TreatLeaf();
      return;
    }
    // sigma_level = sum_{k>level} mu[k][level] t_k
    double sig = 0.0;
    for (int k = level + 1; k < d; k++)
      if (t[k] != 0)
        sig += mu[k][level] * static_cast<double>(t[k]);
    sigma[level] = sig;
    // Every solution has squared norm exactly Fd, so cs lands on the
    // bound and rounding can push it barely over: the same tolerance as
    // the pruning test is needed here, or solutions are lost.
    double remain = Fd * (1.0 + opt.eps) - cs[level + 1];
    if (remain < 0.0)
      return;
    double radius = std::sqrt(remain / cnorm[level]) * (1.0 + opt.eps);
    double centre = -sig;
    long lo = static_cast<long>(std::ceil(centre - radius - opt.eps));
    long hi = static_cast<long>(std::floor(centre + radius + opt.eps));
    if (opt.use_dual_bounds) {
      long db = static_cast<long>(std::floor(dual_bound[level]));
      if (lo < -db) {
        n_prune_dual++;
        lo = -db;
      }
      if (hi > db) {
        n_prune_dual++;
        hi = db;
      }
    }
    // Enumerate outwards from the centre, so that the likely values
    // come first. The number of steps has to reach the farther of the
    // two ends, the centre having no reason to sit in the middle.
    long centre_i = static_cast<long>(std::llround(centre));
    if (centre_i < lo)
      centre_i = lo;
    if (centre_i > hi)
      centre_i = hi;
    long span = std::max(hi - centre_i, centre_i - lo);
    for (long step = 0; step <= 2 * span; step++) {
      long cand =
          (step % 2 == 0) ? centre_i + step / 2 : centre_i - (step + 1) / 2;
      if (cand < lo || cand > hi)
        continue;
      n_node++;
      if (n_node > opt.max_node) {
        budget_hit = true;
        return;
      }
      t[level] = cand;
      double y = static_cast<double>(cand) + sig;
      double cs_new = cs[level + 1] + cnorm[level] * y * y;
      if (cs_new > Fd * (1.0 + opt.eps))
        continue;
      // w_level = w_{level+1} + y b*_level, and its 1-norm
      double norm1 = 0.0;
      std::vector<double> &wl = w[level];
      std::vector<double> const &wp = w[level + 1];
      for (int l = 0; l < n_amb; l++) {
        wl[l] = wp[l] + y * bstar[level][l];
        norm1 += std::fabs(wl[l]);
      }
      if (opt.use_hoelder) {
        if (cs_new > opt.Fq * norm1 * (1.0 + opt.eps)) {
          n_prune_hoelder++;
          continue;
        }
      }
      if (opt.use_only_zeros) {
        bool bad = false;
        for (auto &l : settled_at[level]) {
          if (std::fabs(std::fabs(wl[l]) - opt.Fq) > opt.eps) {
            bad = true;
            break;
          }
        }
        if (bad) {
          n_prune_only_zeros++;
          continue;
        }
      }
      cs[level] = cs_new;
      Explore(level - 1);
      if (budget_hit)
        return;
    }
    t[level] = 0;
  }

public:
  ZeroOneEnumSearch(MyMatrix<Tint> const &_B, ZeroOneEnumOptions const &_opt,
                    F _f, std::ostream &_os)
      : d(_B.rows()), n_amb(_B.cols()), B(_B), opt(_opt), f(_f), os(_os),
        n_node(0), n_solution(0), n_prune_hoelder(0), n_prune_only_zeros(0),
        n_prune_dual(0), budget_hit(false), has_dead_coordinate(false) {
    Fd = static_cast<double>(n_amb) * opt.Fq * opt.Fq;
    Bd.assign(d, std::vector<double>(n_amb, 0.0));
    for (int i = 0; i < d; i++)
      for (int l = 0; l < n_amb; l++)
        Bd[i][l] = UniversalScalarConversion<double, Tint>(B(i, l));
    GramSchmidt();
    ComputeSettled();
    ComputeDualBounds();
    w.assign(d + 1, std::vector<double>(n_amb, 0.0));
    cs.assign(d + 1, 0.0);
    sigma.assign(d + 1, 0.0);
    t.assign(d + 1, 0);
  }

  ZeroOneEnumResult Run() {
    if (d == 0 || has_dead_coordinate) {
#ifdef DEBUG_ZERO_ONE_ENUM
      os << "ZERO_ONE_ENUM: nothing to enumerate, d=" << d
         << " dead coordinate=" << has_dead_coordinate << "\n";
#endif
      return {true, 0, 0, 0, 0, 0};
    }
#ifdef DEBUG_ZERO_ONE_ENUM
    size_t n_settled_end = settled_at[0].size();
    os << "ZERO_ONE_ENUM: d=" << d << " n_amb=" << n_amb << " Fd=" << Fd
       << " coordinates settled only at the last level=" << n_settled_end
       << "\n";
#endif
    Explore(d - 1);
#ifdef DEBUG_ZERO_ONE_ENUM
    os << "ZERO_ONE_ENUM: n_node=" << n_node << " n_solution=" << n_solution
       << " prune_hoelder=" << n_prune_hoelder
       << " prune_only_zeros=" << n_prune_only_zeros << "\n";
#endif
    return {!budget_hit,      n_node,          n_solution,
            n_prune_hoelder,  n_prune_only_zeros, n_prune_dual};
  }
};

// Enumerate the v = t B, t in Z^d, all of whose coordinates are +-Fq,
// calling f on each of them.
template <typename Tint, typename F>
ZeroOneEnumResult EnumeratePlusMinusVectors(MyMatrix<Tint> const &B,
                                            ZeroOneEnumOptions const &opt, F f,
                                            std::ostream &os) {
#ifdef TIMINGS_ZERO_ONE_ENUM
  MicrosecondTime time;
#endif
  ZeroOneEnumSearch<Tint, F> search(B, opt, f, os);
  ZeroOneEnumResult result = search.Run();
#ifdef TIMINGS_ZERO_ONE_ENUM
  os << "ZERO_ONE_ENUM: EnumeratePlusMinusVectors took " << time << "\n";
#endif
  return result;
}

// clang-format off
#endif  // SRC_MILP_ZERO_ONE_ENUM_H_
// clang-format on
