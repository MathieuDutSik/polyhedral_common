// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SHVEC_SIEVE_H_
#define SRC_LATT_SHVEC_SIEVE_H_

// clang-format off
#include "ClassicLLL.h"
#include "LatticeDefinitions.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "QuoIntFcts.h"
#include "Shvec_exact.h"
#include <algorithm>
#include <cmath>
#include <optional>
#include <random>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_SHVEC_SIEVE
#endif

#ifdef DISABLE_DEBUG_SHVEC_SIEVE
#undef DEBUG_SHVEC_SIEVE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_SHVEC_SIEVE
#endif

#ifdef TIMINGS
#define TIMINGS_SHVEC_SIEVE
#endif

/*
  Sieving for short vectors: the Gauss sieve of Micciancio and Voulgaris.

  The enumeration of Shvec_exact.h walks the tree of partial coordinates and
  costs 2^O(n log n). A sieve instead keeps a list L of lattice vectors and
  makes them shorter by subtracting them from one another. In the Gauss sieve
  the list is kept PAIRWISE REDUCED: for any v, w in L neither of v - k w,
  k in Z, is shorter than v, that is

      2 |<v, w>| <= min(|v|^2, |w|^2).

  A new vector, sampled at random or taken from a stack, is first reduced
  against every element of L until no reduction applies; then every element
  of L that the new vector reduces is taken out of L, reduced, and pushed on
  the stack to be processed again; then the new vector enters L. A vector
  reduced to zero is a COLLISION: the sample carried no information that L did
  not already have. The sieve stops once the number of collisions exceeds
  collision_base + collision_ratio * |L|.

  Two pairwise reduced vectors make an angle of at least 60 degrees, so |L| is
  bounded by the kissing-type bound 2^(0.401 n). Heuristically the time is
  2^(0.52 n) and the list ends up containing a shortest vector. Nothing of
  that is a proof: the result of the sieve alone is an UPPER BOUND on the
  minimum, correct as a lattice vector and as a norm, but not certified to be
  the minimum. T_ShortestVectorSieve adds the certification.

  EXACTNESS. The Gram matrix is scaled to an integral one, and every vector
  carries its integral coordinates x, the integral vector G x and the integral
  norm x^T G x. A reduction v <- v - k w, with k the nearest integer to
  <v, w> / |w|^2, updates the three by integer operations, and is performed
  only when 2 |<v, w>| > |w|^2, which makes the norm drop strictly. Every
  decision is therefore taken in exact arithmetic, and strict decrease of an
  integral norm is what rules out a reduction loop.

  THE DOUBLE PRECISION FILTER. Almost every pair tested does not reduce, and
  proving that exactly would cost n multiprecision products per pair. The
  vectors also carry double copies of x and G x, and the inner product is
  first computed in double together with the sum of the absolute values of
  its terms, which bounds its rounding error. The exact test runs only when
  the double one cannot rule the reduction out, so the filter never loses a
  reduction: it is a certified filter, not an approximation. This is where the
  double precision speed of the usual sieves is recovered without giving up
  exactness.

  THE SAMPLER is Klein's randomized nearest plane on the LLL reduced basis, in
  double precision: going from the last Gram-Schmidt direction to the first,
  the coordinate x_i is the rounding of its Babai center plus a Gaussian of
  standard deviation s / |b_i^*|. The sampler only has to produce lattice
  vectors of moderate length spread over the whole lattice; the vectors are
  integral whatever the rounding, so its floating point is harmless.

  References:
  * D. Micciancio, P. Voulgaris, Faster exponential time algorithms for the
    shortest vector problem, SODA 2010, 1468--1480.
  * P. Klein, Finding the closest lattice vector when it's unusually close,
    SODA 2000, 937--941.
  * M. R. Albrecht, L. Ducas, G. Herold, E. Kirshanova, E. W. Postlethwaite,
    M. Stevens, The general sieve kernel and new records in lattice reduction,
    EUROCRYPT 2019, for what a production sieve looks like.
 */

template <typename T> struct SieveOptions {
  // The seed of the sampler, so that a run can be replayed exactly.
  unsigned long seed = 1;
  // The sieve stops after collision_base + collision_ratio * |L| collisions.
  size_t collision_base = 200;
  double collision_ratio = 0.1;
  // Stop as soon as a vector of norm at most target_norm has been found.
  std::optional<T> target_norm;
};

template <typename T, typename Tint> struct SieveResult {
  // The list L at the end of the sieve, by increasing norm, in the
  // coordinates of the input Gram matrix. One vector per antipodal pair.
  std::vector<MyVector<Tint>> list;
  std::vector<T> norms;
  size_t n_samples;
  size_t n_collisions;
};

/*
  A vector of the sieve, in the coordinates of the integral LLL reduced Gram
  matrix: its coordinates, their image by G and the norm, exactly and in
  double.
 */
template <typename Tring> struct SieveVector {
  MyVector<Tring> x;
  MyVector<Tring> Gx;
  Tring norm;
  MyVector<double> x_d;
  MyVector<double> Gx_d;
  double norm_d = 0;
};

template <typename Tring> void SieveRefreshDouble(SieveVector<Tring> &v) {
  int n = v.x.size();
  for (int i = 0; i < n; i++) {
    v.x_d(i) = UniversalScalarConversion<double, Tring>(v.x(i));
    v.Gx_d(i) = UniversalScalarConversion<double, Tring>(v.Gx(i));
  }
  v.norm_d = UniversalScalarConversion<double, Tring>(v.norm);
}

template <typename Tring>
SieveVector<Tring> SieveMakeVector(MyMatrix<Tring> const &G,
                                   MyVector<Tring> const &x) {
  int n = x.size();
  SieveVector<Tring> v{x,
                       MyVector<Tring>(n),
                       Tring(0),
                       MyVector<double>(n),
                       MyVector<double>(n),
                       0.0};
  for (int i = 0; i < n; i++) {
    Tring sum(0);
    for (int j = 0; j < n; j++) {
      sum += G(i, j) * x(j);
    }
    v.Gx(i) = sum;
    v.norm += x(i) * sum;
  }
  SieveRefreshDouble(v);
  return v;
}

/*
  Reduce v by w when that makes v strictly shorter; returns whether it did.
  The double evaluation of <v, w> has an error below n 2^-53 times the sum of
  the absolute values of its terms, and |w|^2 one below 2^-53 |w|^2. The
  margin taken, 1e-12 of those, is several orders of magnitude above both for
  any dimension a sieve can reach, so a skipped pair is one that provably does
  not reduce.
 */
template <typename Tring>
bool SieveReducePair(SieveVector<Tring> &v, SieveVector<Tring> const &w) {
  int n = v.x.size();
  double ip_d = 0;
  double abs_sum = 0;
  for (int i = 0; i < n; i++) {
    double term = v.x_d(i) * w.Gx_d(i);
    ip_d += term;
    abs_sum += std::abs(term);
  }
  double margin = 1.0e-12 * (2 * abs_sum + w.norm_d);
  if (2 * std::abs(ip_d) + margin < w.norm_d) {
    return false;
  }
  Tring ip(0);
  for (int i = 0; i < n; i++) {
    ip += v.x(i) * w.Gx(i);
  }
  Tring ip2 = ip + ip;
  if (ip2 <= w.norm && -ip2 <= w.norm) {
    return false;
  }
  // The nearest integer to ip / |w|^2, as floor((2 ip + |w|^2) / (2 |w|^2)).
  Tring num = ip2 + w.norm;
  Tring den = w.norm + w.norm;
  Tring k = QuoInt(num, den);
#ifdef SANITY_CHECK_SHVEC_SIEVE
  Tring norm_before = v.norm;
#endif
  for (int i = 0; i < n; i++) {
    v.x(i) -= k * w.x(i);
    v.Gx(i) -= k * w.Gx(i);
  }
  v.norm += k * (k * w.norm - ip2);
#ifdef SANITY_CHECK_SHVEC_SIEVE
  if (v.norm >= norm_before || v.norm < 0) {
    std::cerr << "SHVEC_SIEVE: the reduction did not decrease the norm, "
              << norm_before << " -> " << v.norm << "\n";
    throw TerminalException{1};
  }
  Tring norm_check(0);
  for (int i = 0; i < n; i++) {
    norm_check += v.x(i) * v.Gx(i);
  }
  if (norm_check != v.norm) {
    std::cerr << "SHVEC_SIEVE: incremental norm " << v.norm
              << " differs from the recomputed " << norm_check << "\n";
    throw TerminalException{1};
  }
#endif
  SieveRefreshDouble(v);
  return true;
}

/*
  Klein's sampler on the basis of G, in double precision. mu and r are the
  Gram-Schmidt data of G: b_j = b_j^* + sum_{i<j} mu(j,i) b_i^* and
  r(i) = |b_i^*|^2.
 */
struct SieveSampler {
  int n;
  MyMatrix<double> mu;
  std::vector<double> sigma;
  std::mt19937_64 rng;
  std::normal_distribution<double> gauss;

  template <typename Tring>
  SieveSampler(MyMatrix<Tring> const &G, unsigned long seed)
      : n(G.rows()), mu(n, n), sigma(n), rng(seed), gauss(0.0, 1.0) {
    std::vector<double> r(n);
    for (int j = 0; j < n; j++) {
      for (int i = 0; i <= j; i++) {
        double val = UniversalScalarConversion<double, Tring>(G(j, i));
        for (int k = 0; k < i; k++) {
          val -= mu(j, k) * mu(i, k) * r[k];
        }
        if (i < j) {
          mu(j, i) = val / r[i];
        } else {
          r[j] = val;
        }
      }
    }
    // The width s^2 = max_i |b_i^*|^2: every coordinate direction then
    // contributes about s^2 to the norm, which gives vectors a few times
    // longer than the first basis vector, spread over the whole lattice.
    double s2 = 0;
    for (int i = 0; i < n; i++) {
      s2 = std::max(s2, r[i]);
    }
    for (int i = 0; i < n; i++) {
      sigma[i] = std::sqrt(s2 / r[i]);
    }
  }

  template <typename Tring> MyVector<Tring> sample() {
    MyVector<Tring> x(n);
    std::vector<double> xd(n);
    while (true) {
      bool is_zero = true;
      for (int i = n - 1; i >= 0; i--) {
        double center = 0;
        for (int j = i + 1; j < n; j++) {
          center -= mu(j, i) * xd[j];
        }
        xd[i] = std::round(center + sigma[i] * gauss(rng));
        if (xd[i] != 0) {
          is_zero = false;
        }
      }
      if (!is_zero) {
        break;
      }
    }
    for (int i = 0; i < n; i++) {
      x(i) = UniversalScalarConversion<Tring, int64_t>(
          static_cast<int64_t>(xd[i]));
    }
    return x;
  }
};

/*
  The Gauss sieve itself, on an integral Gram matrix whose basis should be
  LLL reduced for the sampler to be effective. f_stop is called on each vector
  entering the list and returns true to end the sieve early. Returns the list
  L, in the coordinates of G.
 */
template <typename Tring, typename Fstop>
std::vector<SieveVector<Tring>>
GaussSieveKernel(MyMatrix<Tring> const &G, size_t const &collision_base,
                 double const &collision_ratio, unsigned long const &seed,
                 Fstop f_stop, size_t &n_samples, size_t &n_collisions,
                 [[maybe_unused]] std::ostream &os) {
  int n = G.rows();
  SieveSampler sampler(G, seed);
  std::vector<SieveVector<Tring>> L;
  std::vector<SieveVector<Tring>> S;
  // The basis vectors come first: they are short already after LLL, and
  // every later sample is reduced against them.
  for (int i = n - 1; i >= 0; i--) {
    MyVector<Tring> e = ZeroVector<Tring>(n);
    e(i) = 1;
    S.emplace_back(SieveMakeVector(G, e));
  }
  n_samples = 0;
  n_collisions = 0;
  auto collision_limit = [&]() -> double {
    return static_cast<double>(collision_base) +
           collision_ratio * static_cast<double>(L.size());
  };
  while (static_cast<double>(n_collisions) < collision_limit()) {
    SieveVector<Tring> v;
    if (S.empty()) {
      v = SieveMakeVector(G, sampler.sample<Tring>());
      n_samples++;
    } else {
      v = std::move(S.back());
      S.pop_back();
    }
    // Reduce v against the list until it is stable. A reduction changes v,
    // after which a vector already passed may reduce it again.
    bool changed = true;
    while (changed && v.norm != 0) {
      changed = false;
      for (auto const &w : L) {
        if (SieveReducePair(v, w)) {
          changed = true;
          if (v.norm == 0) {
            break;
          }
        }
      }
    }
    if (v.norm == 0) {
      n_collisions++;
      continue;
    }
    // Take out of L every vector that v reduces. A vector of L shorter than
    // v cannot be one: 2 |<v, w>| > |v|^2 >= |w|^2 would have let w reduce v.
    size_t pos = 0;
    while (pos < L.size()) {
      if (L[pos].norm > v.norm && SieveReducePair(L[pos], v)) {
        SieveVector<Tring> w = std::move(L[pos]);
        L[pos] = std::move(L.back());
        L.pop_back();
        if (w.norm == 0) {
          n_collisions++;
        } else {
          S.emplace_back(std::move(w));
        }
      } else {
        pos++;
      }
    }
    bool stop = f_stop(v);
    L.emplace_back(std::move(v));
    if (stop) {
      break;
    }
  }
#ifdef DEBUG_SHVEC_SIEVE
  os << "SHVEC_SIEVE: n=" << n << " |L|=" << L.size()
     << " n_samples=" << n_samples << " n_collisions=" << n_collisions << "\n";
#endif
#ifdef SANITY_CHECK_SHVEC_SIEVE
  // The pairwise reducedness that the sieve maintains, checked exactly.
  for (size_t i = 0; i < L.size(); i++) {
    for (size_t j = 0; j < L.size(); j++) {
      if (i == j) {
        continue;
      }
      Tring ip(0);
      for (int u = 0; u < n; u++) {
        ip += L[i].x(u) * L[j].Gx(u);
      }
      Tring ip2 = ip + ip;
      if (ip2 > L[j].norm || -ip2 > L[j].norm) {
        std::cerr << "SHVEC_SIEVE: the list is not pairwise reduced at i=" << i
                  << " j=" << j << "\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  std::sort(L.begin(), L.end(),
            [](SieveVector<Tring> const &a, SieveVector<Tring> const &b) {
              return a.norm < b.norm;
            });
  return L;
}

/*
  The Gauss sieve for a positive definite Gram matrix. The result is a list of
  short vectors of the lattice, by increasing norm; its first element is a
  heuristic shortest vector, whose norm is an upper bound on the minimum.
 */
template <typename T, typename Tint>
SieveResult<T, Tint> GaussSieve(MyMatrix<T> const &GramMat,
                                SieveOptions<T> const &opts, std::ostream &os) {
  using Tring = typename underlying_ring<T>::ring_type;
  [[maybe_unused]] int n = GramMat.rows();
#ifdef TIMINGS_SHVEC_SIEVE
  MicrosecondTime time;
#endif
  LLLreduction<T, Tint> eRec = LLLreducedBasis<T, Tint>(GramMat, os);
#ifdef TIMINGS_SHVEC_SIEVE
  os << "|SHVEC_SIEVE: LLLreducedBasis|=" << time << "\n";
#endif
  // A positive rescaling changes no comparison of norms.
  MyMatrix<Tring> G = RemoveFractionMatrixPlusCoeffRing(eRec.GramMatRed).TheMat;
  auto to_input = [&](MyVector<Tring> const &x) -> MyVector<Tint> {
    MyVector<Tint> x_tint = UniversalVectorConversion<Tint, Tring>(x);
    return eRec.Pmat.transpose() * x_tint;
  };
  // Norms are compared in the scaled integral form; the input norm is only
  // formed for a vector that improves on the best one so far.
  std::optional<Tring> best_norm;
  auto f_stop = [&](SieveVector<Tring> const &v) -> bool {
    if (!opts.target_norm) {
      return false;
    }
    if (best_norm && *best_norm <= v.norm) {
      return false;
    }
    best_norm = v.norm;
    T norm = EvaluationQuadForm<T, Tint>(GramMat, to_input(v.x));
    return norm <= *opts.target_norm;
  };
  size_t n_samples, n_collisions;
  std::vector<SieveVector<Tring>> L =
      GaussSieveKernel<Tring>(G, opts.collision_base, opts.collision_ratio,
                              opts.seed, f_stop, n_samples, n_collisions, os);
#ifdef TIMINGS_SHVEC_SIEVE
  os << "|SHVEC_SIEVE: GaussSieveKernel|=" << time << "\n";
#endif
  std::vector<MyVector<Tint>> list;
  std::vector<T> norms;
  for (auto const &v : L) {
    MyVector<Tint> x = to_input(v.x);
    norms.push_back(EvaluationQuadForm<T, Tint>(GramMat, x));
    list.emplace_back(std::move(x));
  }
#ifdef SANITY_CHECK_SHVEC_SIEVE
  for (size_t i = 1; i < norms.size(); i++) {
    if (norms[i - 1] > norms[i]) {
      std::cerr << "SHVEC_SIEVE: the list is not sorted by norm in the input "
                   "coordinates\n";
      throw TerminalException{1};
    }
  }
  if (n > 0 && list.empty()) {
    std::cerr << "SHVEC_SIEVE: the sieve returned an empty list\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(list), std::move(norms), n_samples, n_collisions};
}

/*
  The shortest vectors, certified. The sieve provides an upper bound B on the
  minimum and the enumeration of Shvec_exact.h then lists every vector of norm
  at most B, from which the minimal ones are kept. The result is exact
  whatever the sieve returned; what the sieve decides is only the cost of the
  enumeration, which is smallest when B is the minimum itself, as it usually
  is. The output has the convention of T_ShortestVector: the rows come in
  pairs v, -v.
 */
template <typename T, typename Tint>
Tshortest<T, Tint> T_ShortestVectorSieve(MyMatrix<T> const &GramMat,
                                         SieveOptions<T> const &opts,
                                         std::ostream &os) {
  int n = GramMat.rows();
  SieveResult<T, Tint> res = GaussSieve<T, Tint>(GramMat, opts, os);
  T const &bound = res.norms[0];
#ifdef TIMINGS_SHVEC_SIEVE
  MicrosecondTime time;
#endif
  CVPSolver<T, Tint> solver(GramMat, os);
  std::vector<MyVector<Tint>> l_vect = solver.at_most_norm_vectors(bound);
#ifdef TIMINGS_SHVEC_SIEVE
  os << "|SHVEC_SIEVE: at_most_norm_vectors|=" << time << "\n";
#endif
  T min = bound;
  std::vector<MyVector<Tint>> l_min;
  for (auto &x : l_vect) {
    T norm = EvaluationQuadForm<T, Tint>(GramMat, x);
    if (norm < min) {
      min = norm;
      l_min.clear();
    }
    if (norm == min) {
      l_min.emplace_back(std::move(x));
    }
  }
#ifdef DEBUG_SHVEC_SIEVE
  os << "SHVEC_SIEVE: sieve bound=" << bound << " minimum=" << min
     << " |l_vect|=" << l_vect.size() << " |l_min|=" << l_min.size() << "\n";
#endif
#ifdef SANITY_CHECK_SHVEC_SIEVE
  if (l_min.empty()) {
    std::cerr << "SHVEC_SIEVE: the enumeration found no vector of norm at "
                 "most the sieve bound "
              << bound << "\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> SHV(2 * l_min.size(), n);
  for (size_t i = 0; i < l_min.size(); i++) {
    for (int j = 0; j < n; j++) {
      SHV(2 * i, j) = l_min[i](j);
      SHV(2 * i + 1, j) = -l_min[i](j);
    }
  }
  return {min, std::move(SHV)};
}

// clang-format off
#endif  // SRC_LATT_SHVEC_SIEVE_H_
// clang-format on
