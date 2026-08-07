// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_PERIODICDELAUNAY_H_
#define SRC_DELAUNAY_PERIODICDELAUNAY_H_

// clang-format off
#include "PeriodicStructures.h"
#include "Shvec_exact.h"
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PERIODIC_DELAUNAY
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_PERIODIC_DELAUNAY
#endif

/*
  Closest vectors of a periodic point set Z^n + {c_1, ..., c_m} to a point.

  The points are expressed by their numerators over the denominator N of
  the point set: the point of numerator u is u / N, and u is congruent
  modulo N to the numerator of the coset it belongs to. Keeping the scaled
  form means the vertex matrices of the periodic Delaunay polytopes are
  integral, exactly like in the lattice case.
 */
template <typename T, typename Tint> struct PeriodicCVPResult {
  T TheNorm;
  MyMatrix<Tint> ListVectScaled;
};

/*
  The closest points of Z^n + c_k to v are c_k plus the closest points of
  Z^n to v - c_k, so the periodic problem is m lattice problems whose
  results are merged on the smallest norm. That merge is the only periodic
  part: the enumeration itself is the existing CVP solver.
 */
template <typename T, typename Tint>
PeriodicCVPResult<T, Tint>
PeriodicClosestVectors(CVPSolver<T, Tint> const &solver,
                       PeriodicPointSet<Tint> const &pps,
                       MyVector<T> const &eV) {
  int n = eV.size();
  int m = pps.cosets_num.rows();
  T N_T = UniversalScalarConversion<T, Tint>(pps.N);
  std::optional<T> min_norm;
  std::vector<MyVector<Tint>> l_vect;
  for (int k = 0; k < m; k++) {
    MyVector<T> eV_shift(n);
    for (int j = 0; j < n; j++) {
      T num = UniversalScalarConversion<T, Tint>(pps.cosets_num(k, j));
      eV_shift(j) = eV(j) - num / N_T;
    }
    resultCVP<T, Tint> res = solver.nearest_vectors(eV_shift);
    if (min_norm && res.TheNorm > *min_norm) {
      continue;
    }
    if (!min_norm || res.TheNorm < *min_norm) {
      min_norm = res.TheNorm;
      l_vect.clear();
    }
    int n_vect = res.ListVect.rows();
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      // The point is c_k + z, of numerator cosets_num(k) + N z.
      MyVector<Tint> u(n);
      for (int j = 0; j < n; j++) {
        u(j) = pps.cosets_num(k, j) + pps.N * res.ListVect(i_vect, j);
      }
      l_vect.push_back(std::move(u));
    }
  }
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  if (!min_norm) {
    std::cerr << "PERIODIC_DELAUNAY: PeriodicClosestVectors: the point set "
                 "has no coset at all\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> ListVectScaled = MatrixFromVectorFamily(l_vect);
  return {*min_norm, std::move(ListVectScaled)};
}

/*
  Whether an affine transformation preserves a periodic point set, and the
  permutation it induces on the cosets when it does.

  With p = c_k + z and the transformation p -> p A + w / d, the image is
  (c_k A + w / d) + z A. It stays in the point set for every z if and only
  if A is unimodular (the lattice Z^n is preserved) and c_k A + w / d is a
  coset modulo Z^n. Multiplied by the denominator N of the point set, that
  last condition is the integral congruence

      cosets_num(k) A + N w / d = cosets_num(sigma(k))   modulo N

  so the whole test is integer arithmetic, provided N w / d is integral --
  which is the invariant of the periodic setting, every transformation of
  the point set having its translation in (1 / N) Z^n.
 */
template <typename Tring>
std::optional<std::vector<size_t>>
PeriodicCosetPermutation(PeriodicPointSet<Tring> const &pps,
                         PeriodicAffineTransform<Tring> const &x) {
  int n = pps.cosets_num.cols();
  int m = pps.cosets_num.rows();
  Tring det = DeterminantMat(x.A);
  if (det != 1 && det != -1) {
    return {};
  }
  // The scaled translation N w / d, integral exactly when the translation
  // lies in (1 / N) Z^n.
  MyVector<Tring> w_scal(n);
  for (int j = 0; j < n; j++) {
    Tring prod = pps.N * x.w(j);
    Tring res = ResInt(prod, x.d);
    if (res != 0) {
      return {};
    }
    w_scal(j) = QuoInt(prod, x.d);
  }
  std::vector<size_t> sigma(m);
  Face f_hit(m);
  for (int k = 0; k < m; k++) {
    MyVector<Tring> img(n);
    for (int j = 0; j < n; j++) {
      Tring eSum = w_scal(j);
      for (int i = 0; i < n; i++) {
        AddMul(eSum, pps.cosets_num(k, i), x.A(i, j));
      }
      img(j) = eSum;
    }
    std::optional<size_t> opt = GetCosetIndex(pps, img);
    if (!opt) {
      return {};
    }
    sigma[k] = *opt;
    f_hit[*opt] = 1;
  }
#ifdef SANITY_CHECK_PERIODIC_DELAUNAY
  // The induced map on the cosets is injective, being the restriction of a
  // bijection of the point set, so all the cosets must be reached.
  if (static_cast<int>(f_hit.count()) != m) {
    std::cerr << "PERIODIC_DELAUNAY: PeriodicCosetPermutation: the induced "
                 "map on the cosets is not a bijection\n";
    throw TerminalException{1};
  }
#endif
  return sigma;
}

template <typename Tring>
bool IsPeriodicPointSetPreserved(PeriodicPointSet<Tring> const &pps,
                                 PeriodicAffineTransform<Tring> const &x) {
  return PeriodicCosetPermutation(pps, x).has_value();
}

/*
  Norm of the difference between the point of numerator u and the point eV,
  for the quadratic form GramMat.
 */
template <typename T, typename Tint>
T PeriodicNormDiff(MyMatrix<T> const &GramMat, Tint const &N,
                   MyVector<Tint> const &u, MyVector<T> const &eV) {
  int n = eV.size();
  T N_T = UniversalScalarConversion<T, Tint>(N);
  MyVector<T> eDiff(n);
  for (int j = 0; j < n; j++) {
    T num = UniversalScalarConversion<T, Tint>(u(j));
    eDiff(j) = num / N_T - eV(j);
  }
  return EvaluationQuadForm<T, T>(GramMat, eDiff);
}

// clang-format off
#endif  // SRC_DELAUNAY_PERIODICDELAUNAY_H_
// clang-format on
