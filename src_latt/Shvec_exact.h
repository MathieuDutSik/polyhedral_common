// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SHVEC_EXACT_H_
#define SRC_LATT_SHVEC_EXACT_H_

// clang-format off
#include "LatticeDefinitions.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "POLY_LinearProgramming.h"
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_SHVEC
#define SANITY_CHECK_SHVEC_EXACT_POLYTOPE
#endif

#ifdef DEBUG
#define DEBUG_SHVEC
#define DEBUG_SHVEC_EXACT_POLYTOPE
// #define DEBUG_SHVEC_VECTOR
// #define DEBUG_SHVEC_MATRIX
#endif

#ifdef DISABLE_DEBUG_SHVEC
#undef DEBUG_SHVEC
#endif

#ifdef DISABLE_DEBUG_SHVEC_EXACT_POLYTOPE
#undef DEBUG_SHVEC_EXACT_POLYTOPE
#endif

#ifdef TIMINGS
#define TIMINGS_SHVEC
#endif

/*
  The native-integer enumeration (int32_t / int64_t when a priori bounds
  certify that every intermediate fits) is OPT-IN: define ENABLE_SHVEC_FAST_PATH
  to use it. It makes the enumeration itself about twice as fast, but each
  computeIt call first pays a preparation whose cost scales with the matrix
  rather than with the number of vectors enumerated: the coset denominator
  scan, the master bound, and the conversion of the Gram matrix and of its
  Bareiss data to the native type. On the small shells of the iso-Delaunay
  work (dimension <= 10, 10 to 500 vectors) that preparation costs more than
  the enumeration saves: measured +1.1% to +4.5% on the five T-space cases,
  worst on the smallest shells. The path is correct and fully tested; it is
  the size regime that does not suit it. See TODO.md.
 */
#ifdef ENABLE_SHVEC_FAST_PATH
#define SHVEC_FAST_PATH
#endif

/*
  The Gram matrix together with its Bareiss decomposition: d(i) the leading
  principal minors and Nmat(i,j) the elimination entries, related to the
  classical Cholesky coefficients by q(i,i) = d(i)/d(i-1) and
  q(i,j) = Nmat(i,j)/d(i). The enumeration reads them but never changes
  them, so they are a property of the matrix rather than of the call: a
  solver enumerating many shells of one matrix computes them once.
 */
template <typename T> struct FullGramInfo {
  int dim;
  MyMatrix<T> gram_matrix;
  MyVector<T> d;
  MyMatrix<T> Nmat;
  FullGramInfo() : dim(0) {}
  /*
    Fraction-free replacement of the Cholesky-style normalization. Keeping d
    and Nmat instead of q removes every quotient from the enumeration: the
    only division here is the Bareiss one, exact over any integral domain.
    The decomposition is a constructor rather than a free function so that
    an instance without it cannot be built by mistake.
   */
  FullGramInfo(MyMatrix<T> const &gram)
      : dim(gram.rows()), gram_matrix(gram), d(gram.rows()), Nmat(gram) {
    T prev(1);
    for (int i = 0; i < dim; i++) {
      d(i) = Nmat(i, i);
#ifdef SANITY_CHECK_SHVEC
      if (d(i) <= 0) {
        std::cerr << "SHVEC: leading principal minor " << i << " is " << d(i)
                  << ", the Gram matrix is not positive definite\n";
        throw TerminalException{1};
      }
#endif
      for (int i2 = i + 1; i2 < dim; i2++) {
        for (int j2 = i + 1; j2 < dim; j2++) {
          T val = Nmat(i2, j2) * Nmat(i, i) - Nmat(i2, i) * Nmat(i, j2);
          // exact division guaranteed by Bareiss's theorem
          T quot = val / prev;
#ifdef SANITY_CHECK_SHVEC
          if (quot * prev != val) {
            std::cerr << "SHVEC: non-exact Bareiss division, T is not an "
                         "integral domain\n";
            throw TerminalException{1};
          }
#endif
          Nmat(i2, j2) = quot;
        }
      }
      prev = Nmat(i, i);
    }
  }
  // For the native-integer path, whose decomposition is converted from an
  // exact one already computed rather than recomputed.
  FullGramInfo(MyMatrix<T> _gram, MyVector<T> _d, MyMatrix<T> _Nmat)
      : dim(_gram.rows()), gram_matrix(std::move(_gram)), d(std::move(_d)),
        Nmat(std::move(_Nmat)) {}
};

template <typename T, typename Tint> struct T_shvec_info {
  std::vector<MyVector<Tint>> short_vectors;
  T minimum;
};

template <typename T, typename Tint>
std::pair<MyVector<Tint>, MyVector<T>>
ReductionMod1vector(MyVector<T> const &V) {
  int n = V.size();
  MyVector<T> v_T(n);
  MyVector<Tint> v_Tint(n);
  for (int i = 0; i < n; i++) {
    T val = V(i);
    Tint red_i = UniversalNearestScalarInteger<Tint, T>(val);
    T red_T = UniversalScalarConversion<T, Tint>(red_i);
    T val_red = val - red_T;
    v_T(i) = val_red;
    v_Tint(i) = red_i;
  }
  return {std::move(v_Tint), std::move(v_T)};
}

/*
  Whether the native-integer fast path can be attempted for the type: the
  preparation needs the denominators, the conversions to and from int64_t
  and the exact square root. The detection is by capability, so a type
  lacking one of them (an algebraic field for example) silently keeps the
  exact path.
 */
template <typename T, typename = void>
struct shvec_has_fast_path : std::false_type {};
#ifdef SHVEC_FAST_PATH
template <typename T>
struct shvec_has_fast_path<
    T, std::void_t<
           decltype(TYPE_CONVERSION(std::declval<stc<T> const &>(),
                                    std::declval<int64_t &>())),
           decltype(TYPE_CONVERSION(std::declval<stc<int64_t> const &>(),
                                    std::declval<T &>())),
           decltype(GetDenominator(std::declval<T const &>())),
           decltype(GcdPair(std::declval<T const &>(),
                            std::declval<T const &>()))>> : std::true_type {};
#endif

/*
  The two bounds of the enumeration in fraction-free form. The admissible
  range of a coordinate n is |den * n + B| <= sqrt(K) with den > 0 and
  K >= 0: no quotient is formed, only products and comparisons.

  Bound_Floor returns the largest n with den * n + B <= sqrt(K). The double
  computation only seeds the search, the exact predicate then walks to the
  right answer.
 */
template <typename T, typename Tint>
Tint Bound_Floor(T const &K, T const &B, T const &den) {
#ifdef SANITY_CHECK_SHVEC
  if (K < 0 || den <= 0) {
    std::cerr << "SHVEC: Bound_Floor, K=" << K << " den=" << den
              << " but K >= 0 and den > 0 are required\n";
    throw TerminalException{1};
  }
#endif
  double K_doubl = UniversalScalarConversion<double, T>(K);
  double B_doubl = UniversalScalarConversion<double, T>(B);
  double den_doubl = UniversalScalarConversion<double, T>(den);
  double alpha = (sqrt(K_doubl) - B_doubl) / den_doubl;
  Tint eReturn = static_cast<Tint>(lround(floor(alpha)));
  auto f = [&](Tint const &n) -> bool {
    T val = den * UniversalScalarConversion<T, Tint>(n) + B;
    if (val <= 0) {
      return true;
    }
    return val * val <= K;
  };
  bool test1 = f(eReturn);
  bool test2 = f(eReturn + 1);
  while (true) {
    if (test1 && !test2) {
      break;
    }
    if (!test1) {
      test2 = test1;
      test1 = f(eReturn - 1);
      eReturn--;
    } else {
      test1 = test2;
      test2 = f(eReturn + 2);
      eReturn++;
    }
  }
  return eReturn;
}

// Smallest n with den * n + B >= -sqrt(K). Same contract as Bound_Floor.
template <typename T, typename Tint>
Tint Bound_Ceil(T const &K, T const &B, T const &den) {
#ifdef SANITY_CHECK_SHVEC
  if (K < 0 || den <= 0) {
    std::cerr << "SHVEC: Bound_Ceil, K=" << K << " den=" << den
              << " but K >= 0 and den > 0 are required\n";
    throw TerminalException{1};
  }
#endif
  double K_doubl = UniversalScalarConversion<double, T>(K);
  double B_doubl = UniversalScalarConversion<double, T>(B);
  double den_doubl = UniversalScalarConversion<double, T>(den);
  double alpha = (-sqrt(K_doubl) - B_doubl) / den_doubl;
  Tint eReturn = static_cast<Tint>(lround(ceil(alpha)));
  auto f = [&](Tint const &n) -> bool {
    T val = den * UniversalScalarConversion<T, Tint>(n) + B;
    if (val >= 0) {
      return true;
    }
    return val * val <= K;
  };
  bool test1 = f(eReturn - 1);
  bool test2 = f(eReturn);
  while (true) {
    if (!test1 && test2) {
      break;
    }
    if (test1) {
      test2 = test1;
      test1 = f(eReturn - 2);
      eReturn--;
    } else {
      test1 = test2;
      test2 = f(eReturn + 1);
      eReturn++;
    }
  }
  return eReturn;
}

/*
  The coset is passed as integral numerators Cnum with the common denominator
  cden: the true coset is Cnum / cden and the bound argument is the scaled
  cden^2 * bound. All state below is scaled accordingly (x enters through
  x_S = cden * x), which keeps the whole enumeration free of fractions even
  for a fractional coset. With cden = 1 the recursion is literally the
  unscaled one. The norm handed to f_insert is the scaled cden^2 * norm.
 */
template <typename T, typename Tint, typename Finsert, typename Fsetbound>
bool computeIt_Gen_Kernel(const FullGramInfo<T> &request,
                          MyVector<T> const &Cnum, T const &cden,
                          bool const &central, const T &bound, Finsert f_insert,
                          Fsetbound f_set_bound) {
  // T may be a ring: the enumeration is fraction-free, the only divisions
  // being the exact Bareiss ones.
  int i, j;
  int dim = request.dim;
  MyVector<Tint> Upper(dim);
  // Fraction-free state, see the derivation above the Bareiss loop below.
  // Ttil(i) = d(i) * (the remaining scaled norm budget at level i) and
  // Bvec(i) = d(i) * Cnum(i) + sum_{j>i} N(i,j) * (x_S(j) + Cnum(j)).
  MyVector<T> Ttil(dim);
  MyVector<T> Bvec(dim);
  MyVector<Tint> x(dim);
  auto is_x_zero = [&]() -> bool {
    for (int u = 0; u < dim; u++) {
      if (x(u) != 0) {
        return false;
      }
    }
    return true;
  };
  MyVector<T> x_S(dim);
#if defined SANITY_CHECK_SHVEC || defined DEBUG_SHVEC
  const MyMatrix<T> &g = request.gram_matrix;
#endif
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: g=\n";
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++)
      std::cerr << " " << g(i, j);
    std::cerr << "\n";
  }
#endif
  // The Bareiss decomposition travels with the matrix, see FullGramInfo.
  MyVector<T> const &d = request.d;
  MyMatrix<T> const &N = request.Nmat;
#ifdef DEBUG_SHVEC_MATRIX
  for (i = 0; i < dim; i++) {
    std::cerr << "SHVEC: d=" << d(i) << "\n";
    for (int j = i + 1; j < dim; j++)
      std::cerr << "   j=" << j << " N=" << N(i, j) << "\n";
  }
#endif
  // The coefficient of x(i) in the scaled offset d(i) * x_S(i) + Bvec(i),
  // handed to f_set_bound.
  MyVector<T> den(dim);
  for (i = 0; i < dim; i++) {
    den(i) = d(i) * cden;
  }
  bool needs_new_bound = true;
  i = dim - 1;
  if (bound < 0) {
    return true;
  }
  Ttil(i) = d(i) * bound;
  Bvec(i) = d(i) * Cnum(i);
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: Before while loop\n";
#endif
#ifdef DEBUG_SHVEC_VECTOR
  size_t n_vector = 0;
#endif
  // Declared once: a fresh T inside the loop means an allocation and a
  // release per iteration for the multiprecision instantiations.
  T Kval, Hval, eNorm, num, quot;
  while (true) {
    if (needs_new_bound) {
      // The admissible range is |d(i) * x(i) + Bvec(i)| <= sqrt(Kval) with
      // Kval = d(i-1) * Ttil(i), which is the fraction-free form of the
      // classical q(i,i) * (x(i) + C(i) + U(i))^2 <= Trem(i).
      Kval = Ttil(i);
      if (i > 0) {
        Kval *= d(i - 1);
      }
      f_set_bound(Kval, Bvec(i), den(i), x, i, Upper(i), x(i));
      x_S(i) = cden * UniversalScalarConversion<T, Tint>(x(i));
      needs_new_bound = false;
    } else {
      x(i) += 1;
      x_S(i) += cden;
    }
    if (x(i) <= Upper(i)) {
      if (i == 0) {
        if (central) {
          if (is_x_zero()) {
#ifdef DEBUG_SHVEC
            std::cerr << "SHVEC: Exiting because x=0 and central run\n";
#endif
            return true;
          }
        }
        Hval = d(0) * x_S(0) + Bvec(0);
        // bound - Trem(-1), the exact norm; the division is exact because
        // Trem(-1) = bound - Q(x + C) is a value of the original form.
        num = Ttil(0) - Hval * Hval;
        quot = num / d(0);
#ifdef SANITY_CHECK_SHVEC
        if (quot * d(0) != num) {
          std::cerr << "SHVEC: non-exact division in the norm evaluation\n";
          throw TerminalException{1};
        }
#endif
        eNorm = bound - quot;
#ifdef SANITY_CHECK_SHVEC
        T norm(0);
        for (int i2 = 0; i2 < dim; i2++)
          for (int j2 = 0; j2 < dim; j2++)
            norm += g(i2, j2) * (x_S(i2) + Cnum(i2)) * (x_S(j2) + Cnum(j2));
        if (norm != eNorm) {
          std::cerr << "Norm inconsistency\n";
          std::cerr << "norm=" << norm << "\n";
          std::cerr << "eNorm=" << eNorm << "\n";
          throw TerminalException{1};
        }
        if (eNorm > bound) {
          std::cerr << "eNorm is too large\n";
          T eDiff = eNorm - bound;
          double bound_doubl = UniversalScalarConversion<double, T>(bound);
          double eNorm_doubl = UniversalScalarConversion<double, T>(eNorm);
          double eDiff_doubl = UniversalScalarConversion<double, T>(eDiff);
          std::cerr << "bound_doubl=" << bound_doubl << "\n";
          std::cerr << "eNorm_doubl=" << eNorm_doubl << "\n";
          std::cerr << "eDiff_doubl=" << eDiff_doubl << "\n";
          std::cerr << "bound=" << bound << "\n";
          std::cerr << "eNorm=" << eNorm << "\n";
          throw TerminalException{1};
        }
#endif
#ifdef DEBUG_SHVEC_VECTOR
        std::cerr << "SHVEC: n_vector=" << n_vector;
        std::cerr << " x=";
        for (int i = 0; i < dim; i++)
          std::cerr << " " << x(i);
        std::cerr << "\n";
        n_vector++;
#endif
        bool ret_val = f_insert(x, eNorm);
        if (!ret_val) {
          return false;
        }
      } else {
        Hval = d(i) * x_S(i) + Bvec(i);
        i--;
        Bvec(i) = d(i) * Cnum(i);
        for (j = i + 1; j < dim; j++)
          Bvec(i) += N(i, j) * (x_S(j) + Cnum(j));
        // Ttil(i) = (d(i) * Ttil(i+1) - Hval^2) / d(i+1), the fraction-free
        // form of Trem(i) = Trem(i+1) - q(i+1,i+1) * hVal^2. The division is
        // exact: d(i) * Trem(i) is a value of the Schur complement of the
        // leading block, which Bareiss keeps over the ring.
        num = d(i) * Ttil(i + 1) - Hval * Hval;
        quot = num / d(i + 1);
#ifdef SANITY_CHECK_SHVEC
        if (quot * d(i + 1) != num) {
          std::cerr << "SHVEC: non-exact division in the remainder update\n";
          throw TerminalException{1};
        }
#endif
        Ttil(i) = quot;
        needs_new_bound = true;
      }
    } else {
      i++;
      if (i == dim) {
        return true;
      }
    }
  }
}

/*
  The enumeration runs in T itself, ring or field alike: computeIt_Gen_Kernel
  forms no quotient, so a ring input no longer has to be lifted to the
  overlying field.
 */
template <typename T, typename Tint, typename Finsert, typename Fsetbound>
inline bool computeIt_Gen(const FullGramInfo<T> &request,
                          MyVector<T> const &coset, bool const &central,
                          const T &bound, Finsert f_insert,
                          Fsetbound f_set_bound) {
  T cden(1);
  return computeIt_Gen_Kernel<T, Tint, Finsert, Fsetbound>(
      request, coset, cden, central, bound, f_insert, f_set_bound);
}

template <typename T, typename Tint, typename Finsert>
bool computeIt_polytope(const FullGramInfo<T> &request,
                        MyVector<T> const &coset, bool const &central,
                        const T &bound, const MyMatrix<T> &FAC,
                        Finsert f_insert, std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n_rows = FAC.rows();
  int n_col = FAC.cols();
#ifdef SANITY_CHECK_SHVEC_EXACT_POLYTOPE
  if (n_col != request.dim + 1) {
    std::cerr << "SHVEC: request.dim=" << request.dim
              << " |request.gram_matrix|=" << request.gram_matrix.rows()
              << "\n";
    std::cerr << "SHVEC: n_col=" << n_col
              << " request.dim + 1=" << (request.dim + 1) << "\n";
    std::cerr << "SHVEC: Error in the size of FAC\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_SHVEC_EXACT_POLYTOPE
  std::cerr << "Beginning of computeIt_polytope\n";
#endif
  auto f_set_bound = [&](const T &K, const T &B, const T &den,
                         const MyVector<Tint> &x, const int &i, Tint &upper,
                         Tint &lower) -> void {
    upper = Bound_Floor<T, Tint>(K, B, den);
    lower = Bound_Ceil<T, Tint>(K, B, den);
    int len = 2 + i;
    MyMatrix<T> FACwork(n_rows, len);
    for (int i_row = 0; i_row < n_rows; i_row++) {
      for (int i_col = 0; i_col < len; i_col++) {
        FACwork(i_row, i_col) = FAC(i_row, i_col);
      }
      for (int i_col = len; i_col < n_col; i_col++) {
        FACwork(i_row, 0) += x(i_col - 1) * FAC(i_row, i_col);
      }
    }
    LpSolution<T> eSol;
    MyVector<T> Vminimize = ZeroVector<T>(len);
    //
    Vminimize(1 + i) = 1;
    eSol = SIMPLEX_LinearProgramming(FACwork, Vminimize, os);
    if (eSol.DualSolution && eSol.DirectSolution) {
      // Well defined so we get a potential lower bound
      Tint eLow = UniversalCeilScalarInteger<Tint, T>(eSol.OptimalValue);
      if (eLow > lower)
        lower = eLow;
    }
    if (!eSol.DualSolution && eSol.DirectSolution) {
      // Infinite direction. Therefore no better bound available
    }
    if (!eSol.DirectSolution) {
      // No feasible solution. Therefore not feasible.
      // This will lead to a backtrack operation
      upper = lower - 1;
      return;
    }
    //
    Vminimize(1 + i) = -1;
    eSol = SIMPLEX_LinearProgramming(FACwork, Vminimize, os);
    if (eSol.DualSolution && eSol.DirectSolution) {
      // Well defined so we get a potential upper bound
      Tint eUpp = UniversalFloorScalarInteger<Tint, T>(-eSol.OptimalValue);
      if (eUpp < upper)
        upper = eUpp;
    }
    if (!eSol.DualSolution &&
        eSol.DirectSolution) { // Infinite direction. Therefore no bound
                               // available
    }
    if (!eSol.DirectSolution) { // No feasible solution. Therefore not feasible.
      upper = lower - 1;       // This will lead to a backtrack operation
      return;
    }
  };
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(
      request, coset, central, bound, f_insert, f_set_bound);
}

/*
  The part of the fast-path preparation that depends only on the Gram
  matrix: the denominator clearing, the Bareiss data and the largest
  intermediate its computation forms. A CVPSolver computes it once and
  reuses it across all its enumerations, since the preparation would
  otherwise dominate the many small per-shell calls.
 */
template <typename T> struct ShvecFastPrep {
  bool usable;
  T D;
  MyMatrix<T> Gscal;
  MyVector<T> d;
  MyMatrix<T> Nmat;
  T maxBareiss;
};

template <typename T>
ShvecFastPrep<T> ComputeShvecFastPrep(MyMatrix<T> const &gram) {
  int dim = gram.rows();
  T D(1);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      D = LCMpair(D, GetDenominator(gram(i, j)));
    }
  }
  MyMatrix<T> Gscal(dim, dim);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      Gscal(i, j) = D * gram(i, j);
    }
  }
  T M(0);
  auto upd = [&](T const &v) -> void {
    T a = T_abs(v);
    if (a > M) {
      M = a;
    }
  };
  // Bareiss with tracking of every intermediate product.
  MyMatrix<T> N = Gscal;
  MyVector<T> d(dim);
  T prev(1);
  for (int i = 0; i < dim; i++) {
    d(i) = N(i, i);
    if (d(i) <= 0) {
      return {false, std::move(D), std::move(Gscal), std::move(d),
              std::move(N), std::move(M)};
    }
    upd(d(i));
    for (int i2 = i + 1; i2 < dim; i2++) {
      for (int j2 = i + 1; j2 < dim; j2++) {
        T prod1 = N(i2, j2) * N(i, i);
        T prod2 = N(i2, i) * N(i, j2);
        upd(prod1);
        upd(prod2);
        T val = prod1 - prod2;
        upd(val);
        T quot = val / prev;
        upd(quot);
        N(i2, j2) = quot;
      }
    }
    prev = N(i, i);
  }
  return {true, std::move(D), std::move(Gscal), std::move(d), std::move(N),
          std::move(M)};
}

/*
  A priori magnitude bound of every quantity the scaled fraction-free
  enumeration forms: the Bareiss data and its intermediates, the state
  Ttil / Bvec / Hval, the products of the range tests including the seeding
  slack of Bound_Floor / Bound_Ceil, and the norm recomputation of the
  sanity block. All of it is simple linear algebra evaluated in exact
  arithmetic, so the bounding itself cannot overflow.
 */
template <typename T>
T shvec_fast_master_bound(ShvecFastPrep<T> const &prep, MyVector<T> const &Cnum,
                          T const &cden, T const &bound_scal) {
  int dim = prep.Gscal.rows();
  MyVector<T> const &d = prep.d;
  MyMatrix<T> const &N = prep.Nmat;
  T M = prep.maxBareiss;
  auto upd = [&](T const &v) -> void {
    T a = T_abs(v);
    if (a > M) {
      M = a;
    }
  };
  // The two searches below produce BOUNDS, so overshooting by a factor of
  // two is harmless: galloping keeps them logarithmic and robust for any
  // magnitude, where a double-seeded search would overflow its seed on the
  // values met before the fits-test has rejected them.
  // Some s with sqrt(v) <= s <= 2 * sqrt(v) + 2.
  auto isqrt_ceil = [&](T const &v) -> T {
    T s(1);
    while (s * s < v) {
      s += s;
    }
    return s;
  };
  // Some t with a / b <= t <= 2 * a / b + 2, for a >= 0, b > 0.
  auto ceil_quotient = [&](T const &a, T const &b) -> T {
    T t(1);
    while (t * b < a) {
      t += t;
    }
    return t;
  };
  // Backward coordinate bounds: Ytot(i) bounds |x_S(i) + Cnum(i)| on any
  // in-range branch, from d(i) * (x_S + Cnum)(i) = Hval(i) - sum_j N(i,j)
  // * (x_S + Cnum)(j) and |Hval(i)| <= sqrt(d(i-1) * d(i) * bound_scal).
  MyVector<T> Ytot(dim);
  for (int i = dim - 1; i >= 0; i--) {
    T dm1 = (i > 0) ? d(i - 1) : T(1);
    T K = dm1 * d(i) * bound_scal;
    upd(K);
    T acc = isqrt_ceil(K);
    for (int j = i + 1; j < dim; j++) {
      acc += T_abs(N(i, j)) * Ytot(j);
    }
    Ytot(i) = ceil_quotient(acc, d(i));
    upd(Ytot(i));
  }
  // The runtime magnitudes of the enumeration itself.
  MyVector<T> HMax(dim);
  for (int i = 0; i < dim; i++) {
    T aCnum = T_abs(Cnum(i));
    T BvMax = d(i) * aCnum;
    for (int j = i + 1; j < dim; j++) {
      BvMax += T_abs(N(i, j)) * Ytot(j);
    }
    upd(BvMax);
    HMax(i) = d(i) * (Ytot(i) + aCnum) + BvMax;
    upd(HMax(i));
    // The double-seeded search of Bound_Floor / Bound_Ceil can probe a
    // bounded number of steps beyond the true range before correcting.
    T probe = HMax(i) + d(i) * cden * T(4096);
    upd(probe * probe);
    upd(d(i) * bound_scal);
  }
  for (int i = 0; i + 1 < dim; i++) {
    upd(d(i) * d(i + 1) * bound_scal + HMax(i + 1) * HMax(i + 1));
  }
  // The norm recomputation of the SANITY_CHECK_SHVEC block.
  T maxG(0), maxY(0);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      T a = T_abs(prep.Gscal(i, j));
      if (a > maxG) {
        maxG = a;
      }
    }
    T y = Ytot(i) + T_abs(Cnum(i));
    if (y > maxY) {
      maxY = y;
    }
  }
  upd(maxG * maxY * maxY * T(dim) * T(dim));
  return M;
}

/*
  Whether every intermediate magnitude of the enumeration fits the native
  type. For int32 the products are computed in the 32-bit int of the integer
  promotion, so the limit is 2^29 with margin; for int64 it is 2^61.
 */
template <typename Tfast, typename T>
bool shvec_fast_fits(T const &Mbound) {
  int64_t limit =
      (sizeof(Tfast) == 4) ? (int64_t(1) << 29) : (int64_t(1) << 61);
  return Mbound < UniversalScalarConversion<T, int64_t>(limit);
}

/*
  The enumeration at a native integer type. The caller has already proven
  through shvec_fast_master_bound and shvec_fast_fits that every
  intermediate fits the type.
 */
template <typename Tfast, typename T, typename Tint, typename Finsert>
bool computeIt_fast(ShvecFastPrep<T> const &prep, MyVector<T> const &Cnum,
                    T const &cden, T const &bound_scal, T const &scal,
                    bool const &central, Finsert f_insert) {
  MyMatrix<T> const &Gscal = prep.Gscal;
  int dim = Gscal.rows();
  MyMatrix<Tfast> G_fast(dim, dim);
  MyMatrix<Tfast> N_fast(dim, dim);
  MyVector<Tfast> d_fast(dim);
  // The Bareiss data is the preparation's, converted: recomputing it at the
  // native type would repeat work already done exactly.
  for (int i = 0; i < dim; i++) {
    d_fast(i) = UniversalScalarConversion<Tfast, T>(prep.d(i));
    for (int j = 0; j < dim; j++) {
      G_fast(i, j) = UniversalScalarConversion<Tfast, T>(Gscal(i, j));
      N_fast(i, j) = UniversalScalarConversion<Tfast, T>(prep.Nmat(i, j));
    }
  }
  FullGramInfo<Tfast> request_fast(std::move(G_fast), std::move(d_fast),
                                   std::move(N_fast));
  MyVector<Tfast> Cnum_fast(dim);
  for (int i = 0; i < dim; i++) {
    Cnum_fast(i) = UniversalScalarConversion<Tfast, T>(Cnum(i));
  }
  Tfast cden_fast = UniversalScalarConversion<Tfast, T>(cden);
  Tfast bound_fast = UniversalScalarConversion<Tfast, T>(bound_scal);
  bool descale = (scal != T(1));
  auto f_insert_fast = [&](MyVector<Tfast> const &x,
                           Tfast const &norm_s) -> bool {
    int n = x.size();
    MyVector<Tint> x_ret(n);
    for (int u = 0; u < n; u++) {
      x_ret(u) = UniversalScalarConversion<Tint, Tfast>(x(u));
    }
    T norm = UniversalScalarConversion<T, Tfast>(norm_s);
    if (descale) {
      norm = norm / scal;
    }
    return f_insert(x_ret, norm);
  };
  auto f_set_bound = [&](const Tfast &K, const Tfast &B, const Tfast &den,
                         [[maybe_unused]] const MyVector<Tfast> &x,
                         [[maybe_unused]] const int &i, Tfast &upper,
                         Tfast &lower) -> void {
    upper = Bound_Floor<Tfast, Tfast>(K, B, den);
    lower = Bound_Ceil<Tfast, Tfast>(K, B, den);
  };
  return computeIt_Gen_Kernel<Tfast, Tfast, decltype(f_insert_fast),
                              decltype(f_set_bound)>(
      request_fast, Cnum_fast, cden_fast, central, bound_fast, f_insert_fast,
      f_set_bound);
}

/*
  The main enumeration entry: attempts the native-integer fast path and
  falls back on the exact one. The optional prep pointer lets a caller
  amortize the Gram-only preparation across many enumerations of the same
  matrix, which is what CVPSolver does; without it the preparation is
  computed for the call.
 */
template <typename T, typename Tint, typename Finsert>
inline bool computeIt(const FullGramInfo<T> &request, MyVector<T> const &coset,
                      bool const &central, const T &bound, Finsert f_insert,
                      ShvecFastPrep<T> const *prep_ptr = nullptr) {
  if constexpr (shvec_has_fast_path<T>::value) {
    if (bound >= 0) {
      int dim = request.dim;
      std::optional<ShvecFastPrep<T>> prep_local;
      ShvecFastPrep<T> const *prep = prep_ptr;
      if (!prep) {
        prep_local = ComputeShvecFastPrep(request.gram_matrix);
        prep = &*prep_local;
      }
      if (prep->usable) {
        T cden(1);
        for (int i = 0; i < dim; i++) {
          cden = LCMpair(cden, GetDenominator(coset(i)));
        }
        T scal = prep->D * cden * cden;
        T bprod = scal * bound;
        // Early exit on a hopeless bound: the master bound dominates the
        // scaled bound (M >= d(dim-1) * bound_scal >= bound_scal), so a
        // scaled bound at or beyond 2^61 (~2.3e18) can never pass the int64
        // fits-test. The double image is compared against 1.0e18, whose
        // margin absorbs the rounding of the conversion.
        if (UniversalScalarConversion<double, T>(bprod) < 1.0e18) {
          // Largest integer-valued bound below bprod: flooring to the norm
          // grid changes no enumerated vector. Over a ring the value is
          // already integral.
          T bound_scal = bprod;
          if constexpr (is_ring_field<T>::value) {
            bound_scal = UniversalFloorScalarInteger<T, T>(bprod);
          }
          MyVector<T> Cnum(dim);
          for (int i = 0; i < dim; i++) {
            Cnum(i) = cden * coset(i);
          }
          T Mbound = shvec_fast_master_bound(*prep, Cnum, cden, bound_scal);
          if (shvec_fast_fits<int32_t, T>(Mbound)) {
            return computeIt_fast<int32_t, T, Tint, Finsert>(
                *prep, Cnum, cden, bound_scal, scal, central, f_insert);
          }
          if (shvec_fast_fits<int64_t, T>(Mbound)) {
            return computeIt_fast<int64_t, T, Tint, Finsert>(
                *prep, Cnum, cden, bound_scal, scal, central, f_insert);
          }
        }
      }
    }
  }
  auto f_set_bound = [&](const T &K, const T &B, const T &den,
                         [[maybe_unused]] const MyVector<Tint> &x,
                         [[maybe_unused]] const int &i, Tint &upper,
                         Tint &lower) -> void {
    upper = Bound_Floor<T, Tint>(K, B, den);
    lower = Bound_Ceil<T, Tint>(K, B, den);
  };
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(
      request, coset, central, bound, f_insert, f_set_bound);
}

template <typename T>
T get_initial_minimum(const FullGramInfo<T> &request, MyVector<T> const &C,
                      bool const &central) {
  int dim = request.dim;
  if (!central) {
    T eNorm(0);
    for (int i = 0; i < dim; i++)
      for (int j = 0; j < dim; j++)
        eNorm += request.gram_matrix(i, j) * C(i) * C(j);
    return eNorm;
  } else {
    T eMin = request.gram_matrix(0, 0);
    for (int i = 1; i < dim; i++) {
      T diag_val = request.gram_matrix(i, i);
      if (eMin > diag_val) {
        eMin = diag_val;
      }
    }
    return eMin;
  }
}

template <typename T, typename Tint>
T_shvec_info<T, Tint> compute_minimum(const FullGramInfo<T> &request,
                                      MyVector<T> const &coset,
                                      bool const &central,
                                      ShvecFastPrep<T> const *prep = nullptr) {
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: compute_minimum, begin\n";
#endif
  std::vector<MyVector<Tint>> short_vectors;
  T minimum = get_initial_minimum(request, coset, central);
  while (true) {
#ifdef DEBUG_SHVEC
    std::cerr << "SHVEC: Before computeIt (in compute_minimum while loop)\n";
#endif
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == minimum) {
        short_vectors.push_back(V);
        if (central) {
          short_vectors.push_back(-V);
        }
        return true;
      } else {
        short_vectors.clear();
        minimum = min;
        return false;
      }
    };
    bool result = computeIt<T, Tint, decltype(f_insert)>(
        request, coset, central, minimum, f_insert, prep);
    if (result) {
      break;
    }
  }
  return {std::move(short_vectors), std::move(minimum)};
}

template <typename T, typename Tint>
T_shvec_info<T, Tint>
compute_minimum_limit(const FullGramInfo<T> &request, MyVector<T> const &coset,
                      bool const &central, std::optional<size_t> const &limit,
                      ShvecFastPrep<T> const *prep = nullptr) {
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: compute_minimum_limit, begin\n";
#endif
  std::vector<MyVector<Tint>> short_vectors;
  T minimum = get_initial_minimum(request, coset, central);
  while (true) {
#ifdef DEBUG_SHVEC
    std::cerr
        << "SHVEC: Before computeIt (in compute_minimum_limit while loop)\n";
#endif
    size_t n_iter = 0;
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == minimum) {
        short_vectors.push_back(V);
        if (central) {
          short_vectors.push_back(-V);
        }
        if (limit) {
          size_t const &limit_val = *limit;
          if (limit_val <= n_iter) {
            return false;
          }
        }
        return true;
      } else {
        short_vectors.clear();
        minimum = min;
        return false;
      }
    };
    bool result = computeIt<T, Tint, decltype(f_insert)>(
        request, coset, central, minimum, f_insert, prep);
    if (result) {
      break;
    }
  }
  return {std::move(short_vectors), std::move(minimum)};
}

template <typename Tint> struct ResultShortest {
  std::vector<MyVector<Tint>> shortest;
  std::optional<MyVector<Tint>> better_vector;
};

template <typename T, typename Tint>
ResultShortest<Tint>
compute_test_shortest(const FullGramInfo<T> &request, MyVector<T> const &coset,
                      bool const &central, T const &bound) {
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: compute_test_shortest, begin\n";
#endif
  std::vector<MyVector<Tint>> shortest;
  std::optional<MyVector<Tint>> better_vector;
  auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
    if (min == bound) {
      shortest.push_back(V);
      if (central) {
        shortest.push_back(-V);
      }
      return true;
    } else {
      shortest.clear();
      better_vector = V;
      return false;
    }
  };
  (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central, bound,
                                               f_insert);
  return {std::move(shortest), std::move(better_vector)};
}

template <typename T, typename Tint> struct CVPSolver {
public:
  MyMatrix<T> GramMat;

private:
  int dim;
  std::ostream &os;
  LLLreduction<T, Tint> eRec;
  MyMatrix<T> Q_T;
  FullGramInfo<T> request;
  // The Gram-only part of the fast-path preparation, computed once and
  // shared by every enumeration of this solver.
  std::optional<ShvecFastPrep<T>> fast_prep;

public:
  CVPSolver(MyMatrix<T> const &_GramMat, std::ostream &_os)
      : GramMat(_GramMat), dim(GramMat.rows()), os(_os),
        eRec(LLLreducedBasisDual<T, Tint>(GramMat, os)) {
    MyMatrix<Tint> Q_i = Inverse(eRec.Pmat);
    Q_T = UniversalMatrixConversion<T, Tint>(Q_i);
    request = FullGramInfo<T>(eRec.GramMatRed);
    if constexpr (shvec_has_fast_path<T>::value) {
      fast_prep = ComputeShvecFastPrep(eRec.GramMatRed);
    }
  }
  ShvecFastPrep<T> const *get_prep() const {
    if (fast_prep) {
      return &*fast_prep;
    }
    return nullptr;
  }
  T comp_norm_vect(MyVector<Tint> const &x) const {
    MyVector<T> eDiff(dim);
    for (int i = 0; i < dim; i++) {
      eDiff(i) = UniversalScalarConversion<T, Tint>(x(i));
    }
    return EvaluationQuadForm<T, T>(GramMat, eDiff);
  }
  T comp_norm_diff(MyVector<Tint> const &x, MyVector<T> const &v) const {
    MyVector<T> eDiff(dim);
    for (int i = 0; i < dim; i++) {
      eDiff(i) = UniversalScalarConversion<T, Tint>(x(i)) - v(i);
    }
#ifdef DEBUG_SHVEC
    std::cerr << "SHVEC: Debug, dim=" << dim << "\n";
    std::cerr << "SHVEC: Debug, |GramMat|=" << GramMat.rows() << " / " << GramMat.cols() << "\n";
    std::cerr << "SHVEC: Debug, eDiff=\n";
    WriteVector(std::cerr, eDiff);
    std::cerr << "SHVEC: Debug, GramMat=\n";
    WriteMatrix(std::cerr, GramMat);
#endif
    return EvaluationQuadForm<T, T>(GramMat, eDiff);
  }
  ResultShortest<Tint> ShortestAtDistance(MyVector<T> const &eV,
                                          T const &TheNorm) const {
    if (IsIntegralVector(eV)) {
      MyVector<Tint> eV_tint = UniversalVectorConversion<Tint, T>(eV);
      return {{}, eV_tint};
    }
    MyVector<T> cosetRed = -Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const &coset = ePair.second;
    bool central = false;
    ResultShortest<Tint> res =
      compute_test_shortest<T, Tint>(request, coset, central, TheNorm);
    if (res.better_vector) {
      MyVector<Tint> const &short_vec = *res.better_vector;
      MyVector<Tint> x = eRec.Pmat.transpose() * (short_vec - ePair.first);
#ifdef SANITY_CHECK_SHVEC
      if (TheNorm <= comp_norm_diff(x, eV)) {
        std::cerr << "The vector should be shorter\n";
        throw TerminalException{1};
      }
#endif
      return {{}, x};
    }
    std::vector<MyVector<Tint>> shortest;
    for (auto &short_vec : res.shortest) {
      MyVector<Tint> x = eRec.Pmat.transpose() * (short_vec - ePair.first);
      shortest.push_back(x);
#ifdef SANITY_CHECK_SHVEC
      if (TheNorm != comp_norm_diff(x, eV)) {
        std::cerr << "Inconsistecy error in the norms\n";
        throw TerminalException{1};
      }
#endif
    }
    return {std::move(shortest), {}};
  }
  // Differences between points of the point set, spanning the space,
  // used to seed the search for a Delaunay polytope: they are the moves
  // of the lattice itself here. The Delaunay geometry asks for them
  // rather than assuming them, since for a periodic point set the moves
  // are those of the period lattice, not the unit vectors.
  std::vector<MyVector<T>> get_seed_differences() const {
    std::vector<MyVector<T>> ret;
    for (int i = 0; i < dim; i++) {
      MyVector<T> V1 = ZeroVector<T>(dim);
      V1(i) = 1;
      ret.emplace_back(std::move(V1));
      MyVector<T> V2 = ZeroVector<T>(dim);
      V2(i) = -1;
      ret.emplace_back(std::move(V2));
    }
    return ret;
  }
  resultCVP<T, Tint> nearest_vectors(MyVector<T> const &eV) const {
    if (IsIntegralVector(eV)) {
      T TheNorm(0);
      MyMatrix<Tint> ListVect(1, dim);
      for (int i = 0; i < dim; i++)
        ListVect(0, i) = UniversalScalarConversion<Tint, T>(eV(i));
      return {std::move(TheNorm), std::move(ListVect)};
    }
    MyVector<T> cosetRed = -Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const &coset = ePair.second;
    bool central = false;
    T_shvec_info<T, Tint> info =
        compute_minimum<T, Tint>(request, coset, central, get_prep());
    T TheNorm = info.minimum;
    int nbVect = info.short_vectors.size();
    MyMatrix<Tint> ListClos(nbVect, dim);
    for (int iVect = 0; iVect < nbVect; iVect++) {
      MyVector<Tint> x =
          eRec.Pmat.transpose() * (info.short_vectors[iVect] - ePair.first);
#ifdef SANITY_CHECK_SHVEC
      if (TheNorm != comp_norm_diff(x, eV)) {
        std::cerr << "Inconsistecy error in the norms\n";
        throw TerminalException{1};
      }
#endif
      for (int i = 0; i < dim; i++) {
        ListClos(iVect, i) = x(i);
      }
    }
    return {TheNorm, std::move(ListClos)};
  }
  resultCVP<T, Tint>
  nearest_vectors_limit(MyVector<T> const &eV,
                        std::optional<size_t> const &limit) const {
    if (IsIntegralVector(eV)) {
      T TheNorm(0);
      MyMatrix<Tint> ListVect(1, dim);
      for (int i = 0; i < dim; i++)
        ListVect(0, i) = UniversalScalarConversion<Tint, T>(eV(i));
      return {std::move(TheNorm), std::move(ListVect)};
    }
    MyVector<T> cosetRed = -Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const &coset = ePair.second;
    bool central = false;
    T_shvec_info<T, Tint> info =
        compute_minimum_limit<T, Tint>(request, coset, central, limit,
                                       get_prep());
    T TheNorm = info.minimum;
    int nbVect = info.short_vectors.size();
    MyMatrix<Tint> ListClos(nbVect, dim);
    for (int iVect = 0; iVect < nbVect; iVect++) {
      MyVector<Tint> x =
          eRec.Pmat.transpose() * (info.short_vectors[iVect] - ePair.first);
#ifdef SANITY_CHECK_SHVEC
      if (TheNorm != comp_norm_diff(x, eV)) {
        std::cerr << "Inconsistecy error in the norms\n";
        throw TerminalException{1};
      }
#endif
      for (int i = 0; i < dim; i++) {
        ListClos(iVect, i) = x(i);
      }
    }
    return {TheNorm, std::move(ListClos)};
  }
  Tshortest<T, Tint> shortest_vectors() const {
    MyVector<T> coset = ZeroVector<T>(dim);
    bool central = true;
    T_shvec_info<T, Tint> info =
        compute_minimum<T, Tint>(request, coset, central, get_prep());
    T TheNorm = info.minimum;
    int nbVect = info.short_vectors.size();
    MyMatrix<Tint> ListClos(nbVect, dim);
    for (int iVect = 0; iVect < nbVect; iVect++) {
      MyVector<Tint> x = eRec.Pmat.transpose() * info.short_vectors[iVect];
#ifdef SANITY_CHECK_SHVEC
      if (TheNorm != comp_norm_vect(x)) {
        std::cerr << "Inconsistecy error in the norms\n";
        throw TerminalException{1};
      }
#endif
      for (int i = 0; i < dim; i++) {
        ListClos(iVect, i) = x(i);
      }
    }
    return {TheNorm, std::move(ListClos)};
  }
  template <typename Fins>
  void fixed_dist_vectors_f(MyVector<T> const &eV, Fins f_ins,
                            T const &TheNorm) const {
    MyVector<T> cosetRed = -Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const &coset = ePair.second;
    bool central = false;
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == TheNorm) {
        MyVector<Tint> x = eRec.Pmat.transpose() * (V - ePair.first);
#ifdef SANITY_CHECK_SHVEC
        if (TheNorm != comp_norm_diff(x, eV)) {
          std::cerr << "Inconsistecy error in the norms\n";
          throw TerminalException{1};
        }
#endif
        f_ins(x);
      }
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central,
                                                 TheNorm, f_insert,
                                                 get_prep());
  }
  std::vector<MyVector<Tint>> fixed_dist_vectors(MyVector<T> const &eV,
                                                 T const &TheNorm) const {
    std::vector<MyVector<Tint>> ListVect;
    auto f_ins = [&](MyVector<Tint> const &x) -> void {
      ListVect.push_back(x);
    };
    fixed_dist_vectors_f<decltype(f_ins)>(eV, f_ins, TheNorm);
    return ListVect;
  }
  std::vector<MyVector<Tint>> fixed_norm_vectors(T const &TheNorm) const {
    MyVector<T> coset = ZeroVector<T>(dim);
    bool central = true;
    std::vector<MyVector<Tint>> ListVect;
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == TheNorm) {
        MyVector<Tint> x = eRec.Pmat.transpose() * V;
#ifdef SANITY_CHECK_SHVEC
        if (TheNorm != comp_norm_vect(x)) {
          std::cerr << "Inconsistecy error in the norms\n";
          throw TerminalException{1};
        }
#endif
        ListVect.emplace_back(std::move(x));
      }
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central,
                                                 TheNorm, f_insert,
                                                 get_prep());
    return ListVect;
  }
  std::vector<MyVector<Tint>> at_most_dist_vectors(MyVector<T> const &eV,
                                                   T const &MaxNorm) const {
    MyVector<T> cosetRed = -Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const &coset = ePair.second;
    bool central = false;
    std::vector<MyVector<Tint>> ListVect;
    auto f_insert = [&](const MyVector<Tint> &V,
                        [[maybe_unused]] const T &min) -> bool {
      MyVector<Tint> x = eRec.Pmat.transpose() * (V - ePair.first);
#ifdef SANITY_CHECK_SHVEC
      if (MaxNorm < comp_norm_diff(x, eV)) {
        std::cerr << "Inconsistecy error in the norms\n";
        throw TerminalException{1};
      }
#endif
      ListVect.emplace_back(std::move(x));
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central,
                                                 MaxNorm, f_insert,
                                                 get_prep());
    return ListVect;
  }
  resultCVP<T, Tint>
  increase_distance_vectors(MyVector<T> const &eV,
                            std::optional<T> const &opt_norm) const {
    if (!opt_norm) {
      return nearest_vectors(eV);
    }
    int n = eV.size();
    T n_T = T(n);
    T factor = T(1) + T(1) / n_T;
    T previous_norm = *opt_norm;
    T eff_norm = previous_norm * factor + T(1);
    while (true) {
      std::vector<MyVector<Tint>> ListV = at_most_dist_vectors(eV, eff_norm);
      std::vector<MyVector<Tint>> list_previous;
      std::vector<MyVector<Tint>> list_above;
      T above_norm(0);
      for (auto &fV : ListV) {
        T norm = comp_norm_diff(fV, eV);
        if (norm <= previous_norm) {
          list_previous.push_back(fV);
        } else {
          if (above_norm == 0 || norm < above_norm) {
            list_above.clear();
            list_above.push_back(fV);
            above_norm = norm;
          } else {
            if (norm == above_norm) {
              list_above.push_back(fV);
            }
          }
        }
      }
      if (above_norm != T(0)) {
        int n_vect = list_previous.size() + list_above.size();
        MyMatrix<Tint> ListVect(n_vect, n);
        int pos = 0;
        auto f_insert = [&](MyVector<Tint> const &v) -> void {
          for (int i = 0; i < n; i++) {
            ListVect(pos, i) = v(i);
          }
          pos++;
        };
        for (auto &fV : list_previous) {
          f_insert(fV);
        }
        for (auto &fV : list_above) {
          f_insert(fV);
        }
        return {std::move(above_norm), std::move(ListVect)};
      }
      eff_norm = eff_norm * factor;
    }
  }
  std::vector<MyVector<Tint>> at_most_norm_vectors(T const &MaxNorm) const {
    MyVector<T> coset = ZeroVector<T>(dim);
    bool central = true;
    std::vector<MyVector<Tint>> ListVect;
    auto f_insert = [&](const MyVector<Tint> &V,
                        [[maybe_unused]] const T &min) -> bool {
      MyVector<Tint> x = eRec.Pmat.transpose() * V;
#ifdef SANITY_CHECK_SHVEC
      if (MaxNorm < comp_norm_vect(x)) {
        std::cerr << "Inconsistecy error in the norms\n";
        throw TerminalException{1};
      }
#endif
      ListVect.emplace_back(std::move(x));
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central,
                                                 MaxNorm, f_insert,
                                                 get_prep());
    return ListVect;
  }
};

/*
  The wisdom of applying LLL all the time can be discussed.
  However, it usually improves the situation, so we do it.
  ---
  If the computation repeats itself with the same Gram matrix,
  then building the CVPSolver should be the right approach.
  ---
  If we have plenty of Gram matrices to consider and computing
  the DualLLL is not a good thing, then you have to write
  specific code.
 */
template <typename T, typename Tint>
resultCVP<T, Tint> NearestVectors(MyMatrix<T> const &GramMat,
                                  MyVector<T> const &eV, std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  return solver.nearest_vectors(eV);
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
FindFixedDistVectors(const MyMatrix<T> &GramMat, const MyVector<T> &eV,
                     const T &norm, std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  return solver.fixed_dist_vectors(eV, norm);
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>>
FindAtMostDistVectors(const MyMatrix<T> &GramMat, const MyVector<T> &eV,
                      const T &norm, std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  return solver.at_most_dist_vectors(eV, norm);
}

// Returns half the vector below a specific bound.
template <typename T, typename Tint>
std::vector<MyVector<Tint>> computeLevel_GramMat(MyMatrix<T> const &GramMat,
                                                 T const &bound,
                                                 std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  std::vector<MyVector<Tint>> full_list = solver.at_most_norm_vectors(bound);
  int n_vect = full_list.size() / 2;
  std::vector<MyVector<Tint>> short_vectors;
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    int pos = 2 * i_vect;
    MyVector<Tint> V = full_list[pos];
    short_vectors.push_back(V);
  }
  return short_vectors;
}

template <typename T, typename Tint>
MyMatrix<Tint> T_ShortVector_fixed(MyMatrix<T> const &GramMat,
                                   T const &SpecNorm, std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  std::vector<MyVector<Tint>> ListVect = solver.fixed_norm_vectors(SpecNorm);
  int dim = GramMat.rows();
  return MatrixFromVectorFamilyDim(dim, ListVect);
}

template <typename T, typename Tint>
Tshortest<T, Tint> T_ShortestVector(MyMatrix<T> const &GramMat,
                                    std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  return solver.shortest_vectors();
}

template <typename T, typename Tint>
Tshortest<T, Tint> T_ShortestVectorHalf(MyMatrix<T> const &GramMat,
                                        std::ostream &os) {
  Tshortest<T, Tint> rec_shv = T_ShortestVector<T, Tint>(GramMat, os);
  return shortest_get_half(rec_shv);
}

// clang-format off
#endif  // SRC_LATT_SHVEC_EXACT_H_
// clang-format on
