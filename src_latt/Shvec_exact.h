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

#ifdef TIMINGS
#define TIMINGS_SHVEC
#endif

template <typename T> struct FullGramInfo {
  int dim;
  MyMatrix<T> gram_matrix;
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


// We return floor(sqrt(A) + epsilon + B)
template <typename T> int Infinitesimal_Floor_V1(T const &a, T const &b) {
  double epsilon = 0.000000001;
#ifdef SANITY_CHECK_SHVEC
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Floor_V1\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  double a_doubl = UniversalScalarConversion<double, T>(a);
  double b_doubl = UniversalScalarConversion<double, T>(b);
  double alpha = sqrt(a_doubl) + epsilon + b_doubl;
  double eD1 = floor(alpha);
  long int eD2 = lround(eD1);
  int eD3 = eD2;
  return eD3;
}

template <typename T> int Infinitesimal_Ceil_V1(T const &a, T const &b) {
  double epsilon = 0.000000001;
#ifdef SANITY_CHECK_SHVEC
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Ceil_V1\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  double a_doubl = UniversalScalarConversion<double, T>(a);
  double b_doubl = UniversalScalarConversion<double, T>(b);
  double alpha = - sqrt(a_doubl) - epsilon + b_doubl;
  double eD1 = ceil(alpha);
  long int eD2 = lround(eD1);
  int eD3 = eD2;
  return eD3;
}

// We return floor(sqrt(a) + b)
// n=floor(sqrt(a) + b) is equivalent to
// n<= sqrt(a) + b < n+1
// And so to n - b <= sqrt(a) (and opposite for n+1)
// And so to (n-b)^2 <= a
template <typename T, typename Tint>
Tint Infinitesimal_Floor(T const &a, T const &b) {
#ifdef SANITY_CHECK_SHVEC
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Floor\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  double a_doubl = UniversalScalarConversion<double, T>(a);
  double b_doubl = UniversalScalarConversion<double, T>(b);
  double alpha = sqrt(a_doubl) + b_doubl;
  double eD1 = floor(alpha);
  long int eD2 = lround(eD1);
  Tint eReturn = eD2;
  auto f = [&](Tint const &x) -> bool {
    T eDiff = UniversalScalarConversion<T, Tint>(x) - b;
    if (eDiff <= 0)
      return true;
    if (eDiff * eDiff <= a)
      return true;
    return false;
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
    }
    if (test2) {
      test1 = test2;
      test2 = f(eReturn + 2);
      eReturn++;
    }
  }
  return eReturn;
}

// We return floor(sqrt(a) + b)
// n=ceil(-sqrt(a) + b) is equivalent to
// n-1 < -sqrt(a) + b <= n
// And so to -sqrt(a) <= n - b  (and opposite for n-1)
// And so to (n-b)^2 <= a (and opposite for n-1)
template <typename T, typename Tint>
Tint Infinitesimal_Ceil(T const &a, T const &b) {
#ifdef SANITY_CHECK_SHVEC
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Ceil\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  double a_doubl = UniversalScalarConversion<double, T>(a);
  double b_doubl = UniversalScalarConversion<double, T>(b);
  double alpha = -sqrt(a_doubl) + b_doubl;
  double eD1 = ceil(alpha);
  long int eD2 = lround(eD1);
  Tint eReturn = eD2;
  auto f = [&](Tint const &x) -> bool {
    T eDiff = UniversalScalarConversion<T, Tint>(x) - b;
    if (eDiff >= 0)
      return true;
    if (eDiff * eDiff <= a)
      return true;
    return false;
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
    }
    if (!test2) {
      test1 = test2;
      test2 = f(eReturn + 1);
      eReturn++;
    }
  }
  return eReturn;
}

template <typename T, typename Tint, typename Finsert, typename Fsetbound>
bool computeIt_Gen_Kernel(const FullGramInfo<T> &request,
                          MyVector<T> const& C,
                          bool const& central,
                          const T &bound,
                          Finsert f_insert, Fsetbound f_set_bound) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int i, j;
  int dim = request.dim;
  MyVector<Tint> Upper(dim);
  MyVector<T> Trem(dim);
  MyVector<T> U(dim);
  MyVector<Tint> x(dim);
  auto is_x_zero=[&]() -> bool {
    for (int u=0; u<dim; u++) {
      if (x(u) != 0) {
        return false;
      }
    }
    return true;
  };
  MyVector<T> x_T(dim);
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
  MyMatrix<T> q = request.gram_matrix;
  for (i = 0; i < dim; i++) {
    for (j = i + 1; j < dim; j++) {
      q(i, j) = q(i, j) / q(i, i);
      q(j, i) = q(j, i) / q(i, i);
    }
    for (int i2 = i + 1; i2 < dim; i2++)
      for (int j2 = i + 1; j2 < dim; j2++)
        q(i2, j2) -= q(i, i) * q(i, i2) * q(i, j2);
#ifdef DEBUG_SHVEC_MATRIX
    std::cerr << "SHVEC: diag q=" << q(i, i) << "\n";
    for (int j = i + 1; j < dim; j++)
      std::cerr << "   j=" << j << " q=" << q(i, j) << "\n";
#endif
  }
  bool needs_new_bound = true;
  i = dim - 1;
  if (bound < 0) {
    return true;
  }
  Trem(i) = bound;
  U(i) = 0;
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: Before while loop\n";
#endif
#ifdef DEBUG_SHVEC_VECTOR
  size_t n_vector = 0;
#endif
  T eQuot, eSum, hVal, eNorm;
  while (true) {
    if (needs_new_bound) {
      eQuot = Trem(i) / q(i, i);
      eSum = -U(i) - C(i);
      f_set_bound(eQuot, eSum, q, x, i, Upper(i), x(i));
      x_T(i) = UniversalScalarConversion<T, Tint>(x(i));
      needs_new_bound = false;
    } else {
      x(i) += 1;
      x_T(i) += 1;
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
        hVal = x_T(0) + C(0) + U(0);
        eNorm = bound - Trem(0) + q(0, 0) * hVal * hVal;
#ifdef SANITY_CHECK_SHVEC
        T norm(0);
        for (int i2 = 0; i2 < dim; i2++)
          for (int j2 = 0; j2 < dim; j2++)
            norm += g(i2, j2) * (x_T(i2) + C(i2)) * (x_T(j2) + C(j2));
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
        i--;
        U(i) = 0;
        for (j = i + 1; j < dim; j++)
          U(i) += q(i, j) * (x_T(j) + C(j));
        hVal = x_T(i + 1) + C(i + 1) + U(i + 1);
        Trem(i) = Trem(i + 1) - q(i + 1, i + 1) * hVal * hVal;
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

template <typename T, typename Tint, typename Finsert, typename Fsetbound>
inline typename std::enable_if<is_ring_field<T>::value, bool>::type
computeIt_Gen(const FullGramInfo<T> &request,
              MyVector<T> const& coset,
              bool const& central,
              const T &bound,
              Finsert f_insert, Fsetbound f_set_bound) {
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: computeIt (field case)\n";
#endif
  return computeIt_Gen_Kernel<T, Tint, Finsert, Fsetbound>(request, coset, central, bound, f_insert, f_set_bound);
}

template <typename T, typename Tint, typename Finsert, typename Fsetbound>
inline typename std::enable_if<!is_ring_field<T>::value, bool>::type
computeIt_Gen(const FullGramInfo<T> &request,
              MyVector<T> const& coset,
              bool const& central,
              const T &bound,
              Finsert f_insert, Fsetbound f_set_bound) {
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: computeIt (ring case)\n";
#endif
  using Tfield = typename overlying_field<T>::field_type;
  //
  Tfield bound_field = UniversalScalarConversion<Tfield, T>(bound);
  FullGramInfo<Tfield> request_field{
      request.dim,
      UniversalMatrixConversion<Tfield, T>(request.gram_matrix)};
  MyVector<Tfield> coset_field = UniversalVectorConversion<Tfield, T>(coset);
  //
  auto f_insert_field = [&](const MyVector<Tint> &V,
                            const Tfield &min_Tfield) -> bool {
    T min_T = UniversalScalarConversion<T, Tfield>(min_Tfield);
    return f_insert(V, min_T);
  };
  bool retVal =
    computeIt_Gen_Kernel<Tfield, Tint, decltype(f_insert_field), Fsetbound>(request_field, coset_field, central, bound_field, f_insert_field, f_set_bound);
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: computeIt (ring case) exit\n";
#endif
  return retVal;
}

template <typename T, typename Tint, typename Finsert>
bool computeIt_polytope(const FullGramInfo<T> &request,
                        MyVector<T> const& coset,
                        bool const& central,
                        const T &bound,
                        const MyMatrix<T> &FAC, Finsert f_insert,
                        std::ostream &os) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n_rows = FAC.rows();
  int n_col = FAC.cols();
#ifdef SANITY_CHECK_SHVEC_EXACT_POLYTOPE
  if (n_col != request.dim + 1) {
    std::cerr << "SHVEC: request.dim=" << request.dim << " |request.gram_matrix|=" << request.gram_matrix.rows() << "\n";
    std::cerr << "SHVEC: n_col=" << n_col << " request.dim + 1=" << (request.dim + 1) << "\n";
    std::cerr << "SHVEC: Error in the size of FAC\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_SHVEC_EXACT_POLYTOPE
  std::cerr << "Beginning of computeIt_polytope\n";
#endif
  auto f_set_bound = [&](const T &eQuot, const T &eSum,
                         [[maybe_unused]] const MyMatrix<T> &q,
                         const MyVector<Tint> &x, const int &i, Tint &upper,
                         Tint &lower) -> void {
    upper = Infinitesimal_Floor<T, Tint>(eQuot, eSum);
    lower = Infinitesimal_Ceil<T, Tint>(eQuot, eSum);
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
    eSol = CDD_LinearProgramming(FACwork, Vminimize, os);
    if (eSol.DualDefined && eSol.PrimalDefined) {
      // Well defined so we get a potential lower bound
      Tint eLow = UniversalCeilScalarInteger<Tint, T>(eSol.OptimalValue);
      if (eLow > lower)
        lower = eLow;
    }
    if (!eSol.DualDefined && eSol.PrimalDefined) {
      // Infinite direction. Therefore no better bound available
    }
    if (!eSol.PrimalDefined) {
      // No feasible solution. Therefore not feasible.
      // This will lead to a backtrack operation
      upper = lower - 1;
      return;
    }
    //
    Vminimize(1 + i) = -1;
    eSol = CDD_LinearProgramming(FACwork, Vminimize, os);
    if (eSol.DualDefined && eSol.PrimalDefined) {
      // Well defined so we get a potential upper bound
      Tint eUpp = UniversalFloorScalarInteger<Tint, T>(-eSol.OptimalValue);
      if (eUpp < upper)
        upper = eUpp;
    }
    if (!eSol.DualDefined &&
        eSol.PrimalDefined) { // Infinite direction. Therefore no bound
                              // available
    }
    if (!eSol.PrimalDefined) { // No feasible solution. Therefore not feasible.
      upper = lower - 1;       // This will lead to a backtrack operation
      return;
    }
  };
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(request, coset, central, bound, f_insert, f_set_bound);
}

template <typename T, typename Tint, typename Finsert>
inline typename std::enable_if<is_ring_field<T>::value, bool>::type
computeIt(const FullGramInfo<T> &request,
          MyVector<T> const& coset,
          bool const& central,
          const T &bound,
          Finsert f_insert) {
  auto f_set_bound =
      [&](const T &eQuot, const T &eSum, [[maybe_unused]] const MyMatrix<T> &q,
          [[maybe_unused]] const MyVector<Tint> &x,
          [[maybe_unused]] const int &i, Tint &upper, Tint &lower) -> void {
    upper = Infinitesimal_Floor<T, Tint>(eQuot, eSum);
    lower = Infinitesimal_Ceil<T, Tint>(eQuot, eSum);
  };
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(request, coset, central, bound, f_insert, f_set_bound);
}

template <typename T, typename Tint, typename Finsert>
inline typename std::enable_if<!is_ring_field<T>::value, bool>::type
computeIt(const FullGramInfo<T> &request,
          MyVector<T> const& coset,
          bool const& central,
          const T &bound, Finsert f_insert) {
  using Tfield = typename overlying_field<T>::field_type;
  auto f_set_bound = [&](const Tfield &eQuot, const Tfield &eSum,
                         [[maybe_unused]] const MyMatrix<Tfield> &q,
                         [[maybe_unused]] const MyVector<Tint> &x,
                         [[maybe_unused]] const int &i, Tint &upper,
                         Tint &lower) -> void {
    upper = Infinitesimal_Floor<Tfield, Tint>(eQuot, eSum);
    lower = Infinitesimal_Ceil<Tfield, Tint>(eQuot, eSum);
  };
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(request, coset, central, bound, f_insert, f_set_bound);
}

template <typename T>
T get_initial_minimum(const FullGramInfo<T> &request, MyVector<T> const& C, bool const& central) {
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
T_shvec_info<T, Tint> compute_minimum(const FullGramInfo<T> &request, MyVector<T> const& coset, bool const& central) {
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
    bool result = computeIt<T, Tint, decltype(f_insert)>(request, coset, central, minimum, f_insert);
    if (result) {
      break;
    }
  }
  return {short_vectors, minimum};
}

template <typename T, typename Tint>
T_shvec_info<T, Tint> compute_minimum_limit(const FullGramInfo<T> &request, MyVector<T> const& coset, bool const& central, std::optional<size_t> const& limit) {
#ifdef DEBUG_SHVEC
  std::cerr << "SHVEC: compute_minimum_limit, begin\n";
#endif
  std::vector<MyVector<Tint>> short_vectors;
  T minimum = get_initial_minimum(request, coset, central);
  while (true) {
#ifdef DEBUG_SHVEC
    std::cerr << "SHVEC: Before computeIt (in compute_minimum_limit while loop)\n";
#endif
    size_t n_iter = 0;
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == minimum) {
        short_vectors.push_back(V);
        if (central) {
          short_vectors.push_back(-V);
        }
        if (limit) {
          size_t const& limit_val = *limit;
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
    bool result = computeIt<T, Tint, decltype(f_insert)>(request, coset, central, minimum, f_insert);
    if (result) {
      break;
    }
  }
  return {short_vectors, minimum};
}





template <typename Tint> struct ResultShortest {
  std::vector<MyVector<Tint>> shortest;
  std::optional<MyVector<Tint>> better_vector;
};

template <typename T, typename Tint>
ResultShortest<Tint> compute_test_shortest(const FullGramInfo<T> &request, MyVector<T> const& coset, bool const& central, T const& bound) {
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
  return {shortest, better_vector};
}

template <typename T, typename Tint> struct CVPSolver {
public:
  MyMatrix<T> const &GramMat;
private:
  int dim;
  std::ostream &os;
  LLLreduction<T, Tint> eRec;
  MyMatrix<T> Q_T;
  FullGramInfo<T> request;
public:
  CVPSolver(MyMatrix<T> const &_GramMat, std::ostream &_os)
    : GramMat(_GramMat), dim(GramMat.rows()), os(_os),
      eRec(LLLreducedBasisDual<T, Tint>(GramMat, os)) {
    MyMatrix<Tint> Q_i = Inverse(eRec.Pmat);
    Q_T = UniversalMatrixConversion<T, Tint>(Q_i);
    request = FullGramInfo<T>{dim, eRec.GramMatRed};
  }
  T comp_norm_vect(MyVector<Tint> const& x) const {
    MyVector<T> eDiff(dim);
    for (int i = 0; i < dim; i++) {
      eDiff(i) = UniversalScalarConversion<T, Tint>(x(i));
    }
    return EvaluationQuadForm<T, T>(GramMat, eDiff);
  }
  T comp_norm_diff(MyVector<Tint> const& x, MyVector<T> const& v) const {
    MyVector<T> eDiff(dim);
    for (int i = 0; i < dim; i++) {
      eDiff(i) = UniversalScalarConversion<T, Tint>(x(i)) - v(i);
    }
    return EvaluationQuadForm<T, T>(GramMat, eDiff);
  }
  ResultShortest<Tint> ShortestAtDistance(MyVector<T> const &eV,
                                          T const &TheNorm) const {
    if (IsIntegralVector(eV)) {
      MyVector<Tint> eV_tint = UniversalVectorConversion<Tint, T>(eV);
      return {{}, eV_tint};
    }
    MyVector<T> cosetRed = - Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const& coset = ePair.second;
    ResultShortest<Tint> res = compute_test_shortest<T, Tint>(request, coset, TheNorm);
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
    return {shortest, {}};
  }
  resultCVP<T, Tint> nearest_vectors(MyVector<T> const &eV) const {
    if (IsIntegralVector(eV)) {
      T TheNorm(0);
      MyMatrix<Tint> ListVect(1, dim);
      for (int i = 0; i < dim; i++)
        ListVect(0, i) = UniversalScalarConversion<Tint, T>(eV(i));
      return {TheNorm, ListVect};
    }
    MyVector<T> cosetRed = - Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const& coset = ePair.second;
    bool central = false;
    T_shvec_info<T, Tint> info = compute_minimum<T, Tint>(request, coset, central);
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
  resultCVP<T, Tint> nearest_vectors_limit(MyVector<T> const &eV, std::optional<size_t> const& limit) const {
    if (IsIntegralVector(eV)) {
      T TheNorm(0);
      MyMatrix<Tint> ListVect(1, dim);
      for (int i = 0; i < dim; i++)
        ListVect(0, i) = UniversalScalarConversion<Tint, T>(eV(i));
      return {TheNorm, ListVect};
    }
    MyVector<T> cosetRed = - Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const& coset = ePair.second;
    bool central = false;
    T_shvec_info<T, Tint> info = compute_minimum_limit<T, Tint>(request, coset, central, limit);
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
    T_shvec_info<T, Tint> info = compute_minimum<T, Tint>(request, coset, central);
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
  template<typename Fins>
  void fixed_dist_vectors_f(MyVector<T> const &eV, Fins f_ins,
                                                 T const &TheNorm) const {
    MyVector<T> cosetRed = - Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const& coset = ePair.second;
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
    (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central, TheNorm, f_insert);
  }
  std::vector<MyVector<Tint>> fixed_dist_vectors(MyVector<T> const &eV,
                                                 T const &TheNorm) const {
    std::vector<MyVector<Tint>> ListVect;
    auto f_ins = [&](MyVector<Tint> const& x) -> void {
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
    (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central, TheNorm, f_insert);
    return ListVect;
  }
  std::vector<MyVector<Tint>> at_most_dist_vectors(MyVector<T> const &eV,
                                                   T const &MaxNorm) const {
    MyVector<T> cosetRed = -Q_T.transpose() * eV;
    std::pair<MyVector<Tint>, MyVector<T>> ePair =
        ReductionMod1vector<T, Tint>(cosetRed);
    MyVector<T> const& coset = ePair.second;
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
    (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central, MaxNorm, f_insert);
    return ListVect;
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
    (void)computeIt<T, Tint, decltype(f_insert)>(request, coset, central, MaxNorm, f_insert);
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
resultCVP<T, Tint> NearestVectors(MyMatrix<T> const &GramMat, MyVector<T> const &eV, std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  return solver.nearest_vectors(eV);
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>> FindFixedDistVectors(const MyMatrix<T> &GramMat,
                                                 const MyVector<T> &eV,
                                                 const T &norm, std::ostream& os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  return solver.fixed_dist_vectors(eV, norm);
}

template <typename T, typename Tint>
std::vector<MyVector<Tint>> FindAtMostDistVectors(const MyMatrix<T> &GramMat,
                                                  const MyVector<T> &eV,
                                                  const T &norm, std::ostream& os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  return solver.at_most_dist_vectors(eV, norm);
}

// Returns half the vector below a specific bound.
template <typename T, typename Tint>
std::vector<MyVector<Tint>> computeLevel_GramMat(MyMatrix<T> const &GramMat,
                                                 T const &bound, std::ostream& os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  std::vector<MyVector<Tint>> full_list = solver.at_most_norm_vectors(bound);
  int n_vect = full_list.size() / 2;
  std::vector<MyVector<Tint>> short_vectors;
  for (int i_vect=0; i_vect<n_vect; i_vect++) {
    int pos = 2 * i_vect;
    MyVector<Tint> V = full_list[pos];
    short_vectors.push_back(V);
  }
  return short_vectors;
}






template <typename T, typename Tint>
MyMatrix<Tint> T_ShortVector(MyMatrix<T> const &GramMat, T const &MaxNorm,
                             std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  std::vector<MyVector<Tint>> ListVect = solver.at_most_norm_vectors(MaxNorm);
  int dim = GramMat.rows();
  return MatrixFromVectorFamilyDim(dim, ListVect);
}

template <typename T, typename Tint>
MyMatrix<Tint> T_ShortVector_fixed(MyMatrix<T> const &GramMat,
                                   T const &SpecNorm, std::ostream& os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  std::vector<MyVector<Tint>> ListVect = solver.fixed_norm_vectors(SpecNorm);
  int dim = GramMat.rows();
  return MatrixFromVectorFamilyDim(dim, ListVect);
}

template <typename T, typename Tint>
Tshortest<T, Tint> T_ShortestVector(MyMatrix<T> const &GramMat, std::ostream &os) {
  CVPSolver<T, Tint> solver(GramMat, os);
  return solver.shortest_vectors();
}

template <typename T, typename Tint>
Tshortest<T, Tint> T_ShortestVectorHalf(MyMatrix<T> const &GramMat, std::ostream &os) {
  Tshortest<T, Tint> RecSHV = T_ShortestVector<T,Tint>(GramMat, os);
  return shortest_get_half(RecSHV);
}





// clang-format off
#endif  // SRC_LATT_SHVEC_EXACT_H_
// clang-format on
