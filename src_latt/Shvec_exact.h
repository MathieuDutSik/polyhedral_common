// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SHVEC_EXACT_H_
#define SRC_LATT_SHVEC_EXACT_H_

#include "LatticeDefinitions.h"
#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "POLY_LinearProgramming.h"
#include <utility>
#include <vector>

#define CHECK_BASIC_CONSISTENCY
// #define PRINT_DEBUG_INFO
// #define PRINT_DEBUG_INFO_VECTOR

namespace TempShvec_globals {
const int TEMP_SHVEC_MODE_UNDEF = -1;
const int TEMP_SHVEC_MODE_BOUND = 0;
const int TEMP_SHVEC_MODE_SHORTEST_VECTORS = 1;
const int TEMP_SHVEC_MODE_MINIMUM = 2;
const int TEMP_SHVEC_MODE_THETA_SERIES = 3;
const int TEMP_SHVEC_MODE_VINBERG_ALGO = 4;
const int TEMP_SHVEC_MODE_LORENTZIAN = 5;
const int TEMP_SHVEC_MODE_HAN_TRAN = 6;
const int STOP_COMPUTATION = 666;
const int NORMAL_TERMINATION_COMPUTATION = 555;
// clang-format off
}  // namespace TempShvec_globals
// clang-format on

template <typename T> struct T_shvec_request {
  int dim;
  int mode;
  T bound;
  MyVector<T> coset;
  MyMatrix<T> gram_matrix;
  bool central;
};

template <typename T, typename Tint> struct T_shvec_info {
  std::vector<MyVector<Tint>> short_vectors;
  T minimum;
};

template <typename Tint> struct cvp_reduction_info {
  MyVector<Tint> shift;
  MyMatrix<Tint> P;
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
  We apply the P reduction
  Gred = P G P^T
  Q = P^{-1}
  G = Q Gred Q^T
  looking for G[x - coset] = bound
  That is (x - coset) Q Gred Q^T (x- coset)^T = bound
  Or (y - coset_red) Gred (y - coset_red)^T = bound
  with coset_red = coset Q and y = x Q
  coset_red = V_T + V_i
  Thus we have if z = y - V_i
  (z - V_T) Gred (z - v_T)^T = bound
  So y = z + V_i
  and x = (z + V_i) P
 */
template <typename T, typename Tint>
std::pair<T_shvec_request<T>, cvp_reduction_info<Tint>>
GetReducedShvecRequest(T_shvec_request<T> const &request) {
  int dim = request.dim;
  LLLreduction<T, Tint> eRec =
      LLLreducedBasisDual<T, Tint>(request.gram_matrix);
  MyMatrix<Tint> Q_i = Inverse(eRec.Pmat);
  MyMatrix<T> Q_T = UniversalMatrixConversion<T, Tint>(Q_i);
  MyVector<T> cosetRed = Q_T.transpose() * request.coset;
  std::pair<MyVector<Tint>, MyVector<T>> ePair =
      ReductionMod1vector<T, Tint>(cosetRed);
  T_shvec_request<T> request_ret{
      dim,          request.mode,    request.bound,
      ePair.second, eRec.GramMatRed, request.central};
  cvp_reduction_info<Tint> cvp_red{ePair.first, eRec.Pmat};
  return {std::move(request_ret), std::move(cvp_red)};
}

template <typename T, typename Tint>
T_shvec_info<T, Tint>
ApplyReductionToShvecInfo(T_shvec_info<T, Tint> const &info,
                          cvp_reduction_info<Tint> const &red_info) {
  std::vector<MyVector<Tint>> short_vectors;
  for (auto &z_vec : info.short_vectors) {
    MyVector<Tint> x = red_info.P.transpose() * (z_vec + red_info.shift);
    short_vectors.emplace_back(std::move(x));
  }
  return {std::move(short_vectors), info.minimum};
}

// We return floor(sqrt(A) + epsilon + B)
template <typename T> int Infinitesimal_Floor_V1(T const &a, T const &b) {
  double epsilon = 0.000000001;
#ifdef CHECK_BASIC_CONSISTENCY
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
#ifdef CHECK_BASIC_CONSISTENCY
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Ceil_V1\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  double a_doubl = UniversalScalarConversion<double, T>(a);
  double b_doubl = UniversalScalarConversion<double, T>(b);
  double alpha = -sqrt(a_doubl) - epsilon + b_doubl;
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
#ifdef CHECK_BASIC_CONSISTENCY
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
    T eDiff = T(x) - b;
    if (eDiff <= 0)
      return true;
    if (eDiff * eDiff <= a)
      return true;
    return false;
  };
  while (true) {
    bool test1 = f(eReturn);
    bool test2 = f(eReturn + 1);
    if (test1 && !test2)
      break;
    if (!test1)
      eReturn--;
    if (test2)
      eReturn++;
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
#ifdef CHECK_BASIC_CONSISTENCY
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
    T eDiff = x - b;
    if (eDiff >= 0)
      return true;
    if (eDiff * eDiff <= a)
      return true;
    return false;
  };
  while (true) {
    bool test1 = f(eReturn - 1);
    bool test2 = f(eReturn);
    if (!test1 && test2)
      break;
    if (test1)
      eReturn--;
    if (!test2)
      eReturn++;
  }
  return eReturn;
}

template <typename T, typename Tint, typename Finsert, typename Fsetbound>
int computeIt_Gen_Kernel(const T_shvec_request<T> &request, const T &bound,
                         Finsert f_insert, Fsetbound f_set_bound) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int i, j;
  int dim = request.dim;
  // The value of bound is assumed to be correct.
  // Thus the Trem values should be strictly positive.
  MyVector<Tint> Upper(dim);
  MyVector<T> Trem(dim);
  MyVector<T> U(dim);
  MyVector<Tint> x(dim);
#if defined CHECK_BASIC_CONSISTENCY || defined PRINT_DEBUG_INFO
  const MyMatrix<T> &g = request.gram_matrix;
#endif
#ifdef PRINT_DEBUG_INFO
  std::cerr << "g=\n";
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
#ifdef PRINT_DEBUG_INFO
    std::cerr << "diag q=" << q(i, i) << "\n";
    for (int j = i + 1; j < dim; j++)
      std::cerr << "   j=" << j << " q=" << q(i, j) << "\n";
#endif
  }
  const bool &central = request.central;
  const MyVector<T> &C = request.coset;
  bool needs_new_bound = true;
  i = dim - 1;
  if (bound < 0) {
    return TempShvec_globals::NORMAL_TERMINATION_COMPUTATION;
  }
  Trem(i) = bound;
  U(i) = 0;
#ifdef PRINT_DEBUG_INFO
  std::cerr << "Before while loop\n";
#endif
#ifdef PRINT_DEBUG_INFO_VECTOR
  size_t n_vector = 0;
#endif
  T eQuot, eSum, hVal, eNorm;
  while (true) {
    if (needs_new_bound) {
      eQuot = Trem(i) / q(i, i);
      eSum = -U(i) - C(i);
      f_set_bound(eQuot, eSum, q, x, i, Upper(i), x(i));
      needs_new_bound = false;
    } else {
      x(i) += 1;
    }
    if (x(i) <= Upper(i)) {
      if (i == 0) {
        if (central) {
          j = dim - 1;
          bool not_finished = false;
          while (j >= 0 && !not_finished) {
            not_finished = (x(j) != 0);
            j--;
          }
          if (!not_finished) {
#ifdef PRINT_DEBUG_INFO
            std::cerr << "Exiting because x=0 and central run\n";
#endif
            return TempShvec_globals::NORMAL_TERMINATION_COMPUTATION;
          }
        }
        hVal = x(0) + C(0) + U(0);
        eNorm = bound - Trem(0) + q(0, 0) * hVal * hVal;
#ifdef CHECK_BASIC_CONSISTENCY
        T norm = 0;
        for (int i2 = 0; i2 < dim; i2++)
          for (int j2 = 0; j2 < dim; j2++)
            norm += g(i2, j2) * (x(i2) + C(i2)) * (x(j2) + C(j2));
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
#ifdef PRINT_DEBUG_INFO_VECTOR
        std::cerr << "n_vector=" << n_vector;
        std::cerr << " x=";
        for (int i = 0; i < dim; i++)
          std::cerr << " " << x(i);
        std::cerr << "\n";
        n_vector++;
#endif
        bool ret_val = f_insert(x, eNorm);
        if (!ret_val)
          return TempShvec_globals::STOP_COMPUTATION;
      } else {
        i--;
        U(i) = 0;
        for (j = i + 1; j < dim; j++)
          U(i) += q(i, j) * (x(j) + C(j));
        hVal = x(i + 1) + C(i + 1) + U(i + 1);
        Trem(i) = Trem(i + 1) - q(i + 1, i + 1) * hVal * hVal;
        needs_new_bound = true;
      }
    } else {
      i++;
      if (i == dim) {
        return TempShvec_globals::NORMAL_TERMINATION_COMPUTATION;
      }
    }
  }
}

template <typename T, typename Tint, typename Finsert, typename Fsetbound>
inline typename std::enable_if<is_ring_field<T>::value, int>::type
computeIt_Gen(const T_shvec_request<T> &request, const T &bound,
              Finsert f_insert, Fsetbound f_set_bound) {
#ifdef PRINT_DEBUG_INFO
  std::cerr << "computeIt (field case)\n";
#endif
  return computeIt_Gen_Kernel<T, Tint, Finsert, Fsetbound>(
      request, bound, f_insert, f_set_bound);
}

template <typename T, typename Tint, typename Finsert, typename Fsetbound>
inline typename std::enable_if<!is_ring_field<T>::value, int>::type
computeIt_Gen(const T_shvec_request<T> &request, const T &bound,
              Finsert f_insert, Fsetbound f_set_bound) {
#ifdef PRINT_DEBUG_INFO
  std::cerr << "computeIt (ring case)\n";
#endif
  using Tfield = typename overlying_field<T>::field_type;
  //
  Tfield bound_field = UniversalScalarConversion<Tfield, T>(bound);
  T_shvec_request<Tfield> request_field{
      request.dim,
      request.mode,
      UniversalScalarConversion<Tfield, T>(request.bound),
      UniversalVectorConversion<Tfield, T>(request.coset),
      UniversalMatrixConversion<Tfield, T>(request.gram_matrix),
      request.central};
  //
  auto f_insert_field = [&](const MyVector<Tint> &V,
                            const Tfield &min_Tfield) -> bool {
    T min_T = UniversalScalarConversion<T, Tfield>(min_Tfield);
    return f_insert(V, min_T);
  };
  int retVal =
      computeIt_Gen_Kernel<Tfield, Tint, decltype(f_insert_field), Fsetbound>(
          request_field, bound_field, f_insert_field, f_set_bound);
#ifdef PRINT_DEBUG_INFO
  std::cerr << "computeIt (ring case) exit\n";
#endif
  return retVal;
}

template <typename T, typename Tint, typename Finsert>
inline typename std::enable_if<is_ring_field<T>::value, int>::type
computeIt(const T_shvec_request<T> &request, const T &bound, Finsert f_insert) {
  auto f_set_bound =
      [&](const T &eQuot, const T &eSum, [[maybe_unused]] const MyMatrix<T> &q,
          [[maybe_unused]] const MyVector<Tint> &x,
          [[maybe_unused]] const int &i, Tint &upper, Tint &lower) -> void {
    upper = Infinitesimal_Floor<T, Tint>(eQuot, eSum);
    lower = Infinitesimal_Ceil<T, Tint>(eQuot, eSum);
  };
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(
      request, bound, f_insert, f_set_bound);
}

template <typename T, typename Tint, typename Finsert>
inline typename std::enable_if<!is_ring_field<T>::value, int>::type
computeIt(const T_shvec_request<T> &request, const T &bound, Finsert f_insert) {
  using Tfield = typename overlying_field<T>::field_type;
  auto f_set_bound = [&](const Tfield &eQuot, const Tfield &eSum,
                         [[maybe_unused]] const MyMatrix<Tfield> &q,
                         [[maybe_unused]] const MyVector<Tint> &x,
                         [[maybe_unused]] const int &i, Tint &upper,
                         Tint &lower) -> void {
    upper = Infinitesimal_Floor<Tfield, Tint>(eQuot, eSum);
    lower = Infinitesimal_Ceil<Tfield, Tint>(eQuot, eSum);
  };
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(
      request, bound, f_insert, f_set_bound);
}

template <typename T, typename Tint, typename Finsert>
int computeIt_polytope(const T_shvec_request<T> &request, const T &bound,
                       const MyMatrix<T> &FAC, Finsert f_insert) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int n_rows = FAC.rows();
  int n_col = FAC.cols();
  if (n_col != request.dim + 1) {
    std::cerr << "Error in the size of FAC\n";
    throw TerminalException{1};
  }
  std::cerr << "Beginning of computeIt_polytope\n";
  auto f_set_bound = [&](const T &eQuot, const T &eSum,
                         [[maybe_unused]] const MyMatrix<T> &q,
                         const MyVector<Tint> &x, const int &i, Tint &upper,
                         Tint &lower) -> void {
    upper = Infinitesimal_Floor<T, Tint>(eQuot, eSum);
    lower = Infinitesimal_Ceil<T, Tint>(eQuot, eSum);
    int len = 2 + i;
    MyMatrix<T> FACwork(n_rows, len);
    for (int i_row = 0; i_row < n_rows; i_row++) {
      for (int i_col = 0; i_col < len; i_col++)
        FACwork(i_row, i_col) = FAC(i_row, i_col);
      for (int i_col = len; i_col < n_col; i_col++) {
        FACwork(i_row, 0) += x(i_col - 1) * FAC(i_row, i_col);
      }
    }
    LpSolution<T> eSol;
    MyVector<T> Vminimize = ZeroVector<T>(len);
    //
    Vminimize(1 + i) = 1;
    eSol = CDD_LinearProgramming(FACwork, Vminimize);
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
    eSol = CDD_LinearProgramming(FACwork, Vminimize);
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
  return computeIt_Gen<T, Tint, Finsert, decltype(f_set_bound)>(
      request, bound, f_insert, f_set_bound);
}

template <typename T, typename Tint>
T_shvec_info<T, Tint> computeMinimum(const T_shvec_request<T> &request) {
#ifdef PRINT_DEBUG_INFO
  std::cerr << "computeMinimum, begin\n";
#endif
  int i, j;
  int dim = request.dim;
  const MyVector<T> &C = request.coset;
  const bool &central = request.central;
  auto get_minimum_atp = [&]() -> T {
    if (!central) {
      T eNorm = 0;
      for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
          eNorm += request.gram_matrix(i, j) * C(i) * C(j);
      return eNorm;
    }
    T eMin = request.gram_matrix(0, 0);
    for (i = 1; i < dim; i++)
      if (eMin > request.gram_matrix(i, i))
        eMin = request.gram_matrix(i, i);
    return eMin;
  };
  T_shvec_info<T, Tint> info;
  info.minimum = get_minimum_atp();
  while (true) {
#ifdef PRINT_DEBUG_INFO
    std::cerr << "Before computeIt (in computeMinimum while loop)\n";
#endif
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == info.minimum) {
        info.short_vectors.push_back(V);
        if (central) {
          info.short_vectors.push_back(-V);
        }
        return true;
      } else {
        info.short_vectors.clear();
        info.minimum = min;
        return false;
      }
    };
    int result =
        computeIt<T, Tint, decltype(f_insert)>(request, info.minimum, f_insert);
    if (result == TempShvec_globals::NORMAL_TERMINATION_COMPUTATION) {
      break;
    }
  }
  return info;
}

template <typename T>
bool get_central(const MyVector<T> &coset, const int &mode) {
  int dim = coset.size();
  for (int i = 0; i < dim; i++) {
    if (coset(i) != 0)
      return false;
  }
  if (mode == TempShvec_globals::TEMP_SHVEC_MODE_VINBERG_ALGO)
    return false;
  return true;
}

template <typename T>
T_shvec_request<T> initShvecReq(const MyMatrix<T> &gram_matrix,
                                const MyVector<T> &coset, const T &bound,
                                int mode) {
  int dim = gram_matrix.rows();
  /*
  if (dim < 2) {
    std::cerr << "dim=" << dim << " while it should be at least 2\n";
    throw TerminalException{1};
  }
  */
  T_shvec_request<T> request;
  request.dim = dim;
  request.coset = coset;
  request.gram_matrix = gram_matrix;
  request.mode = mode;
  request.bound = bound;
  request.central = get_central(coset, mode);
  return request;
}

template <typename T, typename Tint>
T_shvec_info<T, Tint> T_computeShvec_Kernel(const T_shvec_request<T> &request) {
  if (request.mode == TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS) {
    return computeMinimum<T, Tint>(request);
  }
  if (request.mode == TempShvec_globals::TEMP_SHVEC_MODE_MINIMUM) {
    return computeMinimum<T, Tint>(request);
  }
  T_shvec_info<T, Tint> info;
  if (request.mode == TempShvec_globals::TEMP_SHVEC_MODE_BOUND) {
    info.minimum = request.bound;
    auto f_insert = [&](const MyVector<Tint> &V,
                        [[maybe_unused]] const T &min) -> bool {
      info.short_vectors.push_back(V);
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, request.bound,
                                                 f_insert);
    return info;
  }
  if (request.mode == TempShvec_globals::TEMP_SHVEC_MODE_HAN_TRAN) {
    info.minimum = request.bound;
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == request.bound) {
        info.short_vectors.push_back(V);
        return false;
      }
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, request.bound,
                                                 f_insert);
    return info;
  }
  if (request.mode == TempShvec_globals::TEMP_SHVEC_MODE_VINBERG_ALGO) {
    info.minimum = request.bound;
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == request.bound)
        info.short_vectors.push_back(V);
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, request.bound,
                                                 f_insert);
    return info;
  }
  if (request.mode == TempShvec_globals::TEMP_SHVEC_MODE_LORENTZIAN) {
    info.minimum = request.bound;
    auto f_insert = [&](const MyVector<Tint> &V,
                        [[maybe_unused]] const T &min) -> bool {
      info.short_vectors.push_back(V);
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, request.bound,
                                                 f_insert);
    return info;
  }
  std::cerr << "mode=" << request.mode << "\n";
  std::cerr << "T_compiteShvec: Failed to match an entry\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
T_shvec_info<T, Tint> T_computeShvec(const T_shvec_request<T> &request) {
  std::pair<T_shvec_request<T>, cvp_reduction_info<Tint>> ePair =
      GetReducedShvecRequest<T, Tint>(request);
  T_shvec_info<T, Tint> info1 = T_computeShvec_Kernel<T, Tint>(ePair.first);
  T_shvec_info<T, Tint> info2 = ApplyReductionToShvecInfo(info1, ePair.second);
  return info2;
}

template <typename T, typename Tint>
resultCVP<T, Tint> CVPVallentinProgram_exact(MyMatrix<T> const &GramMat,
                                             MyVector<T> const &eV) {
  int dim = GramMat.rows();
  MyVector<T> cosetVect = -eV;
  if (IsIntegralVector(eV)) {
    T TheNorm = 0;
    MyMatrix<Tint> ListVect(1, dim);
    for (int i = 0; i < dim; i++)
      ListVect(0, i) = UniversalScalarConversion<Tint, T>(eV(i));
    return {TheNorm, ListVect};
  }
  // Following value of bound = 0 should not be used
  T bound = 0;
  int mode = TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS;
  T_shvec_request<T> request = initShvecReq(GramMat, cosetVect, bound, mode);
  T_shvec_info<T, Tint> info = T_computeShvec<T, Tint>(request);
  int nbVect = info.short_vectors.size();
  MyMatrix<Tint> ListClos(nbVect, dim);
  for (int iVect = 0; iVect < nbVect; iVect++)
    for (int i = 0; i < dim; i++)
      ListClos(iVect, i) = info.short_vectors[iVect](i);
  MyVector<T> eDiff(dim);
  for (int i = 0; i < dim; i++)
    eDiff(i) = ListClos(0, i) - eV(i);
  T TheNorm = EvaluationQuadForm<T, T>(GramMat, eDiff);
  return {TheNorm, std::move(ListClos)};
}

template <typename T, typename Tint>
MyMatrix<Tint> T_ShortVector_exact(MyMatrix<T> const &GramMat,
                                   T const &MaxNorm) {
  int dim = GramMat.rows();
  if (dim == 1) {
    std::vector<MyVector<Tint>> ListVect;
    int idx = 1;
    while (true) {
      T norm = idx * idx * GramMat(0, 0);
      if (norm > MaxNorm)
        break;
      MyVector<Tint> eVect1(1);
      eVect1(0) = idx;
      ListVect.push_back(eVect1);
      idx++;
    }
    return MatrixFromVectorFamily(ListVect);
  }
  T bound = MaxNorm;
  int mode = TempShvec_globals::TEMP_SHVEC_MODE_BOUND;
  MyVector<T> cosetVect = ZeroVector<T>(dim);
  T_shvec_request<T> request = initShvecReq(GramMat, cosetVect, bound, mode);
  request.central = true;
  //
  T_shvec_info<T, Tint> info = T_computeShvec<T, Tint>(request);
  //
  return MatrixFromVectorFamily(info.short_vectors);
}

template <typename T, typename Tint>
MyMatrix<Tint> T_ShortVector_fixed(MyMatrix<T> const &GramMat,
                                   T const &SpecNorm) {
  int dim = GramMat.rows();
  std::vector<MyVector<Tint>> ListVect;
  if (dim == 1) {
    int idx = 1;
    while (true) {
      T norm = idx * idx * GramMat(0, 0);
      if (norm == SpecNorm) {
        MyVector<Tint> V(1);
        V(0) = idx;
        ListVect.push_back(V);
        break;
      }
      if (norm > SpecNorm)
        break;
      idx++;
    }
  } else {
    int mode = TempShvec_globals::TEMP_SHVEC_MODE_BOUND;
    MyVector<T> cosetVect = ZeroVector<T>(dim);
    T_shvec_request<T> request =
        initShvecReq(GramMat, cosetVect, SpecNorm, mode);
    request.central = true;
    //
    auto f_insert = [&](const MyVector<Tint> &V, const T &min) -> bool {
      if (min == SpecNorm) {
        ListVect.push_back(V);
      }
      return true;
    };
    (void)computeIt<T, Tint, decltype(f_insert)>(request, SpecNorm, f_insert);
  }
  return MatrixFromVectorFamilyDim(dim, ListVect);
}

// clang-format off
#endif  // SRC_LATT_SHVEC_EXACT_H_
// clang-format on
