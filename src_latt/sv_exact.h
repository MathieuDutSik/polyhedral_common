// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SV_EXACT_H_
#define SRC_LATT_SV_EXACT_H_

// clang-format off
#include "Shvec_exact.h"
// clang-format on

template <typename Tint> struct cvp_reduction_info {
  MyVector<Tint> shift;
  MyMatrix<Tint> P;
};

namespace TempShvec_globals {
const int TEMP_SHVEC_MODE_UNDEF = -1;
const int TEMP_SHVEC_MODE_BOUND = 0;
const int TEMP_SHVEC_MODE_SHORTEST_VECTORS = 1;
const int TEMP_SHVEC_MODE_MINIMUM = 2;
const int TEMP_SHVEC_MODE_THETA_SERIES = 3;
const int TEMP_SHVEC_MODE_VINBERG_ALGO = 4;
const int TEMP_SHVEC_MODE_LORENTZIAN = 5;
const int TEMP_SHVEC_MODE_HAN_TRAN = 6;
// clang-format off
}  // namespace TempShvec_globals
// clang-format on

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
FullGramInfo<T> initShvecReq(const MyMatrix<T> &gram_matrix,
                                const MyVector<T> &coset, const T &bound,
                                int mode) {
  int dim = gram_matrix.rows();
  FullGramInfo<T> request;
  request.dim = dim;
  request.coset = coset;
  request.gram_matrix = gram_matrix;
  request.bound = bound;
  request.central = get_central(coset, mode);
  return request;
}

/*
  We apply the P reduction
  Gred = P G P^T
  Q = P^{-1}
  G = Q Gred Q^T
  looking for G[x + coset] = bound
  That is (x + coset) Q Gred Q^T (x + coset)^T = bound
  Or (y + coset_red) Gred (y + coset_red)^T = bound
  with coset_red = coset Q and y = x Q
  coset_red = V_T + V_i
  Thus we have if z = y + V_i
  (z + V_T) Gred (z + v_T)^T = bound
  So y = z - V_i
  and x = (z - V_i) P
 */
template <typename T, typename Tint>
std::pair<FullGramInfo<T>, cvp_reduction_info<Tint>>
GetReducedShvecRequest(FullGramInfo<T> const &request, std::ostream& os) {
  int dim = request.dim;
  LLLreduction<T, Tint> eRec =
    LLLreducedBasisDual<T, Tint>(request.gram_matrix, os);
  MyMatrix<Tint> Q_i = Inverse(eRec.Pmat);
  MyMatrix<T> Q_T = UniversalMatrixConversion<T, Tint>(Q_i);
  MyVector<T> cosetRed = Q_T.transpose() * request.coset;
  std::pair<MyVector<Tint>, MyVector<T>> ePair =
      ReductionMod1vector<T, Tint>(cosetRed);
  FullGramInfo<T> request_ret{dim, request.bound, ePair.second,
                                 eRec.GramMatRed, request.central};
  cvp_reduction_info<Tint> cvp_red{ePair.first, eRec.Pmat};
  return {std::move(request_ret), std::move(cvp_red)};
}

template <typename T, typename Tint>
T_shvec_info<T, Tint>
ApplyReductionToShvecInfo(T_shvec_info<T, Tint> const &info,
                          cvp_reduction_info<Tint> const &red_info) {
  std::vector<MyVector<Tint>> short_vectors;
  for (auto &z_vec : info.short_vectors) {
    MyVector<Tint> x = red_info.P.transpose() * (z_vec - red_info.shift);
    short_vectors.emplace_back(std::move(x));
  }
  return {std::move(short_vectors), info.minimum};
}

template <typename T, typename Tint>
T_shvec_info<T, Tint> T_computeShvec_Kernel(const FullGramInfo<T> &request,
                                            int mode) {
  if (mode == TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS) {
    return computeMinimum<T, Tint>(request);
  }
  if (mode == TempShvec_globals::TEMP_SHVEC_MODE_MINIMUM) {
    return computeMinimum<T, Tint>(request);
  }
  T_shvec_info<T, Tint> info;
  if (mode == TempShvec_globals::TEMP_SHVEC_MODE_BOUND) {
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
  if (mode == TempShvec_globals::TEMP_SHVEC_MODE_HAN_TRAN) {
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
  if (mode == TempShvec_globals::TEMP_SHVEC_MODE_VINBERG_ALGO) {
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
  if (mode == TempShvec_globals::TEMP_SHVEC_MODE_LORENTZIAN) {
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
  std::cerr << "mode=" << mode << "\n";
  std::cerr << "T_computeShvec: Failed to match an entry\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
T_shvec_info<T, Tint> T_computeShvec(const FullGramInfo<T> &request,
                                     int mode,
                                     std::ostream &os) {
#ifdef TIMINGS_SHVEC
  MicrosecondTime time;
#endif
  std::pair<FullGramInfo<T>, cvp_reduction_info<Tint>> ePair =
    GetReducedShvecRequest<T, Tint>(request, os);
#ifdef TIMINGS_SHVEC
  os << "|SHVEC: GetReducedShvecRequest|=" << time << "\n";
#endif
  T_shvec_info<T, Tint> info1 =
      T_computeShvec_Kernel<T, Tint>(ePair.first, mode);
#ifdef TIMINGS_SHVEC
  os << "|SHVEC: T_computeShvec_Kernel|=" << time << "\n";
#endif
  T_shvec_info<T, Tint> info2 = ApplyReductionToShvecInfo(info1, ePair.second);
#ifdef TIMINGS_SHVEC
  os << "|SHVEC: ApplyReductionToShvecInfo|=" << time << "\n";
#endif
  return info2;
}

// clang-format off
#endif  // SRC_LATT_SV_EXACT_H_
// clang-format on
