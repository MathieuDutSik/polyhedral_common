// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_INDEFAPPROXCANINICAL_H_
#define SRC_INDEFINITE_INDEFAPPROXCANINICAL_H_

// clang-format off
#include "Indefinite_LLL.h"
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "WeightMatrix.h"
// clang-format on

#ifdef DEBUG
#define DEBUG_INDEX_APPROX_CANONICAL
#endif

template <typename T, typename Tint>
ResultReduction<T, Tint> CanonicalizationPermutationSigns(MyMatrix<T> const &M,
                                                          std::ostream &os) {
  using Tidx_value = uint16_t;
  using Tidx = uint16_t;
  using Tgr = GraphListAdj;
  Tidx n = M.rows();
  if (n <= 1) {
    MyMatrix<Tint> B = IdentityMat<Tint>(n);
    MyMatrix<T> Mred = M;
    return {B, Mred};
  }
  MyMatrix<T> Mabs(n, n);
  for (Tidx i_row = 0; i_row < n; i_row++) {
    for (Tidx i_col = 0; i_col < n; i_col++) {
      Mabs(i_row, i_col) = T_abs(M(i_row, i_col));
    }
  }
  WeightMatrix<true, T, Tidx_value> WMat =
      WeightedMatrixFromMyMatrix<true, T, Tidx_value>(Mabs, os);
  WMat.ReorderingSetWeight();
  std::vector<Tidx> CanonicOrd =
      GetGroupCanonicalizationVector_Kernel<T, Tgr, Tidx, Tidx_value>(WMat, os)
          .first;
  MyMatrix<T> Mreord(n, n);
  MyMatrix<T> MreordAbs(n, n);
  for (Tidx i_row = 0; i_row < n; i_row++) {
    Tidx j_row = CanonicOrd[i_row];
    for (Tidx i_col = 0; i_col < n; i_col++) {
      Tidx j_col = CanonicOrd[i_col];
      Mreord(i_row, i_col) = M(j_row, j_col);
      MreordAbs(i_row, i_col) = T_abs(M(j_row, j_col));
    }
  }
  MyMatrix<Tint> Mtrans1 = ZeroMatrix<Tint>(n, n);
  for (Tidx i_row = 0; i_row < n; i_row++) {
    Tidx j_row = CanonicOrd[i_row];
    Mtrans1(i_row, j_row) = 1;
  }
#ifdef DEBUG_INDEFINITE_LLL
  MyMatrix<T> Mtrans1_T = UniversalMatrixConversion<T, Tint>(Mtrans1);
  MyMatrix<T> eProd = Mtrans1_T * M * Mtrans1_T.transpose();
  if (eProd != Mreord) {
    std::cerr << "The matrix product does not work as expected\n";
    std::cerr << "eProd=\n";
    WriteMatrix(std::cerr, eProd);
    throw TerminalException{1};
  }
#endif
  GraphBitset GR(n);
  for (Tidx i_row = 0; i_row < n; i_row++) {
    for (Tidx i_col = 0; i_col < n; i_col++) {
      if (Mreord(i_row, i_col) != 0) {
        GR.AddAdjacent(i_row, i_col);
        GR.AddAdjacent(i_col, i_row);
      }
    }
  }
  MyMatrix<Tint> Mtrans2 = IdentityMat<Tint>(n);
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(GR);
  for (auto &eConn : LConn) {
    size_t len = eConn.size();
    std::vector<size_t> Status(len, 0);
    std::vector<size_t> eConnRev(n, std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < len; i++)
      eConnRev[eConn[i]] = i;
    Status[0] = 1;
    while (true) {
      size_t n_done = 0;
      for (size_t i = 0; i < len; i++) {
        if (Status[i] > 0) {
          n_done++;
        }
      }
      if (n_done == len) {
        break;
      }
      for (size_t i = 0; i < len; i++) {
        if (Status[i] > 0) {
          size_t iImg = eConn[i];
          for (auto &eAdjImg : GR.Adjacency(iImg)) {
            size_t eAdj = eConnRev[eAdjImg];
            if (Status[eAdj] == 0) {
              int sign = 1;
              if (Mreord(iImg, eAdjImg) < 0)
                sign = -1;
              Mtrans2(eAdjImg, eAdjImg) = sign * Mtrans2(iImg, iImg);
              Status[eAdj] = 1;
            }
          }
        }
      }
    }
  }
  MyMatrix<Tint> eP = Mtrans2 * Mtrans1;
  MyMatrix<T> eP_T = UniversalMatrixConversion<T, Tint>(eP);
  MyMatrix<T> M_red = eP_T * M * eP_T.transpose();
  return {std::move(eP), std::move(M_red)};
}

template<typename T>
bool IsEichlerTwoBlock(MyMatrix<T> const& M) {
  if (M.rows() != 2) {
    return false;
  }
  for (int i=0; i<2; i++) {
    if (M(i,i) != 0) {
      return false;
    }
  }
  return true;
}

template<typename T>
struct OrderingInfoIndefForm {
  bool is_eichler_two;
  int dim;
  T sum_abs_coeff;
  T max_abs_coeff;
};

template<typename T>
OrderingInfoIndefForm<T> get_ordering_info_indef_form(MyMatrix<T> const& M) {
  bool is_eichler_two = IsEichlerTwoBlock(M);
  int dim = M.rows();
  T sum_abs_coeff(0);
  T max_abs_coeff(0);
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      T val = T_abs(M(i,j));
      sum_abs_coeff += val;
      if (val > max_abs_coeff) {
        max_abs_coeff = val;
      }
    }
  }
  return {is_eichler_two, dim, sum_abs_coeff, max_abs_coeff};
}

/*
  We order the blocks in the way we want it to be.
  The first value is the result of the comparison.
  The second is whether the ordering info are equal or not.
 */
template<typename T>
std::pair<bool, bool> compare_ordering_info(OrderingInfoIndefForm<T> const& ord1, OrderingInfoIndefForm<T> const& ord2) {
  // First Eichler Block of 2x2 are preferred.
  if (ord1.is_eichler_two && !ord2.is_eichler_two) {
    return {true, false};
  }
  if (!ord1.is_eichler_two && ord2.is_eichler_two) {
    return {false, false};
  }
  // Then the dimension is used for the comparison
  if (ord1.dim < ord2.dim) {
    return {true, false};
  }
  if (ord2.dim > ord2.dim) {
    return {false, false};
  }
  // Then the sum of the absolute coefficients
  if (ord1.sum_abs_coeff < ord2.sum_abs_coeff) {
    return {true, false};
  }
  if (ord2.sum_abs_coeff > ord2.sum_abs_coeff) {
    return {false, false};
  }
  // Then the maximum of the absolute value of the coefficients
  if (ord1.max_abs_coeff < ord2.max_abs_coeff) {
    return {true, false};
  }
  if (ord2.max_abs_coeff > ord2.max_abs_coeff) {
    return {false, false};
  }
  // Nothing works, reporting as such.
  return {true, true};
}

template<typename T, typename Tint>
ResultReduction<T, Tint> order_blocks_by_signature(MyMatrix<T> const& M, [[maybe_unused]] std::ostream& os) {
  int n = M.rows();
  std::vector<std::vector<size_t>> LConn = MatrixConnectedComponents(M);
  size_t n_conn = LConn.size();
  std::vector<OrderingInfoIndefForm<T>> l_ord_infos;
  std::vector<MyMatrix<T>> l_block;
  for (auto & eConn: LConn) {
    int len = eConn.size();
    MyMatrix<T> M_conn(len, len);
    for (int i=0; i<len; i++) {
      for (int j=0; j<len; j++) {
        size_t i_big = eConn[i];
        size_t j_big = eConn[j];
        M_conn(i, j) = Min(i_big, j_big);
      }
    }
    l_block.push_back(M_conn);
    OrderingInfoIndefForm<T> eInfo = get_ordering_info_indef_form(Mconn);
    l_ord_infos.push_back(eInfo);
  }
  std::vector<size_t> ListIdx(n_conn);
  for (size_t iConn=0; iConn<n_conn; i_conn++) {
    ListIdx[iConn] = iConn;
  }
  //
#ifdef DEBUG_INDEX_APPROX_CANONICAL
  size_t n_false_equality = 0;
#endif
  std::stable_sort(ListIdx.begin(), ListIdx.end(), [&](size_t idx1, size_t idx2) -> bool {
    std::pair<bool, bool> pair = compare_ordering_info(l_ord_infos[idx1], l_ord_infos[idx2]);
    if (!pair.second) {
      // if the comparison works then use it.
      return pair.first;
    }
#ifdef DEBUG_INDEX_APPROX_CANONICAL
    if (l_block[idx1] != l_block[idx1]) {
      n_false_equality += 1;
    }
#endif
    return false;
  });
#ifdef DEBUG_INDEX_APPROX_CANONICAL
  os << "IAC: order_blocks_by_signature, n_false_equality=" << n_false_equality << "\n";
#endif
  MyMatrix<T> Mret = ZeroMatrix<T>(n, n);
  MyMatrix<Tint> Bret = ZeroMatrix<Tint>(n, n);
  size_t pos = 0;
  for (auto & idx : ListIdx) {
    std::vector<size_t> const& eConn = LConn[idx];
    MyMatrix<T> const& M_block = l_block[idx];
    int len = eConn.size();
    for (int i=0; i<len; i++) {
      for (int j=0; j<len; j++) {
        Mret(pos + i, pos + j) = M_block(i, j);
      }
    }
    for (int i=0; i<len; i++) {
      int j = eConn[i];
      Bret(pos + i, j) = 1;
    }
    pos += len;
  }
#ifdef DEBUG_INDEX_APPROX_CANONICAL
  MyMatrix<T> Bret_T = UniversalMatrixConversion<T,Tint>(Bret);
  MyMatrix<T> prod = Bret_T * M * Bret_T.transpose();
  if (prod != Mret) {
    std::cerr << "The matrix is not an equivalence\n";
    throw TerminalException{1};
  }
#endif
  return {Bret, Mret};
}


template <typename T, typename Tint>
ResultReduction<T, Tint>
get_individual_reduction(MyMatrix<T> const& M, std::ostream &os) {
  std::pair<int, int> signature = GetSignature(M);
  int n_plus = signature.first;
  int n_minus = signature.second;
  int n_zero = n - n_plus - n_minus;
  if (n_zero > 0) {
#ifdef DEBUG_INDEX_APPROX_CANONICAL
    if (n_plus > 0 || n_minus > 0) {
      std::cerr << "We should have a zero matrix\n";
      throw TerminalException{1};
    }
    if (n_zero != 1) {
      std::cerr << "The dimension of the zero matrix should be 1\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<Tint> B = IdentityMat<Tint>(1);
    return {B, M};
  }
  if (n_plus > 0 && n_minus > 0) {
    return CanonicalizationPermutationSigns<T, Tint>(M, os);
  }
  if (n_plus == 0 || n_minus == 0) {
    auto get_posdef=[&]() -> std::pair<int, MyMatrix<T>> {
      if (n_minus > 0) {
        MyMatrix<T> ret = -M;
        return {-1, ret};
      }
      return {1, M};
    };
    std::pair<int, MyMatrix<T>> pair = get_posdef();
    MyMatrix<T> const& Mpos = pair.second;
    Canonic_PosDef<T, Tint> cpd = ComputeCanonicalForm<T, Tint>(Mpos, os);
    MyMatrix<Tint> const& B = cpd.Basis;
    MyMatrix<T> Mred = T(pair.first) * cpd.Mat;
#ifdef DEBUG_INDEX_APPROX_CANONICAL
    MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
    MyMatrix<T> prod = B_T * Mred * B_T.transpose();
    if (prod != M) {
      std::cerr << "The reduction did not work out correctly\n";
      throw TerminalException{1};
    }
#endif
    return {B, Mred};
  }
  std::cerr << "IAC: Could not find the right reduction method\n";
  throw TerminalException{1};
}


template <typename T, typename Tint>
ResultReduction<T, Tint>
apply_reduction_on_blocks(MyMatrix<T> const& M, std::ostream &os) {
  std::vector<std::vector<size_t>> LConn = MatrixConnectedComponents(M);
  MyMatrix<T> Mred = ZeroMatrix<T>(n,n);
  MyMatrix<Tint> B = ZeroMatrix<Tint>(n,n);
  for (auto & eConn : LConn) {
    int len = eConn.size();
    MyMatrix<T> M_block(len,len);
    for (int i=0; i<len; i++) {
      for (int j=0; j<len; j++) {
        int i_big = eConn[i];
        int j_big = eConn[j];
        M_block(i, j) = M(i_big, j_big);
      }
    }
    ResultReduction<T, Tint> res = get_individual_reduction<T,Tint>(M_block, os);
    for (int i=0; i<len; i++) {
      for (int j=0; j<len; j++) {
        int i_big = eConn[i];
        int j_big = eConn[j];
        Mred(i_big, j_big) = res.Mred(i, j);
        B(i_big, j_big) = res.B(i, j);
      }
    }
  }
  return {B, Mred};
}


/*
  This code attempts to find a canonical form for a form.
  We only require the input matrix to be symmetric.
 */
template <typename T, typename Tint>
ResultReduction<T, Tint>
ApproxCanonicalIndefiniteForm(MyMatrix<T> const &M, std::ostream &os) {
  ResultReduction<T, Tint> RRI_A = IndefiniteReduction<T,Tint>(M, os);
  ResultReduction<T, Tint> RRI_B = apply_reduction_on_blocks<T,Tint>(RRI_A.Mred, os);
  ResultReduction<T, Tint> RRI_C = order_blocks_by_signature<T,Tint>(RRI_B.Mred, os);
  MyMatrix<Tint> B = RRI_C.B * RRI_B.B * RRI_A.B;
#ifdef DEBUG_INDEX_APPROX_CANONICAL
  MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
  MyMatrix<T> prod = B_T * M * B_T.transpose();
  if (prod != RRI_C.Mred) {
    std::cerr << "IAC: The matrix B of the reduction did not work out well\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(B), std::move(RRI_B.Mred)};
}

// clang-format off
#endif  //  SRC_INDEFINITE_INDEFAPPROXCANINICAL_H_
// clang-format on
