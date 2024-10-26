// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_INDEFAPPROXCANINICAL_H_
#define SRC_INDEFINITE_INDEFAPPROXCANINICAL_H_

// clang-format off
#include "Indefinite_LLL.h"
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "WeightMatrix.h"
// clang-format on

template <typename T, typename Tint>
ResultReduction<T, Tint> CanonicalizationPermutationSigns(MyMatrix<T> const &M,
                                                          std::ostream &os) {
  using Tidx_value = uint16_t;
  using Tidx = uint16_t;
  using Tgr = GraphListAdj;
  Tidx n = M.rows();
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

template <typename T, typename Tint>
ResultReduction<T, Tint>
ComputeReductionIndefinitePermSign(MyMatrix<T> const &M, std::ostream &os) {
  if (M.rows() == 1) {
    MyMatrix<Tint> eP = IdentityMat<Tint>(1);
    return {std::move(eP), M};
  }
  ResultIndefiniteLLL<T, Tint> RRI_A =
      ComputeReductionIndefinite<T, Tint>(M, os);
  ResultReduction<T, Tint> RRI_B =
      CanonicalizationPermutationSigns<T, Tint>(RRI_A.Mred, os);
  MyMatrix<Tint> eP = RRI_B.B * RRI_A.B;
  return {std::move(eP), std::move(RRI_B.Mred)};
}

// clang-format off
#endif  //  SRC_INDEFINITE_INDEFAPPROXCANINICAL_H_
// clang-format on
