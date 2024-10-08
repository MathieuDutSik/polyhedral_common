// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_INDEFINITE_LLL_H_
#define SRC_INDEFINITE_INDEFINITE_LLL_H_

#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "MAT_Matrix.h"
#include "WeightMatrix.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#ifdef DEBUG
#define DEBUG_INDEFINITE_LLL
#endif

#ifdef TIMINGS
#define TIMINGS_INDEFINITE_LLL
#endif

template <typename T> struct ResultGramSchmidt_Indefinite {
  // true means we have a basis. False that we have an isotropic vector
  bool success;
  std::vector<MyVector<T>> Bstar;
  std::vector<T> Bstar_norms;
  MyMatrix<T> mu;
  MyVector<T> Xisotrop;
};

template <typename T, typename Tint>
ResultGramSchmidt_Indefinite<T>
GramSchmidtOrthonormalization(MyMatrix<T> const &M, MyMatrix<Tint> const &B) {
  int n = M.rows();
  MyMatrix<T> mu(n, n);
  struct PairInf {
    MyVector<T> Bistar_M;
    T Bistar_norm;
  };
  std::vector<PairInf> l_inf;
  std::vector<MyVector<T>> Bstar;
  std::vector<T> Bstar_norms;
#ifdef DEBUG_INDEFINITE_LLL
  T det1 = DeterminantMat(M);
  T det2 = 1;
#endif
  for (int i = 0; i < n; i++) {
    MyVector<T> Bistar = UniversalVectorConversion<T, Tint>(GetMatrixRow(B, i));
    for (int j = 0; j < i; j++) {
      T muij = (Bistar.dot(l_inf[j].Bistar_M)) / (l_inf[j].Bistar_norm);
      mu(i, j) = muij;
      Bistar -= muij * Bstar[j];
    }
    MyVector<T> Bistar_M = M * Bistar;
    T scal = Bistar_M.dot(Bistar);
    if (scal == 0) {
      return {false, {}, {}, {}, Bistar};
    }
#ifdef DEBUG_INDEFINITE_LLL
    det2 = det2 * scal;
#endif
    PairInf epair{std::move(Bistar_M), scal};
    Bstar.push_back(Bistar);
    Bstar_norms.push_back(scal);
    l_inf.push_back(epair);
  }
#ifdef DEBUG_INDEFINITE_LLL
  std::cerr << "ILLL: det1=" << det1 << " det2=" << det2 << "\n";
#endif
  return {true, Bstar, Bstar_norms, mu, {}};
}

template <typename T, typename Tint> struct ResultIndefiniteLLL {
  // true if we obtained the reduced matrix. false if we found an
  // isotropic vector
  MyMatrix<Tint> B;
  MyMatrix<T> Mred;
  std::optional<MyVector<T>> Xisotrop;
};

// Adapted from Denis Simon, Solving Quadratic Equations Using Reduced
// Unimodular Quadratic Forms, Math. Comp. 74(251) 1531--1543
template <typename T, typename Tint>
ResultIndefiniteLLL<T, Tint> Indefinite_LLL(MyMatrix<T> const &M) {
  int n = M.rows();
  // The c constant of the LLL algorithm
  T c = T(7) / T(8);
  MyMatrix<Tint> B = IdentityMat<Tint>(n);
  auto get_matrix = [&]() -> MyMatrix<T> {
    MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
    MyMatrix<T> Mred = B_T * M * B_T.transpose();
    return Mred;
  };
#ifdef DEBUG_INDEFINITE_LLL
  T det = DeterminantMat(M);
#endif
  int k = 1;
  while (true) {
#ifdef DEBUG_INDEFINITE_LLL
    std::cerr << "ILLL: Passing in Indefinite_LLL det=" << det << " k=" << k
              << "\n";
#endif
    ResultGramSchmidt_Indefinite<T> ResGS = GramSchmidtOrthonormalization(M, B);
#ifdef DEBUG_INDEFINITE_LLL
    std::cerr << "ILLL: Bstar_norms =";
    for (auto &eN : ResGS.Bstar_norms)
      std::cerr << " " << eN;
    std::cerr << "\n";
#endif
    if (!ResGS.success) {
      return {B, get_matrix(), ResGS.Xisotrop};
    }
    for (int i = n - 1; i >= 0; i--) {
      for (int j = 0; j < i; j++) {
        T val = ResGS.mu(i, j);
        Tint q = UniversalNearestScalarInteger<Tint, T>(val);
        B.row(i) -= q * B.row(j);
      }
    }
    T mu = ResGS.mu(k, k - 1);
    T sum1_pre = ResGS.Bstar_norms[k] + mu * mu * ResGS.Bstar_norms[k - 1];
    T sum1 = T_abs(sum1_pre);
    T sum2 = c * T_abs(ResGS.Bstar_norms[k - 1]);
    if (sum1 < sum2) {
#ifdef DEBUG_INDEFINITE_LLL
      std::cerr << "ILLL: Swapping k=" << k << " and " << (k - 1) << "\n";
#endif
      for (int i = 0; i < n; i++) {
        std::swap(B(k, i), B(k - 1, i));
      }
      k = std::max(k - 1, 1);
    } else {
      k++;
    }
    if (k >= n)
      break;
  }
  return {B, get_matrix(), {}};
}

template <typename Tint> MyMatrix<Tint> get_random_int_matrix(int const &n) {
  std::vector<int> LPos(n);
  for (int i = 0; i < n; i++)
    LPos[i] = i;
  for (int iter = 0; iter < 4 * n; iter++) {
    int i = random() % n;
    int j = random() % n;
    if (i != j)
      std::swap(LPos[i], LPos[j]);
  }
  std::vector<int> LDiag(n);
  for (int i = 0; i < n; i++) {
    int val = random() % 2;
    LDiag[i] = -1 + 2 * val;
  }
  MyMatrix<Tint> Unit = ZeroMatrix<Tint>(n, n);
  for (int i = 0; i < n; i++)
    Unit(i, LPos[i]) = LDiag[i];
  return Unit;
}

template <typename T, typename Tint>
ResultIndefiniteLLL<T, Tint>
ComputeReductionIndefinite(MyMatrix<T> const &M,
                           [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_INDEFINITE_LLL
  MicrosecondTime time;
#endif
#ifdef DEBUG_INDEFINITE_LLL
  os << "ILLL: Beginning of ComputeReductionIndefinite\n";
#endif
  int n = M.rows();
  auto get_norm = [&](MyMatrix<T> const &mat) -> T {
    T sum = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        sum += T_abs(mat(i, j));
    return sum;
  };
  MyMatrix<Tint> B = IdentityMat<Tint>(n);
  MyMatrix<T> B_T = IdentityMat<T>(n);
  MyMatrix<T> Mwork = M;
  T norm_work = get_norm(Mwork);
  size_t iter_no_improv = 0;
  size_t limit_iter = 2 * n;
  while (true) {
    ResultIndefiniteLLL<T, Tint> res = Indefinite_LLL<T, Tint>(Mwork);
#ifdef DEBUG_INDEFINITE_LLL
    os << "ILLL: We have computed res\n";
#endif
    // Terminating if we find an isotropic vector
    if (res.Xisotrop) {
      MyVector<T> const &Xisotrop = *res.Xisotrop;
      MyVector<T> V = B_T.transpose() * Xisotrop;
#ifdef DEBUG_INDEFINITE_LLL
      T sum = EvaluationQuadForm<T, T>(M, V);
      if (sum != 0) {
        std::cerr << "ILLL: Error in ComputeReductionIndefinite, M=\n";
        WriteMatrix(std::cerr, M);
        std::cerr << "ILLL: V should be an isotropic vector\n";
        throw TerminalException{1};
      }
#endif
      MyMatrix<Tint> Bret = res.B * B;
      Mwork = res.Mred;
#ifdef DEBUG_INDEFINITE_LLL
      MyMatrix<T> Bret_T = UniversalMatrixConversion<T, Tint>(Bret);
      MyMatrix<T> prod = Bret_T * M * Bret_T.transpose();
      if (prod != Mwork) {
        std::cerr << "ILLL: Bret is not the correct reduction matrix\n";
        throw TerminalException{1};
      }
#endif
#ifdef TIMINGS_INDEFINITE_LLL
      os << "ILLL: |ComputeReductionIndefinite(Iso)|=" << time << "\n";
#endif
      return {Bret, Mwork, V};
    }
    // Applying the reduction
    B = res.B * B;
    B_T = UniversalMatrixConversion<T, Tint>(B);
    Mwork = res.Mred;
    T norm = get_norm(res.Mred);
#ifdef DEBUG_INDEFINITE_LLL
    MyMatrix<T> prod = B_T * M * B_T.transpose();
    if (prod != Mwork) {
      std::cerr << "ILLL: B is not the correct reduction matrix\n";
      throw TerminalException{1};
    }
#endif
    // Updating metric
    if (norm >= norm_work) {
      iter_no_improv++;
      if (limit_iter == iter_no_improv) {
#ifdef TIMINGS_INDEFINITE_LLL
        os << "ILLL: |ComputeReductionIndefinite(None)|=" << time << "\n";
#endif
        return {B, Mwork, {}};
      }
    } else {
      iter_no_improv = 0;
      norm_work = norm;
    }
    // Applying the random perturbation and iterating
    MyMatrix<Tint> RandUnit = get_random_int_matrix<Tint>(n);
    MyMatrix<T> RandUnit_T = UniversalMatrixConversion<T, Tint>(RandUnit);
    B = RandUnit * B;
    B_T = RandUnit_T * B_T;
    Mwork = RandUnit_T * Mwork * RandUnit_T.transpose();
  }
}

template <typename T, typename Tint> struct ResultReduction {
  MyMatrix<Tint> B;
  MyMatrix<T> Mred;
};

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

template <typename T, typename Tint>
ResultReduction<T, Tint>
ComputeReductionIndefinite_opt(MyMatrix<T> const &M, bool const &ApplyReduction,
                               std::ostream &os) {
  if (ApplyReduction) {
    ResultIndefiniteLLL<T, Tint> res =
        ComputeReductionIndefinite<T, Tint>(M, os);
    return {res.B, res.Mred};
  } else {
    int n = M.rows();
    MyMatrix<Tint> B = IdentityMat<Tint>(n);
    MyMatrix<T> Mred = M;
    return {B, Mred};
  }
}

// clang-format off
#endif  //  SRC_INDEFINITE_INDEFINITE_LLL_H_
// clang-format on
