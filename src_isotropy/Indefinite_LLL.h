// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_INDEFINITE_LLL_H_
#define SRC_INDEFINITE_INDEFINITE_LLL_H_

// clang-format off
#include "MAT_Matrix.h"
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalBasic.h"
#include "SignatureSymmetric.h"
#include "ClassicLLL.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_INDEFINITE_LLL
#define DEBUG_SIMPLE_INDEFINITE_REDUCTION
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
  T det2(1);
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
  if (n == 0) {
    MyMatrix<Tint> B = IdentityMat<Tint>(0);
    MyMatrix<T> Mred = M;
    MyVector<T> Xisotrop = MyVector<T>(0);
    std::optional<MyVector<T>> opt{Xisotrop};
    return {B, Mred, opt};
  }
  if (n == 1) {
    MyMatrix<Tint> B = IdentityMat<Tint>(1);
    MyMatrix<T> Mred = M;
    auto get_iso=[&]() -> std::optional<MyVector<T>> {
      if (M(0,0) == 0) {
        MyVector<T> Xisotrop = MyVector<T>(1);
        Xisotrop(0) = 1;
        return Xisotrop;
      } else {
        return {};
      }
    };
    return {B, Mred, get_iso()};
  }
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
    std::cerr << "ILLL: M=\n";
    WriteMatrix(std::cerr, M);
    std::cerr << "ILLL: B=\n";
    WriteMatrix(std::cerr, B);
#endif
    ResultGramSchmidt_Indefinite<T> ResGS = GramSchmidtOrthonormalization(M, B);
#ifdef DEBUG_INDEFINITE_LLL
    std::cerr << "ILLL: Bstar_norms =";
    for (auto &eN : ResGS.Bstar_norms)
      std::cerr << " " << eN;
    std::cerr << "\n";
#endif
    if (!ResGS.success) {
#ifdef DEBUG_INDEFINITE_LLL
      std::cerr << "ILLL: Returning from Indefinite_LLL (!ResGS.success)\n";
#endif
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
#ifdef DEBUG_INDEFINITE_LLL
    std::cerr << "ILLL: sum1=" << sum1 << " sum2=" << sum2 << "\n";
#endif
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
    if (k >= n) {
      break;
    }
  }
#ifdef DEBUG_INDEFINITE_LLL
  std::cerr << "ILLL: Returning from Indefinite_LLL 2\n";
#endif
  return {B, get_matrix(), {}};
}

/*
  Permute the columns randomly and change the signs randomly.
 */
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

/*
  Compute a reduced form of the quadratic form by a
  pack of heuristics.
  --
  The goal can be to ultimately find an isotropic vector
  in which case the option "look_for_isotropic" has to
  be selected.
 */
template <typename T, typename Tint>
ResultIndefiniteLLL<T, Tint>
ComputeReductionIndefinite(MyMatrix<T> const &M,
                           bool const& look_for_isotropic,
                           [[maybe_unused]] std::ostream &os) {
  int n = M.rows();
  if (n == 0) {
    MyMatrix<Tint> B = IdentityMat<Tint>(0);
    MyMatrix<T> Mred = M;
    MyVector<T> Xisotrop = MyVector<T>(0);
    std::optional<MyVector<T>> opt = {Xisotrop};
    return {std::move(B), std::move(Mred), opt};
  }
#ifdef TIMINGS_INDEFINITE_LLL
  MicrosecondTime time;
#endif
#ifdef DEBUG_INDEFINITE_LLL
  os << "ILLL: Beginning of ComputeReductionIndefinite\n";
#endif
  auto get_norm = [&](MyMatrix<T> const &mat) -> T {
    T sum(0);
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
    if (res.Xisotrop && look_for_isotropic) {
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
      os << "|ILLL: ComputeReductionIndefinite(Iso)|=" << time << "\n";
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
        os << "|ILLL: ComputeReductionIndefinite(None)|=" << time << "\n";
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
ResultReduction<T, Tint>
LLLIndefiniteReduction(MyMatrix<T> const &M, std::ostream &os) {
  bool look_for_isotropic = false;
  ResultIndefiniteLLL<T, Tint> res =
    ComputeReductionIndefinite<T, Tint>(M, look_for_isotropic, os);
  return {std::move(res.B), std::move(res.Mred)};
}

template<typename T>
GraphBitset GetMatrixGraph(MyMatrix<T> const& M) {
  size_t n = M.rows();
  GraphBitset eG(n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if (M(i, j) != 0) {
        eG.AddAdjacent(i, j);
        eG.AddAdjacent(j, i);
      }
    }
  }
  return eG;
}



template<typename T>
std::vector<std::vector<size_t>> MatrixConnectedComponents(MyMatrix<T> const& M) {
  GraphBitset eG = GetMatrixGraph(M);
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(eG);
  return LConn;
}

template<typename T>
bool IsMatrixConnected(MyMatrix<T> const& M) {
  GraphBitset eG = GetMatrixGraph(M);
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(eG);
  return LConn.size() == 1;
}



template<typename T, typename Tint>
ResultReduction<T, Tint>
SimpleIndefiniteReduction(MyMatrix<T> const &M, [[maybe_unused]] std::ostream &os) {
  int n = M.rows();
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
  auto l1_norm=[&](MyMatrix<T> const& U) -> T {
    T sum = 0;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        sum += T_abs(U(i,j));
      }
    }
    return sum;
  };
#endif
  MyMatrix<Tint> B = IdentityMat<Tint>(n);
  MyMatrix<T> Mwork = M;
  // We choose the indices at random in order to avoid always
  // hitting the same indices
  std::vector<int> indices;
  for (int i=0; i<n; i++)
    indices.push_back(i);
  auto f_transpose=[&]() -> void {
    int i = random() % n;
    int j = random() % n;
    if (i != j) {
      int k = indices[i];
      indices[i] = indices[j];
      indices[j] = k;
    }
  };
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
  os << "=======================================================\n";
  os << "ILLL: n=" << n << "\n";
  os << "ILLL: M=\n";
  WriteMatrix(os, M);
  os << "ILLL: L1(M)=" << l1_norm(M) << "\n";
#endif
  struct ResSearch {
    int i;
    int j;
    Tint c;
  };
  // Apply the transformation
  // tilde(M) = (Id + c U_ij) M (Id + c U_ji)
  //          = (M + c Row(M,j) at row i) ( Id + c U_ji)
  //          = M + c Row(M,j) at row i + c Col(M,j) at col i + c^2
  auto eval=[&](ResSearch const& x) -> T {
    int i = x.i;
    int j = x.j;
    Tint c = x.c;
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
    os << "ILLL: eval: i=" << i << " j=" << j << " c=" << c << "\n";
#endif
    //
    T delta_off = 0;
    for (int k=0; k<n; k++) {
      if (k != i) {
        T val1 = T_abs( Mwork(i,k) );
        T val2 = T_abs( T(Mwork(i,k) + c * Mwork(j,k)) );
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
        os << "ILLL: eval, off: k=" << k << " val1=" << val1 << " val2=" << val2 << "\n";
#endif
        delta_off += val1 - val2;
      }
    }
    T val1 = T_abs( Mwork(i,i) );
    T val2 = T_abs( T(Mwork(i,i) + 2 * c * Mwork(j,i) + c * c * Mwork(j,j)) );
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
    os << "ILLL: eval, diag, val1=" << val1 << " val2=" << val2 << "\n";
#endif
    T delta_diag = val1 - val2;
    T delta = 2 * delta_off + delta_diag;
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
    os << "ILLL: eval, summary, delta_off=" << delta_off << " delta_diag=" << delta_diag << " delta=" << delta << "\n";
#endif
    return delta;
  };
  auto direction_search=[&](ResSearch const& x) -> std::optional<ResSearch> {
    T val = eval(x);
    if (val <= 0) {
      return {};
    }
    T best_delta = val;
    ResSearch x_ret = x;
    while(true) {
      x_ret.c += x.c;
      T val = eval(x_ret);
      if (val <= best_delta) {
        x_ret.c -= x.c; // Revert the change that would increase the norm.
        return x_ret;
      }
      best_delta = val;
    }
  };
  auto f_search=[&]() -> std::optional<ResSearch> {
    f_transpose();
    for (int i_s=0; i_s<n; i_s++) {
      int i = indices[i_s];
      for (int j_s=0; j_s<n; j_s++) {
        int j = indices[j_s];
        if (i != j) {
          for (int pre_sign=0; pre_sign<2; pre_sign++) {
            Tint sign = Tint(-1 + 2 * pre_sign);
            ResSearch x{i, j, sign};
            std::optional<ResSearch> opt = direction_search(x);
            if (opt) {
              ResSearch const& x_final = *opt;
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
              os << "ILLL: x_final=(" << x_final.i << " / " << x_final.j << " / " << x_final.c << ")\n";
#endif
              return x_final;
            }
          }
        }
      }
    }
    return {};
  };
  auto update=[&](ResSearch const& x) -> void {
    int i = x.i;
    int j = x.j;
    Tint c = x.c;
    T c_T = UniversalScalarConversion<T,Tint>(c);
    B.row(i) += c * B.row(j);
    Mwork.row(i) += c_T * Mwork.row(j);
    Mwork.col(i) += c_T * Mwork.col(j);
  };
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
  auto f_compute=[&](ResSearch const& x) -> std::pair<MyMatrix<Tint>, MyMatrix<T>> {
    int i = x.i;
    int j = x.j;
    T c = x.c;
    MyMatrix<Tint> eUnit = IdentityMat<Tint>(n);
    eUnit(i,j) = c;
    MyMatrix<Tint> Btest = eUnit * B;
    MyMatrix<T> BT_test = UniversalMatrixConversion<T,Tint>(Btest);
    MyMatrix<T> Mwork_test = BT_test * M * BT_test.transpose();
    return {std::move(Btest), std::move(Mwork_test)};
  };
  size_t n_oper = 0;
#endif
  //
  while(true) {
    std::optional<ResSearch> opt = f_search();
    if (opt) {
      ResSearch const& x = *opt;
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
      n_oper += 1;
      os << "------------ n_oper=" << n_oper << " --------------------\n";
      os << "ILLL: i=" << x.i << " j=" << x.j << " c=" << x.c << "\n";
      T norm_prev = l1_norm(Mwork);
      os << "ILLL: Mwork=\n";
      WriteMatrix(os, Mwork);
      std::pair<MyMatrix<Tint>, MyMatrix<T>> pair = f_compute(x);
      os << "ILLL: pair.first(B)=\n";
      WriteMatrix(os, pair.first);
      os << "ILLL: pair.second(M)=\n";
      WriteMatrix(os, pair.second);
      os << "ILLL: L1(Mwork)=" << l1_norm(Mwork) << " L1(pair.second)=" << l1_norm(pair.second) << "\n";
      T delta = eval(x);
#endif
      update(x);
#ifdef DEBUG_SIMPLE_INDEFINITE_REDUCTION
      T norm_next = l1_norm(Mwork);
      size_t num_error = 0;
      if (norm_prev - delta != norm_next) {
        std::cerr << "ILLL: norm_prev=" << norm_prev << " delta=" << delta << " norm_next=" << norm_next << "\n";
        std::cerr << "ILLL: B=\n";
        WriteMatrix(std::cerr, B);
        std::cerr << "ILLL: Mwork=\n";
        WriteMatrix(std::cerr, Mwork);
        std::cerr << "L1-Norm inconsistency\n";
        num_error += 1;
      }
      if (pair.first != B) {
        std::cerr << "B inconsistency\n";
        num_error += 1;
      }
      if (pair.second != Mwork) {
        std::cerr << "Mwork inconsistency\n";
        num_error += 1;
      }
      MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
      MyMatrix<T> prod = B_T * M * B_T.transpose();
      if (prod != Mwork) {
        std::cerr << "B / Mwork inconsistency\n";
        num_error += 1;
      }
      if (num_error > 0) {
        std::cerr << "ILLL: num_error=" << num_error << "\n";
        throw TerminalException{1};
      }
#endif
    } else {
      return {std::move(B), std::move(Mwork)};
    }
  }
}

template<typename T>
T get_l1_norm(MyMatrix<T> const& U) {
  int n = U.rows();
  T sum(0);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      sum += T_abs(U(i,j));
    }
  }
  return sum;
}



template<typename T, typename Tint>
ResultReduction<T, Tint>
BlockSimpleIndefiniteReduction(MyMatrix<T> const &M, std::ostream &os) {
  int n = M.rows();
  MyMatrix<Tint> B = ZeroMatrix<Tint>(n, n);
  MyMatrix<T> Mret = ZeroMatrix<T>(n, n);
  std::vector<std::vector<size_t>> LConn = MatrixConnectedComponents(M);
  for (auto & eConn : LConn) {
    int len = eConn.size();
    MyMatrix<T> M_conn(len, len);
    for (int i=0; i<len; i++) {
      for (int j=0; j<len; j++) {
        size_t i_big = eConn[i];
        size_t j_big = eConn[j];
        M_conn(i, j) = M(i_big, j_big);
      }
    }
    ResultReduction<T, Tint> res = SimpleIndefiniteReduction<T,Tint>(M_conn, os);
    for (int i=0; i<len; i++) {
      for (int j=0; j<len; j++) {
        size_t i_big = eConn[i];
        size_t j_big = eConn[j];
        B(i_big, j_big) = res.B(i,j);
        Mret(i_big, j_big) = res.Mred(i,j);
      }
    }
  }
  return {std::move(B), std::move(Mret)};
}

template <typename T>
std::pair<int, int> GetSignature(MyMatrix<T> const &M) {
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(M);
  int nbPlus = DiagInfo.nbPlus;
  int nbMinus = DiagInfo.nbMinus;
  return {nbPlus, nbMinus};
}

template <typename T, typename Tint>
ResultReduction<T, Tint>
IndefiniteReduction(MyMatrix<T> const &M, std::ostream &os) {
  int n = M.rows();
  MyMatrix<Tint> B = ZeroMatrix<Tint>(n, n);
  MyMatrix<T> Mwork = M;
  std::pair<int, int> signature = GetSignature(M);
  int n_plus = signature.first;
  int n_minus = signature.second;
  int n_choice = 4;
  if (n_plus > 0 && n_minus > 0) {
    n_choice = 2;
  }
  auto compute_by_block=[&](auto & Min) -> std::optional<ResultReduction<T,Tint>> {
    std::vector<std::vector<size_t>> LConn = MatrixConnectedComponents(Min);
    if (LConn.size() == 1) {
      return {};
    }
    MyMatrix<Tint> B_block = ZeroMatrix<Tint>(n, n);
    MyMatrix<T> M_block = ZeroMatrix<T>(n, n);
    for (auto & eConn : LConn) {
      int len = eConn.size();
      MyMatrix<T> M_conn(len, len);
      for (int i=0; i<len; i++) {
        for (int j=0; j<len; j++) {
          size_t i_big = eConn[i];
          size_t j_big = eConn[j];
          M_conn(i, j) = Min(i_big, j_big);
        }
      }
      ResultReduction<T, Tint> res = IndefiniteReduction<T,Tint>(M_conn, os);
      for (int i=0; i<len; i++) {
        for (int j=0; j<len; j++) {
          size_t i_big = eConn[i];
          size_t j_big = eConn[j];
          B_block(i_big, j_big) = res.B(i,j);
          M_block(i_big, j_big) = res.Mred(i,j);
        }
      }
    }
    MyMatrix<Tint> Bret = B_block * B;
    ResultReduction<T,Tint> res{std::move(Bret), std::move(M_block)};
    return res;
  };
  auto compute_reduction=[&](MyMatrix<T> const& Min, int const& choice) -> ResultReduction<T, Tint> {
    if (choice == 0) {
      return LLLIndefiniteReduction<T,Tint>(Min, os);
    }
    if (choice == 1) {
      return SimpleIndefiniteReduction<T,Tint>(Min, os);
    }
    if (choice == 2 || choice == 3) {
      if (n_plus > 0 && n_minus > 0) {
        MyMatrix<Tint> B = IdentityMat<Tint>(n);
        return {std::move(B), Min};
      }
      auto get_m=[&]() -> MyMatrix<T> {
        if (n_plus == 0) {
          return -Min;
        }
        return Min;
      };
      MyMatrix<T> Mcall = get_m();
      auto get_redmat=[&]() -> MyMatrix<Tint> {
        if (choice == 2) {
          LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(Mcall);
          return rec.Pmat;
        } else {
          LLLreduction<T, Tint> rec = LLLreducedBasisDual<T, Tint>(Mcall);
          return rec.Pmat;
        }
      };
      MyMatrix<Tint> B = get_redmat();
      MyMatrix<T> B_T = UniversalMatrixConversion<T,Tint>(B);
      MyMatrix<T> Mout = B_T * Min * B_T.transpose();
      return {std::move(B), std::move(Mout)};
    }
    std::cerr << "The chosen option is not accepted\n";
    throw TerminalException{1};
  };


  T norm_work = get_l1_norm(M);
#ifdef DEBUG_INDEFINITE_LLL
  size_t n_iter = 0;
#endif
  while(true) {
    size_t n_operation = 0;
    for (int choice=0; choice<n_choice; choice++) {
      std::optional<ResultReduction<T,Tint>> opt = compute_by_block(Mwork);
      if (opt) {
        return *opt;
      }
      ResultReduction<T,Tint> res = compute_reduction(Mwork, choice);
      T norm = get_l1_norm(res.Mred);
      if (norm < norm_work) {
        n_operation += 1;
        B = res.B * B;
        Mwork = res.Mred;
        norm_work = norm;
      }
    }
#ifdef DEBUG_INDEFINITE_LLL
    os << "ILLL: n_iter=" << n_iter << " n_operation=" << n_operation << "\n";
    n_iter += 1;
#endif
    // If did nothing then exit. Necessarily we will reach that stage at some point
    if (n_operation < 2) {
      return {std::move(B), std::move(Mwork)};
    }
  }
}

template <typename T, typename Tint>
ResultReduction<T, Tint>
IndefiniteReduction_opt(MyMatrix<T> const &M, bool const &ApplyReduction,
                        std::ostream &os) {
  if (ApplyReduction) {
    return IndefiniteReduction<T,Tint>(M, os);
  } else {
    int n = M.rows();
    MyMatrix<Tint> B = IdentityMat<Tint>(n);
    MyMatrix<T> Mred = M;
    return {std::move(B), std::move(Mred)};
  }
}

// clang-format off
#endif  //  SRC_INDEFINITE_INDEFINITE_LLL_H_
// clang-format on
