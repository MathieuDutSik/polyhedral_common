// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_INDEFINITE_INDEFINITE_LLL_H_
#define SRC_INDEFINITE_INDEFINITE_LLL_H_

// clang-format off
#include "MAT_Matrix.h"
#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
// clang-format on

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
    return {B, Mred, opt};
  }
#ifdef TIMINGS_INDEFINITE_LLL
  MicrosecondTime time;
#endif
#ifdef DEBUG_INDEFINITE_LLL
  os << "ILLL: Beginning of ComputeReductionIndefinite\n";
#endif
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
IndefiniteReduction(MyMatrix<T> const &M, std::ostream &os) {
  bool look_for_isotropic = false;
  ResultIndefiniteLLL<T, Tint> res =
    ComputeReductionIndefinite<T, Tint>(M, look_for_isotropic, os);
  return {std::move(res.B), std::move(res.Mred)};
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
    return {B, Mred};
  }
}

// clang-format off
#endif  //  SRC_INDEFINITE_INDEFINITE_LLL_H_
// clang-format on
