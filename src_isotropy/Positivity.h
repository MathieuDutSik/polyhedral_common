// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_POSITIVITY_H_
#define SRC_LATT_POSITIVITY_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "COMB_Vectors.h"
#include "GraverBasis.h"
#include "COMB_Combinatorics.h"
#include "SignatureSymmetric.h"
#include "Indefinite_LLL.h"
#include <vector>
#include <unordered_map>
#include <utility>
// clang-format on

// This code does the following:
// * diagonalization of symmetric matrices
// * Finding a vector of positive norm for a symmetric matrix.
// * Find a short positive norm vector for a symmetric matrix.

#ifdef DEBUG
#define DEBUG_POSITIVITY
#endif

template <typename T> T MinimumDiagonal(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  T MinNorm = eMat(0, 0);
  for (int i = 1; i < n; i++) {
    T eVal = eMat(i, i);
    if (eVal < MinNorm)
      MinNorm = eVal;
  }
  return MinNorm;
}

template <typename T> T MaximumDiagonal(MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  T MaxNorm = eMat(0, 0);
  for (int i = 1; i < n; i++) {
    T eVal = eMat(i, i);
    if (eVal > MaxNorm)
      MaxNorm = eVal;
  }
  return MaxNorm;
}

template <typename T> MyMatrix<T> AnLattice(int const &n) {
  MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      eMat(i, j) = 0;
  for (int i = 0; i < n; i++)
    eMat(i, i) = 2;
  for (int i = 1; i < n; i++) {
    eMat(i, i - 1) = -1;
    eMat(i - 1, i) = -1;
  }
  return eMat;
}

template <typename T>
MyVector<T> GetPositiveNormVector(MyMatrix<T> const &SymMat, std::ostream& os) {
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(SymMat, os);
  MyMatrix<T> const &Transform = DiagInfo.Transform;
  MyMatrix<T> const &RedMat = DiagInfo.RedMat;
  int n = SymMat.rows();
  for (int i = 0; i < n; i++) {
    if (RedMat(i, i) > 0) {
      MyVector<T> eVect = GetMatrixRow(Transform, i);
#ifdef DEBUG_POSITIVITY
      T norm = eVect.dot(SymMat * eVect);
      if (norm != RedMat(i, i)) {
        std::cerr << "We have a consistency error\n";
        throw TerminalException{1};
      }
#endif
      return eVect;
    }
  }
  std::cerr << "Failed to find a negative norm vector\n";
  throw TerminalException{1};
}

template <typename T, typename Tint, typename Ttest>
std::optional<MyVector<Tint>>
TestIntegralVector_family(std::vector<MyVector<Ttest>> const &ListVect,
                          MyMatrix<T> const &M, T const &CritNorm,
                          bool const &StrictIneq,
                          [[maybe_unused]] std::ostream &os) {
  int n = M.rows();
  MyVector<Tint> V_ret(n);
  for (auto &eV : ListVect) {
    auto is_ok_vector = [&]() -> bool {
      T eNorm = EvaluationQuadForm<T, Tint>(M, V_ret);
      if (StrictIneq) {
        return eNorm > CritNorm;
      } else {
        return eNorm >= CritNorm;
      }
    };
    bool test = is_ok_vector();
    if (test) {
      return V_ret;
    }
  }
  return {};
}

template <typename T, typename Tint, typename Ttest>
MyVector<Tint>
GetIntegralVector_family(std::vector<MyVector<Ttest>> const &ListVect,
                         MyMatrix<T> const &M, T const &CritNorm,
                         bool const &StrictIneq, std::ostream &os) {
  int n = M.rows();
  int scal = 1;
  MyVector<Tint> V_ret(n);
  while (true) {
    std::optional<MyVector<Tint>> opt =
        TestIntegralVector_family<T, Tint, Ttest>(ListVect, scal, M, CritNorm,
                                                  StrictIneq, os);
    if (opt) {
      return *opt;
    }
#ifdef DEBUG_POSITIVITY
    os << "POS: GetPositiveNormVector_family |ListVect|=" << ListVect.size()
       << " scal=" << scal << "\n";
#endif
    scal += 1;
  }
}

template <typename T>
std::vector<MyVector<T>> GetPositiveDirections_diag(MyMatrix<T> const &M, std::ostream& os) {
  int n = M.rows();
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(M, os);
  MyMatrix<T> const &Transform = DiagInfo.Transform;
  MyMatrix<T> const &RedMat = DiagInfo.RedMat;
  std::vector<MyVector<T>> ListVect;
  for (int i = 0; i < n; i++) {
    if (RedMat(i, i) > 0) {
      MyVector<T> eVect = GetMatrixRow(Transform, i);
      T sum = 0;
      for (int j = 0; j < n; j++) {
        sum += T_abs(eVect(j));
      }
      MyVector<T> fVect = eVect / sum;
      ListVect.push_back(fVect);
    }
  }
  return ListVect;
}

template <typename T, typename Tint>
MyVector<Tint> GetIntegralVector_diag(MyMatrix<T> const &M, T const &CritNorm,
                                      bool const &StrictIneq,
                                      std::ostream &os) {
  std::vector<MyVector<T>> ListVect = GetPositiveDirections_diag(M, os);
  return GetIntegralVector_family<T, Tint, T>(ListVect, M, CritNorm, StrictIneq,
                                              os);
}

template <typename T>
std::vector<MyVector<double>>
GetPositiveDirections_eigen(MyMatrix<T> const &M) {
  int n = M.rows();
  MyMatrix<double> M_double = UniversalMatrixConversion<double, T>(M);
  Eigen::SelfAdjointEigenSolver<MyMatrix<double>> eig(M_double);
  MyVector<double> ListEig = eig.eigenvalues();
  MyMatrix<double> ListVect = eig.eigenvectors();
  std::vector<MyVector<double>> ListEigVect;
  for (int i = 0; i < n; i++) {
    if (ListEig(i) > 0) {
      MyVector<double> V(n);
      for (int j = 0; j < n; j++) {
        V(j) = ListVect(i, j);
      }
      ListEigVect.push_back(V);
    }
  }
  return ListEigVect;
}

template <typename T, typename Tint>
MyVector<Tint> GetIntegralVector_eigen(MyMatrix<T> const &M, T const &CritNorm,
                                       bool const &StrictIneq,
                                       std::ostream &os) {
  std::vector<MyVector<double>> ListEigVect = GetPositiveDirections_eigen(M);
  return GetIntegralVector_family<T, Tint, double>(ListEigVect, M, CritNorm,
                                                   StrictIneq, os);
}

// That code applies several techniques together in order to get a short
// integral vector V satisfying A[V] >= CritNorm
// * if StrictIneq is true then we want actually A[V] > CritNorm
//
// * The diagonalization method should work all the time but maybe not get us
// an optimal vector
// * The eigenvector method is not guaranteed to work. But if it works, it
// should get us a very good vector.
template <typename T, typename Tint>
MyVector<Tint>
GetIntegralVector_allmeth_V1(MyMatrix<T> const &M, T const &CritNorm,
                             bool const &StrictIneq, std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth_V1: Beginning\n";
#endif
  std::vector<MyVector<T>> ListVect = GetPositiveDirections_diag(M, os);
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth_V1: |ListVect|=" << ListVect.size()
     << "\n";
#endif
  std::vector<MyVector<double>> ListEigVect = GetPositiveDirections_eigen(M);
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth_V1: |ListEigVect|="
     << ListEigVect.size() << "\n";
#endif
  int scal = 1;
  while (true) {
#ifdef DEBUG_POSITIVITY
    os << "POS: GetIntegralVector_allmeth_V1: scal=" << scal << "\n";
#endif
    std::optional<MyVector<Tint>> opt1 = TestIntegralVector_family<T, Tint, T>(
        ListVect, scal, M, CritNorm, StrictIneq, os);
    if (opt1) {
      return *opt1;
    }
    std::optional<MyVector<Tint>> opt2 =
        TestIntegralVector_family<T, Tint, double>(ListEigVect, scal, M,
                                                   CritNorm, StrictIneq, os);
    if (opt2) {
      return *opt2;
    }
    scal += 1;
  }
}

// The ApproxIterator allows to find approximation of a vector by
// integral vectors.
// The sequence of approximant is
// V_app(N) = List(Nint( i V(x))
// The call to "increment()" will just go from one approximant to the next.
template <typename T, typename Tint> struct ApproxIterator {
  MyVector<Tint> Vapprox;
  MyVector<T> l_step;
  MyVector<T> l_next;
  std::vector<std::vector<std::pair<int, int>>> ll_coef;
  ApproxIterator(MyVector<T> const &V) {
    int n = V.size();
    Vapprox = ZeroVector<Tint>(n);
    std::unordered_map<T, std::vector<std::pair<int, int>>> map;
    for (int i = 0; i < n; i++) {
      T val = V(i);
      if (val != 0) {
        T val_abs = T_abs(val);
        T inv_val = 1 / val_abs;
        int sign = T_sign(val);
        std::pair<int, int> pair{i, sign};
        map[inv_val].push_back(pair);
      }
    }
    int len = map.size();
    l_step = ZeroVector<T>(len);
    l_next = ZeroVector<T>(len);
    int pos = 0;
    for (auto &kv : map) {
      l_step(pos) = kv.first;
      l_next(pos) = kv.first / 2;
      ll_coef.push_back(kv.second);
      pos += 1;
    }
#ifdef DEBUG_POSITIVITY
    std::cerr << "POS: ApproxIterator len=" << len << "\n";
#endif
  }
  MyVector<Tint> increment() {
    int len = l_step.size();
    std::vector<int> l_idx_inc{0};
    T min_val = l_next(0);
    for (int u = 1; u < len; u++) {
      if (l_next(u) < min_val) {
        min_val = l_next(u);
        l_idx_inc.clear();
        l_idx_inc.push_back(u);
      } else {
        if (l_next(u) == min_val) {
          l_idx_inc.push_back(u);
        }
      }
    }
    for (auto &u : l_idx_inc) {
      l_next(u) += l_step(u);
      for (auto &pair : ll_coef[u]) {
        Vapprox(pair.first) += pair.second;
      }
    }
    return Vapprox;
  }
};

// Use an increment to approximate the vectors.
// It is more efficient than the scaling method.
//
// The possible positive directions are obtained from
// two sources:
// * The exact diagonalization of the matrix. It should
//   terminate all the time.
// * The diagonalization using eigenvectors of the matrix,.
template <typename T, typename Tint>
MyVector<Tint> GetIntegralVector_allmeth_V2(MyMatrix<T> const &M,
                                            T const &CritNorm,
                                            bool const &StrictIneq,
                                            std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth_V2: Beginning\n";
#endif
  std::vector<ApproxIterator<T, Tint>> l_approx_diag;
  for (auto &eVec : GetPositiveDirections_diag(M, os)) {
#ifdef DEBUG_POSITIVITY
    os << "POS: GetIntegralVector_allmeth_V2: diag, eVec="
       << StringVectorGAP(eVec) << "\n";
#endif
    ApproxIterator<T, Tint> ai(eVec);
    l_approx_diag.push_back(ai);
  }
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth_V2: |l_approx_diag|="
     << l_approx_diag.size() << "\n";
#endif
  std::vector<ApproxIterator<double, Tint>> l_approx_eigen;
  for (auto &eVec : GetPositiveDirections_eigen(M)) {
#ifdef DEBUG_POSITIVITY
    os << "POS: GetIntegralVector_allmeth_V2: eigen, eVec="
       << StringVectorGAP(eVec) << "\n";
#endif
    ApproxIterator<double, Tint> ai(eVec);
    l_approx_eigen.push_back(ai);
  }
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth_V2: |l_approx_eigen|="
     << l_approx_eigen.size() << "\n";
#endif
  auto is_vector_correct = [&](MyVector<Tint> const &V) -> bool {
    T eNorm = EvaluationQuadForm<T, Tint>(M, V);
    if (StrictIneq) {
      return eNorm > CritNorm;
    } else {
      return eNorm >= CritNorm;
    }
  };
#ifdef DEBUG_POSITIVITY
  size_t iter = 0;
#endif
  while (true) {
#ifdef DEBUG_POSITIVITY
    os << "POS: GetIntegralVector_allmeth_V2: beginning of loop, iter=" << iter
       << "\n";
#endif
    // First trying the vectors from diagonalization
    for (auto &ai : l_approx_diag) {
      MyVector<Tint> V = ai.increment();
#ifdef DEBUG_POSITIVITY
      os << "POS: GetIntegralVector_allmeth_V2: diag V=" << StringVectorGAP(V)
         << "\n";
#endif
      if (is_vector_correct(V)) {
#ifdef DEBUG_POSITIVITY
        os << "POS: GetIntegralVector_allmeth_V2: returning from "
              "diagonolization vector\n";
#endif
        return V;
      }
    }
    // Second trying the vectors from eigen
    for (auto &ai : l_approx_eigen) {
      MyVector<Tint> V = ai.increment();
#ifdef DEBUG_POSITIVITY
      os << "POS: GetIntegralVector_allmeth_V2: eigen V=" << StringVectorGAP(V)
         << "\n";
#endif
      if (is_vector_correct(V)) {
#ifdef DEBUG_POSITIVITY
        os << "POS: GetIntegralVector_allmeth_V2: returning from eigen "
              "vector\n";
#endif
        return V;
      }
    }
#ifdef DEBUG_POSITIVITY
    iter += 1;
#endif
  }
}

// This is the algorithm for computing a vector V with M[V] >= CritNorm.
// If StrictIneq = true, then the inequality becomes strict.
//
// The algorithm is the following:
// * If finding an isotropic vector would solve the problem then that
//   path is taken. That is if we find an isotropic vector, then return
// * The algorithm is run on the Iterated Indefinite-LLL reduced form.
// * That algorithm is called in allmeth_V2.
template <typename T, typename Tint>
MyVector<Tint>
GetIntegralVector_allmeth(MyMatrix<T> const &M, T const &CritNorm,
                          bool const &StrictIneq, std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth: beginning\n";
#endif
  auto is_isotropic_fine = [&]() -> bool {
    if (StrictIneq) {
      return CritNorm < 0;
    } else {
      return CritNorm <= 0;
    }
  };
  bool test_isotropic_fine = is_isotropic_fine();
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth: CritNorm=" << CritNorm
     << " StrictIneq=" << StrictIneq << "\n";
  os << "POS: GetIntegralVector_allmeth: test_isotropic_fine="
     << test_isotropic_fine << "\n";
  os << "POS: M=\n";
  WriteMatrix(os, M);
#endif
  //
  // No nicer way to proceed right now. We have two scenarios to handle.
  //
  std::optional<MyVector<T>> res_Xisotrop;
  MyMatrix<Tint> res_B;
  MyMatrix<T> res_Mred;
  if (test_isotropic_fine) {
    bool look_for_isotropic = true;
    ResultIndefiniteLLL<T, Tint> res = ComputeReductionIndefinite<T, Tint>(M, look_for_isotropic, os);
    res_Xisotrop = res.Xisotrop;
    res_B = res.B;
    res_Mred = res.Mred;
  } else {
    ResultReduction<T, Tint> res = IndefiniteReduction<T,Tint>(M, os);
    res_B = res.B;
    res_Mred = res.Mred;
  }
  //
  // Now doing the work.
  //
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralVector_allmeth: We have res\n";
#endif
  if (test_isotropic_fine && res_Xisotrop) {
    MyVector<T> const &Xisotrop = *res_Xisotrop;
#ifdef DEBUG_POSITIVITY
    os << "POS: GetIntegralVector_allmeth: Xisotrop="
       << StringVectorGAP(Xisotrop) << "\n";
#endif
    MyVector<T> V1 = RemoveFractionVector(Xisotrop);
#ifdef DEBUG_POSITIVITY
    os << "POS: GetIntegralVector_allmeth: V1=" << StringVectorGAP(V1) << "\n";
#endif
    MyVector<Tint> V2 = UniversalVectorConversion<Tint, T>(V1);
#ifdef DEBUG_POSITIVITY
    os << "POS: GetIntegralVector_allmeth: V2=" << StringVectorGAP(V2) << "\n";
#endif
    return V2;
  }
#ifdef DEBUG_POSITIVITY
  os << "POS: The gambit using isotropic vectors failed, going classically\n";
  os << "POS: GetIntegralVector_allmeth: Before GetIntegralVector_allmeth_V2 "
        "for Mred\n";
#endif
  MyVector<Tint> V1 =
      GetIntegralVector_allmeth_V2<T, Tint>(res_Mred, CritNorm, StrictIneq, os);
#ifdef DEBUG_POSITIVITY
  os << "POS: res_Mred=\n";
  WriteMatrix(os, res_Mred);
  os << "POS: GetIntegralVector_allmeth: We have V1\n";
#endif
  MyVector<Tint> V2 = res_B.transpose() * V1;
#ifdef DEBUG_POSITIVITY
  T norm_red = EvaluationQuadForm(res_Mred, V1);
  T norm = EvaluationQuadForm(M, V2);
  if (norm_red != norm) {
    std::cerr << "norm=" << norm << " norm_red=" << norm_red << "\n";
    std::cerr << "But they should be equal\n";
    throw TerminalException{1};
  }
  os << "POS: GetIntegralVector_allmeth: We have V2=" << StringVectorGAP(V2) << "\n";
#endif
  return V2;
}

template <typename T, typename Tint>
MyVector<Tint> GetIntegralPositiveVector_allmeth(MyMatrix<T> const &M,
                                                 std::ostream &os) {
  T CritNorm(0);
  bool StrictIneq = true;
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralPositiveVector_allmeth: before "
        "GetIntegralVector_allmeth\n";
#endif
  MyVector<Tint> V = GetIntegralVector_allmeth<T, Tint>(M, CritNorm, StrictIneq, os);
#ifdef DEBUG_POSITIVITY
  os << "POS: GetIntegralPositiveVector_allmeth: We have V\n";
#endif
  return V;
}

template <typename T, typename Tint>
MyVector<Tint> GetShortIntegralVector(MyMatrix<T> const &M, T const &CritNorm,
                                      bool const &StrictIneq,
                                      std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: GetShortIntegralVector: beginning\n";
#endif
  MyMatrix<T> Mwork = -M;
  T CritNormWork = -CritNorm;
#ifdef DEBUG_POSITIVITY
  os << "POS: GetShortIntegralVector: before GetIntegralVector_allmeth\n";
#endif
  MyVector<Tint> V = GetIntegralVector_allmeth<T, Tint>(Mwork, CritNormWork, StrictIneq, os);
#ifdef DEBUG_POSITIVITY
  os << "POS: GetShortIntegralVector: We have V\n";
#endif
  return V;
}

template <typename Tint>
MyMatrix<Tint> GetRandomMatrixPerturbation(int const &n) {
  int choice = rand() % 2;
  auto get_rnd_sign = []() -> int {
    int pos = rand() % 2;
    return 2 * pos - 1;
  };
  if (choice == 0) {
    std::vector<int> ePerm = RandomPermutation<int>(n);
    MyMatrix<Tint> eMat = ZeroMatrix<Tint>(n, n);
    for (int iRow = 0; iRow < n; iRow++) {
      int iCol = ePerm[iRow];
      eMat(iRow, iCol) = get_rnd_sign();
    }
    return eMat;
  }
  if (choice == 1) {
    MyMatrix<Tint> eMat = IdentityMat<Tint>(n);
    int i = rand() % n;
    int j = rand() % n;
    eMat(i, j) = get_rnd_sign();
    return eMat;
  }
  std::cerr << "Failed to find a matching entry in GetRandomMatrixPerturbation "
               "choice="
            << choice << "\n";
  throw TerminalException{1};
}

template <typename T, typename Tint>
MyVector<Tint> INDEFINITE_GetShortPositiveVector(MyMatrix<T> const &M,
                                                 std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: INDEFINITE_GetShortPositiveVector: beginning\n";
#endif
  int n = M.rows();
  auto L1_norm = [&](MyMatrix<T> const &M) -> T {
    T sum = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        sum += T_abs(M(i, j));
      }
    }
    return sum;
  };
  std::vector<MyVector<Tint>> GraverBasis = GetGraverKbasis<Tint>(n, 2, os);
  auto GetAlpha = [&](MyVector<Tint> const &TheVect,
                      MyVector<Tint> const &DirVect) -> int {
    T Norm1 = EvaluationQuadForm<T, Tint>(M, TheVect);
    int alpha = 0;
    while (true) {
      MyVector<Tint> NewVect = TheVect + (alpha + 1) * DirVect;
      T Norm2 = EvaluationQuadForm<T, Tint>(M, NewVect);
      if (Norm2 < Norm1 && Norm2 > 0) {
        Norm1 = Norm2;
        alpha += 1;
      } else {
        return alpha;
      }
    }
  };
  auto DirectImprovement =
      [&](MyVector<Tint> const &TheVect) -> MyVector<Tint> {
    MyVector<Tint> ImpVect = TheVect;
    while (true) {
      int n_chg = 0;
      for (auto &eVectBasis : GraverBasis) {
        int alpha = GetAlpha(ImpVect, eVectBasis);
        ImpVect += alpha * eVectBasis;
        n_chg += alpha;
      }
      if (n_chg == 0) {
        return ImpVect;
      }
    }
  };
  MyMatrix<Tint> ePerturb = IdentityMat<Tint>(n);
  int n_no_improv = 0;
  std::optional<T> CurrNorm;
  std::optional<MyVector<Tint>> CurrVect;
  while (true) {
#ifdef DEBUG_POSITIVITY
    os << "POS: INDEFINITE_GetShortPositiveVector: Loop step 1\n";
#endif
    MyMatrix<T> ePerturb_T = UniversalMatrixConversion<T, Tint>(ePerturb);
    MyMatrix<T> M_Pert = ePerturb_T * M * ePerturb_T.transpose();
    MyVector<Tint> uVect_Pert =
        GetIntegralPositiveVector_allmeth<T, Tint>(M_Pert, os);
#ifdef DEBUG_POSITIVITY
    os << "POS: INDEFINITE_GetShortPositiveVector: Loop step 2, we have "
          "uVect_Pert\n";
#endif
    MyVector<Tint> uVect = ePerturb.transpose() * uVect_Pert;
    MyVector<Tint> TheVect = DirectImprovement(uVect);
    T eNorm = EvaluationQuadForm<T, Tint>(M, TheVect);
#ifdef DEBUG_POSITIVITY
    os << "POS: INDEFINITE_GetShortPositiveVector: Loop step 3, eNorm=" << eNorm
       << "\n";
#endif
    auto is_lower = [&]() -> bool {
      if (CurrNorm) {
        return eNorm < *CurrNorm;
      }
      return true;
    };
    bool test_low = is_lower();
#ifdef DEBUG_POSITIVITY
    os << "POS: INDEFINITE_GetShortPositiveVector: Loop step 4, test_low="
       << test_low << "\n";
#endif
    if (test_low) {
      n_no_improv = 0;
      CurrVect = TheVect;
      CurrNorm = eNorm;
    } else {
      n_no_improv += 1;
    }
    if (n_no_improv == 500 || CurrNorm < 10) {
      MyVector<Tint> Vret = *CurrVect;
#ifdef DEBUG_POSITIVITY
      os << "POS: INDEFINITE_GetShortPositiveVector: returning Vret="
         << StringVector(Vret) << "\n";
#endif
      return Vret;
    }
    T norm1 = L1_norm(M_Pert);
    T norm2 = L1_norm(M);
#ifdef DEBUG_POSITIVITY
    os << "POS: INDEFINITE_GetShortPositiveVector: Loop step 5, norm1=" << norm1
       << " norm2=" << norm2 << "\n";
#endif
    if (norm1 > 1000 * norm2) {
      ePerturb = IdentityMat<Tint>(n);
    }
    ePerturb = ePerturb * GetRandomMatrixPerturbation<Tint>(n);
  }
}

template <typename T>
std::vector<MyVector<T>>
GetSetNegativeOrZeroVector(MyMatrix<T> const &SymMat,
                           std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: GetSetNegativeOrZeroVector: beginning\n";
#endif
  DiagSymMat<T> eRecDiag = DiagonalizeSymmetricMatrix(SymMat, os);
  std::vector<MyVector<T>> TheSet;
  int n = SymMat.rows();
  for (int i = 0; i < n; i++)
    if (eRecDiag.RedMat(i, i) <= 0) {
      MyVector<T> eVect = ZeroVector<T>(n);
      eVect(i) = 1;
      MyVector<T> fVect = (eRecDiag.Transform.transpose()) * eVect;
#ifdef DEBUG_POSITIVITY
      T eEval = EvaluationQuadForm(SymMat, fVect);
      if (eEval > 0) {
        std::cerr << "Big bad error\n";
        throw TerminalException{1};
      }
#endif
      T eMax = fVect.maxCoeff();
      MyVector<T> gVect = fVect / eMax;
      TheSet.push_back(gVect);
    }
  return TheSet;
}

template <typename T, typename Tint>
MyVector<Tint> GetShortVectorSpecified(MyMatrix<T> const &M,
                                       std::vector<MyVector<T>> const &ListVect,
                                       T const &MaxNorm,
                                       [[maybe_unused]] std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: GetShortVectorSpecified: beginning\n";
#endif
  int n = M.rows();
  Tint eMult = 1;
  while (true) {
    for (auto &eVect : ListVect) {
      MyVector<Tint> V(n);
      for (int i = 0; i < n; i++) {
        T eNear = UniversalNearestScalarInteger<T, T>(eMult * eVect(i));
        Tint eNear_i = UniversalScalarConversion<Tint, T>(eNear);
        V(i) = eNear_i;
      }
      T eVal = EvaluationQuadForm(M, V);
      if (eVal < MaxNorm)
        return V;
    }
    eMult++;
  }
}

template <typename T, typename Tint>
MyVector<Tint> GetShortVector(MyMatrix<T> const &M, T const &MaxNorm,
                              std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: GetShortVector: beginning\n";
#endif
  std::vector<MyVector<T>> ListNeg = GetSetNegativeOrZeroVector(M, os);
  return GetShortVectorSpecified<T, Tint>(M, ListNeg, MaxNorm, os);
}

template <typename T, typename Tint>
MyVector<Tint> GetShortVectorDegenerate(MyMatrix<T> const &M, T const &MaxNorm,
                                        std::ostream &os) {
#ifdef DEBUG_POSITIVITY
  os << "POS: GetShortVectorDegenerate: beginning\n";
#endif
  MyMatrix<T> NSP = NullspaceMat(M);
  int nbRow = NSP.rows();
  std::vector<MyVector<T>> ListVect(nbRow);
  for (int i = 0; i < nbRow; i++) {
    MyVector<T> eVect = GetMatrixRow(NSP, i);
    T eMax = eVect.maxCoeff();
    ListVect[i] = eVect / eMax;
  }
  return GetShortVectorSpecified<T, Tint>(M, ListVect, MaxNorm, os);
}

template <typename T> std::vector<T> GetLineVector(MyMatrix<T> const &M) {
  int n = M.rows();
  int dim = n * (n + 1) / 2;
  std::vector<T> V(dim);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      if (i == j) {
        V[pos] = M(i, i);
      } else {
        V[pos] = 2 * M(i, j);
      }
      pos++;
    }
  }
  return V;
}

template <typename T>
T EvaluateLineVector(std::vector<T> const &V_mat, MyVector<T> const &V) {
  T sum = 0;
  T pSum = 0;
  int n = V.size();
  size_t pos = 0;
  for (int i = 0; i < n; i++) {
    pSum = 0;
    for (int j = i; j < n; j++) {
      pSum += V_mat[pos] * V[j];
      pos++;
    }
    sum += pSum * V[i];
  }
  return sum;
}

template <typename T>
T EvaluateLineVectorShift(std::vector<T> const &V_mat, MyVector<T> const &V,
                          int shift) {
  T sum = 0;
  T pSum = 0;
  int n = V.size();
  size_t pos = 0;
  for (int i = shift; i < n; i++) {
    pSum = 0;
    for (int j = i; j < n; j++) {
      pSum += V_mat[pos] * V[j];
      pos++;
    }
    sum += pSum * V[i];
  }
  return sum;
}

// clang-format off
#endif  // SRC_LATT_POSITIVITY_H_
// clang-format on
