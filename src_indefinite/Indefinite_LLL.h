#ifndef SRC_INDEFINITE_INDEFINITE_LLL_H_
#define SRC_INDEFINITE_INDEFINITE_LLL_H_

#include "MAT_Matrix.h"
#include <algorithm>
#include <utility>
#include <vector>

//#define DEBUG_INDEFINITE_LLL

template <typename T> struct ResultGramSchmidt_Indefinite {
  bool success; // true means we have a basis. False that we have an isotropic
                // vector
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
  std::cerr << "det1=" << det1 << " det2=" << det2 << "\n";
#endif
  return {true, Bstar, Bstar_norms, mu, {}};
}

template <typename T, typename Tint> struct ResultIndefiniteLLL {
  bool success; // true if we obtained the reduced matrix. false if we found an
                // isotropic vector
  MyMatrix<Tint> B;
  MyMatrix<T> Mred;
  MyVector<T> Xisotrop;
};

// Adapted from Denis Simon, Solving Quadratic Equations Using Reduced
// Unimodular Quadratic Forms, Math. Comp. 74(251) 1531--1543
template <typename T, typename Tint>
ResultIndefiniteLLL<T, Tint> Indefinite_LLL(MyMatrix<T> const &M) {
  int n = M.rows();
  T c = T(7) / T(8); // The c constant of the LLL algorithm
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
    std::cerr << "Passing in Indefinite_LLL det=" << det << " k=" << k << "\n";
#endif
    ResultGramSchmidt_Indefinite<T> ResGS = GramSchmidtOrthonormalization(M, B);
#ifdef DEBUG_INDEFINITE_LLL
    std::cerr << " Bstar_norms =";
    for (auto &eN : ResGS.Bstar_norms)
      std::cerr << " " << eN;
    std::cerr << "\n";
#endif
    if (!ResGS.success) {
      return {false, B, get_matrix(), ResGS.Xisotrop};
    }
    for (int i = n - 1; i >= 0; i--) {
      for (int j = 0; j < i; j++) {
        T val = ResGS.mu(i, j);
        //        double val_d = UniversalScalarConversion<double,T>(val);
        Tint q = UniversalNearestScalarInteger<Tint, T>(val);
        //        std::cerr << " ResGS.mu=" << ResGS.mu(i,j) << " val_d=" <<
        //        val_d << " q=" << q << "\n"; std::cerr << "i=" << i << "
        //        B.row(i)=" << GetMatrixRow(B, i) << " B.row(j)=" <<
        //        GetMatrixRow(B,j) << "\n";
        B.row(i) -= q * B.row(j);
        //        std::cerr << "   q=" << q << " B.row(i)=" << GetMatrixRow(B,i)
        //        << "\n";
      }
    }
    T mu = ResGS.mu(k, k - 1);
    T sum1_pre = ResGS.Bstar_norms[k] + mu * mu * ResGS.Bstar_norms[k - 1];
    T sum1 = T_abs(sum1_pre);
    T sum2 = c * T_abs(ResGS.Bstar_norms[k - 1]);
    if (sum1 < sum2) {
#ifdef DEBUG_INDEFINITE_LLL
      std::cerr << "Swapping k=" << k << " and " << (k - 1) << "\n";
#endif
      for (int i = 0; i < n; i++)
        std::swap(B(k, i), B(k - 1, i));
      k = std::max(k - 1, 1);
    } else {
      k++;
    }
    if (k >= n)
      break;
  }
  return {true, B, get_matrix(), {}};
}

template <typename T, typename Tint> struct ResultReductionIndefinite {
  MyMatrix<Tint> B;
  MyMatrix<T> Mred;
};

template <typename T, typename Tint>
ResultReductionIndefinite<T, Tint>
ComputeReductionIndefinite(MyMatrix<T> const &M) {
  std::cerr << "Beginning of ComputeReductionIndefinite\n";
  int n = M.rows();
  ResultIndefiniteLLL<T, Tint> eRes = Indefinite_LLL<T, Tint>(M);
  std::cerr << "We have computed eRes\n";
  bool early_term = false;
  if (eRes.success && early_term) {
    return {std::move(eRes.B), std::move(eRes.Mred)};
  }
  MyMatrix<Tint> B = eRes.B;
  MyMatrix<T> Mwork = eRes.Mred;
  auto get_norm = [&](MyMatrix<T> const &mat) -> T {
    T sum = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        sum += T_abs(mat(i, j));
    return sum;
  };
  auto get_random_int_matrix = [&]() -> MyMatrix<Tint> {
    std::vector<int> LPos(n);
    for (int i = 0; i < n; i++)
      LPos[i] = i;
    for (int iter = 0; iter < 4 * n; iter++) {
      int i = rand() % n;
      int j = rand() % n;
      if (i != j)
        std::swap(LPos[i], LPos[j]);
    }
    std::vector<int> LDiag(n);
    for (int i = 0; i < n; i++) {
      int val = rand() % 2;
      LDiag[i] = -1 + 2 * val;
    }
    MyMatrix<Tint> Unit = ZeroMatrix<Tint>(n, n);
    for (int i = 0; i < n; i++)
      Unit(i, LPos[i]) = LDiag[i];
    return Unit;
  };
  T norm_work = get_norm(Mwork);
  size_t iter_no_improv = 0;
  size_t limit_iter = 2 * n;
  while (true) {
    MyMatrix<Tint> RandUnit = get_random_int_matrix();
    //    std::cerr << "RandUnit=\n";
    //    WriteMatrix(std::cerr, RandUnit);
    MyMatrix<T> RandUnit_T = UniversalMatrixConversion<T, Tint>(RandUnit);
    B = RandUnit * B;
    Mwork = RandUnit_T * Mwork * RandUnit_T.transpose();
    ResultIndefiniteLLL<T, Tint> eRes = Indefinite_LLL<T, Tint>(Mwork);
    if (eRes.success && early_term) {
      B = eRes.B * B;
      Mwork = eRes.Mred;
      return {std::move(B), std::move(Mwork)};
    }
    T norm = get_norm(eRes.Mred);
    std::cerr << "norm=" << norm << " norm_work=" << norm_work
              << " iter_no_improv=" << iter_no_improv << "\n";
    if (norm >= norm_work) {
      iter_no_improv++;
      if (limit_iter == iter_no_improv)
        return {std::move(B), std::move(Mwork)};
    } else {
      iter_no_improv = 0;
      norm_work = norm;
      B = eRes.B * B;
      Mwork = eRes.Mred;
    }
  }
}

template <typename T, typename Tint>
ResultReductionIndefinite<T, Tint>
ComputeReductionIndefinite_opt(MyMatrix<T> const &M,
                               bool const &ApplyReduction) {
  if (ApplyReduction) {
    return ComputeReductionIndefinite<T, Tint>(M);
  } else {
    int n = M.rows();
    MyMatrix<Tint> B = IdentityMat<Tint>(n);
    MyMatrix<T> Mred = M;
    return {B, Mred};
  }
}

#endif //  SRC_INDEFINITE_INDEFINITE_LLL_H_
