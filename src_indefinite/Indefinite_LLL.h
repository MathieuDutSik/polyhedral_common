#ifndef SRC_INDEFINITE_INDEFINITE_LLL_H_
#define SRC_INDEFINITE_INDEFINITE_LLL_H_

#include "MAT_Matrix.h"
#include "GRAPH_GraphicalBasic.h"
#include "GRAPH_BitsetType.h"
#include "WeightMatrix.h"
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

template<typename T>
std::pair<MyMatrix<T>,MyMatrix<T>> CanonicalizationPermutationSigns(MyMatrix<T> const& M)
{
  using Tidx_value=uint16_t;
  using Tidx=uint16_t;
  using Tgr = GraphListAdj;
  Tidx n = M.rows();
  MyMatrix<T> Mabs(n, n);
  for (Tidx i_row=0; i_row<n; i_row++)
    for (Tidx i_col=0; i_col<n; i_col++)
      Mabs(i_row,i_col) = T_abs(M(i_row,i_col));
  //  std::cerr << "M=" << M << "\n";
  WeightMatrix<true, T, Tidx_value> WMat = WeightedMatrixFromMyMatrix<true, T, Tidx_value>(Mabs);
  WMat.ReorderingSetWeight();
  std::vector<Tidx> CanonicOrd = GetGroupCanonicalizationVector_Kernel<T, Tgr, Tidx, Tidx_value>(WMat).first;
  //  std::cerr << "We have CanonicOrd\n";
  //  for (Tidx i=0; i<n; i++)
  //    std::cerr << "i=" << i << " val=" << CanonicOrd[i] << "\n";
  MyMatrix<T> Mreord(n,n);
  MyMatrix<T> MreordAbs(n,n);
  for (Tidx i_row=0; i_row<n; i_row++) {
    Tidx j_row = CanonicOrd[i_row];
    for (Tidx i_col=0; i_col<n; i_col++) {
      Tidx j_col = CanonicOrd[i_col];
      Mreord(i_row, i_col) = M(j_row, j_col);
      MreordAbs(i_row, i_col) = T_abs(M(j_row, j_col));
    }
  }
  //  std::cerr << "Mreord=\n";
  //  WriteMatrix(std::cerr, Mreord);
  //  std::cerr << "MreordAbs=\n";
  //  WriteMatrix(std::cerr, MreordAbs);
  MyMatrix<T> Mtrans1 = ZeroMatrix<T>(n,n);
  for (Tidx i_row=0; i_row<n; i_row++) {
    Tidx j_row = CanonicOrd[i_row];
    Mtrans1(i_row,j_row) = 1;
  }
  //  std::cerr << "We have Mtrans1, Mreord\n";
  MyMatrix<T> eProd = Mtrans1 * M * Mtrans1.transpose();
  if (eProd != Mreord) {
    std::cerr << "The matrix product does not work as expected\n";
    std::cerr << "eProd=\n";
    WriteMatrix(std::cerr, eProd);
    throw TerminalException{1};
  }

  GraphBitset GR(n);
  for (Tidx i_row=0; i_row<n; i_row++) {
    for (Tidx i_col=0; i_col<n; i_col++) {
      if (Mreord(i_row,i_col) != 0) {
        GR.AddAdjacent(i_row,i_col);
        GR.AddAdjacent(i_col,i_row);
      }
    }
  }
  MyMatrix<T> Mtrans2 = IdentityMat<T>(n);
  std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(GR);
  //  std::cerr << "|LConn|=" << LConn.size() << "\n";
  for (auto & eConn : LConn) {
    /*
    std::cerr << "eConn =";
    for (auto & eVal : eConn)
      std::cerr << " " << eVal;
    std::cerr << "\n";
    */
    size_t len = eConn.size();
    std::vector<size_t> Status(len, 0);
    std::vector<size_t> eConnRev(n, std::numeric_limits<size_t>::max());
    for (size_t i=0; i<len; i++)
      eConnRev[eConn[i]] = i;
    Status[0] = 1;
    while(true) {
      size_t n_done = 0;
      for (size_t i=0; i<len; i++)
        if (Status[i] > 0)
          n_done++;
      if (n_done == len)
        break;
      for (size_t i=0; i<len; i++) {
        if (Status[i] > 0) {
          size_t iImg = eConn[i];
          for (auto & eAdjImg : GR.Adjacency(iImg)) {
            size_t eAdj = eConnRev[eAdjImg];
            if (Status[eAdj] == 0) {
              int sign = 1;
              if (Mreord(iImg,eAdjImg) < 0)
                sign = -1;
              Mtrans2(eAdjImg,eAdjImg) = sign * Mtrans2(iImg,iImg);
              Status[eAdj]=1;
            }
          }
        }
      }
    }
  }
  //  std::cerr << "We have Mtrans2=\n";
  //  WriteMatrix(std::cerr, Mtrans2);
  MyMatrix<T> eP = Mtrans2 * Mtrans1;
  MyMatrix<T> M_red = eP * M * eP.transpose();
  return {eP, M_red};
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
