// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_TSPACE_SEARCH_POSITIVE_MATRICES_H_
#define SRC_TSPACE_SEARCH_POSITIVE_MATRICES_H_

// clang-format off
#include "POLY_Fundamental.h"
#include "Positivity.h"
#include "SHORT_Realizability.h"
#include <vector>
// clang-format on

/*
  Find one positive definite matrix in the space assuming that one exists.
  If one of the matrix of the basis is positive definite then the first one
  is provided.
  ---
  Maybe we should use semidefinite programming? Is that absolutely
  needed? Or just slower?
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GetOnePositiveDefiniteMatrix_ListV(std::vector<MyMatrix<T>> const &ListMat,
                                   std::vector<MyVector<Tint>> const& ListV_init,
                                   std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr
        << "TSPACE: The number of matrices is 0 so we cannot build a positive "
           "definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPFCT: GetOnePositiveDefiniteMatrix n_mat=" << n_mat << " n=" << n << "\n";
#endif
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> const &eMat = ListMat[i_mat];
    if (IsPositiveDefinite(eMat, os)) {
      return eMat;
    }
  }
  std::vector<MyVector<Tint>> ListV = ListV_init;
#ifdef DEBUG_TSPACE_FUNCTIONS
  int iter=0;
#endif
  while (true) {
    int n_vect = ListV.size();
#ifdef DEBUG_TSPACE_FUNCTIONS
    iter += 1;
    os << "TSPFCT: GetOnePositiveDefiniteMatrix iter=" << iter << " n_vect=" << n_vect << "\n";
#endif
    MyMatrix<T> ListIneq = ZeroMatrix<T>(n_vect, 1 + n_mat);
    MyVector<T> ToBeMinimized = ZeroVector<T>(1 + n_mat);
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      ListIneq(i_vect, 0) = -1;
      MyVector<Tint> V = ListV[i_vect];
      MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
      for (int i_mat = 0; i_mat < n_mat; i_mat++) {
        T val = EvaluationQuadForm(ListMat[i_mat], V_T);
        ListIneq(i_vect, 1 + i_mat) = val;
        ToBeMinimized(1 + i_mat) += val;
      }
    }
    //
    // Solving the linear program
    //
    LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, ToBeMinimized, os);
    if (!eSol.PrimalDefined || !eSol.DualDefined) {
      return {};
    }
    MyMatrix<T> TrySuperMat = ZeroMatrix<T>(n, n);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      TrySuperMat += eSol.DirectSolution(i_mat) * ListMat[i_mat];
    }
    if (IsPositiveDefinite(TrySuperMat, os)) {
      return TrySuperMat;
    }
    //
    // Failed, trying to find a vector
    //
    int rnk = RankMat(TrySuperMat);
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPFCT: GetOnePositiveDefiniteMatrix get_one_vect, rnk=" << rnk << " n=" << n << "\n";
    os << "TSPFCT: TrySuperMat=\n";
    WriteMatrix(os, TrySuperMat);
#endif
    if (rnk < n) {
      std::vector<MyVector<Tint>> list_v = GetShortVectorDegenerate<T, Tint>(TrySuperMat, os);
      for (auto & eV: list_v) {
        ListV.push_back(eV);
      }
    } else {
      T CritNorm(0);
      bool StrictIneq = false;
      MyVector<Tint> V = GetShortIntegralVector<T, Tint>(TrySuperMat, CritNorm, StrictIneq, os);
      ListV.push_back(V);
    }
  }
}

template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GetOnePositiveDefiniteMatrix(std::vector<MyMatrix<T>> const &ListMat,
                             std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr
        << "TSPACE: The number of matrices is 0 so we cannot build a positive "
           "definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  std::vector<MyVector<Tint>> ListV_init = get_initial_vector_test_v<Tint>(n, {}, os);
  return GetOnePositiveDefiniteMatrix_ListV(ListMat, ListV_init, os);
}

/*
  Tries to find a semi-definite matrix.
  Unfortunately, we tend to go into infinite loops.
  And that seems to be an intrinsic problem since
  S_{n,>0} is not a polyhedral cone.
  ---
  But it works sometimes.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GetOnePositiveSemiDefiniteMatrix_ListV(std::vector<MyMatrix<T>> const &ListMat,
                                       std::vector<MyVector<Tint>> const& ListV_init,
                                       std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr
        << "TSPACE: The number of matrices is 0 so we cannot build a positive "
           "definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix n_mat=" << n_mat << " n=" << n << "\n";
#endif
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> const &eMat = ListMat[i_mat];
    if (IsPositiveSemiDefinite(eMat, os)) {
      // A useful ansatz
      return eMat;
    }
  }
  std::vector<MyVector<Tint>> ListV = ListV_init;
  int n_vect_init = ListV.size();
#ifdef DEBUG_TSPACE_FUNCTIONS
  int iter=0;
#endif
  while (true) {
    int n_vect = ListV.size();
#ifdef DEBUG_TSPACE_FUNCTIONS
    iter += 1;
    os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix iter=" << iter << " n_vect=" << n_vect << "\n";
    os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix n_mat=" << n_mat << " n=" << n << "\n";
#endif
    MyMatrix<T> ListIneq = ZeroMatrix<T>(n_vect + 1, 1 + n_mat);
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      MyVector<Tint> const& V = ListV[i_vect];
      MyVector<T> V_T = UniversalVectorConversion<T, Tint>(V);
      for (int i_mat = 0; i_mat < n_mat; i_mat++) {
        T val = EvaluationQuadForm(ListMat[i_mat], V_T);
        ListIneq(i_vect, 1 + i_mat) = val;
      }
    }
    ListIneq(n_vect, 0) = -1;
    MyVector<T> ToBeMinimized = ZeroVector<T>(1 + n_mat);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      T sum(0);
      for (int i_vect=0; i_vect<n_vect_init; i_vect++) {
        sum += ListIneq(i_vect, 1 + i_mat);
      }
      ListIneq(n_vect, 1 + i_mat) = sum;
      ToBeMinimized(1 + i_mat) = sum;
    }
    //
    // Solving the linear program
    //
    LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, ToBeMinimized, os);
    if (!eSol.PrimalDefined || !eSol.DualDefined) {
      return {};
    }
    MyMatrix<T> TrySuperMat = ZeroMatrix<T>(n, n);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      TrySuperMat += eSol.DirectSolution(i_mat) * ListMat[i_mat];
    }
    if (IsPositiveSemiDefinite(TrySuperMat, os)) {
      return TrySuperMat;
    }
    //
    // Failed, trying to find a vector
    //
    MyMatrix<T> NSP_T = NullspaceIntMat(TrySuperMat);
    MyMatrix<Tint> NSP = UniversalMatrixConversion<Tint,T>(NSP_T);
    MyMatrix<Tint> Compl = SubspaceCompletionInt<Tint>(NSP, n);
    MyMatrix<T> Compl_T = UniversalMatrixConversion<T,Tint>(Compl);
    MyMatrix<T> M = Compl_T * TrySuperMat * Compl_T.transpose();
    T CritNorm(0);
    bool StrictIneq = true;
    MyVector<Tint> V1 = GetShortIntegralVector<T, Tint>(M, CritNorm, StrictIneq, os);
    MyVector<Tint> V2 = Compl.transpose() * V1;
#ifdef DEBUG_TSPACE_FUNCTIONS
    os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix V2=" << StringVectorGAP(V2) << "\n";
    os << "TSPFCT: GetOnePositiveSemiDefiniteMatrix TrySuperMat=\n";
    WriteMatrixGAP(os, TrySuperMat);
    os << "\n";
#endif
#ifdef SANITY_CHECK_TSPACE_FUNCTIONS
    T norm = EvaluationQuadForm(TrySuperMat, V2);
    if (norm >= 0) {
      std::cerr << "TSPFCT: The norm should be negative since it should be a counter example\n";
      throw TerminalException{1};
    }
#endif
    ListV.push_back(V2);
  }
}

template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GetOnePositiveSemiDefiniteMatrix(std::vector<MyMatrix<T>> const &ListMat,
                                 std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr
        << "TSPACE: The number of matrices is 0 so we cannot build a positive "
           "definite matrix\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  std::vector<MyVector<Tint>> ListV_init = get_initial_vector_test_v<Tint>(n, {}, os);
  return GetOnePositiveSemiDefiniteMatrix_ListV(ListMat, ListV_init, os);
}





// clang-format off
#endif  // SRC_TSPACE_SEARCH_POSITIVE_MATRICES_H_
// clang-format on

