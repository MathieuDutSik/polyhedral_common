// Copyright (C) 2023 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_LISTMATSATURATION_H_
#define SRC_LATT_TSPACE_LISTMATSATURATION_H_

template <typename T>
bool is_integrally_saturated_matrix_space(
    std::vector<MyMatrix<T>> const &ListMat) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr << "TSPACE: We have n_mat=0\n";
    std::cerr
        << "TSPACE: The code could work with n_mat=0 but we are not sure it "
           "makes sense\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  int sym_dim = (n + 1) * n / 2;
  if (n_mat == sym_dim) {
    std::cerr << "TSPACE: We have n_mat=" << n_mat
              << " equal to sym_dim=" << sym_dim << "\n";
    throw TerminalException{1};
  }
  MyMatrix<T> BigMat(n_mat, sym_dim);
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    if (!IsIntegralMatrix(ListMat[i_mat])) {
      return false;
    }
    int pos = 0;
    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        BigMat(i_mat, pos) = ListMat[i_mat](i, j);
        pos++;
      }
    }
  }
  MyMatrix<T> NSP1 = NullspaceIntTrMat(BigMat);
  MyMatrix<T> BigMat_renorm = NullspaceIntTrMat(NSP1);
  int dim_space = BigMat_renorm.rows();
  MyMatrix<T> sol_mat(dim_space, dim_space);
  for (int i = 0; i < dim_space; i++) {
    MyVector<T> V = GetMatrixRow(BigMat, i);
    std::optional<MyVector<T>> opt = SolutionMat(BigMat_renorm, V);
    if (!opt) {
      std::cerr << "TSP_FCT: IntegralSaturationSpace, no solution at i=" << i
                << "\n";
      throw TerminalException{1};
    }
    MyVector<T> const &V2 = *opt;
    AssignMatrixRow(sol_mat, i, V2);
  }
  if (!IsIntegralMatrix(sol_mat)) {
    std::cerr << "TSP_FCT: IntegralSaturationSpace, The matrix sol_mat should "
                 "be integral\n";
    throw TerminalException{1};
  }
  T det = T_abs(DeterminantMat(sol_mat));
  if (det == 1) {
    return true;
  } else {
    return false;
  }
}

/*
  We need the space of matrices to be spanning the integral saturation so that
  the elements of GlStab are integral.
 */
template <typename T>
std::vector<MyMatrix<T>>
IntegralSaturationSpace(std::vector<MyMatrix<T>> const &ListMat,
                        [[maybe_unused]] std::ostream &os) {
  int n_mat = ListMat.size();
  if (n_mat == 0) {
    std::cerr << "TSPACE: We have n_mat=0\n";
    std::cerr
        << "TSPACE: The code could work with n_mat=0 but we are not sure it "
           "makes sense\n";
    throw TerminalException{1};
  }
  int n = ListMat[0].rows();
  int sym_dim = (n + 1) * n / 2;
  if (n_mat == sym_dim) {
    std::cerr << "TSPACE: We have n_mat=" << n_mat
              << " equal to sym_dim=" << sym_dim << "\n";
    throw TerminalException{1};
  }
  MyMatrix<T> BigMat(n_mat, sym_dim);
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    int pos = 0;
    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        BigMat(i_mat, pos) = ListMat[i_mat](i, j);
        pos++;
      }
    }
  }
  MyMatrix<T> NSP1 = NullspaceIntTrMat(BigMat);
  MyMatrix<T> BigMat_renorm = NullspaceIntTrMat(NSP1);
#ifdef SANITY_CHECK_TSPACE_FUNCTIONS
  int dim_space = BigMat_renorm.rows();
  MyMatrix<T> sol_mat(dim_space, dim_space);
  for (int i = 0; i < dim_space; i++) {
    MyVector<T> V = GetMatrixRow(BigMat, i);
    std::optional<MyVector<T>> opt = SolutionMat(BigMat_renorm, V);
    if (!opt) {
      std::cerr << "TSP_FCT: IntegralSaturationSpace, no solution at i=" << i
                << "\n";
      throw TerminalException{1};
    }
    MyVector<T> const &V2 = *opt;
    AssignMatrixRow(sol_mat, i, V2);
  }
  if (!IsIntegralMatrix(sol_mat)) {
    std::cerr << "TSP_FCT: BigMat=\n";
    WriteMatrix(std::cerr, BigMat);
    std::cerr << "TSP_FCT: BigMat_renorm=\n";
    WriteMatrix(std::cerr, BigMat_renorm);
    std::cerr << "TSP_FCT: sol_mat=\n";
    WriteMatrix(std::cerr, sol_mat);
    std::cerr << "TSP_FCT: IntegralSaturationSpace, The matrix sol_mat should "
                 "be integral\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSP_FCT: IntegralSaturationSpace, det=" << DeterminantMat(sol_mat)
     << "\n";
#endif
#endif
  if (BigMat_renorm.rows() != n_mat) {
    std::cerr << "TSPACE: Incoherence in the dimensions\n";
    throw TerminalException{1};
  }
  std::vector<MyMatrix<T>> ListMatRet;
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> eMat(n, n);
    int pos = 0;
    for (int i = 0; i < n; i++) {
      for (int j = i; j < n; j++) {
        T val = BigMat_renorm(i_mat, pos);
        eMat(i, j) = val;
        eMat(j, i) = val;
        pos++;
      }
    }
    ListMatRet.push_back(eMat);
  }
  return ListMatRet;
}

// clang-format off
#endif  // SRC_LATT_TSPACE_LISTMATSATURATION_H_
// clang-format on
