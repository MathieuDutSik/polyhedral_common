// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ERDAHL_ERDAHL_SPACE_H_
#define SRC_ERDAHL_ERDAHL_SPACE_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <optional>
#include <string>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_ERDAHL_SPACE
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_ERDAHL_SPACE
#endif

/*
  Degree 2 functions on Z^n and the linear spaces of them.

  Conventions, used by all of src_erdahl:
  * A point x of Z^n is the row vector e = (1, x) of length n+1, a direction
    z of Z^n is the row vector (0, z).
  * A function of degree at most 2
        f(x) = Cst(f) + 2 Lin(f).x + Quad(f)[x]
    is the symmetric (n+1) x (n+1) matrix
        F = [ Cst(f)  Lin(f)^T ]
            [ Lin(f)  Quad(f)  ]
    so that f(x) = e F e^T.
  * An affine integral transformation is an integral (n+1) x (n+1) matrix g
    of determinant +-1 whose first column is (1, 0, ..., 0)^T, acting on the
    rows by e -> e g. The function f o g is then g F g^T, and the zero set
    of f o g is the image of the zero set of f under g^{-1}.

  The paper (Dutour Sikiric, Erdahl, "Enumeration of perfect Delaunay
  polytopes") works in the full space E_2(n) of such functions, of dimension
  (n+1)(n+2)/2. Here the computations are done in a linear subspace W of
  E_2(n), given by a basis. The cone of interest is then Erdahl(n) cap W,
  whose faces are the zero sets of its elements. The natural choices are:
  * The full space E_2(n).
  * The space of functions of center c: f(x) = a + Q[x - c]. For a
    half-integral c these are the functions invariant under the point
    reflection x -> 2c - x, and their zero sets are the centrally symmetric
    Delaunay polyhedra of center c.
  * More generally, the space of functions invariant under a finite group
    of affine integral transformations.
 */

template <typename T> struct ErdahlFunctionSpace {
  int n;
  // Linearly independent symmetric (n+1) x (n+1) matrices.
  std::vector<MyMatrix<T>> basis;
  // A short human readable description, e.g. "full", "centered".
  std::string name;
  // For the space of functions of center c, the center. The affine
  // integral transformations preserving that space are then exactly those
  // fixing c.
  std::optional<MyVector<T>> center;
};

template <typename T> MyMatrix<T> erdahl_get_quad(MyMatrix<T> const &F) {
  int n = F.rows() - 1;
  MyMatrix<T> Q(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Q(i, j) = F(i + 1, j + 1);
    }
  }
  return Q;
}

template <typename T> MyVector<T> erdahl_get_lin(MyMatrix<T> const &F) {
  int n = F.rows() - 1;
  MyVector<T> V(n);
  for (int i = 0; i < n; i++) {
    V(i) = F(0, i + 1);
  }
  return V;
}

// The function F from its three components.
template <typename T>
MyMatrix<T> erdahl_assemble_function(T const &cst, MyVector<T> const &lin,
                                     MyMatrix<T> const &quad) {
  int n = quad.rows();
  MyMatrix<T> F(n + 1, n + 1);
  F(0, 0) = cst;
  for (int i = 0; i < n; i++) {
    F(0, i + 1) = lin(i);
    F(i + 1, 0) = lin(i);
    for (int j = 0; j < n; j++) {
      F(i + 1, j + 1) = quad(i, j);
    }
  }
  return F;
}

// e1 F e2^T for two rows of length n+1 (points or directions).
template <typename T, typename Tint>
T erdahl_bilinear(MyMatrix<T> const &F, MyVector<Tint> const &e1,
                  MyVector<Tint> const &e2) {
  int len = F.rows();
  T sum(0);
  for (int i = 0; i < len; i++) {
    if (e1(i) == 0) {
      continue;
    }
    T part(0);
    for (int j = 0; j < len; j++) {
      if (e2(j) != 0) {
        T e2_T = UniversalScalarConversion<T, Tint>(e2(j));
        part += F(i, j) * e2_T;
      }
    }
    T e1_T = UniversalScalarConversion<T, Tint>(e1(i));
    sum += e1_T * part;
  }
  return sum;
}

// The action f -> f o g on functions: F -> g F g^T.
template <typename T, typename Tint>
MyMatrix<T> erdahl_transform_function(MyMatrix<T> const &F,
                                      MyMatrix<Tint> const &g) {
  MyMatrix<T> g_T = UniversalMatrixConversion<T, Tint>(g);
  return g_T * F * g_T.transpose();
}

template <typename T>
MyMatrix<T> erdahl_basis_as_rows(ErdahlFunctionSpace<T> const &W) {
  int n_basis = W.basis.size();
  int dim_sym = (W.n + 1) * (W.n + 2) / 2;
  MyMatrix<T> M(n_basis, dim_sym);
  for (int i = 0; i < n_basis; i++) {
    MyVector<T> V = SymmetricMatrixToVector(W.basis[i]);
    AssignMatrixRow(M, i, V);
  }
  return M;
}

template <typename T>
void erdahl_check_space(ErdahlFunctionSpace<T> const &W) {
  int np1 = W.n + 1;
  for (auto &F : W.basis) {
    if (F.rows() != np1 || F.cols() != np1) {
      std::cerr << "ERDAHL: a basis element is not of size n+1\n";
      throw TerminalException{1};
    }
    if (F != F.transpose()) {
      std::cerr << "ERDAHL: a basis element is not symmetric\n";
      throw TerminalException{1};
    }
  }
  MyMatrix<T> M = erdahl_basis_as_rows(W);
  if (RankMat(M) != M.rows()) {
    std::cerr << "ERDAHL: the basis of the space is not linearly "
                 "independent\n";
    throw TerminalException{1};
  }
}

// The space spanned by symmetric matrices, extracting a basis.
template <typename T>
ErdahlFunctionSpace<T>
erdahl_space_from_spanning(int n, std::vector<MyMatrix<T>> const &l_mat,
                           std::string const &name) {
  int dim_sym = (n + 1) * (n + 2) / 2;
  MyMatrix<T> M(l_mat.size(), dim_sym);
  for (size_t i = 0; i < l_mat.size(); i++) {
    AssignMatrixRow(M, i, SymmetricMatrixToVector(l_mat[i]));
  }
  MyMatrix<T> Mred = RowReduction(M);
  std::vector<MyMatrix<T>> basis;
  for (int i = 0; i < Mred.rows(); i++) {
    MyVector<T> V = GetMatrixRow(Mred, i);
    basis.push_back(VectorToSymmetricMatrix(V, n + 1));
  }
  ErdahlFunctionSpace<T> W{n, basis, name, {}};
  erdahl_check_space(W);
  return W;
}

template <typename T> ErdahlFunctionSpace<T> erdahl_full_space(int n) {
  std::vector<MyMatrix<T>> basis;
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= i; j++) {
      MyMatrix<T> F = ZeroMatrix<T>(n + 1, n + 1);
      F(i, j) = 1;
      F(j, i) = 1;
      basis.push_back(F);
    }
  }
  return {n, basis, "full", {}};
}

/*
  The functions f(x) = a + Q[x - c]. Equivalently, those with
  Lin(f) = - Quad(f) c. For a half-integral c, these are exactly the
  functions invariant under the point reflection x -> 2c - x.
 */
template <typename T>
ErdahlFunctionSpace<T> erdahl_centered_space(MyVector<T> const &c) {
  int n = c.size();
  std::vector<MyMatrix<T>> basis;
  MyMatrix<T> Fcst = ZeroMatrix<T>(n + 1, n + 1);
  Fcst(0, 0) = 1;
  basis.push_back(Fcst);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      MyMatrix<T> Q = ZeroMatrix<T>(n, n);
      Q(i, j) = 1;
      Q(j, i) = 1;
      MyVector<T> Qc = Q * c;
      T cst = c.dot(Qc);
      MyVector<T> lin = -Qc;
      basis.push_back(erdahl_assemble_function(cst, lin, Q));
    }
  }
  return {n, basis, "centered", c};
}

template <typename Tint>
void erdahl_check_affine_transformation(MyMatrix<Tint> const &g) {
  int np1 = g.rows();
  if (g.cols() != np1) {
    std::cerr << "ERDAHL: the transformation should be square\n";
    throw TerminalException{1};
  }
  if (g(0, 0) != 1) {
    std::cerr << "ERDAHL: the transformation should have g(0,0) = 1\n";
    throw TerminalException{1};
  }
  for (int i = 1; i < np1; i++) {
    if (g(i, 0) != 0) {
      std::cerr << "ERDAHL: the transformation should have g(i,0) = 0 for "
                   "i > 0\n";
      throw TerminalException{1};
    }
  }
  Tint det = DeterminantMat(g);
  if (T_abs(det) != 1) {
    std::cerr << "ERDAHL: the transformation should be of determinant +-1\n";
    throw TerminalException{1};
  }
}

/*
  The functions F with g F g^T = F for all the transformations g of the
  list, i.e. the functions invariant under the group they generate.
 */
template <typename T, typename Tint>
ErdahlFunctionSpace<T>
erdahl_invariant_space(int n, std::vector<MyMatrix<Tint>> const &l_gen) {
  ErdahlFunctionSpace<T> Wfull = erdahl_full_space<T>(n);
  int dim_full = Wfull.basis.size();
  int dim_sym = dim_full;
  int n_gen = l_gen.size();
  MyMatrix<T> Equa(n_gen * dim_sym, dim_full);
  for (int i_gen = 0; i_gen < n_gen; i_gen++) {
    erdahl_check_affine_transformation(l_gen[i_gen]);
    for (int k = 0; k < dim_full; k++) {
      MyMatrix<T> const &E = Wfull.basis[k];
      MyMatrix<T> diff = erdahl_transform_function(E, l_gen[i_gen]) - E;
      MyVector<T> V = SymmetricMatrixToVector(diff);
      for (int u = 0; u < dim_sym; u++) {
        Equa(i_gen * dim_sym + u, k) = V(u);
      }
    }
  }
  MyMatrix<T> NSP = NullspaceTrMat(Equa);
  std::vector<MyMatrix<T>> basis;
  for (int i = 0; i < NSP.rows(); i++) {
    MyMatrix<T> F = ZeroMatrix<T>(n + 1, n + 1);
    for (int k = 0; k < dim_full; k++) {
      F += NSP(i, k) * Wfull.basis[k];
    }
    basis.push_back(RemoveFractionMatrix(F));
  }
  return {n, basis, "invariant", {}};
}

template <typename T>
std::optional<MyVector<T>>
erdahl_coordinates_in_space(ErdahlFunctionSpace<T> const &W,
                            MyMatrix<T> const &F) {
  MyMatrix<T> M = erdahl_basis_as_rows(W);
  MyVector<T> V = SymmetricMatrixToVector(F);
  return SolutionMat(M, V);
}

template <typename T>
bool erdahl_is_in_space(ErdahlFunctionSpace<T> const &W, MyMatrix<T> const &F) {
  return erdahl_coordinates_in_space(W, F).has_value();
}

template <typename T>
MyMatrix<T> erdahl_function_from_coefficients(ErdahlFunctionSpace<T> const &W,
                                              MyVector<T> const &coeff) {
  int np1 = W.n + 1;
  MyMatrix<T> F = ZeroMatrix<T>(np1, np1);
  for (size_t i = 0; i < W.basis.size(); i++) {
    if (coeff(i) != 0) {
      F += coeff(i) * W.basis[i];
    }
  }
  return F;
}

// Whether f -> f o g maps W onto W.
template <typename T, typename Tint>
bool erdahl_preserves_space(ErdahlFunctionSpace<T> const &W,
                            MyMatrix<Tint> const &g) {
  for (auto &F : W.basis) {
    MyMatrix<T> Fimg = erdahl_transform_function(F, g);
    if (!erdahl_is_in_space(W, Fimg)) {
      return false;
    }
  }
  return true;
}

// The vector (e1 W_i e2^T)_i of the values of the basis of W.
template <typename T, typename Tint>
MyVector<T> erdahl_bilinear_vector(ErdahlFunctionSpace<T> const &W,
                                   MyVector<Tint> const &e1,
                                   MyVector<Tint> const &e2) {
  int n_basis = W.basis.size();
  MyVector<T> V(n_basis);
  for (int i = 0; i < n_basis; i++) {
    V(i) = erdahl_bilinear(W.basis[i], e1, e2);
  }
  return V;
}

// The evaluation vector (f_i(x))_i of the basis of W at the point e = (1,x).
template <typename T, typename Tint>
MyVector<T> erdahl_evaluation_vector(ErdahlFunctionSpace<T> const &W,
                                     MyVector<Tint> const &e) {
  return erdahl_bilinear_vector(W, e, e);
}

// clang-format off
#endif  // SRC_ERDAHL_ERDAHL_SPACE_H_
// clang-format on
