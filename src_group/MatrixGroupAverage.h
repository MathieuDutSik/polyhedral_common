// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPAVERAGE_H_
#define SRC_GROUP_MATRIXGROUPAVERAGE_H_

// clang-format off
#include "MAT_Matrix.h"
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_MATRIX_GROUP_AVERAGE
#endif

/*
  Builds the smallest vector space containing the space TheBasis and invariant
  by the group generated by LGen. Code adapted from DirectSpannEquivariantSpace
  in the GAP version.
  ---
  The MutableSubspaceBelongingRepetitive was designed for that goal and seems
  to work out fine and seem to be optimal.
 */
template <typename T>
MyMatrix<T> DirectSpannEquivariantSpace(MyMatrix<T> const &TheBasis,
                                        std::vector<MyMatrix<T>> const &LGen) {
  int dim = TheBasis.cols();
  MutableSubspaceBelongingRepetitive<T> msbr(dim);
  msbr.InsertMatrix(TheBasis);
  size_t pos_start = 0;
  size_t pos_end = msbr.get_n_vect();
  MyVector<T> V(dim);
#ifdef DEBUG_MATRIX_GROUP_AVERAGE
  int n_iter = 0;
#endif
  while (true) {
    if (pos_start == pos_end) {
      break;
    }
    for (size_t pos = pos_start; pos < pos_end; pos++) {
      msbr.SetVectout(pos, V);
      for (auto &eGen : LGen) {
        MyVector<T> Vimg = eGen.transpose() * V;
        (void)msbr.InsertVector(Vimg);
      }
    }
    pos_start = pos_end;
    pos_end = msbr.get_n_vect();
#ifdef DEBUG_MATRIX_GROUP_AVERAGE
    n_iter += 1;
    std::cerr << "GA n_iter=" << n_iter << " pos_start=" << pos_start
              << " pos_end=" << pos_end << "\n";
#endif
  }
  return msbr.GetBasis();
}

/*
  Compute the average of a vector under the group generated by the generator.
  The code computes the average under the group generated without building the
  full orbit.
  Code adapted from OrbitBarycenter from the GAP version.
 */
template <typename T>
MyVector<T> OrbitBarycenter(MyVector<T> const &a,
                            std::vector<MyMatrix<T>> const &LGen) {
#ifdef DEBUG_MATRIX_GROUP_AVERAGE
  std::cerr << "GA: OrbitBarycenter beginning\n";
#endif
  auto is_invariant = [&](MyVector<T> const &u) -> bool {
    for (auto &eGen : LGen) {
      MyVector<T> u_img = eGen.transpose() * u;
      if (u_img != u) {
        return false;
      }
    }
    return true;
  };
  if (is_invariant(a)) {
    return a;
  }
  // Not invariant, need to build a linear system
  int n_gen = LGen.size();
  int dim = a.size();
  MyMatrix<T> ListSpann(n_gen, dim);
  for (int i_gen = 0; i_gen < n_gen; i_gen++) {
    MyMatrix<T> const &eGen = LGen[i_gen];
    MyVector<T> a_img = eGen.transpose() * a;
    for (int u = 0; u < dim; u++) {
      ListSpann(i_gen, u) = a_img(u) - a(u);
    }
  }
#ifdef DEBUG_MATRIX_GROUP_AVERAGE
  std::cerr << "GA: rnk(ListSpann)=" << RankMat(ListSpann) << " n_gen=" << n_gen
            << " dim=" << dim << "\n";
#endif
  MyMatrix<T> TheBasis2 = DirectSpannEquivariantSpace(ListSpann, LGen);
#ifdef DEBUG_MATRIX_GROUP_AVERAGE
  std::cerr << "GA: |TheBasis2|=" << TheBasis2.rows() << " / "
            << TheBasis2.cols() << "\n";
#endif
  int dim_space = TheBasis2.rows();
  MyMatrix<T> TheBigMat(dim_space, n_gen * dim);
  MyVector<T> Vbig(n_gen * dim);
  for (int i_gen = 0; i_gen < n_gen; i_gen++) {
    MyMatrix<T> const &eGen = LGen[i_gen];
    MyMatrix<T> diff_M = TheBasis2 * eGen - TheBasis2;
    for (int iRow = 0; iRow < dim_space; iRow++) {
      for (int iCol = 0; iCol < dim; iCol++) {
        TheBigMat(iRow, iCol + dim * i_gen) = diff_M(iRow, iCol);
      }
    }
    MyVector<T> diff_V = a - eGen.transpose() * a;
    for (int iCol = 0; iCol < dim; iCol++) {
      Vbig(iCol + dim * i_gen) = diff_V(iCol);
    }
  }
  std::optional<MyVector<T>> opt = SolutionMat(TheBigMat, Vbig);
  MyVector<T> Alpha = unfold_opt(opt, "Failed to solve the linear system");
  MyVector<T> TheSol = a + TheBasis2.transpose() * Alpha;
#ifdef DEBUG_MATRIX_GROUP_AVERAGE
  if (!is_invariant(TheSol)) {
    std::cerr << "The vector TheSol is not invariant\n";
    throw TerminalException{1};
  }
#endif
  return TheSol;
}

/*
  Now doing it over the action of the group over symmetric matrices by
  M -> g M g^T
 */
template <typename T>
MyMatrix<T>
OrbitBarycenterSymmetricMatrix(MyMatrix<T> const &M,
                               std::vector<MyMatrix<T>> const &LGen) {
  int n = M.rows();
  int sym_dim = (n * (n + 1)) / 2;
  MyVector<T> a(sym_dim);
  MyMatrix<int> C(sym_dim, 2);
  MyMatrix<int> Crev(n, n);
  int off_dim = (n * (n - 1)) / 2;
  std::vector<int> LOff(off_dim);
  int pos_sym = 0;
  int pos_off = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      a(pos_sym) = M(i, j);
      C(pos_sym, 0) = i;
      C(pos_sym, 1) = j;
      Crev(i, j) = pos_sym;
      Crev(j, i) = pos_sym;
      if (i != j) {
        LOff[pos_off] = pos_sym;
        pos_off++;
      }
      pos_sym++;
    }
  }
  MyMatrix<int> LHalf(off_dim, 2);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      LHalf(pos, 0) = i;
      LHalf(pos, 1) = j;
      pos++;
    }
  }
  std::vector<MyMatrix<T>> LGenExt;
  for (auto &eGen : LGen) {
    MyMatrix<T> eGenExt = ZeroMatrix<T>(sym_dim, sym_dim);
    for (int i_sym = 0; i_sym < sym_dim; i_sym++) {
      int i = C(i_sym, 0);
      int j = C(i_sym, 1);
      // A matrix entry with 1 in position (i,j) and (j,i)
      //
      // We need to Compute the products
      // P E_{i,i} P^T and
      // P (E_{i,j} + E_{j,i}) P^T
      //
      // Let us compute the initial product
      // P E_{i,j} P^T
      // First the product P E_{i,j}
      // It is a n x n matrix with the i-th column of P in the column j;
      // The product P E_{i,j} P^T is
      // [p_{k,i} p_{l,j}]_{1 \leq k,l \leq n}
      // k is row and l is the column.
      for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
          int pos = Crev(k, l);
          eGenExt(i_sym, pos) += eGen(k, i) * eGen(l, j);
          if (i != j) {
            eGenExt(i_sym, pos) += eGen(k, j) * eGen(l, i);
          }
        }
      }
    }
    for (int i_off = 0; i_off < off_dim; i_off++) {
      int pos_sym = LOff[i_off];
      for (int i_sym = 0; i_sym < sym_dim; i_sym++) {
        eGenExt(i_sym, pos_sym) /= 2;
      }
    }
    LGenExt.push_back(eGenExt);
  }
  MyVector<T> a_bary = OrbitBarycenter(a, LGenExt);
  MyMatrix<T> M_bary(n, n);
  int pos_sym_b = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      M_bary(i, j) = a_bary(pos_sym_b);
      M_bary(j, i) = a_bary(pos_sym_b);
      pos_sym_b++;
    }
  }
#ifdef DEBUG_MATRIX_GROUP_AVERAGE
  for (auto &eGen : LGen) {
    MyMatrix<T> Mimg = eGen * M_bary * eGen.transpose();
    if (Mimg != M_bary) {
      std::cerr << "Failed to find an invariant matrix. Clear bug\n";
      throw TerminalException{1};
    }
  }
#endif
  return M_bary;
}

// clang-format off
#endif  // SRC_GROUP_MATRIXGROUPAVERAGE_H_
// clang-format on
