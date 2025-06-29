// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_CLASSICLLL_H_
#define SRC_LATT_CLASSICLLL_H_

#ifdef DEBUG
#define DEBUG_CLASSIC_LLL
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_CLASSIC_LLL
#endif


template <typename T, typename Tint> struct LLLreduction {
  MyMatrix<T> GramMatRed;
  MyMatrix<Tint> Pmat;
};

template <typename T, typename Tint>
void CheckLLLreduction(LLLreduction<T, Tint> const &res,
                       MyMatrix<T> const &GramMat) {
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(res.Pmat);
  MyMatrix<T> prod = Pmat_T * GramMat * Pmat_T.transpose();
  if (prod != res.GramMatRed) {
    std::cerr << "LLL: The GramMatRed is not the reduced expression\n";
    throw TerminalException{1};
  }
}

//
// Adapted from LLLReducedBasis   in zlattice.gi GAP code
//
template <typename Tmat, typename Tint>
LLLreduction<Tmat, Tint> LLLreducedBasis(MyMatrix<Tmat> const &GramMat, [[maybe_unused]] std::ostream& os) {
  using Tfield = typename overlying_field<Tmat>::field_type;
  MyMatrix<Tmat> gram = GramMat;
  int nbRow = gram.rows();
#ifdef DEBUG_CLASSIC_LLL
  os << "LLL: nbRow=" << nbRow << " nbCol=" << gram.cols() << "\n";
  os << "LLL: GramMat=\n";
  WriteMatrix(os, GramMat);
#endif
#ifdef SANITY_CHECK_CLASSIC_LLL
  if (nbRow != gram.cols()) {
    std::cerr << "LLL: The matrix should be square\n";
    throw TerminalException{1};
  }
#endif
  int n = nbRow;
#ifdef SANITY_CHECK_CLASSIC_LLL
  int rnk = RankMat(GramMat);
  if (rnk != n) {
    std::cerr << "LLL: rnk=" << rnk << " n=" << n << "\n";
    std::cerr << "LLL: The matrix is not of full rank\n";
    throw TerminalException{1};
  }
#endif
  if (nbRow == 1 || nbRow == 0) {
    MyMatrix<Tmat> GramMatRet = GramMat;
    MyMatrix<Tint> H = IdentityMat<Tint>(nbRow);
    LLLreduction<Tmat, Tint> res{GramMatRet, H};
    return res;
  }
#ifdef SANITY_CHECK_CLASSIC_LLL
  if (!IsPositiveDefinite(GramMat)) {
    std::cerr << "LLL: For the LLL reduction, the matrix needs to be positive definite\n";
    throw TerminalException{1};
  }
#endif
  int k = 1;
  int kmax = 0;
  MyMatrix<Tfield> mue = ZeroMatrix<Tfield>(n, n);
  int r = -1;
  MyVector<Tfield> ak(n);
  MyMatrix<Tint> H = IdentityMat<Tint>(n);
  auto RED = [&](int const &l) -> void {
#ifdef DEBUG_CLASSIC_LLL
    os << "LLL: k=" << k << " l=" << l << " mue(k,l)=" << mue(k,l) << "\n";
#endif
    if (1 < mue(k, l) * 2 || mue(k, l) * 2 < -1) {
      Tint q = UniversalNearestScalarInteger<Tint, Tfield>(mue(k, l));
      Tmat q_T = UniversalScalarConversion<Tmat, Tint>(q);
#ifdef DEBUG_CLASSIC_LLL
      os << "LLL: RED, before oper q=" << q << "\n";
      WriteMatrix(os, gram);
#endif
      gram(k, k) -= q_T * gram(k, l);
      for (int i = r + 1; i <= l; i++)
        gram(k, i) -= q_T * gram(l, i);
      for (int i = l + 1; i <= k; i++)
        gram(k, i) -= q_T * gram(i, l);
      for (int i = k + 1; i < n; i++)
        gram(i, k) -= q_T * gram(i, l);
#ifdef DEBUG_CLASSIC_LLL
      os << "LLL: After gram Oper\n";
      WriteMatrix(os, gram);
#endif
      mue(k, l) = mue(k, l) - q_T;
      for (int i = r + 1; i <= l - 1; i++)
        mue(k, i) -= q_T * mue(l, i);
      H.row(k) -= q * H.row(l);
    }
  };
  Tfield y = Tfield(99) / Tfield(100);
#ifdef DEBUG_CLASSIC_LLL
  os << "LLL: y=" << y << "\n";
#endif
  int i = 0;
  while (true) {
    Tmat eVal = gram(i, i);
    if (eVal > 0)
      break;
    i++;
    if (i == n)
      break;
  }
#ifdef DEBUG_CLASSIC_LLL
  os << "LLL: After the while loop i=" << i << "\n";
#endif
  if (i >= n) {
    r = n - 1;
    k = n;
  } else {
    if (i > 0) {
      for (int j = i + 1; j < n; j++) {
        gram(j, 0) = gram(j, i);
        gram(j, i) = 0;
      }
      gram(0, 0) = gram(i, i);
      gram(i, i) = 0;
      H.row(i).swap(H.row(i));
    }
  }
#ifdef DEBUG_CLASSIC_LLL
  os << "LLL: After the if test r=" << r << " k=" << k << " kmax=" << kmax << "\n";
#endif
  MyVector<Tfield> B(n);
  B(0) = UniversalScalarConversion<Tfield, Tmat>(gram(0, 0));
  while (k < n) {
#ifdef DEBUG_CLASSIC_LLL
    os << "LLL: While loop, step 1 k=" << k << " kmax=" << kmax << "\n";
#endif
    if (k > kmax) {
      kmax = k;
      B(k) = UniversalScalarConversion<Tfield, Tmat>(gram(k, k));
      for (int u = 0; u < n; u++)
        mue(k, u) = 0;
      for (int j = r + 1; j <= k - 1; j++) {
        ak(j) = UniversalScalarConversion<Tfield, Tmat>(gram(k, j));
        for (int i = r + 1; i <= j - 1; i++)
          ak(j) -= mue(j, i) * ak(i);
        mue(k, j) = ak(j) / B(j);
        B(k) -= mue(k, j) * ak(j);
      }
    }
    RED(k - 1);
#ifdef DEBUG_CLASSIC_LLL
    bool test=B(k) < ( y - mue(k,k-1) * mue(k,k-1) ) * B(k-1);
    os << "LLL: While loop, step 3 y=" << y << " mue(k,k-1)=" <<
      mue(k,k-1) << " B(k)=" << B(k) << " B(k-1)=" << B(k-1) << " test=" <<
      test << "\n";
#endif
    while (B(k) < (y - mue(k, k - 1) * mue(k, k - 1)) * B(k - 1)) {
      H.row(k).swap(H.row(k - 1));
      for (int j = r + 1; j <= k - 2; j++)
        std::swap(gram(k, j), gram(k - 1, j));
      for (int j = k + 1; j < n; j++)
        std::swap(gram(j, k), gram(j, k - 1));
      std::swap(gram(k - 1, k - 1), gram(k, k));
      for (int j = r + 1; j <= k - 2; j++)
        std::swap(mue(k, j), mue(k - 1, j));
      Tfield mmue = mue(k, k - 1);
      Tfield BB = B(k) + mmue * mmue * B(k - 1);
      if (BB == 0) {
        B(k) = B(k - 1);
        B(k - 1) = 0;
        for (int i = k + 1; k <= kmax; k++)
          std::swap(mue(i, k), mue(i, k - 1));
      } else {
        if (B(k) == 0 && mmue != 0) {
          B(k - 1) = BB;
          mue(k, k - 1) = 1 / mmue;
          for (int i = k + 1; k <= kmax; k++)
            mue(i, k - 1) /= mmue;
        } else {
          Tfield q = B(k - 1) / BB;
          mue(k, k - 1) = mmue * q;
          B(k) *= q;
          B(k - 1) = BB;
          for (int i = k + 1; i <= kmax; i++) {
            Tfield q = mue(i, k);
            mue(i, k) = mue(i, k - 1) - mmue * q;
            mue(i, k - 1) = q + mue(k, k - 1) * mue(i, k);
          }
        }
      }
      if (k > 1)
        k--;
      RED(k - 1);
    }
    if (B(r + 1) == 0)
      r++;
#ifdef DEBUG_CLASSIC_LLL
    os << "LLL: While loop, step 5 k=" << k << " r=" << r << "\n";
#endif
    for (int l = k - 2; l >= r + 1; l--) {
      RED(l);
    }
    k++;
  }
  for (int i = 1; i < n; i++)
    for (int j = 0; j < i; j++)
      gram(j, i) = gram(i, j);
  LLLreduction<Tmat, Tint> res = {std::move(gram), std::move(H)};
#ifdef DEBUG_CLASSIC_LLL
  CheckLLLreduction(res, GramMat);
#endif
  return res;
}




/*
  We reduced the inverse of GramMat using LLL so
  P GramMat^{-1} P^T = Gred
  Inverting we get
  (P^T)^{-1} GramMat P^-1 = Gred^{-1}
  Thus the result is the pair (Gred^{-1} , Q) with Q = (P^T)^{-1}
  This means that in any code, we can substritute LLLreducedBasis with
  LLLreducedBasisDual and it should work just as well.
 */
template <typename Tmat, typename Tint>
LLLreduction<Tmat, Tint> LLLreducedBasisDual(MyMatrix<Tmat> const &GramMat, std::ostream& os) {
  MyMatrix<Tmat> Ginv = Inverse(GramMat);
  LLLreduction<Tmat, Tint> LLLrec = LLLreducedBasis<Tmat, Tint>(Ginv, os);
  MyMatrix<Tmat> const &Gred = LLLrec.GramMatRed;
  MyMatrix<Tint> const &P = LLLrec.Pmat;
  MyMatrix<Tmat> GredInv = Inverse(Gred);
  MyMatrix<Tint> Q = TransposedMat(Inverse(P));
  LLLreduction<Tmat, Tint> res = {std::move(GredInv), std::move(Q)};
#ifdef DEBUG
  CheckLLLreduction(res, GramMat);
#endif
  return res;
}

template <typename Tmat, typename Tint>
LLLreduction<Tmat, Tint> LLLreducedGeneral(MyMatrix<Tmat> const &GramMat,
                                           std::string const &method, std::ostream& os) {
  if (method == "direct")
    return LLLreducedBasis<Tmat, Tint>(GramMat, os);
  if (method == "dual")
    return LLLreducedBasisDual<Tmat, Tint>(GramMat, os);
  std::cerr << "LLL: No matching method\n";
  throw TerminalException{1};
}

// This is for debugging purposes
template <typename Tmat, typename Tint>
LLLreduction<Tmat, Tint> LLLnoreduction(MyMatrix<Tmat> const &GramMat) {
  return {GramMat, IdentityMat<Tint>(GramMat.rows())};
}

/*
  For a family of vectors of an N-dimensional space (possibly in a higher
  dimensional one), returns a smaller basis.
  LattRed is the reduced matrix and Pmat is the reducing matrix
 */
template <typename T, typename Tint> struct LLLbasis {
  MyMatrix<T> LattRed;
  MyMatrix<Tint> Pmat;
};

template <typename T, typename Tint>
LLLbasis<T, Tint> LLLbasisReduction(MyMatrix<T> const &Latt, std::ostream& os) {
  MyMatrix<T> GramMat = Latt * Latt.transpose();
  LLLreduction<T, Tint> pair = LLLreducedBasis<T, Tint>(GramMat, os);
  MyMatrix<T> LattRed = UniversalMatrixConversion<T, Tint>(pair.Pmat) * Latt;
  return {LattRed, pair.Pmat};
}

template <typename T, typename Tint>
LLLbasis<T, Tint> LLLbasisReductionGeneral(MyMatrix<T> const &Latt, std::string const& method, std::ostream& os) {
  MyMatrix<T> GramMat = Latt * Latt.transpose();
  LLLreduction<T, Tint> pair = LLLreducedGeneral<T, Tint>(GramMat, method, os);
  MyMatrix<T> LattRed = UniversalMatrixConversion<T, Tint>(pair.Pmat) * Latt;
  return {LattRed, pair.Pmat};
}

template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>>
ReduceVectorFamily(MyMatrix<T> const &M, std::string const &method, std::ostream& os) {
  using Tint = typename underlying_ring<T>::ring_type;
  int nbRow = M.rows();
  int nbCol = M.cols();
  auto GetGram = [&](MyMatrix<T> const &VectFamily) {
    MyMatrix<T> TheGram = ZeroMatrix<T>(nbCol, nbCol);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      for (int iCol = 0; iCol < nbCol; iCol++) {
        for (int jCol = 0; jCol < nbCol; jCol++) {
          TheGram(iCol, jCol) +=
              VectFamily(iRow, iCol) * VectFamily(iRow, jCol);
        }
      }
    }
    return TheGram;
  };
  MyMatrix<T> TheGram = GetGram(M);
  LLLreduction<T, Tint> res = LLLreducedGeneral<T, Tint>(TheGram, method, os);
  MyMatrix<Tint> Pmat = TransposedMat(res.Pmat);
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  MyMatrix<T> Mred = M * Pmat_T;
#ifdef SANITY_CHECK_CLASSIC_LLL
  if (GetGram(Mred) != res.GramMatRed) {
    std::cerr << "LLL: Matrix error somewhere\n";
    throw TerminalException{1};
  }
#endif
  return {std::move(Mred), std::move(Pmat_T)};
}


/*
  The code for finding a nice basis of a lattice.
  ----
  We use LLL, but there are other alternatives to consider:
  * Exhaustive search as done for indefinite matrices.
 */
template<typename T>
MyMatrix<T> SublatticeBasisReductionKernel(MyMatrix<T> const& Latt, std::ostream& os) {
  using Tint = typename underlying_ring<T>::ring_type;
  int n_row = Latt.rows();
  int n_col = Latt.cols();
  auto f_norm=[&](MyMatrix<T> const& H) -> T {
    T norm(0);
    for (int i_row=0; i_row<n_row; i_row++) {
      for (int i_col=0; i_col<n_col; i_col++) {
        norm += T_abs(H(i_row, i_col));
      }
    }
    return norm;
  };
  MyMatrix<T> Latt_work = Latt;
  T norm_work = f_norm(Latt);
  std::vector<std::string> l_method{"direct", "dual"};
  while(true) {
    int n_success = 0;
    for (auto & method: l_method) {
      LLLbasis<T, Tint> rec = LLLbasisReductionGeneral<T,Tint>(Latt_work, method, os);
      MyMatrix<T> const& Latt_cand = rec.LattRed;
      T norm_cand = f_norm(Latt_cand);
      if (norm_cand < norm_work) {
        n_success += 1;
        Latt_work = Latt_cand;
        norm_work = norm_cand;
      }
    }
    if (n_success == 0) {
      return Latt_work;
    }
  }
}

template<typename T>
inline typename std::enable_if<is_ring_field<T>::value, MyMatrix<T>>::type
SublatticeBasisReduction(MyMatrix<T> const& Latt, std::ostream& os) {
  return SublatticeBasisReductionKernel(Latt, os);
}

template<typename T>
inline typename std::enable_if<!is_ring_field<T>::value, MyMatrix<T>>::type
SublatticeBasisReduction(MyMatrix<T> const& Latt, std::ostream& os) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> Latt_F = UniversalMatrixConversion<Tfield,T>(Latt);
  MyMatrix<Tfield> LattRed_F = SublatticeBasisReductionKernel(Latt_F, os);
  MyMatrix<T> LattRed = UniversalMatrixConversion<T,Tfield>(LattRed_F);
  return LattRed;
}

// clang-format off
#endif  // SRC_LATT_CLASSICLLL_H_
// clang-format on
