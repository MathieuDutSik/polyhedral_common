// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_CLASSICLLL_H_
#define SRC_LATT_CLASSICLLL_H_

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
    std::cerr << "The GramMatRed is not the reduced expression\n";
    throw TerminalException{1};
  }
}

//
// Adapted from LLLReducedBasis   in zlattice.gi GAP code
//
template <typename Tmat, typename Tint>
LLLreduction<Tmat, Tint> LLLreducedBasis(MyMatrix<Tmat> const &GramMat) {
  using Tfield = typename overlying_field<Tmat>::field_type;
  MyMatrix<Tmat> gram = GramMat;
  int nbRow = gram.rows();
  int nbCol = gram.cols();
  //  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  if (nbRow != nbCol) {
    throw TerminalException{1};
  }
  int n = nbRow;
  int k = 1;
  int kmax = 0;
  MyMatrix<Tfield> mue = ZeroMatrix<Tfield>(n, n);
  int r = -1;
  MyVector<Tfield> ak(n);
  MyMatrix<Tint> H = IdentityMat<Tint>(n);
  auto RED = [&](int const &l) -> void {
    //    std::cerr << "k=" << k << " l=" << l << " mue(k,l)=" << mue(k,l) <<
    //    "\n";
    if (1 < mue(k, l) * 2 || mue(k, l) * 2 < -1) {
      Tint q = UniversalNearestScalarInteger<Tint, Tfield>(mue(k, l));
      Tmat q_T = UniversalScalarConversion<Tmat, Tint>(q);
      //      std::cerr << "RED, before oper q=" << q << "\n";
      //      WriteMatrix(std::cerr, gram);
      gram(k, k) -= q_T * gram(k, l);
      //      std::cerr << "RED step 1\n";
      for (int i = r + 1; i <= l; i++)
        gram(k, i) -= q_T * gram(l, i);
      //      std::cerr << "RED step 2\n";
      for (int i = l + 1; i <= k; i++)
        gram(k, i) -= q_T * gram(i, l);
      //      std::cerr << "RED step 3 n=" << n << "\n";
      for (int i = k + 1; i < n; i++)
        gram(i, k) -= q_T * gram(i, l);
      //      std::cerr << "After gram Oper\n";
      //      WriteMatrix(std::cerr, gram);
      //      std::cerr << "RED step 4\n";
      mue(k, l) = mue(k, l) - q_T;
      //      std::cerr << "RED step 5\n";
      for (int i = r + 1; i <= l - 1; i++)
        mue(k, i) -= q_T * mue(l, i);
      //      std::cerr << "RED step 6\n";
      H.row(k) -= q * H.row(l);
      //      std::cerr << "RED step 7\n";
    }
  };
  Tfield y = Tfield(99) / Tfield(100);
  //  std::cerr << " y=" << y << "\n";
  int i = 0;
  while (true) {
    Tmat eVal = gram(i, i);
    //    std::cerr << "i=" << i << " eVal=" << eVal << "\n";
    if (eVal > 0)
      break;
    i++;
    if (i == n)
      break;
  }
  //  std::cerr << "After the while loop i=" << i << "\n";
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
  //  std::cerr << "After the if test r=" << r << " k=" << k << " kmax=" << kmax
  //  << "\n";
  MyVector<Tfield> B(n);
  B(0) = UniversalScalarConversion<Tfield, Tmat>(gram(0, 0));
  //  std::cerr << "Before while loop\n";
  while (k < n) {
    //    std::cerr << "While loop, step 1 k=" << k << " kmax=" << kmax << "\n";
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
    //    std::cerr << "While loop, step 2 k=" << k << "\n";
    RED(k - 1);
    //    bool test=B(k) < ( y - mue(k,k-1) * mue(k,k-1) ) * B(k-1);
    //    std::cerr << "While loop, step 3 y=" << y << " mue(k,k-1)=" <<
    //    mue(k,k-1) << " B(k)=" << B(k) << " B(k-1)=" << B(k-1) << " test=" <<
    //    test << "\n";
    while (B(k) < (y - mue(k, k - 1) * mue(k, k - 1)) * B(k - 1)) {
      //      std::cerr << "Sec While loop, step 1\n";
      H.row(k).swap(H.row(k - 1));
      //      std::cerr << "Sec While loop, step 2\n";
      for (int j = r + 1; j <= k - 2; j++)
        std::swap(gram(k, j), gram(k - 1, j));
      //      std::cerr << "Sec While loop, step 3\n";
      for (int j = k + 1; j < n; j++)
        std::swap(gram(j, k), gram(j, k - 1));
      //      std::cerr << "Sec While loop, step 4\n";
      std::swap(gram(k - 1, k - 1), gram(k, k));
      //      std::cerr << "Sec While loop, step 5\n";
      for (int j = r + 1; j <= k - 2; j++)
        std::swap(mue(k, j), mue(k - 1, j));
      //      std::cerr << "Sec While loop, step 6\n";
      Tfield mmue = mue(k, k - 1);
      Tfield BB = B(k) + mmue * mmue * B(k - 1);
      //      std::cerr << "Sec While loop, step 7 BB=" << BB << "\n";
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
          //  std::cerr << "Sec While loop, step 7.2.1\n";
          Tfield q = B(k - 1) / BB;
          //  std::cerr << "Sec While loop, step 7.2.2\n";
          mue(k, k - 1) = mmue * q;
          //  std::cerr << "Sec While loop, step 7.2.3\n";
          B(k) *= q;
          //  std::cerr << "Sec While loop, step 7.2.4\n";
          B(k - 1) = BB;
          for (int i = k + 1; i <= kmax; i++) {
            //  std::cerr << "Sec While loop, step 7.2.5.1\n";
            Tfield q = mue(i, k);
            //  std::cerr << "Sec While loop, step 7.2.5.2\n";
            mue(i, k) = mue(i, k - 1) - mmue * q;
            //  std::cerr << "Sec While loop, step 7.2.5.3\n";
            mue(i, k - 1) = q + mue(k, k - 1) * mue(i, k);
            //  std::cerr << "Sec While loop, step 7.2.5.4\n";
          }
          //  std::cerr << "Sec While loop, step 7.2.6\n";
        }
      }
      //  std::cerr << "Sec While loop, step 8\n";
      if (k > 1)
        k--;
      RED(k - 1);
      //  std::cerr << "Sec While loop, step 9\n";
    }
    //  std::cerr << "While loop, step 4\n";
    if (B(r + 1) == 0)
      r++;
    //  std::cerr << "While loop, step 5 k=" << k << " r=" << r << "\n";
    for (int l = k - 2; l >= r + 1; l--) {
      //  std::cerr << "  While loop, step 5.1 l=" << l << "\n";
      RED(l);
    }
    //  std::cerr << "While loop, step 6\n";
    k++;
  }
  for (int i = 1; i < n; i++)
    for (int j = 0; j < i; j++)
      gram(j, i) = gram(i, j);
  LLLreduction<Tmat, Tint> res = {std::move(gram), std::move(H)};
#ifdef DEBUG
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
LLLreduction<Tmat, Tint> LLLreducedBasisDual(MyMatrix<Tmat> const &GramMat) {
  MyMatrix<Tmat> Ginv = Inverse(GramMat);
  LLLreduction<Tmat, Tint> LLLrec = LLLreducedBasis<Tmat, Tint>(Ginv);
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
                                           std::string const &method) {
  if (method == "direct")
    return LLLreducedBasis<Tmat, Tint>(GramMat);
  if (method == "dual")
    return LLLreducedBasisDual<Tmat, Tint>(GramMat);
  std::cerr << "No matching method\n";
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
LLLbasis<T, Tint> LLLbasisReduction(MyMatrix<T> const &Latt) {
  //  std::cerr << "LLLbasisReduction, step 1\n";
  //  std::cerr << "|Latt|=" << Latt.rows() << " / " << Latt.cols() << "\n";
  //  std::cerr << "Latt=\n";
  //  WriteMatrix(std::cerr, Latt);
  MyMatrix<T> GramMat = Latt * Latt.transpose();
  //  std::cerr << "LLLbasisReduction, step 2\n";
  LLLreduction<T, Tint> pair = LLLreducedBasis<T, Tint>(GramMat);
  //  std::cerr << "LLLbasisReduction, step 3\n";
  MyMatrix<T> LattRed = UniversalMatrixConversion<T, Tint>(pair.Pmat) * Latt;
  //  std::cerr << "LLLbasisReduction, step 4\n";
  return {LattRed, pair.Pmat};
}

template <typename T>
std::pair<MyMatrix<T>, MyMatrix<T>>
ReduceVectorFamily(MyMatrix<T> const &M, std::string const &method) {
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
  std::cerr << "ReduceVectorFamily, step 1\n";
  LLLreduction<T, Tint> res = LLLreducedGeneral<T, Tint>(TheGram, method);
  std::cerr << "ReduceVectorFamily, step 2\n";
  MyMatrix<Tint> Pmat = TransposedMat(res.Pmat);
  std::cerr << "ReduceVectorFamily, step 3\n";
  MyMatrix<T> Pmat_T = UniversalMatrixConversion<T, Tint>(Pmat);
  std::cerr << "ReduceVectorFamily, step 4\n";
  std::cerr << "|M|=" << M.rows() << "/" << M.cols()
            << " |Pmat_T|=" << Pmat_T.rows() << "/" << Pmat_T.cols() << "\n";
  MyMatrix<T> Mred = M * Pmat_T;
  std::cerr << "ReduceVectorFamily, step 5\n";
  if (GetGram(Mred) != res.GramMatRed) {
    std::cerr << "Matrix error somewhere\n";
    throw TerminalException{1};
  }
  return {std::move(Mred), std::move(Pmat_T)};
}

// clang-format off
#endif  // SRC_LATT_CLASSICLLL_H_
// clang-format on
