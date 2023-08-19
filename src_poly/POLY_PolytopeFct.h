// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_POLYTOPEFCT_H_
#define SRC_POLY_POLY_POLYTOPEFCT_H_

// clang-format off
#include "Boost_bitset.h"
#include "COMB_Stor.h"
#include "MAT_Matrix.h"
#include "rational.h"
#ifndef DISABLE_FP_CLASS
#include "Fp.h"
#endif
#include "NumberTheoryGeneric.h"
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
# define DEBUG_FLIP
#endif

struct GLPKoption {
  bool UseDouble;
  bool UseExact;
  bool UseXcheck;
};

template <typename T> struct LpSolutionSimple {
  bool PrimalDefined;
  T OptimalValue;
  int nbRow;
  int nbCol;
  MyVector<T> DirectSolution;
  MyVector<T> DirectSolutionExt;
  // Value 0 for not assigned.
  // Value 1 for "B"
  // Value 2 for "NF"
  // Value 3 for "NL"
  MyVector<int> RowStatus;
  MyVector<int> ColumnStatus;
};

template <typename T>
MyVector<T> SumMatrixLineSubset(MyMatrix<T> const &eMat, Face const &eList) {
  int nbCol = eMat.cols();
  MyVector<T> eVec = ZeroVector<T>(nbCol);
  int eSize = eList.count();
  //
  boost::dynamic_bitset<>::size_type aRow = eList.find_first();
  for (int i = 0; i < eSize; i++) {
    for (int iCol = 0; iCol < nbCol; iCol++)
      eVec(iCol) += eMat(aRow, iCol);
    aRow = eList.find_next(aRow);
  }
  return eVec;
}

template <typename T>
MyMatrix<T> SelectRow(MyMatrix<T> const &TheMat, Face const &eList) {
  int nbRowRed = eList.count();
  int nbCol = TheMat.cols();
  MyMatrix<T> TheProv(nbRowRed, nbCol);
  boost::dynamic_bitset<>::size_type jRow = eList.find_first();
  for (int iRow = 0; iRow < nbRowRed; iRow++) {
    TheProv.row(iRow) = TheMat.row(jRow);
    jRow = eList.find_next(jRow);
  }
  return TheProv;
}

template <typename T> MyMatrix<T> CyclicPolytope(int n, int k) {
  int i, j, b;
  MyMatrix<T> TheEXT(n, k + 1);
  for (i = 1; i <= n; i++) {
    T a = 1;
    b = i + 1;
    for (j = 0; j <= k; j++) {
      TheEXT(i - 1, j) = a;
      a = a * b;
    }
  }
  return TheEXT;
}

// Compute the maximal determinant
// We have det(A)^2 <= Pi_{i=1}^n (sum_j x_{ij}^2)
template<typename T>
T sqr_estimate_maximal_determinant(MyMatrix<T> const& M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  std::vector<T> ListSqr;
  for (int iRow=0; iRow<nbRow; iRow++) {
    T sum = 0;
    for (int iCol=0; iCol<nbCol; iCol++)
      sum += M(iRow,iCol) * M(iRow,iCol);
    ListSqr.push_back(sum);
  }
  /*
  std::cerr << "ListSqr=";
  for (auto & eVal : ListSqr)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  */
  std::sort(ListSqr.begin(), ListSqr.end());
  T eProd = 1;
  for (int iRow=nbRow-nbCol; iRow<nbRow; iRow++) {
    eProd *= ListSqr[iRow];
  }
  return eProd;
}

template<typename T>
T sqr_estimate_facet_coefficients(MyMatrix<T> const& M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  T max_coeff = 0;
  MyVector<T> SqrCoeffFacet(nbCol);
  for (int jCol=0; jCol<nbCol; jCol++) {
    MyMatrix<T> Mret(nbRow,nbCol-1);
    for (int iRow=0; iRow<nbRow; iRow++) {
      int pos = 0;
      for (int iCol=0; iCol<nbCol; iCol++) {
        if (iCol != jCol) {
          Mret(iRow,pos) = M(iRow,iCol);
          pos++;
        }
      }
    }
    T est = sqr_estimate_maximal_determinant(Mret);
    if (est > max_coeff) {
      max_coeff = est;
    }
  }
  return max_coeff;
}

template <typename T>
MyVector<T> FindFacetInequalityCheck(MyMatrix<T> const &EXT, Face const &eList) {
  int nb = eList.count();
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  boost::dynamic_bitset<>::size_type aRow = eList.find_first();
  auto f = [&](MyMatrix<T> &M, size_t eRank,
               [[maybe_unused]] size_t iRow) -> void {
    M.row(eRank) = EXT.row(aRow);
    aRow = eList.find_next(aRow);
  };
  MyMatrix<T> NSP = NullspaceTrMat_Kernel<T, decltype(f)>(nb, nbCol, f);
  if (NSP.rows() != 1) {
    std::cerr << "Error in rank in Facetness\n";
    std::cerr << "|NSP|=" << NSP.rows() << "\n";
    throw TerminalException{1};
  }
  int nbPlus = 0;
  int nbMinus = 0;
  int nbError = 0;
  MyVector<T> eVect(nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++)
    eVect(iCol) = NSP(0, iCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eScal(0);
    for (int iCol = 0; iCol < nbCol; iCol++)
      eScal += eVect(iCol) * EXT(iRow, iCol);
    if (eScal == 0) {
      if (eList[iRow] != 1) {
        std::cerr << "The vertex iRow has a zero scalar product but does not belong to eList\n";
        nbError++;
      }
    } else {
      if (eList[iRow] != 0) {
        std::cerr << "The vertex iRow has a non-zero scalar product but does belong to eList\n";
        nbError++;
      }
      if (eScal > 0)
        nbPlus++;
      if (eScal < 0)
        nbMinus++;
    }
  }
  if (nbMinus > 0 && nbPlus > 0) {
    std::cerr << "Some plus and minus signs, illegal\n";
    nbError++;
  }
  if (nbError > 0) {
    std::cerr << "nbMinus=" << nbMinus << " nbPlus=" << nbPlus << "\n";
    std::cerr << "EXT(rows/cols)=" << EXT.rows() << " / " << EXT.cols() << "\n";
    std::cerr << "Rank(EXT)=" << RankMat(EXT) << "\n";
    throw TerminalException{1};
  }
  if (nbPlus > 0)
    return eVect;
  return -eVect;
}

template <typename T>
MyVector<T> FindFacetInequality(MyMatrix<T> const &TheEXT, Face const &OneInc) {
  size_t nb = OneInc.count();
  size_t nbRow = TheEXT.rows();
  size_t nbCol = TheEXT.cols();
  boost::dynamic_bitset<>::size_type aRow = OneInc.find_first();
  auto f = [&](MyMatrix<T> &M, size_t eRank,
               [[maybe_unused]] size_t iRow) -> void {
    M.row(eRank) = TheEXT.row(aRow);
    aRow = OneInc.find_next(aRow);
  };
  MyMatrix<T> NSP = NullspaceTrMatTarget_Kernel<T, decltype(f)>(nb, nbCol, 1, f);
  MyVector<T> eVect(nbCol);
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    eVect(iCol) = NSP(0, iCol);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    if (OneInc[iRow])
      continue;
    T eScal(0);
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      eScal += eVect(iCol) * TheEXT(iRow, iCol);
    if (eScal > 0)
      return eVect;
    if (eScal < 0)
      return -eVect;
  }
  std::cerr << "FindFacetInequality: Failed to find the defining inequality\n";
  throw TerminalException{1};
}

template <typename T>
int GetFacetRank(MyMatrix<T> const &TheEXT, Face const &OneInc) {
  size_t nb = OneInc.count();
  size_t nbCol = TheEXT.cols();
  boost::dynamic_bitset<>::size_type aRow = OneInc.find_first();
  auto f = [&](MyMatrix<T> &M, size_t eRank,
               [[maybe_unused]] size_t iRow) -> void {
    M.row(eRank) = TheEXT.row(aRow);
    aRow = OneInc.find_next(aRow);
  };
  MyMatrix<T> NSP = NullspaceTrMat_Kernel<T, decltype(f)>(nb, nbCol, f);
  return NSP.rows();
}

std::pair<std::vector<int>,std::vector<int>> Dynamic_bitset_to_vectorints(Face const &eList) {
  int len = eList.size();
  std::vector<int> V0;
  std::vector<int> V1;
  for (int i=0; i<len; i++) {
    if (eList[i]) {
      V1.push_back(i);
    } else {
      V0.push_back(i);
    }
  }
  return {std::move(V0), std::move(V1)};
}


// This is the FlippingFramework for a given facet F of a polytope.
//
// After the constructor is built, then we provide a function for
// given a facet G of F, find the facet F' such that G = F \cap F'.
//
// The computation requires a number of finding of the kernel of
// matrices which are of rank 1. Thus we use the NullspaceTrMatTarget_Kernel
// function.
//
// There are two fuctions:
// 1) FlipFace for flipping only a face.
// 2) FlipFaceIneq for flipping the face and the inequality if it is already
// known. Normally, the inequality is produce by the various techniques
// used. But the problem is to store them.
template <typename T> struct FlippingFramework {
private:
  MyMatrix<T> EXT_red;
  int nbRow;
  int nbCol;
  Face OneInc;
  std::pair<std::vector<int>,std::vector<int>> PairIncs;
  int e_incd0;
  int e_incd1;
  std::vector<T> ListInvScal;

public:
  MyMatrix<T> EXT_face;
  FlippingFramework(MyMatrix<T> const &EXT, Face const &_OneInc)
    : OneInc(_OneInc), e_incd0(OneInc.size() - OneInc.count()), e_incd1(OneInc.count()), ListInvScal(e_incd0) {
    PairIncs = Dynamic_bitset_to_vectorints(OneInc);
    MyVector<T> FacetIneq = FindFacetInequality(EXT, OneInc);
    //
    // Idx dropping for the projection
    //
    int idx_drop = 0;
    while (true) {
      if (FacetIneq(idx_drop) != 0)
        break;
      idx_drop++;
    }
    nbRow = EXT.rows();
    nbCol = EXT.cols();
    EXT_red = MyMatrix<T>(nbRow, nbCol - 1);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      int pos = 0;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol != idx_drop) {
          EXT_red(iRow, pos) = EXT(iRow, iCol);
          pos++;
        }
      }
    }
    //
    // Inverse scalar products
    //
    for (int pos_row=0; pos_row<e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      T eSum(0);
      for (int iCol = 0; iCol < nbCol; iCol++)
        eSum += FacetIneq(iCol) * EXT(iRow, iCol);
      ListInvScal[pos_row] = -1 / eSum;
    }
    //
    // Now the EXT face that is used by other procedure
    //
    EXT_face = MyMatrix<T>(e_incd1, nbCol - 1);
    for (int i_row=0; i_row<e_incd1; i_row++) {
      int j_row = PairIncs.second[i_row];
      int pos = 0;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol != idx_drop) {
          EXT_face(i_row, pos) = EXT(j_row, iCol);
          pos++;
        }
      }
    }
  }
  Face InternalFlipFaceIneq(Face const &sInc, const T *out) const {
    // We need to compute a vertex in the facet, but not the ridge
    size_t pos_outside = 0;
    while (true) {
      if (sInc[pos_outside] == 0)
        break;
      pos_outside++;
    }
    int outRow = PairIncs.second[pos_outside];
    T eSum(0);
    for (int iCol = 0; iCol < nbCol - 1; iCol++)
      eSum += EXT_red(outRow, iCol) * out[iCol];
    int eSign = 1;
    if (eSum < 0)
      eSign = -1;
    // F0 should be zero on the ridge
    MyVector<T> F0(nbCol - 1);
    for (int iCol = 0; iCol < nbCol - 1; iCol++)
      F0(iCol) = eSign * out[iCol];
    // The sought inequality is expressed as F0 + beta FacetIneq
    // So for all vectors v in EXT we have F0(v) + beta FacetIneq(v) >= 0
    // beta >= -F0(v) ListInvScal(v) = beta(v)
    // beta >= max beta(v)
    T beta_max(0);
    bool isAssigned = false;
    Face f_select(e_incd0);
    for (int pos_row=0; pos_row<e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      T eSum(0);
      for (int iCol = 0; iCol < nbCol - 1; iCol++)
        eSum += EXT_red(iRow, iCol) * F0(iCol);
      T beta = eSum * ListInvScal[pos_row];
      if (!isAssigned || beta > beta_max) {
        for (int k = 0; k < pos_row; k++)
          f_select[k] = 0;
        beta_max = beta;
      }
      if (beta_max == beta) {
        f_select[pos_row] = 1;
      }
      isAssigned = true;
    }
    Face fret(nbRow);
    for (int pos_row=0; pos_row<e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      fret[iRow] = f_select[pos_row];
    }
    // Now adding the points from the ridge
    boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
    while (jRow != boost::dynamic_bitset<>::npos) {
      int aRow = PairIncs.second[jRow];
      fret[aRow] = 1;
      jRow = sInc.find_next(jRow);
    }
    // returning the found facet
    return fret;
  }
  Face FlipFace(Face const &sInc) const {
#ifdef DEBUG_FLIP
    if (OneInc.count() != sInc.size()) {
      std::cerr << "Error in Flip 1\n";
      throw TerminalException{1};
    }
#endif
    size_t nb = sInc.count();
    boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
    auto f = [&](MyMatrix<T> &M, size_t eRank,
                 [[maybe_unused]] size_t iRow) -> void {
      int aRow = PairIncs.second[jRow];
      M.row(eRank) = EXT_red.row(aRow);
      jRow = sInc.find_next(jRow);
    };
    MyMatrix<T> NSP = NullspaceTrMatTarget_Kernel<T, decltype(f)>(nb, nbCol - 1, 1, f);
    return InternalFlipFaceIneq(sInc, NSP.data());
  }
  Face FlipFaceIneq(std::pair<Face, MyVector<T>> const &pair) const {
    return InternalFlipFaceIneq(pair.first, pair.second.data());
  }
};

// This is a special solution for computing the solutions.
//
// It is a specialization for the mpq_class. It uses special
// techniques for the computation of the Kernel.
//
// That is we are reducing the solution to a type Fp,
// that is a finite field computation. That solution
// is lifted and we check for its correctness.
//
// The computation uses the following types:
// * Tfast = Fp: The finite field type
// * mpz_class : The type for doing the flips
// * long: The types used for making the check that
//    the vector is in the kernel.
//
// There is a computation of bits in order to make
// sure that the computation with long will not
// overflow.
//
// Need to find a better template for the solution
#ifndef DISABLE_FP_CLASS
template <> struct FlippingFramework<mpq_class> {
private:
  using T = mpq_class;
  using Tint = mpz_class;
  using Tfast = Fp<long, 2147389441>;
  MyMatrix<T> EXT_redT; // rational type, but scaled to integer
  MyMatrix<mpz_class> EXT_red;
  MyMatrix<Tfast> EXT_fast;
  MyMatrix<long> EXT_long;
  int nbRow;
  int nbCol;
  bool try_int;
  Face OneInc;
  size_t max_bits;
  std::pair<std::vector<int>,std::vector<int>> PairIncs;
  int e_incd0;
  int e_incd1;
  std::vector<mpz_class> ListScal;

public:
  MyMatrix<T> EXT_face;
  size_t get_bit(mpz_class const& v) const {
    return mpz_sizeinbase(v.get_mpz_t(), 2);
  }
  FlippingFramework(MyMatrix<T> const &EXT, Face const &_OneInc)
    : try_int(false), OneInc(_OneInc), e_incd0(OneInc.size() - OneInc.count()), e_incd1(OneInc.count()), ListScal(e_incd0) {
    PairIncs = Dynamic_bitset_to_vectorints(OneInc);

    MyMatrix<Tint> EXT_scaled = RescaleRows(EXT);
    MyVector<Tint> FacetIneq = RescaleVec(FindFacetInequality(EXT, OneInc));
    //
    // Idx dropping for the projection
    //
    int idx_drop = 0;
    while (true) {
      if (FacetIneq(idx_drop) != 0)
        break;
      idx_drop++;
    }
    nbRow = EXT.rows();
    nbCol = EXT.cols();
    EXT_redT = MyMatrix<T>(nbRow, nbCol - 1);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      int pos = 0;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol != idx_drop) {
          EXT_redT(iRow, pos) = EXT_scaled(iRow, iCol);
          pos++;
        }
      }
    }
    EXT_red = UniversalMatrixConversion<Tint, T>(EXT_redT);
    //
    // Faster modular version of EXT_red
    //
    max_bits = 0;
    EXT_fast = MyMatrix<Tfast>(nbRow, nbCol - 1);
    EXT_long = MyMatrix<long>(nbRow, nbCol - 1);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      for (int iCol = 0; iCol < nbCol - 1; iCol++) {
        Tint const& val = EXT_red(iRow, iCol);
        max_bits = std::max(get_bit(val), max_bits);
        EXT_long(iRow, iCol) = val.get_si();
        EXT_fast(iRow, iCol) = Tfast(EXT_long(iRow, iCol));
      }
    }
    try_int = (max_bits <= 30);
    max_bits += get_bit(mpz_class(nbCol));
    //
    // Scalar products
    //
    for (int pos_row=0; pos_row<e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      Tint eSum = 0;
      for (int iCol = 0; iCol < nbCol; iCol++)
        eSum += FacetIneq(iCol) * EXT_scaled(iRow, iCol);
      ListScal[pos_row] = eSum;
    }
    //
    // Now the EXT face that is used by other procedure
    //
    EXT_face = MyMatrix<T>(e_incd1, nbCol - 1);
    for (int i_row=0; i_row<e_incd1; i_row++) {
      int j_row = PairIncs.second[i_row];
      int pos = 0;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol != idx_drop) {
          EXT_face(i_row, pos) = EXT(j_row, iCol);
          pos++;
        }
      }
    }
  }
  MyVector<Tint> RescaleVec(MyVector<T> const &v) const {
    int cols = v.size();
    std::vector<mpz_class> dens(cols, 1);
    MyVector<Tint> vret = MyVector<Tint>(cols);
    for (int iCol = 0; iCol < cols; iCol++) {
      dens[iCol] = v(iCol).get_den();
    }
    mpz_class scale = LCMlist(dens);
    for (int iCol = 0; iCol < cols; iCol++) {
      vret(iCol) = (scale / v(iCol).get_den()) * v(iCol).get_num();
    }
    return vret;
  }
  MyMatrix<Tint> RescaleRows(MyMatrix<T> const &M) const {
    int rows = M.rows();
    int cols = M.cols();
    std::vector<mpz_class> dens(cols, 1);
    MyMatrix<Tint> Mret(rows, cols);
    for (int iRow = 0; iRow < rows; iRow++) {
      for (int iCol = 0; iCol < cols; iCol++) {
        dens[iCol] = M(iRow, iCol).get_den();
      }
      mpz_class scale = LCMlist(dens);
      for (int iCol = 0; iCol < cols; iCol++) {
        Mret(iRow, iCol) =
            (scale / M(iRow, iCol).get_den()) * M(iRow, iCol).get_num();
      }
    }
    return Mret;
  }

  Face InternalFlipFaceIneq(Face const &sInc, const Tint *out) const {
    // We need to compute a vertex in the facet, but not the ridge
    size_t pos_outside = 0;
    while (true) {
      if (sInc[pos_outside] == 0)
        break;
      pos_outside++;
    }
    int outRow = PairIncs.second[pos_outside];
    Tint eSum = 0;
    for (int iCol = 0; iCol < nbCol - 1; iCol++)
      eSum += EXT_red(outRow, iCol) * out[iCol];
    int eSign = 1;
    if (eSum < 0)
      eSign = -1;
    // F0 should be zero on the ridge
    MyVector<Tint> F0(nbCol - 1);
    for (int iCol = 0; iCol < nbCol - 1; iCol++)
      F0(iCol) = eSign * out[iCol];
    Tint beta_max_num = 0;
    Tint beta_max_den = 1;
    bool isAssigned = false;
    Face f_select(e_incd0);
    for (int pos_row=0; pos_row<e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      Tint eSum = 0;
      T const& eScal = ListScal[pos_row];
      for (int iCol = 0; iCol < nbCol - 1; iCol++)
        eSum += EXT_red(iRow, iCol) * F0(iCol);
      if (!isAssigned ||
          eSum * beta_max_den < beta_max_num * eScal) {
        for (int k = 0; k < pos_row; k++)
          f_select[k] = 0;
        beta_max_num = eSum;
        beta_max_den = eScal;
      }
      if (eSum * beta_max_den == beta_max_num * eScal) {
        f_select[pos_row] = 1;
      }
      isAssigned = true;
    }
    // Now putting things together
    Face fret(nbRow);
    for (int pos_row=0; pos_row<e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      fret[iRow] = f_select[pos_row];
    }
    boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
    while (jRow != boost::dynamic_bitset<>::npos) {
      int aRow = PairIncs.second[jRow];
      fret[aRow] = 1;
      jRow = sInc.find_next(jRow);
    }
    return fret;
  }
  Face FlipFace(Face const &sInc) const {
#ifdef DEBUG_FLIP
    if (OneInc.count() != sInc.size()) {
      std::cerr << "Error in Flip 1\n";
      throw TerminalException{1};
    }
#endif
    size_t nb = sInc.count();
    MyMatrix<Tint> NSP = MyMatrix<Tint>(1, nbCol - 1);
    bool failed_int = false;
    if (try_int) {
      boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
      auto f = [&](MyMatrix<Tfast> &M, size_t eRank,
                   [[maybe_unused]] size_t iRow) -> void {
        int aRow = PairIncs.second[jRow];
        M.row(eRank) = EXT_fast.row(aRow);
        jRow = sInc.find_next(jRow);
      };
      MyMatrix<Tfast> NSP_fastT =
        NullspaceTrMatTarget_Kernel<Tfast, decltype(f)>(nb, nbCol - 1, 1, f);
      // check result at full precision in case of overflows
      bool allzero = true;
      for (int iCol = 0; iCol < nbCol - 1; iCol++) {
        if (NSP_fastT(0, iCol) != 0) {
          allzero = false;
          break;
        }
      }
      if (allzero) {
        std::cerr << "NSPint is all zero"
                  << "\n";
        failed_int = true;
      } else {
        MyVector<long> VZ_long(nbCol - 1);

        // reconstruct
        size_t max_bits_NSP = 0;
        std::vector<long> nums(nbCol - 1, 0);
        std::vector<long> dens(nbCol - 1, 1);
        for (int iCol = 0; iCol < nbCol - 1; iCol++) {
          Rational<long> val = NSP_fastT(0, iCol).rational_lift();
          nums[iCol] = val.get_num();
          dens[iCol] = val.get_den();
        }
        long lcm = LCMlist(dens);
        for (int iCol = 0; iCol < nbCol - 1; iCol++) {
          VZ_long(iCol) = nums[iCol] * (lcm / dens[iCol]);
          NSP(0, iCol) = Tint(VZ_long(iCol));
          max_bits_NSP = std::max(max_bits_NSP, get_bit(NSP(0, iCol)));
        }
        // check if elements are small enough to do computation in
        if (max_bits + max_bits_NSP <= 60) {
          // check if part of kernel
          jRow = sInc.find_first();
          for (size_t iRow = 0; iRow < nb; iRow++) {
            int aRow = PairIncs.second[jRow];
            auto row = EXT_long.row(aRow);
            jRow = sInc.find_next(jRow);
            long sm = 0;
            for (int iCol = 0; iCol < nbCol - 1; iCol++) {
              sm += VZ_long(iCol) * row(iCol);
            }
            if (sm != 0) {
              std::cerr << "Not really a kernel vector " << sm << "\n";
              failed_int = true;
              break;
            }
          }
        } else {
          std::cerr << "Precision too low" << max_bits << " " << max_bits_NSP
                    << std::endl;
          failed_int = true;
        }
      }
    }

    if (failed_int || !try_int) {
      std::cerr << "Rational<long> strategy failed , retrying with mpq_class"
                << "\n";
      boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
      auto f = [&](MyMatrix<T> &M, size_t eRank,
                   [[maybe_unused]] size_t iRow) -> void {
        int aRow = PairIncs.second[jRow];
        M.row(eRank) = EXT_redT.row(aRow);
        jRow = sInc.find_next(jRow);
      };
      NSP =
        RescaleRows(NullspaceTrMatTarget_Kernel<T, decltype(f)>(nb, nbCol - 1, 1, f));
    }
    return InternalFlipFaceIneq(sInc, NSP.data());
  }
  Face FlipFaceIneq(std::pair<Face, MyVector<T>> const &pair) const {
    MyVector<Tint> out = RescaleVec(pair.second);
    return InternalFlipFaceIneq(pair.first, out.data());
  }
};
#endif

template <typename T>
Face ComputeFlipping(MyMatrix<T> const &EXT, Face const &OneInc,
                     Face const &sInc) {
  MyMatrix<T> TheEXT = ColumnReduction(EXT);
  FlippingFramework TheFram(TheEXT, OneInc);
  return TheFram.FlipFace(sInc);
}

void PrintListOrbit(std::ostream &os, vectface const &ListOrbit) {
  size_t nbOrbit = ListOrbit.size();
  os << "nbOrbit=" << nbOrbit << "\n";
  for (size_t iOrbit = 0; iOrbit < nbOrbit; iOrbit++) {
    Face eInc = ListOrbit[iOrbit];
    size_t siz = eInc.count();
    os << "O" << iOrbit + 1 << ": inc=" << siz << "list=";
    boost::dynamic_bitset<>::size_type eVal = eInc.find_first();
    for (size_t i = 0; i < siz; i++) {
      os << " " << eVal;
      eVal = eInc.find_next(eVal);
    }
    os << "\n";
  }
}

struct EngelPolyhedralSubordination {
  int n;
  std::vector<CollectedResult<int>> TheSub;
  std::vector<vectface> ListListFace;
};

template <typename T>
EngelPolyhedralSubordination
ComputeEngelPolyhedralSubordination(MyMatrix<T> const &EXT,
                                    MyMatrix<T> const &FAC) {
  int nbFac = FAC.rows();
  int nbExt = EXT.rows();
  vectface FACset(nbExt);
  int n = FAC.cols();
  std::cerr << "nbFac=" << nbFac << " nbExt=" << nbExt << " n=" << n << "\n";
  MyVector<T> eFac;
  MyVector<T> eExt;
  for (int iFac = 0; iFac < nbFac; iFac++) {
    eFac = FAC.row(iFac);
    Face eFace(nbExt);
    //    std::cerr << "iFac=" << iFac << "\n";
    for (int iExt = 0; iExt < nbExt; iExt++) {
      eExt = EXT.row(iExt);
      T eScal = eFac.dot(eExt);
      int eValIns;
      if (eScal == 0)
        eValIns = 1;
      else
        eValIns = 0;
      eFace[iExt] = eValIns;
      //      std::cerr << "  iExt=" << iExt << " val=" << eValIns << "\n";
    }
    FACset.push_back(eFace);
  }
  std::vector<vectface> ListListFace;
  ListListFace.emplace_back(std::move(FACset));
  std::vector<CollectedResult<int>> TheSub;
  for (int eDim = 0; eDim < n - 1; eDim++) {
    int TheRank = n - 2 - eDim;
    std::cerr << "eDim=" << eDim << " TheRank=" << TheRank << "\n";
    std::vector<int> ListSizes;
    std::unordered_set<Face> NewListFace_set;
    std::cerr << "  siz=" << ListListFace[eDim].size() << "\n";
    for (auto &eFace : ListListFace[eDim]) {
      int nb = eFace.count();
      std::vector<int> eList(nb);
      boost::dynamic_bitset<>::size_type aRow = eFace.find_first();
      for (int i = 0; i < nb; i++) {
        eList[i] = static_cast<int>(aRow);
        aRow = eFace.find_next(aRow);
      }
      std::unordered_set<Face> ListSubFace;
      for (auto &fFace : FACset) {
        Face gFace(nbExt);
        std::vector<int> gList;
        int eIncd = 0;
        for (auto &eVal : eList) {
          if (fFace[eVal] == 1) {
            gList.push_back(eVal);
            gFace[eVal] = 1;
            eIncd++;
          }
        }
        bool IsFace;
        if (eIncd < TheRank) {
          IsFace = false;
        } else {
          MyMatrix<T> EXTmat = SelectRow(EXT, gList);
          int rank = RankMat(EXTmat);
          IsFace = rank == TheRank;
        }
        if (IsFace)
          ListSubFace.insert(gFace);
      }
      int eSize = ListSubFace.size();
      for (auto &rgFace : ListSubFace)
        NewListFace_set.insert(rgFace);
      ListSizes.push_back(eSize);
    }
    TheSub.push_back(Collected(ListSizes));
    vectface NewListFace_vect(nbExt);
    for (auto &rFace : NewListFace_set)
      NewListFace_vect.push_back(rFace);
    ListListFace.emplace_back(std::move(NewListFace_vect));
  }
  return {n, std::move(TheSub), std::move(ListListFace)};
}

template <typename T>
void ComputeFileFaceLatticeInfo(std::string const &eFile,
                                MyMatrix<T> const &EXT,
                                MyMatrix<T> const &FAC) {
  EngelPolyhedralSubordination eEngel =
      ComputeEngelPolyhedralSubordination(EXT, FAC);
  std::ofstream os(eFile);
  os << "return [";
  int len = eEngel.ListListFace.size();
  int nbExt = EXT.rows();
  for (int iDim = 0; iDim < len; iDim++) {
    if (iDim > 0)
      os << ",\n";
    int nbFace = eEngel.ListListFace[iDim].size();
    os << "[";
    for (int iFace = 0; iFace < nbFace; iFace++) {
      if (iFace > 0)
        os << ",";
      Face eFace = eEngel.ListListFace[iDim][iFace];
      bool IsFirst = true;
      os << "[";
      for (int iExt = 0; iExt < nbExt; iExt++) {
        if (eFace[iExt] == 1) {
          if (!IsFirst)
            os << ",";
          IsFirst = false;
          int eVal = iExt + 1;
          os << eVal;
        }
      }
      os << "]";
    }
    os << "]";
  }
  os << "];\n";
}

template <typename T>
void ComputeEngelPolyhedralSubordinationFile(std::string const &eFile,
                                             MyMatrix<T> const &EXT,
                                             MyMatrix<T> const &FAC) {
  EngelPolyhedralSubordination eEngel =
      ComputeEngelPolyhedralSubordination(EXT, FAC);
  std::ofstream os(eFile);
  os << "return [";
  int len = eEngel.TheSub.size();
  for (int i = 0; i < len; i++) {
    if (i > 0)
      os << ",\n";
    CollectedResult<int> eColl = eEngel.TheSub[i];
    int nbEnt = eColl.LVal.size();
    os << "[";
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      int eVal = eColl.LVal[iEnt];
      int eMult = eColl.LMult[iEnt];
      if (iEnt > 0)
        os << ",";
      os << "[" << eVal << "," << eMult << "]";
    }
    os << "]";
  }
  os << "];\n";
}

MyMatrix<int> VectfaceAsMatrix(vectface const &vf) {
  size_t n_ent = vf.size();
  size_t n = vf.get_n();
  MyMatrix<int> M(n_ent, n);
  for (size_t i_ent = 0; i_ent < n_ent; i_ent++) {
    Face f = vf[i_ent];
    for (size_t i = 0; i < n; i++)
      M(i_ent, i) = f[i];
  }
  return M;
}

// clang-format off
#endif  // SRC_POLY_POLY_POLYTOPEFCT_H_
// clang-format on
