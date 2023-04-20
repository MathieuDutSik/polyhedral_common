// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_POLYTOPEFCT_H_
#define SRC_POLY_POLY_POLYTOPEFCT_H_

#include "Boost_bitset.h"
#include "COMB_Stor.h"
#include "MAT_Matrix.h"
#include "rational.h"
#include "Fp.h"
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "NumberTheoryGeneric.h"

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

template <typename T>
void TestFacetness(MyMatrix<T> const &EXT, Face const &eList) {
  MyMatrix<T> TheEXT = ColumnReduction(EXT);
  int nb = eList.count();
  int nbRow = TheEXT.rows();
  int nbCol = TheEXT.cols();
  MyMatrix<T> TheProv(nb, nbCol);
  boost::dynamic_bitset<>::size_type aRow = eList.find_first();
  for (int iRow = 0; iRow < nb; iRow++) {
    TheProv.row(iRow) = TheEXT.row(aRow);
    aRow = eList.find_next(aRow);
  }
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(TheProv);
  MyMatrix<T> NSP = eSelect.NSP;
  if (NSP.rows() != 1) {
    std::cerr << "Error in rank in Facetness\n";
    std::cerr << "|NSP|=" << NSP.rows() << "\n";
    throw TerminalException{1};
  }
  int nbZero = 0;
  int nbPlus = 0;
  int nbMinus = 0;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eScal = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      eScal += NSP(0, iCol) * TheEXT(iRow, iCol);
    if (eScal == 0)
      nbZero++;
    if (eScal > 0)
      nbPlus++;
    if (eScal < 0)
      nbMinus++;
  }
  if (nbZero == EXT.rows()) {
    std::cerr << "All vectors seems to be incident. And that is not allowed "
                 "for a facet\n";
    throw TerminalException{1};
  }
  if (nbZero != nb) {
    std::cerr << "Error in computing incidence\n";
    std::cerr << "nbZero=" << nbZero << " nb=" << nb << "\n";
    std::cerr << "nbMinus=" << nbMinus << " nbPlus=" << nbPlus << "\n";
    std::cerr << "EXT(rows/cols)=" << EXT.rows() << " / " << EXT.cols() << "\n";
    std::cerr << "Rank(EXT)=" << RankMat(EXT) << "\n";
    throw TerminalException{1};
  }
  if (nbMinus > 0 && nbPlus > 0) {
    std::cerr << "Some plus and minus signs, illegal\n";
    throw TerminalException{1};
  }
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
  MyMatrix<T> NSP = NullspaceTrMat_Kernel<T, decltype(f)>(nb, nbCol, f);
  if (NSP.rows() != 1) {
    std::cerr << "We should have just one row in NSP\n";
    throw TerminalException{1};
  }
  MyVector<T> eVect(nbCol);
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    eVect(iCol) = NSP(0, iCol);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    if (OneInc[iRow])
      continue;
    T eScal = 0;
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      eScal += eVect(iCol) * TheEXT(iRow, iCol);
    if (eScal > 0)
      return eVect;
    if (eScal < 0)
      return -eVect;
  }
  std::cerr << "Error in FindFacetInequality\n";
  throw TerminalException{1};
}

std::vector<int> Dynamic_bitset_to_vectorint(Face const &eList) {
  int nb = eList.count();
  boost::dynamic_bitset<>::size_type aRow = eList.find_first();
  std::vector<int> retList(nb);
  for (int i = 0; i < nb; i++) {
    retList[i] = static_cast<int>(aRow);
    aRow = eList.find_next(aRow);
  }
  return retList;
}

template <typename T> struct FlippingFramework {
private:
  MyMatrix<T> EXT_red;
  int nbRow;
  int nbCol;
  Face OneInc;
  std::vector<int> OneInc_V;
  std::vector<T> ListInvScal;

public:
  MyMatrix<T> EXT_face;
  FlippingFramework(MyMatrix<T> const &EXT, Face const &_OneInc)
      : OneInc(_OneInc) {
    OneInc_V = Dynamic_bitset_to_vectorint(OneInc);
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
    ListInvScal = std::vector<T>(nbRow, 0);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      if (OneInc[iRow] == 0) {
        T eSum = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          eSum += FacetIneq(iCol) * EXT(iRow, iCol);
        ListInvScal[iRow] = -1 / eSum;
      }
    }
    //
    // Now the EXT face that is used by other procedure
    //
    size_t e_incd = OneInc.count();
    EXT_face = MyMatrix<T>(e_incd, nbCol - 1);
    boost::dynamic_bitset<>::size_type j_row = OneInc.find_first();
    int i_row = 0;
    while (j_row != boost::dynamic_bitset<>::npos) {
      int pos = 0;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol != idx_drop) {
          EXT_face(i_row, pos) = EXT(j_row, iCol);
          pos++;
        }
      }
      j_row = OneInc.find_next(j_row);
      i_row++;
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
    int outRow = OneInc_V[pos_outside];
    T eSum = 0;
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
    // So for all vectors v in EXT we have F0(v) + beta FacetInea(v) >= 0
    // beta >= -F0(v) ListInvScal(v) = beta(v)
    // beta >= max beta(v)
    Face fret(nbRow);
    T beta_max = 0;
    bool isAssigned = false;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      if (OneInc[iRow] == 0) {
        T eSum = 0;
        for (int iCol = 0; iCol < nbCol - 1; iCol++)
          eSum += EXT_red(iRow, iCol) * F0(iCol);
        T beta = eSum * ListInvScal[iRow];
        if (!isAssigned || beta > beta_max) {
          for (int kRow = 0; kRow < iRow; kRow++)
            fret[kRow] = 0;
          beta_max = beta;
        }
        if (beta_max == beta) {
          fret[iRow] = 1;
        }
        isAssigned = true;
      }
    }
    // Now adding the points from the ridge
    boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
    while (jRow != boost::dynamic_bitset<>::npos) {
      int aRow = OneInc_V[jRow];
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
      int aRow = OneInc_V[jRow];
      M.row(eRank) = EXT_red.row(aRow);
      jRow = sInc.find_next(jRow);
    };
    MyMatrix<T> NSP = NullspaceTrMat_Kernel<T, decltype(f)>(nb, nbCol - 1, f);
#ifdef DEBUG_FLIP
    if (NSP.rows() != 1) {
      std::cerr << "Error in Flip 2\n";
      throw TerminalException{1};
    }
#endif
    return InternalFlipFaceIneq(sInc, NSP.data());
  }
  Face FlipFaceIneq(std::pair<Face, MyVector<T>> const &pair) const {
    return InternalFlipFaceIneq(pair.first, pair.second.data());
  }
};

// Need to find a better template for the solution
#ifndef DISABLE_MPQ_CLASS
template <> struct FlippingFramework<mpq_class> {
private:
  using T = mpq_class;
  using Tint = mpz_class;
  using Tfast = Fp<long, 2147389441>;
  MyMatrix<T> EXT_redT; // rational type, but scaled to integer
  MyMatrix<mpz_class> EXT_red;
  MyMatrix<Tfast> EXT_fastT;
  MyMatrix<long> EXT_fast;
  int nbRow;
  int nbCol;
  bool try_int;
  Face OneInc;
  size_t max_bits;
  std::vector<int> OneInc_V;
  std::vector<mpz_class> ListScal;

public:
  MyMatrix<T> EXT_face;
  FlippingFramework(MyMatrix<T> const &EXT, Face const &_OneInc)
      : try_int(false), OneInc(_OneInc) {
    OneInc_V = Dynamic_bitset_to_vectorint(OneInc);

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
    EXT_fastT = MyMatrix<Tfast>(nbRow, nbCol - 1);
    EXT_fast  = MyMatrix<long>(nbRow, nbCol - 1);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      for (int iCol = 0; iCol < nbCol - 1; iCol++) {
        max_bits = std::max(mpz_sizeinbase(EXT_red(iRow, iCol).get_mpz_t(), 2), max_bits);
        EXT_fast(iRow, iCol) = EXT_red(iRow, iCol).get_si();
        EXT_fastT(iRow, iCol) = Tfast(EXT_fast(iRow, iCol));
      }
    }
    try_int = (max_bits <= 30);
    max_bits += mpz_sizeinbase(mpz_class(nbCol).get_mpz_t(), 2);
    //
    // Scalar products
    //
    ListScal = std::vector<Tint>(nbRow, 0);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      if (OneInc[iRow] == 0) {
        Tint eSum = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          eSum += FacetIneq(iCol) * EXT_scaled(iRow, iCol);
        ListScal[iRow] = eSum;
      }
    }
    //
    // Now the EXT face that is used by other procedure
    //
    size_t e_incd = OneInc.count();
    EXT_face = MyMatrix<T>(e_incd, nbCol - 1);
    boost::dynamic_bitset<>::size_type j_row = OneInc.find_first();
    int i_row = 0;
    while (j_row != boost::dynamic_bitset<>::npos) {
      int pos = 0;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol != idx_drop) {
          EXT_face(i_row, pos) = EXT(j_row, iCol);
          pos++;
        }
      }
      j_row = OneInc.find_next(j_row);
      i_row++;
    }
  }
  MyVector<Tint> RescaleVec(MyVector<T> const &v) const {
    int cols = v.size();
    std::vector<mpz_class> dens(cols,1);
    MyVector<Tint> vret = MyVector<Tint>(cols);
    for( int iCol = 0; iCol < cols; iCol++){
      dens[iCol] = v(iCol).get_den();
    }
    mpz_class scale = LCMlist(dens);
    for( int iCol = 0; iCol < cols; iCol++) {
      vret(iCol) = (scale / v(iCol).get_den()) * v(iCol).get_num();
    }
    return vret;
  }
  MyMatrix<Tint> RescaleRows(MyMatrix<T> const &M) const {
    int rows = M.rows();
    int cols = M.cols();
    std::vector<mpz_class> dens(cols,1);
    MyMatrix<Tint> Mret = MyMatrix<Tint>(rows, cols);
    for( int iRow = 0; iRow < rows; iRow++) {
      for( int iCol = 0; iCol < cols; iCol++){
        dens[iCol] = M(iRow, iCol).get_den();
      }
      mpz_class scale = LCMlist(dens);
      for( int iCol = 0; iCol < cols; iCol++) {
        Mret(iRow, iCol) = (scale / M(iRow,iCol).get_den()) * M(iRow,iCol).get_num();
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
    int outRow = OneInc_V[pos_outside];
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
    // The sought inequality is expressed as F0 + beta FacetIneq
    // So for all vectors v in EXT we have F0(v) + beta FacetInea(v) >= 0
    // beta >= -F0(v) ListInvScal(v) = beta(v)
    // beta >= max beta(v)
    Face fret(nbRow);
    Tint beta_max_num = 0;
    Tint beta_max_den = 1;
    bool isAssigned = false;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      if (OneInc[iRow] == 0) {
        Tint eSum = 0;
        for (int iCol = 0; iCol < nbCol - 1; iCol++)
          eSum += EXT_red(iRow, iCol) * F0(iCol);
        if (!isAssigned || eSum * beta_max_den < beta_max_num * ListScal[iRow]) {
          for (int kRow = 0; kRow < iRow; kRow++)
            fret[kRow] = 0;
          beta_max_num = eSum;
          beta_max_den = ListScal[iRow];
        }
        if (eSum * beta_max_den == beta_max_num * ListScal[iRow]) {
          fret[iRow] = 1;
        }
        isAssigned = true;
      }
    }
    // Now adding the points from the ridge
    boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
    while (jRow != boost::dynamic_bitset<>::npos) {
      int aRow = OneInc_V[jRow];
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
    MyMatrix<Tint> NSP = MyMatrix<Tint>(1, nbCol - 1);
    bool failed_int = false;
    if (try_int) {
      boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
      auto f = [&](MyMatrix<Tfast> &M, size_t eRank,
                   [[maybe_unused]] size_t iRow) -> void {
        int aRow = OneInc_V[jRow];
        M.row(eRank) = EXT_fastT.row(aRow);
        jRow = sInc.find_next(jRow);
      };
      MyMatrix<Tfast> NSP_fastT = NullspaceTrMat_Kernel<Tfast, decltype(f)>(nb, nbCol - 1,f);
      // check result at full precision in case of overflows
      if (NSP_fastT.rows() != 1) {
        std::cerr << "NSP_fastT.rows() != 1"
                  << "\n";
        failed_int = true;
      } else {
        bool allzero = true;
        for (int iCol = 0; iCol < nbCol - 1; iCol++) {
          if (NSP_fastT(0, iCol) != 0){
            allzero = false;
            break;
          }
        }
        if (allzero) {
          std::cerr << "NSPint is all zero"
                    << "\n";
          failed_int = true;
        } else {
          MyMatrix<long> NSP_fast = MyMatrix<long>(1, nbCol - 1);

          // reconstruct
          size_t max_bits_NSP = 0;
          std::vector<long> nums(nbCol-1,0);
          std::vector<long> dens(nbCol-1,1);
          for (int iCol = 0; iCol < nbCol - 1; iCol++) {
            Rational<long> val = NSP_fastT(0, iCol).rational_lift();
            nums[iCol] = val.get_num();
            dens[iCol] = val.get_den();
          }
          long lcm = LCMlist(dens);
          for (int iCol = 0; iCol < nbCol - 1; iCol++) {
            NSP_fast(0, iCol) = nums[iCol] * (lcm / dens[iCol]);
            NSP(0, iCol) = Tint(NSP_fast(0,iCol));
            max_bits_NSP = std::max(max_bits_NSP, mpz_sizeinbase(NSP(0,iCol).get_mpz_t(), 2));
          }


          // check if elements are small enough to do computation in
          if (max_bits + max_bits_NSP <= 60) {
            // check if part of kernel
            jRow = sInc.find_first();
            for (size_t iRow = 0; iRow < nb; iRow++) {
              int aRow = OneInc_V[jRow];
              auto row = EXT_fast.row(aRow);
              jRow = sInc.find_next(jRow);
              long sm = 0;
              for (int iCol = 0; iCol < nbCol - 1; iCol++) {
                sm += NSP_fast(0, iCol) * row(iCol);
              }
              if (sm != 0) {
                std::cerr << "Not really a kernel vector " << sm << "\n";
                failed_int = true;
                break;
              }
            }
          } else {
            std::cerr << "Precision too low" << max_bits << " " << max_bits_NSP << std::endl;
            failed_int = true;
          }
        }
      }
    }

    if (failed_int || !try_int) {
      std::cerr << "Rational<long> strategy failed , retrying with mpq_class"
                << "\n";
      boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
      auto f = [&](MyMatrix<T> &M, size_t eRank,
                   [[maybe_unused]] size_t iRow) -> void {
        int aRow = OneInc_V[jRow];
        M.row(eRank) = EXT_redT.row(aRow);
        jRow = sInc.find_next(jRow);
      };
      NSP = RescaleRows(NullspaceTrMat_Kernel<T, decltype(f)>(nb, nbCol - 1, f));
    }
#ifdef DEBUG_FLIP
    if (NSP.rows() != 1) {
      std::cerr << "Error in Flip 2\n";
      throw TerminalException{1};
    }
#endif
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
  int nbRow = TheEXT.rows();
  int nbCol = TheEXT.cols();
  vectface TwoPlanes(nbRow);
  T eVal, prov1, prov2, prov3, prov4;
  T EXT1_1, EXT1_2, EXT2_1, EXT2_2;
  T TheDet, det12, det1N, det2N, prodDet, h;
  int nbForm, IsNonZero, i;
  //  std::cerr << "Begining of ComputeFlipping\n";
  //  std::cerr << "OneInc(count/size)=" << OneInc.count() << "/" <<
  //  OneInc.size() << "\n"; std::cerr << "sInc(count/size)=" << sInc.count() <<
  //  "/" << sInc.size() << "\n"; std::cerr << "RankMat(TheEXT)=" <<
  //  RankMat(TheEXT) << "\n";
  if (OneInc.count() != sInc.size()) {
    std::cerr << "Error in ComputeFlipping\n";
    throw TerminalException{1};
  }
  int nb = sInc.count();
  MyMatrix<T> TheProv(nb, nbCol);
  MyMatrix<T> LV(2, 2);
  TestFacetness(TheEXT, OneInc);
  std::vector<int> OneInc_V = Dynamic_bitset_to_vectorint(OneInc);
  boost::dynamic_bitset<>::size_type jRow = sInc.find_first();
  for (int iRow = 0; iRow < nb; iRow++) {
    int aRow = OneInc_V[jRow];
    TheProv.row(iRow) = TheEXT.row(aRow);
    jRow = sInc.find_next(jRow);
  }
  //  std::cerr << "RankMat(TheProv)=" << RankMat(TheProv) << "\n";
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(TheProv);
  MyMatrix<T> NSP = eSelect.NSP;
  if (NSP.rows() != 2) {
    int nbRowNSP = NSP.rows();
    std::cerr << "NSP.nbRows=" << nbRowNSP << "\n";
    std::cerr << "Deep inconsistency in ComputeFlipping\n";
    throw TerminalException{1};
  }
  MyMatrix<T> NSPtrans = TransposedMat(NSP);
  MyMatrix<T> LProd = TheEXT * NSPtrans;
  nbForm = 0;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    IsNonZero = 0;
    for (i = 0; i < 2; i++) {
      prov1 = LProd(iRow, i);
      if (prov1 != 0)
        IsNonZero = 1;
    }
    if (IsNonZero == 1) {
      prov1 = LProd(iRow, 0);
      prov2 = LProd(iRow, 1);
      if (nbForm == 0) {
        EXT1_1 = prov1;
        EXT1_2 = prov2;
        nbForm++;
      } else {
        TheDet = prov2 * EXT1_1 - prov1 * EXT1_2;
        if (nbForm == 1) {
          if (TheDet != 0) {
            EXT2_1 = prov1;
            EXT2_2 = prov2;
            det12 = TheDet;
            nbForm++;
          }
        } else {
          det1N = TheDet;
          det2N = prov2 * EXT2_1 - prov1 * EXT2_2;
          prodDet = det1N * det2N;
          if (prodDet > 0) {
            prodDet = det12 * det1N;
            if (prodDet > 0) {
              EXT2_1 = prov1;
              EXT2_2 = prov2;
              det12 = det1N;
            } else {
              EXT1_1 = prov1;
              EXT1_2 = prov2;
              det12 = -det2N;
            }
          }
        }
      }
    }
  }
  if (det12 > 0) {
    eVal = -EXT1_2;
    LV(0, 0) = eVal;
    LV(1, 0) = EXT1_1;
    LV(0, 1) = EXT2_2;
    eVal = -EXT2_1;
    LV(1, 1) = eVal;
  } else {
    LV(0, 0) = EXT1_2;
    eVal = -EXT1_1;
    LV(1, 0) = eVal;
    eVal = -EXT2_2;
    LV(0, 1) = eVal;
    LV(1, 1) = EXT2_1;
  }
  MyMatrix<T> PairFac = NSPtrans * LV;
  MyMatrix<T> ListScal = TheEXT * PairFac;
  for (i = 0; i < 2; i++) {
    Face ePlane(nbRow);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      eVal = ListScal(iRow, i);
      if (eVal < 0) {
        std::cerr << "This should never have happened. Please panic\n";
        throw TerminalException{1};
      }
      if (eVal == 0)
        ePlane[iRow] = 1;
    }
    TwoPlanes.push_back(ePlane);
  }
  int idxFound = -1;
  for (int k = 0; k < 2; k++)
    if (TwoPlanes[k] == OneInc)
      idxFound = k;
  if (idxFound == -1) {
    std::cerr << "We did not find the facet\n";
    throw TerminalException{1};
  }
  return TwoPlanes[1 - idxFound];
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

// clang-format off
#endif  // SRC_POLY_POLY_POLYTOPEFCT_H_
// clang-format on
