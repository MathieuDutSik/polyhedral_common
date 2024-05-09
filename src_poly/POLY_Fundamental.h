/// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_FUNDAMENTAL_H_
#define SRC_POLY_POLY_FUNDAMENTAL_H_

// clang-format off
#include "Boost_bitset.h"
#include "COMB_Stor.h"
#include "MAT_Matrix.h"
#include "NumberTheoryGeneric.h"
#include "MAT_Matrix_SubsetSolver.h"
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <algorithm>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLY_FUNDAMENTAL
#endif

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

// Compute an upper bound on the determinant of maximal minor
// The Hadamard bound is
// We have det(A)^2 <= Pi_{i=1}^n (sum_j x_{ij}^2)
template <typename T> T sqr_estimate_maximal_determinant(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  std::vector<T> ListSqr;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T sum = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      sum += M(iRow, iCol) * M(iRow, iCol);
    ListSqr.push_back(sum);
  }
  std::sort(ListSqr.begin(), ListSqr.end());
  T eProd = 1;
  for (int iRow = nbRow - nbCol; iRow < nbRow; iRow++) {
    eProd *= ListSqr[iRow];
  }
  return eProd;
}

template <typename T> T sqr_estimate_facet_coefficients(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  T max_coeff = 0;
  MyVector<T> SqrCoeffFacet(nbCol);
  for (int jCol = 0; jCol < nbCol; jCol++) {
    MyMatrix<T> Mret(nbRow, nbCol - 1);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      int pos = 0;
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (iCol != jCol) {
          Mret(iRow, pos) = M(iRow, iCol);
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
void CheckFacetInequality(MyMatrix<T> const &EXT, Face const &eList, std::string const& context) {
  int siz = eList.size();
  int nb = eList.count();
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  if (nbRow != siz) {
    std::cerr << "nbRow=" << nbRow << " siz=" << siz << " which is inconsistent\n";
    std::cerr << "context=" << context << "\n";
    throw TerminalException{1};
  }
  boost::dynamic_bitset<>::size_type aRow = eList.find_first();
  auto f = [&](MyMatrix<T> &M, size_t eRank,
               [[maybe_unused]] size_t iRow) -> void {
    M.row(eRank) = EXT.row(aRow);
    aRow = eList.find_next(aRow);
  };
  MyMatrix<T> NSP = NullspaceTrMat_Kernel<T, decltype(f)>(nb, nbCol, f);
  if (NSP.rows() != 1) {
    std::cerr << "Error in rank in Facetness\n";
    std::cerr << "|EXT|=" << nbRow << " / " << nbCol << "\n";
    std::cerr << "|eList|=" << eList.size() << " / " << eList.count() << "\n";
    std::cerr << "|NSP|=" << NSP.rows() << " when it should be 1\n";
    std::cerr << "context=" << context << "\n";
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
        std::cerr << "The vertex iRow has a zero scalar product but does not "
                     "belong to eList\n";
        nbError++;
      }
    } else {
      if (eList[iRow] != 0) {
        std::cerr << "The vertex iRow has a non-zero scalar product but does "
                     "belong to eList\n";
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
    std::cerr << "nbError=" << nbError << "\n";
    std::cerr << "context=" << context << "\n";
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
  MyMatrix<T> NSP =
      NullspaceTrMatTarget_Kernel<T, decltype(f)>(nb, nbCol, 1, f);
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
  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  std::cerr << "|f|=" << OneInc.size() << " / " << OneInc.count() << "\n";
  std::cerr << "eVect=" << StringVectorGAP(eVect) << "\n";
  std::cerr << "|NSP|=" << NSP.rows() << " / " << NSP.cols() << "\n";
  std::cerr << "Rank=" << RankMat(TheEXT) << "\n";
  std::cerr << "TheEXT=\n";
  WriteMatrix(std::cerr, TheEXT);
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

template <typename T>
MyMatrix<T> GetEXT_face(MyMatrix<T> const &EXT, int const &idx_drop,
                        std::vector<int> const &l_idx) {
  int nbCol = EXT.cols();
  int e_incd = l_idx.size();
  MyMatrix<T> EXT_face(e_incd, nbCol - 1);
  for (int i_row = 0; i_row < e_incd; i_row++) {
    int j_row = l_idx[i_row];
    int pos = 0;
    for (int iCol = 0; iCol < nbCol; iCol++) {
      if (iCol != idx_drop) {
        EXT_face(i_row, pos) = EXT(j_row, iCol);
        pos++;
      }
    }
  }
  return EXT_face;
}

std::pair<std::vector<int>, std::vector<int>>
Dynamic_bitset_to_vectorints(Face const &eList) {
  int len = eList.size();
  std::vector<int> V0;
  std::vector<int> V1;
  for (int i = 0; i < len; i++) {
    if (eList[i]) {
      V1.push_back(i);
    } else {
      V0.push_back(i);
    }
  }
  return {std::move(V0), std::move(V1)};
}

template <typename T> int get_idx_drop(MyVector<T> const &V) {
  int idx_drop = 0;
  while (true) {
    if (V(idx_drop) != 0)
      break;
    idx_drop++;
  }
  return idx_drop;
}

size_t get_pos_outside(Face const &f) {
  size_t pos_outside = 0;
  while (true) {
    if (f[pos_outside] == 0)
      break;
    pos_outside++;
  }
  return pos_outside;
}

Face get_fret(std::pair<std::vector<int>, std::vector<int>> const &PairIncs,
              int const &nbRow, Face const &sInc, Face const &f_select) {
  Face fret(nbRow);
  for (size_t pos_row = 0; pos_row < PairIncs.first.size(); pos_row++) {
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
template <typename T> struct FlippingFramework_Field {
private:
  using Tint = typename SubsetRankOneSolver<T>::Tint;
  MyMatrix<T> EXT_red;
  Face OneInc;
  int e_incd0;
  int e_incd1;
  int nbRow;
  int nbCol;
  std::pair<std::vector<int>, std::vector<int>> PairIncs;
  std::vector<T> ListInvScal;
  Face f_select;
  std::ostream &os;
#ifdef DEBUG_POLY_FUNDAMENTAL
  MyMatrix<T> EXT_debug;
#endif

public:
  MyMatrix<T> EXT_face;
  MyMatrix<Tint> EXT_face_int;
  FlippingFramework_Field(MyMatrix<T> const &EXT,
                          [[maybe_unused]] MyMatrix<Tint> const &EXT_int,
                          Face const &_OneInc, std::ostream &_os)
      : OneInc(_OneInc), e_incd0(OneInc.size() - OneInc.count()),
        e_incd1(OneInc.count()), nbRow(EXT.rows()), nbCol(EXT.cols()),
        PairIncs(Dynamic_bitset_to_vectorints(OneInc)), ListInvScal(e_incd0),
        f_select(e_incd0), os(_os) {
#ifdef DEBUG_POLY_FUNDAMENTAL
    EXT_debug = EXT;
#endif
    MyVector<T> FacetIneq = FindFacetInequality(EXT, OneInc);
    //
    // Idx dropping for the projection
    //
    int idx_drop = get_idx_drop(FacetIneq);
    EXT_red = DropColumn(EXT, idx_drop);
    //
    // Inverse scalar products
    //
    for (int pos_row = 0; pos_row < e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      T eSum(0);
      for (int iCol = 0; iCol < nbCol; iCol++)
        eSum += FacetIneq(iCol) * EXT(iRow, iCol);
      ListInvScal[pos_row] = -1 / eSum;
    }
    //
    // Now the EXT face that is used by other procedure
    //
    EXT_face = GetEXT_face(EXT, idx_drop, PairIncs.second);
    EXT_face_int = GetEXT_face(EXT_int, idx_drop, PairIncs.second);
  }
  Face InternalFlipFaceIneq(Face const &sInc, const T *out) {
    // We need to compute a vertex in the facet, but not the ridge
    size_t pos_outside = get_pos_outside(sInc);
    int outRow = PairIncs.second[pos_outside];
    T eSum(0);
    for (int iCol = 0; iCol < nbCol - 1; iCol++)
      eSum += EXT_red(outRow, iCol) * out[iCol];
    int eSign = 1;
    if (eSum < 0)
      eSign = -1;
    // The sought inequality is expressed as F0 + beta FacetIneq
    // So for all vectors v in EXT we have F0(v) + beta FacetIneq(v) >= 0
    // beta >= -F0(v) ListInvScal(v) = beta(v)
    // beta >= max beta(v)
    T beta_max(0);
    bool isAssigned = false;
    for (int k = 0; k < e_incd0; k++)
      f_select[k] = 0;
    for (int pos_row = 0; pos_row < e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      T eSum(0);
      for (int iCol = 0; iCol < nbCol - 1; iCol++)
        eSum += EXT_red(iRow, iCol) * out[iCol];
      T beta = eSum * ListInvScal[pos_row];
      auto f_comp = [&]() -> bool {
        if (eSign == 1)
          return beta > beta_max;
        return beta < beta_max;
      };
      if (!isAssigned || f_comp()) {
        for (int k = 0; k < pos_row; k++)
          f_select[k] = 0;
        beta_max = beta;
        f_select[pos_row] = 1;
      } else {
        if (beta_max == beta) {
          f_select[pos_row] = 1;
        }
      }
      isAssigned = true;
    }
    Face fret = get_fret(PairIncs, nbRow, sInc, f_select);
#ifdef DEBUG_POLY_FUNDAMENTAL
    os << "FLIP: FlippingFramework_Field<T> before check\n";
    CheckFacetInequality(EXT_debug, fret, "InternalFlipFaceIneq 1");
#endif
    return fret;
  }
  Face FlipFace(Face const &sInc) {
#ifdef DEBUG_POLY_FUNDAMENTAL
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
    MyMatrix<T> NSP =
        NullspaceTrMatTarget_Kernel<T, decltype(f)>(nb, nbCol - 1, 1, f);
    return InternalFlipFaceIneq(sInc, NSP.data());
  }
  Face FlipFaceIneq(std::pair<Face, MyVector<T>> const &pair) {
    return InternalFlipFaceIneq(pair.first, pair.second.data());
  }
};

// This is a special solution for computing the solutions.
//
// It is using a special scheme for computing the solution
// that uses the Fp class. It essentially works only for the
// rational case and uses a reduction scheme that is specific
// for rational.
//
// The scheme is implemented for the mpq_class and the
// Rational<SafeInt64>. This scheme cannot be used for the
// Quadratic field, algebraic fields and so on.
template <typename T> struct FlippingFramework_Accelerate {
private:
  using Tint = typename SubsetRankOneSolver<T>::Tint;
  Face OneInc;
  int e_incd0;
  int e_incd1;
  int nbRow;
  int nbCol;
  std::vector<Tint> ListScal;
  Face f_select;
  std::pair<std::vector<int>, std::vector<int>> PairIncs;
  MyVector<Tint> FacetIneq;
  int idx_drop;
  MyMatrix<Tint> EXT_red;
  MyMatrix<Tint> EXT_red_sub;
  SubsetRankOneSolver_Acceleration<T> solver;
#ifdef DEBUG_POLY_FUNDAMENTAL
  MyMatrix<T> EXT_debug;
#endif

public:
  MyMatrix<T> EXT_face;
  MyMatrix<Tint> EXT_face_int;
  std::ostream &os;
  FlippingFramework_Accelerate(MyMatrix<T> const &EXT,
                               MyMatrix<Tint> const &EXT_int,
                               Face const &_OneInc, std::ostream &_os)
      : OneInc(_OneInc), e_incd0(OneInc.size() - OneInc.count()),
        e_incd1(OneInc.count()), nbRow(EXT.rows()), nbCol(EXT.cols()),
        ListScal(e_incd0), f_select(e_incd0),
        PairIncs(Dynamic_bitset_to_vectorints(OneInc)),
        FacetIneq(NonUniqueRescaleVecRing(FindFacetInequality(EXT, OneInc))),
        idx_drop(get_idx_drop(FacetIneq)),
        EXT_red(DropColumn(EXT_int, idx_drop)),
        EXT_red_sub(SelectRow(EXT_red, PairIncs.second)),
        solver(SubsetRankOneSolver_Acceleration<T>(EXT_red_sub)),
        EXT_face(GetEXT_face(EXT, idx_drop, PairIncs.second)),
        EXT_face_int(GetEXT_face(EXT_int, idx_drop, PairIncs.second)), os(_os) {
#ifdef DEBUG_POLY_FUNDAMENTAL
    EXT_debug = EXT;
#endif
    //
    // Scalar products
    //
    for (int pos_row = 0; pos_row < e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      Tint eSum = 0;
      for (int iCol = 0; iCol < nbCol; iCol++)
        eSum += FacetIneq(iCol) * EXT_int(iRow, iCol);
      ListScal[pos_row] = eSum;
    }
  }
  Face InternalFlipFaceIneq(Face const &sInc, MyVector<Tint> const &V) {
    // We need to compute a vertex in the facet, but not the ridge
    size_t pos_outside = get_pos_outside(sInc);
    Tint eSum = 0;
    for (int iCol = 0; iCol < nbCol - 1; iCol++)
      eSum += EXT_red_sub(pos_outside, iCol) * V(iCol);
    int eSign = 1;
    if (eSum < 0)
      eSign = -1;
    Tint beta_max_num = 0; // Those values are arbitrary and put
    Tint beta_max_den = 1; // only to avoid compiler warnings.
    bool isAssigned = false;
    Tint delta;
    for (int k = 0; k < e_incd0; k++)
      f_select[k] = 0;
    for (int pos_row = 0; pos_row < e_incd0; pos_row++) {
      int iRow = PairIncs.first[pos_row];
      Tint eSum = 0;
      Tint const &eScal = ListScal[pos_row];
      for (int iCol = 0; iCol < nbCol - 1; iCol++)
        eSum += EXT_red(iRow, iCol) * V(iCol);
      delta = beta_max_num * eScal - eSum * beta_max_den;
      auto f_comp = [&]() -> bool {
        if (eSign == 1)
          return delta > 0;
        return delta < 0;
      };
      if (!isAssigned || f_comp()) {
        for (int k = 0; k < pos_row; k++)
          f_select[k] = 0;
        beta_max_num = eSum;
        beta_max_den = eScal;
        f_select[pos_row] = 1;
      } else {
        if (delta == 0) {
          f_select[pos_row] = 1;
        }
      }
      isAssigned = true;
    }
    Face fret = get_fret(PairIncs, nbRow, sInc, f_select);
#ifdef DEBUG_POLY_FUNDAMENTAL
    os << "FLIP: |PairIncs|=" << PairIncs.first.size() << " / "
       << PairIncs.second.size() << "\n";
    os << "FLIP: OneInc=" << OneInc.size() << " / " << OneInc.count() << "\n";
    os << "FLIP: f_select=" << f_select.size() << " / " << f_select.count() << "\n";
    os << "FLIP: sInc=" << sInc.size() << " / " << sInc.count() << " eSign=" << eSign
       << "\n";
    os << "FLIP: beta_max_num=" << beta_max_num << " / beta_max_den=" << beta_max_den
       << "\n";
    os << "FLIP: FlippingFramework_Accelerate<mpq_class> before check\n";
    CheckFacetInequality(EXT_debug, fret, "InternalFlipFaceIneq 2");
#endif
    return fret;
  }
  Face FlipFace(Face const &sInc) {
#ifdef DEBUG_POLY_FUNDAMENTAL
    if (OneInc.count() != sInc.size()) {
      std::cerr << "Error in Flip 1\n";
      throw TerminalException{1};
    }
#endif
    MyVector<Tint> Vkernel = solver.GetKernelVector(sInc);
    return InternalFlipFaceIneq(sInc, Vkernel);
  }
  Face FlipFaceIneq(std::pair<Face, MyVector<T>> const &pair) {
    MyVector<Tint> V = NonUniqueRescaleVecRing(pair.second);
    return InternalFlipFaceIneq(pair.first, V);
  }
};

template <typename T, typename T2 = void> struct flipping_type;

template <typename T>
struct flipping_type<
    T, typename std::enable_if<has_reduction_subset_solver<T>::value>::type> {
  typedef FlippingFramework_Accelerate<T> type;
};

template <typename T>
struct flipping_type<
    T, typename std::enable_if<!has_reduction_subset_solver<T>::value>::type> {
  typedef FlippingFramework_Field<T> type;
};

template <typename T> class FlippingFramework {
public:
  typename flipping_type<T>::type flipping;
  using Text_int = typename SubsetRankOneSolver<T>::Tint;
  MyMatrix<T> const &EXT_face;
  MyMatrix<Text_int> const &EXT_face_int;
  FlippingFramework(MyMatrix<T> const &EXT, MyMatrix<Text_int> const &EXT_int,
                    Face const &OneInc, std::ostream &os)
      : flipping(EXT, EXT_int, OneInc, os), EXT_face(flipping.EXT_face),
        EXT_face_int(flipping.EXT_face_int) {}
  Face FlipFace(Face const &sInc) { return flipping.FlipFace(sInc); }
  Face FlipFaceIneq(std::pair<Face, MyVector<T>> const &pair) {
    return flipping.FlipFaceIneq(pair);
  }
};

template <typename T>
inline typename std::enable_if<
    !has_reduction_subset_solver<T>::value,
    MyMatrix<typename SubsetRankOneSolver<T>::Tint>>::type
Get_EXT_int(MyMatrix<T> const &EXT) {
  return EXT;
}

template <typename T>
inline typename std::enable_if<
    has_reduction_subset_solver<T>::value,
    MyMatrix<typename SubsetRankOneSolver<T>::Tint>>::type
Get_EXT_int(MyMatrix<T> const &EXT) {
  return UniqueRescaleRowsRing(EXT);
}

template <typename T>
Face ComputeFlipping(MyMatrix<T> const &EXT, Face const &OneInc,
                     Face const &sInc, std::ostream &os) {
  using Tint = typename SubsetRankOneSolver<T>::Tint;
  MyMatrix<T> TheEXT = ColumnReduction(EXT);
  MyMatrix<Tint> TheEXT_int = Get_EXT_int(TheEXT);
  FlippingFramework TheFram(TheEXT, TheEXT_int, OneInc, os);
  return TheFram.FlipFace(sInc);
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
#endif  // SRC_POLY_POLY_FUNDAMENTAL_H_
// clang-format on
