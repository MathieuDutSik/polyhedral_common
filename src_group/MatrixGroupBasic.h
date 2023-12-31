// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPBASIC_H_
#define SRC_GROUP_MATRIXGROUPBASIC_H_

// clang-format off
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_MATRIX_GROUP_BASIC
#endif

template <typename T, typename F>
std::optional<MyMatrix<T>> FindTransformationGeneral_f(MyMatrix<T> const &EXT1,
                                                       MyMatrix<T> const &EXT2,
                                                       F f) {
  if (EXT1.cols() != EXT2.cols())
    return {};
  if (EXT1.rows() != EXT2.rows())
    return {};
  int nbCol = EXT1.cols();
  int nbRow = EXT1.rows();
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(EXT1);
  int eRank = eSelect.TheRank;
  if (eRank != nbCol)
    return {};
  MyMatrix<T> eMat1(nbCol, nbCol);
  MyMatrix<T> eMat2(nbCol, nbCol);
  for (int iRow = 0; iRow < nbCol; iRow++) {
    int iRow1 = eSelect.ListRowSelect[iRow];
    int iRow2 = f(iRow1);
    eMat1.row(iRow) = EXT1.row(iRow1);
    eMat2.row(iRow) = EXT2.row(iRow2);
  }
  MyMatrix<T> eMat1inv = Inverse(eMat1);
  MyMatrix<T> RetMat = eMat1inv * eMat2;
  MyMatrix<T> CheckMat = EXT1 * RetMat;
  for (int iRow1 = 0; iRow1 < nbRow; iRow1++) {
    int iRow2 = f(iRow1);
    for (int iCol = 0; iCol < nbCol; iCol++)
      if (CheckMat(iRow1, iCol) != EXT2(iRow2, iCol))
        return {};
  }
  return RetMat;
}

template <typename T, typename Telt>
std::optional<MyMatrix<T>> FindTransformationGeneral(MyMatrix<T> const &EXT1,
                                                     MyMatrix<T> const &EXT2,
                                                     Telt const &ePerm) {
  auto f=[&](int iRow1) -> int {
    return ePerm.at(iRow1);
  };
  return FindTransformationGeneral_f(EXT1, EXT2, f);
}

template <typename T, typename F>
MyMatrix<T> FindTransformation_f(MyMatrix<T> const &EXT1, MyMatrix<T> const &EXT2,
                                 F f) {
  std::optional<MyMatrix<T>> opt = FindTransformationGeneral_f(EXT1, EXT2, f);
  MyMatrix<T> eMatr = unfold_opt(opt, "FindTransformationGeneral fails");
  return eMatr;
}

template <typename T, typename Telt>
MyMatrix<T> FindTransformation(MyMatrix<T> const &EXT1, MyMatrix<T> const &EXT2,
                               Telt const &ePerm) {
#ifdef DEBUG_MATRIX_GROUP_BASIC
  size_t sizPerm = ePerm.size();
  size_t sizEXT = EXT1.rows();
  if (sizPerm < sizEXT) {
    std::cerr << "sizPerm=" << sizPerm << " sizEXT=" << sizEXT << "\n";
    std::cerr << "ePerm should be of the same size as EXT1\n";
    throw TerminalException{1};
  }
#endif
  std::optional<MyMatrix<T>> opt = FindTransformationGeneral(EXT1, EXT2, ePerm);
  MyMatrix<T> eMatr = unfold_opt(opt, "FindTransformationGeneral fails");
  return eMatr;
}

template <typename T, typename Tidx>
std::optional<MyMatrix<T>> FindTransformationGeneral_vect(MyMatrix<T> const &EXT1,
                                                          MyMatrix<T> const &EXT2,
                                                          std::vector<Tidx> const& v) {
  auto f=[&](int iRow1) -> int {
    return v[iRow1];
  };
  return FindTransformationGeneral_f(EXT1, EXT2, f);
}

template <typename T, typename Tidx>
MyMatrix<T> FindTransformation_vect(MyMatrix<T> const &EXT1,
                                    MyMatrix<T> const &EXT2,
                                    std::vector<Tidx> const& v) {
  auto f=[&](int iRow1) -> int {
    return v[iRow1];
  };
  return FindTransformation_f(EXT1, EXT2, f);
}

template <typename T, typename Tgroup>
bool IsSymmetryGroupOfPolytope(MyMatrix<T> const &EXT, Tgroup const &GRP) {
  using Telt = typename Tgroup::Telt;
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  std::vector<Telt> ListGen = GRP.GeneratorsOfGroup();
  for (auto const &eGen : ListGen) {
    std::optional<MyMatrix<T>> opt =
        FindTransformationGeneral(EXTred, EXTred, eGen);
    if (!opt)
      return false;
  }
  return true;
}

// clang-format off
#endif  // SRC_GROUP_MATRIXGROUPBASIC_H_
// clang-format on
