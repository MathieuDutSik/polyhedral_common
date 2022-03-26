#ifndef INCLUDE_MATRIX_GROUP_BASIC_H
#define INCLUDE_MATRIX_GROUP_BASIC_H

template <typename T, typename Telt>
std::optional<MyMatrix<T>> FindTransformationGeneral(MyMatrix<T> const &EXT1,
                                                     MyMatrix<T> const &EXT2,
                                                     Telt const &ePerm) {
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
    int iRow2 = OnPoints(iRow1, ePerm);
    eMat1.row(iRow) = EXT1.row(iRow1);
    eMat2.row(iRow) = EXT2.row(iRow2);
  }
  MyMatrix<T> eMat1inv = Inverse(eMat1);
  MyMatrix<T> RetMat = eMat1inv * eMat2;
  MyMatrix<T> CheckMat = EXT1 * RetMat;
  for (int iRow1 = 0; iRow1 < nbRow; iRow1++) {
    int iRow2 = ePerm.at(iRow1);
    for (int iCol = 0; iCol < nbCol; iCol++)
      if (CheckMat(iRow1, iCol) != EXT2(iRow2, iCol))
        return {};
  }
  return RetMat;
}

template <typename T, typename Telt>
MyMatrix<T> FindTransformation(MyMatrix<T> const &EXT1, MyMatrix<T> const &EXT2,
                               Telt const &ePerm) {
  size_t siz1 = ePerm.size();
  size_t siz2 = EXT1.rows();
  if (siz1 < siz2) {
    std::cerr << "siz1=" << siz1 << " siz2=" << siz2 << "\n";
    std::cerr << "ePerm should be of the same size as EXT1\n";
    throw TerminalException{1};
  }
  std::optional<MyMatrix<T>> opt = FindTransformationGeneral(EXT1, EXT2, ePerm);
  if (!opt) {
    std::cerr << "FindTransformationGeneral fails\n";
    throw TerminalException{1};
  }
  return *opt;
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

#endif
