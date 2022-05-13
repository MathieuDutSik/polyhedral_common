#ifndef SRC_POLY_PERM_FCT_H_
#define SRC_POLY_PERM_FCT_H_

#include "Timings.h"
#include <algorithm>
#include <utility>
#include <vector>

template <typename T, typename Tidx>
std::vector<Tidx> SortingPerm(std::vector<T> const &ListV) {
  struct PairData {
    Tidx i;
    T x;
  };
  Tidx len = ListV.size();
  std::vector<PairData> ListPair(len);
  for (Tidx i = 0; i < len; i++) {
    PairData ePair{i, ListV[i]};
    ListPair[i] = ePair;
  }
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const &x1, PairData const &x2) -> bool {
         if (x1.x < x2.x)
           return true;
         if (x2.x < x1.x)
           return false;
         return x1.i < x2.i;
       });
  std::vector<Tidx> v(len);
  for (Tidx i = 0; i < len; i++)
    v[i] = ListPair[i].i;
  return v;
}

template <typename T, typename Telt>
MyMatrix<T> RepresentVertexPermutation(MyMatrix<T> const &EXT1,
                                       MyMatrix<T> const &EXT2,
                                       Telt const &ePerm) {
  std::cerr << "Beginning of RepresentVertexPermutation\n";
  std::cerr << "EXT1=\n";
  WriteMatrix(std::cerr, EXT1);
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(EXT1);
  std::vector<int> const &ListRowSelect = eSelect.ListRowSelect;
  std::cerr << "|ListRowSelect|=" << ListRowSelect.size() << " |EXT1|=" << EXT1.rows() << " / " << EXT1.cols() << "\n";
  MyMatrix<T> M1 = SelectRow(EXT1, ListRowSelect);
  MyMatrix<T> M1inv = Inverse(M1);
  size_t nbRow_s = ListRowSelect.size();
  std::vector<int> ListRowSelectImg(nbRow_s);
  for (size_t iRow = 0; iRow < nbRow_s; iRow++)
    ListRowSelectImg[iRow] = ePerm.at(iRow);
  MyMatrix<T> M2 = SelectRow(EXT2, ListRowSelectImg);
  MyMatrix<T> RetMat = M1inv * M2;
#ifdef SANITY_CHECK
  std::cerr << "Doing sanity_checks in RepresentVertexPermutation\n";
  int nbRow=EXT2.rows();
  int nbCol=EXT2.cols();
  MyMatrix<T> EXT1_img = EXT1 * RetMat;
  MyMatrix<T> EXT2_perm(nbRow, nbCol);
  std::cerr << "EXT1_img=\n";
  WriteMatrix(std::cerr, EXT1_img);
  for (int iRow=0; iRow<nbRow; iRow++) {
    int iRowImg = ePerm.at(iRow);
    for (int iCol=0; iCol<nbCol; iCol++) {
      EXT2_perm(iRow,iCol) = EXT2(iRowImg,iCol);
    }
  }
  if (!TestEqualityMatrix(EXT1_img, EXT2_perm)) {
    std::cerr << "The matrices are not equal\n";
    std::cerr << "EXT1_img=\n";
    WriteMatrix(std::cerr, EXT1_img);
    std::cerr << "EXT2_perm=\n";
    WriteMatrix(std::cerr, EXT2_perm);
    throw TerminalException{1};
  }
  if (RankMat(RetMat) != nbCol) {
    std::cerr << "RetMat should be of full rank\n";
    std::cerr << "RetMat=\n";
    WriteMatrix(std::cerr, RetMat);
    throw TerminalException{1};
  }
#endif
  return RetMat;
}

template<typename T, typename Tidx>
std::optional<std::vector<Tidx>> FindPermutationalEquivalence(MyMatrix<T> const& M1, MyMatrix<T> const& M2)
{
  int n_rows = M1.rows();
  if (n_rows != M2.rows()) {
    return {};
  }
  Tidx n_rows_tidx = n_rows;
  std::unordered_map<MyVector<T>,Tidx> map;
  for (Tidx i=0; i<n_rows_tidx; i++) {
    MyVector<T> V = GetMatrixRow(M1, i);
    map[V] = i + 1;
  }
  std::vector<Tidx> eList(n_rows);
  std::vector<Tidx> LMatch(n_rows, 0);
  for (Tidx i=0; i<n_rows_tidx; i++) {
    MyVector<T> V = GetMatrixRow(M2, i);
    Tidx const& val = map[V];
    if (val == 0)
      return {};
    eList[i] = val - 1;
    LMatch[val-1]++;
  }
  for (Tidx i=0; i<n_rows_tidx; i++) {
    if (LMatch[i] != 1) {
      return {};
    }
  }
  return eList;
}


template <typename T, typename Tfield, typename Tidx, typename F1, typename F2>
std::optional<MyMatrix<Tfield>>
FindMatrixTransformationTest_Generic(size_t nbRow, size_t nbCol, F1 f1, F2 f2,
                             std::vector<Tidx> const &eList) {
  static_assert(is_ring_field<Tfield>::value,
                "Requires Tfield to be a field in DivideVector");
  auto f = [&](MyMatrix<Tfield> &M, size_t eRank, size_t iRow) -> void {
    const MyVector<T> &V = f1(iRow);
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      M(eRank, iCol) = UniversalScalarConversion<Tfield, T>(V(iCol));
    }
  };
  SelectionRowCol<Tfield> eSelect =
      TMat_SelectRowCol_Kernel<Tfield>(nbRow, nbCol, f);
  if (eSelect.TheRank != nbCol) {
#ifdef DEBUG_PERM_FCT
    std::cerr << "FindMatrixTransformationTest_Generic, exit false 1\n";
#endif
    return {};
  }
#ifdef DEBUG_PERM_FCT
  std::cerr << "ListRowSelect =";
  for (auto &eVal : eSelect.ListRowSelect)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  std::cerr << "Img(ListRowSelect) =";
  for (auto &eVal : eSelect.ListRowSelect)
    std::cerr << " " << eList[eVal];
  std::cerr << "\n";
#endif
  MyMatrix<Tfield> M1_field(eSelect.TheRank, nbCol);
  for (size_t iRow = 0; iRow < eSelect.TheRank; iRow++) {
    size_t jRow = eSelect.ListRowSelect[iRow];
    const MyVector<T> &V = f1(jRow);
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      M1_field(iRow, iCol) = UniversalScalarConversion<Tfield, T>(V(iCol));
  }
  MyMatrix<Tfield> M1inv_field = Inverse(M1_field);
  MyMatrix<Tfield> M2_field(eSelect.TheRank, nbCol);
  for (size_t iRow = 0; iRow < eSelect.TheRank; iRow++) {
    size_t jRow = eList[eSelect.ListRowSelect[iRow]];
    const MyVector<T> &V = f2(jRow);
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      M2_field(iRow, iCol) = UniversalScalarConversion<Tfield, T>(V(iCol));
  }
  MyMatrix<Tfield> EqMat = M1inv_field * M2_field;
#ifdef DEBUG_PERM_FCT
  // Coef EqMat(0,1) = sum_j M1inv_field(0,j) * M2_field(j,1)
  Tfield eSP = 0;
  for (int j = 0; j < int(nbCol); j++) {
    eSP += M1inv_field(0, j) * M2_field(j, 1);
  }
  std::cerr << "eSP=" << eSP << "\n";
  std::cerr << "M1_field=\n";
  WriteMatrix(std::cerr, M1_field);
  std::cerr << "M1inv_field=\n";
  WriteMatrix(std::cerr, M1inv_field);
  std::cerr << "M2_field=\n";
  WriteMatrix(std::cerr, M2_field);
  std::cerr << "EqMat=\n";
  WriteMatrix(std::cerr, EqMat);
#endif
  // Now testing that we have EXT1 EqMat = EXT2
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    size_t iRowImg = eList[iRow];
#ifdef DEBUG_PERM_FCT
    std::cerr << "iRow=" << iRow << " iRowImg=" << iRowImg << "\n";
#endif
    // We can have f1 = f2 which zould invalidate reference so copy is needed
    MyVector<T> V1 = f1(iRow);
    const MyVector<T> &V2 = f2(iRowImg);
#ifdef DEBUG_PERM_FCT_NOW_OK
    std::cerr << "V1      =";
    WriteVector(std::cerr, V1);
    std::cerr << "V2      =";
    WriteVector(std::cerr, V2);
    //
    MyVector<Tfield> V1_T = UniversalVectorConversion<Tfield, T>(V1);
    std::cerr << "V1_T    =";
    WriteVector(std::cerr, V1_T);
    //
    MyVector<Tfield> V2_T = UniversalVectorConversion<Tfield, T>(V2);
    std::cerr << "V2_T    =";
    WriteVector(std::cerr, V2_T);
    //
    MyVector<Tfield> V1_img = EqMat.transpose() * V1_T;
    std::cerr << "V1_img  =";
    WriteVector(std::cerr, V1_img);
    //
    MyVector<Tfield> V1_expr = M1inv_field.transpose() * V1_T;
    std::cerr << "V1_expr =";
    WriteVector(std::cerr, V1_expr);
    //
    MyVector<Tfield> V1_imgB = M2_field.transpose() * V1_expr;
    std::cerr << "V1_imgB =";
    WriteVector(std::cerr, V1_imgB);
    //
    MyVector<Tfield> V1_imgC =
        M2_field.transpose() * M1inv_field.transpose() * V1_T;
    std::cerr << "V1_imgC =";
    WriteVector(std::cerr, V1_imgC);
    //
    MyMatrix<Tfield> eProd = M1inv_field * M2_field;
    MyVector<Tfield> V1_imgD = eProd.transpose() * V1_T;
    std::cerr << "V1_imgD =";
    WriteVector(std::cerr, V1_imgD);
#endif
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      Tfield eSum = -V2(iCol);
      for (size_t jRow = 0; jRow < nbCol; jRow++)
        eSum += EqMat(jRow, iCol) * V1(jRow);
      if (eSum != 0) {
#ifdef DEBUG_PERM_FCT
        std::cerr << "eSum=" << eSum << " iRow=" << iRow << " iCol=" << iCol
                  << "\n";
        std::cerr << "FindMatrixTransformationTest_Generic, exit false 2\n";
#endif
        return {};
      }
    }
  }
#ifdef DEBUG_PERM_FCT
  std::cerr << "FindMatrixTransformationTest_Generic, exit true\n";
#endif
  return EqMat;
}

template <typename T, typename Tidx>
std::optional<MyMatrix<T>>
FindMatrixTransformationTest(MyMatrix<T> const &EXT1, MyMatrix<T> const &EXT2,
                             std::vector<Tidx> const &eList) {
  size_t nbRow1 = EXT1.rows();
  size_t nbRow2 = EXT2.rows();
  if (nbRow1 != nbRow2)
    return {};
  int nbCol = EXT1.cols();
  MyVector<T> V1(nbCol);
  MyVector<T> V2(nbCol);
  auto f1 = [&](size_t iRow) -> const MyVector<T> & {
    for (int i = 0; i < nbCol; i++)
      V1(i) = EXT1(iRow, i);
    return V1;
  };
  auto f2 = [&](size_t iRow) -> const MyVector<T> & {
    for (int i = 0; i < nbCol; i++)
      V2(i) = EXT2(iRow, i);
    return V2;
  };
  return FindMatrixTransformationTest_Generic<T,T,Tidx>(nbRow1, nbCol, f1, f2, eList);
}

template <typename T, typename Tfield, typename Tidx>
std::optional<MyMatrix<Tfield>>
FindMatrixTransformationTest_Subset(const MyMatrix<T> &EXT,
                                    const std::vector<Tidx> &Vsubset,
                                    const std::vector<Tidx> &Vin) {
#ifdef TIMINGS
  SingletonTime time1;
#endif
  size_t nbCol = EXT.cols();
#ifdef DEBUG_PERM_FCT
  std::cerr << "Vsubset =";
  for (auto &eVal : Vsubset)
    std::cerr << " " << int(eVal);
  std::cerr << "\n";
  std::cerr << "Vin =";
  for (auto &eVal : Vin)
    std::cerr << " " << int(eVal);
  std::cerr << "\n";
#endif
  MyVector<T> V(nbCol);
  auto g1 = [&](size_t iRow) -> const MyVector<T> & {
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      V(iCol) = EXT(Vsubset[iRow], iCol);
    return V;
  };
  std::optional<MyMatrix<Tfield>> test1 =
      FindMatrixTransformationTest_Generic<T, Tfield, Tidx>(Vsubset.size(), nbCol, g1,
                                                    g1, Vin);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|FindMatrixTransformationTest_Subset|=" << ms(time1, time2)
            << "\n";
#endif
  return test1;
}

template <typename T, typename Tfield, typename Tidx>
bool IsSubsetFullRank(const MyMatrix<T> &EXT,
                      const std::vector<Tidx> &Vsubset) {
  size_t nbCol = EXT.cols();
  if (Vsubset.size() < nbCol)
    return false;
#ifdef TIMINGS
  SingletonTime time1;
#endif
  auto f = [&](MyMatrix<Tfield> &M, size_t eRank, size_t iRow) -> void {
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      M(eRank, iCol) =
          UniversalScalarConversion<Tfield, T>(EXT(Vsubset[iRow], iCol));
  };
  SelectionRowCol<Tfield> TheSol =
      TMat_SelectRowCol_Kernel<Tfield>(Vsubset.size(), nbCol, f);
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|IsSubsetFullRank|=" << ms(time1, time2) << "\n";
#endif
  return TheSol.TheRank == nbCol;
}

template <typename T, typename Tfield, typename Tidx>
std::optional<std::vector<Tidx>>
RepresentVertexPermutationTest(MyMatrix<T> const &EXT1, MyMatrix<T> const &EXT2,
                               MyMatrix<Tfield> const &P) {
#ifdef TIMINGS
  SingletonTime time1;
#endif
  size_t n_rows = EXT1.rows();
  size_t n_cols = EXT1.cols();
  MyMatrix<T> VectorContain(1, n_cols);
  ContainerMatrix<T> Cont(EXT2, VectorContain);
  //
  // We are testing if EXT1 P = perm(EXT2)
  std::vector<Tidx> V(n_rows);
  Face f(n_rows);
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    for (size_t i_col = 0; i_col < n_cols; i_col++) {
      Tfield eSum1 = 0;
      for (size_t j_row = 0; j_row < n_cols; j_row++)
        eSum1 += EXT1(i_row, j_row) * P(j_row, i_col);
      std::pair<bool, T> rec_eSum2 =
          UniversalScalarConversionCheck<T, Tfield>(eSum1);
      if (!rec_eSum2.first) {
#ifdef TIMINGS
        SingletonTime time2;
        std::cerr << "ESC1 |RepresentVertexPermutationTest|="
                  << ms(time1, time2) << "\n";
#endif
        return {}; // We fail because the image is not integral.
      }
      VectorContain(0, i_col) = rec_eSum2.second;
    }
    std::pair<bool, size_t> epair = Cont.GetIdx();
    if (!epair.first) {
#ifdef TIMINGS
      SingletonTime time2;
      std::cerr << "ESC2 |RepresentVertexPermutationTest|=" << ms(time1, time2)
                << "\n";
#endif
      return {}; // We fail because the image does not belong to EXT2
    }
    V[i_row] = epair.second;
    f[epair.second] = 1;
  }
#ifdef SANITY_CHECK
  int n_error = 0;
  for (size_t i_row = 0; i_row < n_rows; i_row++)
    if (f[i_row] == 0)
      n_error++;
  if (n_error > 0) {
    std::cerr << "We found several errors. n_error=" << n_error << "\n";
    throw TerminalException{1};
  }
#endif
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|RepresentVertexPermutationTest|=" << ms(time1, time2) << "\n";
#endif
  return V;
}

template <typename Tidx> struct DataMapping {
  bool correct;
  Face block_status;
  std::vector<Tidx> eGen;
};

template <typename T, typename Tfield, typename Tidx>
DataMapping<Tidx> RepresentVertexPermutationTest_Blocks(
    const MyMatrix<T> &EXT, const MyMatrix<Tfield> &P,
    const std::vector<Tidx> &Vsubset, const std::vector<Tidx> &Vin,
    const std::vector<std::vector<Tidx>> &ListBlocks) {
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
  std::vector<Tidx> eGen(n_rows);
  bool correct = true;
  size_t n_block = ListBlocks.size();
  std::vector<Tidx> MapRev(n_rows);
  auto f_insert = [&](const std::vector<Tidx> &eSubset,
                      const std::vector<Tidx> &eMap) -> void {
    size_t len = eSubset.size();
    for (size_t i = 0; i < len; i++) {
      Tidx pos = eMap[i];
      Tidx i_m = eSubset[i];
      Tidx pos_m = eSubset[pos];
      eGen[i_m] = pos_m;
    }
  };
  f_insert(Vsubset, Vin);
  //
  Face block_status(n_block);
  for (size_t i_block = 0; i_block < n_block; i_block++) {
    size_t blk_size = ListBlocks[i_block].size();
    MyMatrix<T> EXTred(blk_size, n_cols);
    for (size_t i_row = 0; i_row < blk_size; i_row++) {
      size_t pos = ListBlocks[i_block][i_row];
      EXTred.row(i_row) = EXT.row(pos);
    }
    std::optional<std::vector<Tidx>> eEquiv =
        RepresentVertexPermutationTest<T, Tfield, Tidx>(EXTred, EXTred, P);
    if (!eEquiv) {
      correct = false;
    } else {
      f_insert(ListBlocks[i_block], *eEquiv);
    }
    block_status[i_block] = bool(eEquiv);
  }
  return {correct, std::move(block_status), std::move(eGen)};
}

template <typename T, typename Tfield, typename Tidx>
std::vector<Tidx>
ExtendPartialCanonicalization(const MyMatrix<T> &EXT,
                              const std::vector<Tidx> &Vsubset,
                              const std::vector<Tidx> &PartOrd) {
#ifdef TIMINGS
  SingletonTime time1;
#endif
  size_t nbRow = EXT.rows();
  size_t nbCol = EXT.cols();
  auto f = [&](MyMatrix<Tfield> &M, size_t eRank, size_t iRow) -> void {
    size_t pos = Vsubset[PartOrd[iRow]];
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      M(eRank, iCol) = UniversalScalarConversion<Tfield, T>(EXT(pos, iCol));
  };
  SelectionRowCol<Tfield> TheSol =
      TMat_SelectRowCol_Kernel<Tfield>(Vsubset.size(), nbCol, f);
  // Selecting the submatrix
  MyMatrix<Tfield> M(nbCol, nbCol);
  for (size_t iRow = 0; iRow < nbCol; iRow++) {
    size_t pos = Vsubset[PartOrd[TheSol.ListRowSelect[iRow]]];
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      M(iRow, iCol) = UniversalScalarConversion<Tfield, T>(EXT(pos, iCol));
  }
  MyMatrix<Tfield> Minv0 = Inverse(M);
  MyMatrix<Tfield> Minv1 = RemoveFractionMatrix(Minv0);
  MyMatrix<T> Minv2 = UniversalMatrixConversion<T, Tfield>(Minv1);
  MyMatrix<T> Mop = EXT * Minv2;
  std::vector<Tidx> ListIdx(nbRow);
  for (size_t iRow = 0; iRow < nbRow; iRow++)
    ListIdx[iRow] = iRow;
  std::sort(ListIdx.begin(), ListIdx.end(),
            [&](size_t iRow, size_t jRow) -> bool {
              for (size_t iCol = 0; iCol < nbCol; iCol++) {
                if (Mop(iRow, iCol) < Mop(jRow, iCol))
                  return true;
                if (Mop(iRow, iCol) > Mop(jRow, iCol))
                  return false;
              }
              return false;
            });
#ifdef TIMINGS
  SingletonTime time2;
  std::cerr << "|ExtendPartialCanonicalization|=" << ms(time1, time2) << "\n";
#endif
  return ListIdx;
}

template <typename T, typename Tfield, typename Tidx>
bool CheckEquivalence(const MyMatrix<T> &EXT1, const MyMatrix<T> &EXT2,
                      const std::vector<Tidx> &ListIdx,
                      const MyMatrix<Tfield> &P) {
  size_t n_rows = EXT1.rows();
  size_t n_cols = EXT1.cols();
  //
  // We are testing if EXT1 P = perm(EXT2)
  //  std::vector<Tidx> V(n_rows);
  /*
  std::cerr << "CheckEquivalence EXT1=\n";
  WriteMatrix(std::cerr, EXT1);
  std::cerr << "CheckEquivalence EXT2=\n";
  WriteMatrix(std::cerr, EXT2);
  std::cerr << "ListIdx=" << ListIdx << "\n";
  std::cerr << "CheckEquivalence P=\n";
  WriteMatrix(std::cerr, P);
  */

  MyVector<T> Vimg(n_cols);
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    size_t i_row_img = ListIdx[i_row];
    for (size_t i_col = 0; i_col < n_cols; i_col++) {
      Tfield eSum1 = 0;
      for (size_t j_row = 0; j_row < n_cols; j_row++)
        eSum1 += EXT1(i_row, j_row) * P(j_row, i_col);
      std::pair<bool, T> rec_eSum2 =
          UniversalScalarConversionCheck<T, Tfield>(eSum1);
      if (!rec_eSum2.first) {
        return false; // We fail because the image is not integral.
      }
      T Img_EXT1 = rec_eSum2.second;
      Vimg[i_col] = Img_EXT1;
      //
      T EXT2_map = EXT2(i_row_img, i_col);
      if (Img_EXT1 != EXT2_map)
        return false;
    }
    //    std::cerr << "i_row=" << i_row << " Vimg=" << StringVectorGAP(Vimg) <<
    //    " i_row_img=" << i_row_img << " EXTé(i_row_img)=" <<
    //    StringVectorGAP(GetMatrixRow(EXT2, i_row_img)) << "\n";
  }
  return true;
}

template <typename T, typename Tfield, typename Tidx>
MyMatrix<Tfield> GetBasisFromOrdering(const MyMatrix<T> &EXT,
                                      const std::vector<Tidx> &Vsubset) {
  size_t nbCol = EXT.cols();
  auto f = [&](MyMatrix<Tfield> &M, size_t eRank, size_t iRow) -> void {
    size_t pos = Vsubset[iRow];
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      M(eRank, iCol) = UniversalScalarConversion<Tfield, T>(EXT(pos, iCol));
  };
  SelectionRowCol<Tfield> TheSol =
      TMat_SelectRowCol_Kernel<Tfield>(Vsubset.size(), nbCol, f);
  // Selecting the submatrix
  MyMatrix<Tfield> M(nbCol, nbCol);
  for (size_t iRow = 0; iRow < nbCol; iRow++) {
    size_t pos = Vsubset[TheSol.ListRowSelect[iRow]];
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      M(iRow, iCol) = UniversalScalarConversion<Tfield, T>(EXT(pos, iCol));
  }
  return M;
}

template <typename T, typename Telt>
Telt GetPermutationOnVectors(MyMatrix<T> const &EXT1, MyMatrix<T> const &EXT2) {
  using Tidx = typename Telt::Tidx;
  size_t nbVect = EXT1.rows();
  std::vector<MyVector<T>> EXTrow1(nbVect), EXTrow2(nbVect);
  for (size_t iVect = 0; iVect < nbVect; iVect++) {
    EXTrow1[iVect] = GetMatrixRow(EXT1, iVect);
    EXTrow2[iVect] = GetMatrixRow(EXT2, iVect);
  }
  Telt ePerm1 = Telt(SortingPerm<MyVector<T>, Tidx>(EXTrow1));
  Telt ePerm2 = Telt(SortingPerm<MyVector<T>, Tidx>(EXTrow2));
  Telt ePermRet = (~ePerm1) * ePerm2;
#if defined DEBUG_PERM_FCT || defined SANITY_CHECK
  for (size_t iVect = 0; iVect < nbVect; iVect++) {
    size_t jVect = ePermRet.at(iVect);
    if (EXTrow2[jVect] != EXTrow1[iVect]) {
      std::cerr << "iVect=" << iVect << " jVect=" << jVect << "\n";
      std::cerr << "Id:";
      for (size_t k = 0; k < nbVect; k++)
        std::cerr << " " << k;
      std::cerr << "\n";
      std::cerr << " p:";
      for (size_t k = 0; k < nbVect; k++)
        std::cerr << " " << ePermRet.at(k);
      std::cerr << "\n";

      std::cerr << "perm1:";
      for (size_t k = 0; k < nbVect; k++)
        std::cerr << " " << ePerm1.at(k);
      std::cerr << "\n";
      std::cerr << "perm2:";
      for (size_t k = 0; k < nbVect; k++)
        std::cerr << " " << ePerm2.at(k);
      std::cerr << "\n";

      std::cerr << "EXTrow1[iVect]=";
      WriteVector(std::cerr, EXTrow1[iVect]);
      std::cerr << "EXTrow2[jVect]=";
      WriteVector(std::cerr, EXTrow2[jVect]);
      std::cerr << "EXTrow1=\n";
      WriteMatrix(std::cerr, EXT1);
      std::cerr << "EXTrow2=\n";
      WriteMatrix(std::cerr, EXT2);
      std::cerr << "Error in GetPermutationOnVectors\n";
      throw TerminalException{1};
    }
  }
#endif
  return ePermRet;
}

template <typename Tgroup>
Tgroup GetGroupListGen(std::vector<std::vector<unsigned int>> const &ListGen,
                       size_t const &nbVert) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> generatorList;
  std::vector<Tidx> gList(nbVert);
  for (auto &eGen : ListGen) {
    for (size_t iVert = 0; iVert < nbVert; iVert++)
      gList[iVert] = eGen[iVert];
    generatorList.push_back(Telt(gList));
  }
  return Tgroup(generatorList, nbVert);
}

template <typename Tgr>
bool CheckListGenerators(std::vector<std::vector<unsigned int>> const &ListGen,
                         Tgr const &eGR) {
  size_t nbVert = eGR.GetNbVert();
  for (auto &eGen : ListGen) {
    for (size_t iVert = 0; iVert < nbVert; iVert++) {
      int eColor = eGR.GetColor(iVert);
      int iVert_img = eGen[iVert];
      int fColor = eGR.GetColor(iVert_img);
      if (eColor != fColor) {
        return false;
      }
      //
      for (auto &jVert : eGR.Adjacency(iVert)) {
        int jVert_img = eGen[jVert];
        bool test = eGR.IsAdjacent(iVert_img, jVert_img);
        if (!test)
          return false;
      }
    }
  }
  return true;
}

// clang-format off
#endif  // SRC_POLY_PERM_FCT_H_
// clang-format on
