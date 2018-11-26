#ifndef MATRIC_CANONICAL_FORM_H
#define MATRIC_CANONICAL_FORM_H



#include "ShortestUniversal.h"
#include "InvariantVectorFamily.h"
#include "Temp_PolytopeEquiStab.h"
#include "MAT_MatrixInt.h"
#include "MatrixGroup.h"



template<typename T,typename Tint>
std::pair<MyMatrix<Tint>,MyMatrix<T>> ComputeCanonicalForm(MyMatrix<T> const& inpMat)
{
  //
  // Computing the Z-basis on which the computation relies.
  //
  std::cerr << "inpMat=\n";
  WriteMatrix(std::cerr, inpMat);
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(inpMat);
  std::cerr << "SHV=\n";
  WriteMatrix(std::cerr, SHV);
  //
  // Computing the scalar product matrix
  //
  int nbRow=SHV.rows();
  int n=SHV.cols();
  MyMatrix<T> eMat(nbRow,nbRow);
  for (int iRow1=0; iRow1<nbRow; iRow1++) {
    MyVector<Tint> V1 = GetMatrixRow(SHV, iRow1);
    for (int iRow2=iRow1; iRow2<nbRow; iRow2++) {
      MyVector<Tint> V2 = GetMatrixRow(SHV, iRow2);
      T eScal = ScalarProductQuadForm(inpMat, V1, V2);
      eMat(iRow1,iRow2) = eScal;
      eMat(iRow2,iRow1) = eScal;
    }
  }
  WeightMatrix<T,T> WMat = T_TranslateToMatrixOrder(eMat);
  std::cerr << "We have WMat\n";
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector<T,T,GraphBitset>(WMat);
  std::cerr << "We have PairCanonic\n";
  std::vector<int> MapVect = PairCanonic.first;
  std::vector<int> MapVectRev = PairCanonic.second;
  //
  // Building the canonical basis
  //
  MyMatrix<T> SHVcan_T(nbRow,n);
  for (int iRowCan=0; iRowCan<nbRow; iRowCan++) {
    int iRowNative = MapVectRev[iRowCan];
    MyVector<Tint> eRow_Tint = GetMatrixRow(SHV, iRowNative);
    MyVector<T> eRow_T = ConvertMatrixUniversal<T,Tint>(eRow_Tint);
    AssignMatrixRow(SHVcan_T, iRowCan, eRow_T);
  }
  MyMatrix<T> BasisCan_T = GetZbasis(SHVcan_T);
  MyMatrix<Tint> BasisCan_Tint = ConvertMatrixUniversal<Tint,T>(BasisCan_T);
#ifdef DEBUG
  T eDet = DeterminantMat(BasisCan_T);
  T eDet_abs = T_abs(eDet);
  if (eDet_abs != 1) {
    std::cerr << "The matrix should be of determinant 1\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<T> RetMat = BasisCan_T * inpMat * TransposedMat(BasisCan_T);
  return {BasisCan_Tint, RetMat};
}





template<typename T,typename Tint>
MyMatrix<Tint> ComputeCanonicalSpanningSet(MyMatrix<T> const& inpMat)
{
  //
  // Computing the Z-basis on which the computation relies.
  //
  std::cerr << "inpMat=\n";
  WriteMatrix(std::cerr, inpMat);
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(inpMat);
  std::cerr << "SHV=\n";
  WriteMatrix(std::cerr, SHV);
  //
  // Computing the scalar product matrix
  //
  int nbRow=SHV.rows();
  int n=SHV.cols();
  MyMatrix<T> eMat(nbRow,nbRow);
  for (int iRow1=0; iRow1<nbRow; iRow1++) {
    MyVector<Tint> V1 = GetMatrixRow(SHV, iRow1);
    for (int iRow2=iRow1; iRow2<nbRow; iRow2++) {
      MyVector<Tint> V2 = GetMatrixRow(SHV, iRow2);
      T eScal = ScalarProductQuadForm(inpMat, V1, V2);
      eMat(iRow1,iRow2) = eScal;
      eMat(iRow2,iRow1) = eScal;
    }
  }
  WeightMatrix<T,T> WMat = T_TranslateToMatrixOrder(eMat);
  std::cerr << "We have WMat\n";
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector<T,T,GraphBitset>(WMat);
  std::cerr << "We have PairCanonic\n";
  std::vector<int> MapVect = PairCanonic.first;
  std::vector<int> MapVectRev = PairCanonic.second;
  //
  // Building the canonical basis
  //
  MyMatrix<Tint> SHVcan_Tint(nbRow,n);
  for (int iRowCan=0; iRowCan<nbRow; iRowCan++) {
    int iRowNative = MapVectRev[iRowCan];
    MyVector<Tint> eRow_Tint = GetMatrixRow(SHV, iRowNative);
    AssignMatrixRow(SHVcan_Tint, iRowCan, eRow_Tint);
  }
  return SHVcan_Tint;
}





template<typename T,typename Tint>
MyMatrix<T> ComputeMatrixScalarProducts(MyMatrix<T> const& inpMat)
{
  //
  // Computing the Z-basis on which the computation relies.
  //
  std::cerr << "inpMat=\n";
  WriteMatrix(std::cerr, inpMat);
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(inpMat);
  std::cerr << "SHV=\n";
  WriteMatrix(std::cerr, SHV);
  //
  // Computing the scalar product matrix
  //
  int nbRow=SHV.rows();
  MyMatrix<T> eMat(nbRow,nbRow);
  for (int iRow1=0; iRow1<nbRow; iRow1++) {
    MyVector<Tint> V1 = GetMatrixRow(SHV, iRow1);
    for (int iRow2=iRow1; iRow2<nbRow; iRow2++) {
      MyVector<Tint> V2 = GetMatrixRow(SHV, iRow2);
      T eScal = ScalarProductQuadForm(inpMat, V1, V2);
      eMat(iRow1,iRow2) = eScal;
      eMat(iRow2,iRow1) = eScal;
    }
  }
  return eMat;
}






#endif
