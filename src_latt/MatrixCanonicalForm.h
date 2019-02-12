#ifndef MATRIC_CANONICAL_FORM_H
#define MATRIC_CANONICAL_FORM_H



#include "ShortestUniversal.h"
#include "InvariantVectorFamily.h"
#include "Temp_PolytopeEquiStab.h"
#include "MAT_MatrixInt.h"
#include "MatrixGroup.h"

//#include <time.h>


template<typename T,typename Tint>
std::pair<MyMatrix<Tint>,MyMatrix<T>> ComputeCanonicalForm(MyMatrix<T> const& inpMat)
{
  //
  // Computing the Z-basis on which the computation relies.
  //
  //  std::cerr << "inpMat=\n";
  //  WriteMatrix(std::cerr, inpMat);
  //
  //clock_t start = clock();

  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(inpMat);
  int nbRow=SHV.rows();
  int n=SHV.cols();
#ifdef DEBUG
  if (!CheckCentralSymmetry(SHV)) {
    std::cerr << "The set of vector does not respect the central symmetry condition\n";
    throw TerminalException{1};
  }
#endif

  //std::cerr << "Generate SHV: " << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
  
  //  std::cerr << "SHV=\n";
  //  WriteMatrix(std::cerr, SHV);
  //
  // Computing the scalar product matrix
  //
  T TheTol=0;
  WeightMatrix<T,T> WMat = T_TranslateToMatrix_QM_SHV(inpMat, SHV, TheTol);
  //std::cerr << "Generate inner products" << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
  ReorderingSetWeight(WMat);

  //std::cerr << "Generate weight matrix: " << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector<T,T,GraphBitset>(WMat);
  //  std::cerr << "We have PairCanonic\n";
  std::vector<int> MapVect = PairCanonic.first;
  std::vector<int> MapVectRev = PairCanonic.second;
  
  //std::cerr << "Get canonical order: " << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
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

  //std::cerr << "Done: " << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
  return {BasisCan_Tint, RetMat};
}

template<typename T,typename Tint>
std::pair<MyMatrix<Tint>,MyMatrix<Tint>> ComputeCanonicalForm(MyMatrix<Tint> const& inpMat)
{
  //
  // Computing the Z-basis on which the computation relies.
  //
  //  std::cerr << "inpMat=\n";
  //  WriteMatrix(std::cerr, inpMat);
  //
  //clock_t start = clock();

  MyMatrix<T> inpMat_T = ConvertMatrixUniversal<T,Tint>(inpMat);

  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(inpMat_T);
  int nbRow=SHV.rows();
  int n=SHV.cols();
#ifdef DEBUG
  if (!CheckCentralSymmetry(SHV)) {
    std::cerr << "The set of vector does not respect the central symmetry condition\n";
    throw TerminalException{1};
  }
#endif

  //std::cerr << "Generate SHV: " << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
  
  //  std::cerr << "SHV=\n";
  //  WriteMatrix(std::cerr, SHV);
  //
  // Computing the scalar product matrix
  //
  Tint TheTol=0;
  WeightMatrix<Tint,Tint> WMat = T_TranslateToMatrix_QM_SHV(inpMat, SHV, TheTol);
  //std::cerr << "Generate inner products" << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
  ReorderingSetWeight(WMat);

  //std::cerr << "Generate weight matrix: " << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector<Tint,Tint,GraphBitset>(WMat);
  //  std::cerr << "We have PairCanonic\n";
  std::vector<int> MapVect = PairCanonic.first;
  std::vector<int> MapVectRev = PairCanonic.second;
  
  //std::cerr << "Get canonical order: " << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
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
  MyMatrix<Tint> RetMat = BasisCan_Tint * inpMat * TransposedMat(BasisCan_Tint);

  //std::cerr << "Done: " << double(clock()-start)/CLOCKS_PER_SEC << std::endl;
  return {BasisCan_Tint, RetMat};
}




#endif
