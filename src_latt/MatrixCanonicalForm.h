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
  //  std::cerr << "inpMat=\n";
  //  WriteMatrix(std::cerr, inpMat);
#ifdef DEBUG_TIME
  std::cerr << "Begining of ComputeCanonicalForm\n";
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(inpMat);
  std::cerr << "|SHV|=" << SHV.rows() << "\n";
#ifdef DEBUG_TIME
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "ExtractInvariantVectorFamilyZbasis : time2 - time1=" << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << "\n";
#endif
  int nbRow=SHV.rows();
  int n=SHV.cols();
#ifdef DEBUG
  if (!CheckCentralSymmetry(SHV)) {
    std::cerr << "The set of vector does not respect the central symmetry condition\n";
    throw TerminalException{1};
  }
#endif
  //  std::cerr << "SHV=\n";
  //  WriteMatrix(std::cerr, SHV);
  //
  // Computing the scalar product matrix
  //
  T TheTol=0;
  WeightMatrix<T,T> WMat = T_TranslateToMatrix_QM_SHV(inpMat, SHV, TheTol);
#ifdef DEBUG_TIME
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "WMat : time3 - time2=" << std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count() << "\n";
#endif
  ReorderingSetWeight(WMat);
#ifdef DEBUG_TIME
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "ReorderingSetWeight : time4 - time3=" << std::chrono::duration_cast<std::chrono::milliseconds>(time4 - time3).count() << "\n";
#endif
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector<T,T,GraphBitset>(WMat);
#ifdef DEBUG_TIME
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
  std::cerr << "GetCanonicalizationVector : time5 - time4=" << std::chrono::duration_cast<std::chrono::milliseconds>(time5 - time4).count() << "\n";
#endif
  //  std::cerr << "We have PairCanonic\n";
  std::vector<int> MapVect = PairCanonic.first;
  std::vector<int> MapVectRev = PairCanonic.second;
  //
  // Building the canonical basis
  //
  MyMatrix<Tint> SHVcan_Tint(n,nbRow);
  for (int iRowCan=0; iRowCan<nbRow; iRowCan++) {
    int iRowNative = MapVectRev[iRowCan];
    MyVector<Tint> eRow_Tint = GetMatrixRow(SHV, iRowNative);
    AssignMatrixCol(SHVcan_Tint, iRowCan, eRow_Tint);
  }
#ifdef DEBUG_TIME
  std::chrono::time_point<std::chrono::system_clock> time6 = std::chrono::system_clock::now();
  std::cerr << "SHVcan : time6 - time5=" << std::chrono::duration_cast<std::chrono::milliseconds>(time6 - time5).count() << "\n";
#endif
  MyMatrix<Tint> BasisCan_Tint_pre = ComputeRowHermiteNormalForm(SHVcan_Tint).first;
  MyMatrix<Tint> BasisCan_Tint = TransposedMat(Inverse(BasisCan_Tint_pre));
  std::cerr << "SHVcan_Tint=\n";
  WriteMatrix(std::cerr, SHVcan_Tint);
#ifdef DEBUG_TIME
  std::chrono::time_point<std::chrono::system_clock> time7 = std::chrono::system_clock::now();
  std::cerr << "ReductionMatrix : time7 - time6=" << std::chrono::duration_cast<std::chrono::milliseconds>(time7 - time6).count() << "\n";
#endif
  MyMatrix<T> BasisCan_T = ConvertMatrixUniversal<T,Tint>(BasisCan_Tint);
#ifdef DEBUG
  T eDet = DeterminantMat(BasisCan_T);
  T eDet_abs = T_abs(eDet);
  if (eDet_abs != 1) {
    std::cerr << "The matrix should be of determinant 1\n";
    throw TerminalException{1};
  }
#endif
  std::cerr << "BasisCan_T=\n";
  WriteMatrix(std::cerr, BasisCan_T);
  std::cerr << "inpMat=\n";
  WriteMatrix(std::cerr, inpMat);


  MyMatrix<T> RetMat = BasisCan_T * inpMat * TransposedMat(BasisCan_T);
  std::cerr << "RetMat=\n";
  WriteMatrix(std::cerr, RetMat);
#ifdef DEBUG_TIME
  std::chrono::time_point<std::chrono::system_clock> time8 = std::chrono::system_clock::now();
  std::cerr << "Matrix products : time8 - time7=" << std::chrono::duration_cast<std::chrono::milliseconds>(time8 - time7).count() << "\n";
#endif
  return {BasisCan_Tint, RetMat};
}


#endif
