#ifndef MATRIC_CANONICAL_FORM_H
#define MATRIC_CANONICAL_FORM_H


#include "Shvec_double.h"
#include "Temp_PolytopeEquiStab.h"
#include "MAT_MatrixInt.h"
#include "MatrixGroup.h"



template<typename T,typename Tint>
MyMatrix<T> ComputeCanonicalForm(MyMatrix<T> const& inpMat)
{
  //
  // Computing the Z-basis on which the computation relies.
  //
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T,Tint>(inpMat);
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
  //
  // Computing the canonicalization of the scalar product matrix
  //
  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector<T,T,GraphBitset>(WMat);
  std::vector<int> MapVect = PairCanonic.first;
  std::vector<int> MapVectRev = PairCanonic.second;
  //
  // Building the canonical basis
  //
  MyMatrix<Tint> SHVred(nbRow,n);
  std::vector<MyVector<T>> OldBasis_T;
  std::vector<MyVector<T>> NewBasis_T;
  MyMatrix<T> OldBasis_MatT(0,n);
  std::vector<MyVector<Tint>> OldBasis_Tint;
  std::vector<MyVector<Tint>> NewBasis_Tint;
  int TheDim=0;
  for (int iRowCan=0; iRowCan<nbRow; iRowCan++) {
    int iRowNative = MapVectRev[iRowCan];
    MyVector<Tint> eRow_Tint = GetMatrixRow(SHV, iRowNative);
    MyVector<T>    eRow_T    = ConvertVectorUniversal<T,Tint>(eRow_Tint);
    SolMatResult<T> TheSol = SolutionMat(OldBasis_MatT, eRow_T);
    if (TheSol.result) {
      MyVector<T> eVect_T = ZeroVector<T>(n);
      for (int u=0; u<TheDim; u++)
	eVect_T += TheSol.eSol(u) * OldBasis_T[u];
      if (!IsIntegralVector(eVect_T)) {
	std::cerr << "Error the vector should be integral\n";
	throw TerminalException{1};
      }
      MyVector<Tint> eVect_Tint = ConvertVectorUniversal<Tint,T>(eVect_T);
      AssignMatrixRow(SHVred, iRowNative, eVect_Tint);
    }
    else {
      OldBasis_T.push_back(eRow_T);
      OldBasis_Tint.push_back(eRow_Tint);
      OldBasis_MatT = MatrixFromVectorFamily(OldBasis_T);
      MyVector<Tint> TheVect = ZeroVector<Tint>(n);
      if (TheDim == 0) {
	FractionVector<Tint> eFr = RemoveFractionVectorPlusCoeff(eRow_Tint);
	TheVect(0) = T_abs(eFr.TheMult);
      }
      else {
	MyMatrix<Tint> OldBasis_MatTint = MatrixFromVectorFamily(OldBasis_Tint);
	MyMatrix<Tint> OrthMat = NullspaceIntTrMat(OldBasis_MatTint);
	int sizVect = OrthMat.rows();
	if (OrthMat.rows() != nbRow - TheDim) {
	  std::cerr << "Rank error, need to retry something else\n";
	  throw TerminalException{1};
	}
	MyVector<Tint> ScalVect(sizVect);
	for (int uS=0; uS<sizVect; uS++) {
	  Tint eScal=0;
	  for (int iCol=0; iCol<n; iCol++)
	    eScal += OrthMat(uS,iCol) * eRow_Tint(iCol);
	  ScalVect(uS)=eScal;
	}
	FractionVector<Tint> eFr = RemoveFractionVectorPlusCoeff(ScalVect);
	TheVect(TheDim) = T_abs(eFr.TheMult);
      }
      NewBasis_Tint.push_back(TheVect);
      AssignMatrixRow(SHVred, iRowNative, TheVect);
      TheDim++;
    }
  }
  MyMatrix<T> SHVred_T = ConvertMatrixUniversal<T,Tint>(SHVred);
  MyMatrix<T> SHV_T    = ConvertMatrixUniversal<T,Tint>(SHV);
  permlib::Permutation ePerm=IdentityPermutation(nbRow);
  MyMatrix<T> MatEquiv_T = FindTransformation(SHV_T, SHVred_T, ePerm);
  if (!IsIntegralMatrix(MatEquiv_T)) {
    std::cerr << "The Matrix MatEquiv_T should be integral\n";
    throw TerminalException{1};
  }
  MyMatrix<T> MatEquivInv_T = Inverse(MatEquiv_T);
  return MatEquivInv_T * inpMat * TransposedMat(MatEquivInv_T);
}






#endif
