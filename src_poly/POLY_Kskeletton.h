#ifndef TEMP_SKELETTON
#define TEMP_SKELETTON

#include "MAT_Matrix.h"
#include "POLY_PolytopeFct.h"
#include "POLY_LinearProgramming.h"
#include "Temp_PolytopeEquiStab.h"
#include "GRP_GroupFct.h"

// We follow here the conventions of SPAN_face_LinearProgramming
// in Kskeleton.g for the computation.
template<typename T, typename Tgroup>
vectface SPAN_face_LinearProgramming(Face const& face, Tgroup const& StabFace, MyMatrix<T> const& FAC, Tgroup const& FullGRP)
{
  MyMatrix<T> FACred=ColumnReduction(FAC);
  int nbRow=FACred.rows();
  vectface TheReturn(nbRow);
  int nbCol=FACred.cols();
  MyMatrix<T> eMatId=IdentityMat<T>(nbCol);
  std::vector<int> Treated(nbRow,0);
  int sizFace=face.count();
  MyMatrix<T> TestMat(sizFace+1, nbCol);
  int ePt=face.find_first();
  for (int iRow=0; iRow<sizFace; iRow++) {
    Treated[ePt]=1;
    MyVector<T> eVect=FACred.row(ePt);
    AssignMatrixRow(TestMat, iRow, eVect);
    ePt=face.find_next(ePt);
  }
  for (int iRow=0; iRow<nbRow; iRow++)
    if (Treated[iRow] == 0) {
      MyVector<T> eVect=FACred.row(iRow);
      AssignMatrixRow(TestMat, sizFace, eVect);
      SelectionRowCol<T> eSelect=TMat_SelectRowCol(TestMat);
      MyMatrix<T> NSP=eSelect.NSP;
      int nbEqua=NSP.rows();
      Face eCand(nbRow);
      Face eCandCompl(nbRow);
      std::vector<int> gList(nbRow);
      for (int jRow=0; jRow<nbRow; jRow++) {
	int test=1;
	for (int iEqua=0; iEqua<nbEqua; iEqua++)
	  if (test == 1) {
	    T eSum=0;
	    for (int iCol=0; iCol<nbCol; iCol++)
	      eSum += FACred(jRow, iCol)*NSP(iEqua, iCol);
	    if (eSum != 0)
	      test=0;
	  }
	gList[jRow]=test;
	if (test == 1)
	  eCand[jRow]=1;
	else
	  eCandCompl[jRow]=1;
      }
      std::vector<int> rList = OrbitUnion(StabFace, gList);
      for (int jRow=0; jRow<nbRow; jRow++)
	if (rList[jRow] == 1)
	  Treated[jRow]=1;
      Tgroup TheStab=FullGRP.Stabilizer_OnSets(eCand);
      vectface ListOrb=DecomposeOrbitPoint(TheStab, eCand);
      int nbOrb=ListOrb.size();
      MyMatrix<T> ListVectSpann(nbOrb, nbCol);
      for (int iOrb=0; iOrb<nbOrb; iOrb++) {
	Face eOrb=ListOrb[iOrb];
	MyVector<T> eVec=SumMatrixLineSubset(FACred, eOrb);
	AssignMatrixRow(ListVectSpann, iOrb, eVec);
      }
      MyMatrix<T> BasisSpann=RowReduction(ListVectSpann);
      int nbRowSpann=BasisSpann.rows();
      int LPdim=nbCol - nbRowSpann;
      MyMatrix<T> PreTheTot=Concatenate(BasisSpann, eMatId);
      MyMatrix<T> TheTot=RowReduction(PreTheTot);
      // the complement
      vectface ListOrbCompl=DecomposeOrbitPoint(TheStab, eCandCompl);
      int nbOrbCompl=ListOrbCompl.size();
      MyMatrix<T> PreListVectors(nbOrbCompl, nbCol);
      for (int iOrb=0; iOrb<nbOrbCompl; iOrb++) {
	Face eOrb=ListOrbCompl[iOrb];
	MyVector<T> eVec=SumMatrixLineSubset(FACred, eOrb);
	AssignMatrixRow(PreListVectors, iOrb, eVec);
      }
      MyMatrix<T> TheTotInv=Inverse(TheTot);
      MyMatrix<T> PreListVectorsB=PreListVectors*TheTotInv;
      MyMatrix<T> ListVectors(nbOrbCompl, LPdim);
      for (int jRow=0; jRow<nbOrbCompl; jRow++)
	for (int iCol=0; iCol<LPdim; iCol++) {
	  int jCol=iCol + nbRowSpann;
	  T eVal=PreListVectorsB(jRow, jCol);
	  ListVectors(jRow, iCol)=eVal;
	}
      bool eTestExist=TestPositiveRelationSimple(ListVectors);
      if (eTestExist == 0)
	TheReturn.push_back(eCand);
    }
  return TheReturn;
}




// We test if eSet is included in a proper face of the polytope
template<typename T>
bool TestInclusionProperFace(std::vector<int> const& eSet, MyMatrix<T> const& FAC)
{
  int nbCol=FAC.cols();
  int nbRow=FAC.rows();
  //  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  //  std::cerr << "eSet=";
  //  PrintVectorInt_GAP(std::cerr, eSet);
  //  std::cerr << "\n";
  std::vector<int> eVectCand(nbRow, 0);
  for (auto eVal : eSet)
    eVectCand[eVal]=1;
  MyMatrix<T> eMatId=IdentityMat<T>(nbCol);
  while(true) {
    int len=0;
    for (int iRow=0; iRow<nbRow; iRow++)
      if (eVectCand[iRow] == 1)
	len++;
    MyMatrix<T> TestMat(len, nbCol);
    int jRow=0;
    for (int iRow=0; iRow<nbRow; iRow++)
      if (eVectCand[iRow] == 1) {
	for (int iCol=0; iCol<nbCol; iCol++) {
	  T eVal=FAC(iRow, iCol);
	  TestMat(jRow, iCol)=eVal;
	}
	jRow++;
      }
    SelectionRowCol<T> eSelect=TMat_SelectRowCol(TestMat);
    MyMatrix<T> NSP=eSelect.NSP;
    int nbEqua=NSP.rows();
    std::vector<int> eCand, eCandCompl;
    for (int kRow=0; kRow<nbRow; kRow++) {
      int test=1;
      for (int iEqua=0; iEqua<nbEqua; iEqua++)
	if (test == 1) {
	  T eSum=0;
	  for (int iCol=0; iCol<nbCol; iCol++) {
	    T eVal1=FAC(kRow, iCol);
	    T eVal2=NSP(iEqua, iCol);
	    eSum += eVal1*eVal2;
	  }
	  if (eSum != 0)
	    test=0;
	}
      if (test == 1)
	eCand.push_back(kRow);
      else
	eCandCompl.push_back(kRow);
    }
    int nbEltCompl=eCandCompl.size();
    if (nbEltCompl == 0)
      return false;
    int nbElt=eCand.size();
    MyMatrix<T> ListVectSpann(nbElt, nbCol);
    for (int iElt=0; iElt<nbElt; iElt++) {
      int iRow=eCand[iElt];
      for (int iCol=0; iCol<nbCol; iCol++) {
	T eVal=FAC(iRow, iCol);
	ListVectSpann(iElt, iCol)=eVal;
      }
    }
    MyMatrix<T> BasisSpann=RowReduction(ListVectSpann);
    int nbRowSpann=BasisSpann->nbRow;
    int LPdim=nbCol - nbRowSpann;
    MyMatrix<T> PreTheTot=Concatenate(BasisSpann, eMatId);
    MyMatrix<T> TheTot=RowReduction(PreTheTot);
    // the complement
    MyMatrix<T> PreListVectors(nbEltCompl, nbCol);
    for (int iElt=0; iElt<nbEltCompl; iElt++) {
      int iRow=eCandCompl[iElt];
      for (int iCol=0; iCol<nbCol; iCol++) {
	T eVal=FAC(iRow, iCol);
	PreListVectors(iElt, iCol)=eVal;
      }
    }
    MyMatrix<T> TheTotInv=Inverse(TheTot);
    MyMatrix<T> PreListVectorsB=PreListVectors*TheTotInv;
    MyMatrix<T> ListVectors(nbEltCompl, LPdim);
    for (int iElt=0; iElt<nbEltCompl; iElt++)
      for (int iCol=0; iCol<LPdim; iCol++) {
	int jCol=iCol + nbRowSpann;
	T eVal=PreListVectorsB(iElt, jCol);
	ListVectors(iElt, iCol)=eVal;
      }
    PosRelRes<T> eResult=SearchPositiveRelationSimple(ListVectors);
    if (eResult.eTestExist) {
      for (int iElt=0; iElt<nbEltCompl; iElt++) {
	T eVal=eResult.TheRelat(iElt);
	if (eVal > 0)
	  eVectCand[eCandCompl[iElt]]=1;
      }
    }
    else
      return true;
  }
}


template<typename T>
bool TestPositiveRelationSimple(MyMatrix<T> const&ListVect)
{
  PosRelRes<T> eResult=SearchPositiveRelationSimple(ListVect);
  return eResult.eTestExist;
}

template<typename Tgroup>
std::vector<std::vector<int>> GetMinimalReprVertices(Tgroup const& TheGRP)
{
  Face eList;
  int n=TheGRP.n_act();
  for (int i=0; i<n; i++)
    eList[i]=1;
  vectface vvO=DecomposeOrbitPoint(TheGRP, eList);
  std::vector<std::vector<int>> RetList;
  for (auto eOrb : vvO) {
    int MinVal=eOrb.find_first();
    std::vector<int> nList{MinVal};
    RetList.push_back(nList);
  }
  return RetList;
}


template<typename T, typename Tgroup>
std::vector<vectface> EnumerationFaces(Tgroup const& TheGRP, MyMatrix<T> const& FAC, int LevSearch)
{
  std::vector<vectface> RetList;
  int n=TheGRP.n_act();
  vectface ListOrb(n);
  Face eList(n);
  MyMatrix<T> FACred=ColumnReduction(FAC);
  WeightMatrix<T,T> WMat=GetWeightMatrix(FACred);
  for (int i=0; i<n; i++)
    eList[i]=1;
  vectface vvO=DecomposeOrbitPoint(TheGRP, eList);
  for (auto &eOrb : vvO) {
    int MinVal=eOrb.find_first();
    Face nList(n);
    nList[MinVal]=1;
    ListOrb.push_back(nList);
  }
  RetList.push_back(ListOrb);
  for (int iLevel=1; iLevel<=LevSearch; iLevel++) {
    vectface NListOrb(n);
    std::vector<std::vector<int>> ListInv;
    for (auto &eOrb : RetList[iLevel-1]) {
      Tgroup StabFace=TheGRP.Stabilizer_OnSets(eOrb);
      vectface TheSpann=SPAN_face_LinearProgramming(eOrb, StabFace, FACred, TheGRP);
      for (Face fOrb : TheSpann) {
	std::vector<int> eInv=GetLocalInvariantWeightMatrix<T,T,int>(WMat, fOrb);
	GROUP_FuncInsertInSet_UseInv(TheGRP, fOrb, eInv, NListOrb, ListInv);
      }
    }
    RetList.push_back(NListOrb);
  }
  return RetList;
}


void PrintListOrb_GAP(std::ostream &os, std::vector<std::vector<int>> const& ListOrb)
{
  int nbOrb, iOrb, len, i;
  std::vector<int> eOrb;
  os << "return ";
  os << "rec(ListRepresentent:=[";
  nbOrb=ListOrb.size();
  for (iOrb=0; iOrb<nbOrb; iOrb++) {
    if (iOrb>0)
      os << ",";
    os << "[";
    eOrb=ListOrb[iOrb];
    len=eOrb.size();
    for (i=0; i<len; i++) {
      if (i>0)
	os << ",";
      os << eOrb[i]+1;
    }
    os << "]";
  }
  os << "]);\n";
}



void PrintListListOrb_IntGAP(std::ostream &os, std::vector<vectface> const& ListListOrb)
{
  int nbLev=ListListOrb.size();
  os << "return [";
  for (int iLev=0; iLev<nbLev; iLev++) {
    if (iLev > 0)
      os << ",\n";
    os << "rec(ListRepresentent:=[";
    int nbOrb=ListListOrb[iLev].size();
    for (int iOrb=0; iOrb<nbOrb; iOrb++) {
      if (iOrb>0)
	os << ",";
      os << "[";
      Face eOrb=ListListOrb[iLev][iOrb];
      int len=eOrb.count();
      int aRow=eOrb.find_first();
      for (int i=0; i<len; i++) {
	if (i>0)
	  os << ",";
	os << aRow+1;
	aRow=eOrb.find_next(aRow);
      }
      os << "]";
    }
    os << "])";
  }
  os << "];\n";
}
#endif

