#ifndef TEMP_POLYTOPE_FCT_INCLUDE
#define TEMP_POLYTOPE_FCT_INCLUDE

#include "MAT_Matrix.h"
#include "Boost_bitset.h"

template<typename T>
MyVector<T> SumMatrixLineSubset(MyMatrix<T> const& eMat, Face const& eList)
{
  int nbCol=eMat.cols();
  MyVector<T> eVec=ZeroVector<T>(nbCol);
  int eSize=eList.count();
  //
  int aRow=eList.find_first();
  for (int i=0; i<eSize; i++) {
    for (int iCol=0; iCol<nbCol; iCol++)
      eVec(iCol) += eMat(aRow, iCol);
    aRow=eList.find_next(aRow);
  }
  return eVec;
}


template<typename T>
MyMatrix<T> SelectRow(MyMatrix<T> const&TheMat, Face const& eList)
{
  int nbRowRed=eList.count();
  int nbCol=TheMat.cols();
  MyMatrix<T> TheProv(nbRowRed, nbCol);
  int jRow=eList.find_first();
  for (int iRow=0; iRow<nbRowRed; iRow++) {
    TheProv.row(iRow)=TheMat.row(jRow);
    jRow=eList.find_next(jRow);
  }
  return TheProv;
}




template<typename T>
MyMatrix<T> CyclicPolytope(int n, int k)
{
  int i, j, b;
  MyMatrix<T> TheEXT(n, k+1);
  for (i=1; i<=n; i++) {
    T a=1;
    b=i+1;
    for (j=0; j<=k; j++) {
      TheEXT(i-1,j)=a;
      a=a*b;
    }
  }
  return TheEXT;
}


template<typename T>
void TestFacetness(MyMatrix<T> const& EXT, Face const& eList)
{
  MyMatrix<T> TheEXT=ColumnReduction(EXT);
  int nb=eList.count();
  int nbRow=TheEXT.rows();
  int nbCol=TheEXT.cols();
  MyMatrix<T> TheProv(nb, nbCol);
  int aRow=eList.find_first();
  for (int iRow=0; iRow<nb; iRow++) {
    TheProv.row(iRow)=TheEXT.row(aRow);
    aRow=eList.find_next(aRow);
  }
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(TheProv);
  MyMatrix<T> NSP=eSelect.NSP;
  if (NSP.rows() != 1) {
    std::cerr << "Error in rank in Facetness\n";
    std::cerr << "|NSP|=" << NSP.rows() << "\n";
    throw TerminalException{1};
  }
  int nbZero=0;
  int nbPlus=0;
  int nbMinus=0;
  for (int iRow=0; iRow<nbRow; iRow++) {
    T eScal=0;
    for (int iCol=0; iCol<nbCol; iCol++)
      eScal += NSP(0, iCol)*TheEXT(iRow, iCol);
    if (eScal == 0)
      nbZero++;
    if (eScal > 0)
      nbPlus++;
    if (eScal < 0)
      nbMinus++;
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


template<typename T>
MyVector<T> FindFacetInequality(MyMatrix<T> const& TheEXT, Face const& OneInc)
{
  int nb, nbRow, nbCol;
  int nbPlus, nbMinus;
  int iRow, iCol;
  T eScal, eVal, prov;
  nb=OneInc.count();
  //  std::cerr << "nb=" << nb << "\n";
  nbRow=TheEXT.rows();
  nbCol=TheEXT.cols();
  MyMatrix<T> TheProv(nb, nbCol);
  MyVector<T> eVect(nbCol);
  int aRow=OneInc.find_first();
  for (iRow=0; iRow<nb; iRow++) {
    //    std::cerr << "iRow=" << iRow << " aRow=" << aRow << "\n";
    TheProv.row(iRow)=TheEXT.row(aRow);
    aRow=OneInc.find_next(aRow);
  }
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(TheProv);
  MyMatrix<T> NSP=eSelect.NSP;
  if (NSP.rows() != 1) {
    std::cerr << " -------------- BUG TO SOLVE ----------\n";
    std::cerr << "NSP.rows=" << NSP.rows() << "\n";
    std::cerr << "TheProv=\n";
    WriteMatrixGAP(std::cerr, TheProv);
    std::cerr << "\n";
    std::cerr << "NSP=\n";
    WriteMatrixGAP(std::cerr, NSP);
    std::cerr << "\n";
    std::cerr << "Error in rank in FindFacetInequality\n";
    throw TerminalException{1};
  }
  nbPlus=0;
  nbMinus=0;
  for (iRow=0; iRow<nbRow; iRow++) {
    eScal=0;
    for (iCol=0; iCol<nbCol; iCol++) {
      eVal=TheEXT(iRow, iCol);
      prov=NSP(0, iCol);
      eScal=eScal + prov*eVal;
    }
    if (eScal > 0)
      nbPlus++;
    if (eScal < 0)
      nbMinus++;
  }
  if (nbPlus > 0) {
    for (iCol=0; iCol<nbCol; iCol++) {
      prov=NSP(0, iCol);
      eVect(iCol)=prov;
    }
  }
  else {
    for (iCol=0; iCol<nbCol; iCol++) {
      prov=NSP(0, iCol);
      eVal=-prov;
      eVect(iCol)=eVal;
    }
  }
  return eVect;
}


std::vector<int> Dynamic_bitset_to_vectorint(Face const& eList)
{
  //  std::cerr << "Begin of Dynamic_bitset_to_vertorint\n";
  //  int siz=eList.size();
  //  for (int i=0; i<siz; i++)
  //    std::cerr << " i=" << i << " v=" << eList[i] << "\n";
  int nb=eList.count();
  int aRow=eList.find_first();
  std::vector<int> retList;
  for (int i=0; i<nb; i++) {
    retList.push_back(aRow);
    //    std::cerr << "i=" << i << " aRow=" << aRow << "\n";
    aRow=eList.find_next(aRow);
  }
  return retList;
}


template<typename T>
Face ComputeFlipping(MyMatrix<T> const& EXT, Face const& OneInc, Face const& sInc)
{
  MyMatrix<T> TheEXT=ColumnReduction(EXT);
  int nbCol, nbRow;
  std::vector<int> hSub;
  std::vector<Face> TwoPlanes;
  T eVal, prov1, prov2, prov3, prov4;
  T EXT1_1, EXT1_2, EXT2_1, EXT2_2;
  T TheDet, det12, det1N, det2N, prodDet, h;
  int nbForm, IsNonZero, i;
  //  std::cerr << "Begining of ComputeFlipping\n";
  //  std::cerr << "OneInc(count/size)=" << OneInc.count() << "/" << OneInc.size() << "\n";
  //  std::cerr << "sInc(count/size)=" << sInc.count() << "/" << sInc.size() << "\n";
  //  std::cerr << "RankMat(TheEXT)=" << RankMat(TheEXT) << "\n";
  if (OneInc.count() != sInc.size()) {
    std::cerr << "Error in ComputeFlipping\n";
    throw TerminalException{1};
  }
  int nb=sInc.count();
  nbRow=TheEXT.rows();
  nbCol=TheEXT.cols();
  MyMatrix<T> TheProv(nb, nbCol);
  MyMatrix<T> LV(2,2);
  TestFacetness(TheEXT, OneInc);
  std::vector<int> OneInc_V=Dynamic_bitset_to_vectorint(OneInc);
  int jRow=sInc.find_first();
  for (int iRow=0; iRow<nb; iRow++) {
    int aRow=OneInc_V[jRow];
    hSub.push_back(aRow);
    TheProv.row(iRow)=TheEXT.row(aRow);
    jRow=sInc.find_next(jRow);
  }
  //  std::cerr << "RankMat(TheProv)=" << RankMat(TheProv) << "\n";
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(TheProv);
  MyMatrix<T> NSP=eSelect.NSP;
  if (NSP.rows() != 2) {
    int nbRowNSP=NSP.rows();
    std::cerr << "NSP.nbRows=" << nbRowNSP << "\n";
    std::cerr << "Deep inconsistency in ComputeFlipping\n";
    throw TerminalException{1};
  }
  MyMatrix<T> NSPtrans=TransposedMat(NSP);
  MyMatrix<T> LProd=TheEXT*NSPtrans;
  nbForm=0;
  for (int iRow=0; iRow<nbRow; iRow++) {
    IsNonZero=0;
    for (i=0; i<2; i++) {
      prov1=LProd(iRow, i);
      if (prov1 != 0)
	IsNonZero=1;
    }
    if (IsNonZero == 1) {
      prov1=LProd(iRow, 0);
      prov2=LProd(iRow, 1);
      if (nbForm == 0) {
	EXT1_1=prov1;
	EXT1_2=prov2;
	nbForm++;
      }
      else {
	TheDet=prov2*EXT1_1 - prov1*EXT1_2;
	if (nbForm == 1) {
	  if (TheDet != 0) {
	    EXT2_1=prov1;
	    EXT2_2=prov2;
	    det12=TheDet;
	    nbForm++;
	  }
	}
	else {
	  det1N=TheDet;
	  det2N=prov2*EXT2_1 - prov1*EXT2_2;
	  prodDet=det1N*det2N;
	  if (prodDet > 0) {
	    prodDet=det12*det1N;
	    if (prodDet > 0) {
	      EXT2_1=prov1;
	      EXT2_2=prov2;
	      det12=det1N;
	    }
	    else {
	      EXT1_1=prov1;
	      EXT1_2=prov2;
	      det12=-det2N;
	    }
	  }
	}
      }
    }
  }
  if (det12 > 0) {
    eVal=-EXT1_2;
    LV(0, 0)=eVal;
    LV(1, 0)=EXT1_1;
    LV(0, 1)=EXT2_2;
    eVal=-EXT2_1;
    LV(1, 1)=eVal;
  }
  else {
    LV(0, 0)=EXT1_2;
    eVal=-EXT1_1; 
    LV(1, 0)=eVal;
    eVal=-EXT2_2; 
    LV(0, 1)=eVal;
    LV(1, 1)=EXT2_1;
  }
  MyMatrix<T> PairFac=NSPtrans*LV;
  MyMatrix<T> ListScal=TheEXT*PairFac;
  for (i=0; i<2; i++) {
    Face ePlane(nbRow);
    for (int iRow=0; iRow<nbRow; iRow++) {
      eVal=ListScal(iRow, i);
      if (eVal < 0) {
	std::cerr << "This should never have happened. Please panic\n";
	throw TerminalException{1};
      }
      if (eVal == 0)
	ePlane[iRow]=1;
    }
    TwoPlanes.push_back(ePlane);
  }
  int idxFound=-1;
  for (int k=0; k<2; k++)
    if (TwoPlanes[k] == OneInc)
      idxFound=k;
  if (idxFound == -1) {
    std::cerr << "We did not find the facet\n";
    throw TerminalException{1};
  }
  return TwoPlanes[1-idxFound];
}



void PrintListOrbit(std::ostream &os, std::vector<Face> const& ListOrbit)
{
  int iOrbit, nbOrbit, siz, i, eVal;
  nbOrbit=ListOrbit.size();
  os << "nbOrbit=" << nbOrbit << "\n";
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
    Face eInc=ListOrbit[iOrbit];
    siz=eInc.count();
    os << "O" << iOrbit+1 << ": inc=" << siz << "list=";
    eVal=eInc.find_first();
    for (i=0; i<siz; i++) {
      os << " " << eVal;
      eVal=eInc.find_next(eVal);
    }
    os << "\n";
  }
}


struct EngelPolyhedralSubordination {
  int n;
  std::vector<CollectedResult<int>> TheSub;
  std::vector<std::vector<Face>> ListListFace;
};

template<typename T>
EngelPolyhedralSubordination ComputeEngelPolyhedralSubordination(MyMatrix<T> const& EXT, MyMatrix<T> const& FAC)
{
  std::vector<Face> FACset;
  int nbFac=FAC.rows();
  int nbExt=EXT.rows();
  int n=FAC.cols();
  std::cerr << "nbFac=" << nbFac << " nbExt=" << nbExt << " n=" << n << "\n";
  MyVector<T> eFac;
  MyVector<T> eExt;
  for (int iFac=0; iFac<nbFac; iFac++) {
    eFac=FAC.row(iFac);
    Face eFace(nbExt);
    //    std::cerr << "iFac=" << iFac << "\n";
    for (int iExt=0; iExt<nbExt; iExt++) {
      eExt=EXT.row(iExt);
      T eScal=eFac.dot(eExt);
      int eValIns;
      if (eScal == 0)
	eValIns=1;
      else
	eValIns=0;
      eFace[iExt]=eValIns;
      //      std::cerr << "  iExt=" << iExt << " val=" << eValIns << "\n";
    }
    FACset.push_back(eFace);
  }
  std::vector<std::vector<Face> > ListListFace;
  ListListFace.push_back(FACset);
  std::vector<CollectedResult<int> > TheSub;
  for (int eDim=0; eDim<n-1; eDim++) {
    int TheRank=n-2-eDim;
    std::cerr << "eDim=" << eDim << " TheRank=" << TheRank << "\n";
    std::vector<int> ListSizes;
    std::unordered_set<Face> NewListFace_set;
    std::cerr << "  siz=" << ListListFace[eDim].size() << "\n";
    for (auto & eFace : ListListFace[eDim]) {
      int nb=eFace.count();
      std::vector<int> eList(nb);
      int aRow=eFace.find_first();
      for (int i=0; i<nb; i++) {
	eList[i]=aRow;
	aRow=eFace.find_next(aRow);
      }
      std::unordered_set<Face> ListSubFace;
      for (auto & fFace : FACset) {
	Face gFace(nbExt);
	std::vector<int> gList;
	int eIncd=0;
	for (auto & eVal : eList) {
	  if (fFace[eVal] == 1) {
	    gList.push_back(eVal);
	    gFace[eVal]=1;
	    eIncd++;
	  }
	}
	bool IsFace;
	if (eIncd < TheRank) {
	  IsFace=false;
	}
	else {
	  MyMatrix<T> EXTmat=SelectRow(EXT, gList);
	  int rank=RankMat(EXTmat);
	  IsFace=rank == TheRank;
	}
	if (IsFace) {
	  //	  std::cerr << "eIncd=" << eIncd << "\n";
	  ListSubFace.insert(gFace);
	}
      }
      int eSize=ListSubFace.size();
      for (auto & rgFace : ListSubFace)
	NewListFace_set.insert(rgFace);
      ListSizes.push_back(eSize);
    }
    TheSub.push_back(Collected(ListSizes));
    std::vector<Face> NewListFace_vect;
    for (auto & rFace : NewListFace_set)
      NewListFace_vect.push_back(rFace);
    ListListFace.push_back(NewListFace_vect);
  }
  return {n, TheSub, ListListFace};
}

template<typename T>
void ComputeFileFaceLatticeInfo(std::string const& eFile, MyMatrix<T> const& EXT, MyMatrix<T> const& FAC)
{
  EngelPolyhedralSubordination eEngel=ComputeEngelPolyhedralSubordination(EXT, FAC);
  std::ofstream os(eFile);
  os << "return [";
  int len=eEngel.ListListFace.size();
  int nbExt=EXT.rows();
  for (int iDim=0; iDim<len; iDim++) {
    if (iDim>0)
      os << ",\n";
    int nbFace=eEngel.ListListFace[iDim].size();
    os << "[";
    for (int iFace=0; iFace<nbFace; iFace++) {
      if (iFace>0)
	os << ",";
      Face eFace=eEngel.ListListFace[iDim][iFace];
      bool IsFirst=true;
      os << "[";
      for (int iExt=0; iExt<nbExt; iExt++) {
	if (eFace[iExt] == 1) {
	  if (!IsFirst)
	    os << ",";
	  IsFirst=false;
	  int eVal=iExt+1;
	  os << eVal;
	}
      }
      os << "]";
    }
    os << "]";
  }
  os << "];\n";
}



template<typename T>
void ComputeEngelPolyhedralSubordinationFile(std::string const& eFile, MyMatrix<T> const& EXT, MyMatrix<T> const& FAC)
{
  EngelPolyhedralSubordination eEngel=ComputeEngelPolyhedralSubordination(EXT, FAC);
  std::ofstream os(eFile);
  os << "return [";
  int len=eEngel.TheSub.size();
  for (int i=0; i<len; i++) {
    if (i>0)
      os << ",\n";
    CollectedResult<int> eColl=eEngel.TheSub[i];
    int nbEnt=eColl.LVal.size();
    os << "[";
    for (int iEnt=0; iEnt<nbEnt; iEnt++) {
      int eVal=eColl.LVal[iEnt];
      int eMult=eColl.LMult[iEnt];
      if (iEnt > 0)
	os << ",";
      os << "[" << eVal << "," << eMult << "]";
    }
    os << "]";
  }
  os << "];\n";
}




#endif
