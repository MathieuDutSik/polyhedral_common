#ifndef INCLUDE_POLYHEDRAL_TESSELATION
#define INCLUDE_POLYHEDRAL_TESSELATION

#include <vector>
#include <list>

template<typename T>
struct TheBoundary {
  int iRank;
  int lenBoundary;
  int *ListSign;
  int *ListOrbit;
  MyMatrix<T> *ListElt;
};

template<typename T>
struct OneCell {
  MyMatrix<T> *EXT;
  MyVector<T> *IsoEXT;
  list<list<int> > ListIncd;
  /* The thing below is computed by the program */
  MyMatrix<T> *TheSpann;
  MyMatrix<T> *EXTexpSpann;
  MyMatrix<T> *FACspann;
};

template<typename T>
struct AdjacencyStructure {
  int nbFac;
  MyMatrix<T> *FAC;
  int *ListAdjacentDomain;
  MyMatrixPtr *ListMatrix;
};

template<typename T>
struct TotalInformationTesselation {
  int TheRank;
  int n;
  int *ListNbOrbit;
  OneCell **ListCells;
  TheBoundary **ListBoundary;
  AdjacencyStructure *ListAdjStruct;
  MyMatrix<T> *QuadForm;
};

template<typename T>
struct SingleConf {
  int TheStatus;
  int iDomain;
  MyMatrix<T> TheMat;
  MyVector<T> TheVect;
};

/* By complex face, we mean a face that is image under the
   Hecke operator */
template<typename T>
struct ImageFace {
  int iRank;
  int iOrbit;
  MyMatrix<T> *eActMat;
};



template<typename T>
T ComputeDistance(MyVector<T> *eVect1, MyVector<T> *eVect2, TotalInformationTesselation<T> *TotInf)
{
  T prov1, prov2, prov3;
  T eDist;
  int i, j, n;
  n=TotInf->n;
  eDist=0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      {
	prov1=eVect1->ListElt[i] - eVect2->ListElt[i];
	prov2=eVect1->ListElt[j] - eVect2->ListElt[j];
	prov3=prov1*prov2;
	prov1=TotInf->QuadForm(i, j);
	eDist=eDist + prov1*prov3;
      }
  return eDist;
}

template<typename T>
int IsVectorInsideDomain(MyVector<T> *eVectTarget, SingleConf &eConf,
			 TotalInformationTesselation<T> *TotInf)
{
  MyVector<T> *eVectImg;
  int i, n, nbFac, iFac, test;
  int iDomain, WeFound, TheReply;
  T eScal, prov1;
  eVectImg=new MyVector<T>;
  iDomain=eConf->iDomain;
  n=TotInf->n;
  MyMatrix<T> TheInv=Inverse(eConf->TheMat);
  TMat_ImageIntVector(eVectTarget, TheInv, eVectImg);
  nbFac=TotInf->ListAdjStruct[iDomain].nbFac;
  WeFound=-1;
  TheReply=1;
  for (iFac=0; iFac<nbFac; iFac++)
    if (WeFound == -1)
      {
	eScal=0;
	for (i=0; i<n; i++)
	  {
	    prov1=TotInf->ListAdjStruct[iDomain].FAC(iFac, i);
	    eScal=eScal + prov1*eVectImg->ListElt[i];
	  }
	if (eScal < 0)
	  {
	    WeFound=0;
	    TheReply=0;
	  }
      }
  TVec_Free(eVectImg);
  delete eVectImg;
  return TheReply;
}

/* We test if the image of domain eConf (call it P) covers non-trivially
   eFace 
   If F denotes the maximal face of P containing the image, then
   we should have dim F >= dim(eFace).
   F.g included in P.h    equivalent to
   F.u included in P with u=g.h^(-1)
 */
template<typename T> 
int TestIntersectionDomain(ImageFace *eFace, SingleConf &eConf,
			   TotalInformationTesselation<T> *TotInf)
{
  int iDomain, iRank, iOrbit;
  int iFac, nbFac, n;
  MyMatrix<T> *SpannImg, *FACspann, *ListIneq;
  MyVector<T> *eVect;
  MyMatrix<T> *uMat, *uMatCgr, *TheMatProv1, *TheInv;
  list<int> ListIncd;
  list<int>::iterator iter;
  int test, iSp, test2, idxFAC, i;
  int iIneq, nbIneqTotal, nbIneqCell;
  T eScal, prov1, prov2, eSum, eVal;
  eVect=new MyVector<T>;
  ListIneq=new MyMatrix<T>;
  FACspann=new MyMatrix<T>;
  uMat=new MyMatrix<T>;
  uMatCgr=new MyMatrix<T>;
  TheInv=new MyMatrix<T>;
  SpannImg=new MyMatrix<T>;
  TheMatProv1=new MyMatrix<T>;
  iDomain=eConf->iDomain;
  iRank=eFace->iRank;
  iOrbit=eFace->iOrbit;
  n=TotInf->n;
  nbFac=TotInf->ListAdjStruct[iDomain].nbFac;
  TMat_Copy(eConf->TheMat, TheMatProv1);
  TMat_Inverse_destroy(TheMatProv1, TheInv);
  eFace->eActMat=TheInv*uMat;
  TotInf->ListCells[iRank][iOrbit].TheSpann=uMat*SpannImg;
  for (iFac=0; iFac<nbFac; iFac++)
    {
      test=1;
      for (iSp=0; iSp<iRank+1; iSp++)
	{
	  eScal=0;
	  for (i=0; i<n; i++)
	    {
	      pov1=TotInf->ListAdjStruct[iDomain].FAC(iFac, i);
	      prov2=SpannImg(iSp, i);
	      eScal=eScal + prov1*prov2;
	    }
	  if (eScal != 0)
	    test=0;
	}
      if (test == 0)
	ListIncd.push_back(iFac);
    }
  nbIneqCell=TotInf->ListCells[iRank][iOrbit].ListIncd.size();
  nbIneqTotal=nbIneqCell + ListIncd.size();
  TMat_Allocate(ListIneq, nbIneqTotal, iRank+2);
  iter=ListIncd.begin();
  idxFAC=0;
  eVal=-1;
  while(iter != ListIncd.end())
    {
      FACspann(idxFAC, 0)=eVal;
      for (iSp=0; iSp<iRank+1; iSp++)
	{
	  eScal=0;
	  for (i=0; i<n; i++)
	    {
	      prov1=TotInf->ListAdjStruct[iDomain].FAC(iFac, i);
	      prov2=SpannImg(iSp, i);
	      eScal=eScal + prov1*prov2;
	    }
	  ListIneq(idxFAC, iSp+1)=eScal;
	}
      iter++;
      idxFAC++;
    }
  TMat_CongrMap(uMat, uMatCgr);
  for (iIneq=0; iIneq<nbIneqCell; iIneq++) {
    FACspann(idxFAC, 0)=eVal;
    for (iSp=0; iSp<iRank+1; iSp++) {
      prov2=TotInf->ListCells[iRank][iOrbit].FACspann(iIneq, iSp);
      ListIneq(idxFAC, iSp+1)=eScal;
    }
    idxFAC++;
  }
  TVec_Allocate(eVect, iRank+2);
  for (iSp=0; iSp<iRank+2; iSp++) {
    eSum=0;
    for (iIneq=0; iIneq<nbIneqCell; iIneq++) {
      prov2=ListIneq(iIneq, iSp);
      eSum += prov2;
    }
    eVect->ListElt[iSp]=eSum;
  }
  LpSolution<T> eSol=CDD_LinearProgramming(ListIneq, eVect);
  if (eSol.PrimalDefined == true)
    return 1;
  return 0;
}

void CONF_GetAdjacent(SingleConf &eConfInput, SingleConf &eConfOutput,
		      int iFac, TotalInformationTesselation *TotInf)
{
  int iDomain, iDomainAdj;
  MyMatrix<T> *TheProdMat;
  int TheRank;
  TheProdMat=new MyMatrix<T>;
  iDomain=eConfInput->iDomain;
  iDomainAdj=TotInf->ListAdjStruct[iDomain].ListAdjacentDomain[iFac];
  MyMatrix<T> RetMat=TotInf->ListAdjStruct[iDomain].ListMatrix[iFac] * eConfOutput->TheMat;
  ImageIntVector(TotInf->ListCells[TheRank][iDomainAdj].IsoEXT, TheProdMat, eConfOutput->TheVect);
  eConfOutput->TheStatus=1;
  eConfOutput->iDomain=iDomainAdj;
}


void ListSingleConf_FuncInsert(vector<SingleConf> &eFamily,
			       int TheStatus, 
			       int iDomain,
			       MyMatrix<T> &eMat,
			       MyVector<T> &eVect)
{
  SingleConf eConf;
  eConf.TheStatus=TheStatus;
  eConf.iDomain=iDomain;
  TMat_SingleCopy(eMat, eConf.TheMat);
  TVec_SingleCopy(eVect, eConf.TheVect);
  eFamily.push_back(eConf);
}


void LISTCONF_FuncInsert(vector<SingleConf> &ListConf, SingleConf &eConf, 
			 TotalInformationTesselation *TotInf)
{
  int test;
  int iCase, nbCase;
  nbCase=ListConf.size();
  for (iCase=0; iCase<nbCase; iCase++)
    if (eConf->TheVect == ListConf[iCase]->TheVect)
      return;
  ListSingleConf_FuncInsert(ListConf,
			    eConf->TheStatus, 
			    eConf->iDomain,
			    eConf->TheMat,
			    eConf->TheVect);
}


			       

template<typename T>
void FindOneContainingCell(MyVector<T> *eVectTarget,
			   TotalInformationTesselation<T> *TotInf,
			   vector<SingleConf> &CoveringFamily,
			   int iDomainFound, MyMatrix<T> *MatFound)
{
  MyMatrix<T> *RetMat;
  vector<SingleConf> SearchFamily;
  SingleConf eConfOutput;
  T eDist, eDistNew;
  MyMatrix<T> *TheProdMat, *TheMatFound;
  MyVector<T> *eVectImg;
  int idxDomainAdj, iFac, nbFac, test, n, TheStatus;
  int TheRank, iDomain, WeFound, iCase, nbCase;
  n=TotInf->n;
  TheRank=TotInf->TheRank;
  mpq_init(eDist);
  mpq_init(eDistNew);
  /* set to initial point */
  iDomain=0;
  eVectImg=new MyVector<T>;
  TheProdMat=new MyMatrix<T>;
  TheMatFound=new MyMatrix<T>;
  RetMat=new MyMatrix<T>;
  RetMat=IdentityMat<T>(n);
  eDist=ComputeDistance(eVectTarget, TotInf->ListCells[TheRank][iDomain].IsoEXT, TotInf);
  /* Now iterating by using distance function */
  while(1) {
    nbFac=TotInf->ListAdjStruct[iDomain].nbFac;
    WeFound=0;
    for (iFac=0; iFac<nbFac; iFac++)
      if (WeFound == 0) {
	idxDomainAdj=TotInf->ListAdjStruct[iDomain].ListAdjacentDomain[iFac];
	RetMat=TotInf->ListAdjStruct[iDomain].ListMatrix[iFac] * TheProdMat;
	TMat_ImageIntVector(TotInf->ListCells[TheRank][idxDomainAdj].IsoEXT, TheProdMat, eVectImg);
	eDistNew=ComputeDistance(eVectTarget, eVectImg, TotInf);
	if (eDist == eDistNew)
	  WeFound=1;
      }
    if (WeFound == 0)
      break;
    TMat_Copy(RetMat, TheProdMat);
    iDomain=idxDomainAdj;
  }
  /* Now taking all neighbors so as to find one cell covering  */
  TMat_ImageIntVector(TotInf->ListCells[TheRank][idxDomainAdj].IsoEXT, RetMat, eVectImg);
  TheStatus=1;
  ListSingleConf_FuncInsert(SearchFamily, TheStatus, 
			    iDomain, RetMat, eVectImg);
  WeFound=0;
  while(1) {
    nbCase=SearchFamily.size();
    for (iCase=0; iCase<nbCase; iCase++)
      if (SearchFamily[iCase]->TheStatus == 1 && WeFound == 0) {
	test=IsVectorInsideDomain(eVectTarget,SearchFamily[iCase],TotInf);
	iDomain=SearchFamily[iCase]->iDomain;
	if (test == 1) {
	  iDomainFound=iDomain;
	  TMat_Copy(TheMatFound, SearchFamily[iCase]->TheMat);
	  WeFound=1;
	}
	if (WeFound == 0) {
	  nbFac=TotInf->ListAdjStruct[iDomain].nbFac;
	  for (iFac=0; iFac<nbFac; iFac++) {
	    CONF_GetAdjacent(SearchFamily[iCase], eConfOutput,
			     iFac, TotInf);
	    LISTCONF_FuncInsert(SearchFamily, eConfOutput, TotInf);
	  }
	}
      }
  }
}

#endif
