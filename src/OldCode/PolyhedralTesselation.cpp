#include "PolyhedralTesselation.h"

template<typename T>
void ComputeDistance(MyVector<T> *eVect1, MyVector<T> *eVect2, TotalInformationTesselation<T> *TotInf, T &eDist)
{
  T prov1, prov2, prov3;
  int i, j, n;
  n=TotInf->n;
  eDist=0;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      {
	prov1=eVect1->ListElt[i] - eVect2->ListElt[i];
	prov2=eVect1->ListElt[j] - eVect2->ListElt[j];
	prov3=prov1*prov2;
	MATRIX_Get(TotInf->QuadForm, i, j, prov1);
	eDist=eDist + prov1*prov3;
      }
}

template<typename T>
int IsVectorInsideDomain(MyVector<T> *eVectTarget, SingleConf *eConf,
			 TotalInformationTesselation<T> *TotInf)
{
  MyMatrix<T> *TheMatCopy, *TheInv;
  MyVector<T> *eVectImg;
  int i, n, nbFac, iFac, test;
  int iDomain, WeFound, TheReply;
  T eScal, prov1;
  if ((eVectImg = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  if ((TheMatCopy = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheInv = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  iDomain=eConf->iDomain;
  n=TotInf->n;
  TMat_Copy(eConf->TheMat, TheMatCopy);
  TMat_Inverse_destroy(TheMatCopy, TheInv);
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
	    MATRIX_Get(TotInf->ListAdjStruct[iDomain].FAC, iFac, i, prov1);
	    eScal=eScal + prov1*eVectImg->ListElt[i];
	  }
	if (eScal < 0)
	  {
	    WeFound=0;
	    TheReply=0;
	  }
      }
  VECTOR_Free(eVectImg);
  MATRIX_Free(TheMatCopy);
  MATRIX_Free(TheInv);
  free(TheMatCopy);
  free(TheInv);
  free(eVectImg);
  mpq_clear(eScal);
  mpq_clear(prov1);
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
int TestIntersectionDomain(ImageFace *eFace, SingleConf *eConf,
			   TotalInformationTesselation<T> *TotInf)
{
  int iDomain, iRank, iOrbit;
  int iFac, nbFac, n;
  MyMatrix<T> *SpannImg, *FACspann, *ListIneq;
  MyVector<T> *eVect;
  MyMatrix<T> *uMat, *uMatCgr, *TheMatProv1, *TheInv;
  list<int> ListIncd;
  list<int>::iterator iter;
  LpSolution *eSol;
  int test, iSp, test2, idxFAC, i;
  int iIneq, nbIneqTotal, nbIneqCell;
  T eScal, prov1, prov2, eSum, eVal;
  if ((eSol = (LpSolution*)malloc(sizeof(LpSolution))) == 0)
    exit (EXIT_FAILURE);
  if ((eVect = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  if ((ListIneq = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((FACspann = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((uMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((uMatCgr = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheInv = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((SpannImg = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheMatProv1 = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  iDomain=eConf->iDomain;
  iRank=eFace->iRank;
  iOrbit=eFace->iOrbit;
  n=TotInf->n;
  nbFac=TotInf->ListAdjStruct[iDomain].nbFac;
  TMat_Copy(eConf->TheMat, TheMatProv1);
  TMat_Inverse_destroy(TheMatProv1, TheInv);
  TMat_Product(eFace->eActMat, TheInv, uMat);
  TMat_Product(TotInf->ListCells[iRank][iOrbit].TheSpann, uMat, SpannImg);
  for (iFac=0; iFac<nbFac; iFac++)
    {
      test=1;
      for (iSp=0; iSp<iRank+1; iSp++)
	{
	  eScal=0;
	  for (i=0; i<n; i++)
	    {
	      pov1=TMat_Get(TotInf->ListAdjStruct[iDomain].FAC, iFac, i);
	      prov2=TMat_Get(SpannImg, iSp, i);
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
      TMat_Assign(FACspann, idxFAC, 0, eVal);
      for (iSp=0; iSp<iRank+1; iSp++)
	{
	  eScal=0;
	  for (i=0; i<n; i++)
	    {
	      prov1=TMat_Get(TotInf->ListAdjStruct[iDomain].FAC, iFac, i);
	      prov2=TMat_Get(SpannImg, iSp, i);
	      eScal=eScal + prov1*prov2;
	    }
	  TMat_Assign(ListIneq, idxFAC, iSp+1, eScal);
	}
      iter++;
      idxFAC++;
    }
  TMat_CongrMap(uMat, uMatCgr);
  for (iIneq=0; iIneq<nbIneqCell; iIneq++)
    {
      TMat_Assign(FACspann, idxFAC, 0, eVal);
      for (iSp=0; iSp<iRank+1; iSp++)
	{
	  prov2=TMat_Get(TotInf->ListCells[iRank][iOrbit].FACspann, iIneq, iSp);
	  TMat_Assign(ListIneq, idxFAC, iSp+1, eScal);
	}
      idxFAC++;
    }
  TVec_Allocate(eVect, iRank+2);
  for (iSp=0; iSp<iRank+2; iSp++)
    {
      eSum=0;
      for (iIneq=0; iIneq<nbIneqCell; iIneq++)
	{
	  prov2=TMat_Get(ListIneq, iIneq, iSp);
	  eSum=eSum + prov2;
	}
      eVect->ListElt[iSp]=eSum;
    }
  CDD_LinearProgramming(ListIneq, eVect, eSol);
  if (eSol->PrimalDefined)
    return 1;
  return 0;
}

void CONF_GetAdjacent(SingleConf *eConfInput, SingleConf *eConfOutput,
		      int iFac, TotalInformationTesselation *TotInf)
{
  int iDomain, iDomainAdj;
  MyMatrix *RetMat, *TheProdMat;
  int TheRank;
  if ((RetMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheProdMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  iDomain=eConfInput->iDomain;
  iDomainAdj=TotInf->ListAdjStruct[iDomain].ListAdjacentDomain[iFac];
  MATRIX_Product(RetMat, TotInf->ListAdjStruct[iDomain].ListMatrix[iFac], eConfOutput->TheMat);
  ImageIntVector(TotInf->ListCells[TheRank][iDomainAdj].IsoEXT, TheProdMat, eConfOutput->TheVect);
  eConfOutput->TheStatus=1;
  eConfOutput->iDomain=iDomainAdj;
}


void ListSingleConf_FuncInsert(vector<SingleConfPtr> &eFamily,
			       int TheStatus, 
			       int iDomain,
			       MyMatrix *eMat,
			       MyVector *eVect)
{
  int siz;
  SingleConfPtr eConfPtrZ;
  eFamily.push_back(eConfPtrZ);
  siz=eFamily.size();
  if ((eFamily[siz-1] = (SingleConf*)malloc(sizeof(SingleConf))) == 0)
    exit (EXIT_FAILURE);
  eFamily[siz-1]->TheStatus=TheStatus;
  eFamily[siz-1]->iDomain=iDomain;
  if ((eFamily[siz-1]->TheMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((eFamily[siz-1]->TheVect = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  MATRIX_AllocateCopy(eMat, eFamily[siz-1]->TheMat);
  VECTOR_AllocateCopy(eVect, eFamily[siz-1]->TheVect);
}


void LISTCONF_FuncInsert(vector<SingleConfPtr> &ListConf, SingleConf *eConf, 
			 TotalInformationTesselation *TotInf)
{
  int test;
  int iCase, nbCase;
  nbCase=ListConf.size();
  for (iCase=0; iCase<nbCase; iCase++)
    {
      test=VECTOR_IsEqual(eConf->TheVect, ListConf[iCase]->TheVect);
      if (test == 1)
	return;
    }
  ListSingleConf_FuncInsert(ListConf,
			    eConf->TheStatus, 
			    eConf->iDomain,
			    eConf->TheMat,
			    eConf->TheVect);
}

void CONF_Free(SingleConf *eConf)
{
  MATRIX_Free(eConf->TheMat);
  VECTOR_Free(eConf->TheVect);
  free(eConf->TheMat);
  free(eConf->TheVect);
}




			       

template<typename T>
void FindOneContainingCell(MyVector<T> *eVectTarget,
			   TotalInformationTesselation<T> *TotInf,
			   vector<SingleConfPtr> &CoveringFamily,
			   int iDomainFound, MyMatrix<T> *MatFound)
{
  MyMatrix<T> *RetMat;
  vector<SingleConfPtr> SearchFamily;
  SingleConf *eConfOutput = 0;
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
  if ((eVectImg = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  if ((TheProdMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheMatFound = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((RetMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheProdMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  TMat_GetIdentity(RetMat, n);
  ComputeDistance(eVectTarget, TotInf->ListCells[TheRank][iDomain].IsoEXT, TotInf, eDist);
  /* Now iterating by using distance function */
  while(1)
    {
      nbFac=TotInf->ListAdjStruct[iDomain].nbFac;
      WeFound=0;
      for (iFac=0; iFac<nbFac; iFac++)
	if (WeFound == 0)
	  {
	    idxDomainAdj=TotInf->ListAdjStruct[iDomain].ListAdjacentDomain[iFac];
	    TMat_Product(RetMat, TotInf->ListAdjStruct[iDomain].ListMatrix[iFac], TheProdMat);
	    TMat_ImageIntVector(TotInf->ListCells[TheRank][idxDomainAdj].IsoEXT, TheProdMat, eVectImg);
	    ComputeDistance(eVectTarget, eVectImg, TotInf, eDistNew);
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
  while(1)
    {
      nbCase=SearchFamily.size();
      for (iCase=0; iCase<nbCase; iCase++)
	if (SearchFamily[iCase]->TheStatus == 1 && WeFound == 0)
	  {
	    test=IsVectorInsideDomain(eVectTarget,SearchFamily[iCase],TotInf);
	    iDomain=SearchFamily[iCase]->iDomain;
	    if (test == 1)
	      {
		iDomainFound=iDomain;
		TMat_Copy(TheMatFound, SearchFamily[iCase]->TheMat);
		WeFound=1;
	      }
	    if (WeFound == 0)
	      {
		nbFac=TotInf->ListAdjStruct[iDomain].nbFac;
		for (iFac=0; iFac<nbFac; iFac++)
		  {
		    CONF_GetAdjacent(SearchFamily[iCase], eConfOutput,
				     iFac, TotInf);
		    LISTCONF_FuncInsert(SearchFamily, eConfOutput, TotInf);
		  }
	      }
	  }
    }
}
