#include "Skeletton.h"

typedef struct {
  vector<int> ListChoice;
  vector<vector<vector<int> > > ListListPoss;
  int iLevel;
  int nbLevel;
} TreeSearchFace;


// We follow here the conventions of SPAN_face_LinearProgramming
// in Kskeleton.g for the computation.
list<vector<int> > SPAN_face_LinearProgramming(vector<int> face,
					       TheGroupFormat StabFace,
					       MyMatrix* FAC,
					       TheGroupFormat FullGRP)
{
  int nbRow, nbCol;
  int* Treated;
  int* gList;
  int* rList;
  set<int>::iterator iter;
  MyMatrix* TestMat;
  MyMatrix *NSP;
  MyMatrix *PreListVectors;
  MyMatrix *PreListVectorsB;
  MyMatrix *ListVectors;
  MyMatrix *PreTheTot;
  MyMatrix *TheTot;
  MyMatrix *TheTotInv;
  MyMatrix *ListVectSpann;
  MyMatrix *BasisSpann;
  MyMatrix *eMatId;
  MyMatrix* FACred;
  MySelection eSelect;
  vector<int> eCand, eCandCompl;
  mpq_t eVal, eVal1, eVal2, eSum;
  list<vector<int> > TheReturn;
  int iRow, jRow, ePt, iCol, jCol, iElt;
  int sizFace, orbSiz, iOrb, nbOrb, nbOrbCompl;
  int LPdim, eTestExist;
  int test, testB, nbEqua, iEqua, i;
  vector<vector<int> > ListOrb, ListOrbCompl;
  vector<int> eOrb;
  TheGroupFormat TheStab;
  int nbRowSpann;
  mpq_init(eVal);
  mpq_init(eVal1);
  mpq_init(eVal2);
  mpq_init(eSum);
  if ((TestMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((NSP = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((PreListVectors = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((PreListVectorsB = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((ListVectors = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((PreTheTot = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheTot = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheTotInv = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((ListVectSpann = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((BasisSpann = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((eMatId = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((FACred = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  MATRIX_ColReduction(FAC, FACred);


  TheReturn.clear();
  nbRow=FACred->nbRow;
  nbCol=FACred->nbCol;
  MATRIX_Allocate(TheTotInv, nbCol, nbCol);
  MATRIX_GetIdentity(eMatId, nbCol);
  if ((Treated = (int*)malloc(nbRow*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((rList = (int*)malloc(nbRow*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((gList = (int*)malloc(nbRow*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (iRow=0; iRow<nbRow; iRow++)
    Treated[iRow]=0;
  sizFace=face.size();
  for (i=0; i<sizFace; i++)
    {
      ePt=face[i];
      Treated[ePt]=1;
    }
  MATRIX_Allocate(TestMat, sizFace+1, nbCol);
  jRow=0;
  for (jRow=0; jRow<sizFace; jRow++)
    {
      iRow=face[jRow];
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(FACred, iRow, iCol, eVal);
	  MATRIX_Assign(TestMat, jRow, iCol, eVal);
	}
    }
  for (iRow=0; iRow<nbRow; iRow++)
    if (Treated[iRow] == 0)
      {
	for (iCol=0; iCol<nbCol; iCol++)
	  {
	    MATRIX_Get(FACred, iRow, iCol, eVal);
	    MATRIX_Assign(TestMat, sizFace, iCol, eVal);
	  }
	MATRIX_SelectColRow(TestMat, &eSelect, NSP);
	SELECT_Free(&eSelect);
	nbEqua=NSP->nbRow;
	eCand.clear();
	eCandCompl.clear();
	for (jRow=0; jRow<nbRow; jRow++)
	  {
	    test=1;
	    for (iEqua=0; iEqua<nbEqua; iEqua++)
	      if (test == 1)
		{
		  mpq_set_si(eSum, 0, 1);
		  for (iCol=0; iCol<nbCol; iCol++)
		    {
		      MATRIX_Get(FACred, jRow, iCol, eVal1);
		      MATRIX_Get(NSP, iEqua, iCol, eVal2);
		      mpq_mul(eVal1, eVal1, eVal2);
		      mpq_add(eSum, eSum, eVal1);
		    }
		  testB=mpq_sgn(eSum);
		  if (testB != 0)
		    test=0;
		}
	    gList[jRow]=test;
	    if (test == 1)
	      eCand.push_back(jRow);
	    else
	      eCandCompl.push_back(jRow);
	  }
	OrbitUnion(StabFace, gList, rList);
	for (jRow=0; jRow<nbRow; jRow++)
	  if (rList[jRow] == 1)
	    Treated[jRow]=1;
	GetStabilizer(FullGRP, eCand, TheStab);
	ListOrb=DecomposeOrbitPoint(TheStab, eCand);
	nbOrb=ListOrb.size();
	MATRIX_Allocate(ListVectSpann, nbOrb, nbCol);
	for (iOrb=0; iOrb<nbOrb; iOrb++)
	  {
	    eOrb=ListOrb[iOrb];
	    orbSiz=eOrb.size();
	    for (iCol=0; iCol<nbCol; iCol++)
	      {
		mpq_set_si(eSum, 0, 1);
		for (iElt=0; iElt<orbSiz; iElt++)
		  {
		    jRow=eOrb[iElt];
		    MATRIX_Get(FACred, jRow, iCol, eVal1);
		    mpq_add(eSum, eSum, eVal1);
		  }
		MATRIX_Assign(ListVectSpann, iOrb, iCol, eSum);
	      }
	  }
	MATRIX_RowReduction(ListVectSpann, BasisSpann);
	MATRIX_Free(ListVectSpann);
	nbRowSpann=BasisSpann->nbRow;
	LPdim=nbCol - nbRowSpann;
	MATRIX_Concatenate(BasisSpann, eMatId, PreTheTot);
	MATRIX_Free(BasisSpann);
	MATRIX_RowReduction(PreTheTot, TheTot);
	MATRIX_Free(PreTheTot);
	// the complement
	ListOrbCompl=DecomposeOrbitPoint(TheStab, eCandCompl);
	nbOrbCompl=ListOrbCompl.size();
	MATRIX_Allocate(PreListVectors, nbOrbCompl, nbCol);
	for (iOrb=0; iOrb<nbOrbCompl; iOrb++)
	  {
	    eOrb=ListOrbCompl[iOrb];
	    orbSiz=eOrb.size();
	    for (iCol=0; iCol<nbCol; iCol++)
	      {
		mpq_set_si(eSum, 0, 1);
		for (iElt=0; iElt<orbSiz; iElt++)
		  {
		    jRow=eOrb[iElt];
		    MATRIX_Get(FACred, jRow, iCol, eVal1);
		    mpq_add(eSum, eSum, eVal1);
		  }
		MATRIX_Assign(PreListVectors, iOrb, iCol, eSum);
	      }
	  }
	MATRIX_Inverse(TheTot, TheTotInv);
	MATRIX_Free(TheTot);
	MATRIX_Product(PreListVectors, TheTotInv, PreListVectorsB);
	MATRIX_Free(PreListVectors);
	MATRIX_Allocate(ListVectors, nbOrbCompl, LPdim);
	for (jRow=0; jRow<nbOrbCompl; jRow++)
	  for (iCol=0; iCol<LPdim; iCol++)
	    {
	      jCol=iCol + nbRowSpann;
	      MATRIX_Get(PreListVectorsB, jRow, jCol, eVal);
	      MATRIX_Assign(ListVectors, jRow, iCol, eVal);
	    }
	MATRIX_Free(PreListVectorsB);
	eTestExist=TestPositiveRelationSimple(ListVectors);
	MATRIX_Free(ListVectors);
	if (eTestExist == 0)
	  TheReturn.push_back(eCand);
	MATRIX_Free(NSP);
      }
  mpq_clear(eVal);
  mpq_clear(eVal1);
  mpq_clear(eVal2);
  mpq_clear(eSum);
  free(Treated);
  free(rList);
  free(gList);

  MATRIX_Free(TestMat);
  free(TestMat);
  free(NSP);
  free(PreListVectors);
  free(PreListVectorsB);
  free(ListVectors);
  free(PreTheTot);
  free(TheTot);
  MATRIX_Free(TheTotInv);
  free(TheTotInv);
  free(ListVectSpann);
  free(BasisSpann);
  MATRIX_Free(eMatId);
  free(eMatId);
  MATRIX_Free(FACred);
  free(FACred);
  return TheReturn;
}

// We test if eSet is included in a proper face of the polytope
int TestInclusionProperFace(vector<int> eSet, MyMatrix* FAC)
{
  int nbCol, nbRow, iRow, len, i, jRow, iCol, jCol;
  int nbEqua, iEqua, nbEltCompl;
  int test, testB, eTestExist;
  int nbRowSpann;
  int iElt, nbElt, LPdim;
  vector<int> eCand, eCandCompl;
  vector<int> eVectCand;
  MyVector *TheRelat;
  MyMatrix *ListVectors;
  MyMatrix *PreListVectors;
  MyMatrix *PreListVectorsB;
  MyMatrix *BasisSpann;
  MyMatrix *ListVectSpann;
  MyMatrix *eMatId;
  MyMatrix *TheTot;
  MyMatrix *PreTheTot;
  MyMatrix *TheTotInv;
  MyMatrix *TestMat;
  MyMatrix *NSP;
  MySelection eSelect;
  mpq_t eVal, eVal1, eVal2, eSum;
  mpq_init(eVal);
  mpq_init(eVal1);
  mpq_init(eVal2);
  mpq_init(eSum);
  nbCol=FAC->nbCol;
  nbRow=FAC->nbRow;
  fprintf(stderr, "nbRow=%d nbCol=%d\n", nbRow, nbCol);
  fprintf(stderr, "eSet=");
  PrintVectorInt_GAP(stderr, eSet);
  fprintf(stderr, "\n");
  for (iRow=0; iRow<nbRow; iRow++)
    eVectCand.push_back(0);
  len=eSet.size();
  for (i=0; i<len; i++)
    eVectCand[eSet[i]]=1;
  if ((PreListVectors = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((PreListVectorsB = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheRelat = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  if ((ListVectors = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((ListVectSpann = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((BasisSpann = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((eMatId = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((PreTheTot = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheTot = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheTotInv = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TestMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((NSP = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  MATRIX_Allocate(TheTotInv, nbCol, nbCol);
  MATRIX_GetIdentity(eMatId, nbCol);
  while(1)
    {
      len=0;
      for (iRow=0; iRow<nbRow; iRow++)
	if (eVectCand[iRow] == 1)
	  len++;
      fprintf(stderr, "TestInclusionProperFace, step 1\n");
      MATRIX_Allocate(TestMat, len, nbCol);
      fprintf(stderr, "len=%d nbCol=%d\n", len, nbCol);
      jRow=0;
      for (iRow=0; iRow<nbRow; iRow++)
	if (eVectCand[iRow] == 1)
	  {
	    for (iCol=0; iCol<nbCol; iCol++)
	      {
		MATRIX_Get(FAC, iRow, iCol, eVal);
		MATRIX_Assign(TestMat, jRow, iCol, eVal);
	      }
	    jRow++;
	  }
      fprintf(stderr, "TestInclusionProperFace, step 2\n");
      MATRIX_SelectColRow(TestMat, &eSelect, NSP);
      MATRIX_Free(TestMat);
      SELECT_Free(&eSelect);
      nbEqua=NSP->nbRow;
      eCand.clear();
      eCandCompl.clear();
      fprintf(stderr, "TestInclusionProperFace, step 3\n");
      for (jRow=0; jRow<nbRow; jRow++)
	{
	  test=1;
	  for (iEqua=0; iEqua<nbEqua; iEqua++)
	    if (test == 1)
	      {
		mpq_set_si(eSum, 0, 1);
		for (iCol=0; iCol<nbCol; iCol++)
		  {
		    MATRIX_Get(FAC, jRow, iCol, eVal1);
		    MATRIX_Get(NSP, iEqua, iCol, eVal2);
		    mpq_mul(eVal1, eVal1, eVal2);
		    mpq_add(eSum, eSum, eVal1);
		  }
		testB=mpq_sgn(eSum);
		if (testB != 0)
		  test=0;
	      }
	  if (test == 1)
	    eCand.push_back(jRow);
	  else
	    eCandCompl.push_back(jRow);
	}
      MATRIX_Free(NSP);
      fprintf(stderr, "TestInclusionProperFace, step 4\n");
      nbEltCompl=eCandCompl.size();
      if (nbEltCompl == 0)
	{
	  free(PreListVectors);
	  free(PreListVectorsB);
	  free(ListVectors);
	  free(ListVectSpann);
	  free(BasisSpann);
	  MATRIX_Free(eMatId);
	  free(eMatId);
	  free(TheTot);
	  free(PreTheTot);
	  MATRIX_Free(TheTotInv);
	  free(TheTotInv);
	  free(TestMat);
	  free(NSP);
	  free(TheRelat);
	  mpq_clear(eVal);
	  mpq_clear(eVal1);
	  mpq_clear(eVal2);
	  mpq_clear(eSum);
	  return 0;
	}
      fprintf(stderr, "TestInclusionProperFace, step 5\n");
      nbElt=eCand.size();
      MATRIX_Allocate(ListVectSpann, nbElt, nbCol);
      for (iElt=0; iElt<nbElt; iElt++)
	{
	  iRow=eCand[iElt];
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      MATRIX_Get(FAC, iRow, iCol, eVal);
	      MATRIX_Assign(ListVectSpann, iElt, iCol, eVal);
	    }
	}
      fprintf(stderr, "TestInclusionProperFace, step 6\n");
      MATRIX_RowReduction(ListVectSpann, BasisSpann);
      fprintf(stderr, "TestInclusionProperFace, step 6.1\n");
      MATRIX_Free(ListVectSpann);
      fprintf(stderr, "TestInclusionProperFace, step 6.2\n");
      nbRowSpann=BasisSpann->nbRow;
      fprintf(stderr, "TestInclusionProperFace, step 6.3\n");
      LPdim=nbCol - nbRowSpann;
      fprintf(stderr, "TestInclusionProperFace, step 6.4\n");
      MATRIX_Concatenate(BasisSpann, eMatId, PreTheTot);
      fprintf(stderr, "TestInclusionProperFace, step 6.5\n");
      MATRIX_Free(BasisSpann);
      fprintf(stderr, "TestInclusionProperFace, step 6.5\n");
      MATRIX_RowReduction(PreTheTot, TheTot);
      fprintf(stderr, "TestInclusionProperFace, step 6.6\n");
      MATRIX_Free(PreTheTot);
      fprintf(stderr, "TestInclusionProperFace, step 6.7\n");
      // the complement
      MATRIX_Allocate(PreListVectors, nbEltCompl, nbCol);
      fprintf(stderr, "TestInclusionProperFace, step 7\n");
      for (iElt=0; iElt<nbEltCompl; iElt++)
	{
	  iRow=eCandCompl[iElt];
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      MATRIX_Get(FAC, iRow, iCol, eVal);
	      MATRIX_Assign(PreListVectors, iElt, iCol, eVal);
	    }
	}
      fprintf(stderr, "TestInclusionProperFace, step 8\n");
      MATRIX_Inverse(TheTot, TheTotInv);
      fprintf(stderr, "TestInclusionProperFace, step 8.1\n");
      MATRIX_Free(TheTot);
      fprintf(stderr, "TestInclusionProperFace, step 8.2\n");
      MATRIX_Product(PreListVectors, TheTotInv, PreListVectorsB);
      fprintf(stderr, "TestInclusionProperFace, step 8.3\n");
      MATRIX_Free(PreListVectors);
      fprintf(stderr, "TestInclusionProperFace, step 8.4\n");
      MATRIX_Allocate(ListVectors, nbEltCompl, LPdim);
      fprintf(stderr, "TestInclusionProperFace, step 9\n");
      for (iElt=0; iElt<nbEltCompl; iElt++)
	for (iCol=0; iCol<LPdim; iCol++)
	  {
	    jCol=iCol + nbRowSpann;
	    MATRIX_Get(PreListVectorsB, iElt, jCol, eVal);
	    MATRIX_Assign(ListVectors, iElt, iCol, eVal);
	  }
      fprintf(stderr, "TestInclusionProperFace, step 10\n");
      MATRIX_Free(PreListVectorsB);
      VECTOR_Allocate(TheRelat, nbEltCompl);
      SearchPositiveRelationSimple(ListVectors, &eTestExist, TheRelat);
      MATRIX_Free(ListVectors);
      if (eTestExist == 1)
	{
	  for (iElt=0; iElt<nbEltCompl; iElt++)
	    {
	      VECTOR_Get(TheRelat, iElt, eVal);
	      test=mpq_sgn(eVal);
	      if (test == 1)
		eVectCand[eCandCompl[iElt]]=1;
	    }
	  VECTOR_Free(TheRelat);
	}
      else
	{
	  VECTOR_Free(TheRelat);
	  free(PreListVectors);
	  free(PreListVectorsB);
	  free(ListVectors);
	  free(ListVectSpann);
	  free(BasisSpann);
	  MATRIX_Free(eMatId);
	  free(eMatId);
	  free(TheTot);
	  free(PreTheTot);
	  MATRIX_Free(TheTotInv);
	  free(TheTotInv);
	  free(TestMat);
	  free(NSP);
	  free(TheRelat);
	  mpq_clear(eVal);
	  mpq_clear(eVal1);
	  mpq_clear(eVal2);
	  mpq_clear(eSum);
	  return 1;
	}
      fprintf(stderr, "TestInclusionProperFace, step 11\n");
    }
}




void SearchPositiveRelationSimple(MyMatrix *ListVect, int *eTestExist, MyVector *TheRelat)
{
  int iVect, iRel;
  int nbRelation, nbVect;
  mpq_t eVal, eSum;
  MyMatrix NSP;
  MyMatrix* ListVectTr;
  MyMatrix* ListInequalities;
  MyVector *ToBeMinimized;
  MyVector *eVectRel;
  MySelection eSelect;
  LpSolution *eSol;
  int test;
  if ((ListVectTr = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((eSol = (LpSolution*)malloc(sizeof(LpSolution))) == 0)
    exit (EXIT_FAILURE);
  if ((ListInequalities = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((ToBeMinimized = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  MATRIX_Transpose(ListVect, ListVectTr);
  mpq_init(eVal);
  mpq_init(eSum);
  MATRIX_SelectColRow(ListVectTr, &eSelect, &NSP);
  nbVect=ListVect->nbRow;
  nbRelation=NSP.nbRow;
  MATRIX_Allocate(ListInequalities, nbVect+1, nbRelation+1);
  for (iVect=0; iVect<nbVect; iVect++)
    {
      mpq_set_si(eVal, 0, 1);
      MATRIX_Assign(ListInequalities, iVect, 0, eVal);
      for (iRel=0; iRel<nbRelation; iRel++)
	{
	  MATRIX_Get(&NSP, iRel, iVect, eVal);
	  MATRIX_Assign(ListInequalities, iVect, iRel+1, eVal);
	}
    }
  mpq_set_si(eVal, -1, 1);
  MATRIX_Assign(ListInequalities, nbVect, 0, eVal);
  VECTOR_Allocate(ToBeMinimized, nbRelation+1);
  mpq_set_si(eVal, 0, 1);
  VECTOR_Assign(ToBeMinimized, 0, eSum);
  for (iRel=0; iRel<nbRelation; iRel++)
    {
      mpq_set_si(eSum, 0, 1);
      for (iVect=0; iVect<nbVect; iVect++)
	{
	  MATRIX_Get(&NSP, iRel, iVect, eVal);
	  mpq_add(eSum, eSum, eVal);
	}
      MATRIX_Assign(ListInequalities, nbVect, iRel+1, eSum);
      VECTOR_Assign(ToBeMinimized, iRel+1, eSum);
    }
  //
  CDD_LinearProgramming(ListInequalities, ToBeMinimized, eSol);
  if (eSol->PrimalDefined && eSol->DualDefined)
    {
      *eTestExist=1;
      if ((eVectRel = (MyVector*)malloc(sizeof(MyVector))) == 0)
	exit (EXIT_FAILURE);
      VECTOR_Allocate(eVectRel, nbRelation);
      for (iRel=0; iRel<nbRelation; iRel++)
	VECTOR_Assign(eVectRel, iRel, eSol->DirectSolutionVal[iRel+1]);
      VECTOR_MATRIX_Product(eVectRel, &NSP, TheRelat);
      VECTOR_Free(eVectRel);
      free(eVectRel);
      for (iVect=0; iVect<nbVect; iVect++)
	{
	  VECTOR_Get(TheRelat, iVect, eVal);
	  test=mpq_sgn(eVal);
	  if (test == -1)
	    {
	      fprintf(stderr, "We have a clear and present error\n");
	      exit(1);
	    }
	}
    }
  else
    {
      *eTestExist=0;
    }
  LPSOL_Free(eSol);
  free(eSol);
  SELECT_Free(&eSelect);
  VECTOR_Free(ToBeMinimized);
  free(ToBeMinimized);
  MATRIX_Free(ListInequalities);
  free(ListInequalities);
  //
  MATRIX_Free(&NSP);
  MATRIX_Free(ListVectTr);
  free(ListVectTr);
  mpq_clear(eVal);
  mpq_clear(eSum);
}


int TestPositiveRelationSimple(MyMatrix *ListVect)
{
  int eTestExist, nbRow;
  MyVector TheRelat;
  nbRow=ListVect->nbRow;
  VECTOR_Allocate(&TheRelat, nbRow);
  SearchPositiveRelationSimple(ListVect, &eTestExist, &TheRelat);
  VECTOR_Free(&TheRelat);
  return eTestExist;
}


// Test if eFace is a lexicographically minimal face in
// the polytope
int TestLexicographicMinimality(vector<int> eFace, MyMatrix* FAC)
{
  int iRow, jRow, nbRow;
  int len, i, MaxVal, ePos, test;
  vector<int> ListStatus, eSet;
  nbRow=FAC->nbRow;
  for (iRow=0; iRow<nbRow; iRow++)
    ListStatus.push_back(1);
  len=eFace.size();
  MaxVal=0;
  for (i=0; i<len; i++)
    {
      ePos=eFace[i];
      ListStatus[ePos]=0;
      if (ePos > MaxVal)
	MaxVal=ePos;
    }
  for (iRow=0; iRow<MaxVal; iRow++)
    if (ListStatus[iRow] == 1)
      {
	eSet.clear();
	for (jRow=0; jRow<iRow; jRow++)
	  if (ListStatus[jRow] == 0)
	    eSet.push_back(jRow);
	eSet.push_back(iRow);
	fprintf(stderr, "Before TestInclusionProperFace\n");
	test=TestInclusionProperFace(eSet, FAC);
	fprintf(stderr, "After TestInclusionProperFace\n");
	if (test == 1)
	  return 0;
      }
  return 1;
}

vector<vector<int> > GetMinimalReprVertices(TheGroupFormat TheGRP)
{
  int nbOrb, iOrb, MinVal;
  int n, i, siz;
  vector<vector<int> > vvO;
  vector<int> nList, eList, eOrb;
  vector<vector<int> > RetList;
  n=TheGRP.n;
  for (i=0; i<n; i++)
    eList.push_back(i);
  vvO=DecomposeOrbitPoint(TheGRP, eList);
  nbOrb=vvO.size();
  for (iOrb=0; iOrb<nbOrb; iOrb++)
    {
      eOrb=vvO[iOrb];
      siz=eOrb.size();
      MinVal=eOrb[0];
      for (i=1; i<siz; i++)
	    if (eOrb[i] < MinVal)
	      MinVal=eOrb[i];
      nList.clear();
      nList.push_back(MinVal);
      RetList.push_back(nList);
    }
  return RetList;
}

vector<vector<int> > GetLexicographicChoices(vector<int> eFace, TheGroupFormat TheGRP, MyMatrix* FAC)
{
  vector<vector<int> > RetList;
  int len, i, test;
  vector<int> eOrb, fOrb;
  int iRow, jRow, iCol, nbRow, nbCol, sizMat;
  vector<int> eVectOrig, eList, eFaceRed;
  mpq_t eVal;
  list<vector<int> >::iterator iter;
  list<vector<int> > TheSpann;
  MyMatrix *FACface;
  TheGroupFormat StabFace;
  len=eFace.size();
  mpq_init(eVal);
  if (len == 0)
    {
      fprintf(stderr, "Clear error BOMB\n");
      exit(1);
    }
  if ((FACface = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  nbCol=FAC->nbCol;
  nbRow=FAC->nbRow;
  for (iRow=0; iRow<nbRow; iRow++)
    eVectOrig.push_back(-1);
  GetStabilizer(TheGRP, eFace, StabFace);
  TheSpann=SPAN_face_LinearProgramming(eFace, StabFace, FAC, TheGRP);
  iter=TheSpann.begin();
  while(iter != TheSpann.end())
    {
      fOrb=*iter;
      sizMat=fOrb.size();
      MATRIX_Allocate(FACface, sizMat, nbCol);
      fprintf(stderr, "sizMat=%d\n", sizMat);
      for (iRow=0; iRow<nbRow; iRow++)
	eVectOrig[iRow]=-1;
      for (iRow=0; iRow<sizMat; iRow++)
	{
	  jRow=fOrb[iRow];
	  eVectOrig[jRow]=iRow;
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      MATRIX_Get(FAC, jRow, iCol, eVal);
	      MATRIX_Assign(FACface, iRow, iCol, eVal);
	    }
	}
      len=eFace.size();
      eFaceRed.clear();
      for (i=0; i<len; i++)
	eFaceRed.push_back(eVectOrig[eFace[i]]);
      fprintf(stderr, "Before TestLexicographicMinimality\n");
      fprintf(stderr, "eFaceRed=");
      PrintVectorInt_GAP(stderr, eFaceRed);
      fprintf(stderr, "\n");
      fprintf(stderr, "FACface nbRow=%d nbCol=%d\n", FACface->nbRow, FACface->nbCol);
      test=TestLexicographicMinimality(eFaceRed, FACface);
      fprintf(stderr, " After TestLexicographicMinimality\n");
      if (test == 1)
	RetList.push_back(fOrb);
      iter++;
      MATRIX_Free(FACface);
    }
  free(FACface);
  mpq_clear(eVal);
  return RetList;
}

int GoUpNextInTree(TreeSearchFace &TheTR, TheGroupFormat TheGRP, MyMatrix* FAC)
{
  int iLevel, NbChoice, idx;
  iLevel=TheTR.iLevel;
  NbChoice=TheTR.ListListPoss[iLevel-1].size();
  idx=TheTR.ListChoice[iLevel-1];
  if (idx == NbChoice)
    {
      if (iLevel == 1)
	return 0;
      TheTR.iLevel=iLevel-1;
      return GoUpNextInTree(TheTR, TheGRP, FAC);
    }
  TheTR.ListChoice[iLevel-1]=idx+1;
  return 1;
}

int NextInTree(TreeSearchFace &TheTR, TheGroupFormat TheGRP, MyMatrix* FAC)
{
  vector<int> eFace;
  int iLevel, nbLevel, nbChoice, idx;
  iLevel=TheTR.iLevel;
  nbLevel=TheTR.nbLevel;
  fprintf(stderr, "NextInTree iLevel=%d nLevel=%d\n", iLevel, nbLevel);
  if (iLevel == nbLevel)
    return GoUpNextInTree(TheTR, TheGRP, FAC);
  fprintf(stderr, "NextInTree step 2\n");
  nbChoice=TheTR.ListListPoss[iLevel-1].size();
  fprintf(stderr, "NextInTree nbChoice=%d\n", nbChoice);
  fprintf(stderr, "NextInTree step 3\n");
  idx=TheTR.ListChoice[iLevel-1];
  fprintf(stderr, "NextInTree idx=%d\n", idx);
  eFace=TheTR.ListListPoss[iLevel-1][idx-1];
  fprintf(stderr, "NextInTree step 3.1\n");
  fprintf(stderr, "NextInTree eFace=");
  PrintVectorInt_GAP(stderr, eFace);
  fprintf(stderr, "\n");
  fprintf(stderr, "NextInTree step 4\n");
  TheTR.ListListPoss[iLevel]=GetLexicographicChoices(eFace, TheGRP, FAC);
  nbChoice=TheTR.ListListPoss[iLevel].size();
  if (nbChoice == 0)
    return GoUpNextInTree(TheTR, TheGRP, FAC);
  TheTR.ListChoice[iLevel]=1;
  TheTR.iLevel=iLevel+1;
  fprintf(stderr, "NextInTree step 5\n");
  return 1;
}

void InitializeTR(TreeSearchFace &TheTR, TheGroupFormat TheGRP, int LevSearch)
{
  vector<vector<int> > eList;
  int eInt, iLevel;
  fprintf(stderr, "InitializeTR, step 1\n");
  eInt=-1;
  fprintf(stderr, "InitializeTR, step 2\n");
  eList.clear();
  fprintf(stderr, "InitializeTR, step 3\n");
  TheTR.nbLevel=LevSearch;
  fprintf(stderr, "InitializeTR, step 4\n");
  TheTR.iLevel=1;
  fprintf(stderr, "InitializeTR, step 5\n");
  for (iLevel=0; iLevel<=LevSearch; iLevel++)
    {
      TheTR.ListChoice.push_back(eInt);
      TheTR.ListListPoss.push_back(eList);
    }
  fprintf(stderr, "InitializeTR, step 6\n");
  TheTR.ListListPoss[0]=GetMinimalReprVertices(TheGRP);
  TheTR.ListChoice[0]=1;
  fprintf(stderr, "InitializeTR, step 7\n");
}



void Generic_DoMemoryEfficientEnum(TheGroupFormat TheGRP, MyMatrix* FAC, 
				   int LevSearch, 
				   void (*hook)(void *user_param, vector<int> eFace),
				   void *hook_user_format)
{
  int result, idx;
  vector<int> eFace;
  TreeSearchFace TheTR;
  int iter;
  InitializeTR(TheTR, TheGRP, LevSearch);
  fprintf(stderr, "Start main iteration\n");
  iter=0;
  while(1)
    {
      iter++;
      fprintf(stderr, "Before NextInTree iter=%d\n", iter);
      result=NextInTree(TheTR, TheGRP, FAC);
      fprintf(stderr, "After NextInTree\n");
      if (TheTR.iLevel == LevSearch)
	{
	  idx=TheTR.ListChoice[LevSearch-1];
	  fprintf(stderr, "After test idx=%d\n", idx);
	  eFace=TheTR.ListListPoss[LevSearch-1][idx-1];
	  (*hook)(hook_user_format, eFace);
	}
      if (result == 0)
	break;
    }
}

void Generic_DoMemoryEfficientAllEnum(TheGroupFormat TheGRP, MyMatrix* FAC, 
				      int LevSearch, 
				      void (*hook)(void *user_param, int iLev, vector<int> eFace),
				      void *hook_user_format)
{
  int result, idx;
  vector<int> eFace;
  TreeSearchFace TheTR;
  int iter, iLevel;
  int len;
  InitializeTR(TheTR, TheGRP, LevSearch);
  fprintf(stderr, "Start main iteration\n");
  iter=0;
  while(1)
    {
      iter++;
      fprintf(stderr, "Before NextInTree iter=%d\n", iter);
      result=NextInTree(TheTR, TheGRP, FAC);
      fprintf(stderr, "After NextInTree\n");
      iLevel=TheTR.iLevel;
      fprintf(stderr, "Now iLevel=%d\n", iLevel);
      idx=TheTR.ListChoice[iLevel-1];
      fprintf(stderr, "Now idx=%d\n", idx);
      len=TheTR.ListListPoss[iLevel-1].size();
      if (idx > len)
	{
	  fprintf(stderr, "len=%d idx=%d\n", len, idx);
	  fprintf(stderr, "Stop here\n");
	  exit(1);
	}
      eFace=TheTR.ListListPoss[iLevel-1][idx-1];
      fprintf(stderr, "eFace=");
      PrintVectorInt_GAP(stderr, eFace);
      fprintf(stderr, "\n");
      (*hook)(hook_user_format, iLevel, eFace);
      fprintf(stderr, "After hook operation\n");
      if (result == 0)
	break;
    }
}



typedef vector<vector<int> > VectVectInt;
typedef vector<vector<vector<int> > > VectVectVectInt;
static void report_face_vectvectint(void* param, vector<int> eVect)
{
  ((VectVectInt*)param)->push_back(eVect);
}

static void report_allface_vectvectint(void* param, int iLevel, vector<int> eVect)
{
  (*((VectVectVectInt*)param))[iLevel-1].push_back(eVect);
}

/*
static void report_face_print(void* param, vector<int> eVect)
{
  int len, i;
  len=eVect.size();
  fprintf((FILE*)param, "[");
  for (i=0; i<len; i++)
    {
      if (i>0)
	fprintf((FILE*)param, ",");
      fprintf((FILE*)param, "%d", eVect[i]);
    }
  fprintf((FILE*)param, "]");
}
*/

vector<vector<int> > DoMemoryEfficientEnum(TheGroupFormat TheGRP, MyMatrix* FAC, 
					   int LevSearch)
{
  vector<vector<int> > RetList;
  vector<vector<int> > *h;
  h=&RetList;
  Generic_DoMemoryEfficientEnum(TheGRP, FAC, LevSearch, 
				&report_face_vectvectint, (void *)h);
  return RetList;
}


vector<vector<vector<int> > > DoMemoryEfficientAllEnum(TheGroupFormat TheGRP, MyMatrix* FAC, 
						       int LevSearch)
{
  vector<vector<vector<int> > > RetList;
  vector<vector<vector<int> > > *h;
  vector<vector<int> > EmptyList;
  int iLevel;
  EmptyList.clear();
  for (iLevel=0; iLevel<LevSearch; iLevel++)
    RetList.push_back(EmptyList);
  h=&RetList;
  Generic_DoMemoryEfficientAllEnum(TheGRP, FAC, LevSearch, 
				   &report_allface_vectvectint, (void *)h);
  return RetList;
}








vector<vector<vector<int> > > EnumerationFaces(TheGroupFormat TheGRP, MyMatrix* FAC, int LevSearch)
{
  MyMatrix FACred;
  int n;
  WeightMatrix WMat;
  vector<vector<vector<int> > > RetList;
  vector<int> eInv, fOrb, eOrb;
  vector<vector<int> > ListOrb, PrevListOrb;
  list<vector<int> >::iterator iter;
  list<vector<int> > TheSpann;
  int iLevel, iOrb, nbOrb, i;
  vector<vector<int> > vvO, ListInv;
  vector<int> nList, eList;
  TheGroupFormat StabFace;
  n=TheGRP.n;
  MATRIX_ColReduction(FAC, &FACred);
  GetWeightMatrix(WMat, &FACred);
  for (i=0; i<n; i++)
    eList.push_back(i);
  vvO=DecomposeOrbitPoint(TheGRP, eList);
  nbOrb=vvO.size();
  ListOrb.clear();
  for (iOrb=0; iOrb<nbOrb; iOrb++)
    {
      nList.clear();
      nList.push_back(vvO[iOrb][0]);
      ListOrb.push_back(nList);
    }
  RetList.push_back(ListOrb);
  for (iLevel=1; iLevel<=LevSearch; iLevel++)
    {
      PrevListOrb=RetList[iLevel-1];
      nbOrb=PrevListOrb.size();
      ListOrb.clear();
      ListInv.clear();
      for (iOrb=0; iOrb<nbOrb; iOrb++)
	{
	  eOrb=PrevListOrb[iOrb];
	  GetStabilizer(TheGRP, eOrb, StabFace);
	  TheSpann=SPAN_face_LinearProgramming(eOrb, 
					       StabFace, 
					       &FACred, TheGRP);
	  iter=TheSpann.begin();
	  while(iter != TheSpann.end())
	    {
	      fOrb=*iter;
	      eInv=GetLocalInvariantWeightMatrix(WMat, fOrb);
	      GROUP_FuncInsertInSet_UseInv(TheGRP,
					   fOrb, eInv, 
					   ListOrb, ListInv);
	      iter++;
	    }
	}
      RetList.push_back(ListOrb);
    }
  MATRIX_Free(&FACred);
  free(WMat.TheMat);
  return RetList;
}


void PrintListOrb_GAP(FILE *f, vector<vector<int> > ListOrb)
{
  int nbOrb, iOrb, len, i;
  vector<int> eOrb;
  fprintf(f, "return ");
  fprintf(f, "rec(ListRepresentent:=[");
  nbOrb=ListOrb.size();
  for (iOrb=0; iOrb<nbOrb; iOrb++)
    {
      if (iOrb>0)
	fprintf(f, ",");
      fprintf(f, "[");
      eOrb=ListOrb[iOrb];
      len=eOrb.size();
      for (i=0; i<len; i++)
	{
	  if (i>0)
	    fprintf(f, ",");
	  fprintf(f, "%d", eOrb[i]+1);
	}
      fprintf(f, "]");
    }
  fprintf(f, "]);");
}



void PrintListListOrb_IntGAP(FILE *f, vector<vector<vector<int> > > ListListOrb)
{
  int nbLev, iLev, nbOrb, iOrb, len, i;
  vector<int> eOrb;
  nbLev=ListListOrb.size();
  fprintf(f, "return [");
  for (iLev=0; iLev<nbLev; iLev++)
    {
      if (iLev > 0)
	fprintf(f, ",\n");
      fprintf(f, "rec(ListRepresentent:=[");
      nbOrb=ListListOrb[iLev].size();
      for (iOrb=0; iOrb<nbOrb; iOrb++)
	{
	  if (iOrb>0)
	    fprintf(f, ",");
	  fprintf(f, "[");
	  eOrb=ListListOrb[iLev][iOrb];
	  len=eOrb.size();
	  for (i=0; i<len; i++)
	    {
	      if (i>0)
		fprintf(f, ",");
	      fprintf(f, "%d", eOrb[i]+1);
	    }
	  fprintf(f, "]");
	}
      fprintf(f, "])");
    }
  fprintf(f, "];\n");
}
