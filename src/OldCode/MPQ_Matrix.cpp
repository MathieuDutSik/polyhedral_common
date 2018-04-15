#include "MPQ_Matrix.h"
#include "MatrixOperation.h"

#define DEBUG
void VECTOR_Free(MyVector *TheVect)
{
  int n, i;
  n=TheVect->n;
  for (i=0; i<n; i++)
    mpq_clear(TheVect->ListElt[i]);
  free(TheVect->ListElt);
}
void VECTOR_Allocate(MyVector *TheVect, int n)
{
  int i;
  TheVect->n=n;
  if ((TheVect->ListElt = (mpq_t*)malloc(n*sizeof(mpq_t))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<n; i++)
    mpq_init(TheVect->ListElt[i]);
}
void VECTOR_ZeroAssignation(MyVector *TheVect)
{
  int i;
  for (i=0; i<TheVect->n; i++)
    mpq_set_si(TheVect->ListElt[i], 0, 1);
}
void VECTOR_Assign(MyVector *TheVect, int i, mpq_t eVal)
{
  mpq_set(TheVect->ListElt[i], eVal);
}
void VECTOR_Get(MyVector *TheVect, int i, mpq_t eVal)
{
  mpq_set(eVal, TheVect->ListElt[i]);
}
void VECTOR_AllocateCopy(MyVector *eVectI, MyVector *eVectO)
{
  int n, i;
  n=eVectI->n;
  VECTOR_Allocate(eVectO, n);
  for (i=0; i<n; i++)
    mpq_set(eVectO->ListElt[i], eVectI->ListElt[i]);
}
void VECTOR_Print(FILE *f, MyVector *eVect)
{
  int n, i;
  n=eVect->n;
  for (i=0; i<n; i++)
    {
      fprintf(f, " ");
      mpq_out_str(f, 10, eVect->ListElt[i]);
    }
  fprintf(f, "\n");
}
int VECTOR_IsEqual(MyVector *eVect1, MyVector *eVect2)
{
  int n1, n2, i, test;
  n1=eVect1->n;
  n2=eVect2->n;
  if (n1 != n2)
    {
      fprintf(stderr, "Error at n1, n2\n");
      exit(1);
    }
  for (i=0; i<n1; i++)
    {
      test=mpq_cmp(eVect1->ListElt[i], eVect2->ListElt[i]);
      if (test != 0)
	return 0;
    }
  return 1;
}


void MATRIX_ReadFernando(FILE *f, MyMatrix *eMat)
{
  int iEnt, nbEnt, iRow, iCol, nbRow, nbCol;
  mpq_t t_q;
  int ret;
  ret=fscanf(f, "%d %d", &nbRow, &nbCol);
  MATRIX_Allocate(eMat, nbRow, nbCol);
  MATRIX_ZeroAssignation(eMat);
  ret=fscanf(f, "%d", &nbEnt);
  mpq_init(t_q);
  for (iEnt=0; iEnt<nbEnt; iEnt++)
    {
      ret=fscanf(f, "%d", &iRow);
      ret=fscanf(f, "%d", &iCol);
      mpq_inp_str(t_q, f, 10);
      MATRIX_Assign(eMat, iRow, iCol, t_q);
    }
  mpq_clear(t_q);
}


/* We copy eMatI to eMatO */
void MATRIX_Copy(MyMatrix *eMatI, MyMatrix *eMatO)
{
  int nbRowI, nbColI, nbRowO, nbColO, n, i;
  nbRowI=eMatI->nbRow;
  nbRowO=eMatO->nbRow;
  nbColI=eMatI->nbCol;
  nbColO=eMatO->nbCol;
  if (nbRowI != nbRowO || nbColI != nbColO)
    {
      fprintf(stderr, "Error in the input\n");
      exit(1);
    }
  n=nbRowI*nbColI;
  for (i=0; i<n; i++)
    mpq_set(eMatO->ListElt[i], eMatI->ListElt[i]);
}

void MATRIX_AllocateCopy(MyMatrix *eMatI, MyMatrix *eMatO)
{
  int nbRow, nbCol, n, i;
  nbRow=eMatI->nbRow;
  nbCol=eMatI->nbCol;
  MATRIX_Allocate(eMatO, nbRow, nbCol);
  n=nbRow*nbCol;
  for (i=0; i<n; i++)
    mpq_set(eMatO->ListElt[i], eMatI->ListElt[i]);
}




void MATRIX_Free(MyMatrix *TheMat)
{
  int nbRow, nbCol, n, i;
  nbRow=TheMat->nbRow;
  nbCol=TheMat->nbCol;
  n=nbRow*nbCol;
  for (i=0; i<n; i++)
    mpq_clear(TheMat->ListElt[i]);
  free(TheMat->ListElt);
}
void MATRIX_Allocate(MyMatrix *TheMat, int nbRow, int nbCol)
{
  int n, i;
  TheMat->nbRow=nbRow;
  TheMat->nbCol=nbCol;
  n=nbRow*nbCol;
  if ((TheMat->ListElt = (mpq_t*)malloc(n*sizeof(mpq_t))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<n; i++)
    mpq_init(TheMat->ListElt[i]);
}
void MATRIX_ZeroAssignation(MyMatrix *TheMat)
{
  int i, n;
  n=TheMat->nbRow*TheMat->nbCol;
  for (i=0; i<n; i++)
    mpq_set_si(TheMat->ListElt[i], 0, 1);
}


void MATRIX_Assign(MyMatrix *TheMat, int iRow, int iCol, mpq_t eVal)
{
  int nbRow, nbCol, idx;
  nbRow=TheMat->nbRow;
  nbCol=TheMat->nbCol;
  idx=iRow+nbRow*iCol;
#ifdef DEBUG
  if (iRow >= nbRow || iCol >= nbCol)
    {
      fprintf(stderr, "iRow=%d nbRow=%d iCol=%d nbCol=%d\n", iRow, nbRow, iCol, nbCol);
      fprintf(stderr, "Out of bound\n");
      exit(1);
    }
#endif
  mpq_set(TheMat->ListElt[idx], eVal);
}


void MATRIX_Get(MyMatrix *TheMat, int iRow, int iCol, mpq_t eVal)
{
  int nbRow, idx;
  nbRow=TheMat->nbRow;
  idx=iRow+nbRow*iCol;
  mpq_set(eVal, TheMat->ListElt[idx]);
}
void MATRIX_Transpose(MyMatrix *TheMat, MyMatrix *TheTrans)
{
  int nbCol, iCol, nbRow, iRow;
  mpq_t eVal;
  mpq_init(eVal);
  nbCol=TheMat->nbCol;
  nbRow=TheMat->nbRow;
  MATRIX_Allocate(TheTrans, nbCol, nbRow);
  for (iCol=0; iCol<nbCol; iCol++)
    for (iRow=0; iRow<nbRow; iRow++)
    {
      MATRIX_Get(TheMat, iRow, iCol, eVal);
      MATRIX_Assign(TheTrans, iCol, iRow, eVal);
    }
  mpq_clear(eVal);
}
/* We allocate M3 and compute the product M3=M1 M2 */
void MATRIX_Product(MyMatrix *M1, MyMatrix *M2, MyMatrix *TheProd)
{
  int nbRow1, nbCol1, nbRow2, nbCol2;
  mpq_t prov1, prov2, prov3;
  int iCol, iRow, i;
  nbCol1=M1->nbCol;
  nbRow1=M1->nbRow;
  nbCol2=M2->nbCol;
  nbRow2=M2->nbRow;
  if (nbCol1 != nbRow2)
    {
      fprintf(stderr, "Error in matrix sizes\n");
      exit(1);
    }
  MATRIX_Allocate(TheProd, nbRow1, nbCol2);
  mpq_init(prov1);
  mpq_init(prov2);
  mpq_init(prov3);
  for (iCol=0; iCol<nbCol2; iCol++)
    for (iRow=0; iRow<nbRow1; iRow++)
      {
	mpq_set_si(prov1, 0,1);
	for (i=0; i<nbCol1; i++)
	  {
	    MATRIX_Get(M1, iRow, i, prov2);
	    MATRIX_Get(M2, i, iCol, prov3);
	    mpq_mul(prov2, prov2, prov3);
	    mpq_add(prov1, prov1, prov2);
	  }
	MATRIX_Assign(TheProd, iRow, iCol, prov1);
      }
  mpq_clear(prov1);
  mpq_clear(prov2);
  mpq_clear(prov3);
}



// We set TheProd=M1*M2
void VECTOR_MATRIX_Product(MyVector *V1, MyMatrix *M2, MyVector *TheProd)
{
  int nbCol1, nbRow2, nbCol2;
  mpq_t prov1, prov2, prov3;
  int iCol, i;
  nbCol1=V1->n;
  nbCol2=M2->nbCol;
  nbRow2=M2->nbRow;
  if (nbCol1 != nbRow2)
    {
      fprintf(stderr, "Error in matrix sizes\n");
      exit(1);
    }
  mpq_init(prov1);
  mpq_init(prov2);
  mpq_init(prov3);
  for (iCol=0; iCol<nbCol2; iCol++)
    {
      mpq_set_si(prov1, 0,1);
      for (i=0; i<nbCol1; i++)
	{
	  VECTOR_Get(V1, i, prov2);
	  MATRIX_Get(M2, i, iCol, prov3);
	  mpq_mul(prov2, prov2, prov3);
	  mpq_add(prov1, prov1, prov2);
	}
      VECTOR_Assign(TheProd, iCol, prov1);
    }
  mpq_clear(prov1);
  mpq_clear(prov2);
  mpq_clear(prov3);
}





void MATRIX_Inverse_destroy(MyMatrix *Input, MyMatrix *Output)
{
  int nbRow, nbCol, idx;
  int iCol, iRow, WeFound;
  int idx1, idx2, iRowFound;
  int idxMult, iColB, test, idxFound;
  nbRow=Input->nbRow;
  nbCol=Input->nbCol;
  mpq_t prov1, prov2;
  mpq_init(prov1);
  mpq_init(prov2);
  if (nbRow != nbCol)
    {
      fprintf(stderr, "Error on nbRow, nbCol in MATRIX_Inverse_destroy\n");
      exit(1);
    }
  for (iRow=0; iRow<nbRow; iRow++)
    for (iCol=0; iCol<nbRow; iCol++)
      {
	if (iRow == iCol)
	  mpq_set_si(prov1, 1, 1);
	else
	  mpq_set_si(prov1, 0, 1);
	idx=iRow+nbRow*iCol;
	mpq_set(Output->ListElt[idx], prov1);
      }
  iRowFound=-400;
  for (iCol=0; iCol<nbCol; iCol++)
    {
      /*      fprintf(stderr, "iCol=%d/%d\n", iCol, nbCol);*/
      WeFound=0;
      for (iRow=iCol; iRow<nbRow; iRow++)
	if (WeFound == 0)
	  {
	    idx=iRow+nbRow*iCol;
	    test=mpq_sgn(Input->ListElt[idx]);
	    if (test != 0)
	      {
		WeFound=1;
		iRowFound=iRow;
		mpq_inv(prov1, Input->ListElt[idx]);
	      }
	  }
      /*fprintf(stderr, "WeFound=%d iRowFound=%d\n", WeFound, iRowFound);*/
      for (iColB=0; iColB<nbCol; iColB++)
	{
	  idx=iRowFound+nbRow*iColB;
	  mpq_mul(Input->ListElt[idx], prov1, Input->ListElt[idx]);
	  mpq_mul(Output->ListElt[idx], prov1, Output->ListElt[idx]);
	}
      for (iRow=0; iRow<nbRow; iRow++)
	if (iRow != iRowFound)
	  {
	    idxMult=iRow+nbRow*iCol;
	    mpq_set(prov2, Input->ListElt[idxMult]);
	    for (iColB=0; iColB<nbCol; iColB++)
	      {
		idxFound=iRowFound+nbRow*iColB;
		idx=iRow+nbRow*iColB;
		mpq_mul(prov1, prov2, Input->ListElt[idxFound]);
		mpq_sub(Input->ListElt[idx], Input->ListElt[idx], prov1);
		mpq_mul(prov1, prov2, Output->ListElt[idxFound]);
		mpq_sub(Output->ListElt[idx], Output->ListElt[idx], prov1);
	      }
	  }
      if (iRowFound != iCol)
	{
	  for (iColB=0; iColB<nbCol; iColB++)
	    {
	      idx1=iRowFound+nbRow*iColB;
	      idx2=iCol+nbRow*iColB;
	      mpq_swap(Input->ListElt[idx1], Input->ListElt[idx2]);
	      mpq_swap(Output->ListElt[idx1], Output->ListElt[idx2]);
	    }
	}
    }
  mpq_clear(prov1);
  mpq_clear(prov2);
}

void MATRIX_Inverse(MyMatrix *Input, MyMatrix *Output)
{
  int nbRow, nbCol;
  MyMatrix *provMat;
  nbRow=Input->nbRow;
  nbCol=Input->nbCol;
  if ((provMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  MATRIX_Allocate(provMat, nbRow, nbCol);
  MATRIX_Copy(Input, provMat);
  MATRIX_Inverse_destroy(provMat, Output);
  MATRIX_Free(provMat);
  free(provMat);
}

void SELECT_Free(MySelection *eSelect)
{
  free(eSelect->ListColSelect);
  free(eSelect->ListRowSelect);
  free(eSelect->ListColSelect01);
  free(eSelect->ListRowSelect01);
}


void MATRIX_SelectColRow(MyMatrix *Input, MySelection *eSelect, MyMatrix *NSP)
{
  int nbRow, nbCol;
  MyMatrix provMat;
  mpq_t eVal, eVal1, eVal2, eVal3;
  int iRank, test;
  int sizMat, nbVect, iRow, iCol;
  int eCol, IsFinished;
  int nbVectZero, maxRank, eRank, FirstNonZeroCol;
  mpq_init(eVal);
  mpq_init(eVal1);
  mpq_init(eVal2);
  mpq_init(eVal3);
  nbRow=Input->nbRow;
  nbCol=Input->nbCol;
  maxRank=nbRow;
  if (nbCol < maxRank)
    maxRank=nbCol;
  sizMat=maxRank+1;
  MATRIX_Allocate(&provMat, sizMat, nbCol);
  eSelect->nbRow=nbRow;
  eSelect->nbCol=nbCol;
  if ((eSelect->ListColSelect = (int*)malloc(nbCol*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eSelect->ListRowSelect = (int*)malloc(nbRow*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eSelect->ListColSelect01 = (int*)malloc(nbCol*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eSelect->ListRowSelect01 = (int*)malloc(nbRow*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (iCol=0; iCol<nbCol; iCol++)
    eSelect->ListColSelect[iCol]=0;
  for (iRow=0; iRow<nbRow; iRow++)
    eSelect->ListRowSelect[iRow]=0;
  for (iCol=0; iCol<nbCol; iCol++)
    eSelect->ListColSelect01[iCol]=0;
  for (iRow=0; iRow<nbRow; iRow++)
    eSelect->ListRowSelect01[iRow]=0;
  eRank=0;
  for (iRow=0; iRow<nbRow; iRow++)
    {
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(Input, iRow, iCol, eVal);
	  MATRIX_Assign(&provMat, eRank, iCol, eVal);
	}
      for (iRank=0; iRank<eRank; iRank++)
	{
	  eCol=eSelect->ListColSelect[iRank];
	  MATRIX_Get(&provMat, eRank, eCol, eVal1);
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      MATRIX_Get(&provMat, iRank, iCol, eVal2);
	      MATRIX_Get(&provMat, eRank, iCol, eVal3);
	      mpq_mul(eVal, eVal1, eVal2);
	      mpq_sub(eVal3, eVal3, eVal);
	      MATRIX_Assign(&provMat, eRank, iCol, eVal3);
	    }
	}
      IsFinished=1;
      FirstNonZeroCol=-1;
      for (iCol=0; iCol<nbCol; iCol++)
	if (IsFinished == 1)
	  {
	    MATRIX_Get(&provMat, eRank, iCol, eVal);
	    test=mpq_sgn(eVal);
	    if (test != 0)
	      {
		FirstNonZeroCol=iCol;
		IsFinished=0;
	      }
	  }
      if (IsFinished == 0)
	{
	  eSelect->ListColSelect[eRank]=FirstNonZeroCol;
	  eSelect->ListRowSelect[eRank]=iRow;
	  eSelect->ListColSelect01[FirstNonZeroCol]=1;
	  eSelect->ListRowSelect01[iRow]=1;
	  MATRIX_Get(&provMat, eRank, FirstNonZeroCol, eVal);
	  mpq_inv(eVal2, eVal);
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      MATRIX_Get(&provMat, eRank, iCol, eVal);
	      mpq_mul(eVal, eVal, eVal2);
	      MATRIX_Assign(&provMat, eRank, iCol, eVal);
	    }
	  for (iRank=0; iRank<eRank; iRank++)
	    {
	      MATRIX_Get(&provMat, iRank, FirstNonZeroCol, eVal1);
	      for (iCol=0; iCol<nbCol; iCol++)
		{
		  MATRIX_Get(&provMat, iRank, iCol, eVal2);
		  MATRIX_Get(&provMat, eRank, iCol, eVal3);
		  mpq_mul(eVal, eVal1, eVal3);
		  mpq_sub(eVal2, eVal2, eVal);
		  MATRIX_Assign(&provMat, iRank, iCol, eVal2);
		}
	    }
	  eRank++;
	}
    }
  eSelect->TheRank=eRank;
  nbVectZero=nbCol-eRank;
  MATRIX_Allocate(NSP, nbVectZero, nbCol);
  MATRIX_ZeroAssignation(NSP);
  nbVect=0;
  for (iCol=0; iCol<nbCol; iCol++)
    if (eSelect->ListColSelect01[iCol] == 0)
      {
	mpq_set_si(eVal, 1, 1);
	MATRIX_Assign(NSP, nbVect, iCol, eVal);
	for (iRank=0; iRank<eRank; iRank++)
	  {
	    eCol=eSelect->ListColSelect[iRank];
	    MATRIX_Get(&provMat, iRank, iCol, eVal);
	    mpq_neg(eVal, eVal);
	    MATRIX_Assign(NSP, nbVect, eCol, eVal);
	  }
	nbVect++;
      }
  mpq_clear(eVal);
  mpq_clear(eVal1);
  mpq_clear(eVal2);
  mpq_clear(eVal3);
  MATRIX_Free(&provMat);
}

void MATRIX_Print(FILE *f, MyMatrix *TheMat)
{
  int iRow, nbRow, iCol, nbCol, idx;
  nbRow=TheMat->nbRow;
  nbCol=TheMat->nbCol;
  fprintf(f, "nbRow=%d nbCol=%d\n", nbRow, nbCol);
  for (iRow=0; iRow<nbRow; iRow++)
    {
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  fprintf(f, " ");
	  idx=iRow+nbRow*iCol;
	  mpq_out_str(f, 10, TheMat->ListElt[idx]);
	}
      fprintf(f, "\n");
    }
}

// In input a matrix
// in output the lines that spann its rank in sequential order
void MATRIX_RowReduction(MyMatrix *eMatIn, MyMatrix *eMatOut)
{
  MyMatrix NSP;
  MySelection eSelect;
  int nbCol, eRank, iLine, jRow, iCol;
  mpq_t eVal;
  mpq_init(eVal);
  MATRIX_SelectColRow(eMatIn, &eSelect, &NSP);
  nbCol=eMatIn->nbCol;
  eRank=eSelect.TheRank;
  MATRIX_Allocate(eMatOut, eRank, nbCol);
  for (iLine=0; iLine<eRank; iLine++)
    {
      jRow=eSelect.ListRowSelect[iLine];
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(eMatIn, jRow, iCol, eVal);
	  MATRIX_Assign(eMatOut, iLine, iCol, eVal);
	}
    }
  SELECT_Free(&eSelect);
  MATRIX_Free(&NSP);
  mpq_clear(eVal);
}


void MATRIX_ColReduction(MyMatrix *eMatIn, MyMatrix *eMatOut)
{
  MyMatrix NSP;
  MySelection eSelect;
  int nbRow, eRank, iCol, jCol, iRow;
  mpq_t eVal;
  mpq_init(eVal);
  MATRIX_SelectColRow(eMatIn, &eSelect, &NSP);
  nbRow=eMatIn->nbRow;
  eRank=eSelect.TheRank;
  MATRIX_Allocate(eMatOut, nbRow, eRank);
  for (iCol=0; iCol<eRank; iCol++)
    {
      jCol=eSelect.ListColSelect[iCol];
      for (iRow=0; iRow<nbRow; iRow++)
	{
	  MATRIX_Get(eMatIn, iRow, jCol, eVal);
	  MATRIX_Assign(eMatOut, iRow, iCol, eVal);
	}
    }
  SELECT_Free(&eSelect);
  MATRIX_Free(&NSP);
  mpq_clear(eVal);
}





void MATRIX_GetIdentity(MyMatrix *TheMat, int n)
{
  int i, j;
  mpq_t t;
  mpq_init(t);
  MATRIX_Allocate(TheMat, n, n);
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      {
	if (i == j)
	  mpq_set_si(t, 1,1);
	else
	  mpq_set_si(t, 0,1);
	MATRIX_Assign(TheMat, i, j, t);
      }
  mpq_clear(t);
}


void MPQ_ImageIntVector(MyVector *eVect, MyMatrix *TheMat, MyVector *eVectImg)
{
  mpq_t t, prov1, prov2;
  int iCol, nbCol, iRow, nbRow, n;
  n=eVect->n;
  nbRow=TheMat->nbRow;
  nbCol=TheMat->nbCol;
  if (n != nbRow)
    {
      fprintf(stderr, "Error in ImageIntVector\n");
      fprintf(stderr, "n=%d nbRow=%d\n", n, nbRow);
      exit(1);
    }
  mpq_init(t);
  mpq_init(prov1);
  mpq_init(prov2);
  for (iCol=0; iCol<nbCol; iCol++)
    {
      mpq_set_si(t, 0,1);
      for (iRow=0; iRow<n; iRow++)
	{
	  MATRIX_Get(TheMat, iRow, iCol, prov1);
	  mpq_mul(prov2, prov1, eVect->ListElt[iRow]);
	  mpq_add(t, t, prov2);
	}
      mpq_set(eVectImg->ListElt[iCol], t);
    }
  mpq_clear(t);
  mpq_clear(prov1);
  mpq_clear(prov2);
}

void MATRIX_Concatenate(MyMatrix *eMat1, MyMatrix *eMat2, MyMatrix *eMatRet)
{
  int nbRow1, nbRow2, nbCol1, nbCol2;
  int nbCol, iRow, iCol;
  mpq_t eVal;
  mpq_init(eVal);
  nbRow1=eMat1->nbRow;
  nbRow2=eMat2->nbRow;
  nbCol1=eMat1->nbCol;
  nbCol2=eMat2->nbCol;
  if (nbCol1 != nbCol2)
    {
      fprintf(stderr, "nbCol1=%d nbCol2=%d\n", nbCol1, nbCol2);
      fprintf(stderr, "Error in the number of columns\n");
      exit(1);
    }
  nbCol=nbCol1;
  MATRIX_Allocate(eMatRet, nbRow1+nbRow2, nbCol);
  for (iRow=0; iRow<nbRow1; iRow++)
    {
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(eMat1, iRow, iCol, eVal);
	  MATRIX_Assign(eMatRet, iRow, iCol, eVal);
	}
    }
  for (iRow=0; iRow<nbRow2; iRow++)
    {
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(eMat2, iRow, iCol, eVal);
	  MATRIX_Assign(eMatRet, iRow+nbRow1, iCol, eVal);
	}
    }
  mpq_clear(eVal);
}


void MATRIX_CongrMap(MyMatrix *eMat, MyMatrix *RetMat)
{
  MyMatrix *eMatCopy, *TheInv;
  int nbRow, nbCol;
  nbRow=eMat->nbRow;
  nbCol=eMat->nbCol;
  if ((eMatCopy = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheInv = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if (nbRow != nbCol)
    {
      fprintf(stderr, "The matrix is not square, PANIC\n");
      exit(1);
    }
  MATRIX_Allocate(eMatCopy, nbRow, nbRow);
  MATRIX_Allocate(TheInv, nbRow, nbRow);
  MATRIX_Copy(eMat, eMatCopy);
  MATRIX_Inverse_destroy(eMatCopy, TheInv);
  MATRIX_Transpose(TheInv, RetMat);
  MATRIX_Free(eMatCopy);
  MATRIX_Free(TheInv);
  free(eMatCopy);
  free(TheInv);
}


void MATRIX_ConvertToInt(MyMatrix *eMat, MatINT &RetMat)
{
  mpq_t eVal_q;
  mpz_t eNum_z, eDen_z;
  int i, j, n, p;
  int eNum, eDen;
  vector<int> eLine;
  mpq_init(eVal_q);
  mpz_init(eNum_z);
  mpz_init(eDen_z);
  RetMat.clear();
  n=eMat->nbRow;
  p=eMat->nbCol;
  for (i=0; i<n; i++)
    {
      for (j=0; j<p; j++)
	{
	  MATRIX_Get(eMat, i, j, eVal_q);
	  mpq_get_num(eNum_z, eVal_q);
	  mpq_get_den(eDen_z, eVal_q);
	  eNum=mpz_get_si(eNum_z);
	  eDen=mpz_get_si(eDen_z);
	  if (eDen != 1)
	    {
	      fprintf(stderr, "The entry is not integral, die!\n");
	      exit(1);
	    }
	  eLine.push_back(eNum);
	}
      RetMat.push_back(eLine);
      eLine.clear();
    }
  mpq_clear(eVal_q);
  mpz_clear(eNum_z);
  mpz_clear(eDen_z);
}


