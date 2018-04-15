#if !defined MATRIX_INCLUDE_SYMPOL
#include "MatrixOperation.h"
#include <stdlib.h>
#include "gmp.h"
typedef struct {
  int n;
  mpq_t *ListElt;
} MyVector;

typedef MyVector* MyVectorPtr;

typedef struct {
  int nbRow;
  int nbCol;
  mpq_t *ListElt;
} MyMatrix;

typedef MyMatrix* MyMatrixPtr;

typedef struct {
  int nbCol;
  int nbRow;
  int TheRank;
  int *ListColSelect;
  int *ListRowSelect;
  int *ListColSelect01;
  int *ListRowSelect01;
} MySelection;


void VECTOR_Free(MyVector *TheVect);
void VECTOR_Allocate(MyVector *TheVect, int n);
void VECTOR_AllocateCopy(MyVector *eVectI, MyVector *eVectO);
void VECTOR_ZeroAssignation(MyVector *TheVect);
void VECTOR_Assign(MyVector *TheVect, int i, mpq_t eVal);
void VECTOR_Get(MyVector *TheVect, int i, mpq_t eVal);
void VECTOR_Print(FILE *f, MyVector *eVect);
void VECTOR_MATRIX_Product(MyVector *M1, MyMatrix *M2, MyVector *TheProd);
void SELECT_Free(MySelection *eSelect);
void MATRIX_Copy(MyMatrix *eMatI, MyMatrix *eMatO);
void MATRIX_AllocateCopy(MyMatrix *eMatI, MyMatrix *eMatO);
void MATRIX_Allocate(MyMatrix *TheMat, int nbRow, int nbCol);
void MATRIX_Free(MyMatrix *TheMat);
void MATRIX_ZeroAssignation(MyMatrix *TheMat);
void MATRIX_Assign(MyMatrix *TheMat, int iRow, int iCol, mpq_t eVal);
void MATRIX_Get(MyMatrix *TheMat, int iRow, int iCol, mpq_t eVal);
void MATRIX_Inverse_destroy(MyMatrix *Intput, MyMatrix *Output);
void MATRIX_Inverse(MyMatrix *Input, MyMatrix *Output);
void MATRIX_Print(FILE *f, MyMatrix *TheMat);
void MATRIX_Transpose(MyMatrix *TheMat, MyMatrix *TheTrans);
void MATRIX_Product(MyMatrix *M1, MyMatrix *M2, MyMatrix *TheProd);
void MATRIX_SelectColRow(MyMatrix *Input, MySelection *eSelect, MyMatrix *NSP);
void MATRIX_CongrMap(MyMatrix *eMat, MyMatrix *RetMat);
void MATRIX_ConvertToInt(MyMatrix *eMat, MatINT &RetMat);
void MPQ_ImageIntVector(MyVector *eVect, MyMatrix *TheMat, MyVector *eVectImg);
int VECTOR_IsEqual(MyVector *eVect1, MyVector *eVect2);
void MATRIX_GetIdentity(MyMatrix *TheMat, int n);
void MATRIX_ReadFernando(FILE *f, MyMatrix *eMat);
void MATRIX_RowReduction(MyMatrix *eMatIn, MyMatrix *eMatOut);
void MATRIX_ColReduction(MyMatrix *eMatIn, MyMatrix *eMatOut);
void MATRIX_Concatenate(MyMatrix *eMat1, MyMatrix *eMat2, MyMatrix *eMatRet);
#endif
#define MATRIX_INCLUDE_SYMPOL
