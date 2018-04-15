#include "MPQ_CDD_LRS.h"
void LPSOL_Free(LpSolution *eSol)
{
  int i, nbDirect, nbDual;
  nbDirect=eSol->nbDirect;
  nbDual=eSol->nbDual;
  for (i=0; i<nbDual; i++)
    mpq_clear(eSol->DualSolutionVal[i]);
  for (i=0; i<nbDirect; i++)
    mpq_clear(eSol->DirectSolutionVal[i]);
  mpq_clear(eSol->OptimalValue);
  free(eSol->DualSolutionPos);
  free(eSol->DualSolutionVal);
  free(eSol->DirectSolutionPos);
  free(eSol->DirectSolutionVal);
}


void LPSOL_Allocate(LpSolution *eSol, int nbDirect, int nbDual)
{
  int i;
  eSol->nbDirect=nbDirect;
  eSol->nbDual=nbDual;
  if ((eSol->DirectSolutionPos = (int*)malloc(nbDirect*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eSol->DualSolutionPos = (int*)malloc(nbDual*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eSol->DirectSolutionVal = (mpq_t*)malloc(nbDirect*sizeof(mpq_t))) == 0)
    exit (EXIT_FAILURE);
  if ((eSol->DualSolutionVal = (mpq_t*)malloc(nbDual*sizeof(mpq_t))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<nbDirect; i++)
    {
      mpq_init(eSol->DirectSolutionVal[i]);
      mpq_set_si(eSol->DirectSolutionVal[i], 0, 1);
      eSol->DirectSolutionPos[i]=0;
    }
  for (i=0; i<nbDual; i++)
    {
      mpq_init(eSol->DualSolutionVal[i]);
      mpq_set_si(eSol->DualSolutionVal[i], 0, 1);
      eSol->DualSolutionPos[i]=0;
    }
  mpq_init(eSol->OptimalValue);
}


dd_MatrixPtr MyMatrix_PolyFile2Matrix (MyMatrix *TheEXT)
{
  dd_MatrixPtr M=NULL;
  dd_rowrange m_input, i;
  dd_colrange d_input, j;
  dd_RepresentationType rep;
  mpq_t value;
  dd_NumberType NT;
  mpq_init(value);
  m_input=TheEXT->nbRow;
  d_input=TheEXT->nbCol;

  /*  NT=dd_GetNumberType("integer");*/
  NT=dd_Rational;
  rep=dd_Generator; /* using dd_Inequality led to horrible bugs */
  /*  NT=dd_dd_Integer;*/
  M=dd_CreateMatrix(m_input, d_input);
  M->representation=rep;
  M->numbtype=NT;

  for (i = 1; i <= m_input; i++)
    for (j = 1; j <= d_input; j++)
      {
	MATRIX_Get(TheEXT, i-1, j-1, value);
	mpq_set(M->matrix[i-1][j - 1],value);
      }
  mpq_clear(value);
  return M;
}



void CDD_LinearProgramming(MyMatrix *TheEXT, MyVector *eVect, LpSolution *eSol)
{
  dd_ErrorType error=dd_NoError;
  dd_MatrixPtr M;
  dd_LPSolverType solver=dd_DualSimplex;  /* either DualSimplex or CrissCross */
  dd_LPPtr lp;   /* pointer to LP data structure that is not visible by user. */
  dd_colrange j;
  mpq_class tcla;
  int d_input;
  int nbRow, nbCol, idx;
  vector<int> DualSolutionPos;
  vector<mpq_class> DualSolutionVal;
  vector<int> DirectSolutionPos;
  vector<mpq_class> DirectSolutionVal;
  M=MyMatrix_PolyFile2Matrix(TheEXT);
  M->representation=dd_Inequality;
  d_input=TheEXT->nbCol;
  for (j = 1; j <= d_input; j++)
    mpq_set(M->rowvec[j - 1], eVect->ListElt[j-1]);
  lp=dd_Matrix2LP(M, &error);
  lp->objective=dd_LPmin;
  nbRow=TheEXT->nbRow;
  nbCol=TheEXT->nbCol;
  dd_LPSolve(lp,solver,&error);
  LPSOL_Allocate(eSol, nbCol+1, nbRow+1);
 switch (lp->LPS){
 case dd_Optimal:
   eSol->PrimalDefined=1;
   eSol->DualDefined=1;
   break;
 case dd_Inconsistent:
   eSol->PrimalDefined=0;
   eSol->DualDefined=1;
   break;
 case dd_DualInconsistent:
   eSol->PrimalDefined=1;
   eSol->DualDefined=0;
   break;
 case dd_StrucDualInconsistent:
   eSol->PrimalDefined=1;
   eSol->DualDefined=0;
   break;
 default:
   eSol->PrimalDefined=0;
   eSol->DualDefined=0;
   break;
 }
 if (eSol->PrimalDefined)
   {
     for (j=1; j<lp->d; j++)
       {
	 eSol->DirectSolutionPos[j]=1;
	 mpq_set(eSol->DirectSolutionVal[j], lp->sol[j]);
       }
   }
 if (eSol->DualDefined)
   {
     for (j=1; j<lp->d; j++)
       {
	 idx=lp->nbindex[j+1];
	 if (idx>0)
	   {
	     eSol->DualSolutionPos[idx]=1;
	     mpq_set(eSol->DualSolutionVal[idx], lp->dsol[j]);
	   }
       }
   }
 if (eSol->PrimalDefined && eSol->DualDefined)
   {
     mpq_set(eSol->OptimalValue, lp->optvalue);
   }
  dd_FreeMatrix(M);
  dd_FreeLPData(lp);
}

void CDD_DualDescription(MyMatrix *TheEXT, MyMatrix *TheFAC)
{
  dd_PolyhedraPtr poly;
  dd_MatrixPtr M;
  int iCol, iRay, nbRay;
  int nbCol;
  dd_RayPtr RayPtr;
  dd_ErrorType err;
  mpq_t prov;
  nbCol=TheEXT->nbCol;
  M=MyMatrix_PolyFile2Matrix(TheEXT);
  poly=dd_DDMatrix2Poly(M, &err);

  RayPtr = poly->child->FirstRay;
  nbRay=0;
  while (RayPtr != NULL)
    {
      if (RayPtr->feasible)
	nbRay++;
      RayPtr = RayPtr->Next;
    }
  fprintf(stderr, "nbRay=%d\n", nbRay);
  MATRIX_Allocate(TheFAC, nbRay, nbCol);
  iRay=0;
  RayPtr = poly->child->FirstRay;
  mpq_init(prov);
  while (RayPtr != NULL)
    {
      if (RayPtr->feasible)
	{
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      dd_set(prov,RayPtr->Ray[iCol]);
	      MATRIX_Assign(TheFAC, iRay, iCol, prov);
	    }
	  iRay++;
	}
      RayPtr = RayPtr->Next;
    }
  mpq_clear(prov);
  dd_FreePolyhedra(poly);
  dd_FreeMatrix(M);
}

#ifdef LRS
void fillModelLRS(MyMatrix *EXT, lrs_dic *P, lrs_dat *Q)
{
  int j;
  int iRow, nbRow, nbCol;
  mpq_t eVal;
  long n;
  long ineq;
  mpq_init(eVal);
  nbRow=EXT->nbRow;
  nbCol=EXT->nbCol;
  n=nbCol;
  
  lrs_mp_vector num = lrs_alloc_mp_vector (n);
  lrs_mp_vector den = lrs_alloc_mp_vector (n);
  
  ineq=1;
  for (iRow=0; iRow<nbRow; iRow++)
    {
      for (j=0; j<nbCol; ++j)
	{
	  MATRIX_Get(EXT, iRow, j, eVal);
	  mpq_get_num(num[j], eVal);
	  mpq_get_den(den[j], eVal);
	}
      /* NOTE lrs_set_row_mp index i starts at 1 */
      lrs_set_row_mp(P, Q, iRow+1, num, den, ineq);
    }
  lrs_clear_mp_vector (num, n);
  lrs_clear_mp_vector (den, n);
  mpq_clear(eVal);
}

void initLRS(MyMatrix *EXT, lrs_dic* & P, lrs_dat* & Q)
{
  lrs_mp_matrix Lin;
  int nbRow, nbCol;
  /* allocate and init structure for static problem data */
  Q = lrs_alloc_dat ("LRS globals");
  if (Q == NULL)
    {
      fprintf(stderr, "Failed rs_alloc_dat\n");
      exit(1);
    }
  nbRow=EXT->nbRow;
  nbCol=EXT->nbCol;
  Q->n = nbCol;
  Q->m = nbRow;
    
  P = lrs_alloc_dic (Q);
  if (P == NULL)
    {
      fprintf(stderr, "We failed allocation, let's die\n");
      exit(1);
    }
  fillModelLRS(EXT, P, Q);
  
  if (!lrs_getfirstbasis (&P, Q, &Lin, FALSE))
    {
      fprintf(stderr, "Another way to die, there are many\n");
      exit(1);
    }
}


void LRS_DualDescription(MyMatrix *TheEXT, MyMatrix *TheFAC)
{
  lrs_dic *P;   /* structure for holding current dictionary and indices  */
  lrs_dat *Q;   /* structure for holding static problem data             */
  lrs_mp_vector output; /* one line of output:ray,vertex,facet,linearity */
  vector<MyVectorPtr> TheList;
  MyVectorPtr eVectPtrZero = NULL;
  int iCol, jCol, nbRow, nbCol, iFac, nbFac;
  mpq_t eVal;
  int col;
  initLRS(TheEXT, P, Q);
  nbRow=TheEXT->nbRow;
  nbCol=TheEXT->nbCol;
  output = lrs_alloc_mp_vector (Q->n);
  nbFac=0;
  mpq_init(eVal);
  do
    {
      for (col = 0; col < nbCol; col++)
	if (lrs_getsolution (P, Q, output, col))
	  {
	    TheList.push_back(eVectPtrZero);
	    nbFac++;
	    TheList[nbFac-1]=(MyVector*)malloc(sizeof(MyVector));
	    TheList[nbFac-1]->n=nbCol;
	    TheList[nbFac-1]->ListElt=(mpq_t*)malloc(nbCol*sizeof(mpq_t));
	    for (jCol=0; jCol<nbCol; jCol++)
	      mpq_init(TheList[nbFac-1]->ListElt[jCol]);
	    for (jCol=0; jCol<nbCol; jCol++)
	      {
		mpq_set_z(eVal, output[jCol]);
		mpq_set(TheList[nbFac-1]->ListElt[jCol], eVal);
	      }
	  }
    }
  while (lrs_getnextbasis (&P, Q, FALSE));
  MATRIX_Allocate(TheFAC, nbFac, nbCol);
  for (iFac=0; iFac<nbFac; iFac++)
    for (iCol=0; iCol<nbCol; iCol++)
      MATRIX_Assign(TheFAC, iFac, iCol, TheList[iFac]->ListElt[iCol]);
  for (iFac=0; iFac<nbFac; iFac++)
    {
      for (iCol=0; iCol<nbCol; iCol++)
	mpq_clear(TheList[iFac]->ListElt[iCol]);
      free(TheList[iFac]->ListElt);
      free(TheList[iFac]);
    }
  lrs_clear_mp_vector (output, Q->n);
  lrs_free_dic (P,Q);
  lrs_free_dat (Q);
  mpq_clear(eVal);
}
#endif

void Polytopization(MyMatrix *EXT, MyMatrix *EXTret)
{
  MyMatrix *nMat, *eBasis;
  MyVector *eVect;
  LpSolution *eSol;
  mpq_t eVal, eSum, prov;
  int nbRow, nbCol, iRow, iCol, iRowWrite;
  int iColSelect, test;
  nbRow=EXT->nbRow;
  nbCol=EXT->nbCol;
  if ((nMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((eBasis = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((eVect = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  if ((eSol = (LpSolution*)malloc(sizeof(LpSolution))) == 0)
    exit (EXIT_FAILURE);
  MATRIX_Allocate(nMat, nbRow, nbCol+1);
  MATRIX_Allocate(eBasis, nbCol, nbCol);
  VECTOR_Allocate(eVect, nbCol+1);
  mpq_init(prov);
  mpq_init(eSum);
  mpq_init(eVal);
  for (iRow=0; iRow<nbRow; iRow++)
    {
      mpq_set_si(eVal, -1,1);
      MATRIX_Assign(nMat, iRow, 0, eVal);
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(EXT, iRow, iCol, eVal);
	  MATRIX_Assign(nMat, iRow, iCol+1, eVal);
	}
    }
  for (iCol=0; iCol<=nbCol; iCol++)
    {
      mpq_set_si(eSum,0,1);
      for (iRow=0; iRow<nbRow; iRow++)
	{
	  MATRIX_Get(nMat, iRow, iCol, eVal);
	  mpq_add(eSum, eSum, eVal);
	}
      VECTOR_Assign(eVect, iCol, eSum);
    }
  CDD_LinearProgramming(nMat, eVect, eSol);
  MATRIX_ZeroAssignation(eBasis);
  for (iCol=0; iCol<nbCol; iCol++)
    {
      MATRIX_Assign(eBasis, iCol, 0, eSol->DirectSolutionVal[iCol+1]);
    }
  iColSelect=-1;
  for (iCol=0; iCol<nbCol; iCol++)
    if (iColSelect == -1)
      {
	MATRIX_Get(eBasis, iCol, 0, eVal);
	test=mpq_sgn(eVal);
	if (test != 0)
	  iColSelect=iCol;
      }
  if (iColSelect == -1)
    {
      fprintf(stderr, "Apparently, we did not find the column\n");
      fprintf(stderr, "That is reason to panic\n");
      exit(1);
    }
  iRowWrite=0;
  mpq_set_si(eVal, 1, 1);
  for (iCol=0; iCol<nbCol; iCol++)
    if (iCol != iColSelect)
      {
	iRowWrite++;
	MATRIX_Assign(eBasis, iCol, iRowWrite, eVal);
      }
  MATRIX_Product(EXT, eBasis, EXTret);
  for (iRow=0; iRow<nbRow; iRow++)
    {
      MATRIX_Get(EXTret, iRow, 0, eVal);
      mpq_inv(eVal, eVal);
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(EXTret, iRow, iCol, prov);
	  mpq_mul(prov, prov, eVal);
	  MATRIX_Assign(EXTret, iRow, iCol, prov);
	}
    }
  mpq_clear(prov);
  mpq_clear(eSum);
  mpq_clear(eVal);
  MATRIX_Free(nMat);
  MATRIX_Free(eBasis);
  VECTOR_Free(eVect);
  LPSOL_Free(eSol);
  free(nMat);
  free(eBasis);
  free(eVect);
  free(eSol);
}



void Initialize_All()
{
#ifdef LRS
  FILE *ms_fIn = NULL;
  FILE *ms_fOut = NULL;
  ms_fIn = fopen("/dev/null","r");
  ms_fOut = fopen("/dev/null","w");
  if (!lrs_mp_init(0, ms_fIn, ms_fOut))
    {
      fprintf(stderr, "We failed to initialize\n");
      exit(1);
    }
#endif
  dd_set_global_constants();
}

void Free_All()
{
  dd_free_global_constants();
#ifdef LRS
  lrs_mp_close();
#endif
}


void FindOneInitialVertexEXT(MyMatrix *EXT, vector<int> &eInc)
{
  int a, b, eVal_i;
  mpq_t eVal, eSum, prov;
  LpSolution *eSol;
  MyMatrix *RnkMat, *nMat, *NSP;
  MyVector *eVect, *TheVert;
  MySelection *eSelect;
  vector<int>::iterator iter;
  int nbRow, nbCol, iCol, iRow, siz, idx, test;
  int TheRank;
  nbRow=EXT->nbRow;
  nbCol=EXT->nbCol;
  mpq_init(eSum);
  mpq_init(eVal);
  mpq_init(prov);
  if ((RnkMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((nMat = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((NSP = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((eVect = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  if ((TheVert = (MyVector*)malloc(sizeof(MyVector))) == 0)
    exit (EXIT_FAILURE);
  if ((eSelect = (MySelection*)malloc(sizeof(MySelection))) == 0)
    exit (EXIT_FAILURE);
  if ((eSol = (LpSolution*)malloc(sizeof(LpSolution))) == 0)
    exit (EXIT_FAILURE);
  VECTOR_Allocate(eVect, nbCol);
  VECTOR_Allocate(TheVert, nbCol);
  for (iCol=0; iCol<nbCol; iCol++)
    {
      mpq_set_si(eSum, 0, 1);
      for (iRow=0; iRow<nbRow; iRow++)
	{
	  MATRIX_Get(EXT, iRow, iCol, eVal);
	  mpq_add(eSum, eSum, eVal);
	}
      mpq_set_si(eVal, 1, nbRow);
      mpq_mul(eSum, eSum, eVal);
      VECTOR_Assign(eVect, iCol, eSum);
    }
  MATRIX_Allocate(nMat, nbRow, nbCol);
  for (iRow=0; iRow<nbRow; iRow++)
    {
      MATRIX_Get(EXT, iRow, 0, eVal);
      MATRIX_Assign(nMat, iRow, 0, eVal);
      for (iCol=1; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(EXT, iRow, iCol, eVal);
	  mpq_sub(eVal, eVal, eVect->ListElt[iCol]);
	  MATRIX_Assign(nMat, iRow, iCol, eVal);
	}
    }
  while(1)
    {
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  a=rand();
	  b=rand();
	  eVal_i=a-b;
	  mpq_set_si(eVal, eVal_i, 1);
	  VECTOR_Assign(eVect, iCol, eVal);
	}
      CDD_LinearProgramming(nMat, eVect, eSol);
      VECTOR_ZeroAssignation(TheVert);
      eInc.clear();
      mpq_set_si(eVal, 1, 1);
      VECTOR_Assign(TheVert, 0, eVal);
      for (iCol=0; iCol<nbCol-1; iCol++)
	VECTOR_Assign(TheVert, iCol+1, eSol->DirectSolutionVal[iCol+1]);
      for (iRow=0; iRow<nbRow; iRow++)
	{
	  mpq_set_si(eSum, 0, 1);
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      MATRIX_Get(nMat, iRow, iCol, prov);
	      mpq_mul(prov, prov, TheVert->ListElt[iCol]);
	      mpq_add(eSum, eSum, prov);
	    }
	  test=mpq_sgn(eSum);
	  if (test == 0)
	    eInc.push_back(iRow);
	}
      siz=eInc.size();
      MATRIX_Allocate(RnkMat, siz, nbCol);
      iter=eInc.begin();
      idx=0;
      while(iter != eInc.end())
	{
	  iRow=*iter;
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      MATRIX_Get(nMat, iRow, iCol, eVal);
	      MATRIX_Assign(RnkMat, idx, iCol, eVal);
	    }
	  idx++;
	  iter++;
	}
      MATRIX_SelectColRow(RnkMat, eSelect, NSP);
      TheRank=eSelect->TheRank;
      MATRIX_Free(RnkMat);
      MATRIX_Free(NSP);
      SELECT_Free(eSelect);
      if (TheRank == nbCol-1)
	{
	  mpq_clear(eSum);
	  mpq_clear(eVal);
	  mpq_clear(prov);
	  MATRIX_Free(nMat);
	  VECTOR_Free(eVect);
	  VECTOR_Free(TheVert);
	  LPSOL_Free(eSol);
	  free(RnkMat);
	  free(nMat);
	  free(NSP);
	  free(eVect);
	  free(TheVert);
	  free(eSelect);
	  free(eSol);
	  return;
	}
    }
}


void FindOneInitialVertex(MyMatrix *EXT, vector<int> &eInc)
{
  MyMatrix *EXTret;
  if ((EXTret = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  Polytopization(EXT, EXTret);
  FindOneInitialVertexEXT(EXTret, eInc);
  MATRIX_Free(EXTret);
  free(EXTret);
}

