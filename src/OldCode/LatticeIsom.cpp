#include "LatticeIsom.h"

double DOUBL_GRAM_GetUpperBound(MatDOUBL TheMat)
{
  double TheLowEst, FudgeFact, MaxDet;
  int n, i;
  TheLowEst=0;
  n=TheMat.size();
  for (i=0; i<n; i++)
    if (TheMat[i][i] > TheLowEst)
      TheLowEst=TheMat[i][i];
  FudgeFact=1.1;
  MaxDet=TheLowEst*FudgeFact;
  return MaxDet;
}


void DOUBL_UpdateListValue(MatDOUBL eMat, double TheTol, vector<double> &ListVal)
{
  int n, i, j, p, WeFound, nb, iVal;
  double fVal, hVal;
  n=eMat.size();
  //  fprintf(stderr, "TheTol=%lg\n", TheTol);
  //  fprintf(stderr, "n=%d\n", n);
  for (i=0; i<n; i++)
    {
      p=eMat[i].size();
      //      fprintf(stderr, " i=%d p=%d\n", i, p);
      for (j=0; j<p; j++)
	{
	  WeFound=0;
	  nb=ListVal.size();
	  fVal=eMat[i][j];
	  //	  fprintf(stderr, "  j=%d nb=%d fVal=%lg\n", j, nb, fVal);
	  for (iVal=0; iVal<nb; iVal++)
	    {
	      hVal=ListVal[iVal];
	      //  fprintf(stderr, "   iVal=%d hVal=%lg\n", iVal, hVal);
	      if (fabs(hVal - fVal) < TheTol)
		WeFound=1;
	      //	      fprintf(stderr, "After the test\n");
	    }
	  //	  fprintf(stderr, "  WeFound=%d\n", WeFound);
	  if (WeFound == 0)
	    ListVal.push_back(fVal);
	}
    }
}

void PrintListValue(FILE *f, vector<double> &ListVal)
{
  int iVal, nbVal;
  nbVal=ListVal.size();
  fprintf(stderr, "nbVal=%d\n", nbVal);
  for (iVal=0; iVal<nbVal; iVal++)
    fprintf(f, "iVal=%d/%d eVal=%lg\n", iVal, nbVal, ListVal[iVal]);
}

void DOUBL_TranslateToMatrix(WeightMatrix& WMat, MatDOUBL eMat, double TheTol, vector<double> &ListVal)
{
  mpq_t eInv;
  mpq_class tcla;
  int nbRow, nbEnt, iRow, jRow, ThePos, idxMat;
  int iEnt;
  nbRow=eMat.size();
  nbEnt=ListVal.size();
  DefineWeightMatrix(WMat, nbRow);
  for (iRow=0; iRow<nbRow; iRow++)
    for (jRow=0; jRow<nbRow; jRow++)
      {
	ThePos=-1;
	for (iEnt=0; iEnt<nbEnt; iEnt++)
	  if (abs(eMat[iRow][jRow] - ListVal[iEnt]) < TheTol)
	    ThePos=iEnt;
	if (ThePos == -1)
	  {
	    fprintf(stderr, "ThePos wrongly assigned\n");
	    exit(1);
	  }
	idxMat=iRow + nbRow*jRow;
	WMat.TheMat[idxMat]=ThePos;
      }
  mpq_init(eInv);
  for (iEnt=0; iEnt<nbEnt; iEnt++)
    {
      mpq_set_si(eInv, iEnt, 1);
      tcla=mpq_class(eInv);
      WMat.ListWeight.push_back(tcla);
    }
  mpq_clear(eInv);
}

void DOUBL_GRAM_GetScalProdMat(MatDOUBL eMat, MatINT ListShort, MatDOUBL &ScalProdMat)
{
  int nbShort, n, iShort, jShort;
  vector<double> eLine;
  double eScal;
  int i, j;
  nbShort=ListShort.size();
  n=ListShort[0].size();
  for (iShort=0; iShort<nbShort; iShort++)
    {
      for (jShort=0; jShort<nbShort; jShort++)
	{
	  eScal=0;
	  for (i=0; i<n; i++)
	    for (j=0; j<n; j++)
	      eScal=eScal + ListShort[iShort][i]*eMat[i][j]*ListShort[jShort][j];
	  eLine.push_back(eScal);
	}
      ScalProdMat.push_back(eLine);
      eLine.clear();
    }
}

void EnumerationShortVectorAntipodal(MatDOUBL TheMat, double MaxNorm, 
				     MatINT & ListShort)
{
  MatINT PreListShort;
  vector<int> eVect, eShort;
  int iShort, n, i, nbPreShort;
  vector<double> ListNorm;
  EnumerationShortVector(TheMat, MaxNorm, PreListShort, ListNorm);
  nbPreShort=PreListShort.size();
  for (iShort=0; iShort<nbPreShort; iShort++)
    {
      eShort=PreListShort[iShort];
      ListShort.push_back(eShort);
      eVect.clear();
      n=eShort.size();
      for (i=0; i<n; i++)
	eVect.push_back(-eShort[i]);
      ListShort.push_back(eVect);
    }
}


void DOUBL_GetGramMatrixAutomorphismGroup(MatDOUBL eMat, double TheTol, TheGroupFormat &GRPperm, vector<MatINT> &ListMatrGens)
{
  MatINT ListShort;
  double MaxDet;
  MatDOUBL ScalProdMat;
  vector<double> ListVal;
  WeightMatrix WMat;
  vector<double> eVect;
  MyMatrix ListShortGMP;
  mpq_t eVal_q;
  MySelection eSelect;
  MyMatrix NSP;
  MyMatrix M1, M2, M3, M1inv;
  MatINT eMatrGen;
  int iShort, jShort, nbShort, i, j, n, eInt;
  list<Permutation::ptr> ListGen;
  list<Permutation::ptr>::iterator iter;
  n=eMat.size();
  MaxDet=DOUBL_GRAM_GetUpperBound(eMat);
  EnumerationShortVectorAntipodal(eMat, MaxDet, ListShort);
  nbShort=ListShort.size();
  MATRIX_Allocate(&ListShortGMP, nbShort, n);
  mpq_init(eVal_q);
  for (iShort=0; iShort<nbShort; iShort++)
    {
      for (i=0; i<n; i++)
	{
	  eInt=ListShort[iShort][i];
	  mpq_set_si(eVal_q, eInt, 1);
	  MATRIX_Assign(&ListShortGMP, iShort, i, eVal_q);
	}
      eVect.clear();
    }
  DOUBL_GRAM_GetScalProdMat(eMat, ListShort, ScalProdMat);
  DOUBL_UpdateListValue(ScalProdMat, TheTol, ListVal);
  DOUBL_TranslateToMatrix(WMat, ScalProdMat, TheTol, ListVal);
  GetStabilizerWeightMatrix(WMat, GRPperm);
  ListGen=GRPperm.group->S;
  iter=ListGen.begin();
  MATRIX_SelectColRow(&ListShortGMP, &eSelect, &NSP);
  MATRIX_Allocate(&M1, n, n);
  MATRIX_Allocate(&M1inv, n, n);
  MATRIX_Allocate(&M2, n, n);
  MATRIX_Allocate(&M3, n, n);
  for (i=0; i<n; i++)
    {
      iShort=eSelect.ListRowSelect[i];
      for (j=0; j<n; j++)
	{
	  MATRIX_Get(&ListShortGMP, iShort, j, eVal_q);
	  MATRIX_Assign(&M1, i, j, eVal_q);
	}
    }
  MATRIX_Inverse(&M1, &M1inv);
  while(iter != ListGen.end())
    {
      for (i=0; i<n; i++)
	{
	  iShort=eSelect.ListRowSelect[i];
	  jShort=(*iter)->at(iShort);
	  for (j=0; j<n; j++)
	    {
	      MATRIX_Get(&ListShortGMP, jShort, j, eVal_q);
	      MATRIX_Assign(&M2, i, j, eVal_q);
	    }
	}
      MATRIX_Product(&M1inv, &M2, &M3);
      MATRIX_ConvertToInt(&M3, eMatrGen);
      ListMatrGens.push_back(eMatrGen);
      iter++;
    }
  mpq_clear(eVal_q);
  MATRIX_Free(&M1);
  MATRIX_Free(&M1inv);
  MATRIX_Free(&M2);
  MATRIX_Free(&M3);
}

void DOUBL_TestGramMatrixEquivalence(MatDOUBL eMat1, MatDOUBL eMat2, double TheTol, int &test)
{
  double MaxDet1, MaxDet2, MaxDet;
  MatINT ListShort1, ListShort2;
  MatDOUBL ScalProdMat1, ScalProdMat2;
  vector<double> ListVal;
  WeightMatrix WMat1, WMat2;
  int TheReply;
  vector<int> eListEquiv;
  int nbShort1, nbShort2;
  MaxDet1=DOUBL_GRAM_GetUpperBound(eMat1);
  MaxDet2=DOUBL_GRAM_GetUpperBound(eMat2);
  if (MaxDet1 > MaxDet2)
    MaxDet=MaxDet2;
  else
    MaxDet=MaxDet1;
  EnumerationShortVectorAntipodal(eMat1, MaxDet, ListShort1);
  EnumerationShortVectorAntipodal(eMat2, MaxDet, ListShort2);
  DOUBL_GRAM_GetScalProdMat(eMat1, ListShort1, ScalProdMat1);
  DOUBL_GRAM_GetScalProdMat(eMat2, ListShort2, ScalProdMat2);
  //  fprintf(stderr, "ScalProdMat1=\n");
  //  PrintMat_Doubl(stderr, ScalProdMat1);
  //  fprintf(stderr, "ScalProdMat2=\n");
  //  PrintMat_Doubl(stderr, ScalProdMat2);

  nbShort1=ListShort1.size();
  nbShort2=ListShort2.size();
  //  fprintf(stderr, "nbShort1=%d nbShort2=%d\n", nbShort1, nbShort2);
  //  fprintf(stderr, "eMat1\n");
  //  PrintMat_Doubl(stderr, eMat1);
  //  fprintf(stderr, "eMat2\n");
  //  PrintMat_Doubl(stderr, eMat2);
  DOUBL_UpdateListValue(ScalProdMat1, TheTol, ListVal);
  DOUBL_UpdateListValue(ScalProdMat2, TheTol, ListVal);
  //  PrintListValue(stderr, ListVal);
  DOUBL_TranslateToMatrix(WMat1, ScalProdMat1, TheTol, ListVal);
  DOUBL_TranslateToMatrix(WMat2, ScalProdMat2, TheTol, ListVal);
  //  fprintf(stderr, "WMat1=\n");
  //  PrintWeightedMatrix(stderr, WMat1);
  //  fprintf(stderr, "WMat2=\n");
  //  PrintWeightedMatrix(stderr, WMat2);
  TestEquivalenceWeightMatrix(WMat1, WMat2, &TheReply, eListEquiv);
  //  fprintf(stderr, "TheReply=%d\n", TheReply);
  test=TheReply;
}

