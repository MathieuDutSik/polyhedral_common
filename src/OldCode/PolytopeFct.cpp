#include "PolytopeFct.h"

void ReadPolytope(FILE *f, MyMatrix *TheMat)
{
  mpq_t t_q;
  int nbRow, nbCol;
  int iRow, iCol, ret;
  mpq_init(t_q);
  ret=fscanf(f, "%d %d", &nbRow, &nbCol);
  MATRIX_Allocate(TheMat, nbRow, nbCol);
  for (iRow=0; iRow<nbRow; iRow++)
    for (iCol=0; iCol<nbCol; iCol++)
      {
	mpq_inp_str(t_q, f, 10);
	MATRIX_Assign(TheMat, iRow, iCol, t_q);
      }
  mpq_clear(t_q);
}


void WritePolytope(FILE *f, MyMatrix *TheMat)
{
  mpq_t t_q;
  int nbRow, nbCol;
  int iRow, iCol;
  nbRow=TheMat->nbRow;
  nbCol=TheMat->nbCol;
  mpq_init(t_q);
  fprintf(f, "%d %d\n", nbRow, nbCol);
  for (iRow=0; iRow<nbRow; iRow++)
    {
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(TheMat, iRow, iCol, t_q);
	  fprintf(f, " ");
	  mpq_out_str(f, 10, t_q);
	}
      fprintf(f, "\n");
    }
  mpq_clear(t_q);
}

void MatrixReduction(MyMatrix *TheMat, vector<int> SelectedRows, MyMatrix *TheMatRed)
{
  MyMatrix TheProv, NSP;
  MySelection pSelect;
  mpq_t eVal;
  vector<int>::iterator iter;
  int nbCol1, nbRow1, iRow, jRow, iCol, jCol, TheRank;
  mpq_init(eVal);
  nbCol1=TheMat->nbCol;
  nbRow1=SelectedRows.size();
  MATRIX_Allocate(&TheProv, nbRow1, nbCol1);
  iter=SelectedRows.begin();
  iRow=0;
  while(iter != SelectedRows.end())
    {
      jRow=*iter;
      for (iCol=0; iCol<nbCol1; iCol++)
	{
	  MATRIX_Get(TheMat, jRow, iCol, eVal);
	  MATRIX_Assign(&TheProv, iRow, iCol, eVal);
	}
      iter++;
      iRow++;
    }
  MATRIX_SelectColRow(&TheProv, &pSelect, &NSP);
  TheRank=pSelect.TheRank;
  MATRIX_Allocate(TheMatRed, nbRow1, TheRank);
  for (iCol=0; iCol<TheRank; iCol++)
    {
      jCol=pSelect.ListColSelect[iCol];
      for (iRow=0; iRow<nbRow1; iRow++)
	{
	  MATRIX_Get(&TheProv, iRow, jCol, eVal);
	  MATRIX_Assign(TheMatRed, iRow, iCol, eVal);
	}
    }
  MATRIX_Free(&NSP);
  MATRIX_Free(&TheProv);
  SELECT_Free(&pSelect);
  mpq_clear(eVal);
}


void TestFacetness(MyMatrix *TheEXT, vector<int> OneInc)
{
  int nb, nbRow, nbCol;
  MyMatrix *TheProv, *NSP;
  MySelection *pSelect;
  int nbZero, nbPlus, nbMinus;
  int iRow, aRow, iCol, test;
  mpq_t eScal, eVal, prov;
  mpq_init(eScal);
  mpq_init(eVal);
  mpq_init(prov);
  nb=OneInc.size();
  nbRow=TheEXT->nbRow;
  nbCol=TheEXT->nbCol;
  if ((TheProv = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((NSP = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((pSelect = (MySelection*)malloc(sizeof(MySelection))) == 0)
    exit (EXIT_FAILURE);
  MATRIX_Allocate(TheProv, nb, nbCol);
  for (iRow=0; iRow<nb; iRow++)
    {
      aRow=OneInc[iRow];
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(TheEXT, aRow, iCol, eVal);
	  MATRIX_Assign(TheProv, iRow, iCol, eVal);
	}
    }
  MATRIX_SelectColRow(TheProv, pSelect, NSP);
  if (NSP->nbRow != 1)
    {
      fprintf(stderr, "Error in rank in Facetness\n");
      exit(1);
    }
  nbZero=0;
  nbPlus=0;
  nbMinus=0;
  for (iRow=0; iRow<nbRow; iRow++)
    {
      mpq_set_ui(eScal, 0, 1);
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(TheEXT, iRow, iCol, eVal);
	  MATRIX_Get(NSP, 0, iCol, prov);
	  mpq_mul(prov, prov, eVal);
	  mpq_add(eScal, eScal, prov);
	}
      test=mpq_sgn(eScal);
      if (test == 0)
	nbZero++;
      if (test == 1)
	nbPlus++;
      if (test == -1)
	nbMinus++;
    }
  if (nbZero != nb)
    {
      fprintf(stderr, "Error in computing incidence\n");
      exit(1);
    }
  if (nbMinus > 0 && nbPlus > 0)
    {
      fprintf(stderr, "Some plus and minus signs, illegal\n");
      exit(1);
    }
  MATRIX_Free(NSP);
  free(NSP);
  MATRIX_Free(TheProv);
  free(TheProv);
  SELECT_Free(pSelect);
  free(pSelect);
  mpq_clear(eScal);
  mpq_clear(eVal);
  mpq_clear(prov);
}



vector<int> ComputeFlipping(MyMatrix *TheEXT, vector<int> OneInc, vector<int> sInc)
{
  int nbCol, nb, nbRow;
  int iRow, iCol, aRow;
  vector<int> hSub;
  vector<vector<int> > TwoPlanes;
  MyMatrix TheProv, NSP, NSPtrans, PairFac, LV, LProd, ListScal;
  mpq_t eVal, prov1, prov2, prov3, prov4;
  mpq_t EXT1_1, EXT1_2, EXT2_1, EXT2_2;
  mpq_t TheDet, det12, det1N, det2N, prodDet, h;
  int test, nbForm, IsNonZero, i, idxFound;
  vector<int> ePlane;
  MySelection pSelect;
  nb=sInc.size();
  nbRow=TheEXT->nbRow;
  nbCol=TheEXT->nbCol;
  MATRIX_Allocate(&TheProv, nb, nbCol);
  MATRIX_Allocate(&LV, 2, 2);
  mpq_init(eVal);
  mpq_init(prov1);
  mpq_init(prov2);
  mpq_init(prov3);
  mpq_init(prov4);
  mpq_init(EXT1_1);
  mpq_init(EXT1_2);
  mpq_init(EXT2_1);
  mpq_init(EXT2_2);
  mpq_init(TheDet);
  mpq_init(det12);
  mpq_init(det1N);
  mpq_init(det2N);
  mpq_init(prodDet);
  mpq_init(h);
  for (iRow=0; iRow<nb; iRow++)
    {
      aRow=OneInc[sInc[iRow]];
      hSub.push_back(aRow);
      for (iCol=0; iCol<nbCol; iCol++)
	{
	  MATRIX_Get(TheEXT, aRow, iCol, eVal);
	  MATRIX_Assign(&TheProv, iRow, iCol, eVal);
	}
    }
  MATRIX_SelectColRow(&TheProv, &pSelect, &NSP);
  if (NSP.nbRow != 2)
    {
      fprintf(stderr, "Deep inconsistency\n");
      exit(1);
    }
  MATRIX_Transpose(&NSP, &NSPtrans);
  MATRIX_Product(TheEXT, &NSPtrans, &LProd);
  nbForm=0;
  for (iRow=0; iRow<nbRow; iRow++)
    {
      IsNonZero=0;
      for (i=0; i<2; i++)
	{
	  MATRIX_Get(&LProd, iRow, i, prov1);
	  test=mpq_sgn(prov1);
	  if (test != 0)
	    IsNonZero=1;
	}
      if (IsNonZero == 1)
	{
	  if (nbForm == 0)
	    {
	      MATRIX_Get(&LProd, iRow, 0, EXT1_1);
	      MATRIX_Get(&LProd, iRow, 1, EXT1_2);
	      nbForm++;
	    }
	  else
	    {
	      MATRIX_Get(&LProd, iRow, 0, prov1);
	      MATRIX_Get(&LProd, iRow, 1, prov2);
	      mpq_mul(prov3, prov2, EXT1_1);
	      mpq_mul(prov4, prov1, EXT1_2);
	      mpq_sub(TheDet, prov3, prov4);
	      test=mpq_sgn(TheDet);
	      if (nbForm == 1)
		{
		  if (test != 0)
		    {
		      MATRIX_Get(&LProd, iRow, 0, EXT2_1);
		      MATRIX_Get(&LProd, iRow, 1, EXT2_2);
		      mpq_set(det12, TheDet);
		      nbForm++;
		    }
		}
	      else
		{
		  mpq_set(det1N, TheDet);
		  mpq_mul(prov3, prov2, EXT2_1);
		  mpq_mul(prov4, prov1, EXT2_2);
		  mpq_sub(det2N, prov3, prov4);
		  mpq_mul(prodDet, det1N, det2N);
		  test=mpq_sgn(prodDet);
		  if (test == 1)
		    {
		      mpq_mul(prodDet, det12, det1N);
		      test=mpq_sgn(prodDet);
		      if (test == 1)
			{
			  mpq_set(EXT2_1, prov1);
			  mpq_set(EXT2_2, prov2);
			  mpq_set(det12, det1N);
			}
		      else
			{
			  mpq_set(EXT1_1, prov1);
			  mpq_set(EXT1_2, prov2);
			  mpq_neg(det12, det2N);
			}
		    }
		}
	    }
	}
    }
  test=mpq_sgn(det12);
  if (test == 1)
    {
      mpq_neg(h, EXT1_2);
      MATRIX_Assign(&LV, 0, 0, h);
      mpq_set(h, EXT1_1);
      MATRIX_Assign(&LV, 1, 0, h);
      mpq_set(h, EXT2_2);
      MATRIX_Assign(&LV, 0, 1, h);
      mpq_neg(h, EXT2_1);
      MATRIX_Assign(&LV, 1, 1, h);
    }
  else
    {
      mpq_set(h, EXT1_2);
      MATRIX_Assign(&LV, 0, 0, h);
      mpq_neg(h, EXT1_1);
      MATRIX_Assign(&LV, 1, 0, h);
      mpq_neg(h, EXT2_2);
      MATRIX_Assign(&LV, 0, 1, h);
      mpq_set(h, EXT2_1);
      MATRIX_Assign(&LV, 1, 1, h);
    }
  MATRIX_Product(&NSPtrans, &LV, &PairFac);
  MATRIX_Product(TheEXT, &PairFac, &ListScal);
  for (i=0; i<2; i++)
    {
      for (iRow=0; iRow<nbRow; iRow++)
	{
	  MATRIX_Get(&ListScal, iRow, i, eVal);
	  test=mpq_sgn(eVal);
	  if (test == -1)
	    {
	      fprintf(stderr, "This should never have happened. Please panic\n");
	      exit(1);
	    }
	  if (test == 0)
	    ePlane.push_back(iRow);
	}
      TwoPlanes.push_back(ePlane);
      ePlane.clear();
    }
  idxFound=-1;
  for (i=0; i<2; i++)
    if (TwoPlanes[i] == OneInc)
      idxFound=i;
  if (idxFound == -1)
    {
      fprintf(stderr, "We did not find the facet\n");
      exit(1);
    }
  MATRIX_Free(&TheProv);
  MATRIX_Free(&NSP);
  MATRIX_Free(&NSPtrans);
  MATRIX_Free(&PairFac);
  MATRIX_Free(&LV);
  MATRIX_Free(&LProd);
  MATRIX_Free(&ListScal);
  SELECT_Free(&pSelect);
  mpq_clear(eVal);
  mpq_clear(prov1);
  mpq_clear(prov2);
  mpq_clear(prov3);
  mpq_clear(prov4);
  mpq_clear(EXT1_1);
  mpq_clear(EXT1_2);
  mpq_clear(EXT2_1);
  mpq_clear(EXT2_2);
  mpq_clear(TheDet);
  mpq_clear(det12);
  mpq_clear(det1N);
  mpq_clear(det2N);
  mpq_clear(prodDet);
  mpq_clear(h);
  return TwoPlanes[1-idxFound];
}

void CompressedToExpanded(vector<long>& TheCompressed, vector<int>& TheExpanded, long& nbExt)
{
  long ThePow, Pow2, TheLong, res;
  long iExt;
  int ThePos;
  ThePow=1;
  ThePos=0;
  TheLong=TheCompressed[ThePos];
  TheExpanded.clear();
  for (iExt=1; iExt<=nbExt; iExt++)
    {
      Pow2=ThePow*2;
      res=TheLong%Pow2;
      if (res == ThePow)
	{
	  TheExpanded.push_back(1);
	  TheLong=TheLong-ThePow;
	}
      else
	TheExpanded.push_back(0);
      res=iExt % 30;
      if (res == 0 && iExt < nbExt)
	{
	  ThePos++;
	  TheLong=TheCompressed[ThePos];
	  ThePow=1;
	}
      else
	ThePow=ThePow*2;
    }
}

void ExpandedToCompressed(vector<long>& TheCompressed, vector<int>& TheExpanded, long& nbExt)
{
  long ThePow, TheLong, res;
  long iExt;
  int NeedToFlush=0;
  ThePow=1;
  TheLong=0;
  TheCompressed.clear();
  for (iExt=1; iExt<=nbExt; iExt++)
    {
      if (TheExpanded[iExt-1] == 1)
	TheLong=TheLong+ThePow;
      NeedToFlush=1;
      res=iExt % 30;
      if (res == 0)
	{
	  TheCompressed.push_back(TheLong);
	  TheLong=0;
	  NeedToFlush=0;
	  ThePow=1;
	}
      else
	ThePow=ThePow*2;
    }
  if (NeedToFlush == 1)
    TheCompressed.push_back(TheLong);
}

void PrintVectInt(FILE *f, vector<int> eList)
{
  int len, i;
  len=eList.size();
  for (i=0; i<len; i++)
    fprintf(f, " %d", eList[i]);
  fprintf(f, "\n");
}


void eEltImage(vector<long>& fSet, vector<long>& eSet, 
	       Permutation::ptr eElt, 
	       long& nbExt)
{
  vector<int> fSetInt, eSetInt;
  long iExt, jExt;
  //  fprintf(stderr, "nbExt=%ld\n", nbExt);
  //  fprintf(stderr, "Before CompressedToExpanded\n");
  CompressedToExpanded(eSet, eSetInt, nbExt);
  //  fprintf(stderr, "After CompressedToExpanded\n");
  //  fprintf(stderr, "eSetInt=\n");
  //  PrintVectInt(stderr, eSetInt);
  for (iExt=0; iExt<nbExt; iExt++)
    {
      jExt=eElt->at(iExt);
      //      fprintf(stderr, "iExt=%ld jExt=%ld eSetInt=%d\n", iExt, jExt, eSetInt[jExt]);
      fSetInt.push_back(eSetInt[jExt]);
    }
  ExpandedToCompressed(fSet, fSetInt, nbExt);
  //  fprintf(stderr, "fSetInt=\n");
  //  PrintVectInt(stderr, fSetInt);
}

void OrbitSplittingSet(set<vector<long> > ListTotal, 
		       TheGroupFormat &TheGRP, 
		       vector<vector<int> > &TheReturn)
{
  list<Permutation::ptr> ListGen;
  list<Permutation::ptr>::iterator iterGen;
  set<vector<long> > SingleOrbit, Additional, NewElts;
  set<vector<long> >::iterator iter, iterTEST1, iterTEST2, iterAddi, iterTotal, iterNew;
  list<vector<long> >::iterator iterGroup;
  vector<long> fSet, eSet, gSet, eVectLong;
  vector<int> eVectInt;
  vector<int> eVectAdd;
  long nbExt, iExt;
  fprintf(stderr, "TheGRP=\n");
  PrintGroup(stderr, TheGRP);
  nbExt=TheGRP.n;
  TheReturn.clear();
  ListGen=TheGRP.group->S;
  while(1)
    {
      iter=ListTotal.begin();
      if (iter == ListTotal.end())
	break;
      eSet=*iter;
      CompressedToExpanded(eSet, eVectInt, nbExt);
      eVectAdd.clear();
      for (iExt=0; iExt<nbExt; iExt++)
	if (eVectInt[iExt] == 1)
	  eVectAdd.push_back(iExt);
      TheReturn.push_back(eVectAdd);
      Additional.insert(eSet);
      ListTotal.erase(eSet);
      while(1)
	{
	  iterAddi=Additional.begin();
	  while(iterAddi != Additional.end())
	    {
	      gSet=*iterAddi;
	      iterGen=ListGen.begin();
	      while(iterGen != ListGen.end())
		{
		  eEltImage(fSet, gSet, *iterGen, nbExt);
		  iterTEST1=SingleOrbit.find(fSet);
		  iterTEST2=Additional.find(fSet);
		  if (iterTEST1 == SingleOrbit.end() && iterTEST2 == Additional.end())
		    {
		      iterNew=NewElts.find(fSet);
		      iterTotal=ListTotal.find(fSet);
		      if (iterNew == NewElts.end())
			{
			  if (iterTotal == ListTotal.end())
			    {
			      fprintf(stderr, "Orbit do not matched, PANIC!!!\n");
			      exit(1);
			    }
			  else
			    {
			      NewElts.insert(fSet);
			      ListTotal.erase(fSet);
			    }
			}
		    }
		  iterGen++;
		}
	      iterAddi++;
	    }
	  iterAddi=Additional.begin();
	  while(iterAddi != Additional.end())
	    {
	      SingleOrbit.insert(*iterAddi);
	      iterAddi++;
	    }
	  Additional.clear();
	  iterNew=NewElts.begin();
	  while(iterNew != NewElts.end())
	    {
	      Additional.insert(*iterNew);
	      iterNew++;
	    }
	  NewElts.clear();
	  if (Additional.size() == 0)
	    break;
	}
      SingleOrbit.clear();
    }
}




set<vector<long> > GetCompressedFromPolytope(MyMatrix *TheEXT, MyMatrix *TheFAC)
{
  list<Permutation::ptr> ListGen;
  list<Permutation::ptr>::iterator iterGen;
  set<vector<long> > SingleOrbit, Additional, NewElts, ListTotal;
  set<vector<long> >::iterator iter, iterTEST1, iterTEST2, iterAddi, iterTotal, iterNew;
  list<vector<long> >::iterator iterGroup;
  vector<long> fSet, eSet, gSet, eVectLong;
  vector<int> eVectInt;
  vector<int> eVectAdd;
  long nbFac, iFac;
  long nbExt, iExt;
  int nbCol, iCol, test;
  mpq_t prov1, prov2, eSum;
  nbFac=TheFAC->nbRow;
  nbCol=TheFAC->nbCol;
  nbExt=TheEXT->nbRow;
  mpq_init(prov1);
  mpq_init(prov2);
  mpq_init(eSum);
  for (iFac=0; iFac<nbFac; iFac++)
    {
      eVectInt.clear();
      for (iExt=0; iExt<nbExt; iExt++)
	{
	  mpq_set_si(eSum, 0, 1);
	  for (iCol=0; iCol<nbCol; iCol++)
	    {
	      MATRIX_Get(TheFAC, iFac, iCol, prov1);
	      MATRIX_Get(TheEXT, iExt, iCol, prov2);
	      mpq_mul(prov1, prov1, prov2);
	      mpq_add(eSum, eSum, prov1);
	    }
	  test=mpq_sgn(eSum);
	  if (test == 0)
	    eVectInt.push_back(1);
	  else
	    eVectInt.push_back(0);
	}
      ExpandedToCompressed(eVectLong, eVectInt, nbExt);
      ListTotal.insert(eVectLong);
    }
  mpq_clear(prov1);
  mpq_clear(prov2);
  mpq_clear(eSum);
  return ListTotal;
}





void OrbitSplitting(MyMatrix *TheEXT, MyMatrix *TheFAC, 
		    TheGroupFormat &TheGRP, 
		    vector<vector<int> > &TheReturn)
{
  set<vector<long> > ListTotal;
  ListTotal=GetCompressedFromPolytope(TheEXT, TheFAC);
  OrbitSplittingSet(ListTotal, TheGRP, TheReturn);
}


void PrintListOrbit(FILE *f, vector<vector<int> > &ListOrbit)
{
  int iOrbit, nbOrbit, siz, i, eVal;
  vector<int> eInc;
  nbOrbit=ListOrbit.size();
  fprintf(f, "nbOrbit=%d\n", nbOrbit);
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    {
      eInc=ListOrbit[iOrbit];
      siz=eInc.size();
      fprintf(f, "O%d: inc=%d list=", iOrbit+1, siz);
      for (i=0; i<siz; i++)
	{
	  eVal=eInc[i];
	  fprintf(f, " %d", eVal);
	}
      fprintf(f, "\n");
    }
}
