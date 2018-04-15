#include "PolytopeEquiStab.h"

Graph* ReadGraphFromFile(FILE *f, unsigned int &nof_vertices)
{
  unsigned int nof_edges;
  unsigned int color;
  int ret, iEdge, i;
  int a, b;
  Graph *g =0;
  ret=fscanf(f, "%u %u\n", &nof_vertices, &nof_edges);
  g = new Graph(nof_vertices);
  for (i=0; i<(int)nof_vertices; i++)
    {
      ret=fscanf(f, "%u\n", &color);
      g->change_color(i, color);
    }
  for (iEdge=0; iEdge<(int)nof_edges; iEdge++)
    {
      ret=fscanf(f, "%u %u\n", &a, &b);
      g->add_edge(a-1, b-1);
    }
  return g;
}

void PrintAllInfoColoredGraph(FILE *f, TheRecGraph *TheRec)
{
  int nbVert, nbEdge, i, j, idxMat, eVal;
  nbVert=TheRec->nbVertices;
  nbEdge=0;
  for (i=0; i<nbVert-1; i++)
    {
      for (j=i+1; j<nbVert; j++)
	{
	  idxMat=i + nbVert*j;
	  eVal=TheRec->AdjMat[idxMat];
	  if (eVal == 1)
	    nbEdge++;
	}
    }
  fprintf(f, " %d %d", nbVert, nbEdge);
  for (i=0; i<nbVert; i++)
    fprintf(f, " %d\n", TheRec->ListColor[i]);
  fprintf(f, "%d", nbEdge);
  for (i=0; i<nbVert-1; i++)
    {
      for (j=i+1; j<nbVert; j++)
	{
	  idxMat=i + nbVert*j;
	  eVal=TheRec->AdjMat[idxMat];
	  if (eVal == 1)
	    fprintf(f, " %d %d", i+1, j+1);
	}
    }
}



void PrintAllInfoBlissFormat(FILE *f, TheRecGraph *TheRec)
{
  int nbVert, nbEdge, i, j, idxMat, eVal;
  nbVert=TheRec->nbVertices;
  nbEdge=0;
  for (i=0; i<nbVert-1; i++)
    {
      for (j=i+1; j<nbVert; j++)
	{
	  idxMat=i + nbVert*j;
	  eVal=TheRec->AdjMat[idxMat];
	  if (eVal == 1)
	    nbEdge++;
	}
    }
  fprintf(f, "p edge %d %d\n", nbVert, nbEdge);
  for (i=0; i<nbVert; i++)
    fprintf(f, "n %d %d\n", i+1, TheRec->ListColor[i]);
  for (i=0; i<nbVert-1; i++)
    {
      for (j=i+1; j<nbVert; j++)
	{
	  idxMat=i + nbVert*j;
	  eVal=TheRec->AdjMat[idxMat];
	  if (eVal == 1)
	    fprintf(f, "e %d %d\n", i+1, j+1);
	}
    }
}









void DefineWeightMatrix(WeightMatrix &WMat, int nbRow)
{
  int nb;
  nb=nbRow*nbRow;
  WMat.nbRow=nbRow;
  if ((WMat.TheMat = (int*)malloc(nb*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
}




void UpdateWeightMatrixAndGetPosition(WeightMatrix &WMat, int iCol, int iRow, mpq_t eVal)
{
  int nbEnt, WeFound, i, test, ThePos, idxMat;
  mpq_class tcla;
  WeFound=0;
  ThePos=-437;
  nbEnt=WMat.ListWeight.size();
  for (i=0; i<nbEnt; i++)
    if (WeFound == 0)
      {
	test=mpq_equal(eVal, WMat.ListWeight[i].get_mpq_t());
	if (test != 0)
	  {
	    WeFound=1;
	    ThePos=i;
	  }
      }
  if (WeFound == 0)
    {
      tcla=mpq_class(eVal);
      WMat.ListWeight.push_back(tcla);
      ThePos=nbEnt;
    }
  idxMat=iRow+WMat.nbRow*iCol;
  WMat.TheMat[idxMat]=ThePos;
}


void GetWeightMatrix(WeightMatrix &WMat, MyMatrix *TheEXT)
{
  MyMatrix QMat, QMatInv;
  mpq_t eSum, eVal, eVal1, eVal2;
  int iRow, jRow, iCol, jCol;
  int nbRow, nbCol;
  mpq_init(eSum);
  mpq_init(eVal);
  mpq_init(eVal1);
  mpq_init(eVal2);
  nbRow=TheEXT->nbRow;
  nbCol=TheEXT->nbCol;
  MATRIX_Allocate(&QMat, nbCol, nbCol);
  MATRIX_Allocate(&QMatInv, nbCol, nbCol);
  DefineWeightMatrix(WMat, nbRow);
  for (iCol=0; iCol<nbCol; iCol++)
    {
      for (jCol=0; jCol<nbCol; jCol++)
	{
	  mpq_set_ui(eSum, 0, 1);
	  for (iRow=0; iRow<nbRow; iRow++)
	    {
	      MATRIX_Get(TheEXT, iRow, iCol, eVal1);
	      MATRIX_Get(TheEXT, iRow, jCol, eVal2);
	      mpq_mul(eVal, eVal1, eVal2);
	      mpq_add(eSum, eSum, eVal);
	    }
	  MATRIX_Assign(&QMat, iCol, jCol, eSum);
	}
    }
  MATRIX_Inverse_destroy(&QMat, &QMatInv);
  for (iRow=0; iRow<nbRow; iRow++)
    for (jRow=0; jRow<nbRow; jRow++)
      {
	mpq_set_ui(eSum, 0, 1);
	for (iCol=0; iCol<nbCol; iCol++)
	  for (jCol=0; jCol<nbCol; jCol++)
	    {
	      MATRIX_Get(TheEXT, iRow, iCol, eVal1);
	      MATRIX_Get(TheEXT, jRow, jCol, eVal2);
	      mpq_mul(eVal, eVal1, eVal2);
	      MATRIX_Get(&QMatInv, iCol, jCol, eVal1);
	      mpq_mul(eVal, eVal, eVal1);
	      mpq_add(eSum, eSum, eVal);
	    }
	UpdateWeightMatrixAndGetPosition(WMat, iRow, jRow, eSum);
      }
  mpq_clear(eSum);
  mpq_clear(eVal);
  mpq_clear(eVal1);
  mpq_clear(eVal2);
  MATRIX_Free(&QMat);
  MATRIX_Free(&QMatInv);
}



typedef vector<vector<int> > VectVectInt;

static void report_aut_vectvectint(void* param, const unsigned int n, const unsigned int* aut)
{
  int i;
  int siz;
  vector<int> eVect;
  siz=n;
  for (i=0; i<siz; i++)
    eVect.push_back(aut[i]);
  ((VectVectInt*)param)->push_back(eVect);
}

static void report_aut_void(void* param, const unsigned int n, const unsigned int* aut)
{
  int siz;
  siz=n;
}


/*
static void report_aut_print(void* param, const unsigned int n, const unsigned int* aut)
{
  assert(param);
  fprintf((FILE*)param, "Generator: ");
  bliss::print_permutation((FILE*)param, n, aut, 1);
  fprintf((FILE*)param, "\n");
}
*/


void GetGraphAutomorphismGroup(Graph *g, unsigned int nof_vertices, TheGroupFormat &GRP)
{
  int *gList;
  int nbGen, iGen, i;
  list<Permutation::ptr> generatorList;
  bliss::Stats stats;
  VectVectInt ListGen;
  VectVectInt *h;
  h=&ListGen;
  //  fprintf(stderr, "Before running find_automorphism\n");
  g->find_automorphisms(stats, &report_aut_vectvectint, (void *)h);
  delete g;
  /* g->find_automorphisms(stats, &report_aut_print, stderr); */
  //  fprintf(stderr, "After running find_automorphism\n");
  nbGen=ListGen.size();
  if ((gList = (int*)malloc(nof_vertices*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  GRP.n=nof_vertices;
  for (iGen=0; iGen<nbGen; iGen++)
    {
      for (i=0; i<(int)nof_vertices; i++)
	{
	  //	  fprintf(stderr, " %d", ListGen[iGen][i]);
	  gList[i]=ListGen[iGen][i];
	}
      fprintf(stderr, "\n");
      vector<dom_int> v(gList, gList+nof_vertices);
      generatorList.push_back(Permutation::ptr(new Permutation(v)));
    }
  free(gList);
  GRP.n=nof_vertices;
  GRP.group=construct(nof_vertices, generatorList.begin(), generatorList.end());
}

int GetNeededPower(int nb)
{
  int h, eExpo;
  h=0;
  eExpo=1;
  while(1)
    {
      if (nb < eExpo)
	return h;
      h++;
      eExpo=eExpo*2;
    }
}

void GetBinaryExpression(int eVal, int *eVect, int h)
{
  int eWork, eExpo, eExpoB, i, res;
  eWork=eVal;
  eExpo=1;
  for (i=0; i<h; i++)
    {
      eExpoB=eExpo*2;
      res=eWork % eExpoB;
      eVect[i]=res/eExpo;
      eExpo=eExpoB;
      eWork=eWork - res;
    }
}

vector<int> GetLocalInvariantWeightMatrix(WeightMatrix &WMat, vector<int> eSet)
{
  int nbVert, iVert, jVert, aVert, bVert, iPos;
  int iWeight, nbWeight;
  vector<int> eInv;
  nbVert=eSet.size();
  nbWeight=WMat.ListWeight.size();
  eInv.clear();
  for (iWeight=0; iWeight<nbWeight; iWeight++)
    eInv.push_back(0);
  for (iVert=0; iVert<nbVert; iVert++)
    {
      aVert=eSet[iVert];
      for (jVert=0; jVert<nbVert; jVert++)
	{
	  bVert=eSet[jVert];
	  iPos=aVert+WMat.nbRow*bVert;
	  iWeight=WMat.TheMat[iPos];
	  eInv[iWeight]++;
	}
    }
  return eInv;
}




void GetInvariantWeightMatrix(WeightMatrix &WMat, vector<int> eSet, mpq_t eMainInv)
{
  int *ListM;
  int *ListAtt;
  mpq_t ePow, prov1, prov2, eInv;
  int nbInv;
  int nbVert, iVert, jVert, aVert, bVert, iInv;
  int iWeight, nbWeight;
  mpq_init(eInv);
  mpq_init(ePow);
  mpq_init(prov1);
  mpq_init(prov2);
  nbInv=3;
  nbVert=eSet.size();
  if ((ListM = (int*)malloc(nbInv*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  ListM[0]=2;
  ListM[1]=5;
  ListM[2]=23;
  nbWeight=WMat.ListWeight.size();
  if ((ListAtt = (int*)malloc(nbWeight*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (iWeight=0; iWeight<nbWeight; iWeight++)
    ListAtt[iWeight]=0;
  for (iVert=0; iVert<nbVert; iVert++)
    {
      aVert=eSet[iVert];
      for (jVert=0; jVert<nbVert; jVert++)
	{
	  bVert=eSet[jVert];
	  iWeight=aVert+WMat.nbRow*bVert;
	  ListAtt[iWeight]++;
	}
    }
  mpq_set_si(eMainInv, 1, 1);
  for (iInv=0; iInv<nbInv; iInv++)
    {
      mpq_set_si(eInv, 0,1);
      mpq_set_si(ePow, ListM[iInv],1);
      for (iWeight=0; iWeight<nbWeight; iWeight++)
	{
	  mpq_set_si(prov1, ListAtt[iWeight], 1);
	  mpq_mul(prov2, prov1, WMat.ListWeight[iWeight].get_mpq_t());
	  mpq_mul(eInv, ePow, eInv);
	  mpq_add(eInv, eInv, prov2);
	}
      mpq_mul(eMainInv, eMainInv, eInv);
    }
  free(ListM);
  mpq_clear(eInv);
  mpq_clear(ePow);
  mpq_clear(prov1);
  mpq_clear(prov2);
}


void ReadWeightedMatrix(FILE *f, WeightMatrix &WMat)
{
  int nbRow, iRow, jRow, idxMat, eVal;
  int iEnt, nbEnt, ret;
  mpq_t eInv;
  mpq_class tcla;
  ret=fscanf(f, "%d", &nbRow);
  DefineWeightMatrix(WMat, nbRow);
  nbEnt=0;
  for (iRow=0; iRow<nbRow; iRow++)
    for (jRow=0; jRow<nbRow; jRow++)
      {
	ret=fscanf(f, "%d", &eVal);
	idxMat=iRow + nbRow*jRow;
	WMat.TheMat[idxMat]=eVal;
	if (eVal> nbEnt)
	  nbEnt=eVal;
      }
  nbEnt++;
  mpq_init(eInv);
  for (iEnt=0; iEnt<nbEnt; iEnt++)
    {
      mpq_set_si(eInv, iEnt, 1);
      tcla=mpq_class(eInv);
      WMat.ListWeight.push_back(tcla);
    }
  mpq_clear(eInv);  
}

void PrintWeightedMatrix(FILE *f, WeightMatrix &WMat)
{
  int i, siz, iRow, iCol, idx, nbRow, eVal;
  siz=WMat.ListWeight.size();
  fprintf(f, "siz=%d\n", siz);
  for (i=0; i<siz; i++)
    {
      fprintf(f, "  i=%d eWeight=", i);
      mpq_out_str(f, 10, WMat.ListWeight[i].get_mpq_t());
      fprintf(f, "\n");
    }
  nbRow=WMat.nbRow;
  fprintf(f, "nbRow=%d\n", WMat.nbRow);
  for (iRow=0; iRow<nbRow; iRow++)
    {
      for (iCol=0; iCol<nbRow; iCol++)
	{
	  idx=iRow+nbRow*iCol;
	  eVal=WMat.TheMat[idx];
	  fprintf(f, " %d", eVal);
	}
      fprintf(f, "\n");
    }
}



Graph* GetGraphFromWeightedMatrix(WeightMatrix &WMat, TheRecGraph *TheRec)
{
  Graph *g = 0;
  unsigned int nof_vertices;
  int nbMult, nbVert;
  int *eVect;
  int aVert, bVert, iVert, jVert;
  int iH, jH, nbRow, eVal;
  int idx, nbEdge, hS;
  nbMult=WMat.ListWeight.size()+2;
  hS=GetNeededPower(nbMult);
  TheRec->hSize=hS;
  nbRow=WMat.nbRow;
  nbVert=nbRow + 2;
  nof_vertices=hS*nbVert;
  nbEdge=nof_vertices*nof_vertices;
  TheRec->nbVertices=nof_vertices;
  if ((TheRec->AdjMat = (int*)malloc(nbEdge*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (idx=0; idx<nbEdge; idx++)
    TheRec->AdjMat[idx]=0;
  if ((TheRec->ListColor = (int*)malloc(nof_vertices*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  g = new Graph(nof_vertices);
  if ((eVect = (int*)malloc(hS*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (iVert=0; iVert<nbVert; iVert++)
    {
      for (iH=0; iH<hS; iH++)
	{
	  aVert=iVert + nbVert*iH;
	  g->change_color(aVert, iH);
	  TheRec->ListColor[aVert]=iH;
	}
    }
  for (iVert=0; iVert<nbVert; iVert++)
    {
      for (iH=0; iH<hS-1; iH++)
	{
	  for (jH=iH+1; jH<hS; jH++)
	    {
	      aVert=iVert + nbVert*iH;
	      bVert=iVert + nbVert*jH;
	      idx=aVert+nof_vertices*bVert;
	      TheRec->AdjMat[idx]=1;
	      idx=bVert+nof_vertices*aVert;
	      TheRec->AdjMat[idx]=1;
	      g->add_edge(aVert, bVert);
	    }
	}
    }
  for (iVert=0; iVert<nbVert-1; iVert++)
    {
      for (jVert=iVert+1; jVert<nbVert; jVert++)
	{
	  if (jVert == nbRow+1)
	    {
	      if (iVert == nbRow)
		eVal=nbMult;
	      else
		eVal=nbMult+1;
	    }
	  else
	    {
	      if (jVert == nbRow)
		{
		  idx=iVert+nbRow*iVert;
		  eVal=WMat.TheMat[idx];
		}
	      else
		{
		  idx=iVert+nbRow*jVert;
		  eVal=WMat.TheMat[idx];
		}
	    }
	  //	  fprintf(stderr, "eVal=%d nbMult=%d iVert=%d jVert=%d nbRow=%d\n", eVal, nbMult, iVert, jVert, nbRow);
	  GetBinaryExpression(eVal, eVect, hS);
	  for (iH=0; iH<hS; iH++)
	    if (eVect[iH] == 1)
	      {
		aVert=iVert + nbVert*iH;
		bVert=jVert + nbVert*iH;
		idx=aVert+nof_vertices*bVert;
		TheRec->AdjMat[idx]=1;
		idx=bVert+nof_vertices*aVert;
		TheRec->AdjMat[idx]=1;
		g->add_edge(aVert, bVert);
	      }
	}
    }
  free(eVect);
  return g;
}

void RenormalizeWeightMatrix(WeightMatrix WMatRef, WeightMatrix &WMat2, int &test)
{
  int *gList;
  int *gListRev;
  int jFound, nbEnt, nbEnt2, i;
  int eValJ, eValI;
  int nbRow, nbRow2, j, nb, idx;
  nbRow=WMatRef.nbRow;
  nbRow2=WMat2.nbRow;
  if (nbRow != nbRow2)
    {
      test=0;
      return;
    }
  if ((gList = (int*)malloc(nbRow*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((gListRev = (int*)malloc(nbRow*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  nbEnt=WMatRef.ListWeight.size();
  nbEnt2=WMat2.ListWeight.size();
  if (nbEnt != nbEnt2)
    {
      test=0;
      return;
    }
  for (i=0; i<nbEnt; i++)
    {
      jFound=-1;
      for (j=0; j<nbEnt; j++)
	{
	  if (WMatRef.ListWeight[i] == WMat2.ListWeight[j])
	    jFound=j;
	}
      if (jFound == -1)
	{
	  test=0;
	  return;
	}
      gList[i]=jFound;
      gListRev[jFound]=i;
    }
  /*
  fprintf(stderr, "   gList=");
  for (i=0; i<nbEnt; i++)
    fprintf(stderr, " %d", gList[i]);
  fprintf(stderr, "\n");
  fprintf(stderr, "gListRev=");
  for (i=0; i<nbEnt; i++)
    fprintf(stderr, " %d", gListRev[i]);
    fprintf(stderr, "\n");*/
  nb=nbRow*nbRow;
  for (idx=0; idx<nb; idx++)
    {
      eValJ=WMat2.TheMat[idx];
      eValI=gListRev[eValJ];
      WMat2.TheMat[idx]=eValI;
    }
  WMat2.ListWeight=WMatRef.ListWeight;
  test=1;
  free(gList);
  free(gListRev);
}




void GetStabilizerWeightMatrix(WeightMatrix WMat, TheGroupFormat &GRP)
{
  Graph *g = 0;
  list<Permutation::ptr> generatorList;
  bliss::Stats stats;
  VectVectInt ListGen;
  VectVectInt *h;
  int *gList;
  int iVert, jVert;
  int nbRow, nbGen, iGen;
  TheRecGraph TheRec;
  int iRow, jRow, iRowI, jRowI, idxMat1, idxMat2;
  nbRow=WMat.nbRow;
  g=GetGraphFromWeightedMatrix(WMat, &TheRec);
  h=&ListGen;
  g->find_automorphisms(stats, &report_aut_vectvectint, (void *)h);
  if ((gList = (int*)malloc(nbRow*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  nbGen=ListGen.size();
  for (iGen=0; iGen<nbGen; iGen++)
    {
      for (iVert=0; iVert<nbRow; iVert++)
	{
	  jVert=ListGen[iGen][iVert];
	  if (jVert >= nbRow)
	    {
	      fprintf(stderr, "jVert is too large\n");
	      exit(1);
	    }
	  gList[iVert]=jVert;
	}
      for (iRow=0; iRow<nbRow; iRow++)
	for (jRow=0; jRow<nbRow; jRow++)
	  {
	    iRowI=gList[iRow];
	    jRowI=gList[jRow];
	    idxMat1=iRow +nbRow*jRow;
	    idxMat2=iRowI+nbRow*jRowI;
	    if (WMat.TheMat[idxMat1] != WMat.TheMat[idxMat2])
	      {
		fprintf(stderr, "Clear error in automorphism computation\n");
		fprintf(stderr, "iRow=%d jRow=%d\n", iRow, jRow);
		exit(1);
	      }
	  }
      vector<dom_int> v(gList, gList+nbRow);
      generatorList.push_back(Permutation::ptr(new Permutation(v)));
    }
  free(gList);
  GRP.n=nbRow;
  GRP.group=construct(nbRow, generatorList.begin(), generatorList.end());
  delete g;
  free(TheRec.AdjMat);
  free(TheRec.ListColor);
}


void PrintAllInfoBlissFormat(FILE *f, TheRecGraph *TheRec, const unsigned int *cl)
{
  int nbVert, nbEdge, i, j, idxMat, eVal;
  nbVert=TheRec->nbVertices;
  nbEdge=0;
  for (i=0; i<nbVert-1; i++)
    {
      for (j=i+1; j<nbVert; j++)
	{
	  idxMat=i + nbVert*j;
	  eVal=TheRec->AdjMat[idxMat];
	  if (eVal == 1)
	    nbEdge++;
	}
    }
  fprintf(f, "p edge %d %d\n", nbVert, nbEdge);
  for (i=0; i<nbVert; i++)
    fprintf(f, "n %d %d\n", i+1, TheRec->ListColor[i]);
  for (i=0; i<nbVert-1; i++)
    {
      for (j=i+1; j<nbVert; j++)
	{
	  idxMat=i + nbVert*j;
	  eVal=TheRec->AdjMat[idxMat];
	  if (eVal == 1)
	    fprintf(f, "e %d %d\n", i+1, j+1);
	}
    }
}



void PrintRecAdjMatCanoc(FILE *f, TheRecGraph *TheRec, const unsigned int *cl)
{
  int nbVert, i, j, iM, jM, idxMat, eVal;
  nbVert=TheRec->nbVertices;
  fprintf(f, "nbVert=%d\n", nbVert);
  set<int> eEdge;
  set<int>::iterator iterEdge;
  set<set<int> > ListEdges;
  set<set<int> >::iterator iterList;
  vector<int> LDegC;
  for (i=0; i<nbVert; i++)
    {
      LDegC.push_back(0);
    }
  for (i=0; i<nbVert; i++)
    for (j=0; j<nbVert; j++)
      {
	iM=cl[i];
	jM=cl[j];
	idxMat=iM + nbVert*jM;
	eVal=TheRec->AdjMat[idxMat];
	if (eVal == 1)
	  {
	    LDegC[i]=LDegC[i]+1;
	    eEdge.insert(i);
	    eEdge.insert(j);
	    ListEdges.insert(eEdge);
	    eEdge.clear();
	  }
      }
  for (i=0; i<nbVert; i++)
    fprintf(f, "i=%d eDegC=%d\n", i, LDegC[i]);
  iterList=ListEdges.begin();
  while(iterList != ListEdges.end())
    {
      eEdge=*iterList;
      iterEdge=eEdge.begin();
      while(iterEdge != eEdge.end())
	{
	  fprintf(f, " %d", *iterEdge);
	  iterEdge++;
	}
      fprintf(f, "\n");
      iterList++;
    }
}

void PrintRecAdjMatCanocR(FILE *f, TheRecGraph *TheRec, const unsigned int *cl)
{
  int nbVert, i;
  unsigned int *clR;
  nbVert=TheRec->nbVertices;
  if ((clR = (unsigned int*)malloc(nbVert*sizeof(unsigned int))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<nbVert; i++)
    clR[cl[i]]=i;
  PrintRecAdjMatCanoc(f, TheRec, clR);
  free(clR);
}



void PrintRecAdjMat(FILE *f, TheRecGraph *TheRec)
{
  int nbVert, i, j, idxMat, eVal;
  nbVert=TheRec->nbVertices;
  fprintf(f, "nbVertices=%d\n", nbVert);
  for (i=0; i<nbVert; i++)
    {
      for (j=0; j<nbVert; j++)
	{
	  if (i != j)
	    {
	      idxMat=i + nbVert*j;
	      eVal=TheRec->AdjMat[idxMat];
	      fprintf(f, " %d", eVal);
	    }
	  else
	    fprintf(f, " %d", TheRec->ListColor[i]);
	}
      fprintf(f, "\n");
    }
}

void FreeRecAdjMat(TheRecGraph *TheRec)
{
  free(TheRec->ListColor);
  free(TheRec->AdjMat);
}


void TestEquivalenceWeightMatrix(WeightMatrix WMat1, WeightMatrix &WMat2, int *TheReply, vector<int> & TheEquiv)
{
  const unsigned int* cl1;
  const unsigned int* cl2;
  int nof_vertices;
  TheRecGraph TheRec1, TheRec2;
  int idx1, idx2, test, i;
  int iVert1, iVert2, jVert1, jVert2, iVert, jVert;
  Graph *g1, *g2;
  bliss::Stats stats;
  unsigned int *clR1;
  unsigned int *clR2;
  vector<int> TheEquivExp;
  int nbRow;
  RenormalizeWeightMatrix(WMat1, WMat2, test);
  if (test == 0)
    {
      fprintf(stderr, "TheReply=0 case 1\n");
      *TheReply=0;
      return;
    }
  fprintf(stderr, "TestEquivalenceWeightMatrix\n");
  fprintf(stderr, "WMat1=\n");
  PrintWeightedMatrix(stderr, WMat1);
  fprintf(stderr, "WMat2=\n");
  PrintWeightedMatrix(stderr, WMat2);
  g1=GetGraphFromWeightedMatrix(WMat1, &TheRec1);
  g2=GetGraphFromWeightedMatrix(WMat2, &TheRec2);
  if ((TheRec1.nbVertices != TheRec2.nbVertices) || (TheRec1.hSize != TheRec2.hSize))
    {
      fprintf(stderr, "Actually at this stage, they should be equal\n");
      fprintf(stderr, "Please debug\n");
      exit(1);
    }
  nof_vertices=TheRec1.nbVertices;
  cl1=g1->canonical_form(stats, &report_aut_void, stderr);
  cl2=g2->canonical_form(stats, &report_aut_void, stderr);
  if ((clR1 = (unsigned int*)malloc(nof_vertices*sizeof(unsigned int))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<nof_vertices; i++)
    clR1[cl1[i]]=i;
  if ((clR2 = (unsigned int*)malloc(nof_vertices*sizeof(unsigned int))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<nof_vertices; i++)
    clR2[cl2[i]]=i;
  for (iVert=0; iVert<nof_vertices; iVert++)
    TheEquivExp.push_back(-1);
  for (iVert=0; iVert<nof_vertices; iVert++)
    fprintf(stderr, "iVert=%d cl1=%d cl2=%d\n", iVert, cl1[iVert], cl2[iVert]);
  for (iVert=0; iVert<nof_vertices; iVert++)
    TheEquivExp[iVert]=clR2[cl1[iVert]];
  free(clR1);
  free(clR2);
  for (iVert=0; iVert<nof_vertices; iVert++)
    {
      jVert=TheEquivExp[iVert];
      if (TheRec1.ListColor[iVert] != TheRec2.ListColor[jVert])
	{
	  FreeRecAdjMat(&TheRec1);
	  FreeRecAdjMat(&TheRec2);
	  delete g1;
	  delete g2;
	  fprintf(stderr, "TheReply=0 case 2\n");
	  *TheReply=0;
	  return;
	}
    }
  for (iVert1=0; iVert1<nof_vertices; iVert1++)
    {
      iVert2=TheEquivExp[iVert1];
      for (jVert1=0; jVert1<nof_vertices; jVert1++)
	{
	  jVert2=TheEquivExp[jVert1];
	  idx1=iVert1 + nof_vertices*jVert1;
	  idx2=iVert2 + nof_vertices*jVert2;
	  if (TheRec1.AdjMat[idx1] != TheRec2.AdjMat[idx2])
	    {
	      FreeRecAdjMat(&TheRec1);
	      FreeRecAdjMat(&TheRec2);
	      delete g1;
	      delete g2;
	      fprintf(stderr, "TheReply=0 case 3\n");
	      *TheReply=0;
	      return;
	    }
	}
    }
  TheEquiv.clear();
  nbRow=WMat1.nbRow;
  for (i=0; i<nbRow; i++)
    TheEquiv.push_back(TheEquivExp[i]);
  for (iVert1=0; iVert1<nbRow; iVert1++)
    {
      iVert2=TheEquiv[iVert1];
      for (jVert1=0; jVert1<nbRow; jVert1++)
	{
	  jVert2=TheEquiv[jVert1];
	  idx1=iVert1 + nbRow*jVert1;
	  idx2=iVert2 + nbRow*jVert2;
	  if (WMat1.TheMat[idx1] != WMat2.TheMat[idx2])
	    {
	      fprintf(stderr, "Our reduction technique is broken\n");
	      fprintf(stderr, "Please panic and debug\n");
	      exit(1);
	    }
	}
    }
  FreeRecAdjMat(&TheRec1);
  FreeRecAdjMat(&TheRec2);
  delete g1;
  delete g2;
  fprintf(stderr, "We find an isomorphism\n");
  *TheReply=1;
}
