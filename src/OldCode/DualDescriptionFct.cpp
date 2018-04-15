#include "DualDescriptionFct.h"

void VectVectInt_Magma_Print(FILE *f, VectVectInt &ListOrbit)
{
  int iOrbit, nbOrbit, siz;
  vector<int> eRepr;
  int i;
  nbOrbit=ListOrbit.size();
  fprintf(f, "[");
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    {
      if (iOrbit>0)
	fprintf(f, ",\n");
      eRepr=ListOrbit[iOrbit];
      siz=eRepr.size();
      fprintf(f, "[");
      for (i=0; i<siz; i++)
	{
	  if (i>0)
	    fprintf(f, ",");
	  fprintf(f, "%d", eRepr[i]);
	}
      fprintf(f, "]");
    }
  fprintf(f, "]\n");
}

void LISTORBIT_Print(FILE *f, EnumerationOrbit &ListOrbit)
{
  int iOrbit, nbOrbit, siz;
  vector<int> eRepr;
  nbOrbit=ListOrbit.ListRepresentative.size();
  fprintf(f, "nbOrbit=%d\n", nbOrbit);
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    {
      eRepr=ListOrbit.ListRepresentative[iOrbit];
      fprintf(f, " iOrbit=%d/%d", iOrbit, nbOrbit);
      siz=eRepr.size();
      fprintf(f, " inc=%d |O|=", siz);
      mpz_out_str(stderr, 10, ListOrbit.OrbitSizes[iOrbit].get_mpz_t());
      fprintf(f, " status=%d splt=%d\n", ListOrbit.ListStatus[iOrbit], ListOrbit.ListSplittable[iOrbit]);
    }
}

void LISTORBIT_FuncInsert(EnumerationOrbit &ListOrbit, 
			  TheGroupFormat &TheGRP, 
			  TheHeuristic *TheHeuSplitting, 
			  vector<int> eList, 
			  int TheLevel)
{
  TheGroupFormat TheStab;
  int sizGRP, i;
  int iOrbit, nbOrbit;
  TheEntries *TheCand;
  mpz_class tcla, sizGRP_z;
  mpz_t q_z;
  int sizInc, IsSplittable;
  vector<int> TheMin1, TheMin2;
  nbOrbit=ListOrbit.ListRepresentative.size();
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    if (ListOrbit.ListRepresentative[iOrbit] == eList)
      return;
  mpz_init(q_z);
  fprintf(stderr, "FuncInsert, step 1\n");
  if ((TheCand = (TheEntries*)malloc(sizeof(TheEntries))) == 0)
    exit (EXIT_FAILURE);
  TheCand->nbEntry=3;
  if ((TheCand->ListValue = (int*)malloc(3*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheCand->ListEntry = (char**)malloc(3*sizeof(char*))) == 0)
    exit (EXIT_FAILURE);
  for (i=0; i<3; i++)
    if ((TheCand->ListEntry[i] = (char*)malloc(2*sizeof(char))) == 0)
      exit (EXIT_FAILURE);
  GetStabilizer(TheGRP, eList, TheStab);
  strcpy(TheCand->ListEntry[0], "g");
  strcpy(TheCand->ListEntry[1], "i");
  strcpy(TheCand->ListEntry[2], "l");
  sizGRP_z=TheStab.group->order<mpz_class>();
  mpz_divexact(q_z, ListOrbit.GRPsize, sizGRP_z.get_mpz_t());
  tcla=mpz_class(q_z);
  sizInc=eList.size();
  sizGRP=TheStab.group->order();
  TheCand->ListValue[0]=sizGRP;
  TheCand->ListValue[1]=sizInc;
  TheCand->ListValue[2]=TheLevel;
  IsSplittable=HeuristicEvaluation(TheCand, TheHeuSplitting);
  fprintf(stderr, "TheLevel=%d sizInc=%d sizGRP_z=", TheLevel, sizInc);
  mpz_out_str(stderr, 10, sizGRP_z.get_mpz_t());
  fprintf(stderr, " splt=%d", IsSplittable);
  fprintf(stderr, "  |O|=");
  mpz_out_str(stderr, 10, q_z);
  fprintf(stderr, "\n");

  ListOrbit.ListRepresentative.push_back(eList);
  ListOrbit.ListStatus.push_back(1);
  ListOrbit.ListSplittable.push_back(IsSplittable);
  ListOrbit.OrbitSizes.push_back(tcla);
  for (i=0; i<3; i++)
    free(TheCand->ListEntry[i]);
  free(TheCand->ListEntry);
  free(TheCand->ListValue);
  free(TheCand);
  mpz_clear(q_z);
  fprintf(stderr, "FuncInsert, step 8\n");
}

void LISTORBIT_Reordering(EnumerationOrbit &ListOrbit)
{
  int iOrbit, jOrbit, nbOrbit;
  mpz_t eOrbSize;
  int sizI, sizJ, WeChange;
  int eStatus, eSplit, test;
  vector<int> eRepr;
  nbOrbit=ListOrbit.ListRepresentative.size();
  mpz_init(eOrbSize);
  for (iOrbit=0; iOrbit<nbOrbit-1; iOrbit++)
    for (jOrbit=iOrbit+1; jOrbit<nbOrbit; jOrbit++)
      {
	sizI=ListOrbit.ListRepresentative[iOrbit].size();
	sizJ=ListOrbit.ListRepresentative[jOrbit].size();
	WeChange=0;
	if (sizI > sizJ)
	  WeChange=1;
	else
	  {
	    if (sizI == sizJ)
	      {
		test=mpz_cmp(ListOrbit.OrbitSizes[iOrbit].get_mpz_t(), ListOrbit.OrbitSizes[jOrbit].get_mpz_t());
		if (test == 1)
		  WeChange=1;
	      }
	  }
	if (WeChange == 1)
	  {
	    eRepr=ListOrbit.ListRepresentative[iOrbit];
	    ListOrbit.ListRepresentative[iOrbit]=ListOrbit.ListRepresentative[jOrbit];
	    ListOrbit.ListRepresentative[jOrbit]=eRepr;

	    eStatus=ListOrbit.ListStatus[iOrbit];
	    ListOrbit.ListStatus[iOrbit]=ListOrbit.ListStatus[jOrbit];
	    ListOrbit.ListStatus[jOrbit]=eStatus;

	    eSplit=ListOrbit.ListSplittable[iOrbit];
	    ListOrbit.ListSplittable[iOrbit]=ListOrbit.ListSplittable[jOrbit];
	    ListOrbit.ListSplittable[jOrbit]=eSplit;
	    
	    mpz_set(eOrbSize, ListOrbit.OrbitSizes[iOrbit].get_mpz_t());
	    ListOrbit.OrbitSizes[iOrbit]=ListOrbit.OrbitSizes[jOrbit];
	    ListOrbit.OrbitSizes[jOrbit]=mpz_class(eOrbSize);
	  }
      }
  mpz_clear(eOrbSize);
}

int LISTORBIT_IsEnumerationFinished(EnumerationOrbit &ListOrbit, 
				    int LastNeeded)
{
  int nbDone, nbOrbit, iOrbit;
  nbOrbit=ListOrbit.ListRepresentative.size();
  for (iOrbit=0; iOrbit<LastNeeded; iOrbit++)
    if (ListOrbit.ListStatus[iOrbit] == 1)
      return 0;
  nbDone=0;
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    if (ListOrbit.ListStatus[iOrbit] == 0)
      nbDone++;
  if (nbDone == 0)
    return 0;
  return 1;
}

void DUALDESC_DetermineComputeNeeded(EnumerationOrbit &ListOrbit, 
				     TheGroupFormat &TheGRP, 
				     int *LastNeeded)
{
  int *gList, *rList;
  int *intList;
  mpz_t nbUndone, MaxAllowed;
  int iOrbit, nbOrbit, iVert, nbVert, i, siz;
  int eSum, nbDone;
  vector<int> eRepr;
  int test;
  nbOrbit=ListOrbit.ListRepresentative.size();
  nbDone=0;
  for (iOrbit=nbOrbit-1; iOrbit>=0; iOrbit--)
    if (ListOrbit.ListStatus[iOrbit] == 0)
      nbDone++;
  if (nbDone == 0)
    {
      *LastNeeded=1;
      return;
    }
  nbVert=ListOrbit.nbVert;
  if ((intList = (int*)malloc(nbVert*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((gList = (int*)malloc(nbVert*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((rList = (int*)malloc(nbVert*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (iVert=0; iVert<nbVert; iVert++)
    intList[iVert]=1;
  mpz_init(nbUndone);
  mpz_init(MaxAllowed);
  mpz_set_si(nbUndone, 0);
  mpz_set_si(MaxAllowed, ListOrbit.MaxAllowedUndone);
  for (iOrbit=nbOrbit-1; iOrbit>=0; iOrbit--)
    if (ListOrbit.ListStatus[iOrbit] == 1)
      {
	for (iVert=0; iVert<nbVert; iVert++)
	  gList[iVert]=0;
	eRepr=ListOrbit.ListRepresentative[iOrbit];
	siz=eRepr.size();
	for (i=0; i<siz; i++)
	  gList[eRepr[i]]=1;
	OrbitIntersection(TheGRP, gList, rList);
	for (iVert=0; iVert<nbVert; iVert++)
	  intList[iVert]=intList[iVert]*rList[iVert];
	eSum=0;
	for (iVert=0; iVert<nbVert; iVert++)
	  eSum = eSum + intList[iVert];
	mpz_add(nbUndone, nbUndone, ListOrbit.OrbitSizes[iOrbit].get_mpz_t());
	test=mpz_cmp(nbUndone, MaxAllowed);
	if (test > 0 && eSum == 0)
	  {
	    mpz_clear(nbUndone);
	    mpz_clear(MaxAllowed);
	    free(intList);
	    free(gList);
	    free(rList);
	    *LastNeeded=iOrbit+1;
	    return;
	  }
      }
  *LastNeeded=-1;
  free(intList);
  free(gList);
  free(rList);
  mpz_clear(nbUndone);
  mpz_clear(MaxAllowed);
}

void DUALDESC_MPI_ProcDirectSplittingOrbit(EnumerationOrbit &ListOrbit, 
					   int LastNeeded, 
					   int NbProc, 
					   vector<VectVectInt> &ThePartition,
					   int *totalOrbitAssigned)
{
  VectVectInt eListOrbit;
  int iProc, iOrbit;
  vector<int> eRepr;
  ThePartition.clear();
  eListOrbit.clear();
  for (iProc=0; iProc<NbProc; iProc++)
    ThePartition.push_back(eListOrbit);
  iProc=0;
  *totalOrbitAssigned=0;
  for (iOrbit=0; iOrbit<LastNeeded; iOrbit++)
    if (ListOrbit.ListStatus[iOrbit] == 1 && ListOrbit.ListSplittable[iOrbit] == 0)
      {
	ListOrbit.ListStatus[iOrbit]=0;
	eRepr=ListOrbit.ListRepresentative[iOrbit];
	ThePartition[iProc].push_back(eRepr);
	iProc++;
	*totalOrbitAssigned=*totalOrbitAssigned+1;
	if (iProc == NbProc)
	  iProc=0;
      }
}

void DUALDESC_MPI_ProcIndirectSplittingOrbit(EnumerationOrbit &ListOrbit, 
					     int LastNeeded, 
					     int NbProc, 
					     int *ListProc,
					     vector<int> &ListOrbitCons,
					     vector<vector<int> > &ListListProc)
{
  int iProc, jProc, iOrbit, jOrbit;
  set<int> iOrbitList;
  set<int>::iterator iter;
  int *ListProcessorAssigned;
  vector<int> eNewList;
  set<int> ListCons;
  iOrbitList.clear();
  if ((ListProcessorAssigned = (int*)malloc(NbProc*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (iProc=0; iProc<NbProc; iProc++)
    ListProcessorAssigned[iProc]=0;
  ListOrbitCons.clear();
  ListListProc.clear();
  iProc=0;
  while(1)
    {
      for (iOrbit=0; iOrbit<LastNeeded; iOrbit++)
	if (ListOrbit.ListStatus[iOrbit] == 1 && ListOrbit.ListSplittable[iOrbit] == 1)
	  {
	    ListCons.insert(iOrbit);
	    ListProcessorAssigned[iProc]=iOrbit;
	    iOrbitList.insert(iOrbit);
	    iProc++;
	    if (iProc == NbProc)
	      {
		iter=ListCons.begin();
		while(iter != ListCons.end())
		  {
		    jOrbit=*iter;
		    ListOrbit.ListStatus[jOrbit]=0;
		    iter++;
		  }
		iter=iOrbitList.begin();
		while(iter != iOrbitList.end())
		  {
		    jOrbit=*iter;
		    eNewList.clear();
		    for (jProc=0; jProc<NbProc; jProc++)
		      if (ListProcessorAssigned[jProc] == jOrbit)
			eNewList.push_back(ListProc[jProc]);
		    ListOrbitCons.push_back(jOrbit);
		    ListListProc.push_back(eNewList);
		    iter++;
		  }
		free(ListProcessorAssigned);
		return;
	      }
	  }
    }
}

void DUALDESC_DirectComputation(MyMatrix *TheEXT,
				TheGroupFormat &GRP, 
				VectVectInt &TheInput,
				VectVectInt &TheOutput)
{
  int siz, i, iOrb, nbOrb, WeFound;
  vector<int> eRepr, eFlip;
  MyMatrix *TheEXTred1, *TheFACred1;
  VectVectInt TheReturnRed1;
  vector<int> eFlipCan;
  
  TheGroupFormat TheStab, TheStabRed;
  if ((TheEXTred1 = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheFACred1 = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  TheOutput.clear();
  siz=TheInput.size();
  for (i=0; i<siz; i++)
    {
      eRepr=TheInput[i];
      MatrixReduction(TheEXT, eRepr, TheEXTred1);
      CDD_DualDescription(TheEXTred1, TheFACred1);
      GetStabilizer(GRP, eRepr, TheStab);
      TheStabRed=ReducedGroupAction(TheStab, eRepr);
      OrbitSplitting(TheEXTred1, TheFACred1, TheStabRed, TheReturnRed1);
      nbOrb=TheReturnRed1.size();
      for (iOrb=0; iOrb<nbOrb; iOrb++)
	{
	  TestFacetness(TheEXTred1, TheReturnRed1[iOrb]);
	  eFlip=ComputeFlipping(TheEXT, eRepr, TheReturnRed1[iOrb]);
	  eFlipCan=GROUP_Canonicalization(GRP, eFlip);
	  WeFound=0;
	  siz=TheOutput.size();
	  for (i=0; i<siz; i++)
	    if (WeFound == 0)
	      if (TheOutput[i] == eFlipCan)
		WeFound=1;
	  if (WeFound == 0)
	    TheOutput.push_back(eFlipCan);
	  TestFacetness(TheEXT, eFlip);
	}
      MATRIX_Free(TheEXTred1);
      MATRIX_Free(TheFACred1);
    }
  free(TheEXTred1);
  free(TheFACred1);
}

void LISTORBIT_InquireBank(EnumerationOrbit &ListOrbit,
			   MyMatrix *EXT, 
			   TheGroupFormat &GRP, 
			   TheHeuristic *TheHeuSplitting, 
			   int TheLevel, 
			   int iOrbit, 
			   MPI_Comm comm)
{
  int rank, dest, ret, iCos, nbCos, iOrb;
  int TheReply, nbOrb;
  vector<int> eListI, eFlip, eRepr;
  VectVectInt ListListSet, ListRepr;
  WeightMatrix WMat;
  TheGroupFormat TheGRPbig, TheStab, TheStabRed;
  int *eVectOrig, *eVectAction, *eVectReply;
  MyMatrix *EXTred;
  vector<int> eFlipCan;
  MPI_Status status;
  ret=MPI_Comm_rank(comm, &rank);

  dest=0;
	    
  if ((EXTred = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectOrig = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectReply = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectAction = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  eVectOrig[0]=rank;
  ret=MPI_Send(eVectOrig, 1, MPI_INT, dest, 1023, comm);
  
  eVectAction[0]=2;
  ret=MPI_Send(eVectAction, 1, MPI_INT, dest, 1492, comm);

  eListI=ListOrbit.ListRepresentative[iOrbit];
  MatrixReduction(EXT, eListI, EXTred);

  MPI_SEND_MyMatrix(EXTred, dest, 242, comm);
  
  MPI_Recv(eVectReply, 1, MPI_INT, dest, 1897, comm, &status);
  TheReply=eVectReply[0];
  if (TheReply == 1)
    {
      ListOrbit.ListStatus[iOrbit]=0;
      MPI_RECV_VectVectInt(ListRepr, dest, 1917, comm);
      GetWeightMatrix(WMat, EXTred);
      GetStabilizerWeightMatrix(WMat, TheGRPbig);
      free(WMat.TheMat);
      WMat.ListWeight.clear();

      GetStabilizer(GRP, eListI, TheStab);
      TheStabRed=ReducedGroupAction(TheStab, eListI);

      nbOrb=ListRepr.size();
      for (iOrb=0; iOrb<nbOrb; iOrb++)
	{
	  eRepr=ListRepr[iOrb];
	  TestFacetness(EXTred, eRepr);
	  DoubleCosetDescription(TheGRPbig, TheStabRed, eRepr, ListListSet);
	  nbCos=ListListSet.size();
	  for (iCos=0; iCos<nbCos; iCos++)
	    {
	      fprintf(stderr, "Before ComputeFlipping\n");
	      eFlip=ComputeFlipping(EXT, eListI, ListListSet[iCos]);
	      fprintf(stderr, " After ComputeFlipping\n");
	      eFlipCan=GROUP_Canonicalization(GRP, eFlip);
	      TestFacetness(EXT, eFlipCan);
	      LISTORBIT_FuncInsert(ListOrbit, GRP, TheHeuSplitting, 
				   eFlipCan, TheLevel);
	    }
	}
    }
  MATRIX_Free(EXTred);
  free(EXTred);
  free(eVectOrig);
  free(eVectAction);
  free(eVectReply);
}


void DUALDESC_MPI_SelectComputation(MyMatrix *EXT, 
				    TheGroupFormat &GRP, 
				    TheHeuristic *TheHeuSplitting, 
				    int TheLevel, 
				    int NbProc, 
				    int *ListProc, 
				    VectVectInt &TheOutput, 
				    MPI_Comm comm)
{
  vector<int> eInc;
  EnumerationOrbit ListOrbit;
  int LastNeeded;
  vector<VectVectInt> ThePartition;
  vector<VectVectInt> TheResult;
  vector<int> eIncCan, eFlipCan;
  int rank, src, dest, ret;
  int *ListProcSel;
  int *eVectOrig, *eVectAction, *eVectLevel, *eVectNbProc;
  int iOrbit, nbOrbit, i, siz, iOrb, iCos, nbCos;
  int iOrbitCons, nbOrbitCons;
  int iProc, NbProcSel, test, nbRow, nbCol;
  vector<int> eFlip, eReprFac, eRepr;
  set<int>::iterator iter;
  VectVectInt ListListSet;
  TheGroupFormat TheGRPbig, TheStab, TheStabRed;
  mpz_class tcla;
  WeightMatrix WMat;
  int totalOrbitAssigned;
  VectVectInt ListListProc, eListListO;
  vector<int> ListOrbitCons;
  int ord1, ord2;
  fprintf(stderr, "SelectComputation, step 1\n");
  ret=MPI_Comm_rank(comm, &rank);
  if (ListProc[0] != rank)
    {
      fprintf(stderr, "We fail one of the condition of the command\n");
      exit(1);
    }
  fprintf(stderr, "SelectComputation, step 2\n");
  GetWeightMatrix(WMat, EXT);
  fprintf(stderr, "SelectComputation, step 3\n");
  GetStabilizerWeightMatrix(WMat, TheGRPbig);
  fprintf(stderr, "SelectComputation, step 4\n");
  free(WMat.TheMat);
  fprintf(stderr, "SelectComputation, step 5\n");
  WMat.ListWeight.clear();
  fprintf(stderr, "SelectComputation, step 6\n");
  nbRow=EXT->nbRow;
  nbCol=EXT->nbCol;
  fprintf(stderr, "nbRow=%d nbCol=%d\n", nbRow, nbCol);
  fprintf(stderr, "SelectComputation, step 7\n");
  ListOrbit.nbVert=nbRow;
  ListOrbit.MaxAllowedUndone=nbCol-2;
  fprintf(stderr, "SelectComputation, step 8\n");
  tcla=TheGRPbig.group->order<mpz_class>();
  ord1=GRP.group->order();
  ord2=TheGRPbig.group->order();
  fprintf(stderr, "ord1=%d ord2=%d\n", ord1, ord2);
  fprintf(stderr, "SelectComputation, step 9\n");
  mpz_init(ListOrbit.GRPsize);
  mpz_set(ListOrbit.GRPsize, tcla.get_mpz_t());
  fprintf(stderr, "SelectComputation, step 10\n");
  FindOneInitialVertex(EXT, eInc);
  fprintf(stderr, "SelectComputation, step 10.1\n");
  eIncCan=GROUP_Canonicalization(TheGRPbig, eInc);
  fprintf(stderr, "SelectComputation, step 11\n");
  if (EXT->nbRow != GRP.n)
    {
      fprintf(stderr, "We have a clear error here\n");
      exit(1);
    }
  fprintf(stderr, "SelectComputation, step 12\n");
  LISTORBIT_FuncInsert(ListOrbit, TheGRPbig, TheHeuSplitting, 
		       eIncCan, TheLevel);
  fprintf(stderr, "SelectComputation, step 13\n");
  if ((eVectOrig = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectLevel = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectNbProc = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectAction = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  while(1)
    {
      DUALDESC_DetermineComputeNeeded(ListOrbit, TheGRPbig, &LastNeeded);
      fprintf(stderr, "After DetermineComputeNeeded\n");
      test=LISTORBIT_IsEnumerationFinished(ListOrbit, LastNeeded);
      fprintf(stderr, "After IsEnumerationFinished\n");
      if (test == 1)
	break;
      for (iOrbit=0; iOrbit<LastNeeded; iOrbit++)
	if (ListOrbit.ListStatus[iOrbit] == 1 && ListOrbit.ListSplittable[iOrbit] == 1)
	  LISTORBIT_InquireBank(ListOrbit, EXT, TheGRPbig, 
				TheHeuSplitting, TheLevel, iOrbit, comm);
      fprintf(stderr, "After LISTORBIT_InquireBank\n");
      DUALDESC_MPI_ProcDirectSplittingOrbit(ListOrbit, LastNeeded, NbProc, 
					    ThePartition, &totalOrbitAssigned);
      fprintf(stderr, "After MPI_ProcDirectSplittingOrbit\n");
      if (totalOrbitAssigned > 0)
	{
	  for (iProc=1; iProc<NbProc; iProc++)
	    {
	      dest=ListProc[iProc];

	      eVectOrig[0]=rank;
	      ret=MPI_Send(eVectOrig, 1, MPI_INT, dest, 1023, comm);

	      eVectAction[0]=2;
	      ret=MPI_Send(eVectAction, 1, MPI_INT, dest, 1492, comm);

	      MPI_SEND_MyMatrix(EXT, dest, 242, comm);
	      MPI_SEND_Group(TheGRPbig, dest, 447, comm);
	      MPI_SEND_VectVectInt(ThePartition[iProc], dest, 2047, comm);
	    }
	  TheResult.clear();
	  DUALDESC_DirectComputation(EXT, TheGRPbig, 
				     ThePartition[0], eListListO);
	  TheResult.push_back(eListListO);
	  for (iProc=1; iProc<NbProc; iProc++)
	    {
	      src=ListProc[iProc];
	      MPI_RECV_VectVectInt(eListListO, src, 1815, comm);
	      TheResult.push_back(eListListO);
	    }
	  for (iProc=0; iProc<NbProc; iProc++)
	    {
	      siz=TheResult[iProc].size();
	      for (iOrb=0; iOrb<siz; iOrb++)
		{
		  eRepr=TheResult[iProc][iOrb];
		  LISTORBIT_FuncInsert(ListOrbit, TheGRPbig, TheHeuSplitting, 
				       eRepr, TheLevel);
		}
	    }
	  fprintf(stderr, "After LISTORBIT_FuncInsert(s)\n");
	}
      else
	{
	  DUALDESC_MPI_ProcIndirectSplittingOrbit(ListOrbit, LastNeeded, 
						  NbProc, ListProc,
						  ListOrbitCons, ListListProc);
	  fprintf(stderr, "After MPI_ProcIndirectSplittingOrbit\n");
	  TheResult.clear();
	  nbOrbitCons=ListOrbitCons.size();
	  fprintf(stderr, "nbOrbitCons=%d\n", nbOrbitCons);
	  for (iOrbitCons=nbOrbitCons-1; iOrbitCons>=0; iOrbitCons--)
	    {
	      fprintf(stderr, "iOrbitCons=%d/%d\n", iOrbitCons, nbOrbitCons);
	      iOrbit=ListOrbitCons[iOrbitCons];
	      NbProcSel=ListListProc[iOrbitCons].size();
	      if ((ListProcSel = (int*)malloc(NbProcSel*sizeof(int))) == 0)
		exit (EXIT_FAILURE);
	      for (i=0;i<NbProcSel; i++)
		ListProcSel[i]=ListListProc[iOrbitCons][i];
	      eRepr=ListOrbit.ListRepresentative[iOrbit];
	      fprintf(stderr, "We have eRepr\n");

	      if (iOrbitCons == 0)
		{
		  fprintf(stderr, "Before Call to SelectAndFlip\n");
		  DUALDESC_MPI_SelectAndFlipComputation(EXT, eRepr, TheGRPbig,
							TheHeuSplitting, 
							TheLevel+1, NbProcSel, 
							ListProcSel, 
							TheOutput, comm);
		  fprintf(stderr, "After Call to SelectAndFlip\n");
		  TheResult.push_back(TheOutput);
		}
	      else
		{
		  fprintf(stderr, "Before Call to SelectAndFlip EXPORT\n");
		  dest=ListProcSel[0];
		  
		  eVectOrig[0]=rank;
		  ret=MPI_Send(eVectOrig, 1, MPI_INT, dest, 1023, comm);
		  
		  eVectAction[0]=3;
		  ret=MPI_Send(eVectAction, 1, MPI_INT, dest, 1492, comm);
		  
		  eVectLevel[0]=TheLevel+1;
		  ret=MPI_Send(eVectLevel, 1, MPI_INT, dest, 142, comm);
		  
		  eVectNbProc[0]=NbProcSel;
		  ret=MPI_Send(eVectNbProc, 1, MPI_INT, dest, 143, comm);
		  
		  ret=MPI_Send(ListProcSel, NbProcSel, MPI_INT, dest, 144, comm);

		  MPI_SEND_VectInt(eRepr, dest, 145, comm);
		  
		  MPI_SEND_MyMatrix(EXT, dest, 243, comm);
		  MPI_SEND_Group(TheGRPbig, dest, 448, comm);
		  fprintf(stderr, "After Call to SelectAndFlip EXPORT\n");
		}
	      free(ListProcSel);
	    }
	  fprintf(stderr, "Now getting back the results\n");
	  for (iOrbitCons=1; iOrbitCons<nbOrbitCons; iOrbitCons++)
	    {
	      src=ListListProc[iOrbitCons][0];
	      MPI_RECV_VectVectInt(eListListO, src, 1830, comm);
	      siz=eListListO.size();
	      TheResult.push_back(eListListO);
	    }
	  for (iOrbitCons=0; iOrbitCons<nbOrbitCons; iOrbitCons++)
	    {
	      siz=TheResult[iOrbitCons].size();
	      for (iOrb=0; iOrb<siz; iOrb++)
		LISTORBIT_FuncInsert(ListOrbit, TheGRPbig,
				     TheHeuSplitting, 
				     TheResult[iOrbitCons][iOrb], TheLevel);
	    }
	}
      fprintf(stderr, "Before reordering\n");
      LISTORBIT_Reordering(ListOrbit);
      fprintf(stderr, "After reordering\n");
    }
  nbOrbit=ListOrbit.ListRepresentative.size();
  fprintf(stderr, "We found nbOrbit=%d\n", nbOrbit);
  TheOutput.clear();
  for (iOrbit=0; iOrbit<nbOrbit; iOrbit++)
    {
      DoubleCosetDescription(TheGRPbig, GRP, ListOrbit.ListRepresentative[iOrbit], ListListSet);
      nbCos=ListListSet.size();
      for (iCos=0; iCos<nbCos; iCos++)
	TheOutput.push_back(ListListSet[iCos]);
    }
  fprintf(stderr, "After DoubleCoset operations\n");
  free(eVectOrig);
  free(eVectAction);
  free(eVectLevel);
  free(eVectNbProc);
  mpz_clear(ListOrbit.GRPsize);
}

void DUALDESC_MPI_SelectAndFlipComputation(MyMatrix *EXT, 
					   vector<int> eSelect, 
					   TheGroupFormat &GRP, 
					   TheHeuristic *TheHeuSplitting, 
					   int TheLevel, 
					   int NbProc, 
					   int *ListProc, 
					   VectVectInt &TheOutput, 
					   MPI_Comm comm)
{
  MyMatrix *TheEXTred1;
  vector<int> eReprFac, eFlip, eFlipCan;
  int len, i, lenO, j, WeFound;
  VectVectInt ThePreOutput;
  TheGroupFormat TheStab, TheStabRed;
  if ((TheEXTred1 = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  fprintf(stderr, "SelectAndFlip, step 1\n");
  MatrixReduction(EXT, eSelect, TheEXTred1);
  fprintf(stderr, "SelectAndFlip, step 2\n");
  GetStabilizer(GRP, eSelect, TheStab);
  fprintf(stderr, "SelectAndFlip, step 3\n");
  TheStabRed=ReducedGroupAction(TheStab, eSelect);
  fprintf(stderr, "SelectAndFlip, step 4\n");
  DUALDESC_MPI_SelectComputation(TheEXTred1, TheStabRed, TheHeuSplitting, 
				 TheLevel, NbProc, ListProc, 
				 ThePreOutput, comm);
  fprintf(stderr, "SelectAndFlip: After MPI_SelectComputation\n");
  TheOutput.clear();
  len=ThePreOutput.size();
  fprintf(stderr, "SelectAndFlip, step 5\n");
  for (i=0; i<len; i++)
    {
      fprintf(stderr, "i=%d/%d\n", i, len);
      eReprFac=ThePreOutput[i];
      fprintf(stderr, "Iterative, step 1\n");
      eFlip=ComputeFlipping(EXT, eSelect, eReprFac);
      fprintf(stderr, "Iterative, step 2\n");
      TestFacetness(EXT, eFlip);
      fprintf(stderr, "Iterative, step 3\n");
      eFlipCan=GROUP_Canonicalization(GRP, eFlip);
      fprintf(stderr, "Iterative, step 4\n");
      lenO=TheOutput.size();
      fprintf(stderr, "Iterative, step 5\n");
      WeFound=0;
      fprintf(stderr, "Iterative, step 6\n");
      for (j=0; j<lenO; j++)
	if (WeFound == 0)
	  if (TheOutput[j] == eFlipCan)
	    WeFound=1;
      fprintf(stderr, "Iterative, step 7\n");
      if (WeFound == 0)
	TheOutput.push_back(eFlipCan);
      fprintf(stderr, "Iterative, step 8\n");
    }
  fprintf(stderr, "SelectAndFlip: After the flipping and canonicalize\n");
  MATRIX_Free(TheEXTred1);
  free(TheEXTred1);
}

void DUALDESC_MPI_AwaitingOrders(TheHeuristic *TheHeuSplitting, 
				 MPI_Comm comm)
{
  int TheSender;
  int *eVectOrig, *eVectAction, *eVectLevel, *eVectNbProc;
  MyMatrix *TheEXTask;
  TheGroupFormat TheGRP;
  int TheAct, TheLevel, NbProc;
  VectVectInt eListListI, eListListO;
  MPI_Status status;
  vector<int> eListWork;
  int *ListProc;
  if ((eVectOrig = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectLevel = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectAction = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((eVectNbProc = (int*)malloc(1*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheEXTask = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  while(1)
    {
      MPI_Recv(eVectOrig, 1, MPI_INT, MPI_ANY_SOURCE, 1023, comm, &status);
      TheSender=eVectOrig[0];

      MPI_Recv(eVectAction, 1, MPI_INT, TheSender, 1492, comm, &status);
      TheAct=eVectAction[0];
      if (TheAct == 1)
	break;
      if (TheAct == 2)
	{
	  MPI_RECV_MyMatrix(TheEXTask, TheSender, 242, comm);
	  MPI_RECV_Group(TheGRP, TheSender, 447, comm);
	  MPI_RECV_VectVectInt(eListListI, TheSender, 2047, comm);
	  DUALDESC_DirectComputation(TheEXTask,
				     TheGRP, 
				     eListListI, 
				     eListListO);
	  fprintf(stderr, "After DirectComputation\n");
	  MPI_SEND_VectVectInt(eListListO, TheSender, 1815, comm);
	  fprintf(stderr, "After MPI_SEND_VectVectInt\n");
	  MATRIX_Free(TheEXTask);
	}
      if (TheAct == 3)
	{
	  MPI_Recv(eVectLevel, 1, MPI_INT, TheSender, 142, comm, &status);
	  TheLevel=eVectLevel[0];

	  MPI_Recv(eVectNbProc, 1, MPI_INT, TheSender, 143, comm, &status);
	  NbProc=eVectNbProc[0];
	  
	  if ((ListProc = (int*)malloc(NbProc*sizeof(int))) == 0)
	    exit (EXIT_FAILURE);
	  MPI_Recv(ListProc, NbProc, MPI_INT, TheSender, 144, comm, &status);

	  MPI_RECV_VectInt(eListWork, TheSender, 145, comm);
	  
	  MPI_RECV_MyMatrix(TheEXTask, TheSender, 243, comm);
	  MPI_RECV_Group(TheGRP, TheSender, 448, comm);

	  DUALDESC_MPI_SelectAndFlipComputation(TheEXTask, eListWork, 
						TheGRP, TheHeuSplitting, 
						TheLevel, NbProc, ListProc, 
						eListListO, comm);
	  fprintf(stderr, "After MPI_SelectAndFlipComputation\n");
	  free(ListProc);
	  MPI_SEND_VectVectInt(eListListO, TheSender, 1830, comm);
	  fprintf(stderr, "After MPI_SEND_VectVectInt\n");
	  MATRIX_Free(TheEXTask);
	}
    }
  free(eVectOrig);
  free(eVectAction);
  free(eVectLevel);
  free(TheEXTask);
}
