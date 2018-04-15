#include "DualDescriptionFct.h"
#include "BankingSystem.h"
#include "MPQ_CDD_LRS.h"
#include "PolytopeFct.h"
#include "PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  TheGroupFormat TheGRP;
  FILE *DATAEXT=NULL;
  FILE *DATAGRP=NULL;
  MyMatrix *TheEXT;
  vector<vector<int> > TheReturn;
  MPI_Comm comm;
  int rank, NbProc;
  int *ListProc;
  int iret;
  BankingDescription TheBank;
  int jProcCall, jProc, iProc, ret, TheLevel;
  int *eVectAction, *eVectOrig;
  VectVectInt TheOutput;
  TheHeuristic *TheHeuSplitting;
  if ((TheHeuSplitting = (TheHeuristic*)malloc(sizeof(TheHeuristic))) == 0)
    exit (EXIT_FAILURE);

  Initialize_All();
  if ((TheEXT = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if (argc != 3)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneFacPolytope_MPI [DATAEXT] [DATAGRP]\n");
      fprintf(stderr, "DATAEXT: The input data of the polytope vertices\n");
      fprintf(stderr, "DATAGRP: The group for which we want orbit splitting\n");
      fprintf(stderr, "n\n");
      return -1;
    }
  fprintf(stderr, "Reading input\n");
  DATAEXT=fopen(argv[1], "r");
  if (DATAEXT == NULL)
    {
      fprintf(stderr,"input: The file %s was not found\n",argv[1]);
      return -1;
    }
  fprintf(stderr, "  input file:%s\n", argv[1]);
  ReadPolytope(DATAEXT, TheEXT);
  fprintf(stderr, "We have the polytope\n");

  DATAGRP=fopen(argv[2], "r");
  if (DATAGRP == NULL)
    {
      fprintf(stderr,"input: The file %s was not found\n",argv[2]);
      return -1;
    }
  fprintf(stderr, "  input file:%s\n", argv[2]);
  ReadGroup(TheGRP, DATAGRP);
  fprintf(stderr, "We have the group\n");

  comm=MPI_COMM_WORLD;
  ret=MPI_Init(&argc, &argv);
  ret=MPI_Comm_rank(comm, &rank);
  ret=MPI_Comm_size(comm, &NbProc);
  StandardHeuristicADM(TheHeuSplitting, comm);

  if (NbProc < 2)
    {
      fprintf(stderr, "The program needs at least two processors\n");
      ret=MPI_Finalize();
      return 1;
    }
  if ((ListProc = (int*)malloc((NbProc-1)*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  for (iProc=1; iProc<NbProc; iProc++)
    ListProc[iProc-1]=iProc;
  
  if (rank == 0)
    BANK_OperatingPolyhedralBank(TheBank, comm);
  if (rank == 1)
    {
      TheLevel=0;
      DUALDESC_MPI_SelectComputation(TheEXT, TheGRP, TheHeuSplitting, 
				     TheLevel, NbProc-1, ListProc, 
				     TheOutput, comm);
      VectVectInt_Magma_Print(stdout, TheOutput);
      if ((eVectOrig = (int*)malloc(1*sizeof(int))) == 0)
	exit (EXIT_FAILURE);
      if ((eVectAction = (int*)malloc(1*sizeof(int))) == 0)
	exit (EXIT_FAILURE);
      eVectOrig[0]=1;
      eVectAction[0]=1;
      for (jProc=1; jProc<NbProc; jProc++)
	{
	  if (jProc == 1)
	    jProcCall=0;
	  else
	    jProcCall=jProc;
	  ret=MPI_Send(eVectOrig, 1, MPI_INT, jProcCall, 1023, comm);
	  ret=MPI_Send(eVectAction, 1, MPI_INT, jProcCall, 1492, comm);
	}
      free(eVectAction);
      free(eVectOrig);
    }
  if (rank > 1)
    DUALDESC_MPI_AwaitingOrders(TheHeuSplitting, comm);
  
  MATRIX_Free(TheEXT);
  free(TheEXT);
  free(ListProc);
  HEURISTIC_Free(TheHeuSplitting);
  free(TheHeuSplitting);
  fprintf(stderr, "Completion of the program\n");
  Free_All();
  ret=MPI_Finalize();
  iret = fclose(DATAEXT);
  iret = fclose(DATAGRP);
}
