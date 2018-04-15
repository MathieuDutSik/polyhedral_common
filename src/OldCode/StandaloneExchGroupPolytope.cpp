#include "MPQ_CDD_LRS.h"
#include "PolytopeFct.h"
#include "PolytopeEquiStab.h"
#include "MPI_functions.h"

int main(int argc, char *argv[])
{
  TheGroupFormat TheGRP;
  FILE *DATAINP=NULL;
  MyMatrix *TheEXT;
  MPI_Comm comm;
  WeightMatrix WMat;
  int tagEXT, tagGRP;
  int iret, ret, rank, NbProc, iProc;
  if ((TheEXT = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "mpirun -np N StandaloneExchGroupPolytope [DATAINP]\n");
      fprintf(stderr, "DATAINP: a polytope");
      fprintf(stderr, "n\n");
      return -1;
    }
  tagEXT=100;
  tagGRP=123;
  comm=MPI_COMM_WORLD;
  ret=MPI_Init(&argc, &argv);
  ret=MPI_Comm_rank(comm, &rank);
  ret=MPI_Comm_size(comm, &NbProc);
  fprintf(stderr, "rank=%d NbProc=%d\n", rank, NbProc);
  if (rank == 0)
    {
      fprintf(stderr, "Reading input\n");
      DATAINP=fopen(argv[1], "r");
      if (DATAINP == NULL)
	{
	  fprintf(stderr,"input: The file %s was not found\n",argv[1]);
	  return -1;
	}
      fprintf(stderr, "  input file:%s\n", argv[1]);
      fprintf(stderr, "Step main 1\n");
      ReadPolytope(DATAINP, TheEXT);
      fprintf(stderr, "Step main 2\n");
      GetWeightMatrix(WMat, TheEXT);
      fprintf(stderr, "Step main 3\n");
      GetStabilizerWeightMatrix(WMat, TheGRP);
      fprintf(stderr, "Step main 4\n");
      free(WMat.TheMat);
      for (iProc=1; iProc<NbProc; iProc++)
	{
	  fprintf(stderr, "Before the Matrix sending\n");
	  MPI_SEND_MyMatrix(TheEXT, iProc, tagEXT, MPI_COMM_WORLD);
	  fprintf(stderr, "After the Matrix sending\n");
	  MPI_SEND_Group(TheGRP, iProc, tagGRP, MPI_COMM_WORLD);
	  fprintf(stderr, "After the Group sending\n");
	}
    }
  else
    {
      fprintf(stderr, "Waiting for receiving the matrix\n");
      MPI_RECV_MyMatrix(TheEXT, 0, tagEXT, MPI_COMM_WORLD);
      fprintf(stderr, "We got it\n");
      MPI_RECV_Group(TheGRP, 0, tagGRP, MPI_COMM_WORLD);
      fprintf(stderr, "We got the group\n");
    }
  fprintf(stderr, "Step main 5\n");
  fprintf(stderr, "Completion of the program\n");
  MATRIX_Free(TheEXT);
  free(TheEXT);
  ret=MPI_Finalize();
  iret = fclose(DATAINP);
}
