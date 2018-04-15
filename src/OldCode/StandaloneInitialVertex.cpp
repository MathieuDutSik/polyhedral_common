#include "MPQ_CDD_LRS.h"
#include "PolytopeFct.h"
#include "PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  FILE *DATAINP=NULL;
  MyMatrix *TheEXT;
  vector<int> eInc;
  int iret;
  Initialize_All();
  if ((TheEXT = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneInitialVertex [DATAINP]\n");
      fprintf(stderr, "DATAINP: The input data of the form\n");
      fprintf(stderr, "n\n");
      return -1;
    }
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
  FindOneInitialVertex(TheEXT, eInc);

  MATRIX_Free(TheEXT);
  free(TheEXT);
  fprintf(stderr, "Completion of the program\n");
  Free_All();
  iret = fclose(DATAINP);
}
