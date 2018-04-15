#include "IterationAlgorithm.h"

int main(int argc, char *argv[])
{
  FILE *DATAINP=NULL;
  MatDOUBL OrigMat, NewMat;
  int eRet;
  int n;
  int i;
  int iret;
  VectDOUBL eVectReal;
  VectDOUBL ListEigVal;
  MatDOUBL ListEigVect;
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneEigenvalue [DATAINP]\n");
      fprintf(stderr, "DATAINP: The input data of the form\n");
      fprintf(stderr, "n\n");
      fprintf(stderr, "a11 .... an1\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, "a1n .... ann\n");
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

  eRet=fscanf(DATAINP, "%d", &n);
  fprintf(stderr, "n=%d\n", n);
  OrigMat=ReadMat_Doubl(DATAINP, n, n);
  jacobi_double(OrigMat, ListEigVal, ListEigVect);
  for (i=0; i<n; i++)
    {
      fprintf(stderr, "i=%d eEig=%lg\n", i, ListEigVal[i]);
    }
  fprintf(stderr, "Completion of the program\n");
  iret = fclose(DATAINP);
}
