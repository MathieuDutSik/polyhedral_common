#include "IterationAlgorithm.h"

int main(int argc, char *argv[])
{
  FILE *DATAINP=NULL;
  FILE *DATAOUT=NULL;
  MatDOUBL OrigMat, DirMat, NewMat;
  int k;
  int eRet;
  int n;
  VectDOUBL eVectReal;
  int iret;
  if (argc != 3)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneFlipping [DATAINP] [DATAOUT]\n");
      fprintf(stderr, "DATAINP: The input data of the form\n");
      fprintf(stderr, "k\n");
      fprintf(stderr, "n\n");
      fprintf(stderr, "a11 .... an1\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, "a1n .... ann\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "b11 .... bn1\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, "b1n .... bnn\n");
      fprintf(stderr, "DATAOUT: output of the generated objects\n");
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

  fprintf(stderr, "Reading input\n");
  DATAOUT=fopen(argv[2], "w");
  if (DATAOUT == NULL)
    {
      fprintf(stderr,"input: The file %s was not found\n",argv[1]);
      return -1;
    }
  fprintf(stderr, "  input file:%s\n", argv[1]);

  eRet=fscanf(DATAINP, "%d", &k);
  fprintf(stderr, "k=%d\n", k);
  eRet=fscanf(DATAINP, "%d", &n);
  fprintf(stderr, "n=%d\n", n);
  OrigMat=ReadMat_Doubl(DATAINP, n, n);
  DirMat=ReadMat_Doubl(DATAINP, n, n);

  NewMat=KFlippingOperation(OrigMat, k, DirMat);
  fprintf(stderr, "NewMat=\n");
  PrintMat_DoublNice(stdout, NewMat);
  PrintMat_Doubl(DATAOUT, NewMat);
  fprintf(stderr, "Completion of the program\n");
  iret = fclose(DATAINP);
  iret = fclose(DATAOUT);
}
