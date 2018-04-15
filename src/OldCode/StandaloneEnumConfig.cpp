#include "MatrixOperation.h"
#include "EnumerationConfig.h"

int main(int argc, char *argv[])
{
  FILE *DATAINP=NULL;
  FILE *DATAOUT=NULL;
  MatDOUBL TheMat;
  int k;
  double MaxDet;
  vector<MatINT> ListSubspace;
  int eRet;
  int iret;
  int n, i, j;
  int nbSubspace, iSub, iK;
  double eVal;
  vector<double> ListDet;
  VectDOUBL eVectReal;
  if (argc != 3)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneEnum [DATAINP] [DATAOUT]\n");
      fprintf(stderr, "DATAINP: The input data of the form\n");
      fprintf(stderr, "k\n");
      fprintf(stderr, "MaxDet\n");
      fprintf(stderr, "n\n");
      fprintf(stderr, "a11 .... an1\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, ".          .\n");
      fprintf(stderr, "a1n .... ann\n");
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
  eRet=fscanf(DATAINP, "%lg", &MaxDet);
  fprintf(stderr, "MaxDet=%lg\n", MaxDet);
  eRet=fscanf(DATAINP, "%d", &n);
  fprintf(stderr, "n=%d\n", n);
  for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)
	{
	  eRet=fscanf(DATAINP, "%lg", &eVal);
	  eVectReal.push_back(eVal);
	}
      TheMat.push_back(eVectReal);
      eVectReal.clear();
    }
  EnumerationKspaces(TheMat, k, MaxDet, ListSubspace, ListDet);
  nbSubspace=ListSubspace.size();
  fprintf(stderr, "We found %d subspaces\n", nbSubspace);
  fprintf(DATAOUT, "We found %d subspaces\n", nbSubspace);
  for (iSub=0; iSub<nbSubspace; iSub++)
    {
      fprintf(DATAOUT, "%d/%d det=%lg\n", iSub, nbSubspace, ListDet[iSub]);
      for (iK=0; iK<k; iK++)
	{
	  for (i=0; i<n; i++)
	    fprintf(DATAOUT, " %d", ListSubspace[iSub][iK][i]);
	  fprintf(DATAOUT, "\n");
	}
    }
  fprintf(stderr, "Completion of the program\n");
  iret=fclose(DATAINP);
  iret=fclose(DATAOUT);
}
