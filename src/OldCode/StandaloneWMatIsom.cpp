#include "MPQ_CDD_LRS.h"
#include "PolytopeFct.h"
#include "PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  FILE *FILE_WMAT1=NULL;
  FILE *FILE_WMAT2=NULL;
  vector<int> TheEquiv;
  int iret;
  WeightMatrix WMat1, WMat2;
  int TheReply;
  if (argc != 3)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneWMatIsom [WMat1] [WMat2]\n");
      fprintf(stderr, "n\n");
      return -1;
    }
  fprintf(stderr, "Reading input\n");
  FILE_WMAT1=fopen(argv[1], "r");
  if (FILE_WMAT1 == NULL)
    {
      fprintf(stderr,"input: The file %s was not found\n",argv[1]);
      return -1;
    }
  ReadWeightedMatrix(FILE_WMAT1, WMat1);
  iret=fclose(FILE_WMAT1);

  FILE_WMAT2=fopen(argv[2], "r");
  if (FILE_WMAT2 == NULL)
    {
      fprintf(stderr,"input: The file %s was not found\n",argv[1]);
      return -1;
    }
  ReadWeightedMatrix(FILE_WMAT2, WMat2);
  iret=fclose(FILE_WMAT2);

  TestEquivalenceWeightMatrix(WMat1, WMat2, &TheReply, TheEquiv);
  fprintf(stderr, "TheReply=%d\n", TheReply);
  free(WMat1.TheMat);
  free(WMat2.TheMat);
}
