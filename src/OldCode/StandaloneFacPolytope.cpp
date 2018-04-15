#include "MPQ_CDD_LRS.h"
#include "PolytopeFct.h"
#include "PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  TheGroupFormat TheGRP;
  TheGroupFormat TheStab;
  TheGroupFormat TheStabRed;
  FILE *DATAINP=NULL;
  MyMatrix *TheEXT;
  MyMatrix *TheEXTred1;
  MyMatrix *TheFAC;
  MyMatrix *TheFACred1;
  WeightMatrix WMat;
  vector<vector<int> > TheReturn, TheReturnRed1;
  vector<int> eO, eFlip;
  int iret;
  Initialize_All();
  if ((TheEXT = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheEXTred1 = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheFAC = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if ((TheFACred1 = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneFacPolytope [DATAINP]\n");
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
  GetWeightMatrix(WMat, TheEXT);
  fprintf(stderr, "Step main 3\n");
  GetStabilizerWeightMatrix(WMat, TheGRP);
  fprintf(stderr, "Step main 4\n");
  PrintGroup(stderr, TheGRP);
  fprintf(stderr, "Step main 5\n");

  /*  CDD_DualDescription(TheEXT, TheFAC);*/
  LRS_DualDescription(TheEXT, TheFAC);
  fprintf(stderr, "Step main 6\n");
  WritePolytope(stderr, TheFAC);
  fprintf(stderr, "Step main 7\n");
  OrbitSplitting(TheEXT, TheFAC, TheGRP, TheReturn);
  fprintf(stderr, "Step main 8\n");
  PrintListOrbit(stderr, TheReturn);
  fprintf(stderr, "Step main 9\n");
  eO=TheReturn[0];
  fprintf(stderr, "Step main 10\n");
  MatrixReduction(TheEXT, eO, TheEXTred1);
  fprintf(stderr, "Step main 11\n");
  GetStabilizer(TheGRP, eO, TheStab);
  fprintf(stderr, "Step main 12\n");
  PrintGroup(stderr, TheStab);
  fprintf(stderr, "Step main 13\n");
  MATRIX_Print(stderr, TheEXTred1);
  fprintf(stderr, "Step main 14\n");
  CDD_DualDescription(TheEXTred1, TheFACred1);
  fprintf(stderr, "Step main 15\n");
  MATRIX_Print(stderr, TheFACred1);
  fprintf(stderr, "Step main 16\n");
  TheStabRed=ReducedGroupAction(TheStab, eO);
  fprintf(stderr, "Step main 17\n");
  PrintGroup(stderr, TheStabRed);
  fprintf(stderr, "Step main 17\n");
  OrbitSplitting(TheEXTred1, TheFACred1, TheStabRed, TheReturnRed1);
  fprintf(stderr, "Step main 18\n");
  fprintf(stderr, "eO=");
  PrintVectorInt(stderr, eO);
  eFlip=ComputeFlipping(TheEXT, eO, TheReturnRed1[0]);
  fprintf(stderr, "eFlip=");
  PrintVectorInt(stderr, eFlip);

  MATRIX_Free(TheEXT);
  free(TheEXT);
  MATRIX_Free(TheFAC);
  free(TheFAC);
  free(WMat.TheMat);
  Free_All();
  MATRIX_Free(TheEXTred1);
  free(TheEXTred1);
  MATRIX_Free(TheFACred1);
  free(TheFACred1);
  fprintf(stderr, "Completion of the program\n");
  iret=fclose(DATAINP);
}
