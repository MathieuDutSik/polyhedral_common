#include "MPQ_CDD_LRS.h"
#include "PolytopeFct.h"
#include "PolytopeEquiStab.h"
#include "Skeletton.h"

int main(int argc, char *argv[])
{
  TheGroupFormat TheGRP;
  FILE *DATAEXT=NULL;
  FILE *DATAGRP=NULL;
  FILE *DATAOUT=NULL;
  MyMatrix *TheEXT;
  int LevSearch;
  int iret;
  vector<vector<int> > ListOrb;
  Initialize_All();
  if ((TheEXT = (MyMatrix*)malloc(sizeof(MyMatrix))) == 0)
    exit (EXIT_FAILURE);
  if (argc != 5)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneGroupPolytope [DATAEXT] [DATAGRP] [LevSearch] [DATAOUT]\n");
      fprintf(stderr, "DATAEXT: The input data of the polytope\n");
      fprintf(stderr, "DATAGRP: The input data of the group\n");
      fprintf(stderr, "LevSearch: The level of the search\n");
      fprintf(stderr, "DATAOUT: The output of the enumeration\n");
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

  fprintf(stderr, "Step main 1\n");
  ReadPolytope(DATAEXT, TheEXT);
  fprintf(stderr, "Step main 2\n");
  iret=fclose(DATAEXT);

  DATAGRP=fopen(argv[2], "r");
  if (DATAGRP == NULL)
    {
      fprintf(stderr,"input: The file %s was not found\n",argv[2]);
      return -1;
    }
  fprintf(stderr, "  input file:%s\n", argv[2]);
  ReadGroup(TheGRP, DATAGRP);
  fprintf(stderr, "We have the group\n");
  iret=fclose(DATAGRP);

  LevSearch=atoi(argv[3]);
  fprintf(stderr, "LevSearch=%d\n", LevSearch);

  /* PrintGroup(stderr, TheGRP);*/
  fprintf(stderr, "Step main 5\n");
  ListOrb=DoMemoryEfficientEnum(TheGRP, TheEXT, LevSearch);

  fprintf(stderr, "Completion of the program\n");
  MATRIX_Free(TheEXT);
  free(TheEXT);

  DATAOUT=fopen(argv[4], "w");
  PrintListOrb_GAP(DATAOUT, ListOrb);
  iret=fclose(DATAOUT);
  Free_All();
}
