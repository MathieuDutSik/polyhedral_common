#include "MPQ_CDD_LRS.h"
#include "PolytopeFct.h"
#include "PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  TheGroupFormat TheGRP;
  FILE *DATAINP=NULL;
  int iret;
  MyMatrix TheEXT, TheMatInv;
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneGroupPolytope [DATAINP]\n");
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
  MATRIX_ReadFernando(DATAINP, &TheEXT);
    /*  ReadPolytope(DATAINP, &TheEXT);*/
  MATRIX_Print(stderr, &TheEXT);
  MATRIX_Allocate(&TheMatInv, TheEXT.nbRow, TheEXT.nbCol);
  MATRIX_Inverse_destroy(&TheEXT, &TheMatInv);
  MATRIX_Print(stderr, &TheMatInv);
  iret=fclose(DATAINP);
}
