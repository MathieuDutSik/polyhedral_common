#include "PolytopeEquiStab.h"

int main(int argc, char *argv[])
{
  unsigned int nof_vertices;
  TheGroupFormat TheGRP;
  FILE *DATAINP=NULL;
  Graph *g;
  int iret;
  if (argc != 2)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "StandaloneGroupGraph [DATAINP]\n");
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
  g=ReadGraphFromFile(DATAINP, nof_vertices);
  fprintf(stderr, "Step main 2\n");
  GetGraphAutomorphismGroup(g, nof_vertices, TheGRP);
  fprintf(stderr, "Step main 3\n");
  PrintGroup(stderr, TheGRP);
  fprintf(stderr, "Step main 4\n");
  fprintf(stderr, "Completion of the program\n");
  iret=fclose(DATAINP);
}
