#include "Temp_ThreadDualDescription.h"
int main(int argc, char *argv[])
{
  DataBank<PolyhedralEntry<mpq_class>> TheBank;
  std::vector<std::vector<int> > TheReturn;
  int NbThr;
  int TheId;
  if (argc != 6) {
    fprintf(stderr, "Number of argument is = %d\n", argc);
    fprintf(stderr, "This program is used as\n");
    fprintf(stderr, "Temp_ThreadedADM [NbThr] [DATAEXT] [DATAGRP] [PrefixSave] [DATAOUT]\n");
    fprintf(stderr, "NbThr: number of threads\n");
    fprintf(stderr, "DATAEXT: The input data of the polytope vertices\n");
    fprintf(stderr, "DATAGRP: The group for which we want orbit splitting\n");
    fprintf(stderr, "PrefixSave: The path where all data are going to be saved\n");
    fprintf(stderr, "DATAOUT: The output of computed data\n");
    fprintf(stderr, "out:  \n");
    fprintf(stderr, "n\n");
    return -1;
  }
  sscanf(argv[1], "%d", &NbThr);
  //
  fprintf(stderr, "Reading input\n");
  std::ifstream EXTfs;
  EXTfs.open(argv[2]);
  MyMatrix<mpq_class> EXT=ReadMatrix<mpq_class>(EXTfs);
  MyMatrix<mpq_class> EXTred=ColumnReduction<mpq_class>(EXT);
  EXTfs.close();
  //
  fprintf(stderr, "Reading group\n");
  std::ifstream GRPfs;
  GRPfs.open(argv[3]);
  TheGroupFormat GRP=ReadGroup(GRPfs);
  GRPfs.close();
  //
  std::string ePrefix=argv[4];
  //
  fprintf(stderr, "Creating MProc\n");
  MainProcessor MProc(NbThr);
  //
  PolyHeuristic AllArr=AllStandardHeuristic();
  TheId=MProc.MPU_GetId();
  int TheLevel=0;
  std::vector<Face> TheOutput=DUALDESC_THR_AdjacencyDecomposition(MProc, 
           TheId,
           TheBank, EXTred, 
	   GRP, 
	   AllArr, ePrefix, TheLevel);
  //
  std::ofstream OUTfs;
  OUTfs.open(argv[5]);
  VectVectInt_Magma_Print(OUTfs, TheOutput);
  OUTfs.close();
  //
  fprintf(stderr, "Completion of the program\n");
}
