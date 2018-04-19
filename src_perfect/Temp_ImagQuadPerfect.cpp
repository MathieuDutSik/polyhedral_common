#include "NumberTheory.h"
#include "Temp_PerfectForm.h"

int main(int argc, char *argv[])
{
  try {
    std::vector<std::vector<int> > TheReturn;
    int NbThr;
    DataBank<PolyhedralEntry<mpq_class>, PolyhedralEquiv> TheBank;
    ListPerfectForm<mpq_class> ListPerf;
    
    int TheId;
    int n, eSum, eProd;
    mpq_class eSum_T, eProd_T;
    
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_ImagQuadPerfect [NbThr] [n] [eSum] [eProd] out\n";
      std::cerr << "n\n";
      return -1;
    }
    sscanf(argv[1], "%d", &NbThr);
    sscanf(argv[2], "%d", &n);
    sscanf(argv[3], "%d", &eSum);
    sscanf(argv[4], "%d", &eProd);
    eSum_T=eSum;
    eProd_T=eProd;
    MainProcessor MProc(NbThr);
    LinSpaceMatrix<mpq_class> LinSpa=ComputeImagQuadraticSpace(n, eSum_T, eProd_T);
    PolyHeuristic AllArr=AllStandardHeuristic();
    TheId=MProc.MPU_GetId();
    VoronoiAlgo_THR_EnumeratePerfectForm<mpq_class>(MProc, TheId, 
						    TheBank,
						    LinSpa, 
						    AllArr, ListPerf);
    std::ofstream OUTfs(argv[5]);
    VoronoiAlgo_PrintListMatrix(OUTfs, LinSpa, ListPerf);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
