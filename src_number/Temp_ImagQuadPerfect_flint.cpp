#include "Temp_PerfectForm.h"
#include "fmpq_xx.h"

int main(int argc, char *argv[])
{
  try {
    std::vector<std::vector<int> > TheReturn;
    int NbThr;
    DataBank<PolyhedralEntry<fmpq_class>, PolyhedralEquiv> TheBank;
    ListPerfectForm<fmpq_class> ListPerf;
    
    int TheId;
    int n, eSum, eProd;
    fmpq_class eSum_T, eProd_T;
    
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "Temp_ImagQuadPerfect [NbThr] [n] [eSum] [eProd] out\n";
      return -1;
    }
    sscanf(argv[1], "%d", &NbThr);
    sscanf(argv[2], "%d", &n);
    sscanf(argv[3], "%d", &eSum);
    sscanf(argv[4], "%d", &eProd);
    eSum_T=eSum;
    eProd_T=eProd;
    MainProcessor MProc(NbThr);
    LinSpaceMatrix<fmpq_class> LinSpa=ComputeImagQuadraticSpace(n, eSum_T, eProd_T);
    PolyHeuristic AllArr=AllStandardHeuristic();
    TheId=MProc.MPU_GetId();
    VoronoiAlgo_THR_EnumeratePerfectForm<fmpq_class>(MProc, TheId, 
						     TheBank, LinSpa, AllArr, ListPerf);
    std::ofstream OUTfs(argv[5]);
    VoronoiAlgo_PrintListMatrix(OUTfs, LinSpa, ListPerf);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
