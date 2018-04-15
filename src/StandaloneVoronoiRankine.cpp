#include "VoronoiRankine.h"

int main(int argc, char *argv[])
{
  try {
    int k, n;
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "StandaloneVoronoiRankine n k [TOLFILE] [OUTFILE] [OUTFILE_GAP]\n";
      return -1;
    }
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &k);
    //
    std::ifstream DATATOL(argv[3]);
    AllTol<double> RecTol=ReadToleranceFile<double>(DATATOL);
    //
    std::vector<MyMatrix<double> > ListPerfect=VoronoiEnumerationRankin<double>(n, k, RecTol);
    //
    std::ofstream DATAOUT(argv[4]);
    PrintFamilyKPerfectMatrices(DATAOUT, n, k, 
				ListPerfect, RecTol);
    //
    std::ofstream DATAOUT_GAP(argv[5]);
    PrintFamilyKPerfectMatricesGAP(DATAOUT_GAP, n, k, 
				   ListPerfect, RecTol);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
