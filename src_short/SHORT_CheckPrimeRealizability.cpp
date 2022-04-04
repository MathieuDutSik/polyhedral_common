//#include "mpreal_related.h"
#include "NumberTheory.h"
#include "SHORT_ShortestConfig.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "SHORT_CheckPrimeRealizability [FilePrime] [FileIn] [FileOut]\n";
      std::cerr << "\n";
      std::cerr
          << "[FilePrime] : The input file of the list of primes considered\n";
      std::cerr << "[FileIn]    : The input file of the system\n";
      std::cerr
          << "[FileOut]   : The output file of the program (GAP readable)\n";
      return -1;
    }
    //
    std::string FilePrime = argv[1];
    std::string FileIn = argv[2];
    std::string FileOut = argv[3];
    //
    std::vector<PrimeListAllowed> ListP;
    std::cerr << "FilePrime=" << FilePrime << "\n";
    std::ifstream is(FilePrime);
    int n;
    is >> n;
    std::cerr << "n=" << n << "\n";
    int nbPrime;
    is >> nbPrime;
    std::cerr << "nbPrime=" << nbPrime << "\n";
    for (int iPrime = 0; iPrime < nbPrime; iPrime++) {
      int ePrime;
      is >> ePrime;
      std::cerr << "iPrime=" << iPrime << " ePrime=" << ePrime << "\n";
      bool DoWeListFeasible;
      is >> DoWeListFeasible;
      std::cerr << "DoWeListFeasible=" << DoWeListFeasible << "\n";
      int nbCase;
      is >> nbCase;
      std::cerr << "nbCase=" << nbCase << "\n";
      std::vector<MyVector<int>> ListCases;
      for (int iCase = 0; iCase < nbCase; iCase++) {
        MyVector<int> eCase(n);
        std::cerr << "iCase=" << iCase << " / " << nbCase << "\n";
        for (int i = 0; i < n; i++) {
          int eVal;
          is >> eVal;
          eCase(i) = eVal;
        }
        ListCases.push_back(eCase);
      }
      PrimeListAllowed eCaseP{ePrime, DoWeListFeasible, ListCases};
      ListP.push_back(eCaseP);
    }
    std::cerr << "We have ListP\n";
    //
    std::ifstream isI(FileIn);
    std::ofstream osO(FileOut);
    osO << "return [";
    bool IsFirst = true;
    int nbCase = 0;
    int nbCaseCorr = 0;
    while (1) {
      int eStatus;
      isI >> eStatus;
      std::cerr << "nbCase=" << nbCase << " nbCaseCorr=" << nbCaseCorr << "\n";
      if (eStatus == 0)
        break;
      MyMatrix<mpq_class> M = ReadMatrix<mpq_class>(isI);
      //      std::cerr << "Before IsMathingListOfPrime\n";
      bool test = IsMatchingListOfPrimes(ListP, M);
      //      std::cerr << "After IsMathingListOfPrime\n";
      if (test) {
        if (IsFirst == false)
          osO << ",\n";
        IsFirst = false;
        WriteMatrixGAP(osO, M);
        nbCaseCorr++;
      }
      nbCase++;
    }
    osO << "];\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
