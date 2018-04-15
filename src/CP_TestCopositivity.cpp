#define PRINT_SPLIT_CONE

#include "Copositivity.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestCopositivity [DATASYMM]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The input data of the symmetric matrix\n";
      std::cerr << "It returns true if the matrix is copositive. If not in returns a vector V with A[V] <0 and V with all coordinates non-negative\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    std::string eFile=argv[1];
    if (!IsExistingFile(eFile)) {
      std::cerr << "File eFile=" << eFile << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream SYMMfs(eFile);
    MyMatrix<mpq_class> eSymmMat=ReadMatrix<mpq_class>(SYMMfs);
    std::cerr << "eSymmMat=\n";
    WriteMatrix(std::cerr, eSymmMat);
    //
    SingleTestResult eResult = TestCopositivity(eSymmMat);
    if (eResult.test) {
      std::cerr << "The matrix is indeed copositive\n";
    }
    else {
      std::cerr << "The matrix is not copositive\n";
      std::cerr << "Nature of violation=" << eResult.strNature << "\n";
      std::cerr << "V=";
      WriteVector(std::cerr, eResult.eVectResult1);
    }
    //
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
