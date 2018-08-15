#include "Copositivity.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveKernelMin [DATASYMM] [DATABASIS] [MaxNorm]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The input data of the symmetric matrix\n";
      std::cerr << "DATABASIS: The input data of the matrix basis\n";
      std::cerr << "It returns the copositive minimum vectors assuming that the basis give non-negative pairwise scalar products\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    std::string eFile1=argv[1];
    if (!IsExistingFile(eFile1)) {
      std::cerr << "The file eFile1=" << eFile1 << " is missing\n";
      throw TerminalException{1};
    }
    using T=mpq_class;
    using Tint=mpz_class;
    //
    std::ifstream SYMMfs(eFile1);
    MyMatrix<T> eSymmMat=ReadMatrix<T>(SYMMfs);
    //
    std::string eFile2=argv[2];
    if (!IsExistingFile(eFile2)) {
      std::cerr << "The file eFile2=" << eFile2 << " is missing\n";
      throw TerminalException{1};
    }
    std::ifstream BASfs(eFile2);
    MyMatrix<Tint> TheBasis=ReadMatrix<Tint>(BASfs);
    //
    int MaxNorm_i;
    sscanf(argv[3], "%d", &MaxNorm_i);
    T MaxNorm=MaxNorm_i;
    //
    std::vector<MyVector<Tint>> LVect = EnumerateShortVectorInCone_UnderPositivityCond<T,Tint>(eSymmMat, TheBasis, MaxNorm);
    std::cerr << "|LVect|=" << LVect.size() << "\n";
    for (auto & eVect : LVect) {
      std::cerr << "eVect=";
      WriteVector(std::cerr, eVect);
    }
    //
    std::cerr << "Completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
