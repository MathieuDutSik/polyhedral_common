#define PRINT_NBCONE
#include "NumberTheory.h"
#include "Copositivity.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2 || argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_ComputeCopositiveMin [DATASYMM]\n";
      std::cerr << "or\n";
      std::cerr << "CP_ComputeCopositiveMin [DATASYMM] [InitialBasois]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The input data of the symmetric matrix\n";
      std::cerr << "It returns true if the matrix is copositive. If not in returns a vector V with A[V] <0 and V with all coordinates non-negative\n";
      std::cerr << "\n";
      std::cerr << "If InitialBasis is not put in argument, then it is the standard basis {e1, ...., en}\n";
      return -1;
    }
    //
    std::cerr << "Reading input\n";
    //
    using T=mpq_class;
    using Tint=int;
    MyMatrix<T> eSymmMat=ReadMatrixFile<T>(argv[1]);
    std::cerr << "eSymmMat=\n";
    WriteMatrix(std::cerr, eSymmMat);
    //
    MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
    if (argc == 3)
      InitialBasis = ReadMatrixFile<Tint>(argv[2]);
    //
    Tshortest<T,Tint> eSh = T_CopositiveShortestVector<T,Tint>(eSymmMat, InitialBasis);
    std::cerr << "eMin=" << eSh.eMin << "\n";
    WriteMatrix(std::cerr, eSh.SHV);
    //
    std::cerr << "Normal completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
