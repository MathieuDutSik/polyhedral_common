#include "StrictPositivity.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_TestCompletePositivity [eMat]\n";
      std::cerr << "\n";
      std::cerr << "eMat: the symmetric matrix which we want to test\n";
      std::cerr << "If completely positive, we return an expression of it using integral vector\n";
      std::cerr << "If not completely positive, we return a copositive matrix having non-negative scalar product with it\n";
      return -1;
    }
    using T=mpq_class;
    using Tint=mpz_class;
    //
    std::cerr << "Reading input\n";
    //
    std::ifstream SYMMfs(argv[1]);
    MyMatrix<T> eSymmMat=ReadMatrix<T>(SYMMfs);
    //
    TestStrictPositivity<T,Tint> StrictPos = TestingAttemptStrictPositivity<T,Tint>(eSymmMat);
    WriteStrictPositivityResult(std::cerr, StrictPos);
    std::cerr << "Normal completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
