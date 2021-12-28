#include "NumberTheory.h"
#include "coxeter_dynkin.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "COXDYN_ComputeSymbol [FileG] [FileRoot]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    std::string FileG = argv[1];
    std::string FileRoot = argv[2];
    //
    MyMatrix<T> G = ReadMatrixFile<T>(FileG);
    MyMatrix<T> MatRoot = ReadMatrixFile<T>(FileRoot);
    std::vector<MyVector<T>> l_root;
    for (int i=0; i<MatRoot.rows(); i++) {
      MyVector<T> eLine = GetMatrixRow(MatRoot, i);
      l_root.push_back(eLine);
    }
    MyMatrix<T> CoxMat = ComputeCoxeterMatrix(G, l_root).first;
    std::cerr << "Symbol of CoxMat=" << coxdyn_matrix_to_string(CoxMat) << "\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
