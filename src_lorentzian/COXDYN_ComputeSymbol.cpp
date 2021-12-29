#include "NumberTheory.h"
#include "coxeter_dynkin.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 3 && argc != 4) {
      std::cerr << "COXDYN_ComputeSymbol [FileG] [FileRoot]\n";
      std::cerr << "or\n";
      std::cerr << "COXDYN_ComputeSymbol [FileG] [FileRoot] [OutFile]\n";
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
    //    std::cerr << "CoxMat=\n";
    //    WriteMatrix(std::cerr, CoxMat);
    std::string symb = coxdyn_matrix_to_string(CoxMat);
    auto prt=[&](std::ostream & os) -> void {
      os << "return rec(CoxMat:=";
      WriteMatrixGAP(os, CoxMat);
      os << ", symb:=\"" << symb << "\");\n";
    };
    if (argc == 3) {
      prt(std::cerr);
    } else {
      std::string out(argv[3]);
      std::ofstream os(out);
      prt(os);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
