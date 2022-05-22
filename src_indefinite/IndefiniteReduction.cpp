#include "Indefinite_LLL.h"
#include "NumberTheory.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3 && argc != 2) {
      std::cerr << "IndefiniteReduction [FileI] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "IndefiniteReduction [FileI]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    using Tint = mpz_class;

    std::string FileI = argv[1];
    //
    MyMatrix<T> M = ReadMatrixFile<T>(argv[1]);
    std::cerr << "We have M\n";
    //
    auto print_result = [&](std::ostream &os) -> void {
      ResultReductionIndefinite<T, Tint> ResRed =
          ComputeReductionIndefinitePermSign<T, Tint>(M);
      MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(ResRed.B);
      MyMatrix<T> M_Control = B_T * M * B_T.transpose();
      if (T_abs(DeterminantMat(B_T)) != 1) {
        std::cerr << "B_T should have determinant 1 or -1\n";
        throw TerminalException{1};
      }
      if (M_Control != ResRed.Mred) {
        std::cerr << "M_Control is not what it should be\n";
        throw TerminalException{1};
      }
      std::cerr << "B=\n";
      WriteMatrix(std::cerr, ResRed.B);
      std::cerr << "Mred=\n";
      WriteMatrix(std::cerr, ResRed.Mred);
      os << "return rec(B:=";
      WriteMatrixGAP(os, ResRed.B);
      os << ", Mred:=";
      WriteMatrixGAP(os, ResRed.Mred);
      os << ");\n";
    };
    if (argc == 2) {
      print_result(std::cerr);
    } else {
      std::string FileO = argv[2];
      std::ofstream os(FileO);
      print_result(os);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
