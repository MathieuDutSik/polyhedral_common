#include "Indefinite_LLL.h"
#include "NumberTheory.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "TestPermutationSignCanonic [FileI]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;

    std::string FileI = argv[1];
    //
    MyMatrix<T> M = ReadMatrixFile<T>(argv[1]);
    if (!IsSymmetricMatrix(M)) {
      std::cerr << "The matrix M should be symmetric\n";
      throw TerminalException{1};
    }
    int nbRow = M.rows();
    //
    MyMatrix<T> Mcan = CanonicalizationPermutationSigns(M).second;
    for (int iter=0; iter<50; iter++) {
      std::cerr << "iter=" << iter << "\n";
      std::vector<int> ePerm = RandomPermutation<int>(nbRow);
      MyMatrix<T> eP = ZeroMatrix<T>(nbRow,nbRow);
      for (int iRow=0; iRow<nbRow; iRow++) {
        eP(iRow,ePerm[iRow]) = 2 * (random() % 2) - 1;
      }
      MyMatrix<T> Mnew = eP * M * eP.transpose();
      MyMatrix<T> Mnew_can = CanonicalizationPermutationSigns(Mnew).second;
      if (Mnew != Mnew_can) {
        std::cerr << "The matrices are not equal\n";
        throw TerminalException{1};
      }
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
