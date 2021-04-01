#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"

#include "NumberTheory.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "POLY_MakeInitialFileGAP [FileIn]\n";
      throw TerminalException{1};
    }
    //
    std::string FileIn  = argv[1];
    //
    using T = mpq_class;
    using Telt = permutalib::DoubleSidedPerm;
    datagap::DataGAP<T,Telt> dataEXT = datagap::ParseGAPFile<T,Telt>(FileIn);
    //
    MyMatrix<T> EXT1 = datagap::ConvertGAPread_MyMatrixT(dataEXT);
    //
    std::ifstream is(FileIn);
    auto get_random_equivalent=[&](MyMatrix<T> const& eMat) -> MyMatrix<T> {
      int nbRow = eMat.rows();
      int n = eMat.cols();
      std::vector<int> ePerm = RandomPermutation(nbRow);
      std::vector<int> AttV(nbRow,0);
      std::cerr << "ePerm=[";
      for (int iRow=0; iRow<nbRow; iRow++) {
        if (iRow > 0)
          std::cerr << ",";
        std::cerr << ePerm[iRow];
        AttV[ePerm[iRow]] += 1;
      }
      std::cerr << "]\n";
      for (int iRow=0; iRow<nbRow; iRow++) {
        if (AttV[iRow] != 1) {
          std::cerr << "Error in ePerm\n";
          throw TerminalException{1};
        }
      }
      MyMatrix<T> eUnimod = RandomUnimodularMatrix<T>(n);
      MyMatrix<T> eProd = eMat * eUnimod;
      MyMatrix<T> RetMat(nbRow, n);
      for (int iRow=0; iRow<nbRow; iRow++) {
        int jRow = ePerm[iRow];
        for (int i=0; i<n; i++)
          RetMat(iRow, i) = eProd(jRow, i);
      }
      return RetMat;
    };
    MyMatrix<T> EXT1_Can = CanonicalizationPolytope(EXT1);
    int n_iter = 10;
    for (int i_iter=0; i_iter<n_iter; i_iter++) {
      std::cerr << "i_iter=" << i_iter << " / " << n_iter << "\n";
      MyMatrix<T> EXT2 = get_random_equivalent(EXT1);
      MyMatrix<T> EXT2_Can = CanonicalizationPolytope(EXT2);
      if (!TestEqualityMatrix(EXT1_Can, EXT2_Can)) {
        std::cerr << "Inconsistency in the canonical code\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal end of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
