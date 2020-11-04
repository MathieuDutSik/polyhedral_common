#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "CTYP_MakeInitialFile [FileIn]\n";
      throw TerminalException{1};
    }
    using Tmat=mpz_class;
    //
    std::string FileIn  = argv[1];
    //
    std::ifstream is(FileIn);
    int nbType;
    is >> nbType;
    auto get_random_equivalent=[&](MyMatrix<Tmat> const& eMat) -> MyMatrix<Tmat> {
      int nbRow = eMat.rows();
      int n = eMat.cols();
      std::vector<int> ePerm = RandomPermutation(nbRow);
      MyMatrix<Tmat> eUnimod = RandomUnimodularMatrix<Tmat>(n);
      MyMatrix<Tmat> eProd = eUnimod * eMat;
      MyMatrix<Tmat> RetMat(nbRow, n);
      for (int iRow=0; iRow<nbRow; iRow++) {
        int jRow = ePerm[iRow];
        for (int i=0; i<n; i++)
          RetMat(iRow, i) = eProd(jRow, i);
      }
      return RetMat;
    };
    auto TestEqual=[&](MyMatrix<Tmat> const& M1, MyMatrix<Tmat> const& M1) -> bool {
      int nbRow = eMat.rows();
      int n = eMat.cols();
      for (int iRow=0; iRow<nbRow; iRow++)
        for (int i=0; i<n; i++)
          if (M1(iRow, i) != M2(iRow, i))
            return false;
      return true;
    };
    for (int iType=0; iType<nbType; iType++) {
      std::cerr << "iType : " << iType << " / " << nbType << "\n";
      MyMatrix<Tmat> eMat1 = ReadMatrix<Tmat>(is);
      MyMatrix<Tmat> eMat1_Can = LinPolytopeIntegral_CanonicForm<Tmat>(eMat1);
      for (int i=0; i<10; i++) {
        MyMatrix<Tmat> eMat2 = get_random_equivalent(eMat1);
        MyMatrix<Tmat> eMat2_Can = LinPolytopeIntegral_CanonicForm<Tmat>(eMat1);
        if (!TestEqual(eMat1_Can, eMat2_Can)) {
          std::cerr << "Inconsistencw in the canonical code\n";
          throw TerminalException{1};
        }
      }
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
