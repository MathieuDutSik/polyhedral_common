// clang-format off
#include "NumberTheory.h"
#include "Copositivity.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 3 || argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveMin [DATASYMM] [MaxNorm]\n";
      std::cerr << "or\n";
      std::cerr << "CP_CopositiveMin [DATASYMM] [MaxNorm] [InitialBasis]\n";
      std::cerr << "\n";
      std::cerr
          << "DATASYMM: The input data of the copositive symmetric matrix A\n";
      std::cerr << "It returns the list of integer vectors v in Z^n_{>= 0} "
                   "such that A[v] <= MaxNorm\n";
      std::cerr << "\n";
      std::cerr << "If InitialBasis is not put in argument, then it is the "
                   "standard basis {e1, ...., en}\n";
      return -1;
    }
    using T = mpq_class;
    using Tint = mpz_class;
    //
    std::cerr << "Reading input\n";
    //
    MyMatrix<T> eSymmMat = ReadMatrixFile<T>(argv[1]);
    std::cerr << "eSymmMat=\n";
    WriteMatrix(std::cerr, eSymmMat);
    //
    MyMatrix<Tint> InitialBasis = IdentityMat<Tint>(eSymmMat.rows());
    if (argc == 4)
      InitialBasis = ReadMatrixFile<Tint>(argv[3]);
    //
    int MaxNorm_i;
    sscanf(argv[2], "%d", &MaxNorm_i);
    T MaxNorm = MaxNorm_i;
    //
    RequestCopositivity<T> CopoReq{MaxNorm, false};
    CopositivityEnumResult<Tint> CopoRes =
        EnumerateCopositiveShortVector<T, Tint>(eSymmMat, InitialBasis,
                                                CopoReq);
    std::cerr << "nbCone=" << CopoRes.nbCone << "\n";
    if (CopoRes.test == false) {
      std::cerr << "Matrix is not Copositive\n";
      std::cerr << "Nature=" << CopoRes.eResult.strNature << "\n";
      std::cerr << "eVect1=";
      WriteVector(std::cerr, CopoRes.eResult.eVectResult1);
    } else {
      int n = eSymmMat.rows();
      std::cerr << "Matrix is Copositive\n";
      int nbVect = CopoRes.TotalListVect.size();
      std::cerr << "nbVect=" << nbVect << "\n";
      for (int iVect = 0; iVect < nbVect; iVect++) {
        MyVector<Tint> eVect = CopoRes.TotalListVect[iVect];
        T eNorm = EvaluationQuadForm(eSymmMat, eVect);
        std::cerr << "iVect=" << iVect << " eVect=[";
        for (int i = 0; i < n; i++) {
          if (i > 0)
            std::cerr << ",";
          std::cerr << eVect[i];
        }
        std::cerr << "] norm=" << eNorm << "\n";
      }
    }
    //
    std::cerr << "Normal completion of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
