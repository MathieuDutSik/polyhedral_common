#include "NumberTheory.h"
#include "Copositivity.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CP_CopositiveListCone [DATASYMM] [MaxNorm] [OUTFILE]\n";
      std::cerr << "\n";
      std::cerr << "DATASYMM: The input data of the copositive symmetric matrix\n";
      std::cerr << "It returns the list of cones that is used to test that the matrix is copositive\n";
      return -1;
    }
    using T=mpq_class;
    using Tint=mpz_class;
    //
    std::cerr << "Reading input\n";
    //
    MyMatrix<T> eSymmMat=ReadMatrixFile<T>(argv[1]);
    //
    int MaxNorm_i;
    sscanf(argv[2], "%d", &MaxNorm_i);
    T MaxNorm=MaxNorm_i;
    //
    RequestCopositivity<T> CopoReq{MaxNorm, true};
    int n=eSymmMat.rows();
    std::cerr << "eSymmMat.rows=" << eSymmMat.rows() << " eSymmMat.cols=" << eSymmMat.cols() << "\n";
    MyMatrix<Tint> TheBasis=IdentityMat<Tint>(n);
    CopositivityEnumResult<Tint> CopoRes=KernelEnumerateShortVectorInCone<T,Tint>(eSymmMat, TheBasis, CopoReq);
    std::cerr << "nbCone=" << CopoRes.nbCone << "\n";
    if (CopoRes.test == false) {
      std::cerr << "Matrix is not Copositive\n";
    }
    else {
      std::ofstream OUTfs(argv[3]);
      for (int iCone=0; iCone<CopoRes.nbCone; iCone++) {
	OUTfs << "iCone=" << iCone << "/" << CopoRes.nbCone << " Basis\n";
	MyMatrix<Tint> eMat=CopoRes.ListBasis[iCone];
	//      std::cerr << "eMat.rows=" << eMat.rows() << " eMat.cols=" << eMat.cols() << "\n";
	WriteMatrix(OUTfs, eMat);
	OUTfs << "Pairwise scalar products:\n";
	for (int i=0; i<n; i++) {
	  //	std::cerr << "i=" << i << "\n";
	  MyVector<Tint> eVect1=eMat.row(i);
	  //	std::cerr << "eVect1.rows=" << eVect1.rows() << " eVect1.cols=" << eVect1.cols() << "\n";
	  for (int j=0; j<=i; j++) {
	    //	  std::cerr << "  j=" << j << "\n";
	    MyVector<Tint> eVect2=eMat.row(j);
	    //	  std::cerr << "eVect2.rows=" << eVect2.rows() << " eVect2.cols=" << eVect2.cols() << "\n";
	    //	  std::cerr << "  Step 1\n";
	    T eScal=ScalarProductQuadForm<T,Tint>(eSymmMat, eVect1, eVect2);
	    //	  std::cerr << "  Step 2\n";
	    OUTfs << " " << eScal;
	  }
	  OUTfs << "\n";
	}
      }
    }
    std::cerr << "Normal completion of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
