#include "GampMatlab.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "ComputeChromaticLowerBounds [inputAdjMat]\n");
      return -1;
    }
    //
    std::string eFile = argv[1];
    auto GetListDegree=[&]() -> std::vector<int> {
      std::ifstream is(eFile);
      int nbVert, eAdj;
      int nnz=0;
      is >> nbVert;
      std::vector<int> ListDeg(nbVert);
      for (int iVert=0; iVert<nbVert; iVert++) {
        is >> nbAdj;
        for (int iAdj=0; iAdj<nbAdj; iAdj++)
          is >> eAdj;
        ListDeg[iVert] = nbAdj;
      }
      return ListDeg;
    };
    //
    std::vector<int> ListDeg = GetListDegree();
    int nnz=0;
    for (auto & eVal : ListDeg)
      nnz += eVal;
    //
    using T2=Eigen::Triplet<double>;
    std::vector<T2> tripletList_A(nnz);
    std::vector<T2> tripletList_D(nbVert);
    std::ifstream is(eFile);
    int nbVert, eAdj;
    is >> nbVert;
    int iNNZ=0;
    for (int iVert=0; iVert<nbVert; iVert++) {
      int nbAdj;
      is >> nbAdj;
      for (int iAdj=0; iAdj<nbAdj; iAdj++) {
        int eAdj;
        is >> eAdj;
        tripletList_A[iNNZ]=T2(iVert,eAdj,double(1));
        iNNZ++;
      }
      tripleList_D[iVert] = T2(iVert, iVert, nbAdj);
    }
    MySparseMatrix<double> A(nbVert, nbVert);
    A.setFromTriplets(tripletList_A.begin(), tripletList_A.end());
    MySparseMatrix<double> D(nbVert, nbVert);
    D.setFromTriplets(tripletList_D.begin(), tripletList_D.end());
    //
    auto GetEigenvalues=[&](MySparseMatrix<double> const& SpMat) -> MyVector<double> {
      Eigen::SelfAdjointEigenSolver<MySparseMatrix<double>> eig(SpMat);
      MyVector<double> ListEig=eig.eigenvalues();
      return ListEig;
    };
    //
    MyVector<double> mu = GetEigenvalues(A);
    MyVector<double> delta = GetEigenvalues(A + D);
    MyVector<double> theta = GetEigenvalues(D - A);

    std::cerr << "Normal termination of ComputeChromaticLowerBounds\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
