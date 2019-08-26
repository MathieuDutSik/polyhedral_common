#include "MAT_Matrix.h"
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
      is >> nbVert;
      std::vector<int> ListDeg(nbVert);
      for (int iVert=0; iVert<nbVert; iVert++) {
        int nbAdj;
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
    std::ifstream is(eFile);
    int nbVert;
    is >> nbVert;
    std::vector<T2> tripletList_A(nnz);
    std::vector<T2> tripletList_D(nbVert);
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
      tripletList_D[iVert] = T2(iVert, iVert, nbAdj);
    }
    MySparseMatrix<double> A(nbVert, nbVert);
    A.setFromTriplets(tripletList_A.begin(), tripletList_A.end());
    std::cerr << "We have A matrix\n";
    MySparseMatrix<double> D(nbVert, nbVert);
    D.setFromTriplets(tripletList_D.begin(), tripletList_D.end());
    std::cerr << "We have D matrix\n";
    //
    auto GetEigenvalues=[&](MySparseMatrix<double> const& SpMat) -> MyVector<double> {
      Eigen::SelfAdjointEigenSolver<MySparseMatrix<double>> eig(SpMat);
      MyVector<double> ListEig=eig.eigenvalues();
      MyVector<double> ListEigRev(nbVert);
      for (int iVert=0; iVert<nbVert; iVert++)
        ListEigRev(iVert) = ListEig(nbVert - 1 - iVert);
      return ListEigRev;
    };
    //
    MyVector<double> mu = GetEigenvalues(A);
    std::cerr << "We have mu eigenvalues\n";
    MyVector<double> delta = GetEigenvalues(A + D);
    std::cerr << "We have delta eigenvalues\n";
    MyVector<double> theta = GetEigenvalues(D - A);
    std::cerr << "We have theta eigenvalues\n";
    int n = nbVert;
    //
    double BestLowerBound = 0;
    std::string RealizingBound = "unset";
    auto UpdateBounds=[&](double const& eBound, std::string const& RealBound) -> void {
      if (eBound > BestLowerBound) {
        BestLowerBound = eBound;
        RealizingBound = RealBound;
      }
    };
    //
    double LowerBoundHoffman = 1 + mu(0) / (- mu(n-1));
    UpdateBounds(LowerBoundHoffman, "Hoffman bound");
    std::cerr << "Hoffman lower bound = " << LowerBoundHoffman << "\n";
    //
    double LowerBoundNifikorov = 1 + mu(0) / (theta(0) - mu(0));
    UpdateBounds(LowerBoundNifikorov, "Nifikorov bound");
    std::cerr << "Nifikorov lower bound = " << LowerBoundNifikorov << "\n";
    //
    double LowerBoundKotolina1 = 1 + mu(0) / (mu(0) - delta(0) + theta(0));
    UpdateBounds(LowerBoundKotolina1, "Kotolina1 bound");
    std::cerr << "Kotolina lower bound 1 = " << LowerBoundKotolina1 << "\n";
    //
    double LowerBoundKotolina2 = 1 + mu(0) / (mu(0) - delta(n-1) + theta(n-1));
    UpdateBounds(LowerBoundKotolina2, "Kotolina2 bound");
    std::cerr << "Kotolina lower bound 2 = " << LowerBoundKotolina2 << "\n";
    //
    for (int m=1; m<=n; m++) {
      double sum1;
      double sum2;
      //
      sum1=0;
      sum2=0;
      for (int i=1; i<=m; i++) {
        sum1 += mu(i - 1);
        sum2 += mu(n+1-i - 1);
      }
      double UnifiedLower1 = 1 + sum1 / (-sum2);
      UpdateBounds(UnifiedLower1, "Unified lower bound 1: m=" + std::to_string(m));
      std::cerr << "Elphick Wocjan unified lower bound 1 = " << UnifiedLower1 << "\n";
      //
      sum1=0;
      sum2=0;
      for (int i=1; i<=m; i++) {
        sum1 += mu(i - 1);
        sum2 += theta(i-1) - mu(i-1);
      }
      double UnifiedLower2 = 1 + sum1 / sum2;
      UpdateBounds(UnifiedLower2, "Unified lower bound 2: m=" + std::to_string(m));
      std::cerr << "Elphick Wocjan unified lower bound 2 = " << UnifiedLower2 << "\n";
      //
      sum1=0;
      sum2=0;
      for (int i=1; i<=m; i++) {
        sum1 += mu(i - 1);
        sum2 += mu(i - 1) - delta(i-1) + theta(i-1);
      }
      double UnifiedLower3 = 1 + sum1 / sum2;
      UpdateBounds(UnifiedLower3, "Unified lower bound 3: m=" + std::to_string(m));
      std::cerr << "Elphick Wocjan unified lower bound 3 = " << UnifiedLower3 << "\n";
      //
      sum1=0;
      sum2=0;
      for (int i=1; i<=m; i++) {
        sum1 += mu(i - 1);
        sum2 += mu(i - 1) - delta(n-i) + theta(n-i);
      }
      double UnifiedLower4 = 1 + sum1 / sum2;
      UpdateBounds(UnifiedLower4, "Unified lower bound 4: m=" + std::to_string(m));
      std::cerr << "Elphick Wocjan unified lower bound 4 = " << UnifiedLower4 << "\n";
    }
    //
    double splus=0;
    double sminus=0;
    for (int i=0; i<n; i++) {
      if (mu(i) > 0)
        splus += mu(i) * mu(i);
      if (mu(i) < 0)
        sminus += mu(i) * mu(i);
    }
    double LowerBoundAndoLin = 1 + splus / sminus;
    UpdateBounds(LowerBoundAndoLin, "Ando/Lin lower bound");
    std::cerr << "Ando Lin lower bound = " << LowerBoundAndoLin << "\n";
    //
    double nPlus = 0;
    double nMinus = 0;
    for (int i=0; i<n; i++) {
      if (mu(i) > 0)
        nPlus++;
      if (mu(i) < 0)
        nMinus++;
    }
    double LowerBoundInertia = 1 + std::max(nMinus / nPlus , nPlus / nMinus);
    UpdateBounds(LowerBoundInertia, "Inertia lower bound");
    std::cerr << "Elphick Wocjan inertial lower bound = " << LowerBoundInertia << "\n";
    //
    std::cerr << "Realizing method=" << RealizingBound << " BestLowerBound=" << BestLowerBound << "\n";
    std::cerr << "Normal termination of ComputeChromaticLowerBounds\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
