#include "NumberTheoryBoostGmpInt.h"
#include "MAT_Matrix.h"
int main()
{
  using T = boost::multiprecision::mpq_rational;
  int n = 10;
  int m = 10;
  MyMatrix<T> M(n,m);
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      M(i,j) = rand() % 2;
    }
  }
  int rnk = RankMat(M);
  std::cerr << "rnk=" << rnk << "\n";
}
