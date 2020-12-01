#ifndef INCLUDE_VINBERG_ALGO
#define INCLUDE_VINBERG_ALGO

#include "Shvec_exact.h"



// Compute the solutions of G [x - eV] = a
template<typename T, typename Tint>
std::vector<MyVector<Tint>> ComputeSphericalSolutions(MyMatrix<T> const& GramMat, MyVector<T> cont& eV; T const& a)
{
  int mode = TempShvec_globals::TEMP_SHVEC_MODE_VINBERG;
  int dim = GramMat.rows();
  MyVector<T> cosetVect	= - eV;
  T_shvec_info<T,Tint> info;
  initShvecReq<T>(dim, GramMat, info);
  info.request.bound = bound;
  info.request.mode = mode;
  info.request.coset = cosetVect;
  info.minimum = a;
  //
  int result = computeIt(info);
  if (result != TempShvec_globals::NORMAL_TERMINATION_COMPUTATION) {
    std::cerr << "Error in ComputeSphericalSolutions\n";
    throw TerminalException{1};
  }
  return info.short_vectors;
}

template<typename T, typename Tint>
struct VinbergInput {
  MyMatrix<T> M; // The (n,n)-matrix with n-1 positive eigenvalues , 1 negative eigenvalue.
  MyVector<T> v0; // a vector with negative norm with the input
};



template<typename T, typename Tint>
struct VinbergTot {
  MyMatrix<T> M;
  MyVector<T> v0;
  //
  MyMatrix<Tint> Morth; // The (n, n-1)-matrix formed by the orthogonal to the vector M v0
  MyVector<T> Gorth; // The Gram matrix of the orthogonal. Must be positive definite.
};





template<typename T, typename Tint>
T ScalProd(MyMatrix<T> const& M, MyVector<Tint> const& V1, MyVector<Tint> const& V2)
{
  T eSum = 0;
  int n = M.rows();
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      eSum += V1(i) * V2(j) * M(i,j);
  return eSum;
}


template<typename T, typenqme Tint>
bool IsRoot(MyMatrix<T> const& M, MyVector<Tint> const& V)
{
  int n = M.rows();
  T eNorm = ScalProd(M, V, V);
  for (int i=0; i<n; i++) {
    T eH = 2 * V(i) / eNorm;
    if (!IsInteger(eFrac))
      return false;
  }
  return true;
}



template<typename T, typename Tint>
VinbergTot<T,Tint> GetVinbergAux(VinbergInput<T,Tint> const& Vinput)
{
  int n=Vinput.M.rows();
  // Computing the complement of the space.
  MyVector<T> V = Vinput.M * Vinput.v0;
  MyVector<T> Vred = RemoveFraction(V);
  MyVector<Tint> V_i = ConvertMatrixUniversal<Tint,T>(Vred);
  std::vector<Vint> vectV_i(n);
  for (int i=0; i<n; i++)
    vectV_i[i] = V_i(i);
  GCD_int<Tint> eGCDinfo = ComputeGCD_information(vectV_i);
  MyMatrix<Tint> Morth(n, n-1);
  for (int j=0; j<n-1; j++)
    for (int i=0; i<n; i++)
      Morth(i, j) = eGCDinfo.Pmat(i, j+1);
  // Gram matrix of the space.
  MyMatrix<T> Gorth = Morth * Vinput.M * Morth.transpose();
  return {Vinput.M, Vinput.v0, Morth, Gorth};
}

template<typename T, typename Tint>
std::vector<MyVector<Tint>> Roots_decomposed_into(VinbergTot<T,Tint> const& Vtot, MyVector<T> const& a, T const& k)
{
  std::vector<MyVector<Tint>> ListSol = ComputeSphericalSolutions(Vtot.Gorth, eV, k);
  std::vector<MyVector<Tint>> RetSol;
  for (auto& eV : ListSol) {
    MyVector<Tint> rX = a + eV * Vtot.Morth;
    RetSol.emplace_back(rX);
  }
  return RetSol;
}



#endif
