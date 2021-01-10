#ifndef INCLUDE_VINBERG_ALGO
#define INCLUDE_VINBERG_ALGO

#include "Shvec_exact.h"
//#include "MAT_MatrixInt.h"



// Compute the solutions of G [x - eV] = a
template<typename T, typename Tint>
std::vector<MyVector<Tint>> ComputeSphericalSolutions(MyMatrix<T> const& GramMat, MyVector<T> const& eV, T const& a)
{
  int mode = TempShvec_globals::TEMP_SHVEC_MODE_VINBERG;
  int dim = GramMat.rows();
  MyVector<T> cosetVect	= - eV;
  T_shvec_info<T,Tint> info;
  initShvecReq<T>(dim, GramMat, info);
  info.request.bound = a;
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
  MyMatrix<T> G; // The (n,n)-matrix with n-1 positive eigenvalues , 1 negative eigenvalue.
  MyVector<T> v0; // a vector with negative norm with the input
};



template<typename T, typename Tint>
struct VinbergTot {
  MyMatrix<T> G;
  MyVector<Tint> V_i;
  MyVector<Tint> Vtrans;
  MyMatrix<Tint> Mbas;
  MyMatrix<Tint> MbasInv;
  //
  MyMatrix<Tint> Morth; // The (n, n-1)-matrix formed by the orthogonal to the vector M v0
  Tint eDet; // The determinant of the matrix.
  MyMatrix<T> Gorth; // The Gram matrix of the orthogonal. Must be positive definite.
  MyMatrix<T> GM_iGorth; // The inverse of the coefficient for the computation.
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


template<typename T, typename Tint>
bool IsRoot(MyMatrix<T> const& M, MyVector<Tint> const& V)
{
  int n = M.rows();
  T eNorm = ScalProd(M, V, V);
  for (int i=0; i<n; i++) {
    T eFrac = 2 * V(i) / eNorm;
    if (!IsInteger(eFrac))
      return false;
  }
  return true;
}



template<typename T, typename Tint>
VinbergTot<T,Tint> GetVinbergAux(VinbergInput<T,Tint> const& Vinput)
{
  int n=Vinput.G.rows();
  // Computing the complement of the space.
  MyVector<T> V = Vinput.G * Vinput.v0;
  MyVector<T> Vred = RemoveFractionVector(V);
  MyVector<Tint> V_i = ConvertMatrixUniversal<Tint,T>(Vred);
  std::vector<Tint> vectV_i(n);
  for (int i=0; i<n; i++)
    vectV_i[i] = V_i(i);
  GCD_int<Tint> eGCDinfo = ComputeGCD_information(vectV_i);
  std::vector<int> ListZer(n-1);
  for (int j=0; j<n-1; j++)
    ListZer[j] = j + 1;
  MyMatrix<Tint> Morth = SelectColumn(eGCDinfo.Pmat, ListZer);
  MyMatrix<Tint> M = ConcatenateMatVec(Morth, V_i);
  // The dterminant. The scalar tell us how much we need to the quotient.
  // We will need to consider the vectors k (V_i / eDet) for k=1, 2, 3, ....
  Tint eDet = T_abs(DeterminantMat(M));
  // We want to find a vector v such that V_i = (det) v + Morth Z^{n-1}
  auto GetVect = [&]() -> MyVector<Tint> {
    for (int i=0; i<n; i++) {
      for (int j=0; j<2; j++) {
        int eps = -1 + 2 * j;
        MyVector<Tint> V = V_i;
        V(i) -= eps;
        SolMatResult<Tint> Solu=SolutionMat(Morth, V);
        if (Solu.result) {
          MyVector<Tint> Vret = ZeroVector<Tint>(n);
          Vret(i) = eps;
          return Vret;
        }
      }
    }
    std::cerr << "Failed to find the right vector\n";
    throw TerminalException{1};
  };
  MyVector<Tint> Vtrans = GetVect();
  MyMatrix<Tint> Mbas = ConcatenateMatVec(Morth, Vtrans);
  MyMatrix<Tint> MbasInv = Inverse(Mbas);

  // Gram matrix of the space.
  MyMatrix<T> Gorth = Morth * Vinput.G * Morth.transpose();
  MyMatrix<T> GorthInv = Inverse(Gorth);
  // Computing the side comput
  MyMatrix<T> GM_iGorth = Vinput.G * Morth * GorthInv;
  return {Vinput.G, V_i, Vtrans, Mbas, MbasInv, Morth, eDet, Gorth, GM_iGorth};
}


/*
  We look for the solutions of (a+v , a+v) = k
  with v in the Morth space.
  (a, a) + 2 (a, v) + (v,v) = n
   v = M w  with  w in Z^{n-1}
  2 a^t G Mw + w^t {M^t G M} w = n - (a,a)
  2 w Gorth sV + w^t Gorth w = n -(a,a)
  (w + sV)^t Gorth (w + sV) = n - (a,a) + sV^t Gorth sV

 */
template<typename T, typename Tint>
std::vector<MyVector<Tint>> Roots_decomposed_into(VinbergTot<T,Tint> const& Vtot, MyVector<T> const& a, T const& n)
{
  MyVector<T> sV = a * Vtot.GM_iGorth;
  T normi = n - a.dot(Vtot.G * a) + sV.dot(Vtot.Gorth * sV);
  MyVector<T> eV = -sV;
  std::vector<MyVector<Tint>> ListSol = ComputeSphericalSolutions<T,Tint>(Vtot.Gorth, eV, normi);
  std::vector<MyVector<Tint>> RetSol;
  for (auto& eV : ListSol) {
    MyVector<Tint> rX = a + eV * Vtot.Morth;
    RetSol.emplace_back(rX);
  }
  return RetSol;
}


template<typename T, typename Tint>
bool is_FundPoly(VinbergTot<T,Tint> const& Vtot, std::vector<MyVector<Tint>> const& ListRoot)
{
  int n_root = ListRoot.size();
  MyMatrix<T> M(n_root, n_root);
  for (int i_root=0; i_root<n_root; i_root++) {
    MyVector<T> eVG = ListRoot[i_root] * Vtot.G;
    for (int j_root=0; j_root<n_root; j_root++) {
      T eScal = eVG.dot(ListRoot[j_root]);
      M(i_root, j_root) = eScal;
    }
  }
  auto weight=[&](int i, int j) -> int {
    T aII = M(i,i);
    T aJJ = M(j,j);
    T aIJ = M(i,j);
    T cos2 = (aIJ * aIJ) / (aII * aJJ);
    if (cos2 == 0)
      return 2;
    if (cos2 == 1/4)
      return 3;
    if (cos2 == 1/2)
      return 4;
    if (cos2 == 3/4)
      return 6;
    if (cos2 == 1)
      return 0;
    if (cos2 > 1)
      return 1;
    std::cerr << "coxiter.py ERROR: cosine " << cos2 << "\n";
    throw TerminalException{1};
  };
  int d = Vtot.G.rows();
  std::string rnd_str = random_string(20);
  std::string FileI = "/tmp/CoxIter_" + rnd_str + ".input";
  std::string FileO = "/tmp/CoxIter_" + rnd_str + ".out";
  {
    std::ofstream os(FileI);
    os << n_root << " " << d << "\n";
    for (int i=0; i<n_root; i++)
      for (int j=0; j<i; j++)
        if (M(i,j) != 0)
          os << (j+1) << " " << (i+1) << " " << weight(i, j) << "\n";
    os << "\n";
  }
  //
  // Running the CoxIter program
  //
  std::string eCommand = "coxiter";
  std::string opt = "-fv";
  eCommand += " " + opt;
  eCommand += " < " + FileI + " > " + FileO;
  std::cerr << "eCommand=" << eCommand << "\n";
  int iret=system(eCommand.c_str());
  if (iret == -1) {
    printf("Oh dear, something went wrong with glpsol! %s\n", strerror(errno));
    throw TerminalException{1};
  }
  //
  // Reading the output
  //
  std::vector<std::string> RESUL;
  {
    std::ifstream INfs(FileO);
    std::string line;
    while (getline(INfs, line))
      RESUL.push_back(line);
  }
  bool IsFiniteCovolume=false;
  std::string question = "Finite covolume";
  std::string answer = "yes";
  for (auto & eLine : RESUL) {
    std::vector<std::string> LStr1 = STRING_Split(eLine, question);
    std::vector<std::string> LStr2 = STRING_Split(eLine, answer);
    if (LStr1.size() > 1 && LStr2.size() > 1)
      IsFiniteCovolume = true;
  }
  return IsFiniteCovolume;
}





#endif
