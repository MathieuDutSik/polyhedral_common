#ifndef INCLUDE_VINBERG_ALGO
#define INCLUDE_VINBERG_ALGO

#include "Shvec_exact.h"
//#include "MAT_MatrixInt.h"




//
// Finding the closest points
//

// Compute the solutions of G [x - eV] = a
template<typename T, typename Tint>
std::vector<MyVector<Tint>> ComputeSphericalSolutions(const MyMatrix<T>& GramMat, const MyVector<T>& eV, const T& a)
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

//
// Small arithmetic computations for the Vinberg computation
//

template<typename Tint>
Tint En_Quantity(const MyMatrix<Tint>& M)
{
  using Tfield=typename overlying_field<Tint>::field_type;
  Tint eDet = DeternminantMat(M);
  MyMatrix<Tfield> M_field = ConvertMatrixUniversal<Tfield,Tint>(M);
  MyMatrix<Tfield> Minv_field = Inverse(M_field);
  FractionMatrix<Tfield> FrMat = RemoveFractionMatrixPlusCoeff(Minv_field);
  Tint eDen_T = UniversalTypeConversion<Tint,Tfield>(FrMat.TheMult);
  Tint TheEn = T_abs(eDet / eDen_T);
  return TheEn;
}


template<typename Tint>
std::vector<Tint> Get_root_lengths(const MyMatrix<Tint>& M)
{
  Tint TheEn = En_Quantity(M);
  Tint limit = 2 * TheEn;
  std::vector<Tint> root_lengths;
  Tint k = 1;
  while (true) {
    Tint res = ResInt(limit, k);
    if (res == 0)
      root_lengths.push_back(k);
    if (k >= limit)
      break;
    k++;
  }
  return root_lengths;
}


//
// The fundamental data type
//


template<typename T, typename Tint>
struct VinbergInput {
  MyMatrix<T> G; // The (n,n)-matrix with n-1 positive eigenvalues , 1 negative eigenvalue.
  MyVector<T> v0; // a vector with negative norm with the input
};



template<typename T, typename Tint>
struct VinbergTot {
  MyMatrix<T> G;
  MyVector<T> v0;
  MyVector<Tint> V_i;
  MyVector<Tint> Vtrans;
  MyMatrix<Tint> Mbas;
  MyMatrix<Tint> MbasInv;
  //
  MyMatrix<Tint> Morth; // The (n, n-1)-matrix formed by the orthogonal to the vector M v0
  Tint eDet; // The determinant of the matrix.
  MyMatrix<T> Gorth; // The Gram matrix of the orthogonal. Must be positive definite.
  MyMatrix<T> GM_iGorth; // The inverse of the coefficient for the computation.
  std::vector<Tint> root_lengths;
};





template<typename T, typename Tint>
T ScalProd(const MyMatrix<T>& M, const MyVector<Tint>& V1, const MyVector<Tint>& V2)
{
  T eSum = 0;
  int n = M.rows();
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      eSum += V1(i) * V2(j) * M(i,j);
  return eSum;
}


template<typename T, typename Tint>
bool IsRoot(const MyMatrix<T>& M, const MyVector<Tint>& V)
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
VinbergTot<T,Tint> GetVinbergAux(const MyMatrix<Tint>& G, const MyVector<Tint>& v0)
{
  int n=G.rows();
  // Computing the complement of the space.
  MyVector<Tint> V = G * v0;
  std::vector<Tint> vectV(n);
  for (int i=0; i<n; i++)
    vectV[i] = V(i);
  GCD_int<Tint> eGCDinfo = ComputeGCD_information(vectV);
  std::vector<int> ListZer(n-1);
  for (int j=0; j<n-1; j++)
    ListZer[j] = j + 1;
  MyMatrix<Tint> Morth = SelectColumn(eGCDinfo.Pmat, ListZer);
  MyMatrix<Tint> M = ConcatenateMatVec(Morth, V);
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
  std::vector<Tint> root_lengths = Get_root_lengths(Vinput.G);
  return {Vinput.G, Vinput.v0, V_i, Vtrans, Mbas, MbasInv, Morth, eDet, Gorth, GM_iGorth, root_lengths};
}


template<typename T, typename Tint>
struct IterateRootDecompositions {
private:
  std::unordered_map<Tint,int> candidates;
  const VinbergTot<T,Tint>& Vtot;
  Tint len_sW;
  MyVector<Tint> cand_a(const int& n) const {
    size_t len_sW = Vtot.W.size();
    int res = n % len_sW;
    int q = n / len_sW;
    return Vtot.W[res] + q * Vtot.v0;
  }
  Tint get_k() const {
    bool we_found = false;
    double minval_d;
    Tint kfind;
    for (auto & k : Vtot.root_lengths) {
      Tint val = - V.dot(cand_a(candidates[k]));
      double k_d = sqrt(UniversalTypeConversion<double,Tint>(val));
      double val_d = UniversalTypeConversion<double,Tint>(val) / k_d;
      if (!we_found) {
        we_found = true;
        minval_d = val_d;
        kfind = k;
      } else {
        if (val_d < minval_d) {
          minval_d = val_d;
          kfind = k;
        }
      }
    }
    return k;
  }
public:
  IterateRootDecompositions(const VinbergTot<T,Tint>& Vtot) : Vtot(Vtot) {
    for (auto & k : Vtot.root_lengths)
      candidates[k] = 1;
  }
  std::pair<MyVector<Tint>, Tint> get_cand() {
    Tint k = get_k();
    std::cerr << "IterateRootDecomposition, k=" << k << "\n";
    int val = candidates[k];
    candidates[k] = val + 1;
    MyVector<Tint> V = cand_a(val);
    std::cerr << "IterateRootDecomposition, cand_a(candidates[k])=" << V << "\n";
    return {V, k};
  }
};










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
std::vector<MyVector<Tint>> Roots_decomposed_into(const VinbergTot<T,Tint>& Vtot, const MyVector<T>& a, const T& n)
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
bool is_FundPoly(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
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


template<typename T, typename Tint>
std::vector<MyVector<Tint>> GetIntegerPoints(const MyMatrix<Tint>& m)
{
  size_t n_rows = m.rows();
  size_t n_cols = m.cols();
  MyVector<Tint> negative = ZeroVector<Tint>(n_cols);
  MyVector<Tint> positive = ZeroVector<Tint>(n_cols);
  for (size_t i_col=0; i_col<n_cols; i_col++) {
    for (size_t i_row=0; i_row<n_rows; i_row++) {
      Tint val = m(i_row,i_col);
      if (val < 0)
        negative(i_col) += val;
      if (val > 0)
        positive(i_col) += val;
    }
  }
  std::vector<int> ListSize(n_cols);
  for (size_t i_col=0; i_col<n_cols; i_col++) {
    int val1 = UniversalTypeConversion<int,Tint>(negative(i_col));
    int val2 = UniversalTypeConversion<int,Tint>(positive(i_col));
    int len = val2 + 1 - val1;
    ListSize[i_col] = len;
  }
  using Tfield=typename overlying_field<Tint>::field_type;
  MyMatrix<Tfield> M_field = ConvertMatrixUniversal<Tfield,Tint>(m);
  MyMatrix<Tfield> Minv_field = Inverse(M_field);
  FractionMatrix<Tfield> FrMat = RemoveFractionMatrixPlusCoeff(Minv_field);
  Tint eDen = UniversalTypeConversion<Tint,Tfield>(FrMat.TheMult);
  if (eDen_T
  MyMatrix<Tint> Comat = ConvertMatrixUniversal<Tint,Tfield>(FrMat.TheMat);
  auto ParallelepipedContains=[&](const MyVector<Tint>& V) -> bool {
    MyVector<Tint> Q = V * Comat;
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      if (Q(i_col) < 0)
        return false;
      if (Q(i_col) >= eDen)
        return false;
    }
    return true;
  };
  std::vector<MyVector<Tint>> ListPoint;
  BlockIterationMultiple BlIter(ListSize);
  for (auto const& eVect : BlIter) {
    MyVector<Tint> ePoint(n_cols);
    for (int i_col=0; i_col<n_cols; i_col++)
      ePoint(i_col) = negative(i_col) + eVect[i_col];
    if (ParallelepipedContains(ePoint))
      ListPoint.push_back(ePoint);
  }
  return ListPoint;
}


template<typename T, typename Tint>
std::vector<MyVector<Tint>> FundCone(const VinbergTot<T,Tint>& Vtot)
{
  std::vector<MyVector<Tint>> V1_roots;
  size_t n = Vtot.G.rows();
  MyVector<T> a = ZeroVector<T>(n);
  for (auto & k : Vtot.root_lengths) {
    for (const MyVector<Tint>& root_cand : Roots_decomposed_into(Vtot, a, k)) {
      if (IsRoot(Vtot.G, root_cand))
        V1_roots.push_back(root_cand);
    }
  }
  return V1_roots;
}


template<typename T, typename Tint>
std::vector<MyVector<Tint>> FindRoots(const VinbergTot<T,Tint>& Vtot)
{
  std::vector<MyVector<Tint>> ListRoot = FundCone(Vtot);

  IterateRootDecompositions<T,Tint> iter(Vtot);
  while (true) {
    const std::pair<MyVector<Tint>,Tint> pair = iter.get_cand();
    const MyVector<Tint>& a = pair.first;
    const Tint& k = pair.second;
    std::cerr << "  NextRoot a=" << a << " k=" << k << "\n";
    for (const MyVector<Tint>& root_cand : Roots_decomposed_into(Vtot, a, n)) {
      if (IsRoot(Vtot.G, root_cand)) {
        ListRoot.push_back(root_cand);
        if (is_FundPoly(Vtot, ListRoot)) {
          return ListRoot;
        }
      }
    }
  }
  std::cerr << "Should never reach that stage\n";
  throw TerminalException{1};
}









#endif
