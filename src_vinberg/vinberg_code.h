#ifndef INCLUDE_VINBERG_ALGO
#define INCLUDE_VINBERG_ALGO

#include "Shvec_exact.h"
//#include "MAT_MatrixInt.h"
#include "POLY_cddlib.h"


#define DEBUG_VINBERG


//
// Finding the closest points
//

// Compute the solutions of G [x - eV] = a
template<typename T, typename Tint>
std::vector<MyVector<Tint>> ComputeSphericalSolutions(const MyMatrix<T>& GramMat, const MyVector<T>& eV, const T& a)
{
  std::cerr << "ComputeSphericalSolutions, step 1\n";
  std::cerr << "GramMat=\n";
  WriteMatrix(std::cerr, GramMat);
  std::cerr << "eV=\n";
  WriteVector(std::cerr, eV);
  std::cerr << "a=" << a << "\n";
  int mode = TempShvec_globals::TEMP_SHVEC_MODE_VINBERG;
  int dim = GramMat.rows();
  MyVector<T> cosetVect	= - eV;
  T_shvec_request<T> request;
  T_shvec_info<T,Tint> info;
  initShvecReq<T>(dim, GramMat, request, info);
  request.bound = a;
  request.mode = mode;
  request.coset = cosetVect;
  request.only_exact_norm = true; // We want only the vectors of norm exactly a.
  info.minimum = a;
  //
  int result = computeIt(request, a, info);
  if (result != TempShvec_globals::NORMAL_TERMINATION_COMPUTATION) {
    std::cerr << "Error in ComputeSphericalSolutions\n";
    throw TerminalException{1};
  }
#ifdef DEBUG_VINBERG
  std::cerr << "|info.short_vectors|=" << info.short_vectors.size() << "\n";
  for (auto & fV : info.short_vectors) {
    T sum=0;
    for (int i=0; i<dim; i++)
      for (int j=0; j<=i; j++) {
        T eMult = 2;
        if (i == j)
          eMult = 1;
        T val1 = eV(i) - fV(i);
        T val2 = eV(j) - fV(j);
        sum += val1 * val2 * GramMat(i,j) * eMult;
      }
    if (sum != a) {
      std::cerr << "We should have vectors of norm exactly a\n";
      std::cerr << "sum=" << sum << " a=" << a << "\n";
      throw TerminalException{1};
    }
  }
#endif
  return info.short_vectors;
}

//
// Small arithmetic computations for the Vinberg computation
//

template<typename Tint>
Tint En_Quantity(const MyMatrix<Tint>& M)
{
  using Tfield=typename overlying_field<Tint>::field_type;
  size_t dim = M.rows();
  Tint eDet = DeterminantMat(M);
  MyMatrix<Tfield> M_field = UniversalMatrixConversion<Tfield,Tint>(M);
  MyMatrix<Tfield> Minv_field = Inverse(M_field);
  Tfield eDet_Tfield = DeterminantMat(M_field);
  Tint eDet_Tint = UniversalScalarConversion<Tint,Tfield>(eDet_Tfield);
  MyMatrix<Tint> M_adjoint = UniversalMatrixConversion<Tint,Tfield>(eDet_Tfield * Minv_field);
  std::vector<Tint> ListX;
  ListX.reserve(dim * dim);
  for (size_t i=0; i<dim; i++)
    for (size_t j=0; j<dim; j++)
      ListX.push_back(M_adjoint(i,j));
  Tint gcd = ComputeGCD_information(ListX).gcd;
  Tint TheEn = T_abs(eDet_Tint) / T_abs(gcd);
  return TheEn;
}


template<typename Tint>
std::vector<Tint> Get_root_lengths(const MyMatrix<Tint>& M)
{
  Tint TheEn = En_Quantity(M);
  std::cerr << "TheEn=" << TheEn << "\n";
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
struct VinbergTot {
  MyMatrix<Tint> G;
  MyMatrix<T> G_T;
  MyVector<Tint> v0;
  MyVector<Tint> V_i;
  MyVector<Tint> Vtrans;
  MyMatrix<Tint> Mbas;
  //  MyMatrix<Tint> MbasInv_T;
  //
  MyMatrix<Tint> Morth; // The (n, n-1)-matrix formed by the orthogonal to the vector M v0
  MyMatrix<T> Morth_T; // The (n, n-1)-matrix formed by the orthogonal to the vector M v0
  Tint eDet; // The determinant of the matrix.
  MyMatrix<Tint> Gorth; // The Gram matrix of the orthogonal. Must be positive definite.
  MyMatrix<T> Gorth_T; // The Gram matrix of the orthogonal. Must be positive definite.
  MyMatrix<T> GM_iGorth; // The inverse of the coefficient for the computation.
  std::vector<MyVector<Tint>> W;
  std::vector<Tint> root_lengths;
};





template<typename Tint>
Tint ScalProd(const MyMatrix<Tint>& M, const MyVector<Tint>& V1, const MyVector<Tint>& V2)
{
  Tint eSum = 0;
  int n = M.rows();
  for (int i=0; i<n; i++)
    for (int j=0; j<=i; j++) {
      Tint eMult = 2;
      if (i == j)
        eMult = 1;
      eSum += V1(i) * V2(j) * M(i,j) * eMult;
    }
  return eSum;
}


template<typename T, typename Tint>
bool IsRoot(const MyMatrix<Tint>& M, const MyVector<Tint>& V)
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


template<typename Tint>
std::vector<MyVector<Tint>> GetIntegerPoints_V1(const MyMatrix<Tint>& m)
{
  std::cerr << "GetIntegerPoints m=\n";
  WriteMatrix(std::cerr, m);
  Tint det = DeterminantMat(m);
  std::cerr << "det(m)=" << det << "\n";
  size_t n_rows = m.rows();
  size_t n_cols = m.cols();
  std::cerr << "n_rows=" << n_rows << " n_cols=" << n_cols << "\n";
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
  std::cerr << "negative =";
  for (size_t i_col=0; i_col<n_cols; i_col++)
    std::cerr << " " << negative(i_col);
  std::cerr << "\n";
  std::cerr << "positive =";
  for (size_t i_col=0; i_col<n_cols; i_col++)
    std::cerr << " " << positive(i_col);
  std::cerr << "\n";
  std::vector<int> ListSize(n_cols);
  for (size_t i_col=0; i_col<n_cols; i_col++) {
    int val1 = UniversalScalarConversion<int,Tint>(negative(i_col));
    int val2 = UniversalScalarConversion<int,Tint>(positive(i_col));
    int len = val2 + 1 - val1;
    ListSize[i_col] = len;
  }
  std::cerr << "We have ListSizes\n";
  using Tfield=typename overlying_field<Tint>::field_type;
  MyMatrix<Tfield> M_field = UniversalMatrixConversion<Tfield,Tint>(m);
  MyMatrix<Tfield> Minv_field = Inverse(M_field);
  std::cerr << "We have Minv_field\n";
  FractionMatrix<Tfield> FrMat = RemoveFractionMatrixPlusCoeff(Minv_field);
  std::cerr << "We have FrMat\n";
  Tint eDen = UniversalScalarConversion<Tint,Tfield>(FrMat.TheMult);
  std::cerr << "We have eDen=" << eDen << "\n";
  if (eDen < 0) {
    std::cerr << "We should have eDen > 0. eDen=" << eDen << "\n";
    throw TerminalException{1};
  }
  MyMatrix<Tint> Comat = UniversalMatrixConversion<Tint,Tfield>(FrMat.TheMat);
  std::cerr << "We have |Comat|=" << Comat.rows() << " / " << Comat.cols() << "\n";
  WriteMatrix(std::cerr, Comat);
  auto ParallelepipedContains=[&](const MyVector<Tint>& V) -> bool {
    std::cerr << "ParallelepipedContains, begin\n";
    MyVector<Tint> Q = V.transpose() * Comat;
    std::cerr << "ParallelepipedContains, we have Q\n";
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
  for (const auto& eVect : BlIter) {
    MyVector<Tint> ePoint(n_cols);
    std::cerr << "ePoint =";
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      ePoint(i_col) = negative(i_col) + eVect[i_col];
      std::cerr << " " << ePoint(i_col);
    }
    std::cerr << "\n";
    if (ParallelepipedContains(ePoint))
      ListPoint.push_back(ePoint);
  }
  std::cerr << "GetIntegerPoints |ListRet|=" << ListPoint.size() << "\n";
  std::cerr << "GetIntegerPoints ListRet=\n";
  for (auto & eVect : ListPoint) {
    WriteVector(std::cerr, eVect);
  }
  std::cerr << "Leaving GetIntegerPoints\n";
  return ListPoint;
}


template<typename Tint>
std::vector<MyVector<Tint>> GetIntegerPoints(const MyMatrix<Tint>& m)
{
  return ComputeTranslationClasses<Tint,Tint>(m);
}



template<typename T, typename Tint>
VinbergTot<T,Tint> GetVinbergAux(const MyMatrix<Tint>& G, const MyVector<Tint>& v0)
{
  std::cerr << "GetVinbergAux, step 1\n";
  int n=G.rows();
  std::cerr << "GetVinbergAux, step 1.1\n";
  // Computing the complement of the space.
  MyMatrix<T> G_T = UniversalMatrixConversion<T,Tint>(G);
  std::cerr << "GetVinbergAux, step 1.2\n";
  MyVector<Tint> V = G * v0;
  std::cerr << "GetVinbergAux, step 1.3\n";
  std::vector<Tint> vectV(n);
  for (int i=0; i<n; i++)
    vectV[i] = V(i);
  std::cerr << "GetVinbergAux, step 2\n";
  GCD_int<Tint> eGCDinfo = ComputeGCD_information(vectV);
  std::vector<int> ListZer(n-1);
  for (int j=0; j<n-1; j++)
    ListZer[j] = j + 1;
  std::cerr << "GetVinbergAux, step 2.1\n";
  MyMatrix<Tint> Morth = SelectColumn(eGCDinfo.Pmat, ListZer);
  std::cerr << "GetVinbergAux, step 2.2\n";
  MyMatrix<T> Morth_T = UniversalMatrixConversion<T,Tint>(Morth);
  std::cerr << "GetVinbergAux, step 2.3\n";
  MyMatrix<Tint> M = ConcatenateMatVec_Tr(Morth, V);
  std::cerr << "GetVinbergAux, step 2.4\n";
  MyMatrix<Tint> M2 = ConcatenateMatVec_Tr(Morth, v0);
  std::cerr << "GetVinbergAux, step 2.5\n";
  std::vector<MyVector<Tint>> W = GetIntegerPoints(M2);
  std::cerr << "GetVinbergAux, step 3\n";
  // The dterminant. The scalar tell us how much we need to the quotient.
  // We will need to consider the vectors k (V_i / eDet) for k=1, 2, 3, ....
  Tint eDet = T_abs(DeterminantMat(M));
  // We want to find a vector v such that V_i = (det) v + Morth Z^{n-1}
  std::cerr << "GetVinbergAux, step 4\n";
  auto GetVect = [&]() -> MyVector<Tint> {
    for (int i=0; i<n; i++) {
      for (int j=0; j<2; j++) {
        int eps = -1 + 2 * j;
        MyVector<Tint> Vwork = V;
        Vwork(i) -= eps;
        SolMatResult<Tint> Solu=SolutionMat(Morth, Vwork);
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
  std::cerr << "GetVinbergAux, step 5\n";
  MyMatrix<Tint> Mbas = ConcatenateMatVec_Tr(Morth, Vtrans);
  std::cerr << "GetVinbergAux, step 5.1\n";
  //  MyMatrix<T> MbasInv_T = Inverse(UniversalMatrixConversion<T,Tint>(Mbas));
  std::cerr << "GetVinbergAux, step 6\n";

  // Gram matrix of the space.
  MyMatrix<Tint> Gorth = Morth.transpose() * G * Morth;
  std::cerr << "GetVinbergAux, step 6.1\n";
  MyMatrix<T> Gorth_T = UniversalMatrixConversion<T,Tint>(Gorth);
  std::cerr << "GetVinbergAux, step 6.2\n";
  MyMatrix<T> GorthInv = Inverse(Gorth_T);
  std::cerr << "GetVinbergAux, step 7\n";
  // Computing the side comput
  MyMatrix<T> GM_iGorth = G_T * UniversalMatrixConversion<T,Tint>(Morth) * GorthInv;
  std::vector<Tint> root_lengths = Get_root_lengths(G);
  std::cerr << "s.root_lengths =";
  for (auto & eVal : root_lengths)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  std::cerr << "GetVinbergAux, step 8\n";
  //  return {G, G_T, v0, V, Vtrans, Mbas, MbasInv_T, Morth, Morth_T, eDet, Gorth, Gorth_T, GM_iGorth, W, root_lengths};
  return {G, G_T, v0, V, Vtrans, Mbas, Morth, Morth_T, eDet, Gorth, Gorth_T, GM_iGorth, W, root_lengths};
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
    double minval_d=0;
    Tint kfind;
    for (auto & k : Vtot.root_lengths) {
      MyVector<Tint> V2 = cand_a(candidates.at(k));
      Tint val = - Vtot.v0.dot(V2);
      double k_d = sqrt(UniversalScalarConversion<double,Tint>(val));
      double val_d = UniversalScalarConversion<double,Tint>(val) / k_d;
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
    return kfind;
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
  std::cerr << "Roots_decomposed_into, step 1\n";
  MyVector<T> sV = a.transpose() * Vtot.GM_iGorth;
  std::cerr << "Roots_decomposed_into, step 2\n";
  T normi = n - a.dot(Vtot.G_T * a) + sV.dot(Vtot.Gorth_T * sV);
  std::cerr << "Roots_decomposed_into, step 3\n";
  MyVector<T> eV = -sV;
  std::cerr << "Roots_decomposed_into, step 4\n";
  std::vector<MyVector<Tint>> ListSol = ComputeSphericalSolutions<T,Tint>(Vtot.Gorth_T, eV, normi);
  std::cerr << "|ListSol|=" << ListSol.size() << "\n";
  std::cerr << "Roots_decomposed_into, step 5\n";
  std::vector<MyVector<Tint>> RetSol;
  std::cerr << "Roots_decomposed_into, step 6\n";
  for (auto& eV_Tint : ListSol) {
    //    std::cerr << "Roots_decomposed_into, step 6.1\n";
    MyVector<T> eV_T = UniversalVectorConversion<T,Tint>(eV_Tint);
    //    std::cerr << "Roots_decomposed_into, step 6.2\n";
    std::cerr << "|a|=" << a.size() << " |eV_T|=" << eV_T.size() << " |Morth_T|=" << Vtot.Morth_T.rows() << " / " << Vtot.Morth_T.cols() << "\n";
    MyVector<T> rX_T = a + Vtot.Morth_T * eV_T;
    //    std::cerr << "Roots_decomposed_into, step 6.3\n";
    MyVector<Tint> rX = UniversalVectorConversion<Tint,T>(rX_T);
    //    std::cerr << "Roots_decomposed_into, step 6.4\n";
    RetSol.emplace_back(std::move(rX));
  }
  std::cerr << "Roots_decomposed_into, step 7\n";
  return RetSol;
}








template<typename T, typename Tint>
bool is_FundPoly(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
{
  int n_root = ListRoot.size();
  MyMatrix<T> M(n_root, n_root);
  for (int i_root=0; i_root<n_root; i_root++) {
    MyVector<Tint> eVG = ListRoot[i_root] * Vtot.G;
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


template<typename T>
MyVector<T> SignCanonicalizeVector(const MyVector<T>& V)
{
  int len = V.size();
  for (int u=0; u<len; u++) {
    if (V(u) > 0)
      return V;
    if (V(u) < 0)
      return -V;
  }
  std::cerr << "Error in SignCanonicalizeVector\n";
  throw TerminalException{1};
}


template<typename T, typename Tint>
std::vector<MyVector<Tint>> FundCone(const VinbergTot<T,Tint>& Vtot)
{
  //
  // First building the initial set of roots
  //
  std::cerr << "FundCone, step 1\n";
  std::vector<MyVector<Tint>> V1_roots;
  size_t n = Vtot.G.rows();
  MyVector<T> a = ZeroVector<T>(n);
  std::cerr << "FundCone, step 1.1\n";
  for (auto & k : Vtot.root_lengths) {
    T k_T = k;
    std::cerr << " k=" << k << "\n";
    std::set<MyVector<Tint>> set;
    for (const MyVector<Tint>& root_cand : Roots_decomposed_into<T,Tint>(Vtot, a, k_T)) {
      if (IsRoot<T,Tint>(Vtot.G, root_cand)) {
        MyVector<Tint> root_can = SignCanonicalizeVector(root_cand);
        set.insert(root_can);
      }
    }
    for (auto & eV : set)
      V1_roots.push_back(eV);
    std::cerr << "|set|=" << set.size() << " |V1_roots|=" << V1_roots.size() << "\n";
  }
  std::cerr << "FundCone, step 2\n";
  //
  // Selecting a basis of roots as a starting point
  // (Not sure if that initial family is full dimensional or not. We assume full dimensional))
  //
  auto f=[&](MyMatrix<T> & M, size_t eRank, size_t iRow) -> void {
    for (size_t i=0; i<n; i++)
      M(eRank, i) = UniversalScalarConversion<T,Tint>(V1_roots[iRow](i));
  };
  size_t nbRow = V1_roots.size();
  size_t nbCol = n;
  SelectionRowCol<T> eSelect=TMat_SelectRowCol_Kernel<T>(nbRow, nbCol, f);
  size_t TheRank = eSelect.TheRank;
  std::cerr << "eSelect.TheRank=" << TheRank << "\n";
  std::cerr << "FundCone, step 3\n";
  Face selected(nbRow);
  std::vector<MyVector<Tint>> SelectedRoots;
  for (auto & idx : eSelect.ListRowSelect) {
    selected[idx] = 1;
    const MyVector<Tint>& uRoot = V1_roots[idx];
    MyVector<Tint> Vprod = Vtot.G * uRoot;
    size_t n_plus = 0, n_minus = 0;
    for (auto& eRoot : SelectedRoots) {
      Tint scal = Vprod.dot(eRoot);
      if (scal > 0)
        n_plus++;
      if (scal < 0)
        n_minus++;
    }
    if (n_minus > n_plus)
      SelectedRoots.push_back(uRoot);
    else
      SelectedRoots.push_back(-uRoot);
  }
  std::cerr << "FundCone, step 4\n";
  std::cerr << "Rank(SelectedRoots)=" << RankMat(MatrixFromVectorFamily(SelectedRoots)) << "\n";
  //
  // Now iterating over the roots.
  //
  auto get_facets=[&]() -> MyMatrix<T> {
    size_t n_root = SelectedRoots.size();
    MyMatrix<T> Mroot(n_root, TheRank);
    for (size_t i_root=0; i_root<n_root; i_root++)
      for (size_t i=0; i<TheRank; i++) {
        int iCol = eSelect.ListColSelect[i];
        Mroot(i_root, i) = UniversalScalarConversion<T,Tint>(SelectedRoots[i_root](iCol));
      }
    std::cerr << "Mroot=\n";
    WriteMatrix(std::cerr, Mroot);
    std::cerr << "Before cdd::DualDescription\n";
    return cdd::DualDescription(Mroot); // maybe use another dual description function
  };
  MyMatrix<T> FAC = get_facets();
  std::cerr << "FundCone, step 5\n";
  auto insert_root=[&](const MyVector<Tint>& V) -> void {
    size_t n_plus = 0;
    size_t n_minus = 0;
    size_t n_fac = FAC.rows();
    const MyVector<T> V_T = UniversalVectorConversion<T,Tint>(V);
    for (size_t i_fac=0; i_fac<n_fac; i_fac++) {
      T scal = 0;
      for (size_t i=0; i<TheRank; i++) {
        int iCol = eSelect.ListColSelect[i];
        scal += FAC(i_fac,i) * V_T(iCol);
      }
      if (scal > 0)
        n_plus++;
      if (scal < 0)
        n_minus++;
    }
    if (n_plus == 0 || n_minus == 0) // The inequality is valid. Exiting
      return;
    auto get_root=[&]() -> MyVector<Tint> {
      if (n_plus > n_minus) // We look for the vector that splits most
        return -V;
      else
        return V;
    };
    MyVector<Tint> Vsel = get_root();
    SelectedRoots.push_back(Vsel);
    FAC = get_facets();
    n_fac = FAC.rows();
    std::vector<MyVector<T>> ListRowFAC;
    for (size_t i_fac=0; i_fac<n_fac; i_fac++)
      ListRowFAC.push_back(GetMatrixRow(FAC, i_fac));
    std::vector<MyVector<Tint>> TheSelect;
    for (auto & eRoot : SelectedRoots) {
      MyVector<T> eRootRestr_T(TheRank);
      for (size_t i=0; i<TheRank; i++) {
        int iCol = eSelect.ListColSelect[i];
        eRootRestr_T(i) = UniversalScalarConversion<T,Tint>(eRoot(iCol));
      }
      std::vector<size_t> TheIncd;
      for (size_t i_fac=0; i_fac<n_fac; i_fac++) {
        T scal = eRootRestr_T.dot(ListRowFAC[i_fac]);
        if (scal == 0)
          TheIncd.push_back(i_fac);
      }
      size_t eRank = TMat_SelectRowCol_subset(FAC, TheIncd).TheRank;
      if (eRank == TheRank - 1) {
        TheSelect.push_back(eRoot);
      }
    }
    SelectedRoots = TheSelect;
  };
  std::cerr << "FundCone, step 6\n";
  for (size_t iRow=0; iRow<nbRow; iRow++)
    if (selected[iRow] == 0) {
      const MyVector<Tint>& uRoot = V1_roots[iRow];
      insert_root(uRoot);
    }
  std::cerr << "FundCone, step 7\n";
  return SelectedRoots;
}


template<typename T, typename Tint>
std::vector<MyVector<Tint>> FindRoots(const VinbergTot<T,Tint>& Vtot)
{
  std::cerr << "FindRoots, step 1\n";
  std::vector<MyVector<Tint>> ListRoot = FundCone(Vtot);
  std::cerr << "FindRoots, step 2\n";

  IterateRootDecompositions<T,Tint> iter(Vtot);
  std::cerr << "FindRoots, step 3\n";
  while (true) {
    const std::pair<MyVector<Tint>,Tint> pair = iter.get_cand();
    const MyVector<Tint>& a = pair.first;
    const Tint& k = pair.second;
    const MyVector<T> a_T = UniversalVectorConversion<T,Tint>(a);
    const T k_T = k;
    std::cerr << "  NextRoot a=" << a << " k=" << k << " k_T=" << k_T << "\n";
    for (const MyVector<Tint>& root_cand : Roots_decomposed_into<T,Tint>(Vtot, a_T, k_T)) {
      if (IsRoot<T,Tint>(Vtot.G, root_cand)) {
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
