#ifndef INCLUDE_VINBERG_ALGO
#define INCLUDE_VINBERG_ALGO

#include "Shvec_exact.h"
//#include "MAT_MatrixInt.h"
#include "POLY_cddlib.h"
#include "POLY_RedundancyElimination.h"
#include "POLY_PolytopeInt.h"
#include "POLY_lrslib.h"
#include "coxeter_dynkin.h"
#include "Temp_ShortVectorUndefinite.h"


#define DEBUG_VINBERG


//
// Finding the closest points
//

// Compute the solutions of G [x - eV] = a
template<typename T, typename Tint, typename Fins>
void ComputeSphericalSolutions(const MyMatrix<T>& GramMat, const MyVector<T>& eV, const T& norm, Fins f_ins)
{
  bool PrintInput=false;
  if (PrintInput) {
    std::cerr << "GramMat=\n";
    WriteMatrix(std::cerr, GramMat);
    std::cerr << "eV=\n";
    WriteVector(std::cerr, eV);
    std::cerr << "norm=" << norm << "\n";
  }
  int mode = TempShvec_globals::TEMP_SHVEC_MODE_VINBERG_ALGO;
  T_shvec_request<T> request = initShvecReq<T>(GramMat, eV, norm, mode);
  //
  auto f_insert=[&](const MyVector<Tint>& V, const T& min) -> bool {
    if (min == norm) {
      f_ins(V);
    }
    return true;
  };
  (void)computeIt<T,Tint,decltype(f_insert)>(request, norm, f_insert);
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
  bool is_even = true;
  size_t n = M.rows();
  for (size_t i=0; i<n; i++) {
    Tint res = M(i,i) % Tint(2);
    if (res != 0)
      is_even = false;
  }
  std::cerr << "is_even=" << is_even << "\n";
  auto is_correct=[&](Tint k) -> bool {
    if (is_even) {
      Tint res2 = k % Tint(2);
      if (res2 != 0)
        return false;
    }
    Tint res = ResInt(limit, k);
    return res == 0;
  };
  //
  Tint k = 1;
  while (true) {
    if (is_correct(k)) {
      root_lengths.push_back(k);
    }
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
  return ListPoint;
}


template<typename Tint>
std::vector<MyVector<Tint>> GetIntegerPoints(const MyMatrix<Tint>& m)
{
  return ComputeTranslationClasses<Tint,Tint>(m);
}



template<typename Tint>
std::vector<MyVector<Tint>> ReduceListRoot(const std::vector<MyVector<Tint>>& ListRoot)
{
  MyMatrix<Tint> M_Tint = MatrixFromVectorFamily(ListRoot);
  using Tfield=typename overlying_field<Tint>::field_type;
  MyMatrix<Tfield> M_Tfield = UniversalMatrixConversion<Tfield,Tint>(M_Tint);
  int ChosenMethod=2;
  auto get_listidx=[&]() -> std::vector<int> {
    if (ChosenMethod == 1) {
      return Kernel_GetNonRedundant_CDD(M_Tfield);
    }
    if (ChosenMethod == 2) {
      MyMatrix<Tfield> M2 = lrs::FirstColumnZero(M_Tfield);
      return cdd::RedundancyReductionClarkson(M2);
    }
    if (ChosenMethod == 3) {
      std::vector<int> ListIdx1 = Kernel_GetNonRedundant_CDD(M_Tfield);
      //      MyMatrix<Tfield> M2 = Polytopization(M_Tfield);
      MyMatrix<Tfield> M2 = lrs::FirstColumnZero(M_Tfield);
      std::vector<int> ListIdx2 = cdd::RedundancyReductionClarkson(M2);
      if (ListIdx1 != ListIdx2) {
        std::cerr << "Inequality between ListIdx1 and ListIdx2\n";
        std::cerr << "ListIdx1=" << ListIdx1 << "\n";
        std::cerr << "ListIdx2=" << ListIdx2 << "\n";
        throw TerminalException{1};
      }
      return ListIdx1;
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  std::vector<int> ListIdx = get_listidx();
  std::vector<MyVector<Tint>> ListV;
  for (auto & idx : ListIdx)
    ListV.push_back(ListRoot[idx]);
  return ListV;
}









template<typename T, typename Tint>
VinbergTot<T,Tint> GetVinbergAux(const MyMatrix<Tint>& G, const MyVector<Tint>& v0)
{
  int n=G.rows();
  // Computing the complement of the space.
  MyMatrix<T> G_T = UniversalMatrixConversion<T,Tint>(G);
  MyVector<Tint> V = G * v0;
  std::vector<Tint> vectV(n);
  for (int i=0; i<n; i++)
    vectV[i] = V(i);
  GCD_int<Tint> eGCDinfo = ComputeGCD_information(vectV);
  std::vector<int> ListZer(n-1);
  for (int j=0; j<n-1; j++)
    ListZer[j] = j + 1;
  MyMatrix<Tint> Morth = SelectColumn(eGCDinfo.Pmat, ListZer);
  MyMatrix<T> Morth_T = UniversalMatrixConversion<T,Tint>(Morth);
  MyMatrix<Tint> M = ConcatenateMatVec_Tr(Morth, V);
  MyMatrix<Tint> M2 = ConcatenateMatVec_Tr(Morth, v0);
  MyMatrix<Tint> M2_tr = M2.transpose();
  std::vector<MyVector<Tint>> W = GetIntegerPoints(M2_tr);
  std::cerr << "W=\n";
  for (auto & eVect : W) {
    WriteVector(std::cerr, eVect);
  }
  // The determinant. The scalar tell us how much we need to the quotient.
  // We will need to consider the vectors k (V_i / eDet) for k=1, 2, 3, ....
  Tint eDet = T_abs(DeterminantMat(M));
  std::cerr << "eDet=" << eDet << "\n";
  // We want to find a vector v such that V_i = (det) v + Morth Z^{n-1}
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
  MyMatrix<Tint> Mbas = ConcatenateMatVec_Tr(Morth, Vtrans);

  // Gram matrix of the space.
  MyMatrix<Tint> Gorth = Morth.transpose() * G * Morth;
  MyMatrix<T> Gorth_T = UniversalMatrixConversion<T,Tint>(Gorth);
  MyMatrix<T> GorthInv = Inverse(Gorth_T);
  // Computing the side comput
  MyMatrix<T> GM_iGorth = G_T * UniversalMatrixConversion<T,Tint>(Morth) * GorthInv;
  std::vector<Tint> root_lengths = Get_root_lengths(G);
  std::cerr << "s.root_lengths =";
  for (auto & eVal : root_lengths)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  //  return {G, G_T, v0, V, Vtrans, Mbas, MbasInv_T, Morth, Morth_T, eDet, Gorth, Gorth_T, GM_iGorth, W, root_lengths};
  return {G, G_T, v0, V, Vtrans, Mbas, Morth, Morth_T, eDet, Gorth, Gorth_T, GM_iGorth, W, root_lengths};
}


template<typename T, typename Tint>
VinbergTot<T,Tint> GetVinbergFromG(const MyMatrix<T>& G)
{
  T CritNorm = 0;
  bool StrictIneq = true;
  bool NeedNonZero = true;
  MyVector<Tint> eVect=GetShortVector_unlimited_float<Tint,T>(G, CritNorm, StrictIneq, NeedNonZero);
  MyMatrix<Tint> G_i = UniversalMatrixConversion<Tint,T>(G);
  return GetVinbergAux<T,Tint>(G_i, eVect);
}




template<typename T, typename Tint>
struct IterateRootDecompositions {
private:
  std::unordered_map<Tint,size_t> candidates;
  const VinbergTot<T,Tint>& Vtot;
  Tint len_sW;
  MyVector<Tint> cand_a(const size_t& n) const {
    size_t len_sW = Vtot.W.size();
    size_t res = n % len_sW;
    size_t q = n / len_sW;
    return Vtot.W[res] + q * Vtot.v0;
  }
  Tint get_k() const {
    bool we_found = false;
    double minval_d=std::numeric_limits<double>::max();
    Tint kfind;
    for (auto & k : Vtot.root_lengths) {
      MyVector<Tint> V2 = cand_a(candidates.at(k));
      Tint val = - Vtot.v0.dot(Vtot.G * V2);
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
    size_t val = candidates[k];
    std::cerr << "IterateRootDecomposition, k=" << k << " val=" << val << "\n";
    candidates[k] = val + 1;
    MyVector<Tint> V = cand_a(val);
    std::cerr << "IterateRootDecomposition, cand_a(candidates[k])=" << V << "\n";
    return {V, k};
  }
};


template<typename T, typename Tint>
struct DataMappingVinbergProblem {
  MyMatrix<T> G;
  MyVector<T> V;
  T norm;
  MyVector<Tint> shift_u;
  MyMatrix<Tint> trans_u;
  bool is_feasible;
};




template<typename T, typename Tint, typename Fins_root>
void Solutioner_CVP(const DataMappingVinbergProblem<T,Tint>& data, Fins_root f_ins_root)
{
  size_t n_sol = 0;
  auto f_ins=[&](const MyVector<Tint>& u) -> void {
    MyVector<Tint> x = data.shift_u + data.trans_u * u;
    n_sol++;
    f_ins_root(x);
  };
  ComputeSphericalSolutions<T,Tint>(data.G, data.V, data.norm, f_ins);
  std::cerr << "|ListSol|=" << n_sol << "\n";
}






/*
  We look for the solutions of (a+v , a+v) = k
  with v in the Morth space.
  (a, a) + 2 (a, v) + (v,v) = k
   v = M w  with  w in Z^{n-1}
  2 a^t G Mw + w^t {M^t G M} w = k - (a,a)
  2 w Gorth sV + w^t Gorth w = k -(a,a)
  (w + sV)^t Gorth (w + sV) = k - (a,a) + sV^t Gorth sV
 */
/*
  We want to find integer vectors such that (a+v, a+v) = k
  with v in the Morth space and one vector a.
  We want to also impose the condition that    2(a+v) G / k    \in Z^n
  Thus the equations that we have are
  ----
  We write v = M w with w \in Z^{n-1}
  2 a^t G Mw + w^t {M^t G M} w = k - (a,a)
  We have the additional condition (2/k) G (a + Mw) \in Z^n
  ----
  We can use SolutionIntmat to get one initial solution w_0
  Then we need also to use a SolutionIntTrMat to get the vectspace U
  So, we have w = w0 + U z
  Putting it all together we get
  2a^T G M (w0 + U z) + (w0^T + z^T U^T) M^T GM (w0 + U z) = k - (a,a)
  z^T {U^T M^T G M U} z + 2 { a^T G M U + w0^T M^T G M U } z = k - (a,a) - w0^T M^T G M w0
  z^T {U^T M^T G M U} z + 2 { U^T M^T G (a + M w0) }^T z = k - (a,a) - w0^T M^T G M w0
  Write Gs = U^T M^T G M U    and Vs = Gs^{-1} Ws and Ws = U^T M^T G (a + M w0)
  And so we finally get the equation
  z^T Gs z + 2 Vs^T Gs z = k - (a,a) - w0^T M^T G M w0
  or
  (z + Vs)^T Gs (z + Vs) = k - (a,a) - w0^T M^T G M w0 + Vs Gs Vs
  ---
  Also we have x = a + Mw
                 = a + M (w0 + Uz)
                 = (a + Mw0) + MU z
  ---
  The equation to find the w0 and U are going to be
  2 G (a + M w) = k h    with (w,h) in Z^{2 n -1}
  So, we get
  2 GM w - k h = -2 G a
 */
template<typename T, typename Tint>
DataMappingVinbergProblem<T,Tint> Get_DataMapping(const VinbergTot<T,Tint>& Vtot, const MyVector<Tint>& a, const Tint& k)
{
  if (k == 1 || k == 2) { // No need for some complex linear algebra work, returning directly
    MyVector<T> a_T = UniversalVectorConversion<T,Tint>(a);
    MyVector<T> sV = a_T.transpose() * Vtot.GM_iGorth;
    Tint term1 = k - a.dot(Vtot.G * a);
    T normi = T(term1) + sV.dot(Vtot.Gorth_T * sV);
    return {Vtot.Gorth_T, sV, normi, a, Vtot.Morth, true};
  }
  //
  size_t n = Vtot.G.rows();
  MyMatrix<Tint> GM = Vtot.G * Vtot.Morth;
  MyVector<Tint> m2_Ga = -2 * Vtot.G * a;
  MyMatrix<Tint> Bmat(2*n-1, n);
  for (size_t i=0; i<n-1; i++) {
    for (size_t j=0; j<n; j++)
      Bmat(i,j) = 2 * GM(j,i);
  }
  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<n; j++) {
      if (i == j)
        Bmat(n-1 + i, i) = -k;
      else
        Bmat(n-1 + i, j) = 0;
    }
  }
  ResultSolutionIntMat<Tint> res = SolutionIntMat(Bmat, m2_Ga);
  if (!res.TheRes)
    return {{}, {}, 0, {}, {}, false};
  //
  MyVector<Tint> w0(n-1);
  for (size_t i=0; i<n-1; i++)
    w0(i) = res.eSol(i);
  MyMatrix<Tint> U_block = NullspaceIntMat(Bmat);
  size_t dim = U_block.rows();
#ifdef DEBUG_VINBERG
  if (dim != n-1) {
    std::cerr << "The dimension dim should be equal to n-1\n";
    throw TerminalException{1};
  }
#endif
  MyMatrix<Tint> U(n-1,dim);
  for (size_t i=0; i<dim; i++)
    for (size_t j=0; j<n-1; j++)
      U(j,i) = U_block(i,j);
  //
  MyMatrix<Tint> Gs = U.transpose() * Vtot.Morth.transpose() * Vtot.G * Vtot.Morth * U;
  MyMatrix<Tint> Ws = U.transpose() * Vtot.Morth.transpose() * Vtot.G * ( a + Vtot.Morth * w0);
  MyMatrix<T> Gs_T = UniversalMatrixConversion<T,Tint>(Gs);
  MyMatrix<T> InvGs_T = Inverse(Gs_T);
  MyVector<T> Vs = InvGs_T * UniversalVectorConversion<T,Tint>(Ws);
  MyVector<Tint> Mw0 = Vtot.Morth * w0;
  Tint term1 = k - a.dot(Vtot.G * a) - Mw0.dot(Vtot.G * Mw0);
  T term2 = Vs.dot(Gs_T * Vs);
  T normi = T(term1) + term2;
  //
  MyVector<Tint> apMw0 = a + Vtot.Morth * w0;
  MyMatrix<Tint> MU = Vtot.Morth * U;
  return {Gs_T, Vs, normi, apMw0, MU, true};
}



template<typename T, typename Tint, typename Fins_root>
void Roots_decomposed_into_select(const VinbergTot<T,Tint>& Vtot, const MyVector<Tint>& a, const Tint& k, Fins_root f_ins_root)
{
  DataMappingVinbergProblem<T,Tint> data = Get_DataMapping(Vtot, a, k);
  if (data.is_feasible)
    Solutioner_CVP(data, f_ins_root);
}



template<typename T, typename Tint>
std::vector<MyVector<Tint>> FindRoot_no_filter(const VinbergTot<T,Tint>& Vtot, const MyVector<Tint>& a, const Tint& k)
{
  std::vector<MyVector<Tint>> list_root;
  auto f_ins_root=[&](const MyVector<Tint>& V) -> void {
    list_root.push_back(V);
  };
  Roots_decomposed_into_select(Vtot, a, k, f_ins_root);
  std::cerr << "|list_root|=" << list_root.size() << "\n";
  return list_root;
}



template<typename T, typename Tint>
std::vector<MyVector<Tint>> FindRoot_filter(const VinbergTot<T,Tint>& Vtot, const MyVector<Tint>& a, const Tint& k, const std::vector<MyVector<Tint>>& ListRoot, const MyMatrix<T>& FACfeasible)
{
  std::vector<MyVector<Tint>> list_root;
  DataMappingVinbergProblem<T,Tint> data = Get_DataMapping(Vtot, a, k);
  if (!data.is_feasible) {
    std::cerr << "Conclude that no solution is feasible\n";
    return {};
  }

  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  auto fct_CVP=[&]() -> void {
    std::vector<MyVector<Tint>> list_GV;
    for (auto & e_root : ListRoot) {
      MyVector<Tint> e_gv = Vtot.G * e_root;
      list_GV.push_back(e_gv);
    }
    auto f_ins_root=[&](const MyVector<Tint>& V) -> void {
      for (auto & e_gv : list_GV) {
        T scal = V.dot(e_gv);
        if (scal > 0)
          return;
      }
      list_root.push_back(V);
    };
    Solutioner_CVP(data, f_ins_root);
  };
  //
  auto fct_CVP_Poly=[&]() -> void {
    size_t n_root = ListRoot.size();
    size_t dim = Vtot.G.rows();
    MyMatrix<T> FAC(n_root,dim);
    for (size_t i_root=0; i_root<n_root; i_root++) {
      MyVector<T> eFAC = GetMatrixRow(FACfeasible,i_root);
      MyVector<T> v_T = UniversalVectorConversion<T,Tint>(data.shift_u);
      T scal = eFAC.dot(v_T);
      FAC(i_root,0) = scal;
      for (size_t i=1; i<dim; i++) {
        v_T = UniversalVectorConversion<T,Tint>(GetMatrixColumn(data.trans_u,i-1));
        T scal = eFAC.dot(v_T);
        FAC(i_root,i) = scal;
      }
    }
    //
    int mode = TempShvec_globals::TEMP_SHVEC_MODE_VINBERG_ALGO;
    T_shvec_request<T> request = initShvecReq<T>(data.G, data.V, data.norm, mode);
    //
    size_t n_pass = 0;
    auto f_insert=[&](const MyVector<Tint>& V, const T& min) -> bool {
      n_pass++;
      if (min == data.norm) {
        MyVector<Tint> x = data.shift_u + data.trans_u * V;
        list_root.emplace_back(std::move(x));
      }
      return true;
    };
    (void)computeIt_polytope<T,Tint,decltype(f_insert)>(request, data.norm, FAC, f_insert);
    std::cerr << "n_pass=" << n_pass << " |list_root|=" << list_root.size() << "\n";
  };
  //
  /*
    This code does not work because we can have an unbounded polytope
    and for this kind of polytope the enumeration algorithm will not work.
    --------------
  auto fct_PolyInt=[&]() -> void {
    size_t n_root = ListRoot.size();
    size_t dim = Vtot.G.rows();
    MyMatrix<T> FAC(n_root, dim);
    for (size_t i_root=0; i_root<n_root; i_root++) {
      MyVector<Tint> e_gv = - Vtot.G * ListRoot[i_root];
      Tint scal1 = e_gv.dot(data.shift_u);
      FAC(i_root,0) = UniversalScalarConversion<T,Tint>(scal1);
      for (size_t i=1; i<dim; i++) {
        MyVector<Tint> eCol = GetMatrixColumn(data.trans_u, i-1);
        Tint scal2 = e_gv.dot(eCol);
        FAC(i_root,i) = UniversalScalarConversion<T,Tint>(scal2);
      }
    }
    //    std::cerr << "FAC=\n";
    //    WriteMatrix(std::cerr, FAC);
    std::cerr << "fct_PolyInt m step 1\n";
    MyMatrix<T> EXT = cdd::DualDescription(FAC);
    std::cerr << "EXT=\n";
    WriteMatrix(std::cerr, EXT);
    std::cerr << "fct_PolyInt m step 2\n";
    T max_norm = 0;
    size_t n_vert = EXT.rows();
    std::cerr << "n_vert=" << n_vert << "\n";
    MyVector<T> shift_u_T = UniversalVectorConversion<T,Tint>(data.shift_u);
    MyMatrix<T> trans_u_T = UniversalMatrixConversion<T,Tint>(data.trans_u);
    std::cerr << "fct_PolyInt m step 3\n";
    int n_col = trans_u_T.cols();
    MyVector<T> V(n_col);
    size_t n_greater = 0;
    T k_T = k;
    for (size_t i_vert=0; i_vert<n_vert; i_vert++) {
      T val0 = EXT(i_vert,0);
      if (val0 > 0) {
        for (int i=0; i<n_col; i++)
          V(i) = EXT(i_vert,i + 1) / val0;
        MyVector<T> x_T = shift_u_T + trans_u_T * V;
        T scal_T = x_T.dot(Vtot.G_T * x_T);
        std::cerr << "i_vert=" << i_vert << "/" << n_vert << " scal=" << scal_T << "\n";
        if (scal_T > k_T)
          n_greater++;
        if (scal_T > max_norm)
          max_norm = scal_T;
      } else {
        for (int i=0; i<n_col; i++)
          V(i) = EXT(i_vert,i + 1);
        std::cerr << "Finding an infinite ray V=";
        WriteVector(std::cerr, V);
        Tint prod = 1;
        for (int j=0; j<10; j++) {
          MyVector<T> x_T = shift_u_T + (prod - 1) * trans_u_T * V;
          prod *= 2;
          T scal_T = x_T.dot(Vtot.G_T * x_T);
          std::cerr << "Finding an infinite ray x=";
          WriteVector(std::cerr, x_T);
          std::cerr << "scal_T=" << scal_T << "\n";
        }
        max_norm = k + 1;
      }
    }
    std::cerr << "k=" << k << " max_norm=" << max_norm << " n_greater=" << n_greater << "\n";
    if (max_norm < k_T)
      return;
    std::cerr << "fct_PolyInt m step 4\n";
    size_t n_interior = 0;
    auto f_insert=[&](const MyVector<Tint>& V) -> bool {
      MyVector<Tint> x = data.shift_u + data.trans_u * V;
      Tint scal = x.dot(Vtot.G * x);
      if (scal == k) {
        list_root.push_back(x);
      }
      n_interior++;
      return true;
    };
    //    Kernel_GetListIntegralPoint<T,Tint,decltype(f_insert)>(FAC, EXT, f_insert);
    Kernel_GetListIntegralPoint_LP<T,Tint,decltype(f_insert)>(FAC, f_insert);
    std::cerr << "n_interior=" << n_interior << " |list_root|=" << list_root.size() << "\n";
  };
  */
  //
  MyMatrix<Tint> RootMat = MatrixFromVectorFamily(ListRoot);
  if (RankMat(RootMat) == Vtot.G.rows()) {
    fct_CVP_Poly();
  } else {
    fct_CVP();
  }
  //
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|list_root|=" << list_root.size() << " |FindRoot_filter|=" << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
  return list_root;
}




template<typename T, typename Tint>
MyMatrix<int> GetWeightMatrix(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
{
  size_t n_root = ListRoot.size();
  size_t nbCol = Vtot.G.rows();
  std::cerr << "GetWeightMatrix : begin, n_root=" << n_root << " nbCol=" << nbCol << "\n";
  MyMatrix<T> M(n_root, n_root);
  std::unordered_map<T,int> DiagVal;
  std::unordered_map<T,int> OffDiagVal;
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<Tint> eVG = Vtot.G * ListRoot[i_root];
    for (size_t j_root=0; j_root<n_root; j_root++) {
      T eScal = eVG.dot(ListRoot[j_root]);
      M(i_root, j_root) = eScal;
      if (i_root == j_root) {
        DiagVal[eScal] += 1;
      } else {
        if (i_root < j_root)
          OffDiagVal[eScal] += 1;
      }
    }
  }
  std::cerr << "GetWeightMatrix : Diag =";
  for (auto & kv : DiagVal)
    std::cerr << " [" << kv.first << "," << kv.second << "]";
  std::cerr << "\n";
  std::cerr << "GetWeightMatrix : OffDiag =";
  for (auto & kv : OffDiagVal)
    std::cerr << " [" << kv.first << "," << kv.second << "]";
  std::cerr << "\n";

  std::unordered_map<T,int> CosVal;
  MyMatrix<T> Cos(n_root,n_root);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    for (size_t j_root=0; j_root<n_root; j_root++) {
      T aII = M(i_root,i_root);
      T aJJ = M(j_root,j_root);
      T aIJ = M(i_root,j_root);
      T cos2 = (aIJ * aIJ) / (aII * aJJ);
      Cos(i_root,j_root) = cos2;
      if (i_root < j_root)
        CosVal[cos2] += 1;
    }
  }
  std::cerr << "GetWeightMatrix : Cos =";
  for (auto & kv : CosVal)
    std::cerr << " [" << kv.first << "," << kv.second << "]";
  std::cerr << "\n";
  //
  // Building of the input
  //
  T cst1 = T(1) / T(4); //  1/4
  T cst2 = T(1) / T(2); //  1/2
  T cst3 = T(3) / T(4); //  3/4
  T cst4 = 1; //  1
  auto weight=[&](int i, int j) -> int {
    T aII = M(i,i);
    T aJJ = M(j,j);
    T aIJ = M(i,j);
    T cos2 = (aIJ * aIJ) / (aII * aJJ);
    if (cos2 == 0)
      return 2;
    if (cos2 == cst1)
      return 3;
    if (cos2 == cst2)
      return 4;
    if (cos2 == cst3)
      return 6;
    if (cos2 == cst4)
      return 0;
    if (cos2 > cst4)
      return 1;
    std::cerr << "coxiter.py ERROR: cosine " << cos2 << "\n";
    throw TerminalException{1};
  };
  MyMatrix<int> WeightMatrix(n_root,n_root);
  for (size_t i=0; i<n_root; i++)
    for (size_t j=0; j<n_root; j++)
      WeightMatrix(i,j) = 444;
  //
  for (size_t i=0; i<n_root; i++)
    for (size_t j=0; j<i; j++)
      if (M(i,j) != 0) {
        int w = weight(i, j);
        WeightMatrix(i,j) = w;
        WeightMatrix(j,i) = w;
      }
  return WeightMatrix;
}





/*

template<typename T, typename Tint>
bool is_FundPoly(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
{
  MyMatrix<int> WeightMatrix = GetWeightMatrix(Vtot, ListRoot);
  size_t n_root = ListRoot.size();
  int d = Vtot.G.rows() - 1;
  //
  size_t nb_spherical_subdiagramN_parabolic_subdiagramN1 = 0;
  std::vector<int> ListPos;
  ListPos.reserve(n_root);
  size_t pos=0;
  auto is_extendible_to_spherical=[&](const std::vector<int>& V) -> bool {
    
  };
  auto is_spherical
    
  
  while(true) {
  }
}

*/



template<typename T, typename Tint>
std::optional<MyVector<Tint>> GetOneInteriorVertex(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  size_t n_root = ListRoot.size();
  size_t n_col = Vtot.G.rows();
  MyMatrix<Tint> FAC(n_root,n_col);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<Tint> e_gv = - Vtot.G * ListRoot[i_root];
    for (size_t i_col=0; i_col<n_col; i_col++)
      FAC(i_root, i_col) = e_gv(i_col);
  }
  MyMatrix<Tint> FACwork=lrs::FirstColumnZero(FAC);
  bool IsFirst = true;
  std::optional<MyVector<Tint>> opt;
  size_t n_iter = 0;
  auto f=[&](Tint* out) -> bool {
    if (!IsFirst) {
      n_iter++;
      MyVector<Tint> V(n_col);
      for (size_t i_col=0; i_col<n_col; i_col++)
        V(i_col) = out[i_col+1];
      Tint scal = V.dot(Vtot.G * V);
      if (scal <= 0) {
        opt = V;
        return false;
      }
    }
    IsFirst=false;
    return true;
  };
  lrs::Kernel_DualDescription_cond(FACwork, f);
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  bool test=opt.has_value();
  std::cerr << "has found vertex=" << test << " n_iter=" << n_iter << " |GetOneInteriorVertex|=" << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
  return opt;
}




template<typename T, typename Tint>
bool is_FundPoly_LRS(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
{
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  size_t n_root = ListRoot.size();
  size_t n_col = Vtot.G.rows();
  MyMatrix<Tint> FAC(n_root,n_col);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<Tint> e_gv = - Vtot.G * ListRoot[i_root];
    for (size_t i_col=0; i_col<n_col; i_col++)
      FAC(i_root, i_col) = e_gv(i_col);
  }
  MyMatrix<Tint> FACwork=lrs::FirstColumnZero(FAC);
  bool IsFiniteCovolume = true;
  bool IsFirst = true;
  size_t n_iter = 0;
  auto f=[&](Tint* out) -> bool {
    if (!IsFirst) {
      n_iter++;
      MyVector<Tint> V(n_col);
      for (size_t i_col=0; i_col<n_col; i_col++)
        V(i_col) = out[i_col+1];
      Tint scal = V.dot(Vtot.G * V);
      if (scal > 0) {
        IsFiniteCovolume = false;
        return false;
      }
    }
    IsFirst=false;
    return true;
  };
  lrs::Kernel_DualDescription_cond(FACwork, f);
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "IsFiniteCovolume=" << IsFiniteCovolume << " n_iter=" << n_iter << " |is_FundPoly_LRS|=" << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
  return IsFiniteCovolume;
}






template<typename T, typename Tint>
bool is_FundPoly_Coxiter(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
{
  MyMatrix<int> WeightMatrix = GetWeightMatrix(Vtot, ListRoot);
  size_t n_root = ListRoot.size();
  int d = Vtot.G.rows() - 1;
  std::string rnd_str = random_string(20);
  std::string FileI = "/tmp/CoxIter_" + rnd_str + ".input";
  std::string FileO = "/tmp/CoxIter_" + rnd_str + ".out";
  {
    std::ofstream os(FileI);
    os << n_root << " " << d << "\n";
    for (size_t i=0; i<n_root; i++)
      for (size_t j=0; j<i; j++)
        if (WeightMatrix(i,j) != 444)
          os << (j+1) << " " << (i+1) << " " << WeightMatrix(i, j) << "\n";
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
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
  int iret=system(eCommand.c_str());
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  if (iret == -1) {
    printf("Oh dear, something went wrong with coxiter! %s\n", strerror(errno));
    throw TerminalException{1};
  }

  //
  // Reading the output
  //
  bool IsFiniteCovolume=false;
  std::string line, answer = "yes", question = "Finite covolume: ";
  std::ifstream INfs(FileO);
  while (getline(INfs, line)) {
    std::vector<std::string> LStr1 = STRING_Split(line, question);
    if (LStr1.size() > 1)
      if (LStr1[1] == answer)
        IsFiniteCovolume = true;
  }
  std::cerr << "is_FundPoly IsFiniteCovolume=" << IsFiniteCovolume << "  |coxiter|=" << std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count() << "\n";
  return IsFiniteCovolume;
}



template<typename T, typename Tint>
bool is_FundPoly(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
{
  int ChosenMethod = 2;
  if (ChosenMethod == 1)
    return is_FundPoly_Coxiter(Vtot, ListRoot);
  if (ChosenMethod == 2)
    return is_FundPoly_LRS(Vtot, ListRoot);
  //
  std::cerr << "Failed to find a matching entry\n";
  throw TerminalException{1};
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
struct DataReflectionGroup {
  std::vector<MyVector<Tint>> ListRoot;
  MyMatrix<Tint> G;
  MyMatrix<Tint> M;
  MyMatrix<T> Cos;
};


template<typename T, typename Tint>
DataReflectionGroup<T,Tint> GetDataReflectionGroup(const std::vector<MyVector<Tint>>& ListRoot, const MyMatrix<Tint>& G)
{
  size_t n_root = ListRoot.size();
  MyMatrix<Tint> M(n_root,n_root);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<Tint> eVG = G * ListRoot[i_root];
    for (size_t j_root=0; j_root<n_root; j_root++) {
      T eScal = eVG.dot(ListRoot[j_root]);
      M(i_root, j_root) = eScal;
    }
  }
  //
  MyMatrix<T> Cos(n_root,n_root);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    for (size_t j_root=0; j_root<n_root; j_root++) {
      T aII = M(i_root,i_root);
      T aJJ = M(j_root,j_root);
      T aIJ = M(i_root,j_root);
      T cos2 = (aIJ * aIJ) / (aII * aJJ);
      Cos(i_root,j_root) = cos2;
    }
  }
  return {ListRoot, G, std::move(M), std::move(Cos)};
}



template<typename T, typename Tint>
void Print_DataReflectionGroup_TXT(const DataReflectionGroup<T,Tint>& data, std::ostream& os)
{
  size_t n_root = data.ListRoot.size();
  os << "|ListRoot|=" << n_root << "\n";
  for (size_t i=0; i<n_root; i++)
    WriteVector(os, data.ListRoot[i]);
  os << "G=\n";
  WriteMatrix(os, data.G);
  os << "M=\n";
  WriteMatrix(os, data.M);
  os << "Cos=\n";
  WriteMatrix(os, data.Cos);
}

template<typename T, typename Tint>
void Print_DataReflectionGroup_GAP(const DataReflectionGroup<T,Tint>& data, std::ostream& os)
{
  os << "return rec(ListRoot:=";
  MyMatrix<Tint> MatRoot = MatrixFromVectorFamily(data.ListRoot);
  WriteMatrixGAP(os, MatRoot);
  os << ", G:=";
  WriteMatrixGAP(os, data.G);
  os << ", M:=";
  WriteMatrixGAP(os, data.M);
  os << ", Cos:=";
  WriteMatrixGAP(os, data.Cos);
  os << ");\n";
}




template<typename T, typename Tint>
void Print_DataReflectionGroup(const DataReflectionGroup<T,Tint>& data, const std::string& mode, std::ostream& os)
{
  if (mode == "txt")
    return Print_DataReflectionGroup_TXT(data, os);
  if (mode == "gap")
    return Print_DataReflectionGroup_GAP(data, os);
  std::cerr << "Failed to find matching entry\n";
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
  MyVector<Tint> a = ZeroVector<Tint>(n);
  std::cerr << "FundCone, step 1.1\n";
  for (auto & k : Vtot.root_lengths) {
    std::cerr << " k=" << k << "\n";
    std::set<MyVector<Tint>> set;
    std::vector<MyVector<Tint>> list_root_cand = FindRoot_no_filter<T,Tint>(Vtot, a, k);
    for (const MyVector<Tint>& root_cand : list_root_cand) {
      MyVector<Tint> root_can = SignCanonicalizeVector(root_cand);
      set.insert(root_can);
    }
    for (auto & eV : set)
      V1_roots.push_back(eV);
    std::cerr << "k=" << k << " |set|=" << set.size() << " |V1_roots|=" << V1_roots.size() << " |list_root_cand|=" << list_root_cand.size() << "\n";
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
  std::cerr << "SelectedRoots=\n";
  for (auto & eVect : SelectedRoots)
    WriteVector(std::cerr, eVect);
  return SelectedRoots;
}


// We compute here a polyhedral domain in which we expect to find
// possibly roots.
//
// The computed domain lies inside of the cone and it comes before the
// business of Vinberg of selecting the level at which the roots will
// be found.
//
// All tricks are allowed in that search. We only want to restrict the
// search space
template<typename T, typename Tint>
MyMatrix<T> GetInitial_FACfeasible(const VinbergTot<T,Tint>& Vtot, const std::vector<MyVector<Tint>>& ListRoot)
{
  size_t n_root = ListRoot.size();
  size_t dim = Vtot.G.rows();
  MyMatrix<T> FACfeasible(n_root, dim);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<Tint> e_gv = - Vtot.G * ListRoot[i_root];
    for (size_t i=0; i<dim; i++) {
      T val = UniversalScalarConversion<T,Tint>(e_gv[i]);
      FACfeasible(i_root,i) = val;
    }
  }
  return FACfeasible;
}




template<typename T, typename Tint, typename F>
void FindRoots_Kernel(const VinbergTot<T,Tint>& Vtot, F f_exit)
{
  std::cerr << "FindRoots, step 1\n";
  std::vector<MyVector<Tint>> ListRoot = FundCone(Vtot);
  MyMatrix<T> FACfeasible = GetInitial_FACfeasible(Vtot, ListRoot);
  std::cerr << "FindRoots, step 2\n";

  IterateRootDecompositions<T,Tint> iter(Vtot);
  std::cerr << "FindRoots, step 3\n";
  while (true) {
    if (f_exit(ListRoot, FACfeasible))
      break;
    const std::pair<MyVector<Tint>,Tint> pair = iter.get_cand();
    const MyVector<Tint>& a = pair.first;
    const Tint& k = pair.second;
    std::vector<MyVector<Tint>> list_root_cand = FindRoot_filter<T,Tint>(Vtot, a, k, ListRoot, FACfeasible);
    if (list_root_cand.size() > 0) {
      for (auto & eRoot : list_root_cand) {
        ListRoot.push_back(eRoot);
      }
      std::cerr << "ListRoot=\n";
      for (auto & eRoot : ListRoot) {
        WriteVector(std::cerr, eRoot);
      }
      std::cerr << "After insert |ListRoot|=" << ListRoot.size() << "\n";
      ListRoot = ReduceListRoot(ListRoot);
      std::cerr << "After ReduceListRoot |ListRoot|=" << ListRoot.size() << "\n";
      FACfeasible = GetInitial_FACfeasible(Vtot, ListRoot);
    }
  }
}



template<typename T, typename Tint>
std::vector<MyVector<Tint>> FindRoots(const VinbergTot<T,Tint>& Vtot)
{
  std::vector<MyVector<Tint>> ListRootRet;
  auto f_exit=[&](std::vector<MyVector<Tint>> const& ListRoot, MyMatrix<T> const& FACfeasible) -> bool {
    if (is_FundPoly(Vtot, ListRoot)) {
      ListRootRet = ListRoot;
      return true;
    }
    return false;
  };
  FindRoots_Kernel(Vtot, f_exit);
  return ListRootRet;
}




template<typename T, typename Tint>
std::pair<MyVector<Tint>, std::vector<MyVector<Tint>>> FindOneInitialRay(const VinbergTot<T,Tint>& Vtot)
{
  std::pair<MyVector<Tint>, std::vector<MyVector<Tint>>> epair;
  auto f_exit=[&](std::vector<MyVector<Tint>> const& ListRoot, MyMatrix<T> const& FACfeasible) -> bool {
    std::optional<MyVector<Tint>> opt = GetOneInteriorVertex(Vtot, ListRoot);
    if (opt) {
      MyVector<Tint> v = *opt;
      std::vector<MyVector<Tint>> l_roots;
      for (auto & root : ListRoot) {
        Tint scal = v.dot(Vtot.G * root);
        if (scal == 0)
          l_roots.push_back(root);
      }
      epair = {v, l_roots};
      return true;
    }
    return false;
  };
  FindRoots_Kernel(Vtot, f_exit);
  return epair;
}









#endif
