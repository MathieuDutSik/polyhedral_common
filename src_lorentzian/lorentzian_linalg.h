#ifndef DEFINE_LORENTZIAN_LINALG_H
#define DEFINE_LORENTZIAN_LINALG_H


#include "POLY_cddlib.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "COMB_Combinatorics.h"

/*
  A few linear algebra stuff used for the lorentzian computations
 */


/*
  Given a lattice L and a matrix g, find the smallest exponent m such that g^m preserves L
 */
template<typename T>
size_t GetMatrixExponentSublattice(MyMatrix<T> const& g, MyMatrix<T> const& Latt)
{
  int n = Latt.rows();
  auto is_preserving=[&](MyMatrix<T> const& h) -> bool {
    for (int i=0; i<n; i++) {
      MyVector<T> eV = GetMatrixRow(Latt, i);
      MyVector<T> eVimg = h.transpose() * eV;
      std::optional<MyVector<T>> opt = SolutionIntMat(Latt, eVimg);
      if (!opt)
        return false;
    }
    return true;
  };
  size_t ord = 1;
  MyMatrix<T> h = g;
  while(true) {
    if (is_preserving(h))
      break;
    ord++;
    h = h * g;
  }
  return ord;
}




/*
  Given a lattice L and a matrix g, find the smallest exponent m such that g^m preserves L
  and acts trivially on the set of translation classes.
 */
template<typename T>
size_t GetMatrixExponentSublattice_TrivClass(MyMatrix<T> const& g, MyMatrix<T> const& Latt)
{
  // First compute the power that preserves the lattice L
  int n = Latt.rows();
  auto is_preserving=[&](MyMatrix<T> const& h) -> bool {
    for (int i=0; i<n; i++) {
      MyVector<T> eV = GetMatrixRow(Latt, i);
      MyVector<T> eVimg = h.transpose() * eV;
      std::optional<MyVector<T>> opt = SolutionIntMat(Latt, eVimg);
      if (!opt)
        return false;
    }
    return true;
  };
  size_t ord1 = 1;
  MyMatrix<T> h = g;
  while(true) {
    if (is_preserving(h))
      break;
    ord1++;
    h = h * g;
  }
  // Now computing the power of the action on he classes
  std::vector<MyVector<T>> ListTrans = ComputeTranslationClasses<T,T>(Latt);
  size_t n_class = ListTrans.size();
  std::vector<size_t> cl_h = GetActionOnClasses(ListTrans, h, Latt);
  auto is_identity=[&](std::vector<size_t> const& cl_k) -> bool {
    for (size_t i=0; i<n_class; i++)
      if (cl_k[i] != i)
        return false;
    return true;
  };
  std::vector<size_t> cl_pow = cl_h;
  size_t ord2 = 1;
  while(true) {
    if (is_identity(cl_pow))
      break;
    ord2++;
    std::vector<size_t> W(n_class);
    for (size_t i=0; i<n_class; i++)
      W[i] = cl_pow[cl_h[i]];
    cl_pow = W;
  }
  // Now combining the info
  size_t ord = ord1 * ord2;
  return ord;
}




/*
  We are looking forthe smallest solution c>0 for the equation
  u + c k in Latt
  By selecting an adequate basis for Latt we can reduce the problem to
  u + ck = a1 v1 + a2 v2
  with a1, a2 in Z and v1, v2 in Latt.
  We write u = u1 v1 + x2 u2   and   k = k1 v1 + k2 v2
  This gets us
  u1 + ck1 = a1
  u2 + ck2 = a2
  The right way to solve the equation is to compute kG = gcd(k1, k2) and a basis of the kernel.
  We thus remap the equation to
  u1 + c k1 = a1
  u2        = a2
  solvability condition becomes u2 in Z.
  c0 = -u1 / k1
  cS = 1/k1
  Solution is c = c0 + h cS
  k = - c0 / cS
 */
template<typename T>
std::optional<MyVector<T>> ResolveLattEquation(MyMatrix<T> const& Latt, MyVector<T> const& u, MyVector<T> const& k)
{
  //  std::cerr << "ResolveLattEquation k="; WriteVector(std::cerr, RemoveFractionVector(k));
  std::vector<MyVector<T>> l_v = {u,k};
  MyMatrix<T> eIndep = MatrixFromVectorFamily(l_v);
  MyMatrix<T> IntBasis = IntersectionLattice_VectorSpace(Latt, eIndep);
  //  std::cerr << "IntBasis=\n";
  //  WriteMatrix(std::cerr, IntBasis);
  std::optional<MyVector<T>> opt_u = SolutionMat(IntBasis, u);
  if (!opt_u) {
    std::cerr << "We failed to find a solution for u\n";
    throw TerminalException{1};
  }
  const MyVector<T>& sol_u = *opt_u;
  T u1 = sol_u(0);
  T u2 = sol_u(1);
  //  std::cerr << "u1=" << u1 << " u2=" << u2 << "\n";
  std::optional<MyVector<T>> opt_k = SolutionMat(IntBasis, k);
  if (!opt_k) {
    std::cerr << "We failed to find a solution for k\n";
    throw TerminalException{1};
  }
  const MyVector<T>& sol_k = *opt_k;
  T k1 = sol_k(0);
  T k2 = sol_k(1);
  //  std::cerr << "k1=" << k1 << " k2=" << k2 << "\n";
  //  std::cerr << "u=";
  //  WriteVector(std::cerr, u);
  //  std::cerr << "k=";
  //  WriteVector(std::cerr, k);
  //
  GCD_int<T> ep = ComputePairGcd(k1, k2);
  T u1_norm = ep.Pmat(0,0) * u1 + ep.Pmat(1,0) * u2;
  T u2_norm = ep.Pmat(0,1) * u1 + ep.Pmat(1,1) * u2;
  T k1_norm = ep.Pmat(0,0) * k1 + ep.Pmat(1,0) * k2;
  T k2_norm = ep.Pmat(0,1) * k1 + ep.Pmat(1,1) * k2;
  //  std::cerr << "norm : u1=" << u1_norm << " u2=" << u2_norm << "\n";
  //  std::cerr << "norm : k1=" << k1_norm << " k2=" << k2_norm << "\n";
  if (k2_norm != 0) {
    std::cerr << "We should have k2_norm = 0. Likely a bug here\n";
    throw TerminalException{1};
  }
  if (!IsInteger(u2_norm)) // No solution then
    return {};
  //
  T c0 = - u1_norm / k1_norm;
  T cS = 1 / k1_norm;
  //  std::cerr << "c0=" << c0 << " cS=" << cS << "\n";
  T hinp = - c0 / cS;
  T h;
  if (cS > 0) {
    h = UniversalCeilScalarInteger<T,T>(hinp);
    //    std::cerr << "1 : hinp=" << hinp << " h=" << h << "\n";
    if (hinp == h)
      h += 1;
  } else {
    h = UniversalFloorScalarInteger<T,T>(hinp);
    //    std::cerr << "2 : hinp=" << hinp << " h=" << h << "\n";
    if (hinp == h)
      h -= 1;
  }
  T c = c0 + h * cS;
  //  std::cerr << "h=" << h << " c=" << c << "\n";
  if (c <= 0) {
    std::cerr << "We should have c>0\n";
    throw TerminalException{1};
  }
  return u + c * k;
}




template<typename T, typename Tint>
MyMatrix<T> Get_Pplane(MyMatrix<T> const& G, std::vector<MyVector<Tint>> const& l_ui)
{
  int n = G.rows();
  size_t n_root = l_ui.size();
  MyMatrix<T> EquaPplane(n_root,n);
  for (size_t i_root=0; i_root<n_root; i_root++) {
    MyVector<T> eV = UniversalVectorConversion<T,Tint>(l_ui[i_root]);
    MyVector<T> eP = G * eV;
    AssignMatrixRow(EquaPplane, i_root, eP);
  }
  MyMatrix<T> Pplane = NullspaceTrMat(EquaPplane);
  if (Pplane.rows() != 2) {
    std::cerr << "The dimension should be exactly 2\n";
    std::cerr << "We have |Pplane|=" << Pplane.rows() << "\n";
    throw TerminalException{1};
  }
  return Pplane;
}


template<typename T>
struct LatticeProjectionFramework {
  MyMatrix<T> ProjP;
  MyMatrix<T> BasisProj;
  MyMatrix<T> ProjFamily;
  MyMatrix<T> Latt;
  LatticeProjectionFramework(MyMatrix<T> const& G, MyMatrix<T> const& Subspace, MyMatrix<T> const& _Latt) : Latt(_Latt)
  {
    //    std::cerr << "LatticeProjectionFramework, step 1\n";
    int n = G.rows();
    int dim = Latt.rows();
    //    std::cerr << "n=" << n << " dim=" << dim << "\n";
    //    std::cerr << "|Subspace|=" << Subspace.rows() << " / " << Subspace.cols() << "\n";
    //    std::cerr << "LatticeProjectionFramework, step 2\n";
    ProjP = GetProjectionMatrix(G, Subspace);
    //    std::cerr << "LatticeProjectionFramework, step 3\n";
    ProjFamily = MyMatrix<T>(dim,n);
    //    std::cerr << "LatticeProjectionFramework, step 4\n";
    for (int i=0; i<dim; i++) {
      MyVector<T> eVect = GetMatrixRow(Latt, i);
      MyVector<T> eVectProj = ProjP * eVect;
      AssignMatrixRow(ProjFamily, i, eVectProj);
    }
    //    std::cerr << "LatticeProjectionFramework, step 5\n";
    BasisProj = GetZbasis(ProjFamily);
    //    std::cerr << "LatticeProjectionFramework, step 6\n";
  }
  std::optional<MyVector<T>> GetOnePreimage(MyVector<T> const& V) const
  {
    std::optional<MyVector<T>> opt = SolutionIntMat(ProjFamily, V);
    if (!opt)
      return {};
    MyVector<T> const& eSol = *opt;
    MyVector<T> preImage = Latt.transpose() * eSol;
    return preImage;
  }
};



template<typename T>
std::vector<size_t> GetFacetOneDomain_ListIdx(std::vector<MyVector<T>> const& l_vect)
{
  using Tfield = typename overlying_field<T>::field_type;
  std::cerr << "|l_vect|=" << l_vect.size() << "\n";
  int dimSpace = l_vect[0].size();
  std::cerr << "dimSpace=" << dimSpace << "\n";
  if (l_vect.size() < size_t(2*dimSpace)) {
    std::cerr << "Number of roots should be at least 2 * dimspace = " << (2 * dimSpace) << "\n";
    std::cerr << "while |l_vect|=" << l_vect.size() << "\n";
    throw TerminalException{1};
  }
  auto is_corr=[&](MyVector<T> const& w) -> bool {
    for (auto & e_root : l_vect) {
      T scal = e_root.dot(w);
      if (scal == 0)
        return false;
    }
    return true;
  };
  auto get_random_vect=[&]() -> MyVector<T> {
    MyVector<T> w(dimSpace);
    int spr = 1000;
    int tot_spr = 2 * spr + 1;
    while (true) {
      for (int i=0; i<dimSpace; i++)
        w(i) = rand() % tot_spr - spr;
      std::cerr << "get_random_vect. Trying w=" << StringVectorGAP(w) << "\n";
      if (is_corr(w))
        return w;
    }
  };
  MyVector<T> selVect = get_random_vect();
  std::cerr << "Random splitting vector selVect=" << StringVectorGAP(selVect) << "\n";
  int n_vect = l_vect.size() / 2;
  MyMatrix<Tfield> EXT(n_vect,1+dimSpace);
  std::vector<size_t> list_idx(n_vect);
  size_t pos=0;
  for (size_t i=0; i<l_vect.size(); i++) {
    T scal = selVect.dot(l_vect[i]);
    if (scal > 0) {
      list_idx[pos] = i;
      MyVector<T> const& eV = l_vect[i];
      EXT(pos,0);
      for (int i=0; i<dimSpace; i++)
        EXT(pos,i+1) = UniversalScalarConversion<Tfield,T>(eV(i));
      pos++;
    }
  }
  std::vector<int> list_red = cdd::RedundancyReductionClarkson(EXT);
  size_t siz = list_red.size();
  std::vector<size_t> l_idx(siz);
  for (size_t i=0; i<siz; i++) {
    size_t pos = list_idx[list_red[i]];
    l_idx[i] = pos;
  }
  return l_idx;
}


template<typename T>
std::vector<MyVector<T>> GetFacetOneDomain(std::vector<MyVector<T>> const& l_vect)
{
  size_t n_vect = l_vect.size();
  MyMatrix<T> Mvect = MatrixFromVectorFamily(l_vect);
  MyMatrix<T> MvectRed = ColumnReduction(Mvect);
  std::vector<MyVector<T>> l_vect_red(n_vect);
  for (size_t i=0; i<n_vect; i++)
    l_vect_red[i] = GetMatrixRow(MvectRed, i);
  std::vector<size_t> l_idx = GetFacetOneDomain_ListIdx(l_vect_red);
  std::vector<MyVector<T>> l_vect_ret;
  for (auto & idx : l_idx)
    l_vect_ret.push_back(l_vect[idx]);
  return l_vect_ret;
}


/*
  We should compute the image of Subspace1 into Subspace2 using the eEquiv.
  Afterwards, we can assume that eEquiv = Id.
  We select another vector v1 in the complement of Subspace1, such that {Subspace1, v1}
  is a Z-basis on Z^n.
  We hus want to find the vector v2 such that G2[v2] = G1[v1] and V_{1i} G1 v1 = V_{2i} G2 v2
  The linear equation V_{1i} G1 v1 = V_{2i} G2 v2 has a solution set v2 = v2_0 + alpha w2_0
  because G2 is non-degenerate.
  Then we have the quadratic equation G2[v2] and two possible solutions.
  ----
  The question is that it is possible that there are 3 possibilities:
  (a) No Solution
  (b) 1 solution
  (c) two solutions
  We can have integer solutions or rational solutions or maybe irrational solutions.
  Are irrational solutions possible? I would think that this is not possible.
  ----
  If we have a rational solutions, then we can use the usual scheme
  for passing from rational to integral.
  ----
  Actually things are quite beautiful in the case of Subspace1 the orthogonal of an
  isotropic case.
  We select a vector eVect1 outside of Subspace1 and look for its image eVect2
  We thus have G1[eVect1] = G2[eVect2]
  and the equalities V1(i) * G1 * eVect1 = V2(i) * G2 * eVect2
  This gets us a solution space of eVect2 = V0 + t V1
  and the vector V1 is isotropic.
  This gets us
  G1[eVect1] = G2[eVect2] = G2[V0] + 2t V0.dot.V1
  and thus the equation system has an unique solution, which all turn out to be unexpected.
 */
template<typename T>
MyMatrix<T> ExtendOrthogonalIsotropicIsomorphism_Basis(MyMatrix<T> const& G1, MyMatrix<T> const& Subspace1, MyMatrix<T> const& G2, MyMatrix<T> const& Subspace2)
{
  int dim = G1.rows();
  if (Subspace1.rows() != dim-1 || Subspace2.rows() != dim-1) {
    std::cerr << "Subspace1 and Subspace2 are not of the right dimension\n";
    throw TerminalException{1};
  }
  MyMatrix<T> Compl1 = SubspaceCompletionRational(Subspace1, dim);
  if (Compl1.rows() != 1) {
    std::cerr << "Compl1 should be of dimension 1\n";
    throw TerminalException{1};
  }
  MyVector<T> eVect1 = GetMatrixRow(Compl1,0);
  T eNorm = eVect1.dot(G1 * eVect1);
  MyMatrix<T> eProd1 = Subspace1 * G1;
  MyMatrix<T> eProd2 = Subspace2 * G2;
  MyVector<T> Vscal = Subspace1 * G1 * eVect1;
  std::optional<MyVector<T>> opt = SolutionMat(TransposedMat(eProd2), Vscal);
  if (!opt) {
    std::cerr << "ExtendOrthogonalIsotropicIsomorphism : The solutioning failed\n";
    throw TerminalException{1};
  }
  MyVector<T> const& V0 = *opt;
  MyMatrix<T> NSP = NullspaceTrMat(eProd2);
  if (NSP.rows() != 1) {
    std::cerr << "NSP should be of dimension 1\n";
    throw TerminalException{1};
  }
  MyVector<T> V1 = GetMatrixRow(NSP,0);
  T eNorm_V1 = V1.dot(G2 * V1);
  if (eNorm_V1 != 0) {
    std::cerr << "The orthogonal space of Subspace2 should be an isotropic vector\n";
    throw TerminalException{1};
  }
  // Expanding we get eVect2 = V0 + t V1
  // This gets us eNorm = G2[eVect2] = G2[V0] + 2 t V0.G2.V1
  // or scal0 = t scal1
  T scal0 = eNorm - V0.dot(G2 * V0);
  T scal1 = 2 * V0.dot(G2 * V1);
  if (scal1 == 0) {
    std::cerr << "The coefficient scal1 should be non-zero\n";
    throw TerminalException{1};
  }
  T t = scal0 / scal1;
  MyVector<T> eVect2 = V0 + t * V1;
  //
  std::vector<MyVector<T>> LV1{eVect1};
  std::vector<MyVector<T>> LV2{eVect2};
  MyMatrix<T> eVect1_mat = MatrixFromVectorFamily(LV1);
  MyMatrix<T> eVect2_mat = MatrixFromVectorFamily(LV2);
  MyMatrix<T> Trans1 = Concatenate(Subspace1, eVect1_mat);
  MyMatrix<T> Trans2 = Concatenate(Subspace2, eVect2_mat);
  MyMatrix<T> eEquiv = Inverse(Trans1) * Trans2;
#ifdef DEBUG_LORENTZIAN_LINALG
  MyMatrix<T> InvEquiv = Inverse(eEquiv);
  MyMatrix<T> G1_tr = InvEquiv * G1 * InvEquiv.transpose();
  if (G1_tr != G2) {
    std::cerr << "G1 has not been transposed into G2\n";
    throw TerminalException{1};
  }
  MyMatrix<T> testProd = Subspace1 * eEquiv;
  if (testProd != Subspace2) {
    std::cerr << "Subspace1 is not mapped to Subspace2\n";
    throw TerminalException{1};
  }
#endif
  return eEquiv;
}


/*
  The vector of Subspace1 / Subspace2 are no longer assumed independent
 */
template<typename T>
std::optional<MyMatrix<T>> ExtendOrthogonalIsotropicIsomorphism(MyMatrix<T> const& G1, MyMatrix<T> const& Subspace1, MyMatrix<T> const& G2, MyMatrix<T> const& Subspace2)
{
  int dim = G1.rows();
  std::vector<int> ListRowSelect=TMat_SelectRowCol(Subspace1).ListRowSelect;
  MyMatrix<T> Subspace1_red = SelectRow(Subspace1, ListRowSelect);
  MyMatrix<T> Subspace2_red = SelectRow(Subspace2, ListRowSelect);
  if (RankMat(Subspace2_red) != dim - 1) {
    return {};
  }
  MyMatrix<T> eEquiv = ExtendOrthogonalIsotropicIsomorphism_Basis(G1, Subspace1_red, G2, Subspace2_red);
  if (Subspace1 * eEquiv != Subspace2)
    return {};
  return eEquiv;
}




/*
  For a dimension N, we want to find all the possible integers k such that there exist
  an integer matrix A of order k. The solution is given in
  https://en.wikipedia.org/wiki/Crystallographic_restriction_theorem
  and involves the Psi function
 */
template<typename T>
std::vector<T> GetIntegralMatricesPossibleOrders(T const& N)
{
  auto is_prime=[](T const& x) -> bool {
    if (x == 1)
      return false;
    if (x == 2)
      return true;
    T div = 2;
    while (true) {
      T res = ResInt(x, div);
      if (res == 0)
        return false;
      div += 1;
      if (div * div > x)
        break;
    }
    return true;
  };
  std::vector<T> ListPrime;
  for (T val=2; val <= N+1; val++) {
    bool test = is_prime(val);
    std::cerr << "val=" << val << " test=" << test << "\n";
    if (test)
      ListPrime.push_back(val);
  }
  //
  struct pair {
    T fact;
    T dim_cost;
  };
  struct Desc {
    T prime;
    std::vector<pair> l_pair;
  };
  auto get_pair=[&](T const& eprime, int const& k) -> pair {
    if (eprime == 2 && k == 1)
      return {2, 0};
    T pow1 = MyPow(eprime, k-1);
    T fact = pow1 * eprime;
    T dim_cost = pow1 * (eprime - 1);
    return {fact, dim_cost};
  };
  auto get_l_pair=[&](T const& eprime) -> std::vector<pair> {
    std::vector<pair> l_pair;
    l_pair.push_back({1,0});
    int k = 1;
    while(true) {
      pair epair = get_pair(eprime, k);
      if (epair.dim_cost > N)
        break;
      l_pair.push_back(epair);
      k++;
    }
    return l_pair;
  };
  std::vector<Desc> l_desc;
  for (auto & ePrime : ListPrime)
    l_desc.push_back({ePrime, get_l_pair(ePrime)});
  size_t n_case = 1;
  std::vector<int> VectSiz;
  for (auto eDesc : l_desc) {
    size_t len = eDesc.l_pair.size();
    n_case *= len;
    VectSiz.push_back(len);
  }
  std::cerr << "n_case=" << n_case << "\n";
  std::vector<T> l_order;
  for (auto & V : BlockIterationMultiple(VectSiz)) {
    T tot_dim = 0;
    T order = 1;
    for (size_t iPrime=0; iPrime<ListPrime.size(); iPrime++) {
      int pos = V[iPrime];
      pair epair = l_desc[iPrime].l_pair[pos];
      tot_dim += epair.dim_cost;
      order *= epair.fact;
    }
    if (tot_dim <= N)
      l_order.push_back(order);
  }
  std::sort(l_order.begin(), l_order.end());
  return l_order;
}


template<typename T>
bool is_infinite_order(MyMatrix<T> const& M, size_t const& max_finite_order)
{
  int n = M.rows();
  auto is_identity=[&](MyMatrix<T> const& A) -> bool {
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++) {
        if (i == j && A(i,j) != 1)
          return false;
        if (i != j && A(i,j) != 0)
          return false;
      }
    return true;
  };
  MyMatrix<T> ThePow = M;
  size_t expo = 1;
  for (size_t u=0; u<=max_finite_order; u++) { // We go a little bit over the needed range
    if (is_identity(ThePow)) {
      std::cerr << "is_infinite_order expo=" << expo << "\n";
      return true;
    }
    ThePow *= M;
    expo++;
  }
  return false;
}



template<typename T, typename Tint>
struct LorentzianFinitenessGroupTester {
  LorentzianFinitenessGroupTester(MyMatrix<T> const& _G) : G(_G)
  {
    int dim = G.rows();
    T dim_T = dim;
    std::vector<T> V = GetIntegralMatricesPossibleOrders<T>(dim_T);
    max_finite_order = UniversalScalarConversion<int,T>(V[V.size() - 1]);
    InvariantBasis = IdentityMat<Tint>(dim);
    is_finite = true;
  }
  void GeneratorUpdate(MyMatrix<Tint> const& eP)
  {
    MyMatrix<T> eP_T = UniversalMatrixConversion<T,Tint>(eP);
    MyMatrix<T> G_img = eP_T * G * eP_T.transpose();
    if (G_img != G) {
      std::cerr << "G="; WriteMatrix(std::cerr, G);
      std::cerr << "eP_T="; WriteMatrix(std::cerr, eP_T);
      std::cerr << "G_img="; WriteMatrix(std::cerr, G_img);
      std::cerr << "The matrix eP should leave the quadratic form invariant\n";
      throw TerminalException{1};
    }
#ifdef TIMINGS
    std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
    bool test = is_infinite_order(eP, max_finite_order);
    if (!test) {
      is_finite = false;
    }
#ifdef TIMINGS
    std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
    std::cerr << "Timing |is_finite_order|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
    MyMatrix<Tint> eDiff = InvariantBasis * eP - InvariantBasis;
    if (!IsZeroMatrix(eDiff)) {
      MyMatrix<Tint> NSP = NullspaceIntMat(eDiff);
      if (NSP.rows() == 0) {
        is_finite = false;
        InvariantBasis = MyMatrix<Tint>(0,G.rows());
      } else {
        InvariantBasis = NSP * InvariantBasis;
        MyMatrix<T> InvariantBasis_T = UniversalMatrixConversion<T,Tint>(InvariantBasis);
        MyMatrix<T> Ginv = InvariantBasis_T * G * InvariantBasis_T.transpose();
        DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(Ginv);
        if (DiagInfo.nbMinus == 0) {
          is_finite = false;
        }
      }
    }
#ifdef TIMINGS
    std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
    std::cerr << "Timing |InvariantSpace|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif
  }
  bool get_finiteness_status() const
  {
    return is_finite;
  }
  std::string get_infos() const
  {
    return std::string("(dim=") + std::to_string(InvariantBasis.rows()) + "/" + std::to_string(G.rows()) + ")";
  }
private:
  MyMatrix<T> G;
  MyMatrix<Tint> InvariantBasis;
  size_t max_finite_order;
  bool is_finite;
};




#endif
