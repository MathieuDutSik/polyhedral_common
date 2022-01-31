#ifndef DEFINE_LORENTZIAN_LINALG_H
#define DEFINE_LORENTZIAN_LINALG_H

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




#endif
