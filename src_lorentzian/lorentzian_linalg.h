// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LORENTZIAN_LORENTZIAN_LINALG_H_
#define SRC_LORENTZIAN_LORENTZIAN_LINALG_H_

#include "COMB_Combinatorics.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "POLY_cddlib.h"
#include "Positivity.h"
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#ifdef DEBUG
#define DEBUG_LORENTZIAN_LINALG
#endif

#ifdef TIMINGS
#define TIMINGS_LORENTZIAN_LINALG
#endif

template<typename T>
std::optional<std::string> ReasonNonLorentzian(MyMatrix<T> const& G) {
  DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(G);
  if (DiagInfo.nbZero != 0) {
    std::string reason_non_lorentzian = "matrix has non-zero kernel";
    return reason_non_lorentzian;
  }
  int nbMinus = DiagInfo.nbMinus;
  int nbPlus = DiagInfo.nbPlus;
  if (nbMinus != 1) {
    std::string reason_non_lorentzian = "Signature is (" + std::to_string(nbMinus) + "," + std::to_string(nbPlus) + ") but it should be (1,n)";
    return reason_non_lorentzian;
  }
  return {};
}

/*
  A few linear algebra stuff used for the lorentzian computations
 */
template <typename T> void TestLorentzianity(MyMatrix<T> const &G) {
  std::optional<std::string> opt = ReasonNonLorentzian(G);
  if (opt) {
    std::cerr << "LORLIN: G=\n";
    WriteMatrix(std::cerr, G);
    std::cerr << "LORLIN: Reason for non-lorentzianity=" << *opt << "\n";
    throw TerminalException{1};
  }
}

/*
  Given a lattice L and a matrix g, find the smallest exponent m such that g^m
  preserves L
 */
template <typename T>
std::pair<size_t, MyMatrix<T>>
GetMatrixExponentSublattice_Kernel(MyMatrix<T> const &g,
                                   MyMatrix<T> const &Latt) {
  int n = Latt.rows();
  auto is_preserving = [&](MyMatrix<T> const &h) -> bool {
    for (int i = 0; i < n; i++) {
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
  while (true) {
    if (is_preserving(h))
      break;
    ord++;
    h = h * g;
  }
  return {ord, h};
}

template <typename T>
size_t GetMatrixExponentSublattice(MyMatrix<T> const &g,
                                   MyMatrix<T> const &Latt) {
  return GetMatrixExponentSublattice_Kernel(g, Latt).first;
}

/*
  Given a lattice L and a matrix g, find the smallest exponent m such that g^m
  preserves L and acts trivially on the set of translation classes.
 */
template <typename T>
size_t GetMatrixExponentSublattice_TrivClass(MyMatrix<T> const &g,
                                             MyMatrix<T> const &Latt) {
  // First compute the power that preserves the lattice L
  std::pair<size_t, MyMatrix<T>> epair =
      GetMatrixExponentSublattice_Kernel(g, Latt);
  size_t const &ord1 = epair.first;
  MyMatrix<T> const &h = epair.second;
  // Now computing the power of the action on he classes
  std::vector<MyVector<T>> ListTrans = ComputeTranslationClasses<T, T>(Latt);
  size_t n_class = ListTrans.size();
  std::vector<size_t> cl_h = GetActionOnClasses(ListTrans, h, Latt);
  auto is_identity = [&](std::vector<size_t> const &cl_k) -> bool {
    for (size_t i = 0; i < n_class; i++)
      if (cl_k[i] != i)
        return false;
    return true;
  };
  std::vector<size_t> cl_pow = cl_h;
  size_t ord2 = 1;
  while (true) {
    if (is_identity(cl_pow))
      break;
    ord2++;
    std::vector<size_t> W(n_class);
    for (size_t i = 0; i < n_class; i++)
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
  The right way to solve the equation is to compute kG = gcd(k1, k2) and a basis
  of the kernel. We thus remap the equation to u1 + c k1 = a1 u2        = a2
  solvability condition becomes u2 in Z.
  c0 = -u1 / k1
  cS = 1/k1
  Solution is c = c0 + h cS
  k = - c0 / cS
 */
template <typename T>
std::optional<MyVector<T>> ResolveLattEquation(MyMatrix<T> const &Latt,
                                               MyVector<T> const &u,
                                               MyVector<T> const &k) {
  std::vector<MyVector<T>> l_v = {u, k};
  MyMatrix<T> eIndep = MatrixFromVectorFamily(l_v);
  MyMatrix<T> IntBasis = IntersectionLattice_VectorSpace(Latt, eIndep);
  std::optional<MyVector<T>> opt_u = SolutionMat(IntBasis, u);
  if (!opt_u) {
    std::cerr << "LORLIN: We failed to find a solution for u\n";
    throw TerminalException{1};
  }
  const MyVector<T> &sol_u = *opt_u;
  T u1 = sol_u(0);
  T u2 = sol_u(1);
  std::optional<MyVector<T>> opt_k = SolutionMat(IntBasis, k);
  if (!opt_k) {
    std::cerr << "LORLIN: We failed to find a solution for k\n";
    throw TerminalException{1};
  }
  const MyVector<T> &sol_k = *opt_k;
  T k1 = sol_k(0);
  T k2 = sol_k(1);
  //
  GCD_int<T> ep = ComputePairGcd(k1, k2);
  T u1_norm = ep.Pmat(0, 0) * u1 + ep.Pmat(1, 0) * u2;
  T u2_norm = ep.Pmat(0, 1) * u1 + ep.Pmat(1, 1) * u2;
  T k1_norm = ep.Pmat(0, 0) * k1 + ep.Pmat(1, 0) * k2;
  T k2_norm = ep.Pmat(0, 1) * k1 + ep.Pmat(1, 1) * k2;
  if (k2_norm != 0) {
    std::cerr << "LORLIN: We should have k2_norm = 0. Likely a bug here\n";
    throw TerminalException{1};
  }
  if (!IsInteger(u2_norm)) {
    // No solution then
    return {};
  }
  //
  T c0 = -u1_norm / k1_norm;
  T cS = 1 / k1_norm;
  T hinp = -c0 / cS;
  T h;
  if (cS > 0) {
    h = UniversalCeilScalarInteger<T, T>(hinp);
    if (hinp == h)
      h += 1;
  } else {
    h = UniversalFloorScalarInteger<T, T>(hinp);
    if (hinp == h)
      h -= 1;
  }
  T c = c0 + h * cS;
  if (c <= 0) {
    std::cerr << "LORLIN: We should have c>0\n";
    throw TerminalException{1};
  }
  return u + c * k;
}

template <typename T, typename Tint>
MyMatrix<T> Get_Pplane(MyMatrix<T> const &G,
                       std::vector<MyVector<Tint>> const &l_ui) {
  int n = G.rows();
  size_t n_root = l_ui.size();
  MyMatrix<T> EquaPplane(n_root, n);
  for (size_t i_root = 0; i_root < n_root; i_root++) {
    MyVector<T> eV = UniversalVectorConversion<T, Tint>(l_ui[i_root]);
    MyVector<T> eP = G * eV;
    AssignMatrixRow(EquaPplane, i_root, eP);
  }
  MyMatrix<T> Pplane = NullspaceTrMat(EquaPplane);
  if (Pplane.rows() != 2) {
    std::cerr << "LORLIN: The dimension should be exactly 2\n";
    std::cerr << "LORLIN: We have |Pplane|=" << Pplane.rows() << "\n";
    throw TerminalException{1};
  }
  return Pplane;
}

template <typename T> struct LatticeProjectionFramework {
  MyMatrix<T> ProjP;
  MyMatrix<T> BasisProj;
  MyMatrix<T> ProjFamily;
  MyMatrix<T> Latt;
  LatticeProjectionFramework(MyMatrix<T> const &G, MyMatrix<T> const &Subspace,
                             MyMatrix<T> const &_Latt)
      : Latt(_Latt) {
    int n = G.rows();
    int dim = Latt.rows();
    ProjP = GetProjectionMatrix(G, Subspace);
    ProjFamily = MyMatrix<T>(dim, n);
    for (int i = 0; i < dim; i++) {
      MyVector<T> eVect = GetMatrixRow(Latt, i);
      MyVector<T> eVectProj = ProjP * eVect;
      AssignMatrixRow(ProjFamily, i, eVectProj);
    }
    BasisProj = GetZbasis(ProjFamily);
  }
  std::optional<MyVector<T>> GetOnePreimage(MyVector<T> const &V) const {
    std::optional<MyVector<T>> opt = SolutionIntMat(ProjFamily, V);
    if (!opt)
      return {};
    MyVector<T> const &eSol = *opt;
    MyVector<T> preImage = Latt.transpose() * eSol;
    return preImage;
  }
};

template <typename T>
std::vector<size_t>
GetFacetOneDomain_ListIdx(std::vector<MyVector<T>> const &l_vect,
                          std::ostream &os) {
  using Tfield = typename overlying_field<T>::field_type;
  int dimSpace = l_vect[0].size();
  if (l_vect.size() < size_t(2 * dimSpace)) {
    std::cerr << "LORLIN: Number of roots should be at least 2 * dimspace = "
              << (2 * dimSpace) << "\n";
    std::cerr << "LORLIN: while |l_vect|=" << l_vect.size() << "\n";
    throw TerminalException{1};
  }
  auto is_corr = [&](MyVector<T> const &w) -> bool {
    for (auto &e_root : l_vect) {
      T scal = e_root.dot(w);
      if (scal == 0)
        return false;
    }
    return true;
  };
  auto get_random_vect = [&]() -> MyVector<T> {
    MyVector<T> w(dimSpace);
    int spr = 1000;
    int tot_spr = 2 * spr + 1;
    while (true) {
      for (int i = 0; i < dimSpace; i++)
        w(i) = random() % tot_spr - spr;
#ifdef DEBUG_LORENTZIAN_LINALG
      os << "LORLIN: get_random_vect. Trying w=" << StringVectorGAP(w) << "\n";
#endif
      if (is_corr(w))
        return w;
    }
  };
  MyVector<T> selVect = get_random_vect();
#ifdef DEBUG_LORENTZIAN_LINALG
  os << "LORLIN: Random splitting vector selVect=" << StringVectorGAP(selVect) << "\n";
#endif
  int n_vect = l_vect.size() / 2;
  MyMatrix<Tfield> EXT(n_vect, 1 + dimSpace);
  std::vector<size_t> list_idx(n_vect);
  size_t pos = 0;
  for (size_t i = 0; i < l_vect.size(); i++) {
    T scal = selVect.dot(l_vect[i]);
    if (scal > 0) {
      list_idx[pos] = i;
      MyVector<T> const &eV = l_vect[i];
      EXT(pos, 0);
      for (int i = 0; i < dimSpace; i++)
        EXT(pos, i + 1) = UniversalScalarConversion<Tfield, T>(eV(i));
      pos++;
    }
  }
  std::vector<int> list_red = cdd::RedundancyReductionClarkson(EXT, os);
  size_t siz = list_red.size();
  std::vector<size_t> l_idx(siz);
  for (size_t i = 0; i < siz; i++) {
    size_t pos = list_idx[list_red[i]];
    l_idx[i] = pos;
  }
  return l_idx;
}

template <typename T>
std::vector<MyVector<T>>
GetFacetOneDomain(std::vector<MyVector<T>> const &l_vect, std::ostream &os) {
  size_t n_vect = l_vect.size();
  MyMatrix<T> Mvect = MatrixFromVectorFamily(l_vect);
  MyMatrix<T> MvectRed = ColumnReduction(Mvect);
  std::vector<MyVector<T>> l_vect_red(n_vect);
  for (size_t i = 0; i < n_vect; i++)
    l_vect_red[i] = GetMatrixRow(MvectRed, i);
  std::vector<size_t> l_idx = GetFacetOneDomain_ListIdx(l_vect_red, os);
  std::vector<MyVector<T>> l_vect_ret;
  for (auto &idx : l_idx)
    l_vect_ret.push_back(l_vect[idx]);
  return l_vect_ret;
}

/*
  We should compute the image of Subspace1 into Subspace2 using the eEquiv.
  Afterwards, we can assume that eEquiv = Id.
  We select another vector v1 in the complement of Subspace1, such that
  {Subspace1, v1} is a Z-basis on Z^n. We hus want to find the vector v2 such
  that G2[v2] = G1[v1] and V_{1i} G1 v1 = V_{2i} G2 v2 The linear equation
  V_{1i} G1 v1 = V_{2i} G2 v2 has a solution set v2 = v2_0 + alpha w2_0 because
  G2 is non-degenerate. Then we have the quadratic equation G2[v2] and two
  possible solutions.
  ----
  The question is that it is possible that there are 3 possibilities:
  (a) No Solution
  (b) 1 solution
  (c) two solutions
  We can have integer solutions or rational solutions or maybe irrational
  solutions. Are irrational solutions possible? I would think that this is not
  possible.
  ----
  If we have a rational solutions, then we can use the usual scheme
  for passing from rational to integral.
  ----
  Actually things are quite beautiful in the case of Subspace1 the orthogonal of
  an isotropic case. We select a vector eVect1 outside of Subspace1 and look for
  its image eVect2 We thus have G1[eVect1] = G2[eVect2] and the equalities V1(i)
  * G1 * eVect1 = V2(i) * G2 * eVect2 This gets us a solution space of eVect2 =
  V0 + t V1 and the vector V1 is isotropic. This gets us G1[eVect1] = G2[eVect2]
  = G2[V0] + 2t V0.dot.V1 and thus the equation system has an unique solution,
  which all turn out to be unexpected.
 */
template <typename T>
MyMatrix<T> LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1_Basis(
    MyMatrix<T> const &G1, MyMatrix<T> const &Subspace1, MyMatrix<T> const &G2,
    MyMatrix<T> const &Subspace2) {
  int dim = G1.rows();
#ifdef DEBUG_LORENTZIAN_LINALG
  auto terminate = [&](std::string const &msg) -> void {
    std::cerr << "LORLIN: G1=\n";
    WriteMatrix(std::cerr, G1);
    std::cerr << "LORLIN: G2=\n";
    WriteMatrix(std::cerr, G2);
    std::cerr << "LORLIN: Subspace1=\n";
    WriteMatrix(std::cerr, Subspace1);
    std::cerr << "LORLIN: Subspace2=\n";
    WriteMatrix(std::cerr, Subspace2);
    std::cerr << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Basis: " << msg
              << "\n";
    throw TerminalException{1};
  };
  if (Subspace1.rows() != dim - 1 || Subspace2.rows() != dim - 1) {
    terminate("Subspace1 and Subspace2 are not of the right dimension");
  }
  MyMatrix<T> S1_G1_S1t = Subspace1 * G1 * Subspace1.transpose();
  MyMatrix<T> S2_G2_S2t = Subspace2 * G2 * Subspace2.transpose();
  if (S1_G1_S1t != S2_G2_S2t) {
    std::cerr << "LORLIN: S1_G1_S1t=\n";
    WriteMatrix(std::cerr, S1_G1_S1t);
    std::cerr << "LORLIN: S2_G2_S2t=\n";
    WriteMatrix(std::cerr, S2_G2_S2t);
    terminate("The S1 / G1 / S2 / G2 are not coherent on input");
  }
#endif
  MyMatrix<T> Compl1 = SubspaceCompletionRational(Subspace1, dim);
#ifdef DEBUG_LORENTZIAN_LINALG
  if (Compl1.rows() != 1) {
    terminate("Compl1 should be of dimension 1");
  }
#endif
  MyVector<T> eVect1 = GetMatrixRow(Compl1, 0);
  T eNorm = eVect1.dot(G1 * eVect1);
  MyMatrix<T> eProd1 = Subspace1 * G1;
  MyMatrix<T> eProd2 = Subspace2 * G2;
  MyVector<T> Vscal = Subspace1 * G1 * eVect1;
  std::optional<MyVector<T>> opt = SolutionMat(TransposedMat(eProd2), Vscal);
#ifdef DEBUG_LORENTZIAN_LINALG
  if (!opt) {
    terminate("The solutioning failed");
  }
#endif
  MyVector<T> const &V0 = *opt;
  MyMatrix<T> NSP = NullspaceTrMat(eProd2);
#ifdef DEBUG_LORENTZIAN_LINALG
  if (NSP.rows() != 1) {
    terminate("NSP should be of dimension 1");
  }
#endif
  MyVector<T> V1 = GetMatrixRow(NSP, 0);
  T eNorm_V1 = V1.dot(G2 * V1);
#ifdef DEBUG_LORENTZIAN_LINALG
  if (eNorm_V1 != 0) {
    terminate(
        "The orthogonal space of Subspace2 should be an isotropic vector");
  }
#endif
  // Expanding we get eVect2 = V0 + t V1
  // This gets us eNorm = G2[eVect2] = G2[V0] + 2 t V0.G2.V1
  // or scal0 = t scal1
  T scal0 = eNorm - V0.dot(G2 * V0);
  T scal1 = 2 * V0.dot(G2 * V1);
#ifdef DEBUG_LORENTZIAN_LINALG
  if (scal1 == 0) {
    terminate("The coefficient scal1 should be non-zero");
  }
#endif
  T t = scal0 / scal1;
  MyVector<T> eVect2 = V0 + t * V1;
  //
  MyMatrix<T> Trans1 = ConcatenateMatVec(Subspace1, eVect1);
  MyMatrix<T> Trans2 = ConcatenateMatVec(Subspace2, eVect2);
  MyMatrix<T> eEquiv = Inverse(Trans1) * Trans2;
#ifdef DEBUG_LORENTZIAN_LINALG
  MyMatrix<T> InvEquiv = Inverse(eEquiv);
  MyMatrix<T> G1_tr = InvEquiv * G1 * InvEquiv.transpose();
  if (G1_tr != G2) {
    terminate("G1 has not been transposed into G2");
  }
  MyMatrix<T> testProd = Subspace1 * eEquiv;
  if (testProd != Subspace2) {
    terminate("Subspace1 is not mapped to Subspace2");
  }
#endif
  return eEquiv;
}

/*
  The vector of Subspace1 / Subspace2 are no longer assumed independent
 */
template <typename T>
std::optional<MyMatrix<T>> LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1(
    MyMatrix<T> const &G1, MyMatrix<T> const &Subspace1, MyMatrix<T> const &G2,
    MyMatrix<T> const &Subspace2) {
  int dim = G1.rows();
  std::vector<int> ListRowSelect = TMat_SelectRowCol(Subspace1).ListRowSelect;
  MyMatrix<T> Subspace1_red = SelectRow(Subspace1, ListRowSelect);
  MyMatrix<T> Subspace2_red = SelectRow(Subspace2, ListRowSelect);
  if (RankMat(Subspace2_red) != dim - 1) {
    return {};
  }
  MyMatrix<T> eEquiv = LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1_Basis(
      G1, Subspace1_red, G2, Subspace2_red);
  if (Subspace1 * eEquiv != Subspace2)
    return {};
  return eEquiv;
}

template <typename T> struct SolutionSpecial {
  std::vector<MyMatrix<T>> BasisKernel;
  MyMatrix<T> eSol_mat;
};

// Resolution of B = X A + A^T X^T
// b_{ij} = sum_k x_{ik} a_{kj} + a_{ki} x_{jk}
template <typename T>
SolutionSpecial<T> SpecialEquationSolving(MyMatrix<T> const &Amat,
                                          MyMatrix<T> const &Bmat) {
  int dim = Bmat.rows();
  MyMatrix<T> TheMat = ZeroMatrix<T>(dim * dim, dim * dim);
  MyVector<T> Bvect = ZeroVector<T>(dim * dim);
  auto f = [&](int i, int j) -> int { return i + dim * j; };
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      int u = f(i, j);
      Bvect(u) = Bmat(i, j);
      for (int k = 0; k < dim; k++) {
        // contribudion x_{ik} a_{kj}
        int v1 = f(i, k);
        TheMat(v1, u) += Amat(k, j);
        // contribution a_{ki} x_{jk}
        int v2 = f(j, k);
        TheMat(v2, u) += Amat(k, i);
      }
    }
  }
  auto f_getmat = [&](MyVector<T> const &eVect) -> MyMatrix<T> {
    MyMatrix<T> eMat(dim, dim);
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        eMat(i, j) = eVect(f(i, j));
      }
    }
    return eMat;
  };
  std::optional<MyVector<T>> opt = SolutionMat(TheMat, Bvect);
  MyVector<T> eSol_vect = unfold_opt(opt, "getting eSol_vect");
  MyMatrix<T> eSol_mat = f_getmat(eSol_vect);
#ifdef DEBUG_LORENTZIAN_LINALG
  MyMatrix<T> SumMat =
      eSol_mat * Amat + Amat.transpose() * eSol_mat.transpose();
  if (Bmat != SumMat) {
    std::cerr << "LORLIN: Failed to find the correct solution 1\n";
    throw TerminalException{1};
  }
#endif
  std::vector<MyMatrix<T>> BasisKernel;
  MyMatrix<T> NSP = NullspaceMat(TheMat);
  int dimNSP = NSP.rows();
  for (int u = 0; u < dimNSP; u++) {
    MyVector<T> eVect = GetMatrixRow(NSP, u);
    MyMatrix<T> eMat = f_getmat(eVect);
#ifdef DEBUG_LORENTZIAN_LINALG
    MyMatrix<T> SumMat2 = eMat * Amat + Amat.transpose() * eMat.transpose();
    if (!IsZeroMatrix(SumMat2)) {
      std::cerr << "LORLIN: Failed to find the correct solution 2\n";
      throw TerminalException{1};
    }
#endif
    BasisKernel.push_back(eMat);
  }
  return {BasisKernel, eSol_mat};
}

template <typename T>
std::vector<MyMatrix<T>>
IntegralSpaceSaturation_Matrix(std::vector<MyMatrix<T>> const &ListM,
                               int const &n) {
  int dim = n * n;
  int nMat = ListM.size();
  MyMatrix<T> BigM(nMat, dim);
  for (int u = 0; u < nMat; u++) {
    MyVector<T> eV = MatrixToVector(ListM[u]);
    AssignMatrixRow(BigM, u, eV);
  }
  MyMatrix<T> BigM_sat = IntegralSpaceSaturation(BigM);
  std::vector<MyMatrix<T>> RetListM;
  for (int u = 0; u < nMat; u++) {
    MyVector<T> eV = GetMatrixRow(BigM_sat, u);
    MyMatrix<T> eM = VectorToMatrix(eV, n);
    RetListM.push_back(eM);
  }
  return RetListM;
}

template <typename T> struct LORENTZ_ExtendOrthogonalIsotropicIsomorphism {
public:
  MyMatrix<T> G1;
  MyMatrix<T> Subspace1;
  MyMatrix<T> G2;
  MyMatrix<T> Subspace2;
  int dim;
  int rnk;
  MyMatrix<T> Trans1Inv;
  MyMatrix<T> ListVectCand2;
  T denomSol;
  SolutionSpecial<T> TheRec;
  std::vector<MyMatrix<T>> ListEquiv_terms1;
  MyMatrix<T> NSP1;
  MyMatrix<T> NSP2;
  std::ostream &os;
  LORENTZ_ExtendOrthogonalIsotropicIsomorphism(MyMatrix<T> const &_G1,
                                               MyMatrix<T> const &_Subspace1,
                                               MyMatrix<T> const &_G2,
                                               MyMatrix<T> const &_Subspace2,
                                               std::ostream &_os)
      : G1(_G1), Subspace1(_Subspace1), G2(_G2), Subspace2(_Subspace2),
        dim(G1.rows()), rnk(Subspace1.rows()), os(_os) {
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, beginning\n";
    if (rnk != RankMat(Subspace1)) {
      std::cerr
          << "Inconsistent input: We should have Subspace1 of full rank\n";
      throw TerminalException{1};
    }
#endif
    int n = G2.rows();
    // Checking the input
    MyMatrix<T> eProd1 = Subspace1 * G1;
    MyMatrix<T> eProd2 = Subspace2 * G2;
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "eProd1 / eProd2\n";
#endif
    NSP1 = NullspaceTrMat(eProd1);
    NSP2 = NullspaceTrMat(eProd2);
#ifdef DEBUG_LORENTZIAN_LINALG
    MyMatrix<T> ProdMat1 = NSP1 * G1 * NSP1.transpose();
    if (!IsZeroMatrix(ProdMat1)) {
      std::cerr << "LORLIN: Inconsistent input: ProdMat1 should be equal to zero\n";
      throw TerminalException{1};
    }
    MyMatrix<T> ProdMat2 = NSP2 * G2 * NSP2.transpose();
    if (!IsZeroMatrix(ProdMat2)) {
      std::cerr << "LORLIN: Inconsistent input: ProdMat1 should be equal to zero\n";
      throw TerminalException{1};
    }
    MyMatrix<T> SMat1 = Subspace1 * G1 * Subspace1.transpose();
    MyMatrix<T> SMat2 = Subspace2 * G2 * Subspace2.transpose();
    if (SMat1 != SMat2) {
      std::cerr << "LORLIN: Inconsistent input: SMat1 should be equal to SMat2\n";
      throw TerminalException{1};
    }
#endif
    // Now doing the computation side of the job
    MyMatrix<T> TheCompl1 = SubspaceCompletionRational(Subspace1, dim);
    MyMatrix<T> Trans1 = Concatenate(Subspace1, TheCompl1);
    Trans1Inv = Inverse(Trans1);
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "Trans1Inv\n";
#endif
    MyMatrix<T> MatScal1 = TheCompl1 * G1 * TheCompl1.transpose();
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "MatScal1\n";
#endif
    int dimCompl1 = TheCompl1.rows();
    ListVectCand2 = MyMatrix<T>(dimCompl1, n);
    MyMatrix<T> eProd2tr = eProd2.transpose();
    for (int u = 0; u < dimCompl1; u++) {
      MyVector<T> eVect1 = GetMatrixRow(TheCompl1, u);
#ifdef DEBUG_LORENTZIAN_LINALG
      os << "LORLIN: We have eVect1=" << StringVectorGAP(eVect1) << "\n";
#endif
      MyVector<T> Vscal = eProd1 * eVect1;
      std::optional<MyVector<T>> opt = SolutionMat(eProd2tr, Vscal);
      MyVector<T> eVectCand2 = unfold_opt(opt, "getting eVectCand2");
#ifdef DEBUG_LORENTZIAN_LINALG
      os << "LORLIN: We have eVectCand2=" << StringVectorGAP(eVectCand2)
         << "\n";
      MyVector<T> LScal1 = eProd1 * eVect1;
      MyVector<T> LScal2 = eProd2 * eVectCand2;
      if (LScal1 != LScal2) {
        std::cerr << "LORLIN: Inconsistency for LScal1 = LScal2\n";
        throw TerminalException{1};
      }
#endif
      AssignMatrixRow(ListVectCand2, u, eVectCand2);
    }
#ifdef DEBUG_LORENTZIAN_LINALG_DISABLE
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "ListVectCand2\n";
    os << "LORLIN: G2=\n";
    WriteMatrix(os, G2);
    os << "LORLIN: ListVectCand2=\n";
    WriteMatrix(os, ListVectCand2);
#endif
    // The solutions are written as
    // eVect2 = eVectCand2 + c_vect * NSP2
    // Putting together this gets ListVect2 = TheCompl2 = ListVectCand2 + c_Mat
    // * NSP2 MatScal2 = (ListVectCand2 + c_Mat * NSP2) * G2 * (NSP2^T * c_Mat^T
    // + ListVectCand2^T)
    //          = ListVectCand2 * G2 * ListVectCand2^T + c_Mat * NSP2 * G2 *
    //          ListVectCand2^T + ListVectCand2 * G2 * NSP2^T * c_Mat^T
    // c_vect is a vector of length dimCompl which gets into
    // Put all together, this gets us c_Mat a matrix of size (dimCompl,
    // dimCompl) The equation that we get is thus of the form B = X A + A^T X^T
    // The equation is underdefined. This is because B is symmetric and so we
    // have n(n+1)/2 equations for n^2 unknowns. The space NSP2 is uniquely
    // defined as the set of isotropic vectors in the space. It is defined as a
    // Kernel.
    MyMatrix<T> MatScalCand2 = ListVectCand2 * G2 * ListVectCand2.transpose();
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "MatScalCand2\n";
#endif
    MyMatrix<T> LambdaMat = NSP2 * G2 * ListVectCand2.transpose();
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "LambdaMat\n";
#endif
    denomSol = GetDenominatorQuotientSolution(LambdaMat);
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, denomSol="
       << denomSol << "\n";
#endif
    MyMatrix<T> DiffScal = MatScal1 - MatScalCand2;
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "DiffScal\n";
#endif
    TheRec = SpecialEquationSolving(LambdaMat, DiffScal);
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "TheRec ( SpecialEquationSolving )\n";
#endif
    for (auto &eMat : TheRec.BasisKernel) {
      MyMatrix<T> Null_mat = ZeroMatrix<T>(rnk, dim);
      MyMatrix<T> Prod_mat = eMat * NSP2;
      MyMatrix<T> Cont_mat = Concatenate(Null_mat, Prod_mat);
      MyMatrix<T> Ins_mat = Trans1Inv * Cont_mat;
      ListEquiv_terms1.push_back(Ins_mat);
    }
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "ListEquiv_terms1\n";
#endif
  }
#ifdef DEBUG_LORENTZIAN_LINALG
  void check_transformation(MyMatrix<T> const &eEq) {
    MyMatrix<T> eEqInv = Inverse(eEq);
    MyMatrix<T> G1_tr = eEqInv * G1 * eEqInv.transpose();
    if (G1_tr != G2) {
      std::cerr << "LORLIN: G1 was not mapped to G2\n";
      throw TerminalException{1};
    }
    MyMatrix<T> Subspace1_img = Subspace1 * eEq;
    if (Subspace1_img != Subspace2) {
      std::cerr << "LORLIN: Subspace1 is not mapped to Subspace2\n";
      throw TerminalException{1};
    }
  }
#endif
  MyMatrix<T> get_one_transformation() {
#ifdef DEBUG_LORENTZIAN_LINALG_DISABLE
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, "
          "get_one_transformation beginning\n";
    os << "LORLIN: eSol_mat=\n";
    WriteMatrix(os, TheRec.eSol_mat);
    os << "LORLIN: NSP2=\n";
    WriteMatrix(os, NSP2);
    os << "LORLIN: ListVectCand2=\n";
    WriteMatrix(os, ListVectCand2);
#endif
    MyMatrix<T> TheCompl2 = ListVectCand2 + TheRec.eSol_mat * NSP2;
    MyMatrix<T> Trans2 = Concatenate(Subspace2, TheCompl2);
    MyMatrix<T> eEquiv0 = Trans1Inv * Trans2;
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "eEquiv0\n";
#endif
    // The matrix is expressed as eEquiv0 + alpha1 ListEquiv_terms[1] + ..... +
    // alphaN ListEquiv_terms[N]
    MyMatrix<T> RetSol =
        EliminateSuperfluousPrimeDenominators_Matrix(eEquiv0, ListEquiv_terms1);
#ifdef DEBUG_LORENTZIAN_LINALG
    os << "LORLIN: LORENTZ_ExtendOrthogonalIsotropicIsomorphism, We have "
          "RetSol\n";
#endif
#ifdef DEBUG_LORENTZIAN_LINALG
    check_transformation(RetSol);
#endif
    return RetSol;
  }
  std::vector<MyMatrix<T>> get_kernel_generating_set(T const &d) {
#ifdef DEBUG_LORENTZIAN_LINALG
    if (G1 != G2 || Subspace1 != Subspace2) {
      std::cerr << "LORLIN: We should have G1=G2 and Subspace1=Subspace2 in order for "
                   "kernel to make sense\n";
      throw TerminalException{1};
    }
#endif
    std::vector<MyMatrix<T>> ListEquiv_terms2 =
        IntegralSpaceSaturation_Matrix(ListEquiv_terms1, dim);
    std::vector<MyMatrix<T>> ListEquiv_terms3;
    for (auto &eM : ListEquiv_terms2) {
      MyMatrix<T> eGen = IdentityMat<T>(dim) + eM / d;
#ifdef DEBUG_LORENTZIAN_LINALG
      check_transformation(eGen);
#endif
      ListEquiv_terms3.push_back(eGen);
    }
    return ListEquiv_terms3;
  }
};

/*
  For a dimension N, we want to find all the possible integers k such that there
  exist an integer matrix A of order k. The solution is given in
  https://en.wikipedia.org/wiki/Crystallographic_restriction_theorem
  and involves the Psi function
 */
template <typename T>
std::vector<T> GetIntegralMatricesPossibleOrders(T const &N) {
  auto is_prime = [](T const &x) -> bool {
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
  for (T val = 2; val <= N + 1; val++) {
    bool test = is_prime(val);
#ifdef DEBUG_LORENTZIAN_LINALG
    std::cerr << "LORLIN: val=" << val << " test=" << test << "\n";
#endif
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
  auto get_pair = [&](T const &eprime, int const &k) -> pair {
    if (eprime == 2 && k == 1)
      return {2, 0};
    T pow1 = MyPow(eprime, k - 1);
    T fact = pow1 * eprime;
    T dim_cost = pow1 * (eprime - 1);
    return {fact, dim_cost};
  };
  auto get_l_pair = [&](T const &eprime) -> std::vector<pair> {
    std::vector<pair> l_pair;
    l_pair.push_back({1, 0});
    int k = 1;
    while (true) {
      pair epair = get_pair(eprime, k);
      if (epair.dim_cost > N)
        break;
      l_pair.push_back(epair);
      k++;
    }
    return l_pair;
  };
  std::vector<Desc> l_desc;
  for (auto &ePrime : ListPrime)
    l_desc.push_back({ePrime, get_l_pair(ePrime)});
#ifdef DEBUG_LORENTZIAN_LINALG
  size_t n_case = 1;
#endif
  std::vector<int> VectSiz;
  for (auto eDesc : l_desc) {
    size_t len = eDesc.l_pair.size();
#ifdef DEBUG_LORENTZIAN_LINALG
    n_case *= len;
#endif
    VectSiz.push_back(len);
  }
#ifdef DEBUG_LORENTZIAN_LINALG
  std::cerr << "LORLIN: n_case=" << n_case << "\n";
#endif
  std::vector<T> l_order;
  for (auto &V : BlockIterationMultiple(VectSiz)) {
    T tot_dim = 0;
    T order = 1;
    for (size_t iPrime = 0; iPrime < ListPrime.size(); iPrime++) {
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

template <typename T>
bool is_infinite_order(MyMatrix<T> const &M, size_t const &max_finite_order) {
  int n = M.rows();
  auto is_identity = [&](MyMatrix<T> const &A) -> bool {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        if (i == j && A(i, j) != 1)
          return false;
        if (i != j && A(i, j) != 0)
          return false;
      }
    return true;
  };
  MyMatrix<T> ThePow = M;
#ifdef DEBUG_LORENTZIAN_LINALG
  size_t expo = 1;
#endif
  for (size_t u = 0; u <= max_finite_order; u++) {
    // We go a little bit over the needed range
    if (is_identity(ThePow)) {
#ifdef DEBUG_LORENTZIAN_LINALG
      std::cerr << "LORLIN: is_infinite_order expo=" << expo << "\n";
#endif
      return true;
    }
    ThePow *= M;
#ifdef DEBUG_LORENTZIAN_LINALG
    expo++;
#endif
  }
  return false;
}

template <typename T, typename Tint> struct LorentzianFinitenessGroupTester {
  LorentzianFinitenessGroupTester(MyMatrix<T> const &_G) : G(_G) {
    int dim = G.rows();
    T dim_T = dim;
    std::vector<T> V = GetIntegralMatricesPossibleOrders<T>(dim_T);
    max_finite_order = UniversalScalarConversion<int, T>(V[V.size() - 1]);
    InvariantBasis = IdentityMat<Tint>(dim);
    is_finite = true;
  }
  void GeneratorUpdate(MyMatrix<Tint> const &eP) {
    MyMatrix<T> eP_T = UniversalMatrixConversion<T, Tint>(eP);
    MyMatrix<T> G_img = eP_T * G * eP_T.transpose();
    if (G_img != G) {
      std::cerr << "LORLIN: G=";
      WriteMatrix(std::cerr, G);
      std::cerr << "LORLIN: eP_T=";
      WriteMatrix(std::cerr, eP_T);
      std::cerr << "LORLIN: G_img=";
      WriteMatrix(std::cerr, G_img);
      std::cerr << "LORLIN: The matrix eP should leave the quadratic form invariant\n";
      throw TerminalException{1};
    }
#ifdef TIMINGS_LORENTZIAN_LINALG
    SingletonTime time1;
#endif
    bool test = is_infinite_order(eP, max_finite_order);
    if (!test) {
      is_finite = false;
    }
#ifdef TIMINGS_LORENTZIAN_LINALG
    SingletonTime time2;
    std::cerr << "|LORLIN: is_finite_order|=" << ms(time1, time2) << "\n";
#endif
    MyMatrix<Tint> eDiff = InvariantBasis * eP - InvariantBasis;
    if (!IsZeroMatrix(eDiff)) {
      MyMatrix<Tint> NSP = NullspaceIntMat(eDiff);
      if (NSP.rows() == 0) {
        is_finite = false;
        InvariantBasis = MyMatrix<Tint>(0, G.rows());
      } else {
        InvariantBasis = NSP * InvariantBasis;
        DiagSymMat<T> DiagInfo = get_diag_info();
        if (DiagInfo.nbMinus == 0) {
          is_finite = false;
        }
      }
    }
#ifdef TIMINGS_LORENTZIAN_LINALG
    SingletonTime time3;
    std::cerr << "|LORLIN: InvariantSpace|=" << ms(time2, time3) << "\n";
#endif
  }
  DiagSymMat<T> get_diag_info() const {
    MyMatrix<T> InvariantBasis_T =
        UniversalMatrixConversion<T, Tint>(InvariantBasis);
    MyMatrix<T> Ginv = InvariantBasis_T * G * InvariantBasis_T.transpose();
    DiagSymMat<T> DiagInfo = DiagonalizeSymmetricMatrix(Ginv);
    return DiagInfo;
  }
  size_t get_max_finite_order() const { return max_finite_order; }
  MyMatrix<Tint> const &get_invariant_basis() const { return InvariantBasis; }
  bool get_finiteness_status() const { return is_finite; }
  std::string get_infos() const {
    return std::string("(dim=") + std::to_string(InvariantBasis.rows()) + "/" +
           std::to_string(G.rows()) + ")";
  }

private:
  MyMatrix<T> G;
  MyMatrix<Tint> InvariantBasis;
  size_t max_finite_order;
  bool is_finite;
};

// clang-format off
#endif  // SRC_LORENTZIAN_LORENTZIAN_LINALG_H_
// clang-format on
